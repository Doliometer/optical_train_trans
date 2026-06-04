"""
Growth curves with pad-aware sky annulus.

For each star:
  1. Collect centroids from N_STACK frames to measure actual shift extents.
  2. Derive the safe radius: largest circle from the centroid that avoids
     zero-padded borders introduced by shift-and-add registration.
  3. Place sky annulus just inside that safe radius.
  4. Fit Moffat β=3.0 to the stacked radial profile; flag the star as
     unreliable for T_tel if sky_inner < MIN_FWHM_SCALE × FWHM.
  5. Plot encircled-energy fraction and annular flux.

Usage
-----
  python3 growth_moffat.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import glob
from astropy.io import fits
from pathlib import Path
from survey_transmissivity import find_star_roots, find_filter_dir
from telescope_transmissivity import (
    register_frames, pad_safe_radius, fit_fwhm_moffat,
)

BASE_DIR        = Path('.')
BETA            = 3.0
SKY_MARGIN      = 30    # px inside safe radius for sky inner edge
SKY_WIDTH       = 25    # px annulus width
N_STACK         = 20    # frames to stack
R_MAX           = 450   # max radius for growth curve [px]
R_STEP          = 5
MIN_FWHM_SCALE  = 3.0   # flag star if sky_inner < this × FWHM
PLATE_SCALE_MAS = 18.63

# Stars selected previously (moffat_fit.py, seed=42)
TARGET_LABELS = ['_c1_ref_G2P1068R', '_c1_ref_G2P1060R', '_c1_ref_G2P1110R']




# ---------------------------------------------------------------------------
# Growth curve
# ---------------------------------------------------------------------------

def growth_curve(stack, cx, cy, sky_inner, sky_outer,
                 r_max=R_MAX, r_step=R_STEP):
    ny, nx  = stack.shape
    yy, xx  = np.ogrid[:ny, :nx]
    r2      = (xx - cx)**2 + (yy - cy)**2

    sky_mask    = (r2 > sky_inner**2) & (r2 <= sky_outer**2)
    sky_per_pix = np.median(stack[sky_mask]) if sky_mask.sum() >= 20 else 0.0

    radii   = np.arange(r_step, r_max + r_step, r_step)
    cumflux = np.array([stack[r2 <= r**2].sum() - sky_per_pix * (r2 <= r**2).sum()
                        for r in radii], dtype=float)
    return radii, cumflux, sky_per_pix


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    base      = BASE_DIR.resolve()
    all_roots = dict(find_star_roots(base))

    fig, axes = plt.subplots(2, len(TARGET_LABELS),
                             figsize=(6 * len(TARGET_LABELS), 9))

    for col, label in enumerate(TARGET_LABELS):
        if label not in all_roots:
            print(f'{label}: not found — skipping.')
            continue

        root  = all_roots[label]
        fdir  = find_filter_dir(root, 'Luminance')
        if fdir is None:
            print(f'{label}: no Luminance dir.')
            continue
        fits_files = sorted(glob.glob(str(fdir / '*.fits')))
        if not fits_files:
            print(f'{label}: no FITS files.')
            continue

        print(f'\n{label}: registering {N_STACK} frames …', end=' ', flush=True)
        stack, cx, cy, dx, dy, ny, nx = register_frames(fits_files, N_STACK)
        print(f'done.  centroid=({cx},{cy})')
        print(f'  shift extents: dx=[{dx.min():+.0f},{dx.max():+.0f}]  '
              f'dy=[{dy.min():+.0f},{dy.max():+.0f}]')

        safe_r = pad_safe_radius(cx, cy, dx, dy, ny, nx)
        print(f'  safe radius={safe_r} px')

        sky_inner = safe_r - SKY_MARGIN
        sky_outer = safe_r - SKY_MARGIN + SKY_WIDTH
        # ensure outer doesn't exceed safe_r
        sky_outer = min(sky_outer, safe_r - 2)
        print(f'  Sky annulus: {sky_inner}–{sky_outer} px  '
              f'(~{np.pi*(sky_outer**2 - sky_inner**2):.0f} px²)')

        fwhm = fit_fwhm_moffat(stack, cx, cy, ny, nx)
        if fwhm is None:
            print(f'  Moffat fit failed.')
            fwhm_str = 'fit failed'
            reliable = False
        else:
            fwhm_as  = fwhm * PLATE_SCALE_MAS / 1000.0
            reliable = sky_inner >= MIN_FWHM_SCALE * fwhm
            flag     = '' if reliable else f'  [WARN: sky_inner={sky_inner:.0f} < {MIN_FWHM_SCALE:.0f}×FWHM={MIN_FWHM_SCALE*fwhm:.0f}]'
            fwhm_str = f'{fwhm:.0f}px={fwhm_as:.2f}"'
            print(f'  FWHM = {fwhm_str}  (β={BETA}){flag}')
            print(f'  Reliable for T_tel: {"yes" if reliable else "no — sky may be halo-contaminated"}')

        radii, cumflux, sky_per_pix = growth_curve(
            stack, cx, cy, sky_inner, sky_outer)
        peak = np.nanmax(cumflux)
        frac = cumflux / peak if peak > 0 else cumflux
        print(f'  Sky/px={sky_per_pix:.1f} ADU   peak cumflux={peak:,.0f} ADU')

        above = radii[frac >= 0.99]
        r99   = float(above[0]) if len(above) else np.nan
        print(f'  99% EE radius: {r99:.0f} px' if np.isfinite(r99) else '  99% EE: not reached')

        title_flag = '' if reliable else ' [unreliable]'
        ax0 = axes[0, col]
        ax0.plot(radii, frac, 'b-', lw=1.5)
        ax0.axhline(0.99, color='gray', ls=':', lw=1)
        ax0.axvspan(sky_inner, sky_outer, alpha=0.15, color='red', label=f'sky {sky_inner}–{sky_outer}px')
        if np.isfinite(r99):
            ax0.axvline(r99, color='purple', ls=':', lw=1, label=f'99%EE={r99:.0f}px')
        ax0.set_title(f'{label}\nFWHM={fwhm_str}{title_flag}', fontsize=8)
        ax0.set_xlabel('Aperture radius [px]')
        ax0.set_ylabel('Fraction of total flux')
        ax0.set_ylim(-0.15, 1.1)
        ax0.legend(fontsize=7)

        ax1 = axes[1, col]
        dflux = np.diff(cumflux, prepend=cumflux[0])
        dmax  = np.nanmax(np.abs(dflux))
        ax1.plot(radii, dflux / dmax if dmax > 0 else dflux, 'g-', lw=1.5)
        ax1.axhline(0, color='gray', ls=':', lw=0.8)
        ax1.axvspan(sky_inner, sky_outer, alpha=0.15, color='red', label=f'sky {sky_inner}–{sky_outer}px')
        ax1.set_xlabel('Aperture radius [px]')
        ax1.set_ylabel('Normalised flux per annulus')
        ax1.legend(fontsize=7)

    plt.suptitle('Growth curves: pad-aware sky annulus  (β=3.0)', y=1.01)
    plt.tight_layout()
    outpath = Path('growth_moffat.png')
    plt.savefig(outpath, dpi=120, bbox_inches='tight')
    print(f'\nPlot saved: {outpath}')


if __name__ == '__main__':
    main()
