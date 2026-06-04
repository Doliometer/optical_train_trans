"""
Encircled energy growth curve on shift-and-add registered stack.

For each star / filter, produces:
  - Cumulative ADU vs aperture radius (background-subtracted)
  - The radius at which growth flattens to within CONVERGE_TOL of its peak
  - A text table and an optional plot

Usage
-----
  python3 growth_curve.py                 # all stars, both filters
  python3 growth_curve.py --plot          # also save PNG plots
  python3 growth_curve.py --star G2P1060R # one star only
"""

import argparse
import numpy as np
from astropy.io import fits
from pathlib import Path
import glob
import warnings
warnings.filterwarnings('ignore')

from telescope_transmissivity import (
    R_SKY_INNER, R_SKY_OUTER,
    smooth, find_star_in_frame,
)
from survey_transmissivity import find_star_roots, find_filter_dir, FILTER_CONFIG

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

BASE_DIR     = Path('.')
N_STACK      = None          # frames to stack (None = all)
R_MAX        = 450           # maximum radius to probe [pixels]
R_STEP       = 5             # radii sampled every R_STEP pixels
# Convergence: radius where cumulative flux first reaches CONVERGE_TOL of
# the maximum cumulative flux within R_MAX.
CONVERGE_TOL = 0.99          # 99% encircled energy radius

# ---------------------------------------------------------------------------

def growth_curve_from_stack(fits_files, n_max=None, r_max=R_MAX, r_step=R_STEP,
                             r_sky_inner=R_SKY_INNER, r_sky_outer=R_SKY_OUTER):
    """
    Build a shift-and-add stack then measure the growth curve.

    Returns
    -------
    radii : ndarray   aperture radii sampled [pixels]
    cumflux : ndarray background-subtracted cumulative ADU at each radius
    cx, cy : int      star centroid in the stack
    """
    from scipy.ndimage import shift as ndimage_shift

    files = fits_files[:n_max] if n_max else fits_files
    n = len(files)

    # First pass: collect centroids
    centroids = []
    for f in files:
        img = fits.getdata(f).astype(np.float32)
        cx, cy = find_star_in_frame(img)
        centroids.append((cx, cy))

    ny, nx = fits.getdata(files[0]).shape
    cx_ref, cy_ref = centroids[0]

    # Second pass: register and accumulate
    accum = np.zeros((ny, nx), dtype=np.float64)
    for f, (cx, cy) in zip(files, centroids):
        img = fits.getdata(f).astype(np.float32)
        shifted = ndimage_shift(img, shift=(cy_ref - cy, cx_ref - cx),
                                order=1, mode='constant', cval=0.0)
        accum += shifted
    stack = accum / len(files)

    # Sky background from annulus
    yy, xx = np.ogrid[:ny, :nx]
    r2_grid = (xx - cx_ref)**2 + (yy - cy_ref)**2
    sky_mask = (r2_grid > r_sky_inner**2) & (r2_grid <= r_sky_outer**2)
    sky_per_pix = np.median(stack[sky_mask])

    # Growth curve
    radii = np.arange(r_step, r_max + r_step, r_step)
    cumflux = np.zeros(len(radii))
    for i, r in enumerate(radii):
        ap_mask = r2_grid <= r**2
        n_ap = ap_mask.sum()
        cumflux[i] = stack[ap_mask].sum() - sky_per_pix * n_ap

    return radii, cumflux, cx_ref, cy_ref, sky_per_pix


def convergence_radius(radii, cumflux, tol=CONVERGE_TOL):
    """Return the smallest radius enclosing `tol` fraction of peak flux."""
    peak = np.nanmax(cumflux)
    if peak <= 0:
        return np.nan
    above = radii[cumflux >= tol * peak]
    return float(above[0]) if len(above) else np.nan


def print_growth_table(radii, cumflux, r_conv):
    peak = np.nanmax(cumflux)
    print(f'  {"r(px)":>6}  {"cum ADU":>14}  {"fraction":>8}')
    print(f'  {"------":>6}  {"-------":>14}  {"--------":>8}')
    prev = 0
    for r, f in zip(radii, cumflux):
        frac = f / peak if peak > 0 else 0
        marker = ' <-- 99% EE' if abs(r - r_conv) < 1 else ''
        # Print every other row but always print near the convergence radius
        if r % 20 == 0 or abs(r - r_conv) <= 10:
            print(f'  {r:>6.0f}  {f:>14,.0f}  {frac:>8.4f}{marker}')


def save_plot(radii, cumflux, r_conv, title, outpath):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        peak = np.nanmax(cumflux)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(title)

        ax1.plot(radii, cumflux / peak, 'b-')
        ax1.axvline(r_conv, color='r', linestyle='--',
                    label=f'99% EE  r={r_conv:.0f} px')
        ax1.axhline(0.99, color='gray', linestyle=':')
        ax1.set_xlabel('Aperture radius [pixels]')
        ax1.set_ylabel('Fraction of total flux')
        ax1.set_title('Encircled energy fraction')
        ax1.legend()
        ax1.set_ylim(0, 1.05)

        # Derivative — flux added per annulus (normalised)
        dflux = np.diff(cumflux, prepend=0)
        ax2.plot(radii, dflux / np.nanmax(dflux), 'g-')
        ax2.axvline(r_conv, color='r', linestyle='--',
                    label=f'99% EE  r={r_conv:.0f} px')
        ax2.set_xlabel('Aperture radius [pixels]')
        ax2.set_ylabel('Normalised flux per annulus')
        ax2.set_title('Radial profile (annular flux)')
        ax2.legend()

        plt.tight_layout()
        plt.savefig(outpath, dpi=120)
        plt.close()
        print(f'  Plot saved: {outpath}')
    except ImportError:
        print('  (matplotlib not available — skipping plot)')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot',  action='store_true',
                        help='Save growth curve PNG plots')
    parser.add_argument('--star',  default=None,
                        help='Process only directories containing this string')
    parser.add_argument('--rmax',  type=int, default=R_MAX,
                        help=f'Maximum radius to probe in pixels (default {R_MAX})')
    args = parser.parse_args()

    base = BASE_DIR.resolve()
    star_roots = find_star_roots(base)
    if args.star:
        star_roots = [(l, p) for l, p in star_roots if args.star in l]
    print(f'Processing {len(star_roots)} star(s).\n')

    summary = []   # (label, filter, r_conv)

    for label, root in star_roots:
        for fname, fcfg in FILTER_CONFIG.items():
            fdir = find_filter_dir(root, fcfg['substr'])
            if fdir is None:
                continue
            fits_files = sorted(glob.glob(str(fdir / '*.fits')))
            if not fits_files:
                continue

            print(f'--- {label}  [{fname}]  ({len(fits_files)} frames) ---')
            print(f'  Registering …', end=' ', flush=True)
            radii, cumflux, cx, cy, sky = growth_curve_from_stack(
                fits_files, n_max=N_STACK, r_max=args.rmax)
            print(f'done.  centroid=({cx},{cy})  sky={sky:.1f} ADU/px')

            r_conv = convergence_radius(radii, cumflux)
            peak   = np.nanmax(cumflux)
            print(f'  99% EE radius: {r_conv:.0f} px  |  peak cumulative flux: {peak:,.0f} ADU')
            print_growth_table(radii, cumflux, r_conv)

            if args.plot:
                outpath = Path(f'growth_{label}_{fname}.png')
                save_plot(radii, cumflux, r_conv,
                          f'{label}  {fname}', outpath)
            print()
            summary.append((label, fname, r_conv, peak))

    # Summary table
    if len(summary) > 1:
        print('=' * 60)
        print(f'{"Star":<25} {"Filter":<12} {"r_99%EE":>8}  {"Peak ADU":>14}')
        print('-' * 60)
        for label, fname, r_conv, peak in summary:
            print(f'{label:<25} {fname:<12} {r_conv:>8.0f}  {peak:>14,.0f}')

        r99_vals = [r for _, _, r, _ in summary if np.isfinite(r)]
        if r99_vals:
            print(f'\n  Median 99% EE radius across all: {np.median(r99_vals):.0f} px')
            print(f'  Recommended R_APERTURE: {int(np.ceil(np.median(r99_vals) / 10) * 10)} px')


if __name__ == '__main__':
    main()
