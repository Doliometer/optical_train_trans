"""
Estimate seeing-disk size by fitting Moffat profiles to shift-and-add stacks.

For 3 randomly chosen stars (luminance filter, 20 frames each):
  - Build a registered mean stack
  - Extract the azimuthally averaged radial profile
  - Fit a Moffat profile with β fixed at 2.5, 3.0, 3.5
  - Report FWHM and plot radial profile + fits

Usage
-----
  python3 moffat_fit.py [--seed N] [--nstack N] [--nstars N]
"""

import argparse
import random
import glob
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from astropy.io import fits
from pathlib import Path
from scipy.ndimage import shift as ndimage_shift
from scipy.optimize import curve_fit

# ---------------------------------------------------------------------------
# reuse discovery + registration from existing scripts
# ---------------------------------------------------------------------------
from survey_transmissivity import find_star_roots, find_filter_dir, FILTER_CONFIG
from telescope_transmissivity import find_star_in_frame

BASE_DIR = Path('.')

# ---------------------------------------------------------------------------
# Moffat model
# ---------------------------------------------------------------------------

def moffat_1d(r, I0, alpha, beta):
    """
    1-D Moffat profile (circularly symmetric).
      I(r) = I0 * (1 + (r/alpha)^2)^(-beta)
    FWHM = 2 * alpha * sqrt(2^(1/beta) - 1)
    """
    return I0 * (1.0 + (r / alpha)**2)**(-beta)


def moffat_fwhm(alpha, beta):
    return 2.0 * alpha * np.sqrt(2.0**(1.0 / beta) - 1.0)


# ---------------------------------------------------------------------------
# stack + radial profile
# ---------------------------------------------------------------------------

def make_stack(fits_files, n_max, smooth_sigma=15):
    files = fits_files[:n_max]
    centroids = []
    for f in files:
        img = fits.getdata(f).astype(np.float32)
        cx, cy = find_star_in_frame(img, smooth_sigma)
        centroids.append((cx, cy))

    ny, nx = fits.getdata(files[0]).shape
    cx_ref, cy_ref = centroids[0]
    accum = np.zeros((ny, nx), dtype=np.float64)
    for f, (cx, cy) in zip(files, centroids):
        img = fits.getdata(f).astype(np.float32)
        shifted = ndimage_shift(img, shift=(cy_ref - cy, cx_ref - cx),
                                order=1, mode='constant', cval=0.0)
        accum += shifted
    return accum / len(files), cx_ref, cy_ref


def radial_profile(stack, cx, cy, r_max, sky_inner=350, sky_outer=420):
    """
    Azimuthally averaged, sky-subtracted radial profile.

    Sky annulus deliberately pushed well out (sky_inner/sky_outer)
    to avoid sitting in the SAA halo.
    Returns (r_px, mean_adu_per_px).
    """
    ny, nx = stack.shape
    yy, xx = np.ogrid[:ny, :nx]
    r_grid = np.sqrt((xx - cx)**2 + (yy - cy)**2)

    sky_mask = (r_grid > sky_inner) & (r_grid <= sky_outer)
    if sky_mask.sum() < 20:
        sky = 0.0
    else:
        sky = np.median(stack[sky_mask])

    stack_sub = stack - sky

    # bin into 1-px-wide annuli
    r_int = r_grid.astype(int)
    r_bins = np.arange(0, r_max + 1)
    profile = np.zeros(len(r_bins))
    counts  = np.zeros(len(r_bins))
    for r in r_bins:
        mask = r_int == r
        if mask.sum():
            profile[r] = stack_sub[mask].mean()
            counts[r]  = mask.sum()

    return r_bins.astype(float), profile


# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------

BETAS = [2.5, 3.0, 3.5]
COLORS = ['tab:blue', 'tab:orange', 'tab:green']


def fit_moffat_fixed_beta(r, profile, beta, r_fit_max=200):
    """Fit Moffat with beta fixed; free params are I0 and alpha."""
    mask = (r > 0) & (r <= r_fit_max) & np.isfinite(profile) & (profile > 0)
    r_fit = r[mask]
    p_fit = profile[mask]

    def model(r_, I0, alpha):
        return moffat_1d(r_, I0, alpha, beta)

    try:
        p0    = [p_fit.max(), 30.0]
        popt, _ = curve_fit(model, r_fit, p_fit, p0=p0,
                            bounds=([0, 1], [np.inf, 500]),
                            maxfev=5000)
        return popt  # (I0, alpha)
    except RuntimeError:
        return None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed',   type=int, default=42,
                        help='Random seed for star selection (default 42)')
    parser.add_argument('--nstack', type=int, default=20,
                        help='Frames to stack per star (default 20)')
    parser.add_argument('--nstars', type=int, default=3,
                        help='Number of random stars to sample (default 3)')
    parser.add_argument('--rmax',   type=int, default=450,
                        help='Max radius for radial profile (default 450)')
    parser.add_argument('--rfit',   type=int, default=200,
                        help='Max radius used in Moffat fit (default 200)')
    parser.add_argument('--sigma',  type=int, default=15,
                        help='Gaussian smooth sigma for centroiding (default 15)')
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        have_mpl = True
    except ImportError:
        have_mpl = False

    base = BASE_DIR.resolve()
    all_roots = find_star_roots(base)
    random.seed(args.seed)
    chosen = random.sample(all_roots, min(args.nstars, len(all_roots)))
    print(f'Selected stars: {[l for l,_ in chosen]}\n')

    plate_scale_mas = 18.63   # mas/px

    if have_mpl:
        fig, axes = plt.subplots(1, len(chosen),
                                 figsize=(6 * len(chosen), 5),
                                 sharey=False)
        if len(chosen) == 1:
            axes = [axes]

    summary_rows = []

    for idx, (label, root) in enumerate(chosen):
        fdir = find_filter_dir(root, 'Luminance')
        if fdir is None:
            print(f'{label}: no Luminance directory found — skipping.')
            continue
        fits_files = sorted(glob.glob(str(fdir / '*.fits')))
        if not fits_files:
            print(f'{label}: no FITS files — skipping.')
            continue

        n_use = min(args.nstack, len(fits_files))
        print(f'{label}: stacking {n_use} frames …', end=' ', flush=True)
        stack, cx, cy = make_stack(fits_files, n_use, smooth_sigma=args.sigma)
        print(f'done.  centroid=({cx},{cy})')

        r, profile = radial_profile(stack, cx, cy, r_max=args.rmax)

        peak = profile[1:20].max() if profile[1:20].max() > 0 else 1.0

        print(f'  {"beta":>5}  {"alpha(px)":>10}  {"FWHM(px)":>10}  {"FWHM(\")":>8}  {"fit_ok":>6}')
        print(f'  {"-----":>5}  {"----------":>10}  {"----------":>10}  {"--------":>8}  {"------":>6}')

        fits_results = {}
        for beta, col in zip(BETAS, COLORS):
            popt = fit_moffat_fixed_beta(r, profile, beta, r_fit_max=args.rfit)
            if popt is not None:
                I0, alpha = popt
                fwhm_px  = moffat_fwhm(alpha, beta)
                fwhm_as  = fwhm_px * plate_scale_mas / 1000.0
                print(f'  {beta:>5.1f}  {alpha:>10.1f}  {fwhm_px:>10.1f}  {fwhm_as:>8.3f}  {"yes":>6}')
                fits_results[beta] = (I0, alpha, fwhm_px)
                summary_rows.append((label, beta, fwhm_px, fwhm_as))
            else:
                print(f'  {beta:>5.1f}  {"—":>10}  {"—":>10}  {"—":>8}  {"FAIL":>6}')
        print()

        if have_mpl:
            ax = axes[idx]
            r_plot = r[r <= args.rmax]
            p_plot = profile[r <= args.rmax] / peak
            ax.plot(r_plot, p_plot, 'k-', lw=1, alpha=0.6, label='data')
            r_fine = np.linspace(0, args.rfit, 500)
            for beta, col in zip(BETAS, COLORS):
                if beta in fits_results:
                    I0, alpha, fwhm_px = fits_results[beta]
                    y_fine = moffat_1d(r_fine, I0, alpha, beta) / peak
                    ax.plot(r_fine, y_fine, color=col, lw=1.5,
                            label=f'β={beta}  FWHM={fwhm_px:.0f}px')
            ax.set_xlabel('Radius [px]')
            ax.set_ylabel('Normalised intensity')
            ax.set_title(label)
            ax.set_xlim(0, args.rmax)
            ax.set_ylim(-0.1, 1.05)
            ax.axhline(0, color='gray', lw=0.5, ls=':')
            ax.legend(fontsize=8)

    # summary
    if summary_rows:
        print('=' * 60)
        print(f'{"Star":<25} {"β":>5}  {"FWHM(px)":>10}  {"FWHM(\")":>8}')
        print('-' * 60)
        for label, beta, fwhm_px, fwhm_as in summary_rows:
            print(f'{label:<25} {beta:>5.1f}  {fwhm_px:>10.1f}  {fwhm_as:>8.3f}')

    if have_mpl:
        plt.suptitle('Moffat fits to SAA stack radial profiles (luminance)', y=1.01)
        plt.tight_layout()
        outpath = Path('moffat_fits.png')
        plt.savefig(outpath, dpi=120, bbox_inches='tight')
        print(f'\nPlot saved: {outpath}')


if __name__ == '__main__':
    main()
