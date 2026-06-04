"""
Estimate optical-train transmissivity from a Gaia reference star.

Method
------
  T_tel = (observed e-/s) / (predicted e-/s)

where

  predicted e-/s = ∫ F_λ(λ) · A_tel · T_atm(λ) · QE(λ) · T_filter(λ) dλ

F_λ(λ) is the photon flux density from flux_estimate.photon_flux_density().
T_tel is the wavelength-averaged throughput of the telescope + optics (mirrors,
spider, any relay optics), excluding atmosphere, filter, and detector QE.

Usage
-----
  python telescope_transmissivity.py

Adjust the CONFIGURATION section below as needed.
"""

import numpy as np
from astropy.io import fits
from pathlib import Path
import glob
import warnings
warnings.filterwarnings('ignore')

from flux_estimate import photon_flux_density

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

# Gaia DR3 star parameters
G_MAG   = 9.071853
TEFF    = 6114.846
BP_RP   = 0.6772032

# Telescope: Hooker 100-inch
D_PRIMARY   = 2.54    # m
D_SECONDARY = 0.36    # m  (estimated secondary obstruction; adjust as needed)
A_TEL = np.pi / 4.0 * (D_PRIMARY**2 - D_SECONDARY**2)  # m²

# Optical train: mirror/lens budget for T_tel sanity check
# 6 mirrors total: primary (Al, recently recoated) + 4 more Al + 1 Ag fold
# 1 relay refractor lens
# Set R_AL_PRIMARY = R_AL_AGED until confirmed recoated.
R_AL_PRIMARY = 0.91   # recently recoated Al (optimistic: 0.91, conservative: 0.87)
R_AL_AGED_LO = 0.85   # aged Al mirror, conservative
R_AL_AGED_HI = 0.90   # aged Al mirror, optimistic
R_AG         = 0.98   # Ag fold mirror before camera (manufacturer spec)
T_LENS_LO    = 0.92   # relay lens, uncoated or single-layer AR
T_LENS_HI    = 0.97   # relay lens, good broadband AR coating
N_AL_AGED    = 4      # secondary + tertiary + 2 relay tube mirrors

T_TEL_LO = R_AL_PRIMARY * R_AL_AGED_LO**N_AL_AGED * R_AG * T_LENS_LO
T_TEL_HI = R_AL_PRIMARY * R_AL_AGED_HI**N_AL_AGED * R_AG * T_LENS_HI

# Camera
# EGAIN: header reports 0.5166 e-/ADU, which matches the ZWO chart at gain~195 (0.1dB units).
# The data were taken at GAIN=252 (0.1dB units), which is in the HCG regime.
# Reading from the ZWO ASI585MM spec chart (asi585spec.webp) at gain=252:
#   EGAIN(chart) ≈ 0.30 e-/ADU
#   full-well     ≈ 1050 e-  →  EGAIN = 1050/4095 ≈ 0.256 e-/ADU
# Best estimate: 0.28 e-/ADU (midpoint).  Header value appears to be incorrect.
EGAIN = 0.057  # e-/ADU  (derived from T_tel self-consistency: bias RN=28.4 ADU
               #           × EGAIN → ~1.6 e- RN, consistent with HCG at gain=252;
               #           ZWO spec chart is unreadable in this regime)

# Plate scale
PLATE_SCALE_MAS = 18.63     # mas/pixel

# Aperture photometry radii (pixels)
R_APERTURE  = 100   # signal aperture; chosen to enclose >95% of seeing disk
R_SKY_INNER = 120   # sky background annulus inner radius
R_SKY_OUTER = 150   # sky background annulus outer radius

# Number of frames to stack (None = all)
N_STACK = None

# Data directories
DATA_ROOT = Path('G2P1060R_/_c1_ref_G2P1060R')
FILTERS = {
    'luminance': {
        'dir':      DATA_ROOT / '1_Luminance_0.04s',
        'exptime':  0.04,       # s
    },
    'ir_longpass': {
        'dir':      DATA_ROOT / '1_IR Long-Pass_0.06s',
        'exptime':  0.06,       # s
    },
}

# ---------------------------------------------------------------------------
# SPECTRAL CURVES
# ---------------------------------------------------------------------------

# ZWO ASI585MM Pro QE — digitised from ZWO published chart (zwoimx585qe.webp).
# Chart runs 400–1000 nm; set to zero outside that range.
_QE_WL  = np.array([400, 450, 500, 550, 575, 600, 625, 650, 675,
                     700, 725, 750, 775, 800, 825, 850, 875, 900,
                     925, 950, 975, 1000])  # nm
_QE_VAL = np.array([0.92, 0.91, 0.91, 0.91, 0.90, 0.89, 0.84, 0.80, 0.76,
                     0.74, 0.66, 0.65, 0.57, 0.50, 0.48, 0.46, 0.39, 0.30,
                     0.28, 0.29, 0.22, 0.17])  # fraction

def qe_asi585(wavelength_nm):
    """ZWO ASI585MM Pro QE as a function of wavelength [fraction]."""
    return np.interp(wavelength_nm, _QE_WL, _QE_VAL, left=0.0, right=0.0)


# Filter transmission curves — digitised from Baader manufacturer charts.
#
# Luminance: Baader UV/IR-Cut L-Filter CMOS-optimised (baader-uv-ir-cut*.webp)
#   Cut-on ~390 nm, flat ~95–97% from 400–680 nm, cut-off complete by ~710 nm.
#
# IR long-pass: Baader 685 IR-pass (baader685irpass.webp)
#   50% point ~683 nm, plateau ~96–97% from 700 nm onward.
#
_LUM_WL  = np.array([370, 385, 395, 400, 410, 680, 690, 700, 710])  # nm
_LUM_VAL = np.array([0.00, 0.10, 0.50, 0.90, 0.96, 0.96, 0.50, 0.10, 0.00])

_IRP_WL  = np.array([640, 650, 660, 670, 680, 683, 690, 700, 750, 1050])  # nm
_IRP_VAL = np.array([0.00, 0.02, 0.10, 0.40, 0.75, 0.50, 0.90, 0.96, 0.97, 0.96])

def filter_transmission(wavelength_nm, name):
    """Filter transmission [fraction] from digitised Baader manufacturer curves."""
    lam = np.asarray(wavelength_nm, dtype=float)
    if name == 'luminance':
        return np.interp(lam, _LUM_WL, _LUM_VAL, left=0.0, right=0.0)
    elif name == 'ir_longpass':
        return np.interp(lam, _IRP_WL, _IRP_VAL, left=0.0, right=_IRP_VAL[-1])
    else:
        raise ValueError(f'Unknown filter: {name}')


def atmospheric_transmission(wavelength_nm, airmass=1.003, altitude_m=1740):
    """
    Simple broadband atmospheric extinction model.

    Combines Rayleigh scattering (scaled for site altitude) and a fixed
    aerosol + ozone term.  Good to ~10% for routine optical work.
    """
    lam = np.asarray(wavelength_nm, dtype=float)
    lam_um = lam / 1000.0       # nm → µm

    # Rayleigh scattering: k_R ∝ λ^-4, normalised to sea-level
    # Pressure scale: MWO at 1740 m ≈ 830 hPa / 1013 hPa
    p_ratio  = np.exp(-altitude_m / 8500.0)
    k_rayleigh = 0.0094 * (lam_um / 0.4) ** (-4) * p_ratio

    # Aerosol (grey, weakly wavelength-dependent) + ozone (peaks ~600 nm)
    k_aerosol = 0.04 * (lam_um / 0.5) ** (-1.3)
    k_ozone   = 0.016 * np.exp(-0.5 * ((lam - 600) / 55) ** 2)   # Chappuis band

    k_total = k_rayleigh + k_aerosol + k_ozone
    return 10.0 ** (-0.4 * k_total * airmass)

# ---------------------------------------------------------------------------
# PHOTOMETRY
# ---------------------------------------------------------------------------

def smooth(image, sigma=15):
    """Gaussian smooth to wash out speckles and reveal the PSF envelope."""
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(image.astype(np.float32), sigma=sigma)


def find_star_in_frame(image, smooth_sigma=15):
    """
    Return (cx, cy) centroid of the star in a speckle frame.

    Smooths out the speckle pattern with a kernel ~= speckle size, then
    returns the peak of the smoothed image.
    """
    sm = smooth(image, smooth_sigma)
    cy, cx = np.unravel_index(np.argmax(sm), sm.shape)
    return int(cx), int(cy)


def aperture_photometry_single(image, cx, cy, r_ap, r_sky_inner, r_sky_outer):
    """
    Circular aperture photometry with annular sky subtraction.

    Returns net background-subtracted ADU in the aperture.
    """
    ny, nx = image.shape
    yy, xx = np.ogrid[:ny, :nx]
    r2 = (xx - cx)**2 + (yy - cy)**2

    ap_mask  = r2 <= r_ap**2
    sky_mask = (r2 > r_sky_inner**2) & (r2 <= r_sky_outer**2)

    # Guard: sky annulus must not fall outside the image
    if sky_mask.sum() < 10:
        return np.nan

    sky_per_pix = np.median(image[sky_mask])
    n_ap = ap_mask.sum()
    return float(image[ap_mask].sum() - sky_per_pix * n_ap)


def shift_and_add(fits_files, n_max=None, smooth_sigma=15):
    """
    Register frames by integer-pixel shift-and-add, then return the mean stack.

    For each frame the PSF centroid is located via Gaussian smoothing.
    Each frame is then shifted (integer-pixel roll) so that the centroid
    lands at a common reference position (the centroid of the first frame).
    The registered frames are mean-combined.

    Returns
    -------
    stack : 2-D ndarray
        Mean-combined registered image [ADU].
    cx_ref, cy_ref : int
        Star centroid position in the output stack (centre of the image
        after padding is removed).
    n_used : int
        Number of frames that contributed to the stack.
    """
    from scipy.ndimage import shift as ndimage_shift

    files = fits_files[:n_max] if n_max else fits_files
    n = len(files)
    print(f'  Registering {n} frames …', end=' ', flush=True)

    # First pass: collect all centroids
    centroids = []
    shapes = []
    for f in files:
        img = fits.getdata(f).astype(np.float32)
        shapes.append(img.shape)
        cx, cy = find_star_in_frame(img, smooth_sigma)
        centroids.append((cx, cy))

    ny, nx = shapes[0]
    cx_ref, cy_ref = centroids[0]   # register everything to the first frame

    # Second pass: shift and accumulate
    accum = np.zeros((ny, nx), dtype=np.float64)
    n_used = 0
    for f, (cx, cy) in zip(files, centroids):
        img = fits.getdata(f).astype(np.float32)
        dy = cy_ref - cy
        dx = cx_ref - cx
        # ndimage_shift uses (row, col) = (y, x) convention
        shifted = ndimage_shift(img, shift=(dy, dx), order=1,
                                mode='constant', cval=0.0)
        accum += shifted
        n_used += 1

    stack = accum / n_used
    print(f'done.  ({n_used}/{n} frames stacked, ref centroid ({cx_ref},{cy_ref}))')
    return stack, cx_ref, cy_ref, n_used


def photometry_from_stack(fits_files, n_max=None,
                          r_ap=R_APERTURE,
                          r_sky_inner=R_SKY_INNER,
                          r_sky_outer=R_SKY_OUTER,
                          smooth_sigma=15):
    """
    Shift-and-add stack, then aperture-photometer the registered mean image.

    Returns mean net ADU per *single* frame (stack sum / n_frames).
    """
    stack, cx, cy, n_used = shift_and_add(fits_files, n_max, smooth_sigma)
    # stack is already the per-frame mean
    net_adu = aperture_photometry_single(stack, cx, cy,
                                         r_ap, r_sky_inner, r_sky_outer)
    return net_adu


# ---------------------------------------------------------------------------
# PAD-AWARE REGISTRATION AND PHOTOMETRY
# ---------------------------------------------------------------------------

def register_frames(fits_files, n_max=None, smooth_sigma=15):
    """
    Register frames and record per-frame shifts.

    ndimage_shift with mode='constant' pads borders with zero when a frame
    is shifted.  Recording the actual shift vectors lets callers compute the
    largest circle that never overlaps a zero-padded region.

    Returns
    -------
    stack : 2-D ndarray   mean-combined registered image [ADU]
    cx_ref, cy_ref : int  star centroid in the stack
    dx, dy : ndarray      per-frame shifts (dx = cx_ref - cx_i, same for dy)
    ny, nx : int          frame dimensions
    """
    from scipy.ndimage import shift as ndimage_shift

    files = fits_files[:n_max] if n_max else fits_files
    centroids = []
    for f in files:
        img = fits.getdata(f).astype(np.float32)
        cx, cy = find_star_in_frame(img, smooth_sigma)
        centroids.append((cx, cy))

    ny, nx   = fits.getdata(files[0]).shape
    cx_ref, cy_ref = centroids[0]

    dx = np.array([cx_ref - cx for cx, cy in centroids])
    dy = np.array([cy_ref - cy for cx, cy in centroids])

    accum = np.zeros((ny, nx), dtype=np.float64)
    for f, (cx, cy) in zip(files, centroids):
        img = fits.getdata(f).astype(np.float32)
        shifted = ndimage_shift(img, shift=(cy_ref - cy, cx_ref - cx),
                                order=1, mode='constant', cval=0.0)
        accum += shifted

    return accum / len(files), cx_ref, cy_ref, dx, dy, ny, nx


def pad_safe_radius(cx, cy, dx, dy, ny, nx):
    """
    Largest radius centred on (cx, cy) guaranteed free of zero-padded pixels.

    When a frame is shifted by (dx_i, dy_i):
      dx_i > 0  →  left  border padded (width dx_i)
      dx_i < 0  →  right border padded (width |dx_i|)
      dy_i > 0  →  top   border padded (width dy_i)
      dy_i < 0  →  bottom border padded (width |dy_i|)
    Worst-case pad on each side is the maximum over all frames.
    """
    pad_l = int(np.ceil(max(dx.max(), 0)))
    pad_r = int(np.ceil(max(-dx.min(), 0)))
    pad_t = int(np.ceil(max(dy.max(), 0)))
    pad_b = int(np.ceil(max(-dy.min(), 0)))
    return int(min(cx - pad_l, (nx - cx) - pad_r,
                   cy - pad_t, (ny - cy) - pad_b))


def fit_fwhm_moffat(stack, cx, cy, ny, nx, beta=3.0, r_fit_max=200):
    """
    Fit a Moffat profile (fixed beta) to the sky-subtracted azimuthal mean.

    Returns FWHM in pixels, or None if the fit fails.
    """
    from scipy.optimize import curve_fit

    yy, xx  = np.ogrid[:ny, :nx]
    r_grid  = np.sqrt((xx - cx)**2 + (yy - cy)**2)

    # Sky from annulus just beyond the fit region (avoid zero-pad border)
    sky_lo  = r_fit_max + 10
    sky_hi  = min(r_fit_max + 50, int(min(cx, nx - cx, cy, ny - cy)) - 5)
    sky_mask = (r_grid > sky_lo) & (r_grid <= sky_hi)
    sky = np.median(stack[sky_mask]) if sky_mask.sum() >= 20 else 0.0

    r_int   = r_grid.ravel().astype(int)
    vals    = (stack - sky).ravel()
    bins    = np.arange(0, r_fit_max + 1)
    profile = np.array([vals[r_int == r].mean() if (r_int == r).any() else 0.0
                        for r in bins])
    mask = (bins > 0) & (bins <= r_fit_max) & (profile > 0)

    def moffat(r, I0, alpha):
        return I0 * (1.0 + (r / alpha)**2)**(-beta)

    try:
        popt, _ = curve_fit(moffat, bins[mask], profile[mask],
                            p0=[profile[1:5].max(), 30.0],
                            bounds=([0, 1], [np.inf, 500]),
                            maxfev=5000)
        alpha = popt[1]
        return 2.0 * alpha * np.sqrt(2.0**(1.0 / beta) - 1.0)
    except RuntimeError:
        return None


def pad_aware_photometry(fits_files, n_max=None, smooth_sigma=15,
                         sky_margin=30, sky_width=25,
                         ap_gap=10, min_fwhm_scale=3.0):
    """
    Shift-and-add photometry with a sky annulus placed safely clear of
    zero-padded borders.

    The signal aperture extends to sky_inner - ap_gap, capturing the full
    SAA halo rather than just the diffraction core.

    Parameters
    ----------
    sky_margin : int
        Gap [px] between sky annulus inner edge and the pad-safe radius.
    sky_width : int
        Width [px] of the sky annulus.
    ap_gap : int
        Gap [px] between signal aperture and sky annulus inner edge.
    min_fwhm_scale : float
        Flag the result as unreliable if sky_inner < min_fwhm_scale × FWHM.

    Returns
    -------
    net_adu   : float   background-subtracted ADU in signal aperture
    sky_inner : int     inner radius of sky annulus used [px]
    sky_outer : int     outer radius of sky annulus used [px]
    r_ap      : int     signal aperture radius used [px]
    fwhm      : float   Moffat FWHM [px], or None if fit failed
    reliable  : bool    True if sky_inner >= min_fwhm_scale × FWHM
    """
    print(f'  Registering {min(len(fits_files), n_max or len(fits_files))} frames …',
          end=' ', flush=True)
    stack, cx, cy, dx, dy, ny, nx = register_frames(fits_files, n_max, smooth_sigma)
    print(f'done.  centroid=({cx},{cy})')

    safe_r    = pad_safe_radius(cx, cy, dx, dy, ny, nx)
    sky_inner = safe_r - sky_margin
    sky_outer = min(sky_inner + sky_width, safe_r - 2)
    r_ap      = sky_inner - ap_gap

    fwhm     = fit_fwhm_moffat(stack, cx, cy, ny, nx)
    reliable = (fwhm is not None) and (sky_inner >= min_fwhm_scale * fwhm)

    net_adu = aperture_photometry_single(stack, cx, cy,
                                         r_ap, sky_inner, sky_outer)
    return net_adu, sky_inner, sky_outer, r_ap, fwhm, reliable

# ---------------------------------------------------------------------------
# PREDICTED SIGNAL
# ---------------------------------------------------------------------------

def predicted_electrons_per_sec(filter_name, airmass=1.003, altitude_m=1740,
                                 wav_min=330, wav_max=1050, n_wav=1000):
    """
    Predicted photon → electron rate [e-/s] at the detector.

      predicted = ∫ F_λ(λ) · A_tel · T_atm(λ) · QE(λ) · T_filter(λ) dλ

    where F_λ is in ph/s/m²/nm and the integral is over nm.
    """
    lam = np.linspace(wav_min, wav_max, n_wav)   # nm
    dlam = lam[1] - lam[0]

    F    = photon_flux_density(lam, G_MAG, TEFF, BP_RP)  # ph/s/m²/nm
    T_atm  = atmospheric_transmission(lam, airmass, altitude_m)
    qe     = qe_asi585(lam)
    T_filt = filter_transmission(lam, filter_name)

    integrand = F * A_TEL * T_atm * qe * T_filt  # e-/s/nm
    return float(np.trapz(integrand, dx=dlam))

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def analyse_filter(name, cfg):
    print(f'\n=== {name.upper()} (t_exp = {cfg["exptime"]} s) ===')

    fits_files = sorted(glob.glob(str(cfg['dir'] / '*.fits')))
    if not fits_files:
        print('  No FITS files found — skipping.')
        return

    mean_net_adu = photometry_from_stack(fits_files, n_max=N_STACK)

    print(f'  Mean net ADU/frame (r={R_APERTURE} px): {mean_net_adu:,.0f}')

    observed_e_per_s = mean_net_adu * EGAIN / cfg['exptime']
    print(f'  Observed rate: {observed_e_per_s:,.0f} e-/s')

    pred_e_per_s = predicted_electrons_per_sec(name)
    t_tel = observed_e_per_s / pred_e_per_s
    print(f'  Predicted rate (T_tel=1): {pred_e_per_s:,.0f} e-/s')
    print(f'  Telescope transmissivity T_tel = {t_tel:.3f}')
    flag = ''
    if t_tel < T_TEL_LO * 0.7:
        flag = '  *** well below expected range — check QE/filter curves or EGAIN'
    elif t_tel > T_TEL_HI * 1.1:
        flag = '  *** above expected range — check QE/filter curves or EGAIN'
    print(f'  Expected T_tel range: {T_TEL_LO:.3f} – {T_TEL_HI:.3f}'
          f'  (1×Al fresh + {N_AL_AGED}×Al aged + 1×Ag + lens){flag}')
    print(f'  (A_tel = {A_TEL:.3f} m²)')

    return observed_e_per_s, pred_e_per_s


if __name__ == '__main__':
    print(f'Telescope area: {A_TEL:.3f} m²  '
          f'(D_primary={D_PRIMARY} m, D_secondary={D_SECONDARY} m)')
    print(f'Expected T_tel: {T_TEL_LO:.3f} – {T_TEL_HI:.3f}  '
          f'(1×Al fresh + {N_AL_AGED}×Al aged + 1×Ag + lens)')
    for name, cfg in FILTERS.items():
        analyse_filter(name, cfg)
