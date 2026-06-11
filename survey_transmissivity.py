"""
Survey optical-train transmissivity across all reference stars observed on the same night.

For each star directory the script:
  1. Reads the Gaia DR3 source ID from the FITS header (GAIANUM~ keyword).
  2. Queries Gaia DR3 for G magnitude, T_eff, BP-RP, and XP spectrum availability.
  3. Fetches Gaia XP sampled spectra (336-1020 nm) via datalink for all available stars.
     Falls back to a Planck blackbody SED for the two stars without XP spectra.
  4. Performs shift-and-add stacking + aperture photometry for each filter.
  5. Computes T_tel = observed e-/s / predicted e-/s.
  6. Prints a summary table with per-filter statistics.

Usage
-----
  python3 survey_transmissivity.py
"""

import numpy as np
from astropy.io import fits
from pathlib import Path
import glob
import re
import pickle
import warnings
warnings.filterwarnings('ignore')

from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
Gaia.ROW_LIMIT = -1

from flux_estimate import photon_flux_density, xp_photon_flux_density
from telescope_transmissivity import (
    A_TEL, EGAIN, T_TEL_LO, T_TEL_HI, T_TEL_IR_LO, T_TEL_IR_HI, N_AL_AGED,
    qe_asi585, filter_transmission, atmospheric_transmission,
    pad_aware_photometry,
)

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

BASE_DIR = Path('.')        # directory containing all G2P*_ folders
N_STACK  = None             # frames per filter (None = all)
GAIA_CACHE    = Path('gaia_cache.pkl')   # cached Gaia photometry query
XP_CACHE      = Path('xp_cache.pkl')    # cached XP sampled spectra

# Stars excluded from T_tel summary statistics (but still processed and printed):
#   G2P1190R: very red star (bp_rp=1.52, Teff=4776K); XP flux suppressed in
#             luminance band by steep SED slope; unreliable T_tel in both bands.
#   G2P1170R: no XP spectrum; Planck SED at Teff=7594K underpredicts IR flux,
#             giving T_tel > 1 in IR_longpass.
EXCLUDE_FROM_SUMMARY = {'_c1_ref_G2P1190R', '_c1_ref_G2P1170R'}

# Filter name → exposure time; matched by substring in subdir name
FILTER_CONFIG = {
    'luminance':   {'substr': 'Luminance',    'exptime': 0.04},
    'ir_longpass': {'substr': 'IR Long-Pass', 'exptime': 0.06},
}

# ---------------------------------------------------------------------------
# DISCOVERY
# ---------------------------------------------------------------------------

def find_star_roots(base):
    """
    Return list of (star_label, data_root_path) for all star directories.

    Handles two layouts:
      base/_c1_ref_G2P????          (the reference star observed directly)
      base/G2P????_/_c1_ref_G2P???? (the 18 additional stars)
    """
    roots = []
    for p in sorted(base.glob('_c1_ref_*')):
        if p.is_dir():
            roots.append((p.name, p))
    for outer in sorted(base.glob('G2P*_')):
        for inner in sorted(outer.glob('_c1_ref_*')):
            if inner.is_dir():
                roots.append((inner.name, inner))

    seen = set()
    unique = []
    for label, path in roots:
        if label not in seen:
            seen.add(label)
            unique.append((label, path))
    return unique


def find_filter_dir(star_root, substr):
    """Return the subdirectory whose name contains `substr`, or None."""
    for p in star_root.iterdir():
        if p.is_dir() and substr in p.name:
            return p
    return None


def gaia_source_id_from_fits(star_root):
    """Read Gaia DR3 source ID from the first FITS file found."""
    for pattern in ['*_Luminance_*/*.fits', '*_IR*/*.fits', '*/*.fits']:
        files = sorted(glob.glob(str(star_root / pattern)))
        if files:
            h = fits.getheader(files[0])
            raw = h.get('GAIANUM~', '')
            m = re.search(r'DR3\s+(\d+)', raw)
            if m:
                return int(m.group(1))
    return None


def airmass_from_fits(star_root):
    """Read AIRMASS from the first FITS file found."""
    for pattern in ['*_Luminance_*/*.fits', '*_IR*/*.fits', '*/*.fits']:
        files = sorted(glob.glob(str(star_root / pattern)))
        if files:
            h = fits.getheader(files[0])
            return float(h.get('AIRMASS', 1.0))
    return 1.0

# ---------------------------------------------------------------------------
# GAIA QUERIES
# ---------------------------------------------------------------------------

def query_gaia_stars(source_ids):
    """
    Query Gaia DR3 for G magnitude, T_eff, BP-RP, and XP availability.

    Returns dict: source_id → {g_mag, teff, bp_rp, has_xp}
    """
    id_list = ','.join(str(i) for i in source_ids)
    adql = f"""
        SELECT source_id, phot_g_mean_mag, teff_gspphot, bp_rp, has_xp_sampled
        FROM gaiadr3.gaia_source
        WHERE source_id IN ({id_list})
    """
    print('Querying Gaia DR3 …', end=' ', flush=True)
    job = Gaia.launch_job(adql)
    tbl = job.get_results()
    print(f'done ({len(tbl)} rows).')

    result = {}
    sid_col = 'SOURCE_ID' if 'SOURCE_ID' in tbl.colnames else 'source_id'
    for row in tbl:
        sid  = int(row[sid_col])
        teff = float(row['teff_gspphot'])
        if np.ma.is_masked(teff):
            teff = None
        bp_rp = float(row['bp_rp'])
        if np.ma.is_masked(bp_rp):
            bp_rp = None
        result[sid] = {
            'g_mag':  float(row['phot_g_mean_mag']),
            'teff':   teff,
            'bp_rp':  bp_rp,
            'has_xp': bool(row['has_xp_sampled']),
        }
    return result


def fetch_xp_spectra(source_ids):
    """
    Fetch Gaia XP sampled spectra via datalink for a list of source IDs.

    Returns dict: source_id → (wavelength_nm ndarray, flux_W ndarray)
    """
    print(f'Fetching XP spectra for {len(source_ids)} stars …', end=' ', flush=True)
    dl = Gaia.load_data(ids=source_ids,
                        data_release='Gaia DR3',
                        retrieval_type='XP_SAMPLED',
                        format='csv',
                        verbose=False)
    spectra = {}
    for key, val in dl.items():
        m = re.search(r'(\d{10,})', key)
        if m:
            sid = int(m.group(1))
            t   = val[0]
            spectra[sid] = (np.array(t['wavelength'], dtype=float),
                            np.array(t['flux'],       dtype=float))
    print(f'done ({len(spectra)} spectra).')
    return spectra

# ---------------------------------------------------------------------------
# PREDICTED SIGNAL
# ---------------------------------------------------------------------------

def predicted_electrons_per_sec(gp, xp_spectra, filter_name,
                                  airmass=1.0, altitude_m=1740,
                                  wav_min=336, wav_max=1020, n_wav=1000):
    """
    Predicted electron rate [e-/s].

    Uses the Gaia XP sampled spectrum (W/m²/nm → ph/s/m²/nm) if available;
    otherwise falls back to a Planck blackbody normalised to the G magnitude.
    """
    lam  = np.linspace(wav_min, wav_max, n_wav)
    dlam = lam[1] - lam[0]

    sid = gp.get('source_id')
    if sid in xp_spectra:
        wl_xp, fl_xp = xp_spectra[sid]
        F = xp_photon_flux_density(lam, wl_xp, fl_xp)
        sed = 'XP'
    else:
        if gp['teff'] is None or gp['bp_rp'] is None:
            return None, 'no-SED'
        F   = photon_flux_density(lam, gp['g_mag'], gp['teff'], gp['bp_rp'])
        sed = 'Planck'

    T_atm  = atmospheric_transmission(lam, airmass, altitude_m)
    qe     = qe_asi585(lam)
    T_filt = filter_transmission(lam, filter_name)
    pred   = float(np.trapz(F * A_TEL * T_atm * qe * T_filt, dx=dlam))
    return pred, sed

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    base = BASE_DIR.resolve()
    star_roots = find_star_roots(base)
    print(f'Found {len(star_roots)} star directories.')

    # Collect source IDs and airmasses
    star_info  = []
    source_ids = []
    for label, root in star_roots:
        sid = gaia_source_id_from_fits(root)
        am  = airmass_from_fits(root)
        star_info.append({'label': label, 'root': root,
                          'source_id': sid, 'airmass': am})
        if sid and sid not in source_ids:
            source_ids.append(sid)

    # Gaia photometry + XP availability (use cache if present)
    if GAIA_CACHE.exists():
        print(f'Loading Gaia photometry from cache ({GAIA_CACHE}) …')
        with open(GAIA_CACHE, 'rb') as f:
            gaia_params = pickle.load(f)
    else:
        gaia_params = query_gaia_stars(source_ids)
        with open(GAIA_CACHE, 'wb') as f:
            pickle.dump(gaia_params, f)
        print(f'Gaia photometry cached to {GAIA_CACHE}.')
    # Inject source_id into each param dict for lookup in predicted_electrons_per_sec
    for sid, gp in gaia_params.items():
        gp['source_id'] = sid

    # Batch fetch XP spectra (use cache if present)
    xp_ids = [sid for sid, gp in gaia_params.items() if gp['has_xp']]
    if GAIA_CACHE.exists() and XP_CACHE.exists():
        print(f'Loading XP spectra from cache ({XP_CACHE}) …')
        with open(XP_CACHE, 'rb') as f:
            xp_spectra = pickle.load(f)
    else:
        xp_spectra = fetch_xp_spectra(xp_ids) if xp_ids else {}
        with open(XP_CACHE, 'wb') as f:
            pickle.dump(xp_spectra, f)
        print(f'XP spectra cached to {XP_CACHE}.')

    # -----------------------------------------------------------------------
    # Per-star, per-filter analysis
    # -----------------------------------------------------------------------
    print(f'\n{"Star":<25} {"Filter":<12} {"G":>5} {"Teff":>6} {"SED":<6} '
          f'{"Obs e/s":>12} {"Pred e/s":>12} {"T_tel":>7} {"FWHM\"":>6} {"ok?":>4}')
    print('-' * 100)

    results = []
    for si in star_info:
        label = si['label']
        root  = si['root']
        sid   = si['source_id']
        am    = si['airmass']

        if sid not in gaia_params:
            print(f'{label:<25}  Gaia source {sid} not found — skipping.')
            continue

        gp = gaia_params[sid]
        teff_str = f'{gp["teff"]:.0f}' if gp['teff'] else '  —  '

        for fname, fcfg in FILTER_CONFIG.items():
            fdir = find_filter_dir(root, fcfg['substr'])
            if fdir is None:
                continue
            fits_files = sorted(glob.glob(str(fdir / '*.fits')))
            if not fits_files:
                continue

            mean_adu, sky_inner, sky_outer, r_ap, fwhm, reliable = \
                pad_aware_photometry(fits_files, n_max=N_STACK)
            obs_eps  = mean_adu * EGAIN / fcfg['exptime']

            fwhm_as_str = f'{fwhm * 18.63e-3:.2f}' if fwhm else '  — '

            pred_eps, sed = predicted_electrons_per_sec(
                gp, xp_spectra, fname, airmass=am)

            if pred_eps is None:
                print(f'{label:<25} {fname:<12} {gp["g_mag"]:>5.2f} {teff_str:>6} '
                      f'{sed:<6} {obs_eps:>12,.0f} {"—":>12} {"—":>7} '
                      f'{fwhm_as_str:>6} {"—":>4}')
                continue

            t_tel    = obs_eps / pred_eps
            ok_str   = 'yes' if reliable else 'NO'
            print(f'{label:<25} {fname:<12} {gp["g_mag"]:>5.2f} {teff_str:>6} '
                  f'{sed:<6} {obs_eps:>12,.0f} {pred_eps:>12,.0f} {t_tel:>7.3f} '
                  f'{fwhm_as_str:>6} {ok_str:>4}')

            results.append({
                'star': label, 'filter': fname, 'sed': sed,
                'g_mag': gp['g_mag'], 'teff': gp['teff'],
                'source_id': sid, 'airmass': am,
                'obs_e_per_s': obs_eps, 'pred_e_per_s': pred_eps,
                't_tel': t_tel, 'fwhm': fwhm, 'reliable': reliable,
            })

    # -----------------------------------------------------------------------
    # Summary statistics
    # -----------------------------------------------------------------------
    if not results:
        print('No results to summarise.')
        return

    print('-' * 90)
    print(f'\n  Expected T_tel (luminance, 400–710 nm): {T_TEL_LO:.3f} – {T_TEL_HI:.3f}  '
          f'(1×Al fresh + {N_AL_AGED}×Al aged + 1×Ag + lens)')
    print(f'  Expected T_tel (IR longpass, 685–1050 nm): {T_TEL_IR_LO:.3f} – {T_TEL_IR_HI:.3f}  '
          f'(Al reflectivity from Rakić 1995)')

    for fname in FILTER_CONFIG:
        rows_all      = [r for r in results if r['filter'] == fname]
        rows_reliable = [r for r in rows_all
                         if r['reliable'] and r['star'] not in EXCLUDE_FROM_SUMMARY]
        rows_use      = rows_reliable if rows_reliable else rows_all
        tag = '' if rows_reliable else '  [no reliable stars — using all]'

        vals   = np.array([r['t_tel'] for r in rows_use], dtype=float)
        finite = vals[np.isfinite(vals)]
        n_xp   = sum(1 for r in rows_use if r['sed'] == 'XP' and np.isfinite(r['t_tel']))
        n_unreliable = sum(1 for r in rows_all if not r['reliable'])
        n_excluded   = sum(1 for r in rows_all
                           if r['reliable'] and r['star'] in EXCLUDE_FROM_SUMMARY)

        print(f'\n  {fname.upper()}  n={len(finite)} used  '
              f'({n_xp} XP, {len(finite)-n_xp} Planck,  '
              f'{n_unreliable} unreliable, {n_excluded} excluded){tag}')
        if len(finite):
            median_t = np.nanmedian(finite)
            mean_t   = np.nanmean(finite)
            print(f'    median T_tel = {median_t:.3f}')
            print(f'    mean   T_tel = {mean_t:.3f}  ±  {np.nanstd(finite):.3f}')
            print(f'    range        = {finite.min():.3f} – {finite.max():.3f}')
            # Implied EGAIN: the value that would make the median T_tel equal to
            # the midpoint of the expected range for this filter.
            if fname == 'luminance':
                t_mid = (T_TEL_LO + T_TEL_HI) / 2
            else:
                t_mid = (T_TEL_IR_LO + T_TEL_IR_HI) / 2
            egain_implied = EGAIN * median_t / t_mid
            print(f'    Implied EGAIN for T_tel={t_mid:.3f} (mid expected range): '
                  f'{egain_implied:.4f} e-/ADU  (used: {EGAIN:.4f})')
            if fname == 'luminance':
                print(f'    Note: EGAIN={EGAIN:.4f} was set so that luminance T_tel '
                      f'falls within the expected mirror-budget range '
                      f'({T_TEL_LO:.3f}–{T_TEL_HI:.3f}).')


if __name__ == '__main__':
    main()
