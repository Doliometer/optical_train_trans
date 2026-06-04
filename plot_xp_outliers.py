"""
Plot XP spectrum for G2P1190R overlaid with filter transmission curves and QE,
to diagnose why observed flux is far below predicted in both bands.

Usage
-----
  python3 plot_xp_outliers.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
Gaia.ROW_LIMIT = -1

from flux_estimate import xp_photon_flux_density, photon_flux_density
from telescope_transmissivity import (
    qe_asi585, filter_transmission, atmospheric_transmission,
    A_TEL, EGAIN,
)

# ---------------------------------------------------------------------------
# Star parameters (from survey output)
# ---------------------------------------------------------------------------
STARS = {
    'G2P1190R': {'source_id': 170677915834391040, 'g_mag': 9.16,  'teff': 4776,
                 'airmass': 1.0, 'has_xp': True},
    'G2P1170R': {'source_id': 235644312705236992, 'g_mag': 10.95, 'teff': 7594,
                 'airmass': 1.0, 'has_xp': False},
}

# ---------------------------------------------------------------------------
# Fetch stellar parameters and XP spectra by source ID
# ---------------------------------------------------------------------------

def fetch_params(source_ids):
    id_list = ','.join(str(i) for i in source_ids)
    adql = f"""
        SELECT source_id, phot_g_mean_mag, teff_gspphot, bp_rp, has_xp_sampled
        FROM gaiadr3.gaia_source
        WHERE source_id IN ({id_list})
    """
    print('Querying Gaia DR3 …', end=' ', flush=True)
    t = Gaia.launch_job(adql).get_results()
    print(f'{len(t)} rows')
    result = {}
    for row in t:
        sid = int(row['SOURCE_ID'])
        result[sid] = {
            'g_mag':  float(row['phot_g_mean_mag']),
            'teff':   float(row['teff_gspphot']),
            'bp_rp':  float(row['bp_rp']),
            'has_xp': bool(row['has_xp_sampled']),
        }
        print(f"  {sid}  G={result[sid]['g_mag']:.3f}  "
              f"Teff={result[sid]['teff']:.0f}  bp_rp={result[sid]['bp_rp']:.3f}  "
              f"has_xp={result[sid]['has_xp']}")
    return result


def fetch_xp(source_id):
    import re
    dl = Gaia.load_data(ids=[source_id],
                        data_release='Gaia DR3',
                        retrieval_type='XP_SAMPLED',
                        format='csv',
                        verbose=False)
    if not dl:
        return None, None
    key = list(dl.keys())[0]
    t = dl[key][0]
    return np.array(t['wavelength'], dtype=float), np.array(t['flux'], dtype=float)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_star(ax_sed, ax_resp, label, wl_xp, fl_xp,
              g_mag, teff, bp_rp, airmass, sed_type):
    lam = np.linspace(330, 1050, 1000)

    # SED
    if wl_xp is not None:
        F = xp_photon_flux_density(lam, wl_xp, fl_xp)
        sed_label = 'XP spectrum'
    else:
        F = photon_flux_density(lam, g_mag, teff, bp_rp if bp_rp else 0.5)
        sed_label = f'Planck Teff={teff}K'

    T_atm = atmospheric_transmission(lam, airmass)
    qe    = qe_asi585(lam)
    T_lum = filter_transmission(lam, 'luminance')
    T_ir  = filter_transmission(lam, 'ir_longpass')

    # normalise SED to peak=1 for display
    Fn = F / F.max()

    ax_sed.plot(lam, Fn, 'k-', lw=1.5, label=sed_label)
    if wl_xp is not None:
        # also show raw XP points
        ax_sed.plot(wl_xp, fl_xp / fl_xp.max(), 'k.', ms=2, alpha=0.4)

    ax_sed.fill_between(lam, 0, T_lum * 0.5, alpha=0.2, color='steelblue',
                        label='Lum filter (×0.5)')
    ax_sed.fill_between(lam, 0, T_ir  * 0.5, alpha=0.2, color='tomato',
                        label='IR filter (×0.5)')
    ax_sed.plot(lam, qe, 'g--', lw=1, alpha=0.7, label='QE')
    ax_sed.plot(lam, T_atm, color='gray', lw=1, ls=':', alpha=0.7, label='T_atm')
    ax_sed.set_ylabel('Normalised')
    ax_sed.set_title(f'{label}  (G={g_mag:.2f}, Teff={teff}K, {sed_label})')
    ax_sed.legend(fontsize=7)
    ax_sed.set_xlim(330, 1050)
    ax_sed.set_ylim(0, 1.1)

    # Integrand: F × T_atm × QE × T_filter  (what actually hits the detector)
    integrand_lum = F * T_atm * qe * T_lum
    integrand_ir  = F * T_atm * qe * T_ir
    norm = max(integrand_lum.max(), integrand_ir.max())

    ax_resp.plot(lam, integrand_lum / norm, color='steelblue', lw=1.5,
                 label=f'Lum  ∫={np.trapz(integrand_lum, lam):.2e}')
    ax_resp.plot(lam, integrand_ir  / norm, color='tomato',    lw=1.5,
                 label=f'IR   ∫={np.trapz(integrand_ir,  lam):.2e}')
    ax_resp.set_xlabel('Wavelength [nm]')
    ax_resp.set_ylabel('Normalised integrand')
    ax_resp.legend(fontsize=7)
    ax_resp.set_xlim(330, 1050)
    ax_resp.set_ylim(0, 1.1)


def main():
    sid_1190 = STARS['G2P1190R']['source_id']
    sid_1170 = STARS['G2P1170R']['source_id']

    params = fetch_params([sid_1190, sid_1170])

    wl_1190, fl_1190, wl_1170, fl_1170 = None, None, None, None

    if params.get(sid_1190, {}).get('has_xp'):
        wl_1190, fl_1190 = fetch_xp(sid_1190)
    print(f'G2P1190R XP: {"fetched" if wl_1190 is not None else "not available"}')

    if params.get(sid_1170, {}).get('has_xp'):
        wl_1170, fl_1170 = fetch_xp(sid_1170)
    print(f'G2P1170R XP: {"fetched" if wl_1170 is not None else "not available"}')

    p1190 = params.get(sid_1190, STARS['G2P1190R'])
    p1170 = params.get(sid_1170, STARS['G2P1170R'])

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle('Outlier star diagnostics: XP spectrum vs filter bandpasses', y=1.01)

    plot_star(axes[0, 0], axes[1, 0],
              'G2P1190R', wl_1190, fl_1190,
              g_mag=p1190['g_mag'], teff=p1190['teff'],
              bp_rp=p1190.get('bp_rp'), airmass=1.0, sed_type='XP')

    plot_star(axes[0, 1], axes[1, 1],
              'G2P1170R', wl_1170, fl_1170,
              g_mag=p1170['g_mag'], teff=p1170['teff'],
              bp_rp=p1170.get('bp_rp'), airmass=1.0, sed_type='Planck')

    plt.tight_layout()
    outpath = 'outlier_diagnostics.png'
    plt.savefig(outpath, dpi=120, bbox_inches='tight')
    print(f'\nPlot saved: {outpath}')


if __name__ == '__main__':
    main()
