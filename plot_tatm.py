"""
Plot atmospheric transmission T_atm(λ) at MWO altitude (1740 m) for
airmass values representative of the 2025-10-16 survey, together with
the two filter transmission curves and the component extinction terms.

Saves: tatm_plot.png
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from telescope_transmissivity import (
    atmospheric_transmission, filter_transmission, qe_asi585
)

lam = np.linspace(330, 1050, 2000)

# --- Atmospheric components at airmass=1, altitude=1740 m ----------------
lam_um   = lam / 1000.0
p_ratio  = np.exp(-1740 / 8500.0)
k_ray    = 0.0094 * (lam_um / 0.4)**(-4) * p_ratio
k_aer    = 0.04   * (lam_um / 0.5)**(-1.3)
k_oz     = 0.016  * np.exp(-0.5 * ((lam - 600) / 55)**2)

T_ray  = 10**(-0.4 * k_ray)
T_aer  = 10**(-0.4 * k_aer)
T_oz   = 10**(-0.4 * k_oz)
T_tot1 = atmospheric_transmission(lam, airmass=1.0)
T_tot15 = atmospheric_transmission(lam, airmass=1.5)

# --- Filter and QE curves ------------------------------------------------
T_lum = filter_transmission(lam, 'luminance')
T_ir  = filter_transmission(lam, 'ir_longpass')
qe    = qe_asi585(lam)

# -------------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

# --- Top panel: T_atm and components -------------------------------------
ax = axes[0]
ax.plot(lam, T_tot1,  'k',   lw=2,   label='$T_{atm}$ total  (X=1.0)')
ax.plot(lam, T_tot15, 'k--', lw=1.5, label='$T_{atm}$ total  (X=1.5)')
ax.plot(lam, T_ray,  color='steelblue',  lw=1, ls=':', label='Rayleigh')
ax.plot(lam, T_aer,  color='sienna',     lw=1, ls=':', label='Aerosol')
ax.plot(lam, T_oz,   color='olive',      lw=1, ls=':', label='Ozone (Chappuis)')

# Shade filter passbands lightly
ax.axvspan(400, 710,  alpha=0.08, color='royalblue', label='Luminance passband')
ax.axvspan(685, 1050, alpha=0.08, color='tomato',    label='IR longpass passband')

ax.set_ylabel('Transmission', fontsize=11)
ax.set_ylim(0.5, 1.02)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
ax.legend(fontsize=8, ncol=2, loc='lower right')
ax.set_title('Atmospheric transmission model — MWO (1740 m)', fontsize=11)
ax.grid(True, alpha=0.3)

# Annotate Chappuis band minimum
i_min = np.argmin(T_oz)
ax.annotate('Chappuis\nO₃ band',
            xy=(lam[i_min], T_oz[i_min]),
            xytext=(560, 0.82),
            fontsize=7, color='olive',
            arrowprops=dict(arrowstyle='->', color='olive', lw=0.8))

# --- Bottom panel: filter × QE × T_atm integrand weights ----------------
ax2 = axes[1]
ax2.plot(lam, T_lum * qe * T_tot1, color='royalblue', lw=2,
         label='Luminance × QE × $T_{atm}$')
ax2.plot(lam, T_ir  * qe * T_tot1, color='tomato',    lw=2,
         label='IR longpass × QE × $T_{atm}$')
ax2.plot(lam, qe,   'k:', lw=1, label='QE (ASI585MM Pro)')
ax2.plot(lam, T_lum, color='royalblue', lw=1, ls='--', alpha=0.5,
         label='Filter transmission')
ax2.plot(lam, T_ir,  color='tomato',    lw=1, ls='--', alpha=0.5)

ax2.set_xlabel('Wavelength (nm)', fontsize=11)
ax2.set_ylabel('Fraction', fontsize=11)
ax2.set_ylim(0, 1.05)
ax2.legend(fontsize=8, loc='upper right')
ax2.set_title('Integrand weights: filter × QE × $T_{atm}$  (X=1.0)', fontsize=11)
ax2.grid(True, alpha=0.3)

# Representative T_atm values as a text table in the top panel
rows = [330, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000]
t1_vals  = atmospheric_transmission(np.array(rows, float), airmass=1.0)
t15_vals = atmospheric_transmission(np.array(rows, float), airmass=1.5)
table_str = 'λ(nm)  X=1.0  X=1.5\n' + '\n'.join(
    f'{r:5d}  {t1:.3f}  {t15:.3f}' for r, t1, t15 in zip(rows, t1_vals, t15_vals))
axes[0].text(0.013, 0.03, table_str, transform=axes[0].transAxes,
             fontsize=6.5, family='monospace', verticalalignment='bottom',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('tatm_plot.png', dpi=150)
print('Saved tatm_plot.png')
print()
print(table_str)
