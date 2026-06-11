"""
Plot the photon transfer curve (PTC) for the ZWO ASI585MM Pro at GAIN=252.

Two panels:
  Top:    shot-noise variance vs mean signal (classic PTC) with linear fit
  Bottom: derived EGAIN vs exposure time

Saves: ptc_plot.png
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from collections import defaultdict

# ---------------------------------------------------------------------------
# Load bias frames
# ---------------------------------------------------------------------------

with fits.open('bias585/SNAPSHOT_0.00_0000.fits') as h:
    bias1 = h[0].data.astype(np.float64)
with fits.open('bias585/SNAPSHOT_0.00_0000(1).fits') as h:
    bias2 = h[0].data.astype(np.float64)

bias_level = (bias1 + bias2).mean() / 2
var_read   = np.var(bias1 - bias2) / 2

# ---------------------------------------------------------------------------
# Group flat pairs by EXPTIME header
# ---------------------------------------------------------------------------

by_exptime = defaultdict(list)
for f in sorted(glob.glob('flat585/*.fits')):
    exptime = float(fits.getheader(f).get('EXPTIME', -1))
    by_exptime[exptime].append(f)

exptimes  = []
signals   = []
varshots  = []
egains    = []

for exptime in sorted(by_exptime):
    files = sorted(by_exptime[exptime])[:2]
    if len(files) < 2:
        continue
    with fits.open(files[0]) as h:
        f1 = h[0].data.astype(np.float64)
    with fits.open(files[1]) as h:
        f2 = h[0].data.astype(np.float64)

    signal   = (f1 + f2).mean() / 2 - bias_level
    var_shot = np.var(f1 - f2) / 2 - var_read
    if var_shot <= 0:
        continue

    exptimes.append(exptime * 1000)   # s → ms
    signals.append(signal)
    varshots.append(var_shot)
    egains.append(signal / var_shot)

exptimes = np.array(exptimes)
signals  = np.array(signals)
varshots = np.array(varshots)
egains   = np.array(egains)

# Linear fit to PTC: var = signal / EGAIN  →  slope = 1/EGAIN
slope, intercept = np.polyfit(signals, varshots, 1)
egain_fit = 1.0 / slope

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(2, 1, figsize=(7, 8))

# --- Top: classic PTC -------------------------------------------------------
ax = axes[0]
ax.scatter(signals, varshots, s=50, zorder=5, label='Flat pairs')
sig_line = np.array([0, signals.max() * 1.05])
ax.plot(sig_line, slope * sig_line + intercept, 'k--', lw=1.5,
        label=f'Linear fit  →  EGAIN = {egain_fit:.4f} e⁻/stored ADU')
ax.set_xlabel('Mean signal (stored ADU above bias)', fontsize=11)
ax.set_ylabel('Shot-noise variance (stored ADU²)', fontsize=11)
ax.set_title('Photon Transfer Curve — ZWO ASI585MM Pro  GAIN=252', fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

# Annotate each point with its exposure time
for x, y, t in zip(signals, varshots, exptimes):
    ax.annotate(f'{t:.0f} ms', (x, y), textcoords='offset points',
                xytext=(5, 4), fontsize=7, color='steelblue')

# --- Bottom: EGAIN vs exposure time ----------------------------------------
ax2 = axes[1]
ax2.scatter(exptimes, egains, s=50, zorder=5)
ax2.axhline(np.median(egains), color='k', ls='--', lw=1.5,
            label=f'Median  {np.median(egains):.4f} e⁻/stored ADU'
                  f'  =  {np.median(egains)*16:.3f} e⁻/true ADU')
ax2.axhspan(np.median(egains) - np.std(egains),
            np.median(egains) + np.std(egains),
            alpha=0.15, color='steelblue', label=f'±1σ  ({np.std(egains):.4f})')
ax2.set_xlabel('Exposure time (ms)', fontsize=11)
ax2.set_ylabel('EGAIN (e⁻/stored ADU)', fontsize=11)
ax2.set_title('EGAIN vs exposure time', fontsize=11)
ax2.set_ylim(0, egains.max() * 1.4)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('ptc_plot.png', dpi=150)
print('Saved ptc_plot.png')
