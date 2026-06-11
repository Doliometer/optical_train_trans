# Optical Train Transmissivity

Pipeline for measuring the end-to-end optical transmissivity of the
Mount Wilson Observatory Hooker 100-inch telescope speckle imaging train,
using Gaia DR3 flux-calibrated spectra as photometric references.

## Method

For each reference star:

```
T_tel = (observed e-/s) / (predicted e-/s)
```

where

```
predicted e-/s = ∫ F_λ(λ) · A_tel · T_atm(λ) · QE(λ) · T_filter(λ) dλ
```

| Term | Units | Description |
|---|---|---|
| `F_λ(λ)` | ph/s/m²/nm | Photon flux density at the top of atmosphere. Taken from the Gaia DR3 XP sampled spectrum (336–1020 nm, 2 nm steps, directly flux-calibrated) where available; otherwise a Planck blackbody normalised to the Gaia G magnitude. |
| `A_tel` | m² | Effective collecting area of the primary mirror, minus the secondary obstruction (π/4 × (2.54² − 0.36²) ≈ 4.96 m²). |
| `T_atm(λ)` | — | Atmospheric transmission. Analytic model combining Rayleigh scattering (scaled for MWO altitude 1740 m), grey aerosol, and the Chappuis ozone band (~600 nm). |
| `QE(λ)` | — | Detector quantum efficiency of the ZWO ASI585MM Pro, digitised from the manufacturer's published curve. |
| `T_filter(λ)` | — | Filter transmission, digitised from Baader manufacturer curves (UV/IR-cut L or 685 IR-pass). |

The integral yields a predicted rate in photo-electrons per second (e-/s):
each photon that reaches the detector has a wavelength-dependent probability
QE(λ) of liberating one electron from the silicon, and the integral sums
those probabilities weighted by the incoming photon flux. The observed e-/s
is obtained from the measured pixel counts — in ADU (Analog-to-Digital Units,
the integer values stored in the FITS file) — by multiplying by EGAIN
(electrons per ADU) and dividing by the exposure time.

`T_tel` absorbs all throughput not otherwise accounted for: mirror
reflectivities, spider obscuration losses, relay optics, and any vignetting.

## Optical train

The Hooker speckle train passes light through 6 mirrors and 1 relay lens:

| Element | Coating | Notes |
|---|---|---|
| Primary (2.54 m) | Al, recently recoated | |
| Secondary (~0.36 m obstruction) | Al, aged | |
| Tertiary fold | Al, aged | redirects beam into side bench |
| Relay tube mirror 1 | Al, aged | |
| Relay tube mirror 2 | Al, aged | |
| Relay refractor lens | — | broadband AR coated |
| Final fold before camera | Ag | manufacturer spec ~98% |

Expected T_tel:

| Passband | Range | Basis |
|---|---|---|
| Luminance (400–710 nm) | **0.43 – 0.57** | Al reflectivity ~91% (fresh), ~85–90% (aged) |
| IR longpass (685–1050 nm) | **0.28 – 0.38** | Al reflectivity ~84% (fresh), ~78–83% (aged) at ~850 nm effective wavelength (Rakić 1995) |

Al reflectivity declines significantly in the near-IR, driving the lower expected range for the IR passband.

## Camera calibration — EGAIN

ADU (Analog-to-Digital Unit) is the integer count reported by the camera.
The on-chip ADC produces a 12-bit integer (true ADU, range 0–4095). The ZWO
driver left-shifts this by 4 bits before writing to the FITS file, so the
stored value is always a multiple of 16 (stored ADU, range 0–65520). One
true ADU therefore equals 16 stored ADU. The bias pedestal is ~942 stored ADU
(~59 true ADU). EGAIN is expressed in stored units throughout the pipeline.

**Photon transfer curve (PTC):** 10 flat pairs at GAIN=252 in darkness, battery-powered
LED, exposure times 0.005–0.050 s:

| | |
|---|---|
| EGAIN | **0.0303 ± 0.0001 e-/stored ADU** (std/mean < 0.4%) |
| True EGAIN | 0.0303 × 16 = **0.485 e-/true 12-bit ADU** |
| Lou Jackson (priv. comm.) | 0.493–0.495 e-/true ADU at GAIN=252, HCG, 0–15°C |

FITS header EGAIN values are unreliable (ZWO SDK inconsistency) and should not be used.

## Results (2025-10-16)

Camera: ZWO ASI585MM Pro · GAIN = 252 (HCG mode) · EGAIN = 0.0303 e-/stored ADU
Plate scale: 18.63 mas/px · Frame size: 1024 × 1024 · 1000 frames/star/filter

| Filter | Bandpass | Expected T_tel | n stars | T_tel (median) | T_tel (mean ± std) |
|---|---|---|---|---|---|
| Baader UV/IR-cut L | 400–710 nm | 0.43 – 0.57 | 4 | 0.251 | 0.251 ± 0.004 |
| Baader 685 IR-pass | 685–1050 nm | 0.28 – 0.38 | 7 | 0.166 | 0.166 ± 0.004 |

Both results fall well below their respective expected ranges — roughly a factor of
2. The implied EGAIN to bring T_tel up to the mirror-budget midpoint is ~0.015
e-/stored ADU for both bands, inconsistent with the PTC. The discrepancy points to
unaccounted throughput losses in the optical train (spider vignetting, relay optics,
actual mirror reflectivities lower than budget values). Full output in
`survey_results.txt`. Runtime: ~81 minutes.

Stars are flagged reliable if the sky annulus inner radius ≥ 3 × Moffat FWHM
(β = 3.0); unreliable stars are excluded from summary statistics. Two additional
stars were excluded by name:

- **G2P1190R** (Teff = 4776 K, bp_rp = 1.52): very red star; luminance band
  samples the suppressed blue side of the SED, giving unreliable T_tel in
  both filters.
- **G2P1170R** (Teff = 7594 K): no XP spectrum available; Planck SED
  underpredicts IR flux, yielding T_tel > 1 in the IR longpass band.

## Scripts

| Script | Purpose |
|---|---|
| `flux_estimate.py` | Gaia XP and Planck SED → photon flux density [ph/s/m²/nm] |
| `telescope_transmissivity.py` | Single-star T_tel; pad-aware shift-and-add photometry |
| `survey_transmissivity.py` | Full multi-star survey with Gaia queries and reliability flagging |
| `growth_curve.py` | Encircled energy growth curves for all stars |
| `growth_moffat.py` | Moffat FWHM fits; pad-safe sky annulus placement |
| `moffat_fit.py` | Seeing disk characterisation across a random sample of stars |
| `plot_xp_outliers.py` | Diagnostic XP spectrum plots for excluded outlier stars |
| `plot_tatm.py` | Atmospheric transmission model: component breakdown and filter × QE × T_atm integrand weights |

## Dependencies

```
astropy
astroquery
scipy
numpy
matplotlib
```

Install with:

```bash
pip install astropy astroquery scipy numpy matplotlib
```

## Usage

Run the full 19-star survey:

```bash
python3 survey_transmissivity.py
```

Characterise the seeing disk (3 random stars, luminance filter):

```bash
python3 moffat_fit.py --seed 42
```

Growth curves with pad-safe sky annulus:

```bash
python3 growth_moffat.py
```

## Data

FITS frames (1000 per star per filter, 19 stars, 2 filters) are held on the
Stelar-group data server and are not included in this repository.

## Next steps

- Investigate factor-of-2 shortfall in T_tel relative to mirror budget:
  spider vignetting, relay optics losses, actual mirror reflectivities
- Cross-check Baader filter curves against manufacturer data
- Extend pipeline to the MWO Hale 60-inch telescope
