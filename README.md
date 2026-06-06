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

Expected T_tel: **0.43 – 0.57**

## Results (2025-10-16)

Camera: ZWO ASI585MM Pro · GAIN = 252 (HCG mode) · EGAIN = 0.057 e-/ADU\*
Plate scale: 18.63 mas/px · Frame size: 1024 × 1024

| Filter | Bandpass | n stars | T_tel (median) | T_tel (mean ± std) |
|---|---|---|---|---|
| Baader UV/IR-cut L | 400–710 nm | 8 | 0.480 | 0.479 ± 0.013 |
| Baader 685 IR-pass | 685–1050 nm | 10 | 0.314 | 0.316 ± 0.008 |

\* The camera stores 12-bit ADC values left-shifted by 4 bits into 16-bit FITS
pixels (all stored values are multiples of 16; bias pedestal ~80 stored ADU =
5 true ADU).  EGAIN = 0.057 e-/stored ADU corresponds to **0.91 e-/true 12-bit
ADU**, consistent with the HCG plateau on the ZWO spec chart.  Derived from
T_tel self-consistency (luminance median → 0.480 in expected range).  Pending
confirmation by photon transfer curve from flat-field pairs.

The lower IR transmissivity is expected: aged aluminium reflectivity declines
in the near-IR.

Two stars were excluded from the summary statistics:

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
project server and are not included in this repository.

## Next steps

- Photon transfer curve from flat-field pairs to confirm EGAIN independently
- Extend pipeline to the MWO Hale 60-inch telescope
- Cross-check Baader 685 IR filter curve against manufacturer data
  (implied EGAIN from IR band is ~60% higher than from luminance)
