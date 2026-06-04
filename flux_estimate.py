import numpy as np

_h  = 6.626e-34   # J s
_c  = 2.998e8     # m/s


def xp_photon_flux_density(wavelength_nm, xp_wavelength_nm, xp_flux_W):
    """
    Convert a Gaia XP sampled spectrum to photon flux density and interpolate
    onto the requested wavelength grid.

    Parameters
    ----------
    wavelength_nm : array-like
        Target wavelength grid [nm].
    xp_wavelength_nm : array-like
        Wavelengths from the XP table [nm], typically 336–1020 nm in 2 nm steps.
    xp_flux_W : array-like
        Spectral flux density from the XP table [W/m²/nm].

    Returns
    -------
    F : ndarray
        Photon flux density [ph/s/m²/nm] on wavelength_nm, zero outside
        the XP wavelength range.
    """
    lam_nm  = np.asarray(wavelength_nm, dtype=float)
    lam_m   = lam_nm * 1e-9
    # Energy per photon at each target wavelength [J]
    E_ph    = _h * _c / lam_m                         # J/photon
    # Interpolate W/m²/nm onto target grid; zero outside XP range
    F_W     = np.interp(lam_nm, xp_wavelength_nm, xp_flux_W, left=0.0, right=0.0)
    # ph/s/m²/nm = (W/m²/nm) / (J/photon)
    return F_W / E_ph


def photon_flux_density(wavelength_nm, phot_g_mean_mag, teff_gspphot, bp_rp,
                        F0_ref=1.015e8, lam_ref_nm=550.0, T_A0V=9600.0):
    """
    Estimate photon flux density [ph/s/m²/nm] at the top of atmosphere
    from Gaia DR3 photometry.

    The spectrum is modelled as a Planck function at T_eff, normalised so
    that the flux at the Gaia G effective wavelength matches the G-band
    magnitude relative to a zero-magnitude A0V (Vega) reference.

    Parameters
    ----------
    wavelength_nm : array-like
        Wavelength(s) in nm.
    phot_g_mean_mag : float
        Gaia G-band mean magnitude.
    teff_gspphot : float
        Effective temperature from Gaia GSP-Phot [K].
    bp_rp : float
        Gaia BP-RP colour index [mag].
    F0_ref : float
        Photon flux density of a 0-mag A0V star at lam_ref_nm [ph/s/m²/nm].
        Default 1.015e8 at 550 nm.
    lam_ref_nm : float
        Reference wavelength for F0_ref [nm].
    T_A0V : float
        Planck temperature of the A0V reference star [K].  Default 9600 K.

    Returns
    -------
    F : ndarray
        Photon flux density [ph/s/m²/nm].

    Notes
    -----
    Accuracy is limited by the blackbody approximation (~10-30% across the
    optical).  Stars with has_xp_sampled=True in Gaia DR3 have a directly
    calibrated low-resolution spectrum available via the gaiaxpy package.
    """
    h = 6.626e-34  # J s
    c = 2.998e8    # m/s
    k = 1.381e-23  # J/K

    lam = np.asarray(wavelength_nm, dtype=float) * 1e-9  # nm -> m
    lam_ref = lam_ref_nm * 1e-9

    def B_ph(T, lam_m):
        """Planck photon spectral density [ph/s/m³] per metre of wavelength."""
        x = h * c / (lam_m * k * T)
        return 2.0 * c / lam_m**4 / np.expm1(x)

    # G-band effective wavelength varies with star colour (Weiler 2018, A&A 617, A138)
    lam_G_eff = (621.0 + 23.0 * bp_rp) * 1e-9  # m

    # Zero-point photon flux at the G effective wavelength, derived from the
    # A0V reference value at lam_ref by scaling along the A0V Planck spectrum
    F0_G = F0_ref * B_ph(T_A0V, lam_G_eff) / B_ph(T_A0V, lam_ref)

    # This star's flux at the G effective wavelength
    F_G = F0_G * 10.0**(-0.4 * phot_g_mean_mag)

    # Full spectrum: Planck at T_eff, normalised to F_G at lam_G_eff
    return F_G * B_ph(teff_gspphot, lam) / B_ph(teff_gspphot, lam_G_eff)
