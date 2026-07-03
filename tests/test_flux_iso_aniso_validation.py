"""Validate isotropic analytic flux against isotropic numeric (alpha=0) integration.

Mirrors ``examples/flux_aniso_example.py``: compare thin-target bremsstrahlung flux from

* **Analytic:** ``cspax`` (e-i); ``eebapx`` below 150 keV / ``eebls`` at/above (e-e)
* **Numeric:** ``bremsstrahlung_thin_target_aniso`` / ``_eeb_aniso`` at ``alpha_rad=0``
  with isotropic pitch-angle ADF.  Angle averaging uses split Gauss-Legendre quadrature
  on ``mu = cos(psi)`` (see ``int_cs1_mu`` / ``int_eeba_mu`` in ``flux_aniso.f95``).
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aniso_brem import (
    ADF_ISOTROPIC,
    AU,
    bremsstrahlung_thin_target,
    bremsstrahlung_thin_target_aniso,
    bremsstrahlung_thin_target_eeb,
    bremsstrahlung_thin_target_eeb_aniso,
)

# Defaults shared with examples/flux_aniso_example.py
P = 3.0
Q = 3.0
BREAK_ENERGY = 1000.0
LOW_E_CUTOFF = 10.0
HIGH_E_CUTOFF = 100_000.0
Z = 1.2
EEB_SWITCH_KEV = 150.0
NPTS = 32
XPARTS = 32
N_THETA = 128
N_PHI = 128

KW_ISO = dict(efd=True, z=Z, npts=NPTS, distance_cm=AU)
KW_EE = dict(**KW_ISO, ep_switch=EEB_SWITCH_KEV, xparts=XPARTS)
KW_ANISO = dict(
    efd=True,
    z=Z,
    npts=NPTS,
    distance_cm=AU,
    n_theta=N_THETA,
    n_phi=N_PHI,
    dist_type=ADF_ISOTROPIC,
    sigma=0.5,
    alpha_rad=0.0,
)


def _iso_analytic(eph, high_e_cutoff=HIGH_E_CUTOFF):
    params = (P, BREAK_ENERGY, Q, LOW_E_CUTOFF, high_e_cutoff)
    ei = bremsstrahlung_thin_target(eph, *params, **KW_ISO)
    ee = bremsstrahlung_thin_target_eeb(eph, *params, **KW_EE)
    return ei, ee, ei + ee


def _iso_numeric(eph, high_e_cutoff=HIGH_E_CUTOFF):
    params = (P, BREAK_ENERGY, Q, LOW_E_CUTOFF, high_e_cutoff)
    ei = bremsstrahlung_thin_target_aniso(eph, *params, **KW_ANISO)
    ee = bremsstrahlung_thin_target_eeb_aniso(eph, *params, **KW_ANISO)
    return ei, ee, ei + ee


@pytest.mark.parametrize(
    "high_e_cutoff, emin, emax, nsteps, rtol, band_label",
    [
        (1500.0, 10.0, 1400.0, 30, 0.13, "10-1400 keV, 1.5 MeV electron cutoff"),
        (HIGH_E_CUTOFF, 10.0, 150.0, 20, 0.15, "10-150 keV photons"),
        (
            HIGH_E_CUTOFF,
            EEB_SWITCH_KEV,
            30_000.0,
            24,
            0.05,
            "150 keV-30 MeV photons, 100 MeV electron cutoff",
        ),
    ],
)
def test_isotropic_analytic_vs_numeric_total_flux(
    high_e_cutoff, emin, emax, nsteps, rtol, band_label
):
    eph = np.logspace(np.log10(emin), np.log10(emax), nsteps)
    _, _, flux_a = _iso_analytic(eph, high_e_cutoff)
    _, _, flux_n = _iso_numeric(eph, high_e_cutoff)
    mask = flux_a > 0
    assert np.any(mask), f"no analytic flux in band: {band_label}"
    assert_allclose(
        flux_n[mask],
        flux_a[mask],
        rtol=rtol,
        atol=0.0,
        err_msg=f"total flux mismatch ({band_label})",
    )


def test_example_default_per_component_flux():
    """Same parameters and photon grid as flux_aniso_example.py defaults."""
    eph = np.logspace(np.log10(LOW_E_CUTOFF), np.log10(30_000.0), 40)
    ei_a, ee_a, tot_a = _iso_analytic(eph)
    ei_n, ee_n, tot_n = _iso_numeric(eph)

    mask_ei = ei_a > 0
    mask_ee = ee_a > 0
    mask_tot = tot_a > 0

    assert_allclose(ei_n[mask_ei], ei_a[mask_ei], rtol=0.13, atol=0.0)
    assert_allclose(ee_n[mask_ee], ee_a[mask_ee], rtol=0.18, atol=0.0)
    assert_allclose(tot_n[mask_tot], tot_a[mask_tot], rtol=0.13, atol=0.0)

    # Above the e-e cross-section switch, agreement is tighter except near the
    # eebapx -> eebls transition (~150-300 keV photon energy).
    above = eph >= EEB_SWITCH_KEV
    assert_allclose(
        ei_n[above & mask_ei],
        ei_a[above & mask_ei],
        rtol=0.03,
        atol=0.0,
    )
    assert_allclose(
        ee_n[above & mask_ee],
        ee_a[above & mask_ee],
        rtol=0.18,
        atol=0.0,
    )
    assert_allclose(
        tot_n[above & mask_tot],
        tot_a[above & mask_tot],
        rtol=0.05,
        atol=0.0,
    )


def test_ee_contributes_near_half_of_ei_in_mev_range():
    """e-e should be a substantial fraction of e-i above ~1 MeV (order 0.3-0.7)."""
    high_e_cutoff = 32_000.0
    eph = np.array([1000.0, 2000.0, 3000.0, 5000.0])
    ei_a, ee_a, _ = _iso_analytic(eph, high_e_cutoff)
    ratio = ee_a / ei_a
    assert np.all(ratio > 0.15), f"e-e/e-i too small: {ratio}"
    assert np.all(ratio < 1.5), f"e-e/e-i too large: {ratio}"
    assert 0.25 < np.median(ratio) < 0.85, f"median e-e/e-i = {np.median(ratio):.3f}"


def test_analytic_eeb_flux_uses_eebls_above_150_kev():
    """Fortran e-e analytic flux ignores ep_switch and always switches at 150 keV."""
    eph = np.array([1000.0, 2000.0, 5000.0])
    params = (P, BREAK_ENERGY, Q, LOW_E_CUTOFF, HIGH_E_CUTOFF)
    flux_ok = bremsstrahlung_thin_target_eeb(eph, *params, **KW_EE)
    kw_bad = dict(KW_EE, ep_switch=10_000.0)
    flux_fixed = bremsstrahlung_thin_target_eeb(eph, *params, **kw_bad)
    assert_allclose(flux_fixed, flux_ok, rtol=1e-10, atol=0.0)
