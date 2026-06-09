"""Regression tests for ``eebapx`` against Haug (1998) Table I.

``eebapx(ee, ep)`` implements Equation (2), the long-wavelength expansion of
Equation (1). Comparisons are restricted to photon energies <= 300 keV where
the approximation remains meaningful (cf. Table II in Haug 1998).

At E = 1 MeV use d sigma / d(h nu) = eebapx / me. At E = 10 MeV the same
``2 * eebapx`` scaling that matches Table I for ``eebls`` is used.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from table_ref_haug98 import (
    CASE1_ELECTRON_ENERGY_KEV,
    CASE1_PHOTON_ENERGIES_KEV,
    CASE1_REFERENCE_CM2_KEV,
    CASE2_ELECTRON_ENERGY_KEV,
    CASE2_PHOTON_ENERGIES_KEV,
    CASE2_REFERENCE_CM2_KEV,
    EEBAPX_MAX_PHOTON_KEV,
    RTOL,
    dsigma_d_hnu_eebapx,
    dsigma_d_hnu_eebls,
    relative_errors,
    values_for_photons,
)
from test_eebls import ATOL_CASE2, CASE2_EEBLS_TABLE_SCALE, eebls_case2_table_i

ATOL_EEBAPX_CASE1 = 1e-28


def eebapx_case1_table_i(electron_energy_keV, photon_energy_keV):
    return dsigma_d_hnu_eebapx(electron_energy_keV, photon_energy_keV)


def eebapx_case2_table_i(electron_energy_keV, photon_energy_keV):
    from aniso_brem import eebapx

    return CASE2_EEBLS_TABLE_SCALE * eebapx(electron_energy_keV, photon_energy_keV)


def _long_wavelength_mask(photon_energies_keV):
    return photon_energies_keV <= EEBAPX_MAX_PHOTON_KEV


def test_eebapx_case1_1mev_electron_below_300kev():
    """``eebapx`` vs Table I at E = 1 MeV for h nu <= 300 keV."""
    mask = _long_wavelength_mask(CASE1_PHOTON_ENERGIES_KEV)
    calculated = values_for_photons(
        CASE1_PHOTON_ENERGIES_KEV[mask],
        CASE1_ELECTRON_ENERGY_KEV,
        eebapx_case1_table_i,
    )
    assert_allclose(
        calculated,
        CASE1_REFERENCE_CM2_KEV[mask],
        rtol=RTOL,
        atol=ATOL_EEBAPX_CASE1,
    )


def test_eebapx_case2_10mev_electron_below_300kev():
    """``eebapx`` vs Table I at E = 10 MeV for h nu <= 300 keV."""
    mask = _long_wavelength_mask(CASE2_PHOTON_ENERGIES_KEV)
    calculated = values_for_photons(
        CASE2_PHOTON_ENERGIES_KEV[mask],
        CASE2_ELECTRON_ENERGY_KEV,
        eebapx_case2_table_i,
    )
    assert_allclose(
        calculated,
        CASE2_REFERENCE_CM2_KEV[mask],
        rtol=RTOL,
        atol=ATOL_CASE2,
    )


def test_eebapx_errors_vs_eebls_below_300kev():
    """Below 300 keV, ``eebapx`` tracks ``eebls`` and Table I with similar errors."""
    mask1 = _long_wavelength_mask(CASE1_PHOTON_ENERGIES_KEV)
    mask2 = _long_wavelength_mask(CASE2_PHOTON_ENERGIES_KEV)

    eebls1 = values_for_photons(
        CASE1_PHOTON_ENERGIES_KEV[mask1],
        CASE1_ELECTRON_ENERGY_KEV,
        dsigma_d_hnu_eebls,
    )
    apx1 = values_for_photons(
        CASE1_PHOTON_ENERGIES_KEV[mask1],
        CASE1_ELECTRON_ENERGY_KEV,
        eebapx_case1_table_i,
    )
    eebls2 = values_for_photons(
        CASE2_PHOTON_ENERGIES_KEV[mask2],
        CASE2_ELECTRON_ENERGY_KEV,
        eebls_case2_table_i,
    )
    apx2 = values_for_photons(
        CASE2_PHOTON_ENERGIES_KEV[mask2],
        CASE2_ELECTRON_ENERGY_KEV,
        eebapx_case2_table_i,
    )

    err_eebls1 = relative_errors(eebls1, CASE1_REFERENCE_CM2_KEV[mask1])
    err_apx1 = relative_errors(apx1, CASE1_REFERENCE_CM2_KEV[mask1])
    err_eebls2 = relative_errors(eebls2, CASE2_REFERENCE_CM2_KEV[mask2])
    err_apx2 = relative_errors(apx2, CASE2_REFERENCE_CM2_KEV[mask2])

    # E = 1 MeV, h nu <= 300 keV
    assert err_apx1.max() == pytest.approx(0.041504, rel=1e-3)
    assert err_eebls1.max() == pytest.approx(0.001156, rel=1e-3)
    assert err_apx1.max() < 40 * err_eebls1.max()
    assert np.all((apx1 / eebls1) > 0.95)
    assert np.all((apx1 / eebls1) < 1.05)

    # E = 10 MeV, h nu <= 300 keV
    assert err_apx2.max() == pytest.approx(0.022765, rel=1e-3)
    assert err_eebls2.max() == pytest.approx(0.022221, rel=1e-3)
    assert err_apx2.max() < 2 * err_eebls2.max()
    assert np.all((apx2 / eebls2) > 0.99)
    assert np.all((apx2 / eebls2) < 1.01)


def test_eebapx_diverges_above_300kev_case1():
    """Above 300 keV at 1 MeV, ``eebapx`` no longer tracks ``eebls`` or Table I."""
    high_photon = CASE1_PHOTON_ENERGIES_KEV > EEBAPX_MAX_PHOTON_KEV
    eebls_vals = values_for_photons(
        CASE1_PHOTON_ENERGIES_KEV[high_photon],
        CASE1_ELECTRON_ENERGY_KEV,
        dsigma_d_hnu_eebls,
    )
    apx_vals = values_for_photons(
        CASE1_PHOTON_ENERGIES_KEV[high_photon],
        CASE1_ELECTRON_ENERGY_KEV,
        eebapx_case1_table_i,
    )
    reference = CASE1_REFERENCE_CM2_KEV[high_photon]

    assert relative_errors(apx_vals, reference).max() > 10 * relative_errors(
        eebls_vals, reference
    ).max()
    assert np.max(apx_vals / eebls_vals) > 10
