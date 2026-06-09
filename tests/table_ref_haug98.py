"""Haug (1998, Solar Physics 178, Table I) reference cross sections.

Table I lists d sigma_ee / d(h nu) in millibarns per keV for kinetic electron
energies E = 1 MeV and E = 10 MeV. Fortran routines ``eebls`` and ``eebapx``
return d sigma / dk with k = h nu / (m_e c^2); photon and electron energies
passed to the code are in keV.

"""

import numpy as np

import aniso_brem

# Haug (1998) Table I: d sigma_ee / d(h nu) in mb / keV
MB_PER_KEV_TO_CM2_PER_KEV = 1e-27

CASE1_ELECTRON_ENERGY_KEV = 1_000.0
CASE1_PHOTON_ENERGIES_KEV = np.array(
    [0.1, 1, 5, 10, 20, 50, 100, 200, 400, 600, 800]
)
CASE1_TABLE_I_MB_KEV = np.array(
    [391.8, 30.3, 4.78, 2.10, 0.895, 0.268, 0.0952, 0.0264, 0.00585, 0.00176, 0.00012]
)

CASE2_ELECTRON_ENERGY_KEV = 10_000.0
CASE2_PHOTON_ENERGIES_KEV = np.array(
    [1, 5, 10, 50, 100, 500, 1000, 2000, 4000, 6000, 8000]
)
CASE2_TABLE_I_MB_KEV = np.array(
    [63539, 10815, 4986, 784, 339.5, 39.67, 15.10, 5.48, 1.752, 0.790, 0.340]
)

CASE1_REFERENCE_CM2_KEV = CASE1_TABLE_I_MB_KEV * MB_PER_KEV_TO_CM2_PER_KEV
CASE2_REFERENCE_CM2_KEV = CASE2_TABLE_I_MB_KEV * MB_PER_KEV_TO_CM2_PER_KEV

XPARTS = 5000
RTOL = 1e-5

# Valid domain for eebapx (Haug 1998 eq. 2): photon energies below this limit.
EEBAPX_MAX_PHOTON_KEV = 300.0


def dsigma_d_hnu_from_dsigma_dk(dsigma_dk):
    """Convert d sigma / dk (Fortran output) to d sigma / d(h nu) in cm^2 / keV."""
    return dsigma_dk / aniso_brem.me


def dsigma_d_hnu_eebls(electron_energy_keV, photon_energy_keV, xparts=XPARTS):
    """Table I d sigma / d(h nu) from ``eebls`` (Haug 1998 eq. 1)."""
    return dsigma_d_hnu_from_dsigma_dk(
        aniso_brem.eebls(electron_energy_keV, photon_energy_keV, xparts)
    )


def dsigma_d_hnu_eebapx(electron_energy_keV, photon_energy_keV):
    """Table I d sigma / d(h nu) from ``eebapx`` (Haug 1998 eq. 2)."""
    return dsigma_d_hnu_from_dsigma_dk(
        aniso_brem.eebapx(electron_energy_keV, photon_energy_keV)
    )


def table_i_values(electron_energy_keV, photon_energy_keV):
    """Reference d sigma / d(h nu) in cm^2 / keV from Table I."""
    if electron_energy_keV == CASE1_ELECTRON_ENERGY_KEV:
        photons = CASE1_PHOTON_ENERGIES_KEV
        reference = CASE1_REFERENCE_CM2_KEV
    elif electron_energy_keV == CASE2_ELECTRON_ENERGY_KEV:
        photons = CASE2_PHOTON_ENERGIES_KEV
        reference = CASE2_REFERENCE_CM2_KEV
    else:
        raise ValueError(f"unsupported electron energy {electron_energy_keV} keV")

    index = np.where(photons == photon_energy_keV)[0]
    if index.size != 1:
        raise ValueError(f"photon energy {photon_energy_keV} keV not in Table I")
    return reference[index[0]]


def values_for_photons(photon_energies_keV, electron_energy_keV, cross_section):
    return np.array(
        [cross_section(electron_energy_keV, photon_energy) for photon_energy in photon_energies_keV]
    )


def relative_errors(calculated, reference):
    return np.abs((calculated - reference) / reference)
