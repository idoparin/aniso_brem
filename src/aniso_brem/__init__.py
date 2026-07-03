"""Python interface to the cseqs/flux Fortran bremsstrahlung routines."""

import numpy as np

from .cseqs import cseqs as _f
from . import _flux as _flux_ext
from ._flux import flux as _flux

__all__ = [
    "eebls",
    "eeblsnr",
    "eeblser",
    "eebapx",
    "eeba",
    "cspe",
    "cspax",
    "cs1",
    "cs1_integrated",
    "cs1_integrated_array",
    "eeba_integrated",
    "eeba_integrated_array",
    "integrate",
    "electron_flux_density",
    "electron_velocity_rel",
    "eeb_emin",
    "eeb_cross",
    "bremsstrahlung_thin_target",
    "bremsstrahlung_thin_target_eeb",
    "bremsstrahlung_thin_target_aniso",
    "bremsstrahlung_thin_target_eeb_aniso",
    "bremsstrahlung_thin_target_eeb_numeric",
    "bremsstrahlung_thin_target_ei_numeric",
    "pitch_angle_distribution_array",
    "ADF_ISOTROPIC",
    "ADF_GAUSSIAN",
    "pi",
    "r0",
    "c",
    "ech",
    "me",
    "afs",
    "AU",
]

AU = 1.496e13  # astronomical unit, cm (SSW IDL / sunxspex)

eebls = _f.eebls
eeblsnr = _f.eeblsnr
eeblser = _f.eeblser
eebapx = _f.eebapx
eeba = _f.eeba
cspe = _f.cspe
cspax = _f.cspax
cs1 = _f.cs1
integrate = _f.integrate
cs1_integrated = _flux_ext.flux_aniso.cs1_integrated
cs1_integrated_array = _flux_ext.flux_aniso.cs1_integrated_array
eeba_integrated = _flux_ext.flux_aniso.eeba_integrated
eeba_integrated_array = _flux_ext.flux_aniso.eeba_integrated_array

_aniso = _flux_ext.flux_aniso
ADF_ISOTROPIC = int(_aniso.adf_isotropic)
ADF_GAUSSIAN = int(_aniso.adf_gaussian)

pi = _f.pi
r0 = _f.r0
c = _f.c
ech = _f.ech
me = _f.me
afs = _f.afs

# Map a user-facing integrator selector to the Fortran integer code.
_INTEGRATOR_CODES = {
    "trapezoid": 0,
    "trap": 0,
    "trapz": 0,
    "gl": 1,
    "gauss": 1,
    "gauss_legendre": 1,
    "gauss-legendre": 1,
}


def _integrator_code(integrator):
    if integrator is None:
        return 1  # Gauss-Legendre is the default for this isotropic case
    if isinstance(integrator, str):
        try:
            return _INTEGRATOR_CODES[integrator.lower()]
        except KeyError:
            raise ValueError(
                f"Unknown integrator {integrator!r}; "
                f"choose one of {sorted(_INTEGRATOR_CODES)}"
            )
    raise TypeError("integrator must be None or one of "
                    f"{sorted(_INTEGRATOR_CODES)}")


def electron_flux_density(ee, p, eebrk, q, elow, ehigh):
    """Normalized broken power-law electron flux density.

    Accepts either a scalar or a numpy array for ``ee``. For arrays the loop
    runs in compiled Fortran (a single call), so passing a whole array is far
    faster than calling the scalar form element by element from Python.
    """
    arr = np.asarray(ee, dtype=np.float64)
    if arr.ndim == 0:
        return _flux.electron_flux_density(float(arr), p, eebrk, q, elow, ehigh)
    return _flux.electron_flux_density_array(
        np.ascontiguousarray(arr.ravel()), p, eebrk, q, elow, ehigh
    ).reshape(arr.shape)


def electron_velocity_rel(energy_kev):
    """Relativistic electron speed (cm/s) from kinetic energy (keV).

    ``v = c * sqrt(1 - 1/gamma^2)`` with ``gamma = 1 + E/mc2``. Accepts a scalar
    or numpy array (vectorized in numpy)."""
    g = 1.0 + np.asarray(energy_kev, dtype=np.float64) / me
    return np.sqrt(1.0 - (1.0 / g) ** 2) * c


def eeb_emin(photon_energy):
    """Minimum electron kinetic energy (keV) that can emit a photon of energy
    ``photon_energy`` in isotropic e-e bremsstrahlung (Haug 1975, eq. 2.6).

    Accepts a scalar or numpy array (vectorized in numpy)."""
    ep = np.asarray(photon_energy, dtype=np.float64)
    k = ep / me
    s = np.sqrt(k * k + 4.0 * k)
    return ep * (2.0 + 3.0 * k - s) / (1.0 - k * k + k * s)


def eeb_cross(electron_energy, photon_energy, z=1.2, ep_switch=150.0, xparts=32):
    """Isotropic e-e bremsstrahlung d(sigma)/dk [cm^2], scaled by z.

    Uses ``eebapx`` for photon energy below 150 keV and ``eebls`` at/above 150 keV
    (hard switch; ``eebapx`` diverges at higher photon energy). The ``ep_switch``
    argument is accepted for API compatibility but is not used.
    """
    return _flux.eeb_cross(electron_energy, photon_energy, z, ep_switch, xparts)


def bremsstrahlung_thin_target(
    photon_energies,
    p,
    break_energy,
    q,
    low_e_cutoff,
    high_e_cutoff,
    efd=True,
    integrator=None,
    *,
    z=1.2,
    npts=None,
    distance_cm=None,
):
    """Isotropic thin-target e-i bremsstrahlung flux (Fortran-backed).

    Same leading positional signature as ``sunkit_spex`` /
    ``sunxspex.emission.bremsstrahlung_thin_target``. Includes geometric
    spreading ``1 / (4 pi R^2)`` with ``R = distance_cm`` (default 1 AU), as in
    the sunxspex fork (sunkit-spex currently omits this factor).

    Parameters
    ----------
    photon_energies : array_like
        Photon energies to evaluate the flux at (keV).
    p, break_energy, q, low_e_cutoff, high_e_cutoff : float
        Broken power-law electron-distribution parameters (energies in keV).
    efd : bool
        ``True`` (default) input is an electron *flux* density distribution;
        ``False`` input is an electron *number* density distribution (the
        integrand is then multiplied by the relativistic velocity v(E) to
        convert density -> flux, matching the sunxspex ``efd`` convention).
    integrator : {None, "gl", "trapezoid"}, optional
        Quadrature rule. ``None``/``"gl"`` (default) uses fixed-order
        Gauss-Legendre in log-energy; ``"trapezoid"`` uses the trapezoidal rule
        on a log-spaced grid (call only if needed).
    z : float, optional
        Mean ion charge (sunkit-spex uses 1.2).
    npts : int, optional
        Quadrature points per energy segment. Defaults to 32 for Gauss-Legendre
        and 300 for the trapezoid.
    distance_cm : float, optional
        Source-observer distance in cm (default ``AU`` = 1.496e13).

    Returns
    -------
    numpy.ndarray
        Photon flux (photons s^-1 keV^-1 cm^-2) at the observer, per unit
        normalization coefficient (multiply by a_0 or n_ion n_nth V).
    """
    code = _integrator_code(integrator)
    if npts is None:
        npts = 32 if code == 1 else 300
    if distance_cm is None:
        distance_cm = AU
    photon_energies = np.ascontiguousarray(photon_energies, dtype=np.float64)
    return _flux.bremsstrahlung_thin_target(
        photon_energies,
        p,
        break_energy,
        q,
        low_e_cutoff,
        high_e_cutoff,
        z,
        1 if efd else 0,
        code,
        npts,
        distance_cm,
    )


def bremsstrahlung_thin_target_eeb(
    photon_energies,
    p,
    break_energy,
    q,
    low_e_cutoff,
    high_e_cutoff,
    efd=True,
    integrator=None,
    *,
    z=1.2,
    npts=None,
    ep_switch=150.0,
    xparts=32,
    distance_cm=None,
):
    """Isotropic thin-target electron-electron (e-e) bremsstrahlung flux.

    Same machinery, units and normalization as :func:`bremsstrahlung_thin_target`
    (so the two can be summed). Uses ``eebapx`` strictly below 150 keV and ``eebls``
    at/above 150 keV (hard switch fixed in Fortran; ``eebapx`` diverges above
    ~150 keV). Both are d sigma / dk (k = epsilon_gamma / m_e).

    Parameters
    ----------
    photon_energies : array_like
        Photon energies to evaluate the flux at (keV).
    p, break_energy, q, low_e_cutoff, high_e_cutoff : float
        Broken power-law electron-distribution parameters (energies in keV).
    efd : bool
        ``True`` (default) electron flux density input, ``False`` electron number
        density input (integrand multiplied by v(E) to convert density -> flux).
    integrator : {None, "gl", "trapezoid"}, optional
        Quadrature rule (default Gauss-Legendre).
    z : float, optional
        Mean ion charge. The e-e cross section is scaled linearly by ``z`` (Z
        electrons per ion), so e-e ~ Z while e-i ~ Z**2 on a per-ion basis.
    npts : int, optional
        Quadrature points per segment (default 32 for GL, 300 for trapezoid).
    ep_switch : float, optional
        Accepted for API compatibility; the Fortran switch is fixed at 150 keV.
    xparts : int, optional
        Internal integration sub-steps used by ``eebls`` (default 32).
    distance_cm : float, optional
        Source-observer distance in cm (default ``AU`` = 1.496e13).

    Returns
    -------
    numpy.ndarray
        e-e photon flux (photons s^-1 keV^-1 cm^-2) at the observer, per unit
        normalization coefficient.
    """
    code = _integrator_code(integrator)
    if npts is None:
        npts = 32 if code == 1 else 300
    if distance_cm is None:
        distance_cm = AU
    photon_energies = np.ascontiguousarray(photon_energies, dtype=np.float64)
    return _flux.bremsstrahlung_thin_target_eeb(
        photon_energies,
        p,
        break_energy,
        q,
        low_e_cutoff,
        high_e_cutoff,
        z,
        1 if efd else 0,
        code,
        npts,
        ep_switch,
        xparts,
        distance_cm,
    )


def pitch_angle_distribution_array(n_theta, dist_type=ADF_ISOTROPIC, sigma=0.5):
    """Normalized pitch-angle distribution g(theta) on [0, pi].

    Returns ``theta`` [rad], ``g_norm``, and ``solid_angle_norm``. The
    distribution satisfies
    ``integral_0^2pi d phi integral_0^pi g_norm(theta) sin(theta) d theta = 1``.
    Phi-symmetric distributions only.
    """
    return _aniso.pitch_angle_distribution_array(n_theta, dist_type, sigma)


def bremsstrahlung_thin_target_aniso(
    photon_energies,
    p,
    break_energy,
    q,
    low_e_cutoff,
    high_e_cutoff,
    efd=True,
    *,
    z=1.2,
    npts=32,
    alpha_rad=0.0,
    n_theta=32,
    n_phi=32,
    dist_type=ADF_GAUSSIAN,
    sigma=0.5,
    distance_cm=None,
):
    """Anisotropic thin-target e-i bremsstrahlung flux.

    Integrates the double-differential cross section ``cs1`` over electron
    pitch angles weighted by a normalized distribution ``g(theta)``, for an
    observer at angle ``alpha_rad`` from the solar-surface normal. Electron
    energy is integrated with Gauss-Legendre quadrature (``npts`` nodes).
    """
    if distance_cm is None:
        distance_cm = AU
    photon_energies = np.ascontiguousarray(photon_energies, dtype=np.float64)
    return _aniso.bremsstrahlung_thin_target_aniso(
        photon_energies,
        p,
        break_energy,
        q,
        low_e_cutoff,
        high_e_cutoff,
        z,
        1 if efd else 0,
        npts,
        alpha_rad,
        n_theta,
        n_phi,
        dist_type,
        sigma,
        distance_cm,
    )


def bremsstrahlung_thin_target_eeb_aniso(
    photon_energies,
    p,
    break_energy,
    q,
    low_e_cutoff,
    high_e_cutoff,
    efd=True,
    *,
    z=1.2,
    npts=32,
    alpha_rad=0.0,
    n_theta=32,
    n_phi=32,
    dist_type=ADF_GAUSSIAN,
    sigma=0.5,
    distance_cm=None,
):
    """Anisotropic thin-target e-e bremsstrahlung flux.

    Same as :func:`bremsstrahlung_thin_target_aniso` but uses ``eeba`` with the
    direction-dependent lower electron-energy limit (Haug 1975 eq. 3.3).
    """
    if distance_cm is None:
        distance_cm = AU
    photon_energies = np.ascontiguousarray(photon_energies, dtype=np.float64)
    return _aniso.bremsstrahlung_thin_target_eeb_aniso(
        photon_energies,
        p,
        break_energy,
        q,
        low_e_cutoff,
        high_e_cutoff,
        z,
        1 if efd else 0,
        npts,
        alpha_rad,
        n_theta,
        n_phi,
        dist_type,
        sigma,
        distance_cm,
    )


def bremsstrahlung_thin_target_eeb_numeric(
    photon_energies,
    p,
    break_energy,
    q,
    low_e_cutoff,
    high_e_cutoff,
    efd=True,
    *,
    z=1.2,
    npts=32,
    n_mu=64,
    distance_cm=None,
):
    """Isotropic e-e thin-target flux via ``eeba_integrated`` (numeric angle average).

    Validation reference for :func:`bremsstrahlung_thin_target_eeb_aniso` with
    ``ADF_ISOTROPIC`` at ``alpha_rad = 0``. Prefer this over
    :func:`bremsstrahlung_thin_target_eeb` (``eebls``) when checking the
    theta-phi integration of ``eeba``.
    """
    if distance_cm is None:
        distance_cm = AU
    photon_energies = np.ascontiguousarray(photon_energies, dtype=np.float64)
    return _aniso.bremsstrahlung_thin_target_eeb_numeric(
        photon_energies,
        p,
        break_energy,
        q,
        low_e_cutoff,
        high_e_cutoff,
        z,
        1 if efd else 0,
        npts,
        n_mu,
        distance_cm,
    )


def bremsstrahlung_thin_target_ei_numeric(
    photon_energies,
    p,
    break_energy,
    q,
    low_e_cutoff,
    high_e_cutoff,
    efd=True,
    *,
    z=1.2,
    npts=32,
    n_mu=64,
    distance_cm=None,
):
    """Isotropic e-i thin-target flux via ``cs1_integrated`` (numeric angle average).

    Validation reference for :func:`bremsstrahlung_thin_target_aniso` with
    ``ADF_ISOTROPIC`` at ``alpha_rad = 0``. Prefer this over
    :func:`bremsstrahlung_thin_target` (``cspax``) when checking the
    theta-phi integration of ``cs1``.
    """
    if distance_cm is None:
        distance_cm = AU
    photon_energies = np.ascontiguousarray(photon_energies, dtype=np.float64)
    return _aniso.bremsstrahlung_thin_target_ei_numeric(
        photon_energies,
        p,
        break_energy,
        q,
        low_e_cutoff,
        high_e_cutoff,
        z,
        1 if efd else 0,
        npts,
        n_mu,
        distance_cm,
    )
