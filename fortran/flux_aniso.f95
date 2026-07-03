! flux_aniso.f95 -- anisotropic thin-target bremsstrahlung flux.
!
! Integrates double-differential cross sections (cs1 for e-i, eeba for e-e)
! over electron pitch-angle directions (theta, phi) weighted by a normalized
! pitch-angle distribution g(theta), for a fixed observer angle alpha between
! the line of sight and the solar-surface normal.


module flux_aniso

  use cseqs, only: cs1, eeba, integrate, pi, me, c
  use flux, only: electron_flux_density, electron_velocity_rel

  implicit none
  private
  public :: cs1_integrated, cs1_integrated_array
  public :: eeba_integrated, eeba_integrated_array
  public :: pitch_angle_distribution_array
  public :: cs1_aniso, eeba_aniso
  public :: bremsstrahlung_thin_target_aniso, bremsstrahlung_thin_target_eeb_aniso
  public :: bremsstrahlung_thin_target_eeb_numeric
  public :: bremsstrahlung_thin_target_ei_numeric

  real(8), parameter :: au_default = 1.496d13
  ! dist_type codes for pitch_angle_distribution_array / aniso flux routines
  integer, parameter, public :: adf_isotropic = 0, adf_gaussian = 1
  ! Split cos(psi) quadrature: cs1 / eeba peak sharply at mu -> 1; uniform
  ! trapezoid rules on [-1, 1] including mu = 1 over-estimate the integral.
  real(8), parameter :: mu_split = 0.99d0

contains

  real(8) function int_cs1_mu(ep, ee, z, nn, dist_type, sigma) result(cs)
    ! 2 pi int_{-1}^{1} cs1(ep, ee, mu, z) g(theta(mu)) d mu, times m_e  -> d(sigma)/dk.
    ! Gauss-Legendre on [-1, mu_split] plus mu = 1 - t^2 on [mu_split, 1].
    integer, intent(in) :: nn, dist_type
    real(8), intent(in) :: ep, ee, z, sigma
    real(8) :: gx(nn), gw(nn), s, val, mu, t, tmax, half, mid, gwt
    integer :: i

    cs = 0d0
    if (nn < 2 .or. ep <= 0d0 .or. ep >= ee) return

    call gauleg(nn, gx, gw)
    half = 0.5d0*(mu_split + 1d0)
    mid = 0.5d0*(mu_split - 1d0)
    s = 0d0
    do i = 1, nn
      mu = half*gx(i) + mid
      gwt = pitch_angle_distribution_raw(dacos(mu), dist_type, sigma)
      val = cs1(ep, ee, mu, z)
      if (val < 0d0 .or. val /= val) val = 0d0
      s = s + gw(i)*val*gwt*half
    end do

    tmax = dsqrt(1d0 - mu_split)
    half = 0.5d0*tmax
    do i = 1, nn
      t = half*(gx(i) + 1d0)
      mu = 1d0 - t*t
      gwt = pitch_angle_distribution_raw(dacos(mu), dist_type, sigma)
      val = cs1(ep, ee, mu, z)
      if (val < 0d0 .or. val /= val) val = 0d0
      s = s + gw(i)*val*gwt*2d0*t*half
    end do

    cs = 2d0*pi*s*me
  end function int_cs1_mu

  real(8) function int_eeba_mu(ep, ee, nn, dist_type, sigma) result(cs)
    ! 2 pi int eeba d mu (times m_e), with per-direction e-e thresholds.
    integer, intent(in) :: nn, dist_type
    real(8), intent(in) :: ep, ee, sigma
    real(8) :: gx(nn), gw(nn), s, val, mu, t, tmax, half, mid, gwt, emin
    integer :: i

    cs = 0d0
    if (nn < 2 .or. ep <= 0d0 .or. ep >= ee) return

    call gauleg(nn, gx, gw)
    half = 0.5d0*(mu_split + 1d0)
    mid = 0.5d0*(mu_split - 1d0)
    s = 0d0
    do i = 1, nn
      mu = half*gx(i) + mid
      gwt = pitch_angle_distribution_raw(dacos(mu), dist_type, sigma)
      emin = eeb_emin_aniso(ep, mu)
      if (ee < emin) then
        val = 0d0
      else
        val = eeba(ee, ep, mu)
        if (val < 0d0 .or. val /= val) val = 0d0
      end if
      s = s + gw(i)*val*gwt*half
    end do

    tmax = dsqrt(1d0 - mu_split)
    half = 0.5d0*tmax
    do i = 1, nn
      t = half*(gx(i) + 1d0)
      mu = 1d0 - t*t
      gwt = pitch_angle_distribution_raw(dacos(mu), dist_type, sigma)
      emin = eeb_emin_aniso(ep, mu)
      if (ee < emin) then
        val = 0d0
      else
        val = eeba(ee, ep, mu)
        if (val < 0d0 .or. val /= val) val = 0d0
      end if
      s = s + gw(i)*val*gwt*2d0*t*half
    end do

    cs = 2d0*pi*s*me
  end function int_eeba_mu

  pure real(8) function cos_psi(alpha, theta, phi) result(cpsi)
    ! Scattering-angle cosine (Bai & Ramaty 1978 eq. 16).
    real(8), intent(in) :: alpha, theta, phi
    cpsi = dcos(alpha)*dcos(theta) + dsin(alpha)*dsin(theta)*dcos(phi)
  end function cos_psi

  pure real(8) function pitch_angle_distribution_raw(theta, dist_type, sigma) result(g)
    ! Unnormalized pitch-angle weight g(theta); phi-symmetric (no phi factor).
    real(8), intent(in) :: theta, sigma
    integer, intent(in) :: dist_type
    select case (dist_type)
    case (adf_isotropic)
      g = 1d0
    case (adf_gaussian)
      if (sigma > 0d0) then
        g = dexp(-0.5d0*(theta/sigma)**2)
      else
        g = 0d0
      end if
    case default
      g = 1d0
    end select
  end function pitch_angle_distribution_raw

  subroutine pitch_angle_distribution_array(n_theta, dist_type, sigma, theta, g_norm, &
                                            solid_angle_norm)
    ! Normalized pitch-angle distribution on [0, pi] (for plotting / diagnostics).
    !
    ! g_norm = g_raw / solid_angle_norm with
    !   solid_angle_norm = 2*pi * integral g_raw(theta)*sin(theta) d theta
    ! so that integral_0^2pi d phi integral_0^pi g_norm sin(theta) d theta = 1.
    !
    ! Flux integrals in cs1_aniso / eeba_aniso use g_raw (legacy convention:
    ! unnormalized weight, isotropic => g = 1), not g_norm.
    integer, intent(in) :: n_theta, dist_type
    real(8), intent(in) :: sigma
    real(8), intent(out) :: theta(n_theta), g_norm(n_theta), solid_angle_norm
    integer :: i
    real(8), allocatable :: g_raw(:), integrand(:)

    if (n_theta < 2) then
      solid_angle_norm = 0d0
      theta = 0d0
      g_norm = 0d0
      return
    end if

    theta(1) = 0d0
    do i = 2, n_theta
      theta(i) = dble(i)*pi/dble(n_theta)
    end do

    allocate(g_raw(n_theta), integrand(n_theta))
    do i = 1, n_theta
      g_raw(i) = pitch_angle_distribution_raw(theta(i), dist_type, sigma)
      integrand(i) = g_raw(i)*dsin(theta(i))
    end do
    solid_angle_norm = 2d0*pi*integrate(theta, integrand)
    if (solid_angle_norm <= 0d0) then
      g_norm = 0d0
    else
      do i = 1, n_theta
        g_norm(i) = g_raw(i)/solid_angle_norm
      end do
    end if
    deallocate(g_raw, integrand)
  end subroutine pitch_angle_distribution_array

  pure real(8) function eeb_emin_iso(ep) result(emin)
    ! Isotropic e-e lower electron-energy limit (Haug 1975 eq. 2.6).
    real(8), intent(in) :: ep
    real(8) :: k, s
    k = ep/me
    s = dsqrt(k*k + 4d0*k)
    emin = ep*(2d0 + 3d0*k - s)/(1d0 - k*k + k*s)
  end function eeb_emin_iso

  pure real(8) function eeb_emin_aniso(ep, cpsi) result(emin)
    ! Direction-dependent e-e lower electron-energy limit (Haug 1975 eq. 3.3
    ! with theta -> psi; Oparin et al. 2020 p. 891).
    real(8), intent(in) :: ep, cpsi
    real(8) :: k, s, num, den
    k = ep/me
    s = dsqrt(k*k*cpsi*cpsi + 4d0*k)
    num = 2d0*(1d0 + k) + cpsi*(k*cpsi - s)
    den = 1d0 - k*k + k*cpsi*s
    if (den <= 0d0) then
      emin = ep
    else
      emin = ep*num/den
    end if
  end function eeb_emin_aniso

  real(8) function cs1_aniso(ee, ep, z, alpha, n_theta, n_phi, dist_type, sigma) result(cs)
    ! Angle- and observer-dependent e-i cross section d(sigma)/d(ep) (cm^2 keV^-1).
    real(8), intent(in) :: ee, ep, z, alpha, sigma
    integer, intent(in) :: n_theta, n_phi, dist_type
    integer :: i, j, n_mu
    real(8), allocatable :: theta(:), phi(:), cskew(:), csket(:), thint(:)
    real(8) :: cpsi, val, gwt

    cs = 0d0
    if (n_theta < 2 .or. n_phi < 2 .or. ep <= 0d0 .or. ep >= ee) return

    ! At alpha = 0, cos(psi) = cos(theta): split Gauss-Legendre on mu = cos(theta).
    if (dabs(alpha) <= 1d-12) then
      n_mu = max(n_theta, n_phi)
      cs = int_cs1_mu(ep, ee, z, n_mu, dist_type, sigma)
      return
    end if

    allocate(theta(n_theta), phi(n_phi), cskew(n_phi), csket(n_theta))
    theta(1) = 0d0
    do i = 2, n_theta
      theta(i) = dble(i)*pi/dble(n_theta)
    end do

    phi(1) = 0d0
    do j = 2, n_phi
      phi(j) = dble(j)*2d0*pi/dble(n_phi)
    end do

    do i = 1, n_theta
      gwt = pitch_angle_distribution_raw(theta(i), dist_type, sigma)
      do j = 1, n_phi
        cpsi = cos_psi(alpha, theta(i), phi(j))
        val = cs1(ep, ee, cpsi, z)
        if (val < 0d0 .or. val /= val) val = 0d0
        cskew(j) = val*gwt
      end do
      csket(i) = integrate(phi, cskew)
    end do
    allocate(thint(n_theta))
    do i = 1, n_theta
      thint(i) = dsin(theta(i))*csket(i)
    end do
    cs = integrate(theta, thint)
    cs = cs*me
    deallocate(theta, phi, cskew, csket, thint)
  end function cs1_aniso

  real(8) function eeba_aniso(ee, ep, alpha, n_theta, n_phi, dist_type, sigma) result(cs)
    ! Angle- and observer-dependent e-e cross section d(sigma)/dk (cm^2), matching
    ! ``eebapx`` / ``eebls`` after 2 pi int eeba d mu and multiplication by me.
    real(8), intent(in) :: ee, ep, alpha, sigma
    integer, intent(in) :: n_theta, n_phi, dist_type
    integer :: i, j, n_mu
    real(8), allocatable :: theta(:), phi(:), cskew(:), csket(:), thint(:)
    real(8) :: cpsi, val, gwt, emin

    cs = 0d0
    if (n_theta < 2 .or. n_phi < 2 .or. ep <= 0d0 .or. ep >= ee) return

    if (dabs(alpha) <= 1d-12) then
      n_mu = max(n_theta, n_phi)
      cs = int_eeba_mu(ep, ee, n_mu, dist_type, sigma)
      return
    end if

    allocate(theta(n_theta), phi(n_phi), cskew(n_phi), csket(n_theta))
    theta(1) = 0d0
    do i = 2, n_theta
      theta(i) = dble(i)*pi/dble(n_theta)
    end do

    phi(1) = 0d0
    do j = 2, n_phi
      phi(j) = dble(j)*2d0*pi/dble(n_phi)
    end do

    do i = 1, n_theta
      gwt = pitch_angle_distribution_raw(theta(i), dist_type, sigma)
      do j = 1, n_phi
        cpsi = cos_psi(alpha, theta(i), phi(j))
        emin = eeb_emin_aniso(ep, cpsi)
        if (ee < emin) then
          val = 0d0
        else
          val = eeba(ee, ep, cpsi)
          if (val < 0d0 .or. val /= val) val = 0d0
        end if
        cskew(j) = val*gwt
      end do
      csket(i) = integrate(phi, cskew)
    end do
    allocate(thint(n_theta))
    do i = 1, n_theta
      thint(i) = dsin(theta(i))*csket(i)
    end do
    cs = integrate(theta, thint)*me
    deallocate(theta, phi, cskew, csket, thint)
  end function eeba_aniso

  ! ---------------------------------------------------------------------------
  ! Isotropic angle averages (retained for cross-section validation).
  ! ---------------------------------------------------------------------------

  real(8) function cs1_integrated(ee, ep, z, n_mu) result(cs)
    real(8), intent(in) :: ee, ep, z
    integer, intent(in) :: n_mu

    cs = int_cs1_mu(ep, ee, z, n_mu, adf_isotropic, 0d0)
  end function cs1_integrated

  subroutine cs1_integrated_array(photon_energies, n, ee, z, n_mu, cs_out)
    integer, intent(in) :: n, n_mu
    real(8), intent(in) :: photon_energies(n), ee, z
    real(8), intent(out) :: cs_out(n)
    integer :: i
    do i = 1, n
      cs_out(i) = cs1_integrated(ee, photon_energies(i), z, n_mu)
    end do
  end subroutine cs1_integrated_array

  real(8) function eeba_integrated(ee, ep, n_mu) result(cs)
    ! Angle-averaged e-e cross section d(sigma)/dk (cm^2), analogous to cs1_integrated.
    real(8), intent(in) :: ee, ep
    integer, intent(in) :: n_mu

    cs = 0d0
    if (n_mu < 2 .or. ep <= 0d0 .or. ep >= ee) return
    if (ee < eeb_emin_aniso(ep, 1d0)) return
    cs = int_eeba_mu(ep, ee, n_mu, adf_isotropic, 0d0)
  end function eeba_integrated

  subroutine eeba_integrated_array(photon_energies, n, ee, n_mu, cs_out)
    integer, intent(in) :: n, n_mu
    real(8), intent(in) :: photon_energies(n), ee
    real(8), intent(out) :: cs_out(n)
    integer :: i
    do i = 1, n
      cs_out(i) = eeba_integrated(ee, photon_energies(i), n_mu)
    end do
  end subroutine eeba_integrated_array

  ! ---------------------------------------------------------------------------
  ! Gauss-Legendre (copied from module flux for the electron-energy integral).
  ! ---------------------------------------------------------------------------

  pure subroutine gauleg(nn, x, w)
    integer, intent(in) :: nn
    real(8), intent(out) :: x(nn), w(nn)
    integer :: i, j, m
    real(8) :: z, z1, p1, p2, p3, pp
    real(8), parameter :: tol = 1d-15

    m = (nn + 1)/2
    do i = 1, m
      z = dcos(pi*(dble(i) - 0.25d0)/(dble(nn) + 0.5d0))
      do
        p1 = 1d0
        p2 = 0d0
        do j = 1, nn
          p3 = p2
          p2 = p1
          p1 = ((2d0*dble(j) - 1d0)*z*p2 - (dble(j) - 1d0)*p3)/dble(j)
        end do
        pp = dble(nn)*(z*p1 - p2)/(z*z - 1d0)
        z1 = z
        z = z1 - p1/pp
        if (dabs(z - z1) <= tol) exit
      end do
      x(i) = -z
      x(nn + 1 - i) = z
      w(i) = 2d0/((1d0 - z*z)*pp*pp)
      w(nn + 1 - i) = w(i)
    end do
  end subroutine gauleg

  real(8) function aniso_flux_integrand_ei(ee, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                           alpha, n_theta, n_phi, dist_type, sigma) result(g)
    real(8), intent(in) :: ee, ep, p, eebrk, q, elow, ehigh, z, alpha, sigma
    integer, intent(in) :: efd, n_theta, n_phi, dist_type
    real(8) :: w, cs

    if (efd /= 0) then
      w = me/c
    else
      w = (me/c)*electron_velocity_rel(ee)
    end if
    cs = cs1_aniso(ee, ep, z, alpha, n_theta, n_phi, dist_type, sigma)
    g = electron_flux_density(ee, p, eebrk, q, elow, ehigh)*cs*w
  end function aniso_flux_integrand_ei

  real(8) function aniso_flux_integrand_ee(ee, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                             alpha, n_theta, n_phi, dist_type, sigma) result(g)
    real(8), intent(in) :: ee, ep, p, eebrk, q, elow, ehigh, z, alpha, sigma
    integer, intent(in) :: efd, n_theta, n_phi, dist_type
    real(8) :: w, cs

    if (efd /= 0) then
      w = me/c
    else
      w = (me/c)*electron_velocity_rel(ee)
    end if
    cs = z*eeba_aniso(ee, ep, alpha, n_theta, n_phi, dist_type, sigma)
    g = electron_flux_density(ee, p, eebrk, q, elow, ehigh)*cs*w
  end function aniso_flux_integrand_ee

  real(8) function segment_gl_aniso(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                    emission_kind, alpha, n_theta, n_phi, dist_type, sigma, &
                                    gx, gw, nn) result(s)
    integer, intent(in) :: efd, emission_kind, nn, n_theta, n_phi, dist_type
    real(8), intent(in) :: elo, ehi, ep, p, eebrk, q, elow, ehigh, z, alpha, sigma
    real(8), intent(in) :: gx(nn), gw(nn)
    integer :: j
    real(8) :: a, b, xm, xl, ee_j

    s = 0d0
    if (nn < 1 .or. ehi <= elo) return

    a = dlog(elo)
    b = dlog(ehi)
    xm = 0.5d0*(b + a)
    xl = 0.5d0*(b - a)
    do j = 1, nn
      ee_j = dexp(xm + xl*gx(j))
      if (emission_kind == 1) then
        s = s + gw(j)*aniso_flux_integrand_ee(ee_j, ep, p, eebrk, q, elow, ehigh, z, &
                                              efd, alpha, n_theta, n_phi, dist_type, sigma)*ee_j
      else
        s = s + gw(j)*aniso_flux_integrand_ei(ee_j, ep, p, eebrk, q, elow, ehigh, z, &
                                              efd, alpha, n_theta, n_phi, dist_type, sigma)*ee_j
      end if
    end do
    s = s*xl
  end function segment_gl_aniso

  subroutine bremsstrahlung_thin_target_aniso(photon_energies, n, p, break_energy, q, &
                                             low_e_cutoff, high_e_cutoff, z, efd, npts, &
                                             alpha, n_theta, n_phi, dist_type, sigma, &
                                             distance_cm, flux_out)
    ! Anisotropic thin-target e-i bremsstrahlung photon flux.
    integer, intent(in) :: n, efd, npts, n_theta, n_phi, dist_type
    real(8), intent(in) :: photon_energies(n)
    real(8), intent(in) :: p, break_energy, q, low_e_cutoff, high_e_cutoff, z
    real(8), intent(in) :: alpha, sigma, distance_cm
    real(8), intent(out) :: flux_out(n)
    integer :: i
    real(8) :: ep, lo, integral, fcoeff, dist
    real(8) :: gx(npts), gw(npts)

    if (distance_cm > 0d0) then
      dist = distance_cm
    else
      dist = au_default
    end if
    fcoeff = c/(4d0*pi*dist*dist*me**2)
    call gauleg(npts, gx, gw)

    do i = 1, n
      ep = photon_energies(i)
      flux_out(i) = 0d0
      if (ep <= 0d0 .or. ep >= high_e_cutoff) cycle

      integral = 0d0
      if (break_energy > low_e_cutoff) then
        lo = max(ep, low_e_cutoff)
        if (break_energy > lo) integral = integral + &
          segment_gl_aniso(lo, break_energy, ep, p, break_energy, q, low_e_cutoff, &
                           high_e_cutoff, z, efd, 0, alpha, n_theta, n_phi, dist_type, &
                           sigma, gx, gw, npts)
      end if
      if (high_e_cutoff > break_energy) then
        lo = max(ep, break_energy)
        if (high_e_cutoff > lo) integral = integral + &
          segment_gl_aniso(lo, high_e_cutoff, ep, p, break_energy, q, low_e_cutoff, &
                           high_e_cutoff, z, efd, 0, alpha, n_theta, n_phi, dist_type, &
                           sigma, gx, gw, npts)
      end if
      flux_out(i) = fcoeff*integral
    end do
  end subroutine bremsstrahlung_thin_target_aniso

  subroutine bremsstrahlung_thin_target_eeb_aniso(photon_energies, n, p, break_energy, q, &
                                                  low_e_cutoff, high_e_cutoff, z, efd, npts, &
                                                  alpha, n_theta, n_phi, dist_type, sigma, &
                                                  distance_cm, flux_out)
    ! Anisotropic thin-target e-e bremsstrahlung photon flux.
    ! Per-direction lower electron-energy limits are applied inside eeba_aniso.
    integer, intent(in) :: n, efd, npts, n_theta, n_phi, dist_type
    real(8), intent(in) :: photon_energies(n)
    real(8), intent(in) :: p, break_energy, q, low_e_cutoff, high_e_cutoff, z
    real(8), intent(in) :: alpha, sigma, distance_cm
    real(8), intent(out) :: flux_out(n)
    integer :: i
    real(8) :: ep, lo, emin, integral, fcoeff, dist
    real(8) :: gx(npts), gw(npts)

    if (distance_cm > 0d0) then
      dist = distance_cm
    else
      dist = au_default
    end if
    fcoeff = c/(4d0*pi*dist*dist*me**2)
    call gauleg(npts, gx, gw)

    do i = 1, n
      ep = photon_energies(i)
      flux_out(i) = 0d0
      if (ep <= 0d0 .or. ep >= high_e_cutoff) cycle

      emin = eeb_emin_iso(ep)
      integral = 0d0
      if (break_energy > low_e_cutoff) then
        lo = max(emin, low_e_cutoff)
        if (break_energy > lo) integral = integral + &
          segment_gl_aniso(lo, break_energy, ep, p, break_energy, q, low_e_cutoff, &
                           high_e_cutoff, z, efd, 1, alpha, n_theta, n_phi, dist_type, &
                           sigma, gx, gw, npts)
      end if
      if (high_e_cutoff > break_energy) then
        lo = max(emin, break_energy)
        if (high_e_cutoff > lo) integral = integral + &
          segment_gl_aniso(lo, high_e_cutoff, ep, p, break_energy, q, low_e_cutoff, &
                           high_e_cutoff, z, efd, 1, alpha, n_theta, n_phi, dist_type, &
                           sigma, gx, gw, npts)
      end if
      flux_out(i) = fcoeff*integral
    end do
  end subroutine bremsstrahlung_thin_target_eeb_aniso

  real(8) function aniso_flux_integrand_ee_int(ee, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                                 n_mu) result(g)
    ! Isotropic e-e integrand using angle-averaged eeba_integrated (validation reference).
    real(8), intent(in) :: ee, ep, p, eebrk, q, elow, ehigh, z
    integer, intent(in) :: efd, n_mu
    real(8) :: w, cs

    if (efd /= 0) then
      w = me/c
    else
      w = (me/c)*electron_velocity_rel(ee)
    end if
    cs = z*eeba_integrated(ee, ep, n_mu)
    g = electron_flux_density(ee, p, eebrk, q, elow, ehigh)*cs*w
  end function aniso_flux_integrand_ee_int

  real(8) function aniso_flux_integrand_ei_int(ee, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                                 n_mu) result(g)
    ! Isotropic e-i integrand using cs1_integrated (validation reference).
    real(8), intent(in) :: ee, ep, p, eebrk, q, elow, ehigh, z
    integer, intent(in) :: efd, n_mu
    real(8) :: w, cs

    if (efd /= 0) then
      w = me/c
    else
      w = (me/c)*electron_velocity_rel(ee)
    end if
    cs = cs1_integrated(ee, ep, z, n_mu)
    g = electron_flux_density(ee, p, eebrk, q, elow, ehigh)*cs*w
  end function aniso_flux_integrand_ei_int

  real(8) function segment_gl_ei_int(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, n_mu, &
                                     gx, gw, nn) result(s)
    integer, intent(in) :: efd, n_mu, nn
    real(8), intent(in) :: elo, ehi, ep, p, eebrk, q, elow, ehigh, z
    real(8), intent(in) :: gx(nn), gw(nn)
    integer :: j
    real(8) :: a, b, xm, xl, ee_j

    s = 0d0
    if (nn < 1 .or. ehi <= elo) return

    a = dlog(elo)
    b = dlog(ehi)
    xm = 0.5d0*(b + a)
    xl = 0.5d0*(b - a)
    do j = 1, nn
      ee_j = dexp(xm + xl*gx(j))
      s = s + gw(j)*aniso_flux_integrand_ei_int(ee_j, ep, p, eebrk, q, elow, ehigh, z, &
                                                  efd, n_mu)*ee_j
    end do
    s = s*xl
  end function segment_gl_ei_int

  subroutine bremsstrahlung_thin_target_ei_numeric(photon_energies, n, p, break_energy, q, &
                                                   low_e_cutoff, high_e_cutoff, z, efd, npts, &
                                                   n_mu, distance_cm, flux_out)
    ! Isotropic e-i thin-target flux using cs1_integrated (2 pi int cs1 d mu).
    integer, intent(in) :: n, efd, npts, n_mu
    real(8), intent(in) :: photon_energies(n)
    real(8), intent(in) :: p, break_energy, q, low_e_cutoff, high_e_cutoff, z, distance_cm
    real(8), intent(out) :: flux_out(n)
    integer :: i
    real(8) :: ep, lo, integral, fcoeff, dist
    real(8) :: gx(npts), gw(npts)

    if (distance_cm > 0d0) then
      dist = distance_cm
    else
      dist = au_default
    end if
    fcoeff = c/(4d0*pi*dist*dist*me**2)
    call gauleg(npts, gx, gw)

    do i = 1, n
      ep = photon_energies(i)
      flux_out(i) = 0d0
      if (ep <= 0d0 .or. ep >= high_e_cutoff) cycle

      integral = 0d0
      if (break_energy > low_e_cutoff) then
        lo = max(ep, low_e_cutoff)
        if (break_energy > lo) integral = integral + &
          segment_gl_ei_int(lo, break_energy, ep, p, break_energy, q, low_e_cutoff, &
                            high_e_cutoff, z, efd, n_mu, gx, gw, npts)
      end if
      if (high_e_cutoff > break_energy) then
        lo = max(ep, break_energy)
        if (high_e_cutoff > lo) integral = integral + &
          segment_gl_ei_int(lo, high_e_cutoff, ep, p, break_energy, q, low_e_cutoff, &
                            high_e_cutoff, z, efd, n_mu, gx, gw, npts)
      end if
      flux_out(i) = fcoeff*integral
    end do
  end subroutine bremsstrahlung_thin_target_ei_numeric

  real(8) function segment_gl_ee_int(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, n_mu, &
                                     gx, gw, nn) result(s)
    integer, intent(in) :: efd, n_mu, nn
    real(8), intent(in) :: elo, ehi, ep, p, eebrk, q, elow, ehigh, z
    real(8), intent(in) :: gx(nn), gw(nn)
    integer :: j
    real(8) :: a, b, xm, xl, ee_j

    s = 0d0
    if (nn < 1 .or. ehi <= elo) return

    a = dlog(elo)
    b = dlog(ehi)
    xm = 0.5d0*(b + a)
    xl = 0.5d0*(b - a)
    do j = 1, nn
      ee_j = dexp(xm + xl*gx(j))
      s = s + gw(j)*aniso_flux_integrand_ee_int(ee_j, ep, p, eebrk, q, elow, ehigh, z, &
                                                  efd, n_mu)*ee_j
    end do
    s = s*xl
  end function segment_gl_ee_int

  subroutine bremsstrahlung_thin_target_eeb_numeric(photon_energies, n, p, break_energy, q, &
                                                   low_e_cutoff, high_e_cutoff, z, efd, npts, &
                                                   n_mu, distance_cm, flux_out)
    ! Isotropic e-e thin-target flux using eeba_integrated (2 pi int eeba d mu).
    ! Matches bremsstrahlung_thin_target_eeb_aniso with ADF_ISOTROPIC at alpha = 0;
    ! use this instead of eebls when validating the theta-phi integration.
    integer, intent(in) :: n, efd, npts, n_mu
    real(8), intent(in) :: photon_energies(n)
    real(8), intent(in) :: p, break_energy, q, low_e_cutoff, high_e_cutoff, z, distance_cm
    real(8), intent(out) :: flux_out(n)
    integer :: i
    real(8) :: ep, lo, emin, integral, fcoeff, dist
    real(8) :: gx(npts), gw(npts)

    if (distance_cm > 0d0) then
      dist = distance_cm
    else
      dist = au_default
    end if
    fcoeff = c/(4d0*pi*dist*dist*me**2)
    call gauleg(npts, gx, gw)

    do i = 1, n
      ep = photon_energies(i)
      flux_out(i) = 0d0
      if (ep <= 0d0 .or. ep >= high_e_cutoff) cycle

      emin = eeb_emin_iso(ep)
      integral = 0d0
      if (break_energy > low_e_cutoff) then
        lo = max(emin, low_e_cutoff)
        if (break_energy > lo) integral = integral + &
          segment_gl_ee_int(lo, break_energy, ep, p, break_energy, q, low_e_cutoff, &
                            high_e_cutoff, z, efd, n_mu, gx, gw, npts)
      end if
      if (high_e_cutoff > break_energy) then
        lo = max(emin, break_energy)
        if (high_e_cutoff > lo) integral = integral + &
          segment_gl_ee_int(lo, high_e_cutoff, ep, p, break_energy, q, low_e_cutoff, &
                            high_e_cutoff, z, efd, n_mu, gx, gw, npts)
      end if
      flux_out(i) = fcoeff*integral
    end do
  end subroutine bremsstrahlung_thin_target_eeb_numeric

end module flux_aniso
