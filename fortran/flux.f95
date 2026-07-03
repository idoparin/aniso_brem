! flux.f95 -- compute photon flux from electron-ion (e-i) and electron-electron
! (e-e) bremsstrahlung emitted by a broken power-law distribution of electrons.
!
! This is a fast, compiled drop-in replacement for the isotropic thin-target
! model in sunkit-spex (sunkit_spex/models/physical/nonthermal.py). It reuses
! the analytic cross sections from the `cseqs` module and performs the
! electron-energy integral over each smooth power-law segment (split at the
! distribution break energy, exactly as sunkit-spex splits its integral) with a
! selectable quadrature rule:
!   integrator = 0 : trapezoidal rule on a log-spaced grid (cseqs `integrate`)
!   integrator = 1 : fixed-order Gauss-Legendre in log-energy (like sunkit-spex)
!
! Two emission processes share the same integration machinery (selected by the
! `kind` flag of the integrand):
!   kind = 0 : e-i bremsstrahlung, cross section cspax  (bremsstrahlung_thin_target)
!   kind = 1 : e-e bremsstrahlung, cross section eebapx/eebls with the higher
!              e-e kinematic lower energy limit  (bremsstrahlung_thin_target_eeb)


module flux

  use cseqs, only: cspax, eebapx, eebls, integrate, me, c

  implicit none
  private
  public :: electron_flux_density, electron_flux_density_array, &
            electron_velocity_rel, eeb_emin, eeb_cross, &
            bremsstrahlung_thin_target, bremsstrahlung_thin_target_eeb

  real(8), parameter :: pi_loc = 3.141592653589793d0
  ! Fixed photon-energy switch (keV) for isotropic e-e: eebapx below, eebls at/above.
  ! eebapx diverges above ~150 keV (Haug 1998 eq. 2); do not use it at higher photon energy.
  real(8), parameter :: eeb_photon_switch_keV = 150d0
  ! Default observer distance: 1 AU in cm (SSW IDL / sunxspex constant).
  real(8), parameter :: au_default = 1.496d13

contains

  pure real(8) function electron_flux_density(ee, p, eebrk, q, elow, ehigh) result(f)
    ! Normalized broken (double) power-law electron flux density, matching
    ! BrokenPowerLawElectronDistribution.flux in sunkit-spex.
    ! The distribution integrates to unity over [elow, ehigh].
    !   ee    : electron kinetic energy (keV)
    !   p     : spectral index below the break
    !   eebrk : break energy (keV)
    !   q     : spectral index above the break
    !   elow  : low-energy cutoff (keV)
    !   ehigh : high-energy cutoff (keV)
    real(8), intent(in) :: ee, p, eebrk, q, elow, ehigh
    real(8) :: n0, n1, n2, norm

    n0 = (q - 1d0)/(p - 1d0) * eebrk**(p - 1d0) * elow**(1d0 - p)
    n1 = n0 - (q - 1d0)/(p - 1d0)
    n2 = 1d0 - eebrk**(q - 1d0) * ehigh**(1d0 - q)
    norm = 1d0/(n1 + n2)

    if (ee < elow) then
      f = 0d0
    else if (ee < eebrk) then
      f = norm * n0 * (p - 1d0) * ee**(-p) * elow**(p - 1d0)
    else if (ee <= ehigh) then
      f = norm * (q - 1d0) * ee**(-q) * eebrk**(q - 1d0)
    else
      f = 0d0
    end if
  end function electron_flux_density

  pure subroutine electron_flux_density_array(ee, n, p, eebrk, q, elow, ehigh, f)
    ! Vectorized electron_flux_density: evaluate the distribution at every energy
    ! in ee(:) with the loop running in compiled Fortran (so a whole numpy array
    ! can be passed from Python in a single call instead of element-by-element).
    integer, intent(in) :: n
    real(8), intent(in) :: ee(n), p, eebrk, q, elow, ehigh
    real(8), intent(out) :: f(n)
    integer :: i
    do i = 1, n
      f(i) = electron_flux_density(ee(i), p, eebrk, q, elow, ehigh)
    end do
  end subroutine electron_flux_density_array

  pure real(8) function electron_velocity_rel(ee) result(v)
    ! Relativistic electron speed (cm/s) from kinetic energy ee (keV):
    !   v = c * sqrt(1 - 1/gamma^2),  gamma = 1 + ee/mc2.
    real(8), intent(in) :: ee
    v = dsqrt(1d0 - (1d0/(1d0 + ee/me))**2)*c
  end function electron_velocity_rel

  pure real(8) function eeb_emin(ep) result(emin)
    ! Minimum electron kinetic energy (keV) that can emit a photon of energy ep
    ! in isotropic e-e bremsstrahlung -- Haug (1975), Solar Phys. 45, 453,
    ! eq. (2.6); https://ui.adsabs.harvard.edu/abs/1975SoPh...45..453H .
    !   k = ep/mc2
    !   emin = ep * (2 + 3k - sqrt(k^2 + 4k)) / (1 - k^2 + k*sqrt(k^2 + 4k))
    ! The factor multiplying ep decreases from 2 (k->0) to 1 (k->inf), so the
    ! e-e threshold always exceeds the e-i threshold (which is just ep).
    real(8), intent(in) :: ep
    real(8) :: k, s
    k = ep/me
    s = dsqrt(k*k + 4d0*k)
    emin = ep*(2d0 + 3d0*k - s)/(1d0 - k*k + k*s)
  end function eeb_emin

  real(8) function eeb_cross(ee, ep, z, ep_switch, xparts) result(cs)
    ! Isotropic e-e bremsstrahlung d(sigma)/dk (cm^2), scaled by z, for the thin-target
    ! flux integrand (same d(sigma)/dk convention as ``cspax`` / ``cs1`` after angle
    ! integration).  ``eebapx`` and ``eebls`` return d(sigma)/dk with k = ep/me.
    !
    ! Sharp switch at 150 keV (``eeb_photon_switch_keV``):
    !   ep < 150 keV  -> eebapx only
    !   ep >= 150 keV -> eebls only  (eebapx is never called at/above the switch)
    ! The ``ep_switch`` argument is retained for API compatibility but is not used.
    ! NaN / negative values are clamped to 0.
    real(8), intent(in) :: ee, ep, z, ep_switch
    integer, intent(in) :: xparts
    if (ep >= ee) then
      cs = 0d0
      return
    end if
    if (ep < eeb_photon_switch_keV) then
      cs = eebapx(ee, ep)
    else
      cs = eebls(ee, ep, xparts)
    end if
    if (cs /= cs .or. cs < 0d0) cs = 0d0
    cs = z*cs
  end function eeb_cross

  real(8) function flux_integrand(ee, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                  kind, ep_switch, xparts) result(g)
    ! The dE integrand of the thin-target flux:
    !   spectrum(E) * dsigma/dk(E,ep) * w(E)
    ! where spectrum(E) is the broken power law (electron_flux_density) and the
    ! cross section is selected by `kind`:
    !   kind = 0 : e-i, cspax(E,ep,z)
    !   kind = 1 : e-e, eeb_cross(E,ep,z,ep_switch,xparts)
    ! w(E) converts the distribution to the photon-flux integrand:
    !   efd /= 0 : input is the electron FLUX density -> w = mc2/clight
    !   efd == 0 : input is the electron NUMBER density -> multiply by the
    !              relativistic velocity to turn density into flux:
    !              w = (mc2/clight) * v(E)  ==  pc/gamma  (sunxspex fork).
    ! Single source of truth shared by both quadrature rules and both processes.
    real(8), intent(in) :: ee, ep, p, eebrk, q, elow, ehigh, z, ep_switch
    integer, intent(in) :: efd, kind, xparts
    real(8) :: w, cs
    if (efd /= 0) then
      w = me/c
    else
      w = (me/c)*electron_velocity_rel(ee)
    end if
    if (kind == 1) then
      cs = eeb_cross(ee, ep, z, ep_switch, xparts)
    else
      cs = cspax(ee, ep, z)
    end if
    g = electron_flux_density(ee, p, eebrk, q, elow, ehigh)*cs*w
  end function flux_integrand

  pure subroutine gauleg(nn, x, w)
    ! Gauss-Legendre nodes x(:) and weights w(:) on [-1, 1].
    !
    ! Algorithm: the nodes are the roots of the Legendre polynomial P_nn, found
    ! by Newton-Raphson; P_nn and its derivative P_nn' are evaluated with the
    ! standard three-term recurrence, the weights from w_i = 2/((1-x_i^2)P_nn'^2).
    ! The cosine starting guess and recurrence are the classic textbook routine.
    !
    ! References:
    !   * Press, Teukolsky, Vetterling & Flannery, "Numerical Recipes in
    !     Fortran 77", 2nd ed. (1992), section 4.5, routine `gauleg`:
    !     https://numerical.recipes/  (see also NR in C, section 4.5).
    !   * Abramowitz & Stegun, "Handbook of Mathematical Functions" (1972),
    !     sections 25.4.29-25.4.30 (Gauss-Legendre nodes/weights).
    !   * Same scheme as SSW IDL `Brm_GauLeg54.pro` used by sunkit-spex/sunxspex:
    !     https://hesperia.gsfc.nasa.gov/ssw/packages/xray/idl/brm/brm_gauleg54.pro
    !   * Produces the same nodes/weights as `scipy.special.roots_legendre`
    !     (numpy.polynomial.legendre.leggauss), which sunkit-spex calls.
    !
    ! The nodes are independent of the integration limits, so this is computed
    ! once per flux evaluation and reused for every photon energy and segment.
    integer, intent(in) :: nn
    real(8), intent(out) :: x(nn), w(nn)
    integer :: i, j, m
    real(8) :: z, z1, p1, p2, p3, pp
    real(8), parameter :: tol = 1d-15

    m = (nn + 1)/2
    do i = 1, m
      z = dcos(pi_loc*(dble(i) - 0.25d0)/(dble(nn) + 0.5d0))
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

  real(8) function segment_trap(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                kind, ep_switch, xparts, npts) result(s)
    ! Trapezoidal integral of flux_integrand over [elo, ehi], sampled on a
    ! log-spaced electron-energy grid and summed with the shared cseqs
    ! `integrate` routine (the same trapezoid used by the cross sections).
    real(8), intent(in) :: elo, ehi, ep, p, eebrk, q, elow, ehigh, z, ep_switch
    integer, intent(in) :: efd, kind, xparts, npts
    integer :: j
    real(8) :: lx0, lx1, dlx
    real(8), dimension(npts) :: ee, yint

    s = 0d0
    if (npts < 2 .or. ehi <= elo) return

    lx0 = dlog(elo)
    lx1 = dlog(ehi)
    dlx = (lx1 - lx0)/dble(npts - 1)
    do j = 1, npts
      ee(j) = dexp(lx0 + dble(j - 1)*dlx)
      yint(j) = flux_integrand(ee(j), ep, p, eebrk, q, elow, ehigh, z, efd, &
                               kind, ep_switch, xparts)
    end do
    s = integrate(ee, yint)
  end function segment_trap

  real(8) function segment_gl(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, &
                              kind, ep_switch, xparts, gx, gw, nn) result(s)
    ! Fixed-order Gauss-Legendre integral of flux_integrand over [elo, ehi],
    ! carried out in log-energy u = ln(E) so a power-law integrand is a smooth
    ! exponential (this is why GL converges with very few nodes here). The
    ! [-1,1] nodes/weights gx, gw are precomputed once by gauleg.
    !   integral_E f dE = integral_u f(E(u)) * E(u) du   (Jacobian dE = E du)
    real(8), intent(in) :: elo, ehi, ep, p, eebrk, q, elow, ehigh, z, ep_switch
    integer, intent(in) :: efd, kind, xparts, nn
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
      s = s + gw(j)*flux_integrand(ee_j, ep, p, eebrk, q, elow, ehigh, z, efd, &
                                   kind, ep_switch, xparts)*ee_j
    end do
    s = s*xl
  end function segment_gl

  subroutine bremsstrahlung_thin_target(photon_energies, n, p, break_energy, q, &
                                        low_e_cutoff, high_e_cutoff, z, efd, &
                                        integrator, npts, distance_cm, flux_out)
    ! Thin-target e-i bremsstrahlung photon flux for a broken power-law
    ! electron distribution. Mirrors sunxspex bremsstrahlung_thin_target,
    ! returning photons s^-1 keV^-1 cm^-2 at the observer per unit
    ! (nth * V * nnth) coefficient, including geometric spreading 1/(4 pi R^2).
    !
    !   distance_cm     : source-observer distance (cm); <=0 uses 1 AU default
    !
    !   photon_energies : photon energies to evaluate at (keV)
    !   n               : number of photon energies (f2py-optional)
    !   p, break_energy, q, low_e_cutoff, high_e_cutoff : distribution params
    !   z               : mean ion charge (sunkit-spex uses 1.2)
    !   efd             : /=0 electron flux density input, ==0 density input
    !   integrator      : 0 = trapezoid, 1 = Gauss-Legendre
    !   npts            : quadrature points per energy segment
    !                     (trapezoid: grid points; GL: number of nodes)
    !   flux_out        : output flux array (n)
    integer, intent(in) :: n
    real(8), intent(in) :: photon_energies(n)
    real(8), intent(in) :: p, break_energy, q, low_e_cutoff, high_e_cutoff, z
    integer, intent(in) :: efd, integrator, npts
    real(8), intent(in) :: distance_cm
    real(8), intent(out) :: flux_out(n)

    integer :: i
    real(8) :: ep, lo, integral, fcoeff, dist
    real(8) :: gx(npts), gw(npts)

    ! sunxspex: fcoeff = clight / (4*pi*au^2) / mc2^2 ; with the (mc2/clight)
    ! weight in the efd integrand this is c / (4*pi*R^2*mc2^2).
    if (distance_cm > 0d0) then
      dist = distance_cm
    else
      dist = au_default
    end if
    fcoeff = c/(4d0*pi_loc*dist*dist*me**2)

    ! Gauss-Legendre nodes are limit-independent: build them once and reuse.
    if (integrator == 1) call gauleg(npts, gx, gw)

    do i = 1, n
      ep = photon_energies(i)
      flux_out(i) = 0d0
      if (ep <= 0d0 .or. ep >= high_e_cutoff) cycle

      integral = 0d0

      ! Part below the break: [max(ep, low_e_cutoff), break_energy]
      if (break_energy > low_e_cutoff) then
        lo = max(ep, low_e_cutoff)
        if (break_energy > lo) integral = integral + &
          segment(lo, break_energy, ep, p, break_energy, q, low_e_cutoff, &
                  high_e_cutoff, z, efd, 0, 0d0, 1, integrator, npts, gx, gw)
      end if

      ! Part above the break: [max(ep, break_energy), high_e_cutoff]
      if (high_e_cutoff > break_energy) then
        lo = max(ep, break_energy)
        if (high_e_cutoff > lo) integral = integral + &
          segment(lo, high_e_cutoff, ep, p, break_energy, q, low_e_cutoff, &
                  high_e_cutoff, z, efd, 0, 0d0, 1, integrator, npts, gx, gw)
      end if

      flux_out(i) = fcoeff*integral
    end do
  end subroutine bremsstrahlung_thin_target

  subroutine bremsstrahlung_thin_target_eeb(photon_energies, n, p, break_energy, q, &
                                            low_e_cutoff, high_e_cutoff, z, efd, &
                                            integrator, npts, ep_switch, xparts, &
                                            distance_cm, flux_out)
    ! Thin-target electron-electron (e-e) bremsstrahlung photon flux for a broken
    ! power-law electron distribution. Same integration machinery, units and
    ! normalization coefficient as the e-i routine bremsstrahlung_thin_target,
    ! but with:
    !   * the e-e differential cross section eeb_cross (eebapx below ep_switch,
    !     eebls above) scaled by z, and
    !   * the higher e-e kinematic lower energy limit eeb_emin(ep), Haug (1975),
    !     applied to both spectral segments (the cross section vanishes below it).
    !
    !   photon_energies : photon energies to evaluate at (keV)
    !   n               : number of photon energies (f2py-optional)
    !   p, break_energy, q, low_e_cutoff, high_e_cutoff : distribution params
    !   z               : mean ion charge (linear weight for the e-e component)
    !   efd             : /=0 electron flux density input, ==0 density input
    !   integrator      : 0 = trapezoid, 1 = Gauss-Legendre
    !   npts            : quadrature points per energy segment
    !   ep_switch       : ignored; switch fixed at 150 keV (eeb_photon_switch_keV)
    !   xparts          : internal integration sub-steps for eebls
    !   distance_cm     : source-observer distance (cm); <=0 uses 1 AU default
    !   flux_out        : output flux array (n)
    integer, intent(in) :: n
    real(8), intent(in) :: photon_energies(n)
    real(8), intent(in) :: p, break_energy, q, low_e_cutoff, high_e_cutoff, z
    real(8), intent(in) :: ep_switch, distance_cm
    integer, intent(in) :: efd, integrator, npts, xparts
    real(8), intent(out) :: flux_out(n)

    integer :: i
    real(8) :: ep, lo, emin, integral, fcoeff, dist
    real(8) :: gx(npts), gw(npts)

    if (distance_cm > 0d0) then
      dist = distance_cm
    else
      dist = au_default
    end if
    fcoeff = c/(4d0*pi_loc*dist*dist*me**2)
    if (integrator == 1) call gauleg(npts, gx, gw)

    do i = 1, n
      ep = photon_energies(i)
      flux_out(i) = 0d0
      if (ep <= 0d0 .or. ep >= high_e_cutoff) cycle

      emin = eeb_emin(ep)
      integral = 0d0

      ! Part below the break: [max(eeb_emin(ep), low_e_cutoff), break_energy]
      if (break_energy > low_e_cutoff) then
        lo = max(emin, low_e_cutoff)
        if (break_energy > lo) integral = integral + &
          segment(lo, break_energy, ep, p, break_energy, q, low_e_cutoff, &
                  high_e_cutoff, z, efd, 1, ep_switch, xparts, integrator, npts, gx, gw)
      end if

      ! Part above the break: [max(eeb_emin(ep), break_energy), high_e_cutoff]
      if (high_e_cutoff > break_energy) then
        lo = max(emin, break_energy)
        if (high_e_cutoff > lo) integral = integral + &
          segment(lo, high_e_cutoff, ep, p, break_energy, q, low_e_cutoff, &
                  high_e_cutoff, z, efd, 1, ep_switch, xparts, integrator, npts, gx, gw)
      end if

      flux_out(i) = fcoeff*integral
    end do
  end subroutine bremsstrahlung_thin_target_eeb

  real(8) function segment(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, &
                           kind, ep_switch, xparts, integrator, npts, gx, gw) result(s)
    ! Dispatch one segment integral to the requested quadrature rule.
    real(8), intent(in) :: elo, ehi, ep, p, eebrk, q, elow, ehigh, z, ep_switch
    integer, intent(in) :: efd, kind, xparts, integrator, npts
    real(8), intent(in) :: gx(npts), gw(npts)
    if (integrator == 1) then
      s = segment_gl(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, &
                     kind, ep_switch, xparts, gx, gw, npts)
    else
      s = segment_trap(elo, ehi, ep, p, eebrk, q, elow, ehigh, z, efd, &
                       kind, ep_switch, xparts, npts)
    end if
  end function segment

end module flux
