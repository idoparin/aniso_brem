! cseqs.f95 -- compute fully relativistic isotropic (dsigma/dk) and 
! anisotropic cross-sections (dsigma/(dk domega)) of electron-ion 
! and electron-electron bremsstrahlung.
! by I. Oparin <ido4@njit.edu> 2019

module cseqs

  implicit none

  private

  public :: eebls,eeblsnr,eeblser,eebapx,eeba,cspe,cspax,cs1,integrate
  public pi,r0,c,ech,me,afs

  real(8), parameter :: pi=3.141592653589793d0, r0=2.817940322719e-13,&
     c=2.99792458e10, ech=4.8e-10, me=510.9989461d0, afs=0.007297352566417d0

     
  real(8) :: dilog
  external dilog ! compile with dilog490.f Comm. ACM 18 200 (c)
  
contains

  real(8) function eebls (ee, ep, xparts)
  
    ! Haug 1998 EEB cross-section in lab. system
    ! https://doi.org/10.1023/A:1005098624121
    ! https://ui.adsabs.harvard.edu/abs/1998SoPh..178..341H/abstractZ
    ! Calculate differential electron-electron bremsstrahlung cross-section in laboratory system
    ! All values expressed in units of electron rest energy mc**2 and momenta in units of mc

    implicit none

    integer, intent(in) :: xparts
    real(8), intent(in) :: ee, ep !input electron and photon energies (keV)
    integer i
    real(8) :: eps, k, k0, p !eps is electron total energy, k is photon energy
    real(8), dimension(xparts) :: x, chiint !dummy variable
    real(8) :: x1, x2 !upper and lower integration limits
    real(8) :: W1, W4

    common W1, W4, p, eps, k, k0

    eps=1d0+(ee/me)
    k=ep/me
    p=dsqrt((eps**2)-1d0)

    k0=(eps-1d0)/(eps+p+1d0)
    ! print *, 'k', k
    x1=k*(eps-p)
    if(k.lt.k0) then
      x2=k*(eps+p)
    else
      x2=eps-k-1d0
    end if

    !******

    W1=dsqrt((p**2)-(2d0*eps*k)+(k**2))
    W4=dsqrt(((eps-k+1)**2)-4d0)

    !******

    forall(i=1:xparts) x(i)=(i-1)*((x2-x1)/(xparts))+x1

    eebls=sigma(x2,2)-sigma(x1,1)
    do i=1, xparts !PERFORM INTEGRATION
      chiint(i)=chi(x(i))
    end do
    eebls=eebls+integrate(x,chiint)
    eebls=eebls*((afs*(r0**2))/(k*(p**2)))

    return
  end function eebls

  real(8) function chi(x)
    use ieee_arithmetic, only: ieee_is_nan
    implicit none
    real(8), intent(in) :: x
    real(8) eps, p, k
    real(8) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, &
      b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27
    real(8) :: W, W1, W4, L1, L2, L3, L4, L, rho, Lh, Wh, L5, L5h, &
        W5, W5h, L6, W6
    common W1, W4, p, eps, k

    !******

    W=dsqrt(((eps-k-x)**2)-1d0)
    L1=dlog(eps-k-x+W)
    L2=dlog(((((p**2)-(2d0*eps*k)+(k**2))+(W1*W))/x)-eps+k)
    L3=dlog((((p**2)+(p*W))/(x+k))-eps)
    L4=dlog(1d0+(((eps-k+1d0)/(((p**2)-(2d0*eps*k)+(k**2))+(2d0*x))) &
      *(((p**2)-(2d0*eps*k)+(k**2))-(x*(eps-k+1d0))+(W4*W))))
    rho=dsqrt(2d0*(eps+1d0-k-x))
    L=dlog((1d0/(2d0*x))*((rho*(eps-k-1d0)) &
      +dsqrt((((p**2)-(2d0*eps*k)+(k**2))+(2d0*x))*((rho**2)-4d0))))
    Wh=dsqrt(((eps-x)**2)+(2d0*k)-1d0)
    Lh=dlog((1d0/(2d0*k))*((rho*(eps-1d0-x)) &
      +(Wh*dsqrt((rho**2)-4d0))))
    W5=dsqrt(((k**2)*(W**2))+(2d0*k*x*(x+k)))
    W5h=dsqrt(((x**2)*(W**2))+(2d0*k*x*(x+k)))
    L5h=dlog(1d0+(((x*W)+W5h)*W/(k*(x+k))))
    L5=dlog(1d0+(((k*W)+W5)*W/(x*(x+k))))
    W6=dsqrt((((eps-1d0)**2)*(W**2))+(2d0*k*x*(eps-1d0)))
    L6=dlog(1d0+((((eps-1d0)*W)+W6)*W/(k*x)))

    chi=0

    !******

    b1=(k*(5d0+(4d0*k)-(9d0*eps)+(((2d0*k)+1.2d1)/(p**2))+(6d0/(p**4)))) &
    +(2.d0*eps*((4d0*eps)-3d0-(2d0/(eps-1d0))))-(8d0*(eps**3)/k) &
    +((eps*((2d0*eps)-(k**2)))/(eps-k-1d0))

    b2=(k*((2d0*(eps+1))+(4d0*((eps**3)-(3d0*(eps**2))+1d0)/(((eps-1d0)**2) &
    -(k**2)))))+((8d0*(eps**3))/k)

    b3=(eps*k)-(p**2)+(3d0*eps)-k+1d0-((2d0*eps)/k)+((2d0-(6d0*k) &
    -(2d0*(eps-1d0)*(k**2)))/(eps+k-1d0))

    b4=(2d0*k)+(5d-1*(9d0-eps))-(3d0*k/(eps-k-1d0))

    b5=(1.5d0*(p**2))-(k*((4d0*eps)-k+8d0))+((k**2)/(eps-k-1d0))

    b6=(k*((4d0*(eps-k))-2d0))/(eps-k-1d0)

    b7=5d-1*((p**2)-(k*(eps-1)))*((p**2)-(k*((3d0*eps)-(2d0*k)+1d0)))

    b8=5d-1*(eps+k-1d0)*((p**2)-(k*((3d0*eps)-(2d0*k)+1d0)))

    b9=(2d0*k)-(4d0*eps)-9d0+(4d0/k)

    b10=(1.5d0*(p**2))+(5d-1*(9d0-eps)*k)

    b11=eps-1d0-(3d0*k)+(2d0/k)

    b12=eps*((eps*k)-(2d0*(p**2))-k)

    b13=5d-1*((5d0*(eps**2))+(eps*k)-(2d0*eps)+(3d0*k)-3d0)

    b14=eps+k-1d0

    b15=5d-1*(p**2)*((p**2)-(k*(eps-1)))

    b16=(6d0*eps*k)+(4d0*eps)-(2d0*k)+(k*((2d0*(k**2))-(eps*k) &
    +(2d0*eps)-4d0)/(eps-k-1d0))

    b17=k*(2d0+(k/(eps-k-1d0)))

    b18=(4d0*(eps**3))-(k*((8d0*(eps**2))-(7d0*(eps*k))+(2d0*(k**2)) &
    +(4d0*eps)-(5d0*k)+(2d0*eps/(eps-k-1d0))))

    b19=k*((2d0*eps)-(3d0*eps*k)-k+(2d0*((eps**2) &
    -(eps*(k**2))+k)/(eps-k-1d0)))-(4d0*(eps**3))

    b20=(4d0*eps*k)-(8d0*(eps**2))-(k**2)-(4d0*eps)-k

    b21=(7d0*eps)-(2d0*k)+5d0

    b22=eps*((4d0*(eps**2))+(2d0*eps*k)+(k**2)+(2d0*k))

    b23=eps*((eps*k)-(k**2)+k-2d0)

    b24=k*((4d0*(eps**2))-(2d0*(eps**3))+(k**2)-(eps*k)-(3d0*k))

    b25=(4d0*(eps**2))-(6d0*(eps*k))+(3d0*(k**2))-(9d0*eps)+(4d0*k)-1d0 &
    +(((k**2)+(4d0*k)-2d0)/(eps-1d0))+(2d0*k/((eps-1d0)**2))

    b26=4d0-(3d0*(eps-k))+((k+2d0)/(eps-1d0))

    b27=(3d0*((eps-k)**2))-((eps-k)**3)-((eps-k)*((eps**2) &
    -(3d0*eps)+(2d0*k/(eps-1))))+(k**2)+k-2d0

    !******

    chi=b23+(b24/(x+k))
    chi=chi*(1d0/(eps-1d0-x))
    chi=chi+b20-(2d0*(x**2))+(b21*x)+(4d0*eps*k/x)+(b22/(x+k))
    chi=chi*(L5h/W5h)

    chi=chi+((L6/W6)*(b25+(b26*x)+(x**2)+(b27/x)))
    chi=chi+((L5/W5)*(b16-(b17*x)+(b18/x)+(b19/(x+k))))
    chi=chi+((Lh/Wh)*(b9+x+(b10/x)+(b11/(eps-1d0-x)) &
    +((b12+(b13*x)-(b14*(x**2))+(b15/x))/(Wh**2))))
    chi=chi+((b4+(b5/x)+(b6/(x**2))+(((b7/x)-b8)/ &
    (((p**2)-(2d0*eps*k)+(k**2))+(2d0*x)))) &
    *L/(dsqrt(((p**2)-(2d0*eps*k)+(k**2))+(2d0*x))))
    chi=chi+((L3/p)*((b1/x)+(b2/(x+k))+(b3/(eps-1d0-x))))
    chi=chi-(2d0*k*L1/(x+k))

    if(ieee_is_nan(rho)) then
      chi=0
    end if

  end function chi

  real(8) function sigma(x,num)
    use ieee_arithmetic, only: ieee_is_nan
    implicit none
    integer, intent(in) :: num
    real(8), intent(in) :: x
    real(8) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14,&
      a15, a16, a17, a18, a19
    real(8) :: W, W1, W4, L1, L2, L3, L4, p, eps, k, k0
    common W1, W4, p, eps, k, k0

    !******
    W=dsqrt(((eps-k-x)**2)-1d0)
    L1=dlog(eps-k-x+W)
    L2=dlog(((((p**2)-(2d0*eps*k)+(k**2))+(W1*W))/x)-eps+k)
    L3=dlog((((p**2)+(p*W))/(x+k))-eps)
    L4=dlog(1d0+(((eps-k+1d0)*(((p**2)-(2d0*eps*k)+(k**2))-(x*(eps-k+1d0))+(W4*W))) &
    /(((p**2)-(2d0*eps*k)+(k**2))+(2d0*x))))

    W1=dsqrt((p**2)-(2d0*eps*k)+(k**2))

    sigma=0
    !******

    a1=(5d0*eps/(12d0*(k**3)))-(((16d0*eps)+13d0)/(12d0*(k**2))) &
    +((1d0/(6d0*k))*((8d0*eps)+35d0+(((15d0*eps)+14d0)/(p**2)))) &
    +(((eps**2)-(eps*k)+(19d0*eps)-(9d0*k))/(6d0*(p**2)))+((eps-(5d-1))/(p**4))

    a2=((2d0*(eps**2))+1d0)/(2d0*(p**4)*(k**3))

    a3=(1d0/((p**2)*(k**2)))*(((1d0/3d0)*((eps*k)+(8d0*eps)-1d0&
    +(eps/(2d0*k))))+(((4d0*eps)-(5d-1))/(p**2)))

    a4=((1d0/((p**2)*k))*((k*((eps/6d0)-(5d-1)))-(((4d0/3d0)*eps)*(eps-1)) &
    -(13d0/3d0)+((eps-(7d0/2d0))/(p**2))-(eps*((2d0*eps)+1d0)/(3d0*k)))) &
    +(1d0/(4d0*(k**3)))

    a5=(k/(p**2))*((3d0*(eps**2))-(eps*k)-(4d0*eps)+(9d0)+(8d0/(p**2)) &
     +((1d0/((p**2)-(2d0*eps*k)+(k**2)))*(((eps-k)*((eps**3)+(2.5d0*(p**2))-(5d0*eps/3d0))) &
     -(3d0*(p**2)/2d0)-(2d0*(k**2)/3d0)+(6d0*k)+((k/(p**2))*(8d0-(2d0*k) &
     -(4d0*eps*k)-(k*((2d0*(eps**2))-(k**2)+2d0)/((p**2)-(2d0*eps*k)+(k**2))))))))

    a6=((k**2)/(p**2))*((4d0*eps)-(1d0/3d0)+(((6d0*eps)-1d0)/(p**2)) &
    +(((eps*k)-(k**2)+2d0)/(3d0*((p**2)-(2d0*eps*k)+(k**2)))))

    a7=((2d0*(k**3))/(3d0*(p**4)))*((2d0*(eps**2))+1d0)

    a8=(2.5d-1*((p**2)-(eps*k)-(2d0*(k**2))))+(4d0*eps)+(5.25d0*k)+7d0 &
    -((k**2)/(eps+1d0))+((k**2)/(eps-k-1d0))-((1d0/(p**2))*((5d0*eps*k) &
    +(2.5d0*eps)-(5d0*k)+(4.5d0)))-((1d0/(p**4))*((2d0*eps*k) &
    +(4d0*eps)-(4d0*k)+(4d0)))-((1d0/k)*((8d0*(eps**2)/3d0) &
    +(2d0*eps)+(1d0/6d0)))+((eps+1d0)*(2d0*eps)/(3d0*(k**2)))&
    -(((2d0*(eps**2))+3d0)/(1.2d1*(k**3)))

    a9=(2d0*k*(eps+2d0))

    a10=(2.5d0*(p**2))-(6d0*(eps*k))+(5d0*(k**2))-(7d0*eps)+(1.5d1*k) &
    -(2d0*eps/(eps-k-1d0))-((1d0/(p**2))*(((k**2)*((eps*k)+(6d0*eps)-k-4d0)) &
    +((1d0/3d0)*eps*(eps+(3.5d1*k)))-(1d1*k)+5d0))+(((p**2)-(2d0*eps*k)+(k**2))*(((4d0 &
    -(2d0*eps))/(p**4))+((p**2)/(2d0*((p**2)+(2d0*k))))))+((k/((p**2)-(2d0*eps*k)+(k**2))) &
    *( ((3.5d0-(2d0*k))*(eps-k))+4.5d0+((1d0/(p**2))*( (4d0*(k**3)) &
    -(1.3d1*(eps*(k**2))/3d0)-(2d0*(eps**3)/3d0)-(2d0*eps*k)+eps &
    -(8d0*k)+4d0))))-((((eps*(eps-k))+1d0)*(k**2)/((p**2)*((((p**2)-(2d0*eps*k)+(k**2))**2)))))

    a11=((k/3d0)*((5d0*k)+7d0))-(6d0*eps*k)+((k/(p**2))*(eps &
    +((k+4d0)/3d0)))+((k/(p**4))*((1.05d1*(eps**3))+1d0))+(((2d0*eps) &
    -(k**3))/(eps-k-1d0))

    a12=(eps*((eps**2)+2d0))/(2d0*(p**4)*(k**3))

    a13=(1d0/(3d0*(p**2)*(k**2)))*((k*((2d0*(eps**2))-1d0))+(2d0*eps &
    *((2d0*eps)-1d0)*(1d0+(3d0/(p**2))))+((2d0*(eps**2))/k))

    a14=(1d0/(k*(p**2)))*(6d0-(2d0*eps)-k-(2d0*(eps**2)*(eps-k-1d0)) &
    +((5d0-(4d0*eps))/(p**2))-((2d0*(eps**3))/k))

    a15=5d0-(8d0*eps)+(5d0*k)+((1d0/(p**2))*((2d0*k)+1.2d1+(6d0/(p**2)))) &
    +((1d0/k)*((4d0*(eps**2))-(4d0*eps)+(2d0/(p**2))))-((k**2)/(eps-k-1d0))

    a16=k*((k*((4d0*eps)-(2d0*k)-4d0-((1.6d1+k)/(p**2))+(((8d0*(eps**3))-1d1) &
    /(p**4))))-(4d0*(eps**2))+(4d0*eps)-(2d0/(p**2))+((2d0*eps)/(eps-k-1d0)))

    a17=(eps*(k**2)/(p**2))*((2d0*(eps**2))-(k*((2d0*eps)-1d0) &
    *(1d0+(3d0/(p**2)))))

    a18=((2d0*eps*(k**3))/(3d0*(p**2)))*(eps-k-(3d0*k/(p**2)))

    a19=(((k**2)*(eps-k))/(2d0*((eps-k-1d0)**2)))-(2.5d-1*(eps+1d0))

    !******

    sigma=a1+(1d1*k/(x+k))-(a2*(x**3))+(a3*(x**2))+(a4*x)+(a5/x) &
    -(a6/(x**2))+(a7/(x**3))

    sigma=sigma*W

    sigma=sigma+(L1*(a8+(4d0*((x/k)-(k/x)))-(a9/(x+k))))

    sigma=sigma+(a10*(L2/W1))


    sigma=sigma+((L3/p)*(a11-(a12*(x**4))+(a13*(x**3))+(a14*(x**2))+(a15*x) &
    +(a16/x)+(a17/(x**2))-(a18/(x**3))))

    sigma=sigma+(a19*W4*L4)+hg(eps,k,k0,x,num)


    if((num.eq.2).and.(k.ge.k0)) then
      W=0
      L1=0
      L2=0
      L4=0
      sigma=((dsqrt((2d0*k)-(k**2)-(2d0*k*dsqrt(abs((2d0*k)-1d0))))/(4d0*k))*dasin(((1d0+((k+dsqrt(abs((2d0*k)-1d0))) &
            *(eps-k-x))))/(eps+dsqrt(abs((2d0*k)-1d0))-x)) &
            *(((1d0/dsqrt(abs((2d0*k)-1d0)))*((3d0*eps*k)-eps-k+1d0)) &
            +eps+(2d0*k)-1d0-(((eps/dsqrt(abs((2d0*k)-1d0)))-1d0) &
            *((p**2)*k)/((p**2)+(2d0*k)))))+(dsqrt((2d0*k)-(k**2) &
            +(2d0*k*dsqrt(abs((2d0*k)-1d0))))/(4d0*k)) &
            *dasin((1d0+((k-dsqrt(abs((2d0*k)-1d0)))*(eps-k-x))) &
            /(eps-dsqrt(abs((2d0*k)-1d0))-x)) &
            *(((1d0/dsqrt(abs((2d0*k)-1d0)))*(eps+k-(3d0*eps*k)-1d0))+eps+(2d0*k)-1d0 &
            +(((eps/dsqrt(abs((2d0*k)-1d0)))+1d0)*((p**2)*k)/((p**2)+(2d0*k))))
      if(ieee_is_nan(sigma)) then
        sigma=0
      end if
    end if
  end function sigma

  real(8) function hg(eps,k,k0,x,num)
    use ieee_arithmetic, only: ieee_is_nan
    implicit none
    integer, intent(in):: num
    real(8), intent(in) :: eps, k, k0, x
    real(8) :: W7, W8, W, L0, ww, r, s, u, v, p

    !******
    p=dsqrt((eps**2)-1d0)

    W=dsqrt(((eps-k-x)**2)-1d0)
    ww=dsqrt(abs((2d0*k)-1d0))
    r=dsqrt((5d-1*k)*((dsqrt((k**2)+(4d0*k)))+k-2d0))
    s=dsqrt((5d-1*k)*((dsqrt((k**2)+(4d0*k)))-k+2d0))
    u=(k*((eps-x)**2))+(((2d0*k)-(k**2)-(r*W))*(eps-x))-(2d0*(k**2))+k-(ww*s*W)
    W7=dsqrt((2d0*k)-(k**2)-(2d0*k*ww))
    W8=dsqrt((2d0*k)-(k**2)+(2d0*k*ww))
    L0=dlog((dsqrt((u**2)+(v**2)))/(((eps-x)**2)+(2d0*k)-1d0))
    v=(ww*((eps-x)**2))-(((2d0*k*ww)+(s*W))*(eps-x))+(ww*((k**2)-1d0+(r*W)))
    hg=0
    if((num.eq.2).and.(k.ge.k0)) then
      L0=0
      v=0
    end if
    !******
    if (k<5d-1) then
      hg=(W7/(4d0*k))*dasin(((1d0+((k+ww)*(eps-k-x))))/(eps+ww-x)) &
      *(((1d0/ww)*((3d0*eps*k)-eps-k+1d0))+eps+(2d0*k)-1d0 &
      -(((eps/ww)-1d0)*((p**2)*k)/((p**2)+(2d0*k))))
      if(ieee_is_nan(hg)) then
        
        hg=0
      end if
      hg=hg+(W8/(4d0*k))*dasin((1d0+((k-ww)*(eps-k-x)))/(eps-ww-x)) &
      *(((1d0/ww)*(eps+k-(3d0*eps*k)-1d0))+eps+(2d0*k)-1d0 &
      +(((eps/ww)+1d0)*((p**2)*k)/((p**2)+(2d0*k))))
      if(ieee_is_nan(hg)) then
        hg=0
      end if
    else if (k==5d-1) then
      hg=(((3d0*eps/2d0)+(1d0/4d0)-(1d0/(2d0*eps)) &
      -(3d0/(4d0*(eps**2))))/(dsqrt(3d0)))
      hg=hg*dasin((5d-1)+(3d0/(4d0*(eps-x))))
      if(ieee_is_nan(hg)) then
        hg=0
      end if
      hg=hg+(W*(eps+1d0)/(2d0*eps*(eps-x)))
    else if (k>5d-1) then
      hg=((r/ww)*(eps+k-(3d0*eps*k)-1d0))+(s*(eps+(2d0*k)-1d0)) &
      +(((eps*r/ww)+s)*(p**2)*k/((p**2)+(2d0*k)))
      hg=hg*(1d0/(2d0*k))*datan(v/u)
      if(ieee_is_nan(hg)) then
        hg=0
      end if
      hg=(L0/(2d0*k))*(((s/ww)*(eps+k-(3d0*eps*k)-1d0)) &
      -(r*(eps+(2d0*k)-1d0))+(((eps*s/ww)-r)*(p**2)*k/((p**2)+(2d0*k))))
    else
      print *, 'k not defined'
    end if
    return
  end function hg

  real(8) function eeblsnr(ee,ep) !nonrelativistic limit
    real(8), intent(in)  :: ee, ep
    real(8) :: eps, p, k

    eps=1+(ee/me)
    k=ep/me
    p=dsqrt((eps**2)-1d0)

    eeblsnr=((1.3d1/8d0)*dlog((p**2)/k))+(1.25d2/5.6d1)
    eeblsnr=eeblsnr*(p**2)*(-1d0)
    eeblsnr=eeblsnr+(1.7d1/6d0)+((7d0/3d0)*dlog((p**2)/k))
    eeblsnr=(-(k**2)/(p**4))*eeblsnr
    eeblsnr=eeblsnr-(((2d0*dlog((p**2)/k))+(2.9d1/6d0)+(((6.2d1* &
    (dlog((p**2)/k))/3d0)-5d0)*(p**2)/7d0))*k/(p**2))
    eeblsnr=eeblsnr+((p**2)*(((1.3d1/2.8d1)*(dlog((p**2)/k)))-(3.91d2/3.36d2)))
    eeblsnr=eeblsnr+(1.7d1/1.2d1)+dlog((p**2)/k)
    eeblsnr=eeblsnr*(1.6d1/5d0)*(afs*(r0**2))/(k)
  end function eeblsnr

  real(8) function eeblser(ee,ep) !extremely relativistic approximation
    real(8), intent(in)  :: ee, ep
    real(8) :: eps, p, k

    eps=1d0+(ee/me)
    k=ep/me
    p=dsqrt((eps**2)-1d0)

    eeblser=(dlog(2d0*(eps**3)/(k**2)))+(1d0/1.2d1) &
    -(k*((dlog(eps/k))+(2.3d1/2.4d1)))
    eeblser=eeblser+((k**2)*((dlog(eps/k))+(8d0/1.5d1)))
    eeblser=eeblser*(1.6d1/3d0)*(afs*(r0**2))/(k)
  end function eeblser

  real(8) function eebapx(ee,ep) 
  
    !Haug 1998 approximation by expansion on the three leading terms of a series in k
    !Long-Wavelength limit eeb approximation
    !Includes the three leading terms of a series in k

    real(8), intent(in)  :: ee, ep
    real(8) :: eps, p, k
    real(8) :: eebk, eebkk

    eps=1d0+(ee/me)
    k=ep/me
    p=dsqrt((eps**2)-1d0)

    eebapx=dlog((2d0*(eps-1d0)*(p**2))/(k**2))
    eebapx=eebapx*(2d0*eps*(2d0+((eps+1d0)/(p**2))))
    eebapx=eebapx-(eps**2)+(6d0*eps)-(4d0/3d0)+(((5d0*eps)+(8d0/3d0))/(p**2))
    eebapx=eebapx*((4d0/p)*dlog(eps+p))
    eebapx=eebapx+(8d0*(eps**2)*((2d0/3d0)-(eps/(p**2))) &
    *dlog((2d0*(eps-1d0)*(p**2))/(k**2)))

    eebapx=eebapx-((8d0*eps/p)*(2d0+(1d0/(p**2))) &
    *((2d0*dlog(eps+p)*dlog(2d0*(p**2)/k))+dilog(-eps-p)-dilog(-eps+p)))
    eebapx=eebapx+((eps**2)*4d0/9d0)+(eps*4d0/3d0)-(3.2d1/3d0) &
    -(((5.6d1*eps/3d0)+1.2d1)/(p**2))

    eebk=-(((8d0*eps)+(6d0/(p**2)))*(dlog((2d0*(eps-1d0)*(p**2)) &
    /(k**2)))/(eps-1d0))
    eebk=eebk-(2d0*p*dlog(eps+p))+(8d0*(eps**3)/3d0)-(8d0*(eps**2)) &
    +(2d1*eps/3d0)-(1d1/3d0)+((1.55d2-(8.9d1*eps))/(3d0*(p**2))) &
    +((3.6d1-(1.6d1*eps))/(p**4))
    eebk=eebk*dlog(eps+p)/p
    eebk=eebk+((6d0-(8d0*eps*(eps+1d0)/3d0)+(4d0*eps/(p**2)) &
    +(((1.2d1*eps)-6d0)/(p**4)))*dlog((2d0*(eps-1d0)*(p**2))/(k**2)))
    eebk=eebk+((2d0/p)*(4d0-(2d0*eps)+(((3d0*eps)-2d0)/(p**2)) &
    +(((6d0*eps)-3d0)/(p**4)))*((2d0*dlog(eps+p)*dlog(2d0*(p**2)/k)) &
    +dilog(-eps-p)-dilog(-eps+p)))
    eebk=eebk-(4.6d1*(eps**2)/9d0)+(2.8d1*eps/9d0)+1.1d1 &
    -(((7.1d1*eps/3d0)-2.1d1)/(p**2))-(((1.8d1*eps)+2d0)/(p**4))

    eebkk=(-8d0+((1d1-(6d0*eps))/(p**2))+((2.6d1+(2.4d1*eps))/(p**4)) &
    +((4d0*(eps+1d0))/(p**6)))
    eebkk=eebkk*dlog((2d0*(eps-1d0)*(p**2))/(k**2))
    eebkk=eebkk-(8d0*(eps**3)/3d0)+(5.2d1*(eps**2)/5d0)-(2d0*eps/3d0) &
    -(1.3d1/1.5d1)+(9.7d1*eps/(3d0*(p**2)))-(2.29d2/(5d0*(p**2))) &
    +(((1.594d3*eps)-1.083d3)/(5d0*(p**4)))+(((2.81d2*eps)-(4.64e2/3d0))/(p**6))
    eebkk=eebkk*dlog(eps+p)/p
    eebkk=eebkk+((dlog((2d0*(eps-1d0)*(p**2))/(k**2)))*((8d0*(eps**2)/3d0) &
    -(6.7d1*eps/1.5d1)+(1d1/3d0)+(4.6d1/(3d0*(p**2)))-(2d0*eps/(p**2)) &
    +((2.68d2-(9.4d1*eps))/(3d0*(p**4)))+((7.6d1-(3.4d1*eps))/(p**6)) ))
    eebkk=eebkk+((1d0/(p**3))*((8d0*(eps**2))+((2.8d1-(8.8d1*eps))/(p**2)) &
    +((3d1-(8d1*eps))/(p**4)))*((2d0*dlog(eps+p)*dlog(2d0*(p**2)/k)) &
    +dilog(-eps-p)-dilog(-eps+p)))
    eebkk=eebkk+(1.28d2*(eps**2)/4.5d1)-(4.49d2*eps/2.25d2)+(9.9d1/5d0) &
    -(((9.3d1*eps)+3.43d2)/(5d0*(p**2)))+(((4.79d2*eps) &
    -(1.2041d4/3d0))/(1.5d1*(p**4)))+(((1.34d2*eps/3d0)-1.71d2)/(p**6))

    eebapx=eebapx+(k*eebk)+((k**2)*eebkk)
    eebapx=eebapx*(afs*(r0**2))/(k*(p**2))
  end function eebapx

  real(8) function cspe(ee,ep) !d(sigma)/dk cross section in photon energy (KOCH-MOTZ)
    real(8) ep, ee
    real(8) p0, p, E0, E, eps, eps0, k, L, me
    integer, parameter :: rel=1
    if(ep.lt.ee) then
      E0=ee/me+1d0
      k=ep/me
      E=E0-k

      p0=dsqrt(E0**2-1d0)
      p=dsqrt(E**2-1d0)

      eps0=dlog((E0+p0)/(E0-p0))
      eps=dlog((E+p)/(E-p))
      L=2*dlog((E0*E+p0*p-(1d0))/(k))
      !print *, E0, E, p0, p, eps0, eps
      cspe=(eps0*((E0*E+p0**2)/(p0**3)))+((2*k*E0*E)/((p**2)*(p0**2))) &
      -(eps*((E0*E+p**2)/(p**3)))
      cspe=(k/(2*p0*p))*cspe
      cspe=(8.)*((E0*E)/(3.*p0*p))+(((k**2)*((E0**2)*(E**2) &
      +(p0**2)*(p**2)))/((p0**3)*(p**3)))+cspe

      cspe=L*cspe

      cspe=(4./3.)-(2*E0*E*((p**2+p0**2)/((p0**2)*(p**2))))+((eps0*E)/(p0**3)) &
      +((eps*E0)/(p**3))-((eps*eps0)/(p0*p))+cspe
      cspe=(p/(p0*k))*cspe
      cspe=(r0**2)*afs*cspe!*(1/(p0**2))*dlog((p0+p)/(p0-p))
      ! cspe=cspe/(511)
      if(rel.eq.0) then
        cspe=(dlog((1d0+dsqrt(1d0-(ep/ee)))/(1d0-dsqrt(1d0-(ep/ee)))))/((ep)*(ee)) !NRBH
      end if
    else
        cspe=0
    end if
  end function cspe

  real(8) function cspax(ee,ep,z)
    ! Approximated equation for e-i bremsstrahlung cross-section (equation 4 in Haug 1997a)
    ! Accurate to semirelativistic energies. 
    ! https://ui.adsabs.harvard.edu/abs/1997A%26A...326..417H/abstract
    ! See IDL version from SSWIDL xray module:
    ! https://hesperia.gsfc.nasa.gov/ssw/packages/xray/idl/brm/brm_bremcross.pro

    real(8), intent(in) :: ee, ep, z
    real(8) :: eps1, eps2, p1, p2, k
    real(8) :: term1, term2, factor1, factor2, factor3
    real(8) :: a_1, a_2, elw_factor

    if(ep.lt.ee) then
      eps1 = ee/me + 1d0
      k = ep/me
      eps2 = eps1 - k
      
      if (eps2 <= 1d0) then
        cspax = 0d0
        return
      end if

      p1 = dsqrt(eps1**2 - 1d0)
      p2 = dsqrt(eps2**2 - 1d0)

      term1 = 1d0 + (1d0/(eps1*eps2)) + (7d0/2d1)*((p1**2 + p2**2)/((eps1*eps2)**3))
      term1 = term1 + (((9d0/28d0)*(k**2) + (263d0/21d1)*(p1**2 * p2**2))/((eps1*eps2)**3))
      term1 = term1 * ((p1 * p2)/(eps1*eps2))
      
      term2 = 2d0 * dlog((eps1*eps2 + p1*p2 - 1d0)/(k))

      factor1 = term2 - term1

      factor2 = (4d0/3d0)*(eps1*eps2) + k**2 - (7d0/15d0)*(k**2 / (eps1*eps2))
      factor2 = factor2 - (11d0/7d1)*((k**2)*(p1**2 + p2**2)/((eps1*eps2)**4))

      factor3 = 2d0 * afs * (z**2) * (r0**2) / (k * (p1**2))

      a_1 = afs * z * eps1 / p1
      a_2 = afs * z * eps2 / p2
      elw_factor = (a_2 / a_1) * (1d0 - dexp(-2d0 * pi * a_1)) / (1d0 - dexp(-2d0 * pi * a_2))
      
      cspax = elw_factor * factor1 * factor2 * factor3
      
    else 
      cspax=0d0
    end if

  end function cspax

  real(8) function cs1(ep,ee,ca,z)
    !GLUCKSTERN-HULL (d sigma / dk d Omega)
    real(8) ep, ee, ca ! ca is pitch angle cosine
    real(8) z  ! average ion charge
    real(8) p0, p1, g0, g1, k, T, L, eps, epst, delta, Qe, b0, b1


    k=ep/me!(me*(c**2))
    g0=1+ee/me!(me*(c**2))
    g1=g0-k
    p0=dsqrt(g0**2-1.)
    p1=dsqrt(g1**2-1.)
    T=dsqrt((p0**2)+(k**2)-2.*p0*k*ca) !cos(psi)=ca
    L=dlog((g0*g1-1.+p0*p1)/(g0*g1-1.-p0*p1))
    eps=dlog((g1+p1)/(g1-p1))
    epst=dlog((T+p1)/(T-p1))
    delta=g0-p0*ca
    b0=p0/g0
    b1=p1/g1

    Qe=(b0*(1-dexp(-(2*pi*z/137)/b0)))/(b1*(1-dexp(-(2*pi*z/137)/b1)))

    !VERSION 3
    cs1=4.*(g0**2)*((g0**2)+(g1**2))
    cs1=cs1-(2.*(7.*(g0**2)-3.*(g0*g1)+(g1**2)))
    cs1=cs1+2.
    cs1=cs1/((p0**2)*(delta**2))
    cs1=cs1+2.*k*((g0**2)+(g0*g1)-1.)/((p0**2)*(delta))
    cs1=cs1+4.*(g0*(1.-(ca**2)))*(3.*k-((p0**2)*g1))/((p0**2)*(delta**4))
    cs1=cs1*(L/(p1*p0))
    cs1=cs1+8.*(1.-(ca**2))*(2.*(g0**2)+1.)/((p0**2)*(delta**4))
    cs1=cs1-2.*(5.*(g0**2)+2.*(g1*g0)+3.)/((p0**2)*(delta**2))
    cs1=cs1-2.*((p0**2)-(k**2))/((T**2)*(delta**2))
    cs1=cs1+4.*g0/(delta*(p0**2))
    cs1=cs1+(epst/(p1*T))*((4./(delta**2))-(6.*(k/delta))-(2.*k*((p0**2)-(k**2)))/((T**2)*delta))
    cs1=cs1-(4.*eps/(p1*delta))

    cs1=cs1*(p1/p0)
    cs1=cs1*(1./ep)
    cs1=cs1*((z**2)*(r0**2)*afs)/((8*pi))!*Qe

    if(ep.ge.ee) then
      cs1=0
    end if
  end function cs1

  real(8) function eeba(ee,ep,ca) 

    ! d2(sigma)/dk dOmega in lab system (Haug 1975a)
    ! Bremsstrahlung and pair production in the field of free electrons, 
    ! Z. Naturforsch., A: Phys. Sci., 1975a, vol. 30, no. 9, pp. 1099-1113
    ! Calculate double differential electron-electron bremsstrahlung cross-section in lab. system
    ! Input: electron kinetic energy, photon energy in keV and angle cosine
    
    use ieee_arithmetic, only: ieee_is_nan
    double precision ee,ep,ca
    double precision eps1,p1,k !normalized units
    double precision x1,x2,omega2,rho2,omega,rho
    double precision a1,a2,beta1,beta2,Fee !Coulomb Correction

    eps1=(ee/me)+1d0
    p1=dsqrt((eps1**2)-1d0)
    k=(ep/me)
    omega2=2d0*(eps1+1d0)
    rho2=2d0*(eps1+1d0-(k*(eps1+1d0-(p1*ca))))
    omega=dsqrt(omega2)
    rho=dsqrt(rho2)
    x1=k*(eps1-(p1*ca))
    x2=k

    beta1=omega*dsqrt(omega2-4d0)/(omega2-2d0)
    beta2=rho*dsqrt(rho2-4d0)/(rho2-2d0)
    a1=afs/beta1
    a2=afs/beta2
    Fee=(a2*(exp(2d0*pi*a1)-1d0))/(a1*(exp(2d0*pi*a2)-1d0))

    eeba=alpha(omega2,rho2,x1,x2)+alpha(omega2,rho2,x2,x1) !dimensionless
    if(ieee_is_nan(eeba)) then
      eeba=0
    else
      eeba=eeba*afs*(r0**2)*k/(pi*omega*rho*dsqrt(omega2-4d0))
      eeba=eeba*Fee
      eeba=eeba/me !d(sigma)/ d ep d(Omega) (Multiply by two???)
    end if
  end function eeba

  real(8) function alpha(omega2,rho2,x1,x2)

    double precision rho,omega
    double precision x,R1,R2,W2,W4,L,L1,L2,L3,L4
    double precision , intent(in)  :: omega2,rho2,x1,x2
    double precision eeba1,eeba2,eeba3,eeba4,eeba5

    omega=dsqrt(omega2)
    rho=dsqrt(rho2)
    x=(omega2-rho2)/(2d0)

    R1=rho2-4d0+(4d0*x1)+(4d0*(x1**2)/rho2)
    R2=rho2-4d0+(2d0*x1)
    W2=dsqrt(x2*((x2*(rho2-4d0)/4d0)+(2d0*x*x1/rho2)))
    W4=dsqrt((omega2-4d0)*((2.5d-1*(omega2-4d0)*(rho2-4d0))+(4d0*x1*x2/rho2)))
    L=dlog(rho*(R2+dsqrt((rho2-4d0)*R1))/(4d0*x1))
    L1=dlog(5d-1*(rho+dsqrt(rho2-4d0)))
    L2=dlog(1d0+((rho2/(4d0*x*x1))*((x2*(rho2-4d0))+(2d0*dsqrt(rho2-4d0)*W2))))
    L3=dlog((((omega*dsqrt(rho2-4d0))+(rho*dsqrt(omega2-4d0)))**2)/(4d0*(omega2-rho2)))
    L4=dlog(1d0+((rho2/(8d0*x1*x2))*(((omega2-4d0)*(rho2-4d0))+(2d0*dsqrt(rho**2-4d0)*W4))))

    eeba1=omega2-(5d-1*rho2*(omega2-2d0))
    eeba1=eeba1*(1d0/(x1**2))
    eeba1=eeba1+((omega2-(4d0*x2))/(rho2))-((omega2*(omega2-4d0))/(4d0*x1*x2))-(4d0/x1)
    eeba1=eeba1/R1
    eeba1=eeba1+((((x1-x2)/x)**2)*(omega2+rho2)/(4d0*x1*x2))
    eeba1=eeba1-(2.5d-1*(((1d0/x1)-(1d0/x2))**2))-(rho2/(2d0*(x**2)))
    eeba1=eeba1+((1d0+(1d0/(omega2-4d0)))*2d0*rho2/((omega2-4d0)*x1*x2))
    eeba1=eeba1+(4d0*rho2*((3d0*x1*x2*((omega2-2d0)**2)/(omega2-4d0))-(2d0*(x**2)*(1d0 &
    +(6d0/(omega2*(omega2-4d0))))))/(omega2*(omega2-4d0)*(x1**4)))
    eeba1=eeba1-(4d0*rho2*x/((omega2-4d0)*(x1**3)))
    eeba1=eeba1+(rho2*((4d0/omega2)-(3d0/2d0)-(8d0/(omega2*((omega2-4d0)**2))) &
    +(x*(omega2-2d0)/(omega2*(omega2-4d0))))/(x1**2))
    eeba1=eeba1-(rho2/(x1*(omega2-4d0)))
    eeba1=eeba1*dsqrt(rho2-4d0)

    eeba2=(omega2*R2/(4d0*x2))-(rho2-2d0)-(2d0*x1)
    eeba2=eeba2*(omega2-4d0-(4d0*x1*x2/rho2))/(x1*R1)
    eeba2=eeba2+(4d0*(x2-(3d0*x1)+4d0+(2d0*(rho2-3d0)/x1))/(x1*R2))
    eeba2=eeba2+(3d0*omega2*(omega2-4d0)/(4d0*x1*x2))
    eeba2=eeba2+4d0+((1d1-(omega2/2d0))/x2)-(2d0*((2d0*omega2)-x2+4d0)/x1)
    eeba2=eeba2*L/(dsqrt(R1))
    eeba2=eeba2+(rho*L1*(((rho2+2d0)/(x**2))+(8d0/(x1**2))))

    eeba3=(x2*(((omega2+rho2)/2d0)-2d0))-(omega2-2d0)
    eeba3=eeba3*(rho2-2d0)
    eeba3=eeba3-(2d0*(x2**2))-(4d0*x2)
    eeba3=eeba3/(R2*x1)
    eeba3=eeba3+(((5d-1*(rho2-2d0)*((3d0*(rho2-4d0))-(omega2*(rho2-5d0))))+(2d0*(x1**2)) &
    -(6d0*x1)+(x2*(omega2-2d0)))/(R2*x))
    eeba3=eeba3+(2d0*(omega2-2d0)/x2)-((rho2-2d0-x2)/x1)+((rho2-2d0)*(rho2+x1)/(2d0*x)) &
    +((omega2-4d0)/R2)-2d0+((rho2-2d0)*((omega2+rho2-4d0)**2)/(8d0*x*x1))
    eeba3=eeba3*L2/W2

    eeba4=(1.2d1*x/(omega2*(omega2-4d0)))-(5d-1*(rho2-2d0))
    eeba4=eeba4*((4d0*(omega2-2d0)*x)/((x1**2)*(omega2-4d0)))
    eeba4=eeba4*(1d0-(x/(omega2*x1)))
    eeba4=eeba4+(2d0*((omega2*(omega2-2d0)*(rho2-4d0))-(2d0*(rho2-2d0)) &
    +(4d0*rho2/omega2))/(x1*((omega2-4d0)**2)))
    eeba4=eeba4+(((5d-1*(rho2-2d0)*((3d0*(rho2-4d0))-(omega2*(rho2-5d0)))) &
    -(x1*(omega2-(2d0*x1)+4d0)))/(x*R2))
    eeba4=eeba4-(((2d0*x2*(rho2-2d0))+(x*(omega2-4d0))-(4d0*(omega2-2d0)) &
    +((8d0-rho2)*(omega2-2d0)/(2d0*x1)))/R2)
    eeba4=eeba4+(4d0*(x**2)/(omega2*(omega2-4d0)*x1))
    eeba4=eeba4+(x*(((omega2-2d0)**2)+((rho2-2d0)*(rho2-4d0)) &
    -(8d0*(rho2-2d0)/(omega2-4d0)))/(4d0*x1*x2))
    eeba4=eeba4+4d0+(8d0*(omega2-2d0)/((omega2-4d0)**2))-(3d0*(omega2-2d0)/(2d0*x2)) &
    +(((x2**2)-(2d0*omega2)-(x2*(omega2-1d0))+(5d-1*((omega2-2d0)**2)))/x1)
    eeba4=eeba4+(rho2*(rho2-2d0)/(2d0*x))-(x2*(omega2-2d0)/(2d0*x)) &
    -(((omega2-2d0)**2)*(omega2+rho2-4d0)/(4d0*x*x2))
    eeba4=eeba4*(2d0*rho*L3)/(omega*x1*dsqrt(omega2-4d0))

    eeba5=((omega2-2d0)**2)+((rho2-2d0)**2)-(6d0*(omega2+rho2-4d0))+(1.6d1*x/(omega2-4d0))
    eeba5=eeba5*((rho2-2d0)/(8d0*x1*x2))
    eeba5=eeba5+1d0+(2d0*(1d0-x1-(x1**2))/(x1*x2))+((rho2-4d0-(8d0/(omega2-4d0)))/(omega2-4d0))
    eeba5=eeba5*L4/W4

    alpha=eeba1+eeba2+eeba3+eeba4-eeba5
  end function alpha

  real(8) function eecskl(ee,ep)
    integer i,n
    real(8) :: ee,ep,eps,p,k,eecskl1,eecskl2
    eps=1.+ee/me
    k=ep/me
    p=dsqrt((eps**2)-1.)

    eecskl=(1.+((p**2)*13./(7)))*(16./(5.))
    eecskl=eecskl*dlog(4.*((p**2)/(k)))
    eecskl=eecskl+(68./(15.))-((p**3)*444./(35.))

    eecskl1=(8./(5.))*(1.+((p**2)*73./(14.)))*dlog((4.*(p**2))/k)
    eecskl1=eecskl1+(58./(15.))-((p**2)*33./(7.))
    eecskl1=eecskl1*(k/(p**2))

    eecskl2=(7./(15.))-((p**2)*46./(15.))
    eecskl2=eecskl2*dlog((4.*(p**2))/k)
    eecskl2=eecskl2+(17./(30.))-(97./(42.))*(p**2)
    eecskl2=eecskl2*((k**2)/(p**4))

    eecskl=eecskl-eecskl1-eecskl2
    eecskl=eecskl*(afs*(r0**2))/(k)
  end function eecskl

  real(8) function eecsk(ee,ep)
    real(8) eecsk1,eecsk2
    real(8) :: ee,ep,eps,p,k,jjvalue
    eps=1.+ee/me
    k=ep/me
    p=dsqrt((eps**2)-1.)

    jjvalue=jj(eps)

    eecsk1=(-32./3.)+(1./(6.*(p**2)))+(1./((eps**2)*(p**4)))
    eecsk1=eecsk1-((3.*(eps**4))/(4*(p**6)))-((3.*(p**2))/(4.*(eps**4)))
    eecsk1=(1./eps)*eecsk1
    eecsk1=eecsk1+(1./(4.*p))*dlog(eps+p)*((1./(p**4))-(1./(eps**4))+(3./(p**6))+(3./(eps**6)))
    eecsk1=eecsk1*dlog(4.*eps*(p**2)/(k))
    eecsk1=eecsk1-(16./3.)*((p**2)/(eps**3))-(2.)/(eps*(p**2))+(2.)/(3.*eps*(p**4))
    eecsk1=eecsk1-(1./(2.*(eps**5)))-(1./(eps*(p**6)))
    eecsk1=eecsk1+(1./((eps**2)*p))*dlog(eps+p)*(-2.+(3./(eps**2))-(1./((eps**2)*(p**2)))&
    -(1./(2.*(eps**4)))+(1./(2.*(eps**2)*(p**4)))+(1./(p**6)))
    eecsk1=eecsk1-(jjvalue/(2.*(eps**2)*(p**4)))*(1.+(3./(4.*(eps**2)))+(3./(4.*(p**2))))
    eecsk1=eecsk1*k

    eecsk2=(7./(4.*(eps**2)))-(1./(2.*(eps**4)))+(5./(8.*(eps**6)))&
    -(39./(4.*(p**2)))+(1./(4.*(p**4)))+(5./(8.*(p**6)))-(7./(8.*(p**8)))
    eecsk2=eecsk2*(1./(eps*p))*dlog(eps+p)
    eecsk2=eecsk2+(49./(4.*(p**2)))+(1./(6.*(p**4)))-(29./(24.*(p**6)))+(7./(8.*(p**8)))&
    -(17./(4.*(eps**2)))-(37./(12.*(eps**4)))+(5./(8.*(eps**6)))
    eecsk2=eecsk2*dlog(4.*eps*(p**2)/k)
    eecsk2=eecsk2+(4./(3.*(eps**2)))+(7./(2.*(eps**4)))-(9./(16.*(eps**6)))&
    -(29./(48.*(eps**2)*(p**4)))-(5./(4.*(eps**4)*(p**4)))-(6./(p**6))-(31./(16.*(eps**2)*(p**8)))
    eecsk2=eecsk2+(1./(eps*p))*dlog(eps+p)*((2.*(eps**4)/(p**8))-((eps**2)/(16.*(p**8)))&
    +(19./(8.*(p**6)))+(4.*(eps**2)/(p**4))+(9./(4.*(eps**4)*(p**4)))-(2./(eps**2))-(9./(16.*(eps**6))))
    eecsk2=eecsk2+(jjvalue/(eps*(p**2)))*(2.-(33./(32.*(p**2)))-(9./(4.*(p**4)))&
    -(23./(32.*(p**6)))-(15./(32.*(eps**2)))+(9./(32.*(eps**4))))
    eecsk2=eecsk2*(k**2)

    eecsk=(32./3.)+(8./(3.*(eps**2)*(p**2)))-((((eps**2)+(p**2))**3)/((eps**4)*(p**4)))&
    +((((eps**2)+(p**2))**2)/((eps**5)*(p**5)))*dlog(eps+p)
    eecsk=eecsk*dlog(4.*eps*(p**2)/k)
    eecsk=eecsk-(16./3.)+(4./(eps**2))-(2.*(eps**2)/(p**4))
    eecsk=eecsk+(2./(eps*p))*(1.+(3./(eps**2))+((eps**4)/(p**4)))*dlog(eps+p)&
    +(jjvalue/(eps*(p**2)))*((1./(eps**2))- 4.-(1./(2.*(eps**2)*(p**2))))

    eecsk=eecsk+eecsk1+eecsk2
    eecsk=eecsk*(afs*(r0**2)/k)
  end function eecsk

  real(8) function eecs(ee,ep, thetaparts) 
  
    ! PRECISE e-e bremsstrahlung Haug 1989 (center-of-mass system)
    ! https://ui.adsabs.harvard.edu/abs/1989A%26A...218..330H/abstract
    ! Equation (A1)
    
    implicit none
    integer i
    integer, intent(in) :: thetaparts
    real(8), intent(in) :: ee, ep
    real(8), dimension(thetaparts) :: theta,ctheta,k1,k2,W2,W4,L,L2,L4,&
    R1,eecs3t,eecs3t1,eecs3t2,eecs3t3
    real(8) :: eps,p,k,eecs1,eecs2,eecs3
    real(8) :: L1,L3
    
    eecs1=1
    eecs2=1
    eecs3=1
    eecs3t1=0
    eecs3t2=0
    eecs3t3=0
    eps=1.+ee/me
    k=ep/me
    p=dsqrt((eps**2)-1.)
    L1=dlog(dsqrt((eps**2)-(eps*k))+dsqrt((p**2)-(eps*k)))
    L3=dlog((((eps*dsqrt((p**2)-(eps*k)))+(p*dsqrt((eps**2)-(eps*k))))**2)/(eps*k))
    
    do i=1,thetaparts
      theta(i)=pi-(i)*((pi)/thetaparts)
      ctheta(i)=dcos(theta(i))
      k1(i)=k*(eps-(p*ctheta(i)))
      k2(i)=k*(eps+(p*ctheta(i)))
      W2(i)=dsqrt(((eps*k2(i))/(eps-k))*((2*(k**2))+(k2(i)*((p**2)+(k**2)&
      -(2.*eps*k)))))
      W4(i)=(2.*p)*dsqrt((4*(p**2)*((p**2)-(eps*k)))+((k1(i)*k2(i))&
      /((eps**2)-(eps*k))))
      R1(i)=(4.*((p**2)-(eps*k)))+(4.*k1(i))+(((k1(i))**2)/((eps**2)-(eps*k)))
      L(i)=dlog(((dsqrt((eps**2)-(eps*k)))/(k1(i)))*((2.*(p**2))-(k2(i))&
      +dsqrt(R1(i)*((p**2)-(eps*k)))))
      L2(i)=dlog(1.+(2.*((eps-k)/(k*k1(i))))*((k2(i)*((p**2)-(eps*k)))&
      +(W2(i)*dsqrt((p**2)-(eps*k)))))
      L4(i)=dlog(1.+(2.*((eps**2)-(eps*k))/(k1(i)*k2(i)))*(((4.*(p**2))&
      *((p**2)-(eps*k)))+(W4(i)*dsqrt((p**2)-(eps*k)))))


  	  !*************************************************************************

      eecs3t1(i)=(2.*((eps**2)+(p**2)-(2.*eps*k))*((eps*k)+ 1.-(4.*(p**2)*((p**2)&
      -(eps*k)))))+(k1(i)**2)-(3.*(k1(i)))+(((eps**2)+(p**2))*(k2(i)))
      eecs3t1(i)=eecs3t1(i)*(1./(4.*(eps*k)*((2.*(p**2))-(k2(i)))))
      eecs3t1(i)=eecs3t1(i)+(1./((k1(i))*((2*(p**2))-k2(i))))*(((eps**2)+(p**2)&
      -(2.*eps*k))*((k2(i)*((eps**2)+(p**2)-(eps*k)))-((eps**2)&
      +(p**2)))-k2(i)-((k2(i)**2)/(2.)))
      eecs3t1(i)=eecs3t1(i)+((p**2)/((2.*(p**2))-(k2(i))))+(2.*((eps**2)&
      +(p**2))/(k2(i)))-(3./(2.))
      eecs3t1(i)=eecs3t1(i)+((eps**2)-(eps*k)+((k1(i))/(4.)))*(((eps**2)&
      +(p**2)-(2.*eps*k))/(eps*k))
      eecs3t1(i)=eecs3t1(i)*((L2(i))/(W2(i)))

      !*************************************************************************

      eecs3t2(i)=(eps/(k2(i)*(R1(i))))*((eps*k1(i))-((eps-k)*(k2(i))))*((2.*(p**2)&
      /(k1(i)))+(((eps**2)+1.)/((eps**2)-(eps*k))))
      eecs3t2(i)=eecs3t2(i)+(((4.*(p**2))-(3.*eps*k)+(1./(k1(i))))/(k1(i)&
      *((2.*(p**2))-k2(i))))
      eecs3t2(i)=eecs3t2(i)+((2.-(eps**2))/(k2(i)))+(2./((k1(i))**2))&
      +((((4.*eps)-k)*(eps*(p**2)))/(k1(i)*k2(i)))
      eecs3t2=eecs3t2(i)*((2.*L(i))/(dsqrt(R1(i))))

  	  ! *************************************************************************

      eecs3t3(i)=(4*(p**2)*((p**2)-(eps*k)))+((eps**2)*((2*(k**2))-2.+((eps*k)/(p**2))))
      eecs3t3(i)=eecs3t3(i)*(((eps**2)+(p**2)-(2*eps*k))/(k1(i)*k2(i)))
      eecs3t3(i)=eecs3t3(i)+1.-((eps*k)/(2.*(p**2)))-(1./(4.*(p**4)))&
      +((1.-k1(i)+(k1(i)**2))/(k1(i)*k2(i)))
      eecs3t3(i)=eecs3t3(i)*((L4(i))/(W4(i)))

  	  !*************************************************************************
  	  !*************************************************************************

      eecs3t(i)=eecs3t1(i)+eecs3t2(i)-eecs3t3(i)
    end do

    eecs3=integrate(ctheta,eecs3t)
    eecs3=eecs3*((k**2)/(dsqrt((eps**2)-(eps*k))))

    eecs1=(eps-k)*(1.+((4.-(2.*eps*k))/(p**2))+(1./(p**4))+(k/((eps*(p**2))+k)))
    eecs1=eecs1+(1./(eps))+(2.*eps)
    eecs1=eecs1+((p**2)/((p**2)-(eps*k)))*((4.*k)-eps+(2.*k/((p**2)-(eps*k))))
    eecs1=eecs1*(1./p)*(dlog(eps+p))
    eecs1=eecs1-(eps-k)*((eps*16./(3.))+(2.*((2.*eps)-k)/(p**2))+(1./(eps*(p**4))))
    eecs1=eecs1-(4.*(k**2))-(k/(eps))-(2.*(p**2)/((p**2)-(eps*k)))
    eecs1=eecs1*(dsqrt(((p**2)-(eps*k))/((eps**2)-(eps*k))))

    eecs2=((eps**3)*(k**2))-((k*((eps**2)+(p**2)))*3./(4.))&
    -((p**2)/(4.))*((eps+k)/((p**2)-(eps*k)))
    eecs2=eecs2*(1./((p**2)-(eps*k)))
    eecs2=eecs2+(2.*(eps**3))-(4.*eps)+k+(3./(4.*eps))
    eecs2=eecs2*(1./(p))*dlog(((2.*(p**2))-(k*(eps-p)))/((2.*(p**2))-(k*(eps+p))))
    eecs2=eecs2+(1./(p))*dlog(eps+p)*((16.*eps*(p**2))-(10.*(p**2)*k)-(4.*eps)-k&
    +(11./(2*eps))+(((4.*k)-eps)/(p**2))+(k/(p**4))+(k*((3*((eps**2)+(p**2)))&
    -(4.*(eps**3)*k))/(2.*((p**2)-(eps*k))))+(((p**2)/2.)*(eps+k)/(((p**2)&
    -(eps*k))**2))-(2.*(((eps**2)+(p**2))**3)/((eps**2)*k)))
    eecs2=eecs2+(((eps**3)*(eps-k))*32./(3))+(8.*(eps**2)*(k**2))+(k**2)&
    -(14.*(eps**2))-((p**2)*14./(3))+(eps*k*19./(3))-(k/eps)-((eps-k)*k/(p**2))&
    -((k)/(2*eps*(p**4)))-(((eps**2)+(k**2))/((p**2)-(eps*k)))
    eecs2=L3*eecs2*(1./(eps*p))

    eecs=(2.*(eps-k)/((p**2)-(eps*k)))*dlog(((2.*(eps*p))+(k*(eps-p)))&
    /((2.*(eps*p))-(k*(eps+p))))
    eecs=(-1.)*eecs*(((eps**2)*(k**2))/((2*p)*dsqrt(((p**2)-(eps*k))*((eps**2)-(eps*k)))))
    eecs=eecs+((k*dsqrt((p**2)-(eps*k)))/(2.*p))*(dlog(((2.*p*(dsqrt((eps**2)-(eps*k))))+k)&
    /((2.*p*(dsqrt((eps**2)-(eps*k))))-k)))*((eps/((eps*(p**2))+k))-(((eps**2)+(p**2))&
    /((p**2)-(eps*k)))-((2.*eps*k)/(((p**2)-(eps*k))**2)))
    eecs=eecs+((18.-(2.*k/(eps))+(1./(eps**2)))*L1)
    eecs=eecs+(((2.*L1*L3)/(p*(dsqrt(((p**2)-(eps*k))*((eps**2)-(eps*k))))))*((6.*(p**2)*k)&
    -(16.*eps*(p**2))-(eps*(k**2))-(2.*eps)-(3./eps)+((((eps**2)+(p**2))**3)/((eps**2)*k))))

    eecs=eecs+eecs1+eecs2+eecs3

    eecs=eecs*(afs*(r0**2))/(eps*p*k)
    
  end function eecs

  real(8) function jj(eps)
    implicit none
    integer i,n
    integer, parameter :: limit=10000
    real (kind=8) jiter, it1 !it1 iteration value
    real (kind=8) x, eps, p
    p=dsqrt((eps**2)-1.)
    it1=1
    n=1
    jiter=0
    jj=0
    do n=1,1000!while (it1.ge.1e-5)
      jj=jj+jiter
      jiter=log(eps+p)*((((eps+p)/(2.*eps))**n)+(((eps-p)/(2.*eps))**n))
      jiter=jiter-(1./n)*((((eps+p)/(2.*eps))**n)-(((eps-p)/(2.*eps))**n))
      jiter=(1./p)*(1./n)*jiter
      it1=jiter
      ! n=n+1
      print *, 'n',n,'jj', jj
    end do
    print *, '***', n
    return
  end

  pure function integrate(x,y) result(r)
      !! Calculates the integral of an array y with respect to x using the trapezoid
      !! approximation. Note that the mesh spacing of x does not have to be uniform.
      !! Source: https://fortranwiki.org/fortran/show/integration
      !! Created on April 22, 2016 17:31:53 by jabirali (46.9.153.214) (2538 characters / 1.0 pages) 
      integer, parameter :: wp=8
      real(wp), intent(in)  :: x(:)         !! Variable x
      real(wp), intent(in)  :: y(size(x))   !! Function y(x)
      real(wp)  :: r            !! Integral (y(x)*dx)
      ! Integrate using the trapezoidal rule
      associate(n => size(x))
        r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
      end associate
  end function integrate
  
end module cseqs
