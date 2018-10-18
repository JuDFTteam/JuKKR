!------------------------------------------------------------------------------------
!> Summary: PBE exchange correlation functional
!> Author: M. Ogura, K. Burke, E. Engel 
!> Exchange correlation potential making use of GGA, the parametrization is given 
!> by the PBE functional
!------------------------------------------------------------------------------------
module mod_mkxcpe2
  use :: mod_datatypes, only: dp
  use :: constants, only : pi
  private :: dp

contains

  !-------------------------------------------------------------------------------   
  !> Summary: PBE exchange correlation functional
  !> Author: M. Ogura
  !> Category: xc-potential, potential, KKRhost
  !> Deprecated: False
  !> Exchange correlation potential making use of GGA, the parametrization is given  
  !> by the PBE functional
  !-------------------------------------------------------------------------------   
  subroutine mkxcpe2(ir,np,rv,rholm,vxcp,excp,ylm,dylmt1,dylmf1,dylmf2,dylmtf,drrl, &
    ddrrl,drrul,ddrrul,irmd,lmpotd,lmmax,use_sol)
    ! ------------------------------------------------------------------
    ! Calculation of the exchange-correlation potential.
    ! coded by M. Ogura, Apr. 2015, Munich
    ! ------------------------------------------------------------------
    implicit none
    integer :: ijd
    parameter (ijd=434)
    ! .. Input variables
    integer, intent(in) :: ir
    integer, intent(in) :: np 
    integer, intent(in) :: irmd    !! Maximum number of radial points                
    integer, intent(in) :: lmmax   !! (LMAX+1)^2
    integer, intent(in) :: lmpotd  !! (lpot+1)**2 
    real (kind=dp), intent(in) :: rv
    logical, intent(in) :: use_sol !! use_sol=0 -> PBE, use_sol=1 -> PBEsol
    real (kind=dp), dimension(ijd, lmpotd), intent(in)  :: ylm !! real spherical harmonic to a given l,m
    real (kind=dp), dimension(irmd, lmpotd), intent(in) :: drrl                      
    real (kind=dp), dimension(lmpotd, 2), intent(in)    :: rholm !! l,m decomposed charge density
    real (kind=dp), dimension(irmd, lmpotd), intent(in) :: drrul                     
    real (kind=dp), dimension(irmd, lmpotd), intent(in) :: ddrrl                     
    real (kind=dp), dimension(irmd, lmpotd), intent(in) :: ddrrul                    
    real (kind=dp), dimension(ijd, lmpotd), intent(in)  :: dylmf1                    
    real (kind=dp), dimension(ijd, lmpotd), intent(in)  :: dylmf2                    
    real (kind=dp), dimension(ijd, lmpotd), intent(in)  :: dylmt1                    
    real (kind=dp), dimension(ijd, lmpotd), intent(in)  :: dylmtf    
    ! .. In/Out variables                                                            
    real (kind=dp), dimension(ijd), intent(inout) :: excp !! XC-energy               
    real (kind=dp), dimension(ijd,2), intent(inout) :: vxcp !! XC-potential   
    ! .. Local variables
    integer :: ip, ispin, l1, lm, n
    real (kind=dp) :: c, s, um, bet
    real (kind=dp), dimension(2) :: d, dl
    real (kind=dp), dimension(3,2) :: d1
    real (kind=dp), dimension(5,2) :: d2

    ! Set 'UM' for subroutine 'excpbex' and 'BET' for subroutine 'excpbec' for
    ! PBE or PBEsol
    if (use_sol) then
      um = 0.123456790123456e0_dp
      bet = 0.046e0_dp
    else
      um = 0.2195149727645171e0_dp
      bet = 0.06672455060314922e0_dp
    end if

    ! --- surface integration
    do ip = 1, np
      do ispin = 1, 2
        d(ispin) = 0e0_dp
        dl(ispin) = 0e0_dp
        do n = 1, 3
          d1(n, ispin) = 0e0_dp
        end do
        do n = 1, 5
          d2(n, ispin) = 0e0_dp
        end do
      end do
      do lm = 1, lmmax
        l1 = nint(sqrt(real(lm,kind=dp)-5e-1_dp))
        d(1) = d(1) + rholm(lm, 1)*ylm(ip, lm)
        d(2) = d(2) + rholm(lm, 2)*ylm(ip, lm)
        dl(1) = dl(1) + real(l1*(l1+1), kind=dp)*rholm(lm, 1)*ylm(ip, lm)
        dl(2) = dl(2) + real(l1*(l1+1), kind=dp)*rholm(lm, 2)*ylm(ip, lm)
        d1(1, 2) = d1(1, 2) + drrul(ir, lm)*ylm(ip, lm)
        d1(1, 1) = d1(1, 1) + (drrl(ir,lm)-drrul(ir,lm))*ylm(ip, lm)
        d1(2, 1) = d1(2, 1) + rholm(lm, 1)*dylmt1(ip, lm)
        d1(2, 2) = d1(2, 2) + rholm(lm, 2)*dylmt1(ip, lm)
        d1(3, 1) = d1(3, 1) + rholm(lm, 1)*dylmf1(ip, lm)
        d1(3, 2) = d1(3, 2) + rholm(lm, 2)*dylmf1(ip, lm)
        d2(1, 2) = d2(1, 2) + ddrrul(ir, lm)*ylm(ip, lm)
        d2(1, 1) = d2(1, 1) + (ddrrl(ir,lm)-ddrrul(ir,lm))*ylm(ip, lm)
        d2(2, 2) = d2(2, 2) + drrul(ir, lm)*dylmt1(ip, lm)
        d2(2, 1) = d2(2, 1) + (drrl(ir,lm)-drrul(ir,lm))*dylmt1(ip, lm)
        d2(3, 2) = d2(3, 2) + drrul(ir, lm)*dylmf1(ip, lm)
        d2(3, 1) = d2(3, 1) + (drrl(ir,lm)-drrul(ir,lm))*dylmf1(ip, lm)
        d2(4, 1) = d2(4, 1) + rholm(lm, 1)*dylmtf(ip, lm)
        d2(4, 2) = d2(4, 2) + rholm(lm, 2)*dylmtf(ip, lm)
        d2(5, 1) = d2(5, 1) + rholm(lm, 1)*dylmf2(ip, lm)
        d2(5, 2) = d2(5, 2) + rholm(lm, 2)*dylmf2(ip, lm)
      end do
      c = ylm(ip, 3)*sqrt(4e0_dp*pi/3e0_dp)
      s = sqrt(1e0_dp-c**2)
      do ispin = 1, 2
        dl(ispin) = dl(ispin)/rv**2
        d1(2, ispin) = d1(2, ispin)/rv
        d2(2, ispin) = d2(2, ispin)/rv
        d2(4, ispin) = d2(4, ispin)/rv**2
        if (s>1e-8_dp) then
          d1(3, ispin) = d1(3, ispin)/rv/s
          d2(3, ispin) = d2(3, ispin)/rv/s
          d2(5, ispin) = d2(5, ispin)/rv**2/s
          call fpexcpbe(d, dl, d1, d2, rv, s, c, vxcp(ip,1), vxcp(ip,2), excp(ip), um, bet)
        else
          d1(3, ispin) = 0e0_dp
          d2(3, ispin) = 0e0_dp
          d2(5, ispin) = 0e0_dp
          vxcp(ip, 1) = 0e0_dp
          vxcp(ip, 2) = 0e0_dp
          excp(ip) = 0e0_dp
        end if
      end do
    end do

    return
  end subroutine mkxcpe2

  !-------------------------------------------------------------------------------    
  !> Summary: Driver routine for PBE GGA subroutines 
  !> Author: M. Ogura
  !> Category: xc-potential, potential, KKRhost
  !> Deprecated: False
  !> Driver routine for the PBE GGA subroutines. It is the place where the actual
  !> exchange and correlation potentials are calculated.
  !-------------------------------------------------------------------------------   
  subroutine fpexcpbe(ro, rol, ro1, ro2, xr, s, c, v1, v2, exc, um, bet)
    ! ----------------------------------------------------------------------
    ! driver routine for PBE GGA subroutines.
    ! based on excpbe.f in Munich code (version on 20 Dec 2009)
    ! coded by M. Ogura, Jun. 2011, Munich
    ! ----------------------------------------------------------------------
    implicit none

    ! .. Input variables
    real (kind=dp), intent(in) :: c
    real (kind=dp), intent(in) :: s
    real (kind=dp), intent(in) :: xr
    real (kind=dp), intent(in) :: um
    real (kind=dp), intent(in) :: bet
    real (kind=dp), dimension(2), intent(in) :: ro
    real (kind=dp), dimension(2), intent(in) :: rol
    real (kind=dp), dimension(3,2), intent(in) :: ro1
    real (kind=dp), dimension(5,2), intent(in) :: ro2
    ! .. Output variables
    real (kind=dp), intent(out) :: v1 !! xc-potential up
    real (kind=dp), intent(out) :: v2 !! xc-potential down
    real (kind=dp), intent(out) :: exc !! xc-energy
    ! .. Local variables
    integer :: jsp, llda
    real (kind=dp) :: conf, conrs, d, drv1, drv2, drv2s, drv3, drv4, ec, ex, fk, g 
    real (kind=dp) :: rs, sk, ss, thrd, thrd2, tt, uu, vcdn, vcup, vv, vx
    real (kind=dp) :: vxcdn, vxcup, ww, x, xd, xu, y, z, zet

    thrd  = 1e0_dp/3e0_dp
    thrd2 = 2e0_dp/3e0_dp
    conf  = (3e0_dp*pi**2)**thrd
    conrs = (3e0_dp/(4e0_dp*pi))**thrd
    llda  = 0
    exc   = 0e0_dp
    vxcup = 0e0_dp
    vxcdn = 0e0_dp
    if (ro(1)>1e-12_dp .and. ro(2)>1e-12_dp) then

      ! ---begin the spin loop for exchange
      if (ro(1)<=1e-6_dp .or. ro(2)<=1e-6_dp) llda = 1
      do jsp = 1, 2
        d  = 2e0_dp*ro(jsp)
        fk = conf*d**thrd
        x  = ro1(1, jsp)
        y  = ro1(2, jsp)
        z  = ro1(3, jsp)
        drv1 = sqrt(x**2+y**2+z**2)*2e0_dp
        if (abs(drv1)<1e-8_dp) then
          drv2 = 0e0_dp
        else
          drv2s = 2e0_dp*y*z*ro2(4, jsp) + (z**2-y**2)*ro2(5, jsp) - c/xr*y*(z**2+y**2)
          drv2 = x**2*ro2(1, jsp) + 2e0_dp*x*y*ro2(2, jsp) + 2e0_dp*x*z*ro2(3, jsp) - x*y**2/xr - x*z**2/xr - y**2*rol(jsp)
          if (abs(drv2s)>=1e-10_dp) drv2 = drv2 + drv2s/s
          drv2 = drv2/drv1*8e0_dp
        end if
        drv3 = ro2(1, jsp) + 2e0_dp/xr*x - rol(jsp)
        drv3 = drv3*2e0_dp
        ss = drv1/(d*2e0_dp*fk)
        uu = drv2/(d**2*(2e0_dp*fk)**3)
        vv = drv3/(d*(2e0_dp*fk)**2)
        call excpbex(d, ss, uu, vv, ex, vx, llda, um)
        exc = exc + ex*(d/2e0_dp)/(ro(1)+ro(2))
        if (jsp==1) vxcup = vx
        if (jsp==2) vxcdn = vx
      end do

      ! ---correlation
      d = ro(1) + ro(2)
      zet = (ro(1)-ro(2))/d
      rs = conrs/d**thrd
      fk = 1.91915829e0_dp/rs
      sk = sqrt(4e0_dp*fk/pi)
      g = ((1e0_dp+zet)**thrd2+(1e0_dp-zet)**thrd2)/2e0_dp
      x = ro1(1, 1) + ro1(1, 2)
      y = ro1(2, 1) + ro1(2, 2)
      z = ro1(3, 1) + ro1(3, 2)
      drv1 = sqrt(x**2+y**2+z**2)
      if (drv1<1e-8_dp) then
        drv2 = 0e0_dp
      else
        drv2s = 2e0_dp*y*z*(ro2(4,1)+ro2(4,2)) + (z**2-y**2)*(ro2(5,1)+ro2(5,2)) - c/xr*y*(z**2+y**2)
        drv2 = x**2*(ro2(1,1)+ro2(1,2)) + 2e0_dp*x*y*(ro2(2,1)+ro2(2,2)) + 2e0_dp*x*z*(ro2(3,1)+ro2(3,2)) - x*y**2/xr - x*z**2/xr - y**2*(rol(1)+rol(2))
        if (abs(drv2s)>=1e-10_dp) drv2 = drv2 + drv2s/s
        drv2 = drv2/drv1
      end if
      drv3 = ro2(1, 1) + ro2(1, 2) + 2e0_dp/xr*x - rol(1) - rol(2)
      drv4 = x*(ro1(1,1)-ro1(1,2)-zet*x) + y*(ro1(2,1)-ro1(2,2)-zet*y) + z*(ro1(3,1)-ro1(3,2)-zet*z)
      tt = drv1/(d*2e0_dp*sk*g)
      uu = drv2/(d**2*(2e0_dp*sk*g)**3)
      vv = drv3/(d*(2e0_dp*sk*g)**2)
      ww = drv4/(d**2*(2e0_dp*sk*g)**2)
      call excpbec(rs, zet, tt, uu, vv, ww, ec, vcup, vcdn, llda, bet)
      exc = exc + ec
      vxcup = vxcup + vcup
      vxcdn = vxcdn + vcdn

    end if
    ! ---convert from h to ry
    exc = 2e0_dp*exc
    xu  = 2e0_dp*vxcup
    xd  = 2e0_dp*vxcdn

    v1 = xu
    v2 = xd
  end subroutine fpexcpbe

  !-------------------------------------------------------------------------------    
  !> Summary: PBE exchange for a spin-**unpolarized** electronic system
  !> Author: K. Burke, E. Engel
  !> Category: xc-potential, potential, KKRhost
  !> Deprecated: False
  !> PBE exchange for a spin-**unpolarized** electronic system. This is calculated
  !> using the fact that the exchange in LDA is given by
  !> \begin{equation}
  !> e_x[unif]=a_x\rho^{\frac{4}{3}}
  !> \end{equation}
  !> with
  !> \begin{equation}
  !> a_x = -\frac{3}{4}\left(\frac{3}{\pi}\right)^{\frac{1}{3}}
  !> \end{equation}
  !> where the \(e_x[PBE]\) is given by 
  !> \begin{equation}
  !> e_x[PBE] =e_x[unif]*F_x^{PBE}(s)
  !> \end{equation}
  !> where
  !> \begin{equation}
  !> F_x^{PBE}(s)=1+u_k-\frac{u_k}{1+u_l s^2}
  !> \end{equation}
  !> with \(u_k\) and \(u_l\) being defined in Eq.13 of [a]
  !> @note
  !> [a] J.P. Perdew, K. Burke, and M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996).
  !> [b] J.P. Perdew and Y. Wang, Phys. Rev. B33, 8800 (1986); B40, 3399 (1989)(E).
  !> 
  !> K Burke's modification of PW91 codes, May 14, 1996
  !> Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
  !> 
  !> **All input and output is in atomic units**
  !> @endnote
  !-------------------------------------------------------------------------------   
  subroutine excpbex(rho, s, u, v, ex, vx, llda, um)
    ! ----------------------------------------------------------------------
    ! PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
    ! K Burke's modification of PW91 codes, May 14, 1996
    ! Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
    ! ----------------------------------------------------------------------
    ! INPUT rho : DENSITY
    ! INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
    ! INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
    ! INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
    ! (for U,V, see PW86(24))
    ! OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
    ! ----------------------------------------------------------------------

    implicit none

    ! PARAMETER definitions
    real (kind=dp), parameter :: thrd  = 1.e0_dp/3.e0_dp
    real (kind=dp), parameter :: thrd4 = 4.e0_dp/3.e0_dp
    real (kind=dp), parameter :: ax = -0.738558766382022405884230032680836e0_dp
    real (kind=dp), parameter :: uk = 0.8040e0_dp

    ! Dummy arguments
    integer, intent(in) :: llda
    real (kind=dp), intent(in) :: rho !! density
    real (kind=dp), intent(in) :: s !! \(\frac{\| \nabla \rho\|}{2K_F \rho}\) with \(\K_F=\left(3 \pi^2 \rho\right)^\frac{1}{3}\)
    real (kind=dp), intent(in) :: u !! \(\frac{\nabla^2 \rho}{\rho\left(2K_F\right)^2}\)
    real (kind=dp), intent(in) :: v
    real (kind=dp), intent(in) :: um
    real (kind=dp), intent(out) :: vx !! exchange potential
    real (kind=dp), intent(out) :: ex !! exchange energy per electron

    ! Local variables
    real (kind=dp) :: exunif, fs, fss, fxpbe, p0, s2
    real (kind=dp) :: ul

    ! ----------------------------------------------------------------------
    ! Define UL with via UM and UK
    ul = um/uk

    ! ----------------------------------------------------------------------
    ! construct LDA exchange energy density
    exunif = ax*rho**thrd
    if (llda==1) then
      ex = exunif
      vx = ex*thrd4
      return
    end if
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    ! construct PBE enhancement factor
    s2 = s*s
    p0 = 1.e0_dp + ul*s2
    fxpbe = 1e0_dp + uk - uk/p0
    ex = exunif*fxpbe
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    ! ENERGY DONE. NOW THE POTENTIAL:
    ! find first and second derivatives of Fx w.r.t s.
    ! Fs=(1/s)*d FxPBE/ ds
    ! Fss=d Fs/ds
    fs = 2.e0_dp*uk*ul/(p0*p0)
    fss = -4.e0_dp*ul*s*fs/p0
    ! ----------------------------------------------------------------------
    ! ----------------------------------------------------------------------
    ! calculate potential from [b](24)
    vx = exunif*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
  end subroutine excpbex

  !-------------------------------------------------------------------------------    
  !> Summary: This subroutine evaluates the correlation energy per particle for PBE xc-potential
  !> Author: K. Burke, E. Engel
  !> Category: xc-potential, potential, KKRhost
  !> Deprecated: False
  !> This subroutine evaluates the correlation energy per particle and
  !> spin-up and spin-dn correlation potentials within the Perdew-Burke-
  !> Ernzerhof GGA. It is a slightly modified version of K. Burke's
  !> official PBE subroutine.
  !> @note
  !> [a] J.P. Perdew, K. Burke, and M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996).
  !> [b] J. P. Perdew, K. Burke, and Y. Wang, Phys. Rev. B54, 16533 (1996).
  !> [c] J. P. Perdew and Y. Wang, Phys. Rev. B45, 13244 (1992).
  !> 
  !> **All input and output is in atomic units**
  !> @endnote
  !-------------------------------------------------------------------------------   
  subroutine excpbec(rs, zeta, t, uu, vv, ww, ec, vcup, vcdn, llda, bet)

    implicit none

    ! PARAMETER definitions
    real (kind=dp), parameter :: thrd = 1.e0_dp/3.e0_dp
    real (kind=dp), parameter :: thrdm = -thrd
    real (kind=dp), parameter :: thrd2 = 2.e0_dp*thrd
    real (kind=dp), parameter :: sixthm = thrdm/2.e0_dp
    real (kind=dp), parameter :: thrd4 = 4.e0_dp*thrd
    real (kind=dp), parameter :: gam = 0.5198420997897463295344212145565e0_dp
    real (kind=dp), parameter :: fzz = 8.e0_dp/(9.e0_dp*gam)
    real (kind=dp), parameter :: gamma = 0.03109069086965489503494086371273e0_dp
    real (kind=dp), parameter :: eta = 1.e-12_dp

    ! Dummy arguments
    real (kind=dp), intent(in) :: rs !! Wigner-Seitz radius \(\left(\frac{3}{4\pi\left(\rho_{up}+\rho_{down}\right)}\right)^\frac{1}{3}\)
    real (kind=dp), intent(in) :: t  !! \(\frac{\| \nabla D\|}{2S_K G D} \)
    real (kind=dp), intent(in) :: uu !! \(\frac{\nabla D \nabla\|\nabla D \|}{\left(2S_K G\right)^3 D^2}\)
    real (kind=dp), intent(in) :: vv !! \(\frac{\nabla^2 D}{\left(2S_K G\right)^2 D}\)
    real (kind=dp), intent(in) :: ww !! \(\frac{\nabla D \nabla Z}{\left(2S_K G\right)^2 D}\)
    real (kind=dp), intent(in) :: zeta !! \(Z=\frac{\rho_{up}-\rho_{down}}{\rho_{up}+\rho_{down}}\) relative spin polarization
    real (kind=dp), intent(in) :: bet !! coefficient in gradient expansion for correlation, [a](4).
    integer, intent(in) :: llda
    real (kind=dp), intent(out) :: ec !! correlation energy per particle
    real (kind=dp), intent(out) :: vcdn !! spin-up correlation potential
    real (kind=dp), intent(out) :: vcup !! spin-dn correlation potential

    !! Local Feri momentum  \(F_K = \left(3\pi^2\left(\rho_{up}+\rho_{down}\right)\right)^\frac{1}{3}\)
    !! Local screening momentum  \(S_K = \lef(\frac{4*F_K}{\pi}\right)^\frac{1}{2}\)
    ! Local variables
    real (kind=dp) :: alfm, alfrsm, b, b2, bec, bg, comm, ecrs, eczeta, ep, eprs, eu
    real (kind=dp) :: eurs, f, fac, fact0, fact1, fact2, fact3, fact5, fz, g, g3, g4
    real (kind=dp) :: gz, h, hb, hbt, hrs, hrst, ht, htt, hz, hzt, pon, pref, q4, q5
    real (kind=dp) :: q8, q9, rsthrd, rtrs, t2, t4, t6, z4, delt

    ! thrd*=various multiples of 1/3
    ! numbers for use in LSD energy spin-interpolation formula, [c](9).
    ! GAM= 2^(4/3)-2
    ! FZZ=f''(0)= 8/(9*GAM)
    ! numbers for construction of PBE
    ! gamma=(1-log(2))/pi^2
    ! bet=coefficient in gradient expansion for correlation, [a](4).
    ! eta=small number to stop d phi/ dzeta from blowing up at
    ! |zeta|=1.
    ! ----------------------------------------------------------------------
    ! Define DELT via BET and GAMMA
    delt = bet/gamma

    ! ----------------------------------------------------------------------
    ! find LSD energy contributions, using [c](10) and Table I[c].
    ! EU=unpolarized LSD correlation energy
    ! EURS=dEU/drs
    ! EP=fully polarized LSD correlation energy
    ! EPRS=dEP/drs
    ! ALFM=-spin stiffness, [c](3).
    ! ALFRSM=-dalpha/drs
    ! F=spin-scaling factor from [c](9).
    ! construct ec, using [c](8)
    if (rs<3.e5_dp) then
      rtrs = sqrt(rs)
      call excgcor2(0.0310907e0_dp, 0.21370e0_dp, 7.5957e0_dp, 3.5876e0_dp, 1.6382e0_dp, 0.49294e0_dp, rtrs, eu, eurs)
      call excgcor2(0.01554535e0_dp, 0.20548e0_dp, 14.1189e0_dp, 6.1977e0_dp, 3.3662e0_dp, 0.62517e0_dp, rtrs, ep, eprs)
      call excgcor2(0.0168869e0_dp, 0.11125e0_dp, 10.357e0_dp, 3.6231e0_dp, 0.88026e0_dp, 0.49671e0_dp, rtrs, alfm, alfrsm)
      z4 = zeta**4
      f = ((1.e0_dp+zeta)**thrd4+(1.e0_dp-zeta)**thrd4-2.e0_dp)/gam
      ec = eu*(1.e0_dp-f*z4) + ep*f*z4 - alfm*f*(1.e0_dp-z4)/fzz
      ! ----------------------------------------------------------------------
      ! ----------------------------------------------------------------------
      ! LSD potential from [c](A1)
      ! ECRS = dEc/drs [c](A2)
      ! ECZETA=dEc/dzeta [c](A3)
      ! FZ = dF/dzeta [c](A4)
      ecrs = eurs*(1.e0_dp-f*z4) + eprs*f*z4 - alfrsm*f*(1.e0_dp-z4)/fzz
      fz = thrd4*((1.e0_dp+zeta)**thrd-(1.e0_dp-zeta)**thrd)/gam
      eczeta = 4.e0_dp*(zeta**3)*f*(ep-eu+alfm/fzz) + fz*(z4*ep-z4*eu-(1.e0_dp-z4)*alfm/fzz)
      comm = ec - rs*ecrs/3.e0_dp - zeta*eczeta
      vcup = comm + eczeta
      vcdn = comm - eczeta
      if (llda==1) return
      ! ----------------------------------------------------------------------
      ! ----------------------------------------------------------------------
      ! PBE correlation energy
      ! G=phi(zeta), given after [a](3)
      ! DELT=bet/gamma
      ! B=A of [a](8)
      g = ((1.e0_dp+zeta)**thrd2+(1.e0_dp-zeta)**thrd2)/2.e0_dp
      g3 = g**3
      pon = -ec/(g3*gamma)
      b = delt/(exp(pon)-1.e0_dp)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      q4 = 1.e0_dp + b*t2
      q5 = 1.e0_dp + b*t2 + b2*t4
      h = g3*(bet/delt)*log(1.e0_dp+delt*q4*t2/q5)
      ec = ec + h
      ! ----------------------------------------------------------------------
      ! ----------------------------------------------------------------------
      ! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      g4 = g3*g
      t6 = t4*t2
      rsthrd = rs/3.e0_dp
      gz = (((1.e0_dp+zeta)**2+eta)**sixthm-((1.e0_dp-zeta)**2+eta)**sixthm)/3.e0_dp
      fac = delt/b + 1.e0_dp
      bg = -3.e0_dp*b2*ec*fac/(bet*g4)
      bec = b2*fac/(bet*g3)
      q8 = q5*q5 + delt*q4*q5*t2
      q9 = 1.e0_dp + 2.e0_dp*b*t2
      hb = -bet*g3*b*t6*(2.e0_dp+b*t2)/q8
      hrs = -rsthrd*hb*bec*ecrs
      fact0 = 2.e0_dp*delt - 6.e0_dp*b
      fact1 = q5*q9 + q4*q9*q9
      hbt = 2.e0_dp*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
      hrst = rsthrd*t2*hbt*bec*ecrs
      hz = 3.e0_dp*gz*h/g + hb*(bg*gz+bec*eczeta)
      ht = 2.e0_dp*bet*g3*q9/q8
      hzt = 3.e0_dp*gz*ht/g + hbt*(bg*gz+bec*eczeta)
      fact2 = q4*q5 + b*t2*(q4*q9+q5)
      fact3 = 2.e0_dp*b*q5*q9 + delt*fact2
      htt = 4.e0_dp*bet*g3*t*(2.e0_dp*b/q8-(q9*fact3/q8)/q8)
      comm = h + hrs + hrst + t2*ht/6.e0_dp + 7.e0_dp*t2*t*htt/6.e0_dp
      pref = hz - gz*t2*ht/g
      fact5 = gz*(2.e0_dp*ht+t*htt)/g
      comm = comm - pref*zeta - uu*htt - vv*ht - ww*(hzt-fact5)
      vcup = vcup + comm + pref
      vcdn = vcdn + comm - pref
    else
      vcup = 0.e0_dp
      vcdn = 0.e0_dp
    end if
  end subroutine excpbec

  !-------------------------------------------------------------------------------    
  !> Summary: Slimmed down version of GCOR used in PW91 routines, to interpolate LSD correlation energy
  !> Author: K. Burke
  !> Category: xc-potential, potential, KKRhost
  !> Deprecated: False
  !> Slimmed down version of GCOR used in PW91 routines, to interpolate LSD
  !> correlation energy as given by Eq. 10 of [a] 
  !> @note
  !> [a] J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992).
  !> @endnote
  !-------------------------------------------------------------------------------   
  subroutine excgcor2(a, a1, b1, b2, b3, b4, rtrs, gg, ggrs)

    implicit none

    ! Dummy arguments
    real (kind=dp) :: a, a1, b1, b2, b3, b4, gg, ggrs, rtrs

    ! Local variables
    real (kind=dp) :: q0, q1, q2, q3

    q0 = -2.e0_dp*a*(1.e0_dp+a1*rtrs*rtrs)
    q1 = 2.e0_dp*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
    q2 = log(1.e0_dp+1.e0_dp/q1)
    gg = q0*q2
    q3 = a*(b1/rtrs+2.e0_dp*b2+rtrs*(3.e0_dp*b3+4.e0_dp*b4*rtrs))
    ggrs = -2.e0_dp*a*a1*q2 - q0*q3/(q1*(1.e0_dp+q1))
  end subroutine excgcor2

end module mod_mkxcpe2
