module SingleSiteHelpers_mod
  implicit none
  private
  public :: beshan, bessel
  public :: regsol, irwsol
  public :: zgeinv1, vllns
  public :: wftsca, wfmesh
  public :: wfint0, wfint
  public :: csinwd, csout
  public :: cradwf

  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0), ci=(0.d0,1.d0)
  double precision, parameter :: pi=4.d0*atan(1.d0)

  contains

  subroutine beshan(hl, jl, nl, z, lmax)
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders lmin  <=  l  <=  lmax.
!  for |z|  <   1 the taylor expansions of jl and nl are used.
!  for |z|  >=  1 the explicit expressions for hl(+), hl(-) are used.
!
!                            r. zeller   jan. 1990
!-----------------------------------------------------------------------
    double complex, intent(in) :: z
    integer, intent(in) :: lmax
    double complex, intent(out) :: hl(0:lmax), jl(0:lmax), nl(0:lmax)
    
    double complex :: termj, termn, z2, zj, zn
    double precision :: rl, rn, rnm
    integer :: l, m, n

    zj = 1.d0
    zn = 1.d0
    z2 = z*z
    if (abs(z) < lmax+1.d0) then
      do l = 0, lmax
        rl = l + l
        termj = -0.5d0/(rl+3.d0)*z2
        termn =  0.5d0/(rl-1.d0)*z2
        jl(l) = 1.d0
        nl(l) = 1.d0
        do n = 2, 25
          jl(l) = jl(l) + termj
          nl(l) = nl(l) + termn
          rn = n + n
          termj = -termj/(rl+rn+1.d0)/rn*z2
          termn =  termn/(rl-rn+1.d0)/rn*z2
        enddo ! n
        jl(l) = jl(l)*zj
        nl(l) = -nl(l)*zn/z
        hl(l) = jl(l) + nl(l)*ci

        zj = zj*z/(rl+3.d0)
        zn = zn/z*(rl+1.d0)
      enddo ! l
    endif

    do l = 0, lmax
      if (abs(z) >= l+1.d0) then
        hl(l) = 0.d0
        nl(l) = 0.d0
        rnm = 1.d0
        do m = 0, l
          hl(l) = hl(l) + rnm/(-ci*(z+z))**m
          nl(l) = nl(l) + rnm/( ci*(z+z))**m
          rnm = rnm*(l*l+l-m*m-m)/(m+1.d0)
        enddo ! m
        hl(l) = hl(l)*(-ci)**l*exp( ci*z)/( ci*z)
        nl(l) = nl(l)*( ci)**l*exp(-ci*z)/(-ci*z)
        jl(l) = (hl(l)+nl(l))*0.5d0
        nl(l) = (hl(l)-jl(l))/ci
      endif
    enddo ! l

  endsubroutine beshan
    
  subroutine bessel(jl, nl, hl, z, lmx)
!**********************************************************************
!    attention : contrary to abramowitz and stegun and
!                contrary to subroutine beshan
!
!                the bessel functions of third kind ( hankel functions)
!                are defined as:      hl(l) = nl(l) - i * jl(l)
!**********************************************************************
    double complex, intent(in) :: z
    integer, intent(in) :: lmx
    double complex, intent(out) :: hl(0:lmx), jl(0:lmx), nl(0:lmx)
    
    call beshan(hl, jl, nl, z, lmx)
    hl(0:lmx) = -ci*hl(0:lmx) ! scale with -i

  endsubroutine bessel
  
  
  subroutine cradwf(e,ek,nsra,alpha,ipan,ircut,cvlight,rs,s,pz,fz,qz,sz,tmat,vm2z,drdi,r,z,ldau,nldau,lldau,wmldauav,ldaucut, lmaxd, irmd, ipand)
!-----------------------------------------------------------------------
!  subroutine for radial wave functions of spherical potentials
!
!             the generalized phase shifts are calculated by
!             a wronski relation :
!
!                 alpha(z,l) =-sqrt(z)*wronski{hl(r;z),rl(r;z)}; r->0
!
!             where hl is the free hankel function and rl the regular
!             solution . using the analytical behaviour of rl at the
!             origin (rl = alphal * r**(l+1)  ; r->0),
!             the generalized phase shifts can be calculated
!             directly with the renormalization alphal .
!                                           b.drittler nov.1987
!-----------------------------------------------------------------------
    integer, intent(in) :: lmaxd, irmd, ipand
    double complex, intent(in) :: e, ek
    double precision, intent(in) :: cvlight, z
    integer, intent(in) :: ipan,nsra,nldau
    logical, intent(in) :: ldau
    double complex, intent(out) :: alpha(0:lmaxd)
    double complex, intent(out) :: fz(irmd,0:lmaxd), pz(irmd,0:lmaxd), qz(irmd,0:lmaxd), sz(irmd,0:lmaxd)
    double complex, intent(out) :: tmat(0:lmaxd)
    double precision, intent(in) :: drdi(irmd),r(irmd), rs(irmd,0:lmaxd),s(0:lmaxd), vm2z(irmd), ldaucut(irmd), wmldauav(lmaxd+1)
    integer, intent(in) :: ircut(0:ipand), lldau(lmaxd+1)
    
    double complex :: alphal,arg,bl,eklfac,hl,pn,qf,slope,tl,tlsqez, value,w,x,y
    double precision :: rirc,rsirc,s1
    integer :: ir,irc1,l,n, lmaxp1
    double complex :: bessjw(0:lmaxd+1),bessyw(0:lmaxd+1),dlogdp(0:lmaxd)
    double complex :: hamf(irmd,0:lmaxd),hankws(0:lmaxd+1),mass(irmd)
    double precision :: dror(irmd)

    lmaxp1 = lmaxd+1

! initialisations
    dror = 0.d0
    hamf = zero
    mass = zero
! initialisations

    irc1 = ircut(ipan)
    do ir = 2,irc1
      dror(ir) = drdi(ir)/r(ir)
    enddo
    rirc = r(irc1)
    arg = rirc*ek
    call beshan(hankws,bessjw,bessyw,arg,lmaxp1)
!
!    attention : contrary to abramowitz and stegun and
!                the bessel functions of third kind ( hankel functions)
!                are defined as:      hl(l) = nl(l) - i * jl(l)
!
    do l = 0,lmaxp1
      hankws(l) = bessyw(l) - ci*bessjw(l)
    enddo ! l
!
!---> calculate regular wavefunctions
!
    call regsol(cvlight,e,nsra,dlogdp,fz,hamf,mass,pz,dror,r,s,vm2z, z,ipan,ircut,irmd,ipand,lmaxd, ldau,nldau,lldau,wmldauav,ldaucut)
!
    eklfac = ek
!
    do l = 0, lmaxd
      s1 = s(l)
      rsirc = rs(irc1,l)
      eklfac = eklfac/ek*dble(2*l+1)
!
!---> determine t - matrix
!
      pn = pz(irc1,l)*rsirc
      n = l+1
      qf = dble(l)/rirc
      hl = hankws(l)
      bl = bessjw(l)
      x = qf*hl - ek*hankws(n)
      y = qf*bl - ek*bessjw(n)
      w = dlogdp(l)
      tlsqez = (bl*w - y)/(x - hl*w)
      tl = tlsqez/ek
      tmat(l) = tl
!
!---> determine the renormalization
!
      alphal = (bl + hl*tlsqez)*rirc/pn
!
!---> determine the alpha matrix
!
      alpha(l) = alphal*eklfac

      pz(2:irc1,l) = pz(2:irc1,l)*alphal
      fz(2:irc1,l) = fz(2:irc1,l)*alphal
!
      value = hl*rirc*rsirc
      slope = dble(l+1)*hl - rirc*ek*hankws(l+1)
      slope = (slope*rsirc + s1/rirc*value)
      qz(irc1,l) = value
      sz(irc1,l) = (slope*rirc - (s1 + 1.d0)*value)/mass(irc1)*dror(irc1)
    enddo ! l
!
!---> calculate irregular wavefunctions
!
    call irwsol(ek,fz,hamf,mass,pz,qz,sz,dror,s,ipan,ircut,irmd, ipand,lmaxd)

    do l = 0, lmaxd
      if (nsra == 2) then
        pz(2:irc1,l) = pz(2:irc1,l)*rs(2:irc1,l)
        qz(2:irc1,l) = qz(2:irc1,l)/rs(2:irc1,l)
        fz(2:irc1,l) = fz(2:irc1,l)*rs(2:irc1,l)/cvlight
        sz(2:irc1,l) = sz(2:irc1,l)/rs(2:irc1,l)/cvlight
      else
        pz(2:irc1,l) = pz(2:irc1,l)*rs(2:irc1,l)
        qz(2:irc1,l) = qz(2:irc1,l)/rs(2:irc1,l)
        fz(2:irc1,l) = zero
        sz(2:irc1,l) = zero
      endif ! nsra
    enddo ! l 

  endsubroutine cradwf

  
    
      
  subroutine csinwd(f, fint, lmmsqd, irmind, irmd, ipan, ircut)
!-----------------------------------------------------------------------
!     this subroutine does an inwards integration of llmax
!     functions f with an extended 3-point-simpson :
!
!
!                               irmax
!                   fint(:,i) = { f(:,i') di'
!                                ir
!
!  the starting value for this integration at is - 1 is determined by
!    a 4 point lagrangian integration  , coefficients given by
!    m. abramowitz and i.a. stegun, handbook of mathematical functions,
!    nbs applied mathematics series 55 (1968)
!
!  attention in case of radial integration :
!       the weights drdi have to be multiplied before calling this
!       subroutine .
!
!                                     b. drittler mar. 1989
!
!    modified for functions with kinks - at each kink the integration
!      is restarted
!
!    attention : it is supposed that irmin + 3 is less than imt !
!
!
!                                     b. drittler july 1989
!    modified by m. ogura, june 2015
!-----------------------------------------------------------------------
    
    integer, intent(in) :: ipan, irmd, irmind, lmmsqd
    double complex, intent(in) :: f(lmmsqd,irmind:irmd)
    double complex, intent(out) :: fint(lmmsqd,irmind:irmd)
    integer, intent(in) :: ircut(0:ipan)
    
    double precision, parameter :: a1=5.d0/12.d0, a2=8.d0/12.d0, a3=-1.d0/12.d0
    integer :: ir, ie, ip, is
!
!---> loop over kinks
!
    do ip = ipan, 1, -1
      is = ircut(ip)
      ie = ircut(ip-1) + 1
      if (ip == 1) ie = irmind ! first panel
!
      if (ip == ipan) then
        fint(:,is) = 0.d0 ! last panel
      else
        fint(:,is) = fint(:,is+1)
      endif ! last panel
!
!---> calculate fint with an extended 3-point-simpson
!
      do ir = is, ie+2, -2
        fint(:,ir-1) = fint(:,ir-0) + f(:,ir)*a1 + f(:,ir-1)*a2 + f(:,ir-2)*a3
        fint(:,ir-2) = fint(:,ir-1) + f(:,ir)*a3 + f(:,ir-1)*a2 + f(:,ir-2)*a1
      enddo ! ir
      if (mod(is-ie,2) == 1) then
        fint(:,ie  ) = fint(:,ie+1) + f(:,ie)*a1 + f(:,ie+1)*a2 + f(:,ie+2)*a3
      endif
    enddo ! ip
    
  endsubroutine ! csinwd
    
      
      
  subroutine csout(f,fint,lmmsqd,irmind,irmd,ipan,ircut)
!-----------------------------------------------------------------------
!     this subroutine does an outwards integration of llmax
!     functions f with an extended 3-point-simpson :
!
!
!                                ir
!                   fint(ll,i) = { f(ll,i') di'
!                               irmin
!
!  the starting value for this integration at irmin+1 is determined by
!    a 4 point lagrangian integration  , coefficients given by
!    m. abramowitz and i.a. stegun, handbook of mathematical functions,
!    nbs applied mathematics series 55 (1968)
!
!  attention in case of radial integration :
!       the weights drdi have to be multiplied before calling this
!       subroutine .
!
!                                     b. drittler mar. 1989
!
!    modified for functions with kinks - at each kink the integration
!      is restarted
!
!    attention : it is supposed that irmin + 3 is less than imt !
!
!
!                                     b. drittler july 1989
!    modified by m. ogura, june 2015
!-----------------------------------------------------------------------
    integer, intent(in) :: ipan, irmd, irmind, lmmsqd
    double complex, intent(in) ::  f(lmmsqd,irmind:irmd)
    double complex, intent(inout) :: fint(lmmsqd,irmind:irmd)
    integer, intent(in) :: ircut(0:ipan)

    double precision, parameter :: a1=5.d0/12.d0, a2=8.d0/12.d0, a3=-1.d0/12.d0
    integer :: i, ien, ip, ist
!
!---> loop over kinks
!
    do ip = 1, ipan
      ien = ircut(ip)
      ist = ircut(ip-1) + 1

      if (ip == 1) then
        ist = irmind
        fint(:,ist) = 0.d0
      else
        fint(:,ist) = fint(:,ist-1)
      endif
!
!---> calculate fint with an extended 3-point-simpson
!
      do i = ist, ien-2, 2
        fint(:,i+1) = fint(:,i+0) + f(:,i)*a1 + f(:,i+1)*a2 + f(:,i+2)*a3
        fint(:,i+2) = fint(:,i+1) + f(:,i)*a3 + f(:,i+1)*a2 + f(:,i+2)*a1
      enddo ! i
        
      if (mod(ien-ist, 2) == 1) then
        fint(:,ien) = fint(:,ien-1) + f(:,ien)*a1 + f(:,ien-1)*a2 + f(:,ien-2)*a3
      endif
    enddo ! ip

  endsubroutine ! csout
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine zgeinv1(a,u,aux,ipiv,n)
! ------------------------------------------------------------------------
!   - inverts a general double complex matrix a,
!   - the result is return in u,
!   - input matrix a is returned unchanged,
!   - aux is a auxiliary matrix,
!   - a,u and aux are of dimension (n,n),
! ------------------------------------------------------------------------
    integer, intent(in) :: n
    double complex, intent(in) :: a(n,n)
    double complex, intent(out) :: u(n,n)
    double complex, intent(inout) :: aux(n,n)
    integer, intent(inout) :: ipiv(:)
    
    integer :: i, info
    external :: zcopy, zgetrs, zgetrf ! from BLAS

    u(:,1:n) = zero
    do i = 1, n
      u(i,i) = cone
    enddo ! i

    call zcopy(n*n,a,1,aux,1)
    call zgetrf(n,n,aux,n,ipiv,info)
    call zgetrs('n',n,n,aux,n,ipiv,u,n,info)
    
  endsubroutine
      
      
#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

  subroutine regsol(cvlight, e, nsra, dlogdp, fz, hamf, mass, pz, dror, r, &
                 s, vm2z, z, ipan, ircut, irmd, ipand, lmaxd, &
                 ldau, nldau, lldau, wmldauav, ldaucut)
  !-----------------------------------------------------------------------
  !  calculates the regular solution of the schroedinger equation or
  !    in semi relativistic approximation for a spherically averaged
  !    potential and given energy . to archieve greater presion the
  !    leading power r**s ( in schroedinger case s = l , in case of sra
  !    s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
  !    from the wavefunction .
  !
  !  the t - matrix has to be determined at the mt radius in case of
  !    a mt calculation or at the ws radius in case of a ws calcu-
  !    lation . therefore the logarithmic derivative is calculated
  !    at that point (=ircut(ipan) )
  !
  !  the differential equation is solved with a 5 point adams - bashforth
  !    and adams - moulton predictor corrector method integrating
  !    outwards and extended for potentials with kinks
  !
  !                                               b.drittler   nov 1989
  !-----------------------------------------------------------------------
    double precision, intent(in) :: cvlight
    double complex, intent(in) :: e
    integer, intent(in) :: nsra
    double complex, intent(out) :: dlogdp(0:lmaxd)
    double complex, intent(out) :: hamf(irmd,0:lmaxd)
    double complex, intent(out) :: mass(irmd)
    double complex, intent(out) :: fz(irmd,0:lmaxd)
    double complex, intent(out) :: pz(irmd,0:lmaxd)
    double precision, intent(in) :: dror(irmd)
    double precision, intent(in) :: r(irmd)
    double precision, intent(in) :: s(0:lmaxd)
    double precision, intent(in) :: vm2z(irmd)
    double precision, intent(in) :: z
    integer, intent(in) :: ipan
    integer, intent(in) :: ircut(0:ipand)
    integer, intent(in) :: irmd, ipand, lmaxd, nldau
    logical, intent(in) :: ldau
    integer, intent(in) :: lldau(lmaxd+1)
    double precision, intent(in) :: wmldauav(lmaxd+1)
    double precision, intent(in) :: ldaucut(irmd)

    double complex :: dfd0, dpd0, fip0, fip1, hamf1, k1f, k1p, k2f, k2p, k3f, k3p, k4f, k4p, mass1, pip0, pip1, vme, vmefac, vmetr1
    
    double precision :: dror1, drsm1, drsp1, s1, sm1, sp1, srafac, zocsq
    integer :: ip, ir, irc, ire, irs, irsp1, j, k, L, ildau
    double complex :: a(-1:4)
    double complex :: b(0:4)
    double complex :: dfdi(-4:0)
    double complex :: dpdi(-4:0)
    double complex :: hamfldau(irmd)
  
    CHECKASSERT(irmd > 6)

    if (nsra == 2) then
      !
      !---> in case of sra  srafac = 1/c - otherwise srafac = 0
      !
      srafac = 1.d0/cvlight

    else
      srafac = 0.d0
    endif
    !
    irc = ircut(ipan)
    !
    do ir = 2, irc
      vmetr1 = (vm2z(ir) - e)*r(ir) - 2.d0*z
      hamf(ir,0) = vmetr1*dror(ir)
      mass(ir) = r(ir) - srafac*srafac*vmetr1
    enddo ! ir
    !
    do L = 1, lmaxd
      do ir = 7, irc
        hamf(ir,L) = dble(L*L+L)/mass(ir)*dror(ir) + hamf(ir,0)
      enddo ! ir
    enddo ! L
    !
    !-----------------------------------------------------------------------
    ! LDA+U
    !
    !  Account for potential shift in case of LDA+U (averaged over m)
    !  by adding the average WLDAUAV to the spherical part of the
    !  potential.
    !
    if (ldau) then
      do ildau = 1, nldau
        if (lldau(ildau) >= 0) then
          s1 = dble(lldau(ildau)*lldau(ildau) + lldau(ildau))
          do ir = 2, irc
            vmetr1 = (vm2z(ir) - e + wmldauav(ildau)*ldaucut(ir))*r(ir)
            hamfldau(ir) = (vmetr1 - 2.d0*z)*dror(ir)
          enddo ! ir
          do ir = 7, irc
            hamf(ir,lldau(ildau)) = s1/mass(ir)*dror(ir) + hamfldau(ir)
          enddo ! ir
        endif ! lldau(ildau) >= 0
      enddo ! ildau
    endif ! ldau
    !
    ! LDA+U
    !-----------------------------------------------------------------------
    !

    do ir = 2, irc
      mass(ir) = mass(ir)*dror(ir)
    enddo ! ir
    
    do L = 0, lmaxd
      !
      s1 = s(L)
      sm1 = s1 - 1.d0
      sp1 = s1 + 1.d0
      !
      !---> loop over the number of kinks
      !
      do ip = 1, ipan
        !
        if (ip == 1) then
          irs = 2
          ire = ircut(1)

          !
          !---> initial values
          !
          vme = vm2z(2) - e
          vmefac = 1.d0 - vme*srafac*srafac
          
          if (nsra == 2 .and. z > 0.d0) then

            zocsq = -2.d0*z*z/(cvlight*cvlight)
            a(-1) = 0.d0
            a(0) = 1.d0
            b(0) = cmplx(sm1*cvlight*cvlight/(2*z), 0.d0)
            do j = 1,3
              a(j) = (0.d0,0.d0)
              b(j) = (0.d0,0.d0)
            enddo ! j

          else  ! nsra == 2 .and. z > 0.d0

            a(0) = 0.d0
            b(0) = dble(L)/vmefac
            a(1) = 1.d0
            do j = 2, 4
              a(j) = (vme*vmefac*a(j-2) - 2.d0*z*a(j-1))/dble((j-1)*(j+2*L))
              b(j-1) = dble(L+j-1)*a(j)/vmefac
            enddo ! j

          endif ! nsra == 2 .and. z > 0.d0
          !
          k = -4
          !
          !---> power series near origin
          !
          do ir = 2, 6
            pip0 = a(3)
            dpd0 = 3.d0*a(3)
            fip0 = b(3)
            dfd0 = 3.d0*b(3)
            do j = 2, 0, -1
              pip0 = a(j) + pip0*r(ir)
              dpd0 = dble(j)*a(j) + dpd0*r(ir)
              fip0 = b(j) + fip0*r(ir)
              dfd0 = dble(j)*b(j) + dfd0*r(ir)
            enddo ! j
            !
            pz(ir,L) = pip0
            fz(ir,L) = fip0
            dpdi(k) = dpd0*dror(ir)
            dfdi(k) = dfd0*dror(ir)
            !
            k = k + 1
          enddo ! ir

        else  ! ip == 1 ! panels 2, 3, ...
          !
          !---> runge kutta step to restart algorithm
          !
          irs = ircut(ip-1) + 1
          ire = ircut(ip)
          CHECKASSERT((ire - irs) >= 4) ! minimum of 5 points per panel!
          irsp1 = irs + 1
          pip0 = pz(irs,L)
          fip0 = fz(irs,L)
          drsp1 = dror(irs)*sp1
          drsm1 = dror(irs)*sm1
          dpdi(-4) = mass(irs)*  fip0 - drsm1*pip0
          dfdi(-4) = hamf(irs,L)*pip0 - drsp1*fip0
          !
          !---> first step - 4 point runge kutta with interpolation
          !
          k1p = dpdi(-4)
          k1f = dfdi(-4)
          !
          dror1 = (3.d0*dror(irs+3)   - 15.d0*dror(irs+2)   + 45.d0*dror(irsp1)   + 15.d0*dror(irs)  )/48.d0
          mass1 = (3.d0*mass(irs+3)   - 15.d0*mass(irs+2)   + 45.d0*mass(irsp1)   + 15.d0*mass(irs)  )/48.d0
          hamf1 = (3.d0*hamf(irs+3,L) - 15.d0*hamf(irs+2,L) + 45.d0*hamf(irsp1,L) + 15.d0*hamf(irs,L))/48.d0
          drsp1 = dror1*sp1
          drsm1 = dror1*sm1
          k2p = mass1*(fip0 + 0.5d0*k1f) - drsm1*(pip0 + 0.5d0*k1p)
          k2f = hamf1*(pip0 + 0.5d0*k1p) - drsp1*(fip0 + 0.5d0*k1f)
          k3p = mass1*(fip0 + 0.5d0*k2f) - drsm1*(pip0 + 0.5d0*k2p)
          k3f = hamf1*(pip0 + 0.5d0*k2p) - drsp1*(fip0 + 0.5d0*k2f)
          !
          drsp1 = dror(irsp1)*sp1
          drsm1 = dror(irsp1)*sm1
          k4p = mass(irsp1)  *(fip0 + k3f) - drsm1*(pip0 + k3p)
          k4f = hamf(irsp1,L)*(pip0 + k3p) - drsp1*(fip0 + k3f)
          pip0 = pip0 + (k1p + 2.d0*(k2p + k3p) + k4p)/6.d0
          fip0 = fip0 + (k1f + 2.d0*(k2f + k3f) + k4f)/6.d0
          !
          pz(irsp1,L) = pip0
          fz(irsp1,L) = fip0
          dpdi(-3) = mass(irsp1)*  fip0 - drsm1*pip0
          dfdi(-3) = hamf(irsp1,L)*pip0 - drsp1*fip0
          !
          k = -2
          !
          !---> 4 point runge kutta with h = i+2 - i
          !
          do ir = irs + 2, irs + 4
            pip0 = pz(ir-2,L)
            fip0 = fz(ir-2,L)
            k1p = dpdi(k-2)
            k1f = dfdi(k-2)
            k2p = mass(ir-1)* (fip0+k1f) - drsm1* (pip0+k1p)
            k2f = hamf(ir-1,L)* (pip0+k1p) - drsp1* (fip0+k1f)
            k3p = mass(ir-1)* (fip0+k2f) - drsm1* (pip0+k2p)
            k3f = hamf(ir-1,L)* (pip0+k2p) - drsp1* (fip0+k2f)
            !
            drsp1 = dror(ir)*sp1
            drsm1 = dror(ir)*sm1
            !
            k4p = mass(ir)*  (fip0 + 2.d0*k3f) - drsm1*(pip0 + 2.d0*k3p)
            k4f = hamf(ir,L)*(pip0 + 2.d0*k3p) - drsp1*(fip0 + 2.d0*k3f)
            pip0 = pip0 + (k1p + 2.d0*(k2p + k3p) + k4p)/3.d0
            fip0 = fip0 + (k1f + 2.d0*(k2f + k3f) + k4f)/3.d0
            !
            pz(ir,L) = pip0
            fz(ir,L) = fip0
            dpdi(k) = mass(ir)*  fip0 - drsm1*pip0
            dfdi(k) = hamf(ir,L)*pip0 - drsp1*fip0
            k = k + 1
          enddo ! ir
          
        endif ! ip == 1 

        ! if >5 points in panel, calculate those points with 5-point formula
        do ir = irs + 5, ire
          drsp1 = dror(ir)*sp1
          drsm1 = dror(ir)*sm1
          !
          !---> predictor : 5 point adams - bashforth
          !
          pip1 = pip0 + (1901.d0*dpdi(0) - 2774.d0*dpdi(-1) + 2616.d0*dpdi(-2) - 1274.d0*dpdi(-3) + 251.d0*dpdi(-4))/720.d0
          fip1 = fip0 + (1901.d0*dfdi(0) - 2774.d0*dfdi(-1) + 2616.d0*dfdi(-2) - 1274.d0*dfdi(-3) + 251.d0*dfdi(-4))/720.d0
          !
          dpdi(-4) = dpdi(-3)
          dpdi(-3) = dpdi(-2)
          dpdi(-2) = dpdi(-1)
          dpdi(-1) = dpdi( 0)
          dfdi(-4) = dfdi(-3)
          dfdi(-3) = dfdi(-2)
          dfdi(-2) = dfdi(-1)
          dfdi(-1) = dfdi( 0)
          !
          dpdi( 0) = mass(ir)*  fip1 - drsm1*pip1
          dfdi( 0) = hamf(ir,L)*pip1 - drsp1*fip1
          !
          !---> corrector : 5 point adams - moulton
          !
          pip0 = pip0 + (251.d0*dpdi(0) + 646.d0*dpdi(-1) - 264.d0*dpdi(-2) + 106.d0*dpdi(-3) - 19.d0*dpdi(-4))/720.d0
          fip0 = fip0 + (251.d0*dfdi(0) + 646.d0*dfdi(-1) - 264.d0*dfdi(-2) + 106.d0*dfdi(-3) - 19.d0*dfdi(-4))/720.d0
          !
          pz(ir,L) = pip0
          fz(ir,L) = fip0
          dpdi(0) = mass(ir)*  fip0 - drsm1*pip0
          dfdi(0) = hamf(ir,L)*pip0 - drsp1*fip0
        enddo ! ir
        !
        !---> remember that the r - mesh contains the kinks two times
        !     store the values of pz and fz to restart the algorithm
        !
        if (ip /= ipan) then
          pz(ire+1,L) = pip0
          fz(ire+1,L) = fip0
        endif

      enddo ! ip ! end loop over panels
      
      !
      !---> logarithmic derivate of real wavefunction ( r**s *pz / r)
      !
      dlogdp(L) = (dpdi(0)/(pip0*dror(irc)) + sm1)/r(irc)
     
    enddo ! L ! end loop over L

  endsubroutine regsol
 
 
 
  subroutine irwsol(ek,fz,hamf,mass,pz,qz,sz,dror,s,ipan,ircut,irmd,ipand,lmaxd)
!-----------------------------------------------------------------------
!  calculates the irregular solution of the schroedinger equation or
!    in semi relativistic approximation for a spherically averaged
!    potential and given energy . to achieve greater precision the
!    leading power r**-s ( in schroedinger case s = l , in case of sra
!    s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
!    from the wavefunction .
!
!
!  the differential equation is solved with a 5 point adams - bashforth
!    and adams - moulton predictor corrector method integrating
!    inwards and extended for potentials with kinks
!
!                                               b.drittler   nov.1989
!-----------------------------------------------------------------------
    double complex, intent(in) :: ek
    integer, intent(in) :: ipan,ipand,irmd,lmaxd
    double complex, intent(in) :: fz(irmd,0:lmaxd),hamf(irmd,0:lmaxd),mass(irmd), pz(irmd,0:lmaxd)
    double complex, intent(inout) :: qz(irmd,0:lmaxd), sz(irmd,0:lmaxd)
    double precision, intent(in) :: dror(irmd), s(0:lmaxd)
    integer, intent(in) :: ircut(0:ipand)

    double complex :: hamf1,k1f,k1p,k2f,k2p,k3f,k3p,k4f,k4p,mass1,qim0,qim1,sim0,sim1
    double precision :: dror1,drsm1,drsp1,s1,sm1,sp1
    integer :: ip,ir,ire,irs,irwsk,k,l
    double complex :: dqdi(0:4),dsdi(0:4)
    
    irwsk = ircut(1)/9

    do l = 0, lmaxd
      s1 = s(l)
      sm1 = s1 - 1.d0
      sp1 = s1 + 1.d0
!
!---> loop over kinks
!
      do ip = ipan, 1, -1
        irs = ircut(ip)
        ire = max(ircut(ip-1), irwsk) + 1
        drsp1 = dror(irs)*sp1
        drsm1 = dror(irs)*sm1
        qim0 = qz(irs,l)
        sim0 = sz(irs,l)
        dqdi(4) = mass(irs)*  sim0 + drsp1*qim0
        dsdi(4) = hamf(irs,l)*qim0 + drsm1*sim0
!
!---> start algorithm - 4 point runge kutta with interpolation
!
        k1p = dqdi(4)
        k1f = dsdi(4)

        dror1 = (3.d0*dror(irs-3)   - 15.d0*dror(irs-2)   + 45.d0*dror(irs-1)   + 15.d0*dror(irs)  )/48.d0
        mass1 = (3.d0*mass(irs-3)   - 15.d0*mass(irs-2)   + 45.d0*mass(irs-1)   + 15.d0*mass(irs)  )/48.d0
        hamf1 = (3.d0*hamf(irs-3,l) - 15.d0*hamf(irs-2,l) + 45.d0*hamf(irs-1,l) + 15.d0*hamf(irs,l))/48.d0
        drsp1 = dror1*sp1
        drsm1 = dror1*sm1
        k2p = mass1*(sim0 - 0.5d0*k1f) + drsp1*(qim0 - 0.5d0*k1p)
        k2f = hamf1*(qim0 - 0.5d0*k1p) + drsm1*(sim0 - 0.5d0*k1f)
        k3p = mass1*(sim0 - 0.5d0*k2f) + drsp1*(qim0 - 0.5d0*k2p)
        k3f = hamf1*(qim0 - 0.5d0*k2p) + drsm1*(sim0 - 0.5d0*k2f)
        drsp1 = dror(irs-1)*sp1
        drsm1 = dror(irs-1)*sm1
        k4p = mass(irs-1)*  (sim0 - k3f) + drsp1*(qim0 - k3p)
        k4f = hamf(irs-1,l)*(qim0 - k3p) + drsm1*(sim0 - k3f)
        qim0 = qim0 - (k1p + 2.d0*(k2p + k3p) + k4p)/6.d0
        sim0 = sim0 - (k1f + 2.d0*(k2f + k3f) + k4f)/6.d0
        qz(irs-1,l) = qim0
        sz(irs-1,l) = sim0
        dqdi(3) = mass(irs-1)*  sim0 + drsp1*qim0
        dsdi(3) = hamf(irs-1,l)*qim0 + drsm1*sim0

        k = 2
!
!---> 4 point runge kutta with h = i+2 - 1
!
        do ir = irs - 2, irs - 4, -1
          qim0 = qz(ir+2,l)
          sim0 = sz(ir+2,l)
          k1p = dqdi(k+2)
          k1f = dsdi(k+2)
          k2p = mass(ir+1)*  (sim0 - k1f) + drsp1*(qim0 - k1p)
          k2f = hamf(ir+1,l)*(qim0 - k1p) + drsm1*(sim0 - k1f)
          k3p = mass(ir+1)*  (sim0 - k2f) + drsp1*(qim0 - k2p)
          k3f = hamf(ir+1,l)*(qim0 - k2p) + drsm1*(sim0 - k2f)

          drsp1 = dror(ir)*sp1
          drsm1 = dror(ir)*sm1

          k4p = mass(ir)*  (sim0 - 2.d0*k3f) + drsp1* (qim0 - 2.d0*k3p)
          k4f = hamf(ir,l)*(qim0 - 2.d0*k3p) + drsm1* (sim0 - 2.d0*k3f)
          qim0 = qim0 - (k1p + 2.d0*(k2p + k3p) + k4p)/3.d0
          sim0 = sim0 - (k1f + 2.d0*(k2f + k3f) + k4f)/3.d0
          qz(ir,l) = qim0
          sz(ir,l) = sim0
          dqdi(k) = mass(ir)*  sim0 + drsp1*qim0
          dsdi(k) = hamf(ir,l)*qim0 + drsm1*sim0
          k = k - 1
        enddo ! ir

        do ir = irs - 5, ire, -1
!
!---> predictor : 5 point adams - bashforth
!
          qim1 = qim0 - (1901.d0*dqdi(0) - 2774.d0*dqdi(1) + 2616.d0*dqdi(2) - 1274.d0*dqdi(3) + 251.d0*dqdi(4))/720.d0
          sim1 = sim0 - (1901.d0*dsdi(0) - 2774.d0*dsdi(1) + 2616.d0*dsdi(2) - 1274.d0*dsdi(3) + 251.d0*dsdi(4))/720.d0
!
          dqdi(4) = dqdi(3)
          dqdi(3) = dqdi(2)
          dqdi(2) = dqdi(1)
          dqdi(1) = dqdi(0)
          
          dsdi(4) = dsdi(3)
          dsdi(3) = dsdi(2)
          dsdi(2) = dsdi(1)
          dsdi(1) = dsdi(0)
!
          drsp1 = dror(ir)*sp1
          drsm1 = dror(ir)*sm1
!
          dqdi(0) = mass(ir)*  sim1 + drsp1*qim1
          dsdi(0) = hamf(ir,l)*qim1 + drsm1*sim1
!
!---> corrector : 5 point adams - moulton
!
          qim0 = qim0 - (251.d0*dqdi(0) + 646.d0*dqdi(1) - 264.d0*dqdi(2) + 106.d0*dqdi(3) - 19.d0*dqdi(4))/720.d0
          sim0 = sim0 - (251.d0*dsdi(0) + 646.d0*dsdi(1) - 264.d0*dsdi(2) + 106.d0*dsdi(3) - 19.d0*dsdi(4))/720.d0

          dqdi(0) = mass(ir)*  sim0 + drsp1*qim0
          dsdi(0) = hamf(ir,l)*qim0 + drsm1*sim0

          qz(ir,l) = qim0
          sz(ir,l) = sim0
        enddo ! ir

        if (ip /= 1) then
          qz(ire-1,l) = qim0
          sz(ire-1,l) = sim0
        endif ! ip /= 1

      enddo ! ip

!
!---> use wronski relation near origin
!
      do ir = irwsk, 2, -1
!
!---> 2 point corrector - predictor
!
        qim1 = qim0 - 1.5d0*dqdi(0) + 0.5d0*dqdi(1)

        dqdi(1) = dqdi(0)
        drsp1 = dror(ir)*sp1

        dqdi(0) = mass(ir)*(1.d0/ek + qim1*fz(ir,l))/pz(ir,l) + drsp1*qim1

        qim0 = qim0 - 0.5d0*dqdi(0) - 0.5d0*dqdi(1)

        dqdi(0) = mass(ir)*(1.d0/ek + qim0*fz(ir,l))/pz(ir,l) + drsp1*qim0

        qz(ir,l) = qim0
      enddo ! ir 
!
      do ir = irwsk, 2, -1
        sz(ir,l) = (1.d0/ek + qz(ir,l)*fz(ir,l))/pz(ir,l)
      enddo ! ir
      
    enddo ! l

  endsubroutine ! irwsol
      
 
 
  subroutine vllns(vnspll, vins, cleb, icleb, iend, lmax, irmd, irnsd, ncleb)
!-----------------------------------------------------------------------
!     calculates v_ll' from v_l.
!     to determine the non - spherical wavefunctions the potential
!         has to be lm1 and lm2 dependent . the potential is stored
!         only as lm dependent , therefore a transformation in the
!         following way has to be done :
!
!        vnsll(r,lm1,lm2)   =   {  c(lm1,lm2,lm3) *vins(r,lm3)  }
!                                  (summed over lm3 at the right site )
!        where c(lm1,lm2,lm3) are the gaunt coeffients .
!
!             (see notes by b.drittler)
!
!     attention : the gaunt coeffients are stored in an index array
!                  only for lm1 > lm2
!                 (see subroutine gaunt)
!
!                               b.drittler   july 1988
!-----------------------------------------------------------------------
!                          modified by r. zeller sep. 2000
!-----------------------------------------------------------------------
    integer, intent(in) :: lmax, irmd, irnsd, ncleb, iend
    double precision, intent(in) :: cleb(ncleb,2)
    double precision, intent(in) :: vins(irmd-irnsd:irmd,(2*lmax+1)**2)
    double precision, intent(inout) :: vnspll((lmax+1)**2,(lmax+1)**2,irmd-irnsd:irmd)
    integer, intent(in) :: icleb(ncleb,3)

    integer :: ir, j, lm1, lm2, lm3, lmmaxd, irmind
    
    irmind = irmd-irnsd
    lmmaxd = (lmax+1)**2

    do ir = irmind, irmd
      do lm1 = 1, lmmaxd
        do lm2 = 1, lm1
          vnspll(lm1,lm2,ir) = 0.d0 ! clear upper triangular matrix including diagonal
        enddo ! lm2
      enddo ! lm1
    enddo ! ir

    do j = 1, iend
      lm1 = icleb(j,1)
      lm2 = icleb(j,2)
      lm3 = icleb(j,3)
      do ir = irmind, irmd
        vnspll(lm1,lm2,ir) = vnspll(lm1,lm2,ir) + cleb(j,1)*vins(ir,lm3)
      enddo ! ir
    enddo ! j

!
!---> use symmetry of the gaunt coef.
!
    do ir = irmind, irmd
      do lm1 = 1, lmmaxd
        do lm2 = 1, lm1-1
          vnspll(lm2,lm1,ir) = vnspll(lm1,lm2,ir) ! copy lower trinagular matrix
        enddo ! lm2
        vnspll(lm1,lm1,ir) = vnspll(lm1,lm1,ir) + vins(ir,1) ! add diagonal terms
      enddo ! lm1
    enddo ! ir

! todo: check if an out loop over ir for all ops is better (will depend on ratio between iend and #ir)    
  endsubroutine ! vllns
  
  
  
  subroutine wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
!-----------------------------------------------------------------------
!      determines the integrands cder, dder or ader, bder in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinspll.
!
!      r. zeller      aug. 1994
!-----------------------------------------------------------------------
    integer, intent(in) :: irmd,irmind,lmmaxd,nsra
    double complex, intent(out) :: cder(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(out) :: dder(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(in) :: pzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: qns(lmmaxd,lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: qzekdr(lmmaxd,irmind:irmd,2)
    double precision, intent(in) :: vnspll(lmmaxd,lmmaxd,irmind:irmd)

#define _USE_zgemm_in_WFINT_      
#ifndef _USE_zgemm_in_WFINT_     

    external :: dgemm ! from BLAS
    integer :: ir, lm
    double precision :: qnsi(lmmaxd,lmmaxd), qnsr(lmmaxd,lmmaxd)
    double precision :: vtqnsi(lmmaxd,lmmaxd), vtqnsr(lmmaxd,lmmaxd)
    
    do ir = irmind, irmd
    
      qnsr =  dble(qns(:,:,ir,1))
      qnsi = dimag(qns(:,:,ir,1))
      
      call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
      call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
      
      do lm = 1, lmmaxd
        cder(:,lm,ir) = qzekdr(:,ir,1)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
        dder(:,lm,ir) = pzekdr(:,ir,1)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
      enddo ! lm
  
      if (nsra == 2) then
        qnsr =  dble(qns(:,:,ir,2))
        qnsi = dimag(qns(:,:,ir,2))
  
        call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
        call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
        do lm = 1, lmmaxd
          cder(:,lm,ir) = cder(:,lm,ir) + qzekdr(:,ir,2)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
          dder(:,lm,ir) = dder(:,lm,ir) + pzekdr(:,ir,2)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
        enddo ! lm
      endif ! nsra

    enddo ! ir
#else
! new version
    external :: zgemm ! from BLAS
    integer :: ir, lm
    double complex :: cvnspll(lmmaxd,lmmaxd), vtqns(lmmaxd,lmmaxd)
    
    do ir = irmind, irmd
    
      cvnspll = dcmplx(vnspll(:,:,ir), 0.d0) ! convert the potential to complex

      call zgemm('n','n',lmmaxd,lmmaxd,lmmaxd,cone,cvnspll,lmmaxd,qns(1,1,ir,1),lmmaxd,zero,vtqns,lmmaxd)
      ! if the multiplication were the other way round, i.e. call zgemm('n','n',N,N,N,cone,qns(1,1,ir,1),N,cvnspll,N,zero,vtqns,N),
      ! we could call dgemm('n','n',2*N,N,N,cone,qns(1,1,ir,1),N,vnspll(1,1,ir),N,zero,vtqns,N)
      
      do lm = 1, lmmaxd
        cder(:,lm,ir) = qzekdr(:,ir,1)*vtqns(:,lm)
        dder(:,lm,ir) = pzekdr(:,ir,1)*vtqns(:,lm)
      enddo ! lm
  
      if (nsra == 2) then
        call zgemm('n','n',lmmaxd,lmmaxd,lmmaxd,cone,cvnspll,lmmaxd,qns(1,1,ir,2),lmmaxd,zero,vtqns,lmmaxd)
  
        do lm = 1, lmmaxd
          cder(:,lm,ir) = cder(:,lm,ir) + qzekdr(:,ir,2)*vtqns(:,lm)
          dder(:,lm,ir) = dder(:,lm,ir) + pzekdr(:,ir,2)*vtqns(:,lm)
        enddo ! lm
      endif ! nsra

    enddo ! ir
#endif
  endsubroutine wfint

      
  subroutine wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra, irmind,irmd,lmmaxd)
!-----------------------------------------------------------------------
!      determines the integrands cder, dder or ader, bder in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinspll.
!        (this subroutine is used in zeroth order born approximation,
!         otherwise subroutine wfint must be used)
!      r. zeller      aug. 1994
!-----------------------------------------------------------------------
    integer, intent(in) :: irmd,irmind,lmmaxd,nsra
    double complex, intent(out) :: cder(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(out) :: dder(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(in) :: pzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: qzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: qzlm(lmmaxd,irmind:irmd,2)
    double precision, intent(in) :: vnspll(lmmaxd,lmmaxd,irmind:irmd)
    
    integer :: ir, lm

    if (nsra == 2) then
      do ir = irmind, irmd
        do lm = 1, lmmaxd
          cder(:,lm,ir) = qzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1) + qzekdr(:,ir,2)*vnspll(:,lm,ir)*qzlm(lm,ir,2)
          dder(:,lm,ir) = pzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1) + pzekdr(:,ir,2)*vnspll(:,lm,ir)*qzlm(lm,ir,2)
        enddo ! lm
      enddo ! ir
    else
      do ir = irmind, irmd
        do lm = 1, lmmaxd
          cder(:,lm,ir) = qzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1)
          dder(:,lm,ir) = pzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1)
        enddo ! lm
      enddo ! ir
    endif
      
  endsubroutine wfint0
      
      
  subroutine wfmesh(e, ek, cvlight, nsra, z, r, s, rs, irm, irmd, lmaxd)
    double complex, intent(in) :: e
    double complex, intent(out) :: ek
    double precision, intent(in) :: cvlight, z
    integer, intent(in) :: irm, irmd, lmaxd, nsra
    double precision, intent(in) :: r(irmd)
    double precision, intent(out) :: rs(irmd,0:lmaxd), s(0:lmaxd)

    double precision :: s1
    integer :: ir, l

    if (nsra == 2) then
      ek = sqrt(e+e*e/ (cvlight*cvlight))
    else
      ! assume(nsra == 1)
      ek = sqrt(e)
    endif
      
    do l = 0, lmaxd

      if (nsra == 2) then
        s1 = sqrt(dble(l*l+l+1) - 4.d0*z*z/(cvlight*cvlight))
        if (z == 0.d0) s1 = dble(l)
      else
        s1 = dble(l)
      endif
      s(l) = s1
      rs(1,l) = 0.d0
      do ir = 2, irm
        rs(ir,l) = r(ir)**s1
      enddo ! ir
      do ir = irm+1, irmd
          rs(ir,l) = 0.d0
      enddo ! ir

    enddo ! l
    
  endsubroutine wfmesh
      
      
      
!>    @param[out] efac
!>    @param[in]  pz
!>    @param[in]  qz
!>    @param[in]  fz
!>    @param[in]  sz
!>    @param[in]  nsra
!>    @param[out] pzlm
!>    @param[out] qzlm
!>    @param[out] pzekdr
!>    @param[out] qzekdr
!>    @param[in]  ek
  subroutine wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm, pzekdr, qzekdr, ek, loflm, irmind, irmd, lmaxd, lmmaxd)
!-----------------------------------------------------------------------
!                 r. zeller      oct. 1993
!-----------------------------------------------------------------------
    double complex, intent(in) :: ek
    integer, intent(in) :: irmd, irmind, lmaxd, lmmaxd, nsra
    double complex, intent(out) :: efac(lmmaxd)
    double complex, intent(in) :: fz(irmd,0:lmaxd), pz(irmd,0:lmaxd)
    double complex, intent(out) :: pzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(out) :: pzlm(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: qz(irmd,0:lmaxd)
    double complex, intent(out) :: qzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(out) :: qzlm(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: sz(irmd,0:lmaxd)
    double precision, intent(in) :: drdi(*)
    integer, intent(in) :: loflm(*)
    
    double complex :: efac1, v1
    integer :: ir, j, l, l1, lm1, m

!
!---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!
!
    efac(1) = cone
    v1 = cone
    do l = 1, lmaxd
      v1 = v1*ek/dble(2*l-1)
      do m = -l, l
        lm1 = l*(l+1)+m+1
        efac(lm1) = v1
      enddo ! m
    enddo ! l
!
!---> get wfts of same magnitude by scaling with efac
!
    do lm1 = 1, lmmaxd
      l1 = loflm(lm1)
      efac1 = efac(lm1)
      do ir = irmind, irmd
        pzlm(lm1,ir,1) = pz(ir,l1)/efac1
        qzlm(lm1,ir,1) = qz(ir,l1)*efac1
      enddo ! ir
      if (nsra == 2) then
        do ir = irmind, irmd
          pzlm(lm1,ir,2) = fz(ir,l1)/efac1
          qzlm(lm1,ir,2) = sz(ir,l1)*efac1
        enddo ! ir
      endif ! nsra == 2

      do j = 1, nsra
        do ir = irmind, irmd
          pzekdr(lm1,ir,j) = pzlm(lm1,ir,j)*ek*drdi(ir)
          qzekdr(lm1,ir,j) = qzlm(lm1,ir,j)*ek*drdi(ir)
        enddo ! ir
      enddo ! j
  
    enddo ! lm1

  endsubroutine wftsca

endmodule SingleSiteHelpers_mod