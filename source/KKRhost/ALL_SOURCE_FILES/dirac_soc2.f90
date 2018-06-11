subroutine dirabmsoc2(getirrsol, c, socscl, it, e, l, mj, kap1, kap2, pis, &
  cg1, cg2, cg4, cg5, cg8, v, b, z, nucleus, r, drdi, dovr, nmesh, dxp, pr, &
  qr, pi, qi, d_p, dq, nrmax)
!   ********************************************************************
!   *                                                                  *
!   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
!   *                                                                  *
!   *   with the SOC-operator manipulated                              *
!   *                                                                  *
!   *           H_soc = xi [ (sig_x*l_x+sig*l_y) + sig_z*l_z ]         *
!   *                                                                  *
!   *   SOCSCL = -1  ==  IXY = 0 : neglect 1st term of SOC-operator    *
!   *   SOCSCL = -2: ==  IXY = 1 : neglect 2nd term of SOC-operator    *
!   *                                                                  *
!   *   the outward integration is started by a power expansion        *
!   *   and continued by ADAMS-BASHFORTH-MOULTON - pred./corr.-method  *
!   *   NABM = 4(5) selects the 4(5)-point formula                     *
!   *                                                                  *
!   *   the inward integration is started analytically                 *
!   *                                                                  *
!   *   returns the wave functions up to the mesh point NMESH          *
!   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
!   *   and    R/I standing for regular/irregular solution             *
!   *                                                                  *
!   *  12/12/97  HE                                                    *
!   *  21/01/98  HE  finite nucelus                                    *
!   ********************************************************************

  use :: mod_types, only: t_inc
  use mod_DataTypes, only: dp
  implicit none

! PARAMETER definitions

  integer :: mpsmax, npemax, nabm
  parameter (mpsmax=40, npemax=4)
  parameter (nabm=4)
!PARAMETER ( NABM   =      5 )
  complex (kind=dp) :: cz
  parameter (cz=(0.0d0,0.0d0))
  real (kind=dp) :: tol
  parameter (tol=1.0d-9)
  integer :: itmax
  parameter (itmax=50)

! Dummy arguments
  real (kind=dp) :: c, cg1, cg2, cg4, cg5, cg8, mj, socscl
  complex (kind=dp) :: e
  logical :: getirrsol
  integer :: it, kap1, kap2, l, nmesh, nrmax, nucleus, z
  complex (kind=dp) :: pis
  real (kind=dp) :: b(nrmax), dovr(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
  complex (kind=dp) :: d_p(2, 2, nrmax), dq(2, 2, nrmax), dxp(2, 2), &
    pi(2, 2, nrmax), pr(2, 2, nrmax), qi(2, 2, nrmax), qr(2, 2, nrmax)

! Local variables
  complex (kind=dp) :: aa11, aa12, aa21, aa22, arg, bb1, bb2, bpp, bqq, cfac, &
    d14, db14, dbh, dh, diffa, diffb, dv14, dvh, emvpp, emvqq, s0, t0
  real (kind=dp) :: acorr(0:nabm-1), acorr0(0:nabm-1), apred(nabm), &
    apred0(nabm), astep, b14, bc(0:npemax), bh, bhlp(nabm+4), cgd(2), cgmd(2), &
    cgo, cgoz, cgz(2), cm(npemax, npemax), cmi(npemax, npemax), csqr, &
    dbdr(nrmax), dbhlp(nabm+4), dhlp(nabm+4), dvdr(nrmax), dvhlp(nabm+4), &
    gam(2), gpm, hlp(nabm+4), hlp1, kap(2), kpx(2), kpy(2), r14, rh, &
    rhlp(nabm+4), rpwgpm, rr, sk(2), sk1, sk2, so2, so6, srk, tz, v14, &
    vc(0:npemax), vh, vhlp(nabm+4), x14, xh
  complex (kind=dp) :: cjlz
  real (kind=dp) :: dabs, dble, dsqrt
  complex (kind=dp) :: detd, mp1(2, 2), mp2(2, 2), mp3(2, 2), mp4(2, 2), &
    mq1(2, 2), mq2(2, 2), mq3(2, 2), mq4(2, 2), p1(2, 2), p2(2, 2), p3(2, 2), &
    p4(2, 2), pc(2, 2, -npemax:mpsmax), pnew(2, 2), pold(2, 2), q1(2, 2), &
    q2(2, 2), q3(2, 2), q4(2, 2), qc(2, 2, -npemax:mpsmax), qnew(2, 2), &
    qold(2, 2), socpp(2), zz
  integer :: i, ic, ip, irk, isk1, isk2, iv, ixy, j, jcorr, k, lb(2), lb1, &
    lb2, m, mps, n, nacorr, ndiv, nhlp, nm, npe, nsol, ntop
  integer :: int, isign, nint
  real (kind=dp) :: ylag

  data apred0/55.0d0, -59.0d0, +37.0d0, -9.0d0/
  data acorr0/9.0d0, +19.0d0, -5.0d0, +1.0d0/
  data astep/24.0d0/

  csqr = c*c
  cfac = pis*c/(e+csqr)

! find   NPE  expansion coefficients for the potential and b-field
  npe = 4

  tz = dble(2*z)

  do iv = 1, npe
    do n = 1, npe
      cm(n, iv) = r(n)**(iv-1)
    end do
  end do

  call rinvgj(cmi, cm, npemax, npe)

  do iv = 1, npe
    vc(iv-1) = 0.0d0
    do n = 1, npe
      vc(iv-1) = vc(iv-1) + cmi(iv, n)*(v(n)+tz/r(n))
    end do
  end do

  do iv = 1, npe
    bc(iv-1) = 0.0d0
    do n = 1, npe
      bc(iv-1) = bc(iv-1) + cmi(iv, n)*b(n)
    end do
  end do

!*HF* IXY=0 ==> ONLY ZZ-PART, IXY=1 ==> ONLY XY-PART
  if (nint(socscl)==-1) then
    ixy = 0
  else
    ixy = 1
  end if
!*HF*

!    calculate g-coefficients of b-field

  isk1 = isign(1, kap1)
  isk2 = isign(1, kap2)
  sk1 = dble(isk1)
  sk2 = dble(isk2)
  lb1 = l - isk1
  lb2 = l - isk2

  cg1 = -mj/(kap1+0.5d0)
  cg5 = -mj/(-kap1+0.5d0)
  cgd(1) = cg1
  cgmd(1) = cg5
  kap(1) = dble(kap1)
  gam(1) = dsqrt(kap(1)**2-(tz/c)**2)
  if (ixy==0) then
    cgz(1) = kap(1) + 1 + mj*cgd(1) - 0.5d0
  else
    cgz(1) = -(mj*cgd(1)-0.5d0)
  end if
  lb(1) = lb1
  sk(1) = sk1
  if (dabs(mj)>l) then
    cg2 = 0.0d0
    cg4 = 0.0d0
    cg8 = 0.0d0
    nsol = 1
    cgd(2) = 0.0d0
    cgo = 0.0d0
    cgoz = 0.0d0
    cgmd(2) = 0.0d0
    gam(2) = 0.0d0
    kap(2) = 0.0d0
    cgz(2) = 0.0d0
    lb(2) = 0
    sk(2) = 0.0d0
  else
    cg2 = -dsqrt(1.0d0-(mj/(kap1+0.5d0))**2)
    cg4 = -mj/(kap2+0.5d0)
    cg8 = -mj/(-kap2+0.5d0)
    nsol = 2
    cgd(2) = cg4
    cgo = cg2
    cgmd(2) = cg8
    kap(2) = dble(kap2)
    gam(2) = dsqrt(kap(2)**2-(tz/c)**2)
    if (ixy==0) then
      cgz(2) = kap(2) + 1 + mj*cgd(2) - 0.5d0
      cgoz = mj*cgo
    else
      cgz(2) = -(mj*cgd(2)-0.5d0)
      cgoz = -mj*cgo
    end if
    lb(2) = lb2
    sk(2) = sk2
  end if
  do i = 1, nsol
    kpx(i) = kap(i)
    kpy(i) = kap(i)
  end do

  call dvdrspline(v, r, dvdr, nmesh)
  call derspl(nmesh, r, b, dbdr)

  do ip = 1, nabm
    ic = ip - 1
    apred(ip) = apred0(ip)/astep
    acorr(ic) = acorr0(ic)/astep
  end do
  nacorr = nabm - 1

  do i = 1, 2
    do j = 1, 2
      do ip = -npemax, mpsmax
        pc(i, j, ip) = cz
        qc(i, j, ip) = cz
      end do
    end do
  end do

! ======================================================================

  if ((tz>=2) .and. (nucleus==0)) then

    do j = 1, nsol
      i = 3 - j
      pc(j, j, 0) = dsqrt(abs(kpy(j))-gam(j))
      qc(j, j, 0) = (kpy(j)+gam(j))*(csqr/tz)*pc(j, j, 0)
      pc(i, j, 0) = cz
      qc(i, j, 0) = cz
    end do

!  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS

    mps = 40

    aa12 = -tz/csqr
    aa21 = tz
    emvqq = (e-vc(0)+csqr)/csqr
    emvpp = -e + vc(0)
    bqq = bc(0)/csqr
    do i = 1, nsol
      socpp(i) = cz
    end do

    do j = 1, nsol

      do m = 1, mps
        do i = 1, nsol
          k = 3 - i
          bb1 = (emvqq+bqq*cgmd(i))*qc(i, j, m-1)
          bb2 = (emvpp+bc(0)*cgd(i))*pc(i, j, m-1) + bc(0)*cgo*pc(k, j, m-1)
          do ip = 1, npe - 1
            bb1 = bb1 + (-vc(ip)+bc(ip)*cgmd(i))*qc(i, j, m-1-ip)/csqr
            bb2 = bb2 + (+vc(ip)+bc(ip)*cgd(i))*pc(i, j, m-1-ip) + &
              bc(ip)*cgo*pc(k, j, m-1-ip)
          end do

          aa11 = gam(j) + m + kpy(i)
          aa22 = gam(j) + m - kpy(i)
          detd = aa11*aa22 - aa12*aa21
          pc(i, j, m) = (bb1*aa22-aa12*bb2)/detd
          qc(i, j, m) = (aa11*bb2-bb1*aa21)/detd
        end do
      end do

    end do


!  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
!  FOR THE FIRST   NABM   R - MESH - POINTS

    do n = 1, nabm
      rr = r(n)

      do j = 1, nsol
        rpwgpm = rr**gam(j)

        do i = 1, nsol
          pr(i, j, n) = pc(i, j, 0)*rpwgpm
          qr(i, j, n) = qc(i, j, 0)*rpwgpm
          d_p(i, j, n) = pc(i, j, 0)*rpwgpm*gam(j)*dovr(n)
          dq(i, j, n) = qc(i, j, 0)*rpwgpm*gam(j)*dovr(n)
        end do

        do m = 1, mps
          rpwgpm = rpwgpm*rr
          gpm = gam(j) + m

          do i = 1, nsol
            pr(i, j, n) = pr(i, j, n) + pc(i, j, m)*rpwgpm
            qr(i, j, n) = qr(i, j, n) + qc(i, j, m)*rpwgpm
            d_p(i, j, n) = d_p(i, j, n) + pc(i, j, m)*rpwgpm*gpm*dovr(n)
            dq(i, j, n) = dq(i, j, n) + qc(i, j, m)*rpwgpm*gpm*dovr(n)
          end do

        end do
      end do
    end do
! ======================================================================
!                                  == EMPTY SPHERE  or FINITE NUCLEUS ==
  else

!        assume constant pot: V=V(1)   ignore coupling: B=0

    t0 = e - v(1)
    s0 = (e-v(1))/csqr + 1

    do n = 1, nabm

      do j = 1, nsol
        do i = 1, nsol
          pr(i, j, n) = cz
          qr(i, j, n) = cz
          d_p(i, j, n) = cz
          dq(i, j, n) = cz
        end do
      end do

      zz = sqrt(s0*t0)*r(n)

      do j = 1, nsol
        pr(j, j, n) = cjlz(l, zz)*r(n)
        d_p(j, j, n) = (dble(l+1)*cjlz(l,zz)-zz*cjlz(l+1,zz))*drdi(n)

        qr(j, j, n) = (d_p(j,j,n)/drdi(n)+pr(j,j,n)*(kap(j)/r(n)))/s0
        dq(j, j, n) = qr(j, j, n)*(kap(j)/r(n)) - pr(j, j, n)*t0
      end do
    end do

  end if
! ===================================================================

! =============================================================== N ====
!     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-BASHFORTH-MOULTON)

  do n = nabm + 1, nmesh

!    EVALUATE PREDICTOR

    do j = 1, nsol
      do i = 1, nsol
        pnew(i, j) = pr(i, j, n-1)
        qnew(i, j) = qr(i, j, n-1)

        do ip = 1, nabm
          pnew(i, j) = pnew(i, j) + apred(ip)*d_p(i, j, n-ip)
          qnew(i, j) = qnew(i, j) + apred(ip)*dq(i, j, n-ip)
        end do
      end do
    end do

    emvqq = (e-v(n)+csqr)*drdi(n)/csqr
    emvpp = -(e-v(n))*drdi(n)
    bqq = b(n)*drdi(n)/csqr
    bpp = b(n)*drdi(n)
    do i = 1, nsol
      socpp(i) = dovr(n)*csqr*(dvdr(n)-dbdr(n)*cgmd(i))/ &
        (e-v(n)+csqr+b(n)*cgmd(i))**2

    end do

!    EVALUATE CORRECTOR


    do jcorr = 1, itmax

      do j = 1, nsol
        do i = 1, nsol
          k = 3 - i
          pold(i, j) = pnew(i, j)
          qold(i, j) = qnew(i, j)
          d_p(i, j, n) = -kpx(i)*pnew(i, j)*dovr(n) + &
            (emvqq+bqq*cgmd(i))*qnew(i, j)
          dq(i, j, n) = kpx(i)*qnew(i, j)*dovr(n) + (emvpp+bpp*cgd(i))*pnew(i, &
            j) + (bpp*cgo+socpp(k)*cgoz)*pnew(k, j) + &
            socpp(i)*cgz(i)*pnew(i, j)
          pnew(i, j) = pr(i, j, n-1)
          qnew(i, j) = qr(i, j, n-1)
          do ic = 0, nacorr
            pnew(i, j) = pnew(i, j) + acorr(ic)*d_p(i, j, n-ic)
            qnew(i, j) = qnew(i, j) + acorr(ic)*dq(i, j, n-ic)
          end do
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          diffa = pold(i, j) - pnew(i, j)
          if (abs(real(diffa, kind=dp))>(tol*abs(real(pnew(i,j), kind=dp)))) go to 100
          if (abs(aimag(diffa))>(tol*abs(aimag(pnew(i,j))))) go to 100

          diffb = qold(i, j) - qnew(i, j)
          if (abs(real(diffb, kind=dp))>(tol*abs(real(qnew(i,j), kind=dp)))) go to 100
          if (abs(aimag(diffb))>(tol*abs(aimag(qnew(i,j))))) go to 100
        end do
      end do
      go to 110

100 end do
    if (t_inc%i_write>0) write (1337, 140) kap1, n, r(n), diffa, diffb, it, l, &
      int(2*mj), 'REG'

!                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS

110 continue
    do j = 1, nsol
      do i = 1, nsol
        k = 3 - i
        pr(i, j, n) = pnew(i, j)
        qr(i, j, n) = qnew(i, j)
        d_p(i, j, n) = -kpx(i)*pnew(i, j)*dovr(n) + (emvqq+bqq*cgmd(i))*qnew(i, &
          j)
        dq(i, j, n) = kpx(i)*qnew(i, j)*dovr(n) + (emvpp+bpp*cgd(i))*pnew(i, j &
          ) + (bpp*cgo+socpp(k)*cgoz)*pnew(k, j) + socpp(i)*cgz(i)*pnew(i, j)
      end do
    end do


  end do
! =============================================================== N ====

  do i = 1, 2
    do j = 1, 2
      dxp(i, j) = d_p(i, j, nmesh)
    end do
  end do


  call cinit(2*2*nrmax, pi)
  call cinit(2*2*nrmax, qi)

  if (.not. getirrsol) return

! #####################################################################

!             IRREGULAR SOLUTION IRREGULAR SOLUTION  IRREGULAR SOLUTION

!  CALCULATE THE INITIAL VALUES OF THE WAVEFUNCTION AT THE SPHERE
!  BOUNDARY

  do n = nmesh, nmesh + nabm
    arg = pis*r(n)

    do j = 1, nsol
      i = 3 - j
      pi(j, j, n) = cjlz(l, arg)*r(n)
      qi(j, j, n) = cfac*sk(j)*cjlz(lb(j), arg)*r(n)*c
      d_p(j, j, n) = (dble(l+1)*cjlz(l,arg)-arg*cjlz(l+1,arg))*drdi(n)
      m = lb(j)
      dq(j, j, n) = cfac*sk(j)*(dble(m+1)*cjlz(m,arg)-arg*cjlz(m+1,arg))* &
        drdi(n)*c

      pi(i, j, n) = cz
      qi(i, j, n) = cz
      d_p(i, j, n) = cz
      dq(i, j, n) = cz
    end do
  end do
!        ------------------------------------------------------------
!              INITIALIZE INWARD INTEGRATION WITH RUNGE - KUTTA
!        ------------------------------------------------------------
  ndiv = 60
  if (ndiv/=0) then

    srk = 1.0d0/dble(ndiv)
    so2 = srk/2.0d0
    so6 = srk/6.0d0

    n = nmesh

    emvqq = (e-v(n)+csqr)*drdi(n)/csqr
    emvpp = -(e-v(n))*drdi(n)
    bqq = b(n)*drdi(n)/csqr
    bpp = b(n)*drdi(n)
    do i = 1, nsol
      socpp(i) = dovr(n)*csqr*(dvdr(n)-dbdr(n)*cgmd(i))/ &
        (e-v(n)+csqr+b(n)*cgmd(i))**2

    end do

! *** reinitialize Q using only DP and PI

    do j = 1, nsol
      i = 3 - j
      qi(j, j, n) = (d_p(j,j,n)+kpx(j)*pi(j,j,n)*dovr(n))/(emvqq+bqq*cgmd(j))
      qi(i, j, n) = cz
    end do

    do j = 1, nsol
      do i = 1, nsol
        k = 3 - i
        dq(i, j, n) = kpx(i)*qi(i, j, n)*dovr(n) + (emvpp+bpp*cgd(i))*pi(i, j, &
          n) + (bpp*cgo+socpp(k)*cgoz)*pi(k, j, n) + &
          socpp(i)*cgz(i)*pi(i, j, n)
      end do
    end do

    do j = 1, nsol
      do i = 1, nsol
        p1(i, j) = pi(i, j, n)
        q1(i, j) = qi(i, j, n)
        mp1(i, j) = d_p(i, j, n)
        mq1(i, j) = dq(i, j, n)
      end do
    end do

    x14 = dble(n)
    nhlp = nabm + 4
    hlp1 = dble(nmesh-nhlp)
    do i = 1, nhlp
      hlp(i) = dble(i)
      vhlp(i) = v(nmesh-nhlp+i)
      bhlp(i) = b(nmesh-nhlp+i)
      dhlp(i) = drdi(nmesh-nhlp+i)
      rhlp(i) = r(nmesh-nhlp+i)
      dvhlp(i) = dvdr(nmesh-nhlp+i)
      dbhlp(i) = dbdr(nmesh-nhlp+i)
    end do

    do irk = 1, (nabm-1)*ndiv

      xh = x14 - so2
      vh = ylag(xh-hlp1, hlp, vhlp, 0, 3, nhlp)
      bh = ylag(xh-hlp1, hlp, bhlp, 0, 3, nhlp)
      rh = ylag(xh-hlp1, hlp, rhlp, 0, 3, nhlp)
      dh = ylag(xh-hlp1, hlp, dhlp, 0, 3, nhlp)
      dvh = ylag(xh-hlp1, hlp, dvhlp, 0, 3, nhlp)
      dbh = ylag(xh-hlp1, hlp, dbhlp, 0, 3, nhlp)

      emvqq = (e-vh+csqr)*dh/csqr
      emvpp = -(e-vh)*dh
      bqq = bh*dh/csqr
      bpp = bh*dh
      do i = 1, nsol
        socpp(i) = (dh/rh)*csqr*(dvh-dbh*cgmd(i))/(e-vh+csqr+bh*cgmd(i))**2

      end do
      n = nmesh - irk/ndiv - n + nint(xh)

      do j = 1, nsol
        do i = 1, nsol
          p2(i, j) = p1(i, j) - so2*mp1(i, j)
          q2(i, j) = q1(i, j) - so2*mq1(i, j)
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          k = 3 - i
          mp2(i, j) = -kpx(i)*p2(i, j)*dh/rh + (emvqq+bqq*cgmd(i))*q2(i, j)
          mq2(i, j) = kpx(i)*q2(i, j)*dh/rh + (emvpp+bpp*cgd(i))*p2(i, j) + &
            (bpp*cgo+socpp(k)*cgoz)*p2(k, j) + socpp(i)*cgz(i)*p2(i, j)
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          p3(i, j) = p1(i, j) - so2*mp2(i, j)
          q3(i, j) = q1(i, j) - so2*mq2(i, j)
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          k = 3 - i
          mp3(i, j) = -kpx(i)*p3(i, j)*dh/rh + (emvqq+bqq*cgmd(i))*q3(i, j)
          mq3(i, j) = kpx(i)*q3(i, j)*dh/rh + (emvpp+bpp*cgd(i))*p3(i, j) + &
            (bpp*cgo+socpp(k)*cgoz)*p3(k, j) + socpp(i)*cgz(i)*p3(i, j)
        end do
      end do

      x14 = x14 - srk
      v14 = ylag(x14-hlp1, hlp, vhlp, 0, 3, nhlp)
      b14 = ylag(x14-hlp1, hlp, bhlp, 0, 3, nhlp)
      r14 = ylag(x14-hlp1, hlp, rhlp, 0, 3, nhlp)
      d14 = ylag(x14-hlp1, hlp, dhlp, 0, 3, nhlp)
      dv14 = ylag(x14-hlp1, hlp, dvhlp, 0, 3, nhlp)
      db14 = ylag(x14-hlp1, hlp, dbhlp, 0, 3, nhlp)


      emvqq = (e-v14+csqr)*d14/csqr
      emvpp = -(e-v14)*d14
      bqq = b14*d14/csqr
      bpp = b14*d14
      do i = 1, nsol
        socpp(i) = (d14/r14)*csqr*(dv14-db14*cgmd(i))/ &
          (e-v14+csqr+b14*cgmd(i))**2

      end do

      do j = 1, nsol
        do i = 1, nsol
          p4(i, j) = p1(i, j) - srk*mp3(i, j)
          q4(i, j) = q1(i, j) - srk*mq3(i, j)
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          k = 3 - i
          mp4(i, j) = -kpx(i)*p4(i, j)*d14/r14 + (emvqq+bqq*cgmd(i))*q4(i, j)
          mq4(i, j) = kpx(i)*q4(i, j)*d14/r14 + (emvpp+bpp*cgd(i))*p4(i, j) + &
            (bpp*cgo+socpp(k)*cgoz)*p4(k, j) + socpp(i)*cgz(i)*p4(i, j)
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          p1(i, j) = p1(i, j) - so6*(mp1(i,j)+2*(mp2(i,j)+mp3(i,j))+mp4(i,j))
          q1(i, j) = q1(i, j) - so6*(mq1(i,j)+2*(mq2(i,j)+mq3(i,j))+mq4(i,j))
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          k = 3 - i
          mp1(i, j) = -kpx(i)*p1(i, j)*d14/r14 + (emvqq+bqq*cgmd(i))*q1(i, j)
          mq1(i, j) = kpx(i)*q1(i, j)*d14/r14 + (emvpp+bpp*cgd(i))*p1(i, j) + &
            (bpp*cgo+socpp(k)*cgoz)*p1(k, j) + socpp(i)*cgz(i)*p1(i, j)
        end do
      end do

      if (mod(irk,ndiv)==0) then
        n = nmesh - irk/ndiv
        if (abs(x14-dble(n))>1.0d-5) then
          write (*, *) ' <DIRAC> RUNGE-KUTTA: ', irk, ndiv, n, x14
          stop
        end if
        do j = 1, nsol
          do i = 1, nsol
            pi(i, j, n) = p1(i, j)
            qi(i, j, n) = q1(i, j)
            d_p(i, j, n) = mp1(i, j)
            dq(i, j, n) = mq1(i, j)
          end do
        end do

      end if

    end do

  end if

! =============================================================== N ====

!     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-BASHFORTH-MOULTON)

  if (ndiv/=0) then
    ntop = nmesh - nabm
  else
    ntop = nmesh
  end if

  do nm = 1, ntop
    n = 1 + ntop - nm

!    EVALUATE PREDICTOR

    do j = 1, nsol
      do i = 1, nsol
        pnew(i, j) = pi(i, j, n+1)
        qnew(i, j) = qi(i, j, n+1)

        do ip = 1, nabm
          pnew(i, j) = pnew(i, j) - apred(ip)*d_p(i, j, n+ip)
          qnew(i, j) = qnew(i, j) - apred(ip)*dq(i, j, n+ip)
        end do
      end do
    end do

    emvqq = (e-v(n)+csqr)*drdi(n)/csqr
    emvpp = -(e-v(n))*drdi(n)
    bqq = b(n)*drdi(n)/csqr
    bpp = b(n)*drdi(n)
    do i = 1, nsol
      socpp(i) = dovr(n)*csqr*(dvdr(n)-dbdr(n)*cgmd(i))/ &
        (e-v(n)+csqr+b(n)*cgmd(i))**2
    end do

!    EVALUATE CORRECTOR

    do jcorr = 1, itmax
      do j = 1, nsol
        do i = 1, nsol
          k = 3 - i
          pold(i, j) = pnew(i, j)
          qold(i, j) = qnew(i, j)
          d_p(i, j, n) = -kpx(i)*pnew(i, j)*dovr(n) + &
            (emvqq+bqq*cgmd(i))*qnew(i, j)
          dq(i, j, n) = kpx(i)*qnew(i, j)*dovr(n) + (emvpp+bpp*cgd(i))*pnew(i, &
            j) + (bpp*cgo+socpp(k)*cgoz)*pnew(k, j) + &
            socpp(i)*cgz(i)*pnew(i, j)

          pnew(i, j) = pi(i, j, n+1)
          qnew(i, j) = qi(i, j, n+1)
          do ic = 0, nacorr
            pnew(i, j) = pnew(i, j) - acorr(ic)*d_p(i, j, n+ic)
            qnew(i, j) = qnew(i, j) - acorr(ic)*dq(i, j, n+ic)
          end do
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          diffa = pold(i, j) - pnew(i, j)
          if (abs(real(diffa, kind=dp))>(tol*abs(real(pnew(i,j), kind=dp)))) go to 120
          if (abs(aimag(diffa))>(tol*abs(aimag(pnew(i,j))))) go to 120

          diffb = qold(i, j) - qnew(i, j)
          if (abs(real(diffb, kind=dp))>(tol*abs(real(qnew(i,j), kind=dp)))) go to 120
          if (abs(aimag(diffb))>(tol*abs(aimag(qnew(i,j))))) go to 120
        end do
      end do
      go to 130

120 end do
    if (t_inc%i_write>0) write (1337, 140) kap1, n, r(n), diffa, diffb, it, l, &
      int(2*mj), 'IRR'

!                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS

130 continue
    do j = 1, nsol
      do i = 1, nsol
        k = 3 - i
        pi(i, j, n) = pnew(i, j)
        qi(i, j, n) = qnew(i, j)
        d_p(i, j, n) = -kpx(i)*pnew(i, j)*dovr(n) + (emvqq+bqq*cgmd(i))*qnew(i, &
          j)
        dq(i, j, n) = kpx(i)*qnew(i, j)*dovr(n) + (emvpp+bpp*cgd(i))*pnew(i, j &
          ) + (bpp*cgo+socpp(k)*cgoz)*pnew(k, j) + socpp(i)*cgz(i)*pnew(i, j)
      end do
    end do

  end do
! =============================================================== N ====

!     the minor component for the soc-manipulated wf is meaningless
!     =>  set it to zero

  call cinit(2*2*nrmax, qr)
  call cinit(2*2*nrmax, qi)
140 format (' PRE/CORR NOT CONV. IN <DIRABMSOC2> ', 2i4, f10.7, 2x, 4e12.4, &
    3i2, '/2 ', a3)

end subroutine

subroutine dvdrspline(v, r, dvdr, n)
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE DERIVATIVE OF THE POTENTIAL                          *
!   *   DV/DR= ( D(R**2*V)/DR-2.R*V ) / R**2                           *
!   *                                                                  *
!   *   calling J.REDINGER's routine <DERSPL>                          *
!   *                                                                  *
!   * 27/10/94  HE  adopted for SPRKKR-package                         *
!   ********************************************************************
  use mod_DataTypes, only: dp
  implicit none

! Dummy arguments
  integer :: n
  real (kind=dp) :: dvdr(n), r(n), v(n)

! Local variables
  real (kind=dp) :: dr2vdr(n), r2v(n)
  integer :: i

  do i = 1, n
    r2v(i) = v(i)*r(i)**2
  end do

  call derspl(n, r, r2v, dr2vdr)

  do i = 1, n
    dvdr(i) = (dr2vdr(i)-2.d0*r2v(i)/r(i))/r(i)**2
  end do
end subroutine

subroutine derspl(n, x, f, d)
!   ********************************************************************
!   *                                                                  *
!   *   F(I) ARE THE FUNCTION VALUES AT THE POINTS X(I) FOR I=1,N      *
!   *   AND THE SPLINE DERIVATIVES D(I) ARE FOUND.                     *
!   *   THE DIMENSION OF A MUST NOT BE LESS THAN 3*N.                  *
!   *                                                                  *
!   *   sepp redinger 1985                                             *
!   *                                                                  *
!   ********************************************************************
  use :: mod_types, only: t_inc
  use mod_DataTypes, only: dp
  implicit none

! Dummy arguments
  integer :: n
  real (kind=dp) :: d(n), f(n), x(n)

! Local variables
  real (kind=dp) :: a(n*3), h1, h2, p
  integer :: i, j, k

  do i = 2, n
    if (x(i)<=x(i-1)) then
      if (t_inc%i_write>0) write (1337, 110) i
      a(1) = 1.d0
      return
    end if
  end do
  do i = 1, n
    j = 2
    if (i/=1) then
      j = n - 1
      if (i/=n) then
        h1 = 1.d0/(x(i)-x(i-1))
        h2 = 1.d0/(x(i+1)-x(i))
        a(3*i-2) = h1
        a(3*i-1) = 2.d0*(h1+h2)
        a(3*i) = h2
        d(i) = 3*(f(i+1)*h2*h2+f(i)*(h1*h1-h2*h2)-f(i-1)*h1*h1)
        go to 100
      end if
    end if
    h1 = 1.d0/(x(j)-x(j-1))
    h2 = 1.d0/(x(j+1)-x(j))
    a(3*i-2) = h1*h1
    a(3*i-1) = h1*h1 - h2*h2
    a(3*i) = -h2*h2
    d(i) = 2.d0*(f(j)*(h2*h2*h2+h1*h1*h1)-f(j+1)*h2*h2*h2-f(j-1)*h1*h1*h1)
100 end do
  p = a(4)/a(1)
  a(5) = a(5) - p*a(2)
  a(6) = a(6) - p*a(3)
  d(2) = d(2) - p*d(1)
  do i = 3, n
    k = 3*i - 4
    p = a(k+2)/a(k)
    a(k+3) = a(k+3) - p*a(k+1)
    d(i) = d(i) - p*d(i-1)
    if (i==n-1) then
      p = a(k+5)/a(k)
      a(k+5) = a(k+6) - p*a(k+1)
      a(k+6) = a(k+7)
      d(n) = d(n) - p*d(n-2)
    end if
  end do
  d(n) = d(n)/a(3*n-1)
  do i = 3, n
    j = n + 2 - i
    d(j) = (d(j)-a(3*j)*d(j+1))/a(3*j-1)
  end do
  d(1) = (d(1)-d(2)*a(2)-d(3)*a(3))/a(1)
  a(1) = 0.d0
110 format (' RETURN FROM DERSPL  ', i3, ' OUT OF ORDER')
end subroutine
