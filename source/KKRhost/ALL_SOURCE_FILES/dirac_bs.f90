subroutine dirbs(getirrsol, c, e, l, mj, kap1, kap2, pis, cg1, cg2, cg4, cg5, &
  cg8, v, b, z, nucleus, r, drdi, dovr, nmesh, pr, qr, pi, qi, dp, dq)
!   ********************************************************************
!   *                                                                  *
!   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
!   *                                                                  *
!   *   the outward integration is started by a power expansion        *
!   *   and the inward integration is started analytically             *
!   *   the integration itself is done by the BURLISCH-STOER method    *
!   *   see: numerical recipes chapter 15.4                            *
!   *                                                                  *
!   *   returns the wave functions up to the mesh point NMESH          *
!   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
!   *   and    R/I standing for regular/irregular solution             *
!   *                                                                  *
!   *   bug fixed 93/11/24                                             *
!   *  31/10/94  HE  arg. list changed - return P,Q instead of g,f     *
!   *  06/12/94  HE  CM real                                           *
!   *  29/04/95  MB  Adopted for finite nucleus                        *
!   ********************************************************************
  implicit none
  include 'sprkkr_rmesh.dim'

! PARAMETER definitions
  integer :: mpsmax, npemax, nabm
  parameter (mpsmax=40, npemax=4, nabm=5)
  double complex :: c0
  parameter (c0=(0.0d0,0.0d0))
  double precision :: epsbs
  parameter (epsbs=2.0d-7)

! COMMON variables
  double precision :: cgd(2), cgmd(2), cgo, csqr, kap(2)
  double complex :: ebs
  integer :: nradbs, nsolbs
  common /commbs/ebs, csqr, cgd, cgmd, cgo, kap, nsolbs, nradbs

! Dummy arguments
  double precision :: c, cg1, cg2, cg4, cg5, cg8, mj
  double complex :: e
  logical :: getirrsol
  integer :: kap1, kap2, l, nmesh, nucleus, z
  double complex :: pis
  double precision :: b(nrmax), dovr(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
  double complex :: dp(2, 2, nrmax), dq(2, 2, nrmax), pi(2, 2, nrmax), &
    pr(2, 2, nrmax)
  double complex :: qi(2, 2, nrmax), qr(2, 2, nrmax)

! Local variables
  double complex :: a11, a12, a21, a22, aa11, aa12, aa21, aa22, alpha, bb1, &
    bb2, beta
  double complex :: bqq, cfac, emvpp, emvqq, w1, w2, w3, w4, w5, w6, w7
  double precision :: bc(0:npemax), cm(npemax, npemax), cmi(npemax, npemax), &
    dix
  double precision :: gam(2), gpm, hbs, rpwgpm, rr, sk(2), sk1, sk2, tz, &
    vc(0:npemax), x
  double complex :: cjlz

  double complex :: detd, dy(ncfmax), efac, fy(ncfmax), &
    pc(2, 2, -npemax:mpsmax)
  double complex :: qc(2, 2, -npemax:mpsmax), zz

  integer :: i, ip, isk1, isk2, iv, j, k, lb(2), lb1, lb2, m, mps, n, nfy, &
    npe, nsol
  integer :: isign
  save :: a11, a12, a21, a22, aa11, aa12, aa21, aa22, alpha, bb1, bb2, bc, &
    beta, bqq, cfac, cm, cmi, detd, dix, dy, efac, emvpp, emvqq, fy, gam, gpm, &
    hbs, i, ip, isk1, isk2, iv, j, k, lb, lb1, lb2, m, mps, n, nfy, npe, nsol, &
    pc, qc, rpwgpm, rr, sk, sk1, sk2, tz, vc, w1, w2, w3, w4, w5, w6, w7, x, &
    zz

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
      if (nucleus==0) then
        vc(iv-1) = vc(iv-1) + cmi(iv, n)*(v(n)+tz/r(n))
      else
        vc(iv-1) = vc(iv-1) + cmi(iv, n)*v(n)
      end if
    end do
  end do

  do iv = 1, npe
    bc(iv-1) = 0.0d0
    do n = 1, npe
      bc(iv-1) = bc(iv-1) + cmi(iv, n)*b(n)
    end do
  end do



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
! MB
  if (nucleus==0) then
    gam(1) = dsqrt(kap(1)**2-(tz/c)**2)
  else
    gam(1) = dabs(kap(1))
  end if
! MB
  lb(1) = lb1
  sk(1) = sk1
  if (dabs(mj)>l) then
    cg2 = 0.0d0
    cg4 = 0.0d0
    cg8 = 0.0d0
    nsol = 1
    cgd(2) = 0.0d0
    cgo = 0.0d0
    cgmd(2) = 0.0d0
    gam(2) = 0.0d0
    kap(2) = 0.0d0
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
! MB
    if (nucleus==0) then
      gam(2) = dsqrt(kap(2)**2-(tz/c)**2)
    else
      gam(2) = dabs(kap(2))
    end if
! MB
    lb(2) = lb2
    sk(2) = sk2
  end if

  nsolbs = nsol
  ebs = e

  do i = 1, 2
    do j = 1, 2
      do ip = -npemax, mpsmax
        pc(i, j, ip) = c0
        qc(i, j, ip) = c0
      end do
    end do
  end do

! ======================================================================
  if (tz>=2) then

    do j = 1, nsol
      i = 3 - j
      pc(j, j, 0) = dsqrt(abs(kap(j))-gam(j))
      qc(j, j, 0) = (kap(j)+gam(j))*(csqr/tz)*pc(j, j, 0)
      pc(i, j, 0) = c0
      qc(i, j, 0) = c0
    end do

!  determine higher expansion coefficients for the wave functions

    mps = 40

    aa12 = -tz/csqr
    aa21 = tz
    emvqq = (e-vc(0)+csqr)/csqr
    emvpp = -e + vc(0)
    bqq = bc(0)/csqr
!MBA
    if (nucleus==0) then
!MBE

      do j = 1, nsol

        do m = 1, mps
          do i = 1, nsol
            k = 3 - i
            bb1 = (emvqq+bqq*cgmd(i))*qc(i, j, m-1)
            bb2 = (emvpp+bc(0)*cgd(i))*pc(i, j, m-1) + bc(0)*cgo*pc(k, j, m-1)
            do ip = 1, npe - 1
              bb1 = bb1 + (-vc(ip)+bc(ip)*cgmd(i))*qc(i, j, m-1-ip)/csqr
              bb2 = bb2 + (+vc(ip)+bc(ip)*cgd(i))*pc(i, j, m-1-ip) + &
                (+bc(ip)*cgo)*pc(k, j, m-1-ip)
            end do

            aa11 = gam(j) + m + kap(i)
            aa22 = gam(j) + m - kap(i)
            detd = aa11*aa22 - aa12*aa21
            pc(i, j, m) = (bb1*aa22-aa12*bb2)/detd
            qc(i, j, m) = (aa11*bb2-bb1*aa21)/detd
          end do
        end do

      end do
! MBA
    else
! EXPANSION ADAPTED FOR POTENTIALS WITH FINITE NUCLEUS
! EXPANSION OF POTENTIAL UP TO SECOND ORDER: V_O+V_1*R+V_2*R*R
      do j = 1, nsol
        i = 3 - j
        if (kap(j)>0) then
! ARBITRARY STARTING VALUES
          alpha = 0.0d0
          beta = 0.174d0
        else
          beta = 0.0d0
          alpha = 0.174d0
        end if
        pc(j, j, 0) = alpha
        qc(j, j, 0) = beta
        pc(i, j, 0) = 0.0d0
        qc(i, j, 0) = 0.0d0
      end do
      w4 = bc(0)*cgo
      w2 = vc(1)/csqr
      w5 = vc(1)
      w6 = vc(2)/csqr
      w7 = vc(2)
      do j = 1, nsol
        do i = 1, nsol
          w1 = emvqq + bqq*cgmd(i)
          w3 = -emvpp + bc(0)*cgd(i)
          a11 = gam(j) + kap(i) + 1d0
          a12 = gam(j) - kap(i) + 1d0
          if (a11/=0) pc(i, j, 1) = w1/a11*qc(i, j, 0)
          if (a12/=0) qc(i, j, 1) = (-w3*pc(i,j,0)+w4*pc(3-i,j,0))/a12
        end do
      end do
      do j = 1, nsol
        do i = 1, nsol
          w1 = emvqq + bqq*cgmd(i)
          w3 = -emvpp + bc(0)*cgd(i)
          a11 = gam(j) + kap(i) + 2d0
          a12 = gam(j) - kap(i) + 2d0
          if (a11/=0) pc(i, j, 2) = (w1*qc(i,j,1)-w2*qc(i,j,0))/a11
          if (a12/=0) qc(i, j, 2) = (-w3*pc(i,j,1)+w4*pc(3-i,j,1)+w5*pc(i,j,0) &
            )/a12
        end do
      end do
      do j = 1, nsol
        do m = 3, mps
          do i = 1, nsol
            w1 = emvqq + bqq*cgmd(i)
            w3 = -emvpp + bc(0)*cgd(i)
            a21 = gam(j) + kap(i) + dble(m)
            a22 = gam(j) - kap(i) + dble(m)
            if (a21/=0) pc(i, j, m) = (w1*qc(i,j,m-1)-w2*qc(i,j,m-2)-w6*qc(i,j &
              ,m-3))/a21
            if (a22/=0) qc(i, j, m) = (-w3*pc(i,j,m-1)+w4*pc(3-i,j,m-1)+w5*pc( &
              i,j,m-2)+w7*pc(i,j,m-3))/a22
          end do
        end do
      end do
    end if
!MBE



!  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
!  FOR THE FIRST   NABM   R - MESH - POINTS

    do n = 1, nabm
      rr = r(n)

      do j = 1, nsol
        rpwgpm = rr**gam(j)

        do i = 1, nsol
          pr(i, j, n) = pc(i, j, 0)*rpwgpm
          qr(i, j, n) = qc(i, j, 0)*rpwgpm
          dp(i, j, n) = pc(i, j, 0)*rpwgpm*gam(j)*dovr(n)
          dq(i, j, n) = qc(i, j, 0)*rpwgpm*gam(j)*dovr(n)
        end do

        do m = 1, mps
          rpwgpm = rpwgpm*rr
          gpm = gam(j) + m

          do i = 1, nsol
            pr(i, j, n) = pr(i, j, n) + pc(i, j, m)*rpwgpm
            qr(i, j, n) = qr(i, j, n) + qc(i, j, m)*rpwgpm
            dp(i, j, n) = dp(i, j, n) + pc(i, j, m)*rpwgpm*gpm*dovr(n)
            dq(i, j, n) = dq(i, j, n) + qc(i, j, m)*rpwgpm*gpm*dovr(n)
          end do

        end do
      end do
    end do

! ======================================================================
!                                                     == EMPTY SPHERE ==
  else
!                     assume constant pot: V=V(1)   ignore coupling: B=0

    do n = 1, nabm
      zz = sqrt(e-v(1))*r(n)
      efac = (zz/r(n))*c/(e+csqr)

      do j = 1, nsol
        i = 3 - j
        pr(j, j, n) = cjlz(l, zz)*r(n)
        qr(j, j, n) = efac*sk(j)*cjlz(lb(j), zz)*r(n)*c
        dp(j, j, n) = (dble(l+1)*cjlz(l,zz)-zz*cjlz(l+1,zz))*drdi(n)
        m = lb(j)
        dq(j, j, n) = efac*sk(j)*(dble(m+1)*cjlz(m,zz)-zz*cjlz(m+1,zz))* &
          drdi(n)*c

        pr(i, j, n) = c0
        qr(i, j, n) = c0
        dp(i, j, n) = c0
        dq(i, j, n) = c0
      end do
    end do

  end if

! =============================================================== n ====
!     DO 400 J=1,NSOL

  nfy = 0
  do j = 1, nsol
    do i = 1, nsol
      fy(nfy+1) = pr(i, j, 1)
      fy(nfy+2) = qr(i, j, 1)
      dy(nfy+1) = dp(i, j, 1)
      dy(nfy+2) = dq(i, j, 1)
      nfy = nfy + 2
    end do
  end do
  x = 1.0d0
  dix = 1.0d0

  do n = 2, nmesh

    nradbs = n

    call dirbsstp(fy, dy, nfy, x, dix, epsbs, fy, b, v, r, drdi, nmesh)


    nfy = 0
    do j = 1, nsol
      do i = 1, nsol
        pr(i, j, n) = fy(nfy+1)
        qr(i, j, n) = fy(nfy+2)
        nfy = nfy + 2
      end do
    end do

  end do


  if (.not. getirrsol) return


! =============================================================== n ====

! ######################################################################
!                     irregular solution
! ######################################################################

!         calculate the initial values of the wavefunction
!                     at the sphere boundary

  n = nmesh

  zz = pis*r(n)

  do j = 1, nsol
    pi(j, j, n) = cjlz(l, zz)*r(n)
    qi(j, j, n) = cfac*sk(j)*cjlz(lb(j), zz)*r(n)*c
    dp(j, j, n) = (dble(l+1)*cjlz(l,zz)-zz*cjlz(l+1,zz))*drdi(n)

    m = lb(j)
    dq(j, j, n) = cfac*sk(j)*(dble(m+1)*cjlz(m,zz)-zz*cjlz(m+1,zz))*drdi(n)*c

    i = 3 - j
    pi(i, j, n) = c0
    qi(i, j, n) = c0
    dp(i, j, n) = c0
    dq(i, j, n) = c0
  end do

! =============================================================== n ====
  hbs = -1.0d0

  nfy = 0
  do j = 1, nsol
    do i = 1, nsol
      fy(nfy+1) = pi(i, j, nmesh)
      fy(nfy+2) = qi(i, j, nmesh)
      dy(nfy+1) = dp(i, j, nmesh)
      dy(nfy+2) = dq(i, j, nmesh)
      nfy = nfy + 2
    end do
  end do
  x = dble(nmesh)

  nradbs = nmesh
  call dirbsrad(x, fy, dy, drdi, b, v, r, nmesh)

  do n = nmesh - 1, 1, -1
    nradbs = n

    call dirbsstp(fy, dy, nfy, x, hbs, epsbs, fy, b, v, r, drdi, nmesh)

    nfy = 0
    do j = 1, nsol
      do i = 1, nsol
        pi(i, j, n) = fy(nfy+1)
        qi(i, j, n) = fy(nfy+2)
!              DP(I,J,N) = DY(NFY+1)
!              DQ(I,J,N) = DY(NFY+2)
        nfy = nfy + 2
      end do
    end do

  end do

! =============================================================== n ====
end subroutine
