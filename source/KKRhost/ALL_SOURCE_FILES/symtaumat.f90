module mod_symtaumat
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine symtaumat(rotname, rotmat, drot, nsym, isymindex, symunitary, nqmax, nkmmax, nq, nl, krel, iprint, nsymaxd)
    ! ********************************************************************
    ! *                                                                  *
    ! *  Find the symmetry matrices DROT that act on t, tau, ....        *
    ! *  KREL=0: for real spherical harmonics                            *
    ! *  KREL=1: for relativistic represntation                          *
    ! *                                                                  *
    ! *  The NSYM allowed symmetry operations are indicated by ISYMINDEX *
    ! *  in the table  ROTMAT. For KREL=1, SYMUNITARY=T/F indicates a    *
    ! *  unitary/antiunitary symmetry operation.                         *
    ! *                                                                  *
    ! *  The routine determines first the Euler angles correponding      *
    ! *  to a symmetry operation. Reflections are decomposed into        *
    ! *  inversion + rotation for this reason.                           *
    ! *                                                                  *
    ! ********************************************************************

    use :: mod_calcrotmat
    use :: mod_checkrmat
    use :: mod_cinit
    use :: mod_cmatstr
    use :: mod_ddet33
    use :: mod_taustruct
    use :: mod_errortrap
    implicit none

    ! PARAMETER definitions
    complex (kind=dp) :: ci, c1, c0
    parameter (ci=(0.0e0_dp,1.0e0_dp), c1=(1.0e0_dp,0.0e0_dp), c0=(0.0e0_dp,0.0e0_dp))

    ! Dummy arguments
    integer :: iprint, krel, nkmmax, nl, nq, nqmax, nsym, nsymaxd
    complex (kind=dp) :: drot(nkmmax, nkmmax, 48)
    integer :: isymindex(nsymaxd)
    real (kind=dp) :: rotmat(64, 3, 3)
    character (len=10) :: rotname(64)
    logical :: symunitary(48)

    ! Local variables
    real (kind=dp) :: a, b, co1, co2, co3, det, fact(0:100), pi, si1, si2, si3, sk, symeulang(3, 48), tet1, tet2, tet3, rj, rmj
    real (kind=dp) :: dble
    complex (kind=dp) :: dinv(nkmmax, nkmmax), dtim(nkmmax, nkmmax), rc(nkmmax, nkmmax), w1(nkmmax, nkmmax), w2(nkmmax, nkmmax)
    logical :: equal
    integer :: i, i1, i2, ind0q(nqmax), invflag(48), iq, irel, ireleff, isym, itop, j, k, l, loop, m, n, nk, nkeff, nkm, nlm, nok, ns
    integer :: nint
    real (kind=dp) :: rmat(3, 3)
    real (kind=dp) :: w

    equal(a, b) = (abs(a-b)<1e-7_dp)

    write (1337, 110)

    pi = 4e0_dp*atan(1e0_dp)

    irel = krel*3
    nk = (1-krel)*nl + krel*(2*nl-1)
    nkm = (1+krel)*nl**2

    ! -----------------------------------------------------------------------
    fact(0) = 1.0e0_dp
    do i = 1, 100
      fact(i) = fact(i-1)*dble(i)
    end do
    ! -----------------------------------------------------------------------

    ind0q(1) = 0
    do iq = 2, nq
      ind0q(iq) = ind0q(iq-1) + nkm
    end do

    ! ----------------------------------------------------------------------
    ! RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
    ! |LC> = sum[LR] |LR> * RC(LR,LC)
    ! ----------------------------------------------------------------------
    if (krel==0) then
      nlm = nkm

      call cinit(nkmmax*nkmmax, rc)

      w = 1.0e0_dp/sqrt(2.0e0_dp)

      do l = 0, (nl-1)
        do m = -l, l
          i = l*(l+1) + m + 1
          j = l*(l+1) - m + 1

          if (m<0) then
            rc(i, i) = -ci*w
            rc(j, i) = w
          end if
          if (m==0) then
            rc(i, i) = c1
          end if
          if (m>0) then
            rc(i, i) = w*(-1.0e0_dp)**m
            rc(j, i) = ci*w*(-1.0e0_dp)**m
          end if
        end do
      end do
    end if

    ! =======================================================================
    ! The routine determines first the Euler angles correponding
    ! to a symmetry operation. Reflections are decomposed into
    ! inversion + rotation for this reason.
    ! =======================================================================

    do isym = 1, nsym

      do i1 = 1, 3
        do i2 = 1, 3
          rmat(i1, i2) = rotmat(isymindex(isym), i1, i2)
        end do
      end do

      det = ddet33(rmat)

      invflag(isym) = 0
      if (det<0e0_dp) then
        call dscal(9, -1.0e0_dp, rmat, 1)
        invflag(isym) = 1
      end if

      ! ----------------------------------------------------------------------
      co2 = rmat(3, 3)
      tet2 = acos(co2)
      loop = 0
100   continue
      if (loop==1) tet2 = -tet2
      si2 = sin(tet2)

      if (equal(co2,1.0e0_dp)) then
        tet1 = acos(rmat(1,1))
        if (.not. equal(rmat(1,2),sin(tet1))) then
          tet1 = -tet1
          if (.not. equal(rmat(1,2),sin(tet1))) write (1337, *) '>>>>>>>>>>>>>>> STRANGE 1'
        end if
        tet2 = 0e0_dp
        tet3 = 0e0_dp
      else if (equal(co2,-1e0_dp)) then
        tet1 = acos(-rmat(1,1))
        if (.not. equal(rmat(1,2),-sin(tet1))) then
          tet1 = -tet1
          if (.not. equal(rmat(1,2),-sin(tet1))) write (1337, *) '>>>>>>>>>>>>>>> STRANGE 2'
        end if
        tet2 = pi
        tet3 = 0e0_dp
      else
        tet1 = acos(rmat(3,1)/si2)
        if (.not. equal(rmat(3,2),si2*sin(tet1))) then
          tet1 = -tet1
          if (.not. equal(rmat(3,2),si2*sin(tet1))) write (1337, *) '>>>>>>>>>>>>>>> STRANGE 3'
        end if

        tet3 = acos(-rmat(1,3)/si2)
        if (.not. equal(rmat(2,3),si2*sin(tet3))) then
          tet3 = -tet3
          if (.not. equal(rmat(2,3),si2*sin(tet3))) write (1337, *) '>>>>>>>>>>>>>>> STRANGE 4'
        end if

      end if

      co1 = cos(tet1)
      si1 = sin(tet1)
      co2 = cos(tet2)
      si2 = sin(tet2)
      co3 = cos(tet3)
      si3 = sin(tet3)

      nok = 0
      do i1 = 1, 3
        do i2 = 1, 3
          if (checkrmat(rmat,co1,si1,co2,si2,co3,si3,i1,i2)) then
            nok = nok + 1
          else if (loop<1) then
            loop = loop + 1
            go to 100
          end if
        end do
      end do

      symeulang(1, isym) = tet1*(180e0_dp/pi)
      symeulang(2, isym) = tet2*(180e0_dp/pi)
      symeulang(3, isym) = tet3*(180e0_dp/pi)

      if (nok/=9) write (1337, 130) nok
      write (1337, 120) isym, rotname(isymindex(isym)), invflag(isym), (symeulang(i,isym), i=1, 3), symunitary(isym)

    end do
    write (1337, '(8X,57("-"),/)')

    ! -----------------------------------------------------------------------
    ! initialize all rotation matrices
    ! -----------------------------------------------------------------------

    call cinit(nkmmax*nkmmax*nsym, drot)

    ! -----------------------------------------------------------------------
    ! create rotation matrices
    ! -----------------------------------------------------------------------

    if (irel<=2) then
      ireleff = 0
      nkeff = nl
    else
      ireleff = 3
      nkeff = nk
    end if

    do isym = 1, nsym

      call calcrotmat(nkeff, ireleff, symeulang(1,isym), symeulang(2,isym), symeulang(3,isym), drot(1,1,isym), fact, nkmmax)

    end do
    ! -----------------------------------------------------------------------
    ! create matrix for inversion
    ! -----------------------------------------------------------------------
    call cinit(nkmmax*nkmmax, dinv)

    i = 0
    if (irel>2) then
      ns = 2
    else
      ns = 1
    end if
    do l = 0, (nl-1)
      do m = 1, ns*(2*l+1)
        i = i + 1
        dinv(i, i) = (-1.0e0_dp)**l
      end do
    end do
    itop = i

    ! -----------------------------------------------------------------------
    ! include inversion
    ! -----------------------------------------------------------------------
    do isym = 1, nsym
      if (invflag(isym)/=0) then

        call zgemm('N', 'N', nkm, nkm, nkm, c1, drot(1,1,isym), nkmmax, dinv, nkmmax, c0, w2, nkmmax)

        do j = 1, nkm
          call zcopy(nkm, w2(1,j), 1, drot(1,j,isym), 1)
        end do
      end if
    end do

    ! -----------------------------------------------------------------------
    ! add second spin-diagonal block for  IREL=2
    ! spin off-diagonal blocks have been initialized before
    ! -----------------------------------------------------------------------
    if (irel==2) then
      nlm = nkm/2
      if (itop/=nlm) call errortrap('SYMTAUMAT', 11, 1)
      do isym = 1, nsym

        do j = 1, nlm
          call zcopy(nlm, drot(1,j,isym), 1, drot(nlm+1,nlm+j,isym), 1)
        end do
      end do
    end if
    ! -----------------------------------------------------------------------
    ! transform to real spherical representation for  KREL=0
    ! -----------------------------------------------------------------------
    n = nkm
    m = nkmmax
    if (krel==0) then
      do isym = 1, nsym
        call zgemm('N', 'N', n, n, n, c1, rc, m, drot(1,1,isym), m, c0, w1, m)
        call zgemm('N', 'C', n, n, n, c1, w1, m, rc, m, c0, drot(1,1,isym), m)
      end do
    end if
    ! -----------------------------------------------------------------------
    ! create matrix for time reversal
    ! -----------------------------------------------------------------------
    if (irel>1) then

      call cinit(nkmmax*nkmmax, dtim)

      i = 0
      do k = 1, nk
        l = k/2
        if (l*2==k) then
          sk = -1e0_dp
        else
          sk = +1e0_dp
        end if
        rj = l + sk*0.5_dp
        do rmj = -rj, +rj
          i1 = nint(2*l*(rj+0.5e0_dp)+rj+rmj+1)
          i2 = nint(2*l*(rj+0.5e0_dp)+rj-rmj+1)
          dtim(i1, i2) = sk*(-1)**nint(rmj+0.5e0_dp)
        end do
      end do
      if (iprint>0) then
        call cmatstr('Inversion     MATRIX', 20, dinv, nkm, nkmmax, 3, 3, 0, 1e-8_dp, 6)
        call cmatstr('Time reversal MATRIX', 20, dtim, nkm, nkmmax, 3, 3, 0, 1e-8_dp, 6)
      end if

    end if
    ! =======================================================================
    ! set up of transformation matrices completed
    ! =======================================================================

    ! =======================================================================
    ! include time reversal operation for anti-unitary transformations
    ! =======================================================================
    do isym = 1, nsym
      if (.not. symunitary(isym)) then
        if (irel==2) call errortrap('SYMTAUMAT', 14, 1)

        call zgemm('N', 'N', nkm, nkm, nkm, c1, drot(1,1,isym), nkmmax, dtim, nkmmax, c0, w2, nkmmax)
        do j = 1, nkm
          call zcopy(nkm, w2(1,j), 1, drot(1,j,isym), 1)
        end do
      end if
    end do

    ! -----------------------------------------------------------------------
    ! for testing

    ! ccc      write (6,*) ' NUMBER OF SYMMETRIES : ', NSYM
    ! ccc
    ! ccc      do isym = 1,nsym
    ! ccc         write(6,*) ' ISYM = ',isym
    ! ccc         call cmatstr('DROT',4,drot(1,1,isym),nkm,nkmmax,krel*3,krel*3,
    ! ccc     &        0,1d-12,6)
    ! ccc         write(6,*)
    ! ccc      end do

    ! -----------------------------------------------------------------------

    if (iprint==0) return

    ! =======================================================================
    ! find the structure of the site-diagonal TAU - matrices  TAUQ
    ! =======================================================================

    call taustruct(drot, nsym, symunitary, nkm, nq, nqmax, nkmmax, iprint, irel)

    return
110 format (5x, '<SYMTAUMAT> : rotation matrices acting on t/G/tau', /, /, 8x, 57('-'), /, 8x, 'ISYM            INV          Euler angles      Unitarity', /, 8x, 57('-'))
120 format (8x, i2, 3x, a, i3, 3f10.5, 3x, l1)
130 format (50('>'), ' trouble in <SYMTAUMAT>', i3, f10.5)
  end subroutine symtaumat

end module mod_symtaumat
