!------------------------------------------------------------------------------------
!> Summary: Generates the lattice vectors of direct and reciprocal space from basic translation vectors for a 2D system
!> Author: V. Popescu
!> Generates the lattice vectors of direct and reciprocal space from
!> basic translation vectors for a 2D system
!------------------------------------------------------------------------------------
module mod_lattice2d
  use :: mod_datatypes, only: dp
  use :: constants, only: pi
  private :: dp

contains
  !-------------------------------------------------------------------------------
  !> Summary: Generates the lattice vectors of direct and reciprocal space from basic translation vectors for a 2D system
  !> Author: V. Popescu
  !> Category: geometry, k-points, electrostatics, KKRhost 
  !> Deprecated: False 
  !> Generates the lattice vectors of direct and reciprocal space from basic 
  !> translation vectors for a 2D system
  !-------------------------------------------------------------------------------
  !> @note Popescu May 2004: The routine has been brought to a form which is very similar to
  !> `LATTICE2D` -- from which it has been originally derived. Dimension of arrays GN,RM
  !> changed from `(4,*)` to `(2,*)`, the 4th one it is used only locally (GNR/RMR)
  !> -- only `GN/RM(2,*)` are actually needed in EWALD2D
  !> @endnote
  !------------------------------------------------------------------------------- 
  subroutine lattice2d(alat,bravais,recbv,ngmax,nrmax,nshlg,nshlr,nsg,nsr,gn,rm,    &
    rmax,gmax,iprint,nmaxd,ishld)

    implicit none
    ! ..
    ! .. Input variables
    integer, intent (in) :: nmaxd  !! Paremeters for the Ewald summations
    integer, intent (in) :: ishld  !! Paremeters for the Ewald summations
    integer, intent (in) :: iprint !! Printing index control
    real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
    real (kind=dp), dimension (3, 3), intent (in) :: recbv !! Reciprocal basis vectors
    real (kind=dp), dimension (3, 3), intent (in) :: bravais !! Bravais lattice vectors
    ! ..
    ! .. Input/Output variables
    real (kind=dp), intent (inout) :: gmax !! Ewald summation cutoff parameter for reciprocal space summation
    real (kind=dp), intent (inout) :: rmax !! Ewald summation cutoff parameter for real space summation
    integer, dimension (ishld), intent (inout) :: nsg
    integer, dimension (ishld), intent (inout) :: nsr
    real (kind=dp), dimension (2, nmaxd), intent (inout) :: gn !! x,y,z of reciprocal lattice vectors
    real (kind=dp), dimension (2, nmaxd), intent (inout) :: rm !! x,y,z of real space vectors
    ! .. Output variables
    integer, intent (out) :: nshlr !! Shells in real space
    integer, intent (out) :: nshlg !! Shells in reciprocal space
    integer, intent (out) :: nrmax !! Number of real space vectors rr
    integer, intent (out) :: ngmax !! Number of reciprocal space vectors
    ! ..
    ! .. Local scalars ..
    integer :: i, k, l, m, n, n1, ng, nr, nsh, nshl, numg, numgh, numr, numrh
    real (kind=dp) :: rx, ry, vmin
    real (kind=dp) :: a, absgm, absrm, ag, ar, b, da, db, gx, gy
    ! ..
    ! .. Local arrays ..
    real (kind=dp), dimension (nmaxd) :: gnr
    real (kind=dp), dimension (nmaxd) :: rmr
    real (kind=dp), dimension (3) :: absg
    real (kind=dp), dimension (3) :: absr
    real (kind=dp), dimension (nmaxd) :: length
    real (kind=dp), dimension (3, 3) :: bg
    real (kind=dp), dimension (3, 3) :: br
    real (kind=dp), dimension (4, nmaxd) :: cj
    ! ----------------------------------------------------------------------------
    ! OUTPUT
    ! ----------------------------------------------------------------------------
    write (1337, '(5X,2A,/)') '< LATTICE2D > : ', 'generating direct/reciprocal lattice vectors'
    ! ----------------------------------------------------------------------------
    ! OUTPUT
    ! ----------------------------------------------------------------------------
    rmax = rmax*alat
    gmax = gmax/alat
    ! ----------------------------------------------------------------------------
    ! OUTPUT
    ! ----------------------------------------------------------------------------
    write (1337, fmt=100) rmax, gmax
    ! ----------------------------------------------------------------------------
    ! OUTPUT
    ! ----------------------------------------------------------------------------

    ! ----------------------------------------------------------------------------
    ! Basic trans. vectors and basis vectors
    ! ----------------------------------------------------------------------------
    do i = 1, 3
      br(1, i) = bravais(1, i)*alat
      br(2, i) = bravais(2, i)*alat
    end do
    ! ----------------------------------------------------------------------------
    ! Generate primitive vectors BG of reciprocal space
    ! ----------------------------------------------------------------------------
    do i = 1, 3
      bg(1, i) = recbv(1, i)*2e0_dp*pi/alat
      bg(2, i) = recbv(2, i)*2e0_dp*pi/alat
    end do
    ! ----------------------------------------------------------------------------
    ! Estimate no. of lattice vectors
    ! ----------------------------------------------------------------------------
    do i = 1, 3
      absr(i) = sqrt(br(1,i)**2+br(2,i)**2)
      absg(i) = sqrt(bg(1,i)**2+bg(2,i)**2)
    end do

    absrm = max(absr(1), absr(2))
    absgm = max(absg(1), absg(2))
    absrm = 2.0e0_dp*pi/absrm
    absgm = 2.0e0_dp*pi/absgm
    numr = 2*(int(rmax/absgm)+1) + 1
    numg = 2*(int(gmax/absrm)+1) + 1
    numrh = numr/2 + 1
    numgh = numg/2 + 1

    !--------------------------------------------------------------------------------
    ! generate lattice vectors of real space
    !--------------------------------------------------------------------------------

    write (1337, *) 'Real space...'
    nr = 0
    ! ----------------------------------------------------------------------------
    do l = 1, numr
      a = real(l-numrh, kind=dp)
      do m = 1, numr
        b = real(m-numrh, kind=dp)
        ! ----------------------------------------------------------------------
        rx = a*br(1, 1) + b*br(1, 2)
        ry = a*br(2, 1) + b*br(2, 2)
        ar = sqrt(rx*rx+ry*ry)
        ! ----------------------------------------------------------------------
        if (ar<=rmax) then
          nr = nr + 1
          if (nr>nmaxd) then
            write (6, *) 'lattice2d: ERROR: Dimension NMAXD in the inputcard too small', nr, nmaxd
            stop 'lattice2d'
          end if
          cj(1, nr) = rx
          cj(2, nr) = ry
          cj(3, nr) = 0e0_dp
          cj(4, nr) = ar
        end if
      end do
    end do
    ! ----------------------------------------------------------------------------
    nrmax = nr
    ! ----------------------------------------------------------------------------
    ! Sort vectors in order of increasing absolute value
    ! ----------------------------------------------------------------------------
    write (1337, fmt='(A11,I8,A11)') '...sorting ', nrmax, ' vectors...'

    da = 1.e-06_dp
    nsh = 0
    nshl = -1
    ! ----------------------------------------------------------------------------
    do k = 1, nr
      vmin = rmax + 1.0e0_dp
      do n = 1, nr
        if (cj(4,n)-vmin<0e0_dp) then
          vmin = cj(4, n)
          n1 = n
        end if
      end do

      nshl = nshl + 1
      rm(1, k) = cj(1, n1)
      rm(2, k) = cj(2, n1)
      rmr(k) = cj(4, n1)
      db = vmin
      ! -------------------------------------------------------------------------
      if (db>da+1.e-06_dp) then
        nsh = nsh + 1
        if (nsh>ishld) then
          write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', nsh, ishld
          stop 'lattice2d'
        end if

        nsr(nsh) = nshl
        nshl = 0
        da = db
      end if
      ! -------------------------------------------------------------------------
      cj(4, n1) = rmax + 1.0e0_dp
    end do
    ! ----------------------------------------------------------------------------
    nsh = nsh + 1
    nshl = nshl + 1
    if (nsh>ishld) then
      write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', nsh, ishld
      stop 'lattice2d'
    end if

    nsr(nsh) = nshl
    nshlr = nsh
    if (nshlr<=1) stop 'lattice2d: ERROR: cut-off radius RMAX too small '
    write (1337, *) '...done.'
    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! generate lattice vectors of real space
    !--------------------------------------------------------------------------------

    write (1337, *) 'Reciprocal space...'
    ng = 0
    ! ----------------------------------------------------------------------------
    do l = 1, numg
      a = real(l-numgh, kind=dp)
      do m = 1, numg
        b = real(m-numgh, kind=dp)
        ! ----------------------------------------------------------------------
        gx = a*bg(1, 1) + b*bg(1, 2)
        gy = a*bg(2, 1) + b*bg(2, 2)
        ag = sqrt(gx*gx+gy*gy)
        ! ----------------------------------------------------------------------
        if (ag<=gmax) then
          ng = ng + 1
          if (ng>nmaxd) then
            write (6, *) ' ERROR: Dimension NMAXD in the inputcard too small', ng, nmaxd
            stop 'lattice2d'
          end if
          cj(1, ng) = gx
          cj(2, ng) = gy
          cj(3, ng) = 0e0_dp
          cj(4, ng) = ag
        end if
      end do
    end do
    ! ----------------------------------------------------------------------------
    ngmax = ng
    ! ----------------------------------------------------------------------------

    ! ----------------------------------------------------------------------------
    ! Sort vectors in order of increasing abs. value
    ! ----------------------------------------------------------------------------
    write (1337, fmt='(A11,I8,A11)') '...sorting ', ngmax, ' vectors...'
    do n = 1, ng
      length(n) = cj(4, n)
    end do
    da = 1.e-06_dp
    nsh = 0
    nshl = -1
    ! ----------------------------------------------------------------------------
    do k = 1, ng
      vmin = gmax + 1.0e0_dp
      do n = 1, ng
        if (length(n)<vmin) then   ! ( CJ(4,N).LT.VMIN ) THEN
          vmin = length(n)         ! CJ(4,N)
          n1 = n
        end if
      end do

      nshl = nshl + 1
      gn(1, k) = cj(1, n1)
      gn(2, k) = cj(2, n1)
      gnr(k) = length(n1)          ! CJ(4,N1)
      db = vmin
      ! -------------------------------------------------------------------------
      if (db>da+1.e-07_dp) then
        nsh = nsh + 1              ! Number of shells of different length
        if (nsh>ishld) then
          write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', nsh, ishld
          stop 'lattice2d'
        end if

        nsg(nsh) = nshl            ! Number of vectors in shell
        nshl = 0
        da = db
      end if
      ! -------------------------------------------------------------------------
      length(n1) = gmax + 1.0e0_dp ! CJ(4,N1) = GMAX + 1.0D0
    end do
    ! ----------------------------------------------------------------------------
    nsh = nsh + 1
    nshl = nshl + 1
    if (nsh>ishld) then
      write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', nsh, ishld
      stop 'lattice2d'
    end if

    nsg(nsh) = nshl
    nshlg = nsh
    if (nshlg<=1) stop 'lattice2dERROR: cut-off radius GMAX too small '

    write (1337, *) '...done.'
    !--------------------------------------------------------------------------------
    ! OUTPUT
    !--------------------------------------------------------------------------------
    write (1337, fmt=110)
    write (1337, fmt=120) 'Direct  lattice', nrmax, nshlr, rmr(nrmax)
    write (1337, fmt=120) 'Recipr. lattice', ngmax, nshlg, gnr(ngmax)
    write (1337, fmt=130)

    if (iprint<3) return
    !--------------------------------------------------------------------------------
    ! OUTPUT
    !--------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    k = 0
    write (1337, fmt=140) 'real-space'
    do l = 1, nshlr
      write (1337, 150) l, nsr(l), rmr(k+1), (rm(m,k+1), m=1, 2)
      do n = 2, nsr(l)
        write (1337, fmt=160)(rm(m,k+n), m=1, 2)
      end do
      if (l/=nshlr) write (1337, 170)
      k = k + nsr(l)
    end do
    write (1337, 180)
    k = 0
    write (1337, fmt=140) 'reciprocal'
    do l = 1, nshlg
      write (1337, 150) l, nsg(l), gnr(k+1), (gn(m,k+1), m=1, 2)
      do n = 2, nsg(l)
        write (1337, fmt=160)(gn(m,k+n), m=1, 2)
      end do
      if (l/=nshlg) write (1337, 170)
      k = k + nsg(l)
    end do
    write (1337, 180)
    ! ----------------------------------------------------------------------------

100 format (10x, 'R max =', f10.5, ' (a.u.)', /, 10x, 'G max =', f10.5, ' (1/a.u.)', /)
110 format (10x, '               vectors  shells  max. R ', /, 10x, '               ------------------------------')
120 format (10x, a, i7, 2x, i6, 2x, f9.5)
130 format (10x, '               ------------------------------', /)
140 format (10x, 45('+'), /, 13x, 'generated ', a, ' lattice vectors', /, 10x, 45('+'), /, 10x, 'shell Nvec    radius          x         y', /, 10x, 45('-'))
150 format (10x, i5, i5, f12.6, 2x, 2f10.5)
160 format (34x, 2f10.5)
170 format (13x, 42('-'))
180 format (10x, 45('+'), /)
  end subroutine lattice2d

end module mod_lattice2d