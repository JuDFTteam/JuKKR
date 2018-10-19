!------------------------------------------------------------------------------------
!> Summary: Generates the lattice vectors of direct and reciprocal space frombasic translation vectors for a 3D system
!> Author: V. Popescu 
!> Generates the lattice vectors of direct and reciprocal space from
!> basic translation vectors for a 3D system
!------------------------------------------------------------------------------------
module mod_lattice3d
  use :: mod_datatypes, only: dp
  use :: constants, only: pi
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Generates the lattice vectors of direct and reciprocal space from basic translation vectors for a 3D system
  !> Author: V. Popescu
  !> Category: geometry, k-points, electrostatics, KKRhost 
  !> Deprecated: False 
  !> Generates the lattice vectors of direct and reciprocal space from basic 
  !> translation vectors for a 3D system
  !-------------------------------------------------------------------------------
  !> @note Popescu May 2004: The routine has been brought to a form which is very similar to
  !> `LATTICE2D` -- from which it has been originally derived. Dimension of arrays GN,RM
  !> changed from `(4,*)` to `(2,*)`, the 4th one it is used only locally (GNR/RMR)
  !> @endnote
  !------------------------------------------------------------------------------- 
  subroutine lattice3d(alat,bravais,recbv,ngmax,nrmax,nshlg,nshlr,nsg,nsr,gn,rm,    &
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
    real (kind=dp), intent (inout) :: rmax !! Ewald summation cutoff parameter for real space summation
    real (kind=dp), intent (inout) :: gmax !! Ewald summation cutoff parameter for reciprocal space summation
    ! ..
    ! .. Input/Output variables
    integer, dimension (ishld), intent (inout) :: nsg
    integer, dimension (ishld), intent (inout) :: nsr
    real (kind=dp), dimension (3, nmaxd), intent (inout) :: gn !! x,y,z   of
    ! reciprocal
    ! lattice vectors
    real (kind=dp), dimension (3, nmaxd), intent (inout) :: rm !! x,y,z  of
    ! real space
    ! vectors
    ! .. Ouptut Variables
    integer, intent (out) :: ngmax !! Number of reciprocal space vectors
    integer, intent (out) :: nrmax !! Number of real space vectors rr
    integer, intent (out) :: nshlg !! Shells in reciprocal space
    integer, intent (out) :: nshlr !! Shells in real space
    ! ..
    ! .. Local scalars ..
    integer :: i, k, l, m, n, n1, ng, nr, nsh, nshl, numg, numgh, numr, numrh
    real (kind=dp) :: a, absgm, absrm, ag, ar, b, c, da, db, gx, gy, gz, rx, ry, rz, vmin
    ! ..
    ! .. Local arrays ..
    real (kind=dp), dimension (nmaxd) :: gnr
    real (kind=dp), dimension (nmaxd) :: rmr
    real (kind=dp), dimension (3) :: absg
    real (kind=dp), dimension (3) :: absr
    real (kind=dp), dimension (3, 3) :: bg
    real (kind=dp), dimension (3, 3) :: br
    real (kind=dp), dimension (4, nmaxd) :: cj
    ! ----------------------------------------------------------------------------
    ! OUTPUT
    ! ----------------------------------------------------------------------------
    write (1337, '(5X,2A,/)') '< LATTICE3D > : ', 'generating direct/reciprocal lattice vectors'
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
      br(3, i) = bravais(3, i)*alat
    end do
    ! ----------------------------------------------------------------------------
    ! Generate primitive vectors BG of reciprocal space
    ! ----------------------------------------------------------------------------
    do i = 1, 3
      bg(1, i) = recbv(1, i)*2e0_dp*pi/alat
      bg(2, i) = recbv(2, i)*2e0_dp*pi/alat
      bg(3, i) = recbv(3, i)*2e0_dp*pi/alat
    end do
    ! ----------------------------------------------------------------------------
    ! Estimate no. of lattice vectors
    ! ----------------------------------------------------------------------------
    do i = 1, 3
      absr(i) = sqrt(br(1,i)**2+br(2,i)**2+br(3,i)**2)
      absg(i) = sqrt(bg(1,i)**2+bg(2,i)**2+bg(3,i)**2)
    end do

    absrm = max(absr(1), absr(2), absr(3))
    absgm = max(absg(1), absg(2), absg(3))
    absrm = 2.0e0_dp*pi/absrm
    absgm = 2.0e0_dp*pi/absgm
    numr = 2*(int(rmax/absgm)+1) + 1
    numg = 2*(int(gmax/absrm)+1) + 1
    numrh = numr/2 + 1
    numgh = numg/2 + 1

    !--------------------------------------------------------------------------------
    ! generate lattice vectors of real space
    !--------------------------------------------------------------------------------
    nr = 0
    ! ----------------------------------------------------------------------------
    do l = 1, numr
      a = real(l-numrh, kind=dp)
      do m = 1, numr
        b = real(m-numrh, kind=dp)
        do n = 1, numr
          c = real(n-numrh, kind=dp)
          ! -------------------------------------------------------------------
          rx = a*br(1, 1) + b*br(1, 2) + c*br(1, 3)
          ry = a*br(2, 1) + b*br(2, 2) + c*br(2, 3)
          rz = a*br(3, 1) + b*br(3, 2) + c*br(3, 3)
          ar = sqrt(rx*rx+ry*ry+rz*rz)
          ! -------------------------------------------------------------------
          if (ar<=rmax) then
            nr = nr + 1
            if (nr>nmaxd) then
              write (6, *) ' ERROR: Dimension NMAXD in the inputcard is too small', nr, nmaxd
              stop
            end if
            cj(1, nr) = rx
            cj(2, nr) = ry
            cj(3, nr) = rz
            cj(4, nr) = ar
          end if
        end do
      end do
    end do
    ! ----------------------------------------------------------------------------
    nrmax = nr
    ! ----------------------------------------------------------------------------
    ! Sort vectors in order of increasing absolute value
    ! ----------------------------------------------------------------------------
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
      rm(3, k) = cj(3, n1)
      rmr(k) = cj(4, n1)
      db = vmin
      ! -------------------------------------------------------------------------
      if (db>da+1.e-06_dp) then
        nsh = nsh + 1
        if (nsh>ishld) then
          write (6, *) ' ERROR: Dimension ISHLD in the inputcard is too small', nsh, ishld
          stop
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
      write (6, *) ' ERROR: Dimension ISHLD in the inputcard is too small', nsh, ishld
      stop
    end if

    nsr(nsh) = nshl
    nshlr = nsh
    if (nshlr<=1) stop ' ERROR: cut-off radius RMAX too small '
    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! generate lattice vectors of real space
    !--------------------------------------------------------------------------------
    ng = 0
    ! ----------------------------------------------------------------------------
    do l = 1, numg
      a = real(l-numgh, kind=dp)
      do m = 1, numg
        b = real(m-numgh, kind=dp)
        do n = 1, numg
          c = real(n-numgh, kind=dp)
          ! -------------------------------------------------------------------
          gx = a*bg(1, 1) + b*bg(1, 2) + c*bg(1, 3)
          gy = a*bg(2, 1) + b*bg(2, 2) + c*bg(2, 3)
          gz = a*bg(3, 1) + b*bg(3, 2) + c*bg(3, 3)
          ag = sqrt(gx*gx+gy*gy+gz*gz)
          ! -------------------------------------------------------------------
          if (ag<=gmax) then
            ng = ng + 1
            if (ng>nmaxd) then
              write (6, *) ' ERROR: Dimension NMAXD in the inputcard is too small', ng, nmaxd
              stop
            end if
            cj(1, ng) = gx
            cj(2, ng) = gy
            cj(3, ng) = gz
            cj(4, ng) = ag
          end if
        end do
      end do
    end do
    ! ----------------------------------------------------------------------------
    ngmax = ng
    ! ----------------------------------------------------------------------------

    ! ----------------------------------------------------------------------------
    ! Sort vectors in order of increasing abs. value
    ! ----------------------------------------------------------------------------
    da = 1.e-06_dp
    nsh = 0
    nshl = -1
    ! ----------------------------------------------------------------------------
    do k = 1, ng
      vmin = gmax + 1.0e0_dp
      do n = 1, ng
        if (cj(4,n)<vmin) then
          vmin = cj(4, n)
          n1 = n
        end if
      end do

      nshl = nshl + 1
      gn(1, k) = cj(1, n1)
      gn(2, k) = cj(2, n1)
      gn(3, k) = cj(3, n1)
      gnr(k) = cj(4, n1)
      db = vmin
      ! -------------------------------------------------------------------------
      if (db>da+1.e-07_dp) then
        nsh = nsh + 1
        if (nsh>ishld) then
          write (6, *) ' ERROR: Dimension ISHLD in the inputcard is too small', nsh, ishld
          stop
        end if

        nsg(nsh) = nshl
        nshl = 0
        da = db
      end if
      ! -------------------------------------------------------------------------
      cj(4, n1) = gmax + 1.0e0_dp
    end do
    ! ----------------------------------------------------------------------------
    nsh = nsh + 1
    nshl = nshl + 1
    if (nsh>ishld) then
      write (6, *) ' ERROR: Dimension ISHLD in the inputcard is too small', nsh, ishld
      stop
    end if

    nsg(nsh) = nshl
    nshlg = nsh
    if (nshlg<=1) stop ' ERROR: cut-off radius GMAX too small '
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
    k = 0
    write (1337, fmt=140) 'real-space'
    do l = 1, nshlr
      write (1337, 150) l, nsr(l), rmr(k+1), (rm(m,k+1), m=1, 3)
      do n = 2, nsr(l)
        write (1337, fmt=160)(rm(m,k+n), m=1, 3)
      end do
      if (l/=nshlr) write (1337, 170)
      k = k + nsr(l)
    end do
    write (1337, 180)
    k = 0
    write (1337, fmt=140) 'reciprocal'
    do l = 1, nshlg
      write (1337, 150) l, nsg(l), gnr(k+1), (gn(m,k+1), m=1, 3)
      do n = 2, nsg(l)
        write (1337, fmt=160)(gn(m,k+n), m=1, 3)
      end do
      if (l/=nshlg) write (1337, 170)
      k = k + nsg(l)
    end do
    write (1337, 180)
    ! ----------------------------------------------------------------------------

100 format (10x, 'R max =', f9.5, ' (a.u.)', /, 10x, 'G max =', f9.5, ' (1/a.u.)', /)
110 format (10x, '               vectors  shells  max. R ', /, 10x, '               ------------------------------')
120 format (10x, a, i7, 2x, i6, 2x, f9.5)
130 format (10x, '               ------------------------------', /)
140 format (10x, 55('+'), /, 18x, 'generated ', a, ' lattice vectors', /, 10x, 55('+'), /, 10x, 'shell Nvec    radius          x         y         z', /, 10x, 55('-'))
150 format (10x, i5, i5, f12.6, 2x, 3f10.5)
160 format (34x, 3f10.5)
170 format (13x, 52('-'))
180 format (10x, 55('+'), /)
  end subroutine lattice3d

end module mod_lattice3d
