  subroutine baryint_gf2(ne,ze,ia,ja,z0,n0,gf,onsite,struct)
! Barycentric interpolation with n0 points
! Works for any point distribution, as long as n0 is kept low
! See JP Berrut and LN Trefethen, SIAM Review vol. 46, pp. 501
  use global

  implicit none

! Number of energies
  integer(kind=i4b), intent(in)  :: ne
! List of energies available for interpolation
  complex(kind=c8b), intent(in)  :: ze(ne)
! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! Desired point
  complex(kind=c8b), intent(in)  :: z0
! Desired order of interpolation
  integer(kind=i4b), intent(in)  :: n0
! Result of interpolation
  complex(kind=c8b), intent(out) :: gf(nlmsb,nlmsb)
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-8
! -----------------------------------------
  integer(kind=i4b) :: ie, key(ne), info, i, j
  real(kind=r8b)    :: dist(ne), hscale
  complex(kind=c8b) :: dz, we(n0), den, gfi(nlmsb,nlmsb)

! Good morning
  if (n0 > ne) stop 'baryint_gf2: n0 > ne'
  if (n0 < 2)  stop 'baryint_gf2: n0 < 2'
  if (mod(n0,2) /= 0)  stop 'baryint_gf2: make n0 even'
!-----------------------------------------------------------------------
! First find the n0 points closest to z0
  do ie=1,ne
    dist(ie) = abs(z0 - ze(ie))
    key(ie)  = ie
  end do
! This sorts the energies by distance (ScaLAPACK routine)
  call dlasrt2_saved('I',ne,dist,key,info)
  if (info /= 0) stop 'baryint_gf2: fail in dlasrt2'
  write(*,'("z0=",2es16.8,/"ze=")') z0
  do i=1,n0
    write(*,'(3es16.8)') ze(key(i)), dist(i)
  end do
!-----------------------------------------------------------------------
! z0 corresponds to some energy already tabulated
  if (dist(1) < tol) then
    ie = key(1)
    call projected_gf(ie,ia,ja,gf,onsite,struct)
    return 
  end if
!-----------------------------------------------------------------------
  hscale = dist(n0)
! Compute the interpolation weights
  do j=1,n0
    we(j) = 1.d0
    do i=1,n0
      if (i == j) cycle
      dz = (ze(key(j)) - ze(key(i)))/hscale
      if (abs(dz) < tol) stop 'baryint_gf2: repeated energy points!'
      we(j) = dz*we(j)
    end do
!    we(j) = 1.d0/we(j)
    we(j) = (-1.d0)**j
  end do
  we(1)  = 0.5d0*we(1)
  we(n0) = 0.5d0*we(n0)
!-----------------------------------------------------------------------
! Compute the interpolated GF
  gf = 0.d0; den = 0.d0
  do i=1,n0
    ie = key(i)
    call projected_gf(ie,ia,ja,gfi,onsite,struct)
    dz = (z0 - ze(key(i)))/hscale
!   there is no divide by zero here; dist(1) would be zero above
    dz = we(i)/dz
    write(*,'("dz=",2es16.8)') dz
    gf = gf + dz*gfi
    den = den + dz
  end do
  gf = gf/den
! All done!
  end subroutine baryint_gf2
