  subroutine ratint_gf2(ne,ze,ia,ja,z0,n0,gf,onsite,struct)
! Rational interpolation with n0 points
! Works for any point distribution, as long as n0 is kept low
! Extension of standard NR algorithm
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
  real(kind=r8b), parameter :: tol = 1.d-8, small = 1.d-25
! -----------------------------------------
  integer(kind=i4b) :: ie, key(ne), info, i, j, iq, jq
  real(kind=r8b)    :: dist(ne)
  complex(kind=c8b) :: dz, fac, zn0(n0), gfn0(n0), dgf
  complex(kind=c8b), allocatable :: gf0(:,:,:)

! Good morning
  if (n0 > ne) stop 'ratint_gf2: n0 > ne'
  if (n0 < 0)  stop 'ratint_gf2: n0 < 0'
  gf = 0.d0
!-----------------------------------------------------------------------
! First find the n0 points closest to z0
  do ie=1,ne
    dist(ie) = abs(z0 - ze(ie))
    key(ie)  = ie
  end do
! This sorts the energies by distance (ScaLAPACK routine)
  call dlasrt2_saved('I',ne,dist,key,info)
  if (info /= 0) stop 'ratint_gf2: fail in dlasrt2'
!  write(*,'("z0=",2es16.8,/"ze=")') z0
!  do i=1,n0
!    write(*,'(3es16.8)') ze(key(i)), dist(i)
!  end do
!-----------------------------------------------------------------------
! z0 corresponds to some energy already tabulated
  if (dist(1) < tol) then
    ie = key(1)
    call projected_gf(ie,ia,ja,gf,onsite,struct)
    return
  end if
!-----------------------------------------------------------------------
! Initialize
  allocate(gf0(n0,nlmsba(ia),nlmsba(ja)))
  do i=1,n0
    ie = key(i)
    zn0(i) = ze(key(i))
    call projected_gf(ie,ia,ja,gf,onsite,struct)
    do jq=1,nlmsba(ja)
      do iq=1,nlmsba(ia)
        gf0(i,iq,jq) = gf(iq,jq)
      end do
    end do
  end do
! Rational interpolation from NR
  do jq=1,nlmsba(ja)
    do iq=1,nlmsba(ia)
      gfn0 = gf0(:,iq,jq)
      call zratint(n0,zn0,gfn0,z0,gf(iq,jq),dgf)
    end do
  end do
  deallocate(gf0)
! All done!
  end subroutine ratint_gf2
