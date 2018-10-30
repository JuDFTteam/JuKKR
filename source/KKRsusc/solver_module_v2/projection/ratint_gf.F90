  subroutine ratint_gf(ia,ja,onsite,struct,e0,gf)
! Rational function interpolation

  implicit none

! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! Desired complex energy
  complex(kind=c8b), intent(in)  :: e0
! Result of fit
  complex(kind=c8b), intent(out) :: gf(nlmsb,nlmsb)
! -----------------------------------------
! Number of energy points to include
  integer(kind=i4b), parameter :: npts = 10
  real(kind=r8b),    parameter :: tol  = 1.d-8
  complex(kind=c8b) :: ze(npts), gfdata(nlmsb,nlmsb,npts), gfint(npts), dgf
  integer(kind=i4b) :: i, j, ie0(npts)


! Are we in the upper complex plane?
  if (aimag(e0) < 0.d0) stop 'ratint_gf: Im E < 0'
  gf = 0.d0
! Find the closest energy points
  call sort_energies(npts,nesusc,e0,esusc,ie0)
! Desired energy is already stored
  if (abs(esusc(ie0(1)) - e0) < tol) then
    call projected_gf(ie0(i),ia,ja,gf,onsite,struct)
    gf = gfdata(ie0(1),:,:)
    return
  end if
! Fill up the GFs
  ze = esusc(ie0)
  do i=1,npts
    call projected_gf(ie0(i),ia,ja,gfdata(:,:,i),onsite,struct)
  end do
! Get coefficients
  do j=1,nlmsba(ja)
    do i=1,nlmsba(ia)
      gfint = gfdata(i,j,:)
      if (sum(abs(gfint)) > gfilter) then
        call zratint(npts,ze,gfint,e0,gf(i,j),dgf)
!        write(iodb,'("ratint_gf: ia,ja, i,j=",8i4," maxij=",2es16.8)') ia, ja, i2lmsb(:,i,ia), i2lmsb(:,j,ja), dgf
      end if
    end do
  end do
! All done!
  end subroutine ratint_gf

