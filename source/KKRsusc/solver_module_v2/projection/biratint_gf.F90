  subroutine biratint_gf(ze,ia,ja,gfint,onsite,struct)
! Block of the projected GF
! Assumes energies are distributed in panels parallel to real axis
  use global

  implicit none

! which energy
  complex(kind=c8b), intent(in)  :: ze
! which atoms
  integer(kind=i4b), intent(in)  :: ia, ja
! the block of the GF in the projection basis
  complex(kind=c8b), intent(out) :: gfint(nlmsb,nlmsb)
! whether to include the onsite and the structural parts
  logical,           intent(in)  :: onsite, struct
! ----------------------------------------------------------------------
  integer(kind=i4b), parameter :: npan = 5, ninter = 5
  real(kind=r8b),    parameter :: tol = 1.d-12
  integer(kind=i4b) :: npts(npan)
  complex(kind=c8b) :: gfere(nlmsb,nlmsb,ninter), gfeim(nlmsb,nlmsb,npan)
  real(kind=r8b)    :: ere, eim
  integer(kind=i4b) :: ipan, ipan0, ipan1, inear, i, j, ie
  complex(kind=c8b) :: zepan(npan), gf(ninter), zeinter(ninter), gf2(npan), dgf

  ere = real(ze)
  npts(:) = (/11,21,41,81,161/)
! as many energies as total number of points in panels
  if (sum(npts(:)) /= nesusc) stop 'biratint_gf: sum(npts) /= nesusc'
! find interpolated GF in each panel
  do ipan=1,npan
!   where does each panel start and end
    ipan0 = sum(npts(1:ipan-1))
    ipan1 = sum(npts(1:ipan))
!   energy in current panel closest to desired energy
    inear = minloc(abs(esusc(ipan0+1:ipan1)-ze),dim=1)
!    eim = aimag(esusc(ipan0+inear))
    eim = aimag(esusc(ipan0+1))
!   interpolation energy in current panel
    zepan(ipan) = cmplx(ere,eim)
!   energies to use for interpolation:
!   - not enough points
    if (npts(ipan) < ninter) stop 'biratint_gf: npts(ipan) < ninter'
!   - as many as needed
    if (npts(ipan) == ninter) then
      inear = ipan0
!   - more than needed: try to center
    else 
      inear = ipan0 + inear - ninter/2
      if (inear < ipan0) inear = ipan0
      if (inear + ninter > ipan1) inear = ipan1 - ninter
    end if
    write(iodb,'(" ipan, ipan0, ipan1, inear=",4i4)') ipan, ipan0, ipan1, inear
!   projected GF for interpolation energies
    do ie=1,ninter
      zeinter(ie) = esusc(ie+inear)
      gfere(:,:,ie) = 0.d0
      call projected_gf(ie+inear,ia,ja,gfere(:,:,ie),onsite,struct)
      write(iodb,'(2i4,6es16.8)') ipan, ie+inear, ze, zeinter(ie), sum(gfere(:,:,ie))
    end do
!   rational interpolation on current panel
    gfeim(:,:,ipan) = 0.d0
    do j=1,nlmsba(ja)
      do i=1,nlmsba(ia)
        gf(:) = gfere(i,j,:)
        if (sum(abs(gf)) > tol) then
          call zratint(ninter,zeinter,gf,zepan(ipan),gfeim(i,j,ipan),dgf)
!        write(iodb,'("biratint_gf: ia,ja, i,j=",8i4," maxij=",2es16.8)') ia, ja, i2lmsb(:,i,ia), i2lmsb(:,j,ja), dgf
        end if
      end do
    end do
    write(iodb,'(i4,4es16.8)') ipan, zepan(ipan), sum(gfeim(:,:,ipan))
    write(iodb,'("------------------------------------------------------")')
  end do
! rational interpolation across the panels
  gfint(:,:) = 0.d0
  do j=1,nlmsba(ja)
    do i=1,nlmsba(ia)
      gf2(:) = gfeim(i,j,:)
      if (sum(abs(gf2)) > tol) then
        call zratint(npan,zepan,gf2,ze,gfint(i,j),dgf)
!      write(iodb,'("biratint_gf: ia,ja, i,j=",8i4," maxij=",2es16.8)') ia, ja, i2lmsb(:,i,ia), i2lmsb(:,j,ja), dgf
      end if
    end do
  end do
  write(iodb,'("  ia=",i4,4es16.8)') ia, ze, sum(gfint(:,:))
! All done!
  end subroutine biratint_gf
