  subroutine get_gfsum(ia,gfsum,egfsum,tgfsum,onsite,struct)
  use global

  implicit none

! Atom of interest
  integer(kind=i4b), intent(in)  :: ia
! Energy integrated GF
  complex(kind=c8b), intent(out) :: gfsum(nlmsb,nlmsb), egfsum(nlmsb,nlmsb), tgfsum(nlmsb,nlmsb)
! Onsite GF
  logical,           intent(in)  :: onsite
! Structural GF
  logical,           intent(in)  :: struct
! ----------------------------------------------------------------------
! i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi), cone = (1.d0,0.d0), czero = (0.d0,0.d0)
  integer(kind=i4b) :: ie, i, j, k
  real(kind=r8b)    :: zre, zim
  integer(kind=i4b), allocatable :: imap(:,:)
  complex(kind=c8b), allocatable :: gf(:,:)


  if (.not.lfit .and. nesusc < nescf) stop 'get_gfsum: nesusc /= nescf'
  allocate(imap(nlmsb,nlmsb),gf(nlmsb,nlmsb))
! use the energy mesh from SCF calculation here
! construct the energy integrated site-diagonal GF
  gfsum = 0.d0; egfsum = 0.d0; tgfsum = 0.d0
  do ie=1,nescf
!   ****  site diagonal part for local quantities  ****
!   GF at projection energies
    if (.not.lfit) then
      call projected_gf(ie,ia,ia,gf,onsite,struct)
!   GF from interpolation
    else if (ifit == 1) then
!      call baryint_gf2(nesusc,esusc,ia,ia,esusc(ie),numd,gf,onsite,struct)
!      call ratint_gf2(nesusc,esusc,ia,ia,escf(ie),numd,gf,onsite,struct)
      call biratint_gf(escf(ie),ia,ia,gf,onsite,struct)
!   GF from fit
    else if (ifit == 2) then
     call ratval_gf(ia,ia,escf(ie),gf)
    end if
!    do j=1,nlmsb
!      do i=1,nlmsb
!        write(*,'(3i4,2es16.8)') ie, i, j, gfsum(i,j)
!      end do
!    end do
!   (G^dagger - G)/(i 2pi)
    gfsum = gfsum + (conjg(transpose(descf(ie)*gf)) - descf(ie)*gf)/i2pi
!   (z* G^dagger - z G)/(i 2pi)
    if (lebandnoef) then
      egfsum = egfsum + (conjg(transpose(descf(ie)*escf(ie)*gf)) - descf(ie)*escf(ie)*gf)/i2pi
    else
      egfsum = egfsum + (conjg(transpose(descf(ie)*(escf(ie)-efscf)*gf)) - descf(ie)*(escf(ie)-efscf)*gf)/i2pi
    end if
!    egfsum = egfsum + (conjg(transpose(descf(ie)*escf(ie)*gf)) - descf(ie)*escf(ie)*gf)/i2pi
!    do k=1,3
!      tgfsum(:,:,k) = tgfsum(:,:,k) + matmul((conjg(transpose(descf(ie)*gf)) - descf(ie)*gf)/i2pi,torque(:,:,k,ia,ie))
!      tgfsum(:,:,k) = tgfsum(:,:,k) + (conjg(transpose(descf(ie)*matmul(gf,torque(:,:,k,ia,ie)))) - descf(ie)*matmul(gf,torque(:,:,k,ia,ie)))/i2pi
!    end do
!   ***************************************************
  end do
!  write(*,'("G(EF) for ia=",i4)') ia
!  do i=1,nlmsba(ia)
!    write(*,'(1000es10.1)') gf(i,1:nlmsba(ia))
!  end do
! matrix elements for torque
!  if (lsoc .and. isoc(ia) == 1) then
!    do k=1,3
!      call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ia),cone,gfsum,nlmsb,torque(:,:,k,ia,nescf),nlmsb,czero,tgfsum(:,:,k),nlmsb)
!    end do
!  end if
! filter small elements of the GF away
!  imap = 0
  do j=1,nlmsba(ia)
    do i=1,nlmsba(ia)
      zre = real(gfsum(i,j)); zim = aimag(gfsum(i,j))
      if (abs(zre) < gstol) zre = 0.d0
      if (abs(zim) < gstol) zim = 0.d0
      gfsum(i,j) = cmplx(zre,zim)
      if (abs(gfsum(i,j)) > gstol) imap(i,j) = 1
!      write(*,'(2i4,2es16.8)') i, j, gfsum(i,j)
    end do
  end do
!  do i=1,nlmsb
!    write(*,'(1000i2)') imap(i,1:nlmsb)
!  end do
!  do i=1,nlmsb
!    write(*,'("gfsum ",1000es10.1)') gfsum(i,1:nlmsb)
!  end do
! All done!
  end subroutine get_gfsum
