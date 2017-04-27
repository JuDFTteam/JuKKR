  subroutine spinrot_scatt_sol(ie,magdir0,magdir1,rotate_onsite,rotate_regular)
! Rotates the onsite GF and regular solutions according to local spin rotations
  use global

  implicit none

! Energy point 
  integer(kind=i4b), intent(in) :: ie
! Initial spin axes
  real(kind=r8b), intent(in) :: magdir0(3,nasusc)
! Final spin axis
  real(kind=r8b), intent(in) :: magdir1(3,nasusc)
! Rotate the structural GF ?
  logical,        intent(in) :: rotate_onsite, rotate_regular
! -------------------------------------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
  real(kind=r8b),    parameter :: pi = 4.d0*atan(1.d0)
  complex(kind=c8b) :: spinrot(nsmax,nsmax,nasusc), trace, qe, ms(3), mo(3)
  complex(kind=c8b), allocatable :: pzsave(:,:), gfsave(:,:)
  complex(kind=c8b), allocatable :: sr1(:,:,:), sr2(:,:,:)
  integer(kind=i4b) :: i3(3), i2(2), ilms, jlms, ia, ja, jalms, ialms
  integer(kind=i4b) :: i, j, k, l, is, js, ib, jb, ilm, jlm, il, jl
  real(kind=r8b)    :: start, finish, time1, time2
! ------------------------------------------------------------------------------------------------- 

! Allocate only once Energy loop
  allocate(pzsave(nlmsb,nlms),gfsave(nlmsb,nlmsb),sr1(nlms,nlms,nasusc),sr2(nlmsb,nlmsb,nasusc))
 
! Spin rotation matrices
  sr1 = czero; sr2 = czero
  do ia=1,nasusc
    call spin_rotation(magdir0(:,ia),magdir1(:,ia),pauli,spinrot(:,:,ia))
!   for the nlms basis
    do ilm=1,lmmax
      do js=1,nsmax
        j = lms2i(ilm,js)
        do is=1,nsmax
          i = lms2i(ilm,is)
          sr1(i,j,ia) = spinrot(is,js,ia)
        end do
      end do
    end do
!   for the nlmsb basis
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
!      jl = i2lm(2,j)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
!       il = i2lm(2,i)
        if (ilm == jlm .and. ib == jb) sr2(i,j,ia) = spinrot(is,js,ia)
!       if (il == jl .and. ib == jb) sr2(i,j,ia) = spinrot(is,js,ia)
      end do
    end do
  end do
! **************
  do ia=1,nasusc
! **************
    pzsave = czero; gfsave = czero
    if (rotate_onsite) then 
!     ------------------------------------------------------------------
!                      Rotation of the single-site GF
!     ------------------------------------------------------------------
!     gfpqsave = gfpq.U^\dagger
      call zgemm('N','C',nlmsba(ia),nlmsba(ia),nlmsba(ia),cone,gfpq(:,:,ia,ie),nlmsb,sr2(:,:,ia),nlmsb,czero,gfsave,nlmsb)
!     gfpq = U.gfpqsave
      call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ia),cone,sr2(:,:,ia),nlmsb,gfsave,nlmsb,czero,gfpq(:,:,ia,ie),nlmsb)
    end if ! onsite 
    if (rotate_regular) then
!     ------------------------------------------------------------------
!                   Rotation of the scattering solutions
!     ------------------------------------------------------------------
!     pzsave = pzr.U^\dagger
      call zgemm('N','C',nlmsba(ia),nlms,nlms,cone,pzr(:,:,ia,ie),nlmsb,sr1(:,:,ia),nlms,czero,pzsave,nlmsb)
!     pzr = U.pzsave
      call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cone,sr2(:,:,ia),nlmsb,pzsave,nlmsb,czero,pzr(:,:,ia,ie),nlmsb)
!     Careful new version, saving pzl not pzl^T 
!     pzsave = pzl.U^\dagger
!      call zgemm('N','C',nlmsba(ia),nlms,nlms,cone,pzl(:,:,ia,ie),nlmsb,sr1(:,:,ia),nlms,czero,pzsave,nlmsb)
!     pzl = U.pzsave
!      call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cone,sr2(:,:,ia),nlmsb,pzsave,nlmsb,czero,pzl(:,:,ia,ie),nlmsb)
!     **************************** THE PROBLEM WAS HERE NOW I DO NOT SAVE pzl^T but i save pzl from kkrflex*******************
!     pzsave = pzl.U^T
      call zgemm('N','T',nlmsba(ia),nlms,nlms,cone,pzl(:,:,ia,ie),nlmsb,sr1(:,:,ia),nlms,czero,pzsave,nlmsb)
!     pzl = U^*.pzsave
      gfsave = conjg(sr2(:,:,ia))
      call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cone,gfsave,nlmsb,pzsave,nlmsb,czero,pzl(:,:,ia,ie),nlmsb)
!     ************************************************************************************************************************
    end if 
! ******
  end do
! ******
! ----------------------------------------------------------------------
  deallocate(pzsave,gfsave,sr1,sr2)
! All done!
  end subroutine spinrot_scatt_sol
