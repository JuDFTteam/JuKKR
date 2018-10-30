  subroutine lippmann_schwinger(ie,natoms,nlms,nlmsb,nlmsba,deltapot,pzl0,pzr0,dtmat,gf0)
  use global, only: i4b, r8b, c8b , iodb!, i2lmsb, i2lms

  implicit none

! Energy point
  integer(kind=i4b), intent(in)    :: ie  
! Number of atoms
  integer(kind=i4b), intent(in)    :: natoms
! Maximum combined dimension of angular momentum and spin
  integer(kind=i4b), intent(in)    :: nlms
! Maximum combined dimension of angular momentum, spin and basis functions
  integer(kind=i4b), intent(in)    :: nlmsb
! Pointers to storage
  integer(kind=i4b), intent(in)    :: nlmsba(natoms)
! Change in potential
  complex(kind=c8b), intent(in)    :: deltapot(nlmsb,nlmsb,natoms)
! Initial lhs scattering solution (overwritten)
  complex(kind=c8b), intent(inout) :: pzl0(nlmsb,nlms,natoms)
! Initial rhs scattering solution (overwritten)
  complex(kind=c8b), intent(inout) :: pzr0(nlmsb,nlms,natoms)
! Change in t-matrix from change in potential (overwritten)
  complex(kind=c8b), intent(inout) :: dtmat(nlms,nlms,natoms)
! Old single-site GF (overwritten)
  complex(kind=c8b), intent(inout) :: gf0(nlmsb,nlmsb,natoms)
! ----------------------------------------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-4
! complex parameters
  complex(kind=c8b), parameter :: cone   = ( 1.d0, 0.d0)
  complex(kind=c8b), parameter :: cminus = (-1.d0, 0.d0)
  complex(kind=c8b), parameter :: czero  = ( 0.d0, 0.d0)
! permutation array
  integer(kind=i4b), allocatable :: ipiv(:)
! new scattering solutions
  complex(kind=c8b), allocatable :: pzl(:,:), pzr(:,:)
! new single-site GF
  complex(kind=c8b), allocatable :: gf(:,:)
! change in t-matrices
  complex(kind=c8b), allocatable :: tmatl(:,:), tmatr(:,:)
! work array
  complex(kind=c8b), allocatable :: work(:,:)
! misc
  integer(kind=i4b) :: i, j, ia, info
  real(kind=r8b)    :: maxelem


  allocate(pzl(nlmsb,nlms),pzr(nlmsb,nlms),gf(nlmsb,nlmsb),tmatl(nlms,nlms),tmatr(nlms,nlms),work(nlmsb,nlmsb),ipiv(nlmsb))
  do ia=1,natoms
!   --> construct 1 - G0.deltaV
    call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ia),cminus,gf0(:,:,ia),nlmsb,deltapot(:,:,ia),nlmsb,czero,work,nlmsb)
    do i=1,nlmsba(ia)
      work(i,i) = work(i,i) + cone
    end do
!   --> LU decomposition of 1 - G0.deltaV
    call zgetrf(nlmsba(ia),nlmsba(ia),work,nlmsb,ipiv,info)
    if (info /= 0) then
      write(*,*) "LU fail 1 - G0.deltaV at ie,ia=", ie, ia, info
      stop
    end if
!   --> new single-site GF: (1 - G0.deltaV).G = G0
    gf = gf0(:,:,ia)
    call zgetrs('N',nlmsba(ia),nlmsba(ia),work,nlmsb,ipiv,gf,nlmsb,info)
    if (info /= 0) then
      write(*,*) "single-site GF failed at ie,ia=", ie, ia
      stop
    end if
!   --> solve Lipp-Schw for rhs wfn: (1 - G0.deltaV).R^r = R0^r
!     nlms are the boundary conditions
    pzr = pzr0(:,:,ia)  ! the initial rhs solution is needed in t^l
    call zgetrs('N',nlmsba(ia),nlms,work,nlmsb,ipiv,pzr,nlmsb,info)
    if (info /= 0) then
      write(*,*) "rhs Lipp Schw failed at ie,ia=", ie, ia
      stop
    end if
!    write(iodb,*) "non-zero elements of rhs wfn"
!    write(iodb,'("  ib ilm  is jlm  js      pzr0      pzr")')
!    maxelem = maxval(abs(pzr))
!    write(iodb,'("maxelem  ",es14.6)') maxelem
!    do j=1,nlms
!      do i=1,nlmsba(ia)
!        if (abs(pzr(i,j)) > tol*maxelem) then
!          write(iodb,'(5i4,4es14.6)') i2lmsb(:,i,ia), i2lms(:,j), pzr0(i,j,ia), pzr(i,j)
!        else
!          pzr(i,j) = 0.d0
!        end if
!      end do
!    end do
!   --> change in t-matrix: t^r = R0^l.Vsoc.R^r
!   deltaV.R^r --> work
    call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cone,deltapot(:,:,ia),nlmsb,pzr,nlmsb,czero,work,nlmsb)
!   R0^l.work -> t^r
    call zgemm('T','N',nlms,nlms,nlmsba(ia),cone,pzl0(:,:,ia),nlmsb,work,nlmsb,czero,tmatr,nlms)
!    write(iodb,*) "tmatr for ia=", ia
!    do j=1,nlms
!      write(iodb,'(100es10.2)') tmatr(:,j)
!    end do
!   --> construct 1 - G0^T.deltaV^T
    call zgemm('T','T',nlmsba(ia),nlmsba(ia),nlmsba(ia),cminus,gf0(:,:,ia),nlmsb,deltapot(:,:,ia),nlmsb,czero,work,nlmsb)
    do i=1,nlmsba(ia)
      work(i,i) = work(i,i) + cone
    end do
!   --> LU decomposition of 1 - G0^T.deltaV^T
    call zgetrf(nlmsba(ia),nlmsba(ia),work,nlmsb,ipiv,info)
    if (info /= 0) then
      write(*,*) "LU fail 1 - G0^T.deltaV^T at ie,ia=", ie, ia, info
      stop
    end if
!   --> solve Lipp-Schw for left hand side wfn: (1 - G0^T.Vsoc^T).R^l = R0^l
!     nlms are the boundary conditions
    pzl = pzl0(:,:,ia)
    call zgetrs('N',nlmsba(ia),nlms,work,nlmsb,ipiv,pzl,nlmsb,info)
    if (info /= 0) then
      write(*,*) "lhs Lipp Schw failed at ie,ia=", ie, ia
      stop
    end if
!    write(iodb,*) "non-zero elements of lhs wfn"
!    write(iodb,'("  ib ilm  is jlm  js      pzl0      pzl")')
!    maxelem = maxval(abs(pzl))
!    write(iodb,'("maxelem  ",es14.6)') maxelem
!    do j=1,nlms
!      do i=1,nlmsba(ia)
!        if (abs(pzl(i,j)) > tol*maxelem) then
!          write(iodb,'(5i4,4es14.6)') i2lmsb(:,i,ia), i2lms(:,j), pzl0(i,j,ia), pzl(i,j)
!        else
!          pzl(i,j) = 0.d0
!        end if
!      end do
!    end do
!   --> change in t-matrix: t^l = R^l.Vsoc.R0^r
!   R^l.Vsoc -> work
    call zgemm('T','N',nlms,nlmsba(ia),nlmsba(ia),cone,pzl,nlmsb,deltapot(:,:,ia),nlmsb,czero,work,nlmsb)
!   work.R0^r -> t^l
    call zgemm('N','N',nlms,nlms,nlmsba(ia),cone,work,nlmsb,pzr0(:,:,ia),nlmsb,czero,tmatl,nlms)
!    write(iodb,*) "tmatl for ia=", ia
!    do j=1,nlms
!      write(iodb,'(100es10.2)') tmatl(:,j)
!    end do
!    maxelem = maxval(abs(tmatl - tmatr))
!    write(iodb,'("|tmatl - tmatr| for ia=",i4," is ",es14.6)') ia, maxelem
!    do j=1,nlms
!      write(iodb,'(100es10.2)') abs(tmatl(:,j) - tmatr(:,j))
!    end do
!   --> update wfns
    pzl0(:,:,ia) = pzl; pzr0(:,:,ia) = pzr
!   --> update t-matrix
    dtmat(:,:,ia) = tmatr
!   --> update single-site GF
    gf0(:,:,ia) = gf
  end do
  deallocate(pzl,pzr,gf,tmatl,tmatr,work,ipiv)
! All done!
  end subroutine lippmann_schwinger
