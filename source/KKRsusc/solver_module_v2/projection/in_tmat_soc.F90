  subroutine in_tmat_soc(ie)

! Get t-matrix from outsusc_tmat.dat
! Reads the whole thing for all spin channels and a given energy point

  use global

  implicit none
  
  integer(kind=i4b), intent(in)  :: ie
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ja, ilms, jlms, i, j, ic, k, ib 
  character*60      :: header
! 1D Arrays for MPI gathering
  real(kind=c8b), allocatable    :: tmat_1d(:)
! Pointers  
  integer(kind=i4b), allocatable :: i4d_1dt(:,:,:), i1d_4dt(:,:)
  integer(kind=i4b)              :: ndimt, i4(3), irec
  integer(kind=i4b)              :: length_t
! Careful inc.p from JMcode WLENGTH
  integer(kind=i4b), parameter   :: wlength = 4
! ------------------------------------------------------------------

! Size tmat 
  ndimt = nlms**2*nasusc

! length strings
  length_t = wlength*2*ndimt

! Open tmat-file
  open(access='direct',recl=length_t,file='outsusc_tmat.dat',unit=iomain2,form='unformatted',status='old') 

! Tmat
  allocate(i4d_1dt(nlms,nlms,nasusc),i1d_4dt(3,ndimt))
  allocate(tmat_1d(2*ndimt))

! Pointers for 4d <--> 1d 
  k = 0
  do ia = 1, nasusc
    do ilms = 1, nlms
      do jlms = 1, nlms
        k  = k + 1
        i4d_1dt(ilms,jlms,ia) = k 
        i1d_4dt(:,k) = (/ilms,jlms,ia/)
      end do ! j
    end do ! i
  end do ! ia 

! Read in as written from "gather_tgmat_mpi.f90"
  irec = ie
  read(iomain2,rec=irec)(tmat_1d(i), i=1,2*ndimt)

! Gstructural 1d
  ic = 1
  do ib = 1, ndimt
    i4=i1d_4dt(:,ib)    
    ilms=i4(1); jlms=i4(2); ia=i4(3)
    i = alms2i(ilms,ia)
    j = alms2i(jlms,ia)  
!    tmatrix(i,j,ia,ie) = cmplx(tmat_1d(ic),tmat_1d(ic+1))
    ic = ic + 2 
  end do

  deallocate(tmat_1d,i4d_1dt,i1d_4dt)
! All done
  end subroutine in_tmat_soc
