  subroutine in_gmat_soc(ie)

! Get structural GF from outsusc_green.dat
! Reads the whole thing for all spin channels and a given energy point

  use global

  implicit none
 
  integer(kind=i4b), intent(in) :: ie
! ------------------------------------------------------------------
  integer(kind=i4b) :: ia, ja, ilms, jlms, i, j, ic, ib, k
  character*60      :: header
! 1D Arrays for MPI gathering
  real(kind=r8b), allocatable    :: gmat_1d(:)
! Pointers  
  integer(kind=i4b), allocatable :: i4d_1dgs(:,:,:,:), i1d_4dgs(:,:)
  integer(kind=i4b)              :: ndimgs, i4(4), irec
  integer(kind=i4b)              :: length_gs
! Careful inc.p from JMcode WLENGTH
  integer(kind=i4b), parameter   :: wlength = 4
! ------------------------------------------------------------------

! Size gstruct
  ndimgs = (nlms*nasusc)**2

! length strings
  length_gs = wlength*2*ndimgs

! Open gmat-files
  open(access='direct',recl=length_gs,file='outsusc_green.dat',unit=iomain1,form='unformatted',status='old')

! Gstructural 
  allocate(i4d_1dgs(nlms,nlms,nasusc,nasusc),i1d_4dgs(4,ndimgs))
  allocate(gmat_1d(2*ndimgs))

! 1d <-> 4d pointer
  k = 0
  do ia = 1, nasusc
    do ja = 1, nasusc
      do ilms = 1, nlms
        do jlms = 1, nlms
          k  = k + 1
          i4d_1dgs(ilms,jlms,ia,ja) = k 
          i1d_4dgs(:,k) = (/ilms,jlms,ia,ja/)  
        end do ! j
      end do ! i
    end do ! ja
  end do ! ia 

! Read in as written from "gather_tgmat_mpi.f90"
  irec = ie 
  read(iomain1,rec=irec)(gmat_1d(i),i=1,2*ndimgs)

! Gstructural 1d
  ic = 1
  do ib = 1, ndimgs
    i4=i1d_4dgs(:,ib)    
    ilms=i4(1); jlms=i4(2); ia=i4(3); ja=i4(4)
    i = alms2i(ilms,ia)
    j = alms2i(jlms,ja)  
    gstruct(i,j,ie) = cmplx(gmat_1d(ic),gmat_1d(ic+1))
    ic = ic + 2 
  end do

  deallocate(gmat_1d,i4d_1dgs,i1d_4dgs)
! All done
  end subroutine in_gmat_soc
