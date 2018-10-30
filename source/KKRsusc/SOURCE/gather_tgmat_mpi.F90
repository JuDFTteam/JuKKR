  subroutine gather_tgmat_mpi(lmaxd,ielast,ez,wz,my_rank,mpi_size,mpi_iebounds)

  ! --> modules
  use type_config
  use mod_config
  use mod_physic_params, only: cvlight
  use global

  implicit none

  ! --> types 
  type(config_type)             ::  config
  ! --> energy point label, mpi_stuff
  integer(kind=i4b), intent(in) :: lmaxd, ielast, my_rank, mpi_size
  ! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(in) :: ez(ielast), wz(ielast)
  ! --> mpi_bounds
  integer(kind=i4b), intent(in) :: mpi_iebounds(2,0:mpi_size-1) 
  ! ***********************************************************************************
  integer(kind=i4b) :: ia, ja, ih, jh, il, il2, im, im2, ilm, ilm2, ilmsn, ilmsn2
  integer(kind=i4b) :: nb, ib, jb, i, j, k, l, ilms, jlms, irank
  ! ***********************************************************************************
  ! Write for all energies
  logical                        :: write_mpi
  ! Dimensions of 1d Arrays for MPI
  integer(kind=i4b)              :: ndimgs, ndimgs_all, ntot_lmsb, lmsize, ie, ierror
  integer(kind=i4b)              :: ndimt, ndimt_all, ne_rank, ndim_skip, ne_nusd 
  integer(kind=i4b)              :: irec
  complex(kind=c8b)              :: ek_susc 
  ! Pointers  
  integer(kind=i4b), allocatable :: i4d_1dgs(:,:,:,:), i1d_4dgs(:,:)
  integer(kind=i4b), allocatable :: i4d_1dt(:,:,:)   , i1d_4dt(:,:)
  ! 1D Arrays for MPI gathering
  complex(kind=c8b), allocatable :: gmat_1d(:,:), gmat_1d_all(:)
  complex(kind=c8b), allocatable :: tmat_1d(:,:), tmat_1d_all(:)
  ! Testing the arrays
  complex(kind=c8b), allocatable :: gmat_out(:,:), tmat_out(:,:)
  integer(kind=i4b)              :: is, js, jlm
  integer(kind=i4b)              :: length_gs, length_t
  ! Careful inc.p from JMcode WLENGTH
  integer(kind=i4b), parameter   :: wleng = 4 
  ! ***********************************************************************************

  ! Define MPI constants
#ifdef MPI
  INCLUDE "mpif.h"
#endif

  ! Size gstruct
  ndimgs = (nlms*nasusc)**2 
  ! Size tmat
  ndimt = nasusc*nlms**2

  ! length of the strings
  length_gs = wleng*2*ndimgs
  length_t  = wleng*2*ndimt

  ! Open io files 
  open(access='direct',recl=length_gs,file='outsusc_green.dat',unit=iomain1,form='unformatted',status='replace')
  open(access='direct',recl=length_gs,file='outsusc_green_test.dat',unit=1500,form='unformatted',status='replace')
  open(access='direct',recl=length_t,file='outsusc_tmat.dat',  unit=iomain2,form='unformatted',status='replace')

!  write(*,*) my_rank, mpi_iebounds(1,my_rank), mpi_iebounds(2,my_rank)
  ! Allocat e 1d array for mpi gathering
  allocate(gmat_1d(ndimgs,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))
  allocate(tmat_1d(ndimt, mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))

  ! Pointers for 4d <--> 1d 
  ! Gstructural 
  allocate(i4d_1dgs(nlms,nlms,nasusc,nasusc))
  allocate(i1d_4dgs(4,ndimgs))
  k = 0
  do ia = 1, nasusc
    do ja = 1, nasusc
      do i = 1, nlms
        do j = 1, nlms
          k  = k + 1
          i4d_1dgs(i,j,ia,ja) = k 
          i1d_4dgs(:,k) = (/i,j,ia,ja/)  
        end do ! j
      end do ! i
    end do ! ja
  end do ! ia 
  ! Pointers for 4d <--> 1d 
  ! Tmat 
  allocate(i4d_1dt(nlms,nlms,nasusc))
  allocate(i1d_4dt(3,ndimt))
  k = 0
  do ia = 1, nasusc
    do i = 1, nlms
      do j = 1, nlms
        k  = k + 1
        i4d_1dt(i,j,ia) = k 
        i1d_4dt(:,k) = (/i,j,ia/)  
      end do ! j
    end do ! i
  end do ! ia 

  ! Loop over energies
  do ie = mpi_iebounds(1,my_rank), mpi_iebounds(2,my_rank)
 
    !   Gstructural 1d 
    do ia = 1, nasusc
      do ja = 1, nasusc 
        do ilms = 1, nlms
          do jlms = 1, nlms
            i = alms2i(ilms,ia)
            j = alms2i(jlms,ja)
            k = i4d_1dgs(ilms,jlms,ia,ja)         
            gmat_1d(k,ie) = gstruct(i,j,ie)
          end do ! j 
        end do ! i
      end do ! ja
    end do ! ia
    !   Tmat 1d
    do ia = 1, nasusc
      ih = iasusc(ia)
      do i = 1, nlms
        do j = 1, nlms
          k = i4d_1dt(i,j,ia) 
          tmat_1d(k,ie) = tmatrix(i,j,ia,ie)
        end do ! j
      end do ! i
    end do ! ia

  end do ! ie
#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  if (ierror/=MPI_SUCCESS) stop 'gathering tgmat 1'
#endif

  ! Write into a single unformatted file   
  do ie = mpi_iebounds(1,my_rank), mpi_iebounds(2,my_rank)
    irec = ie   
    !   gmat  
    write(iomain1,rec=irec) (gmat_1d(i,ie), i=1,ndimgs)
    if (ie == 1) then
      write(1500,rec=irec) (gmat_1d(i,ie),i=1,ndimgs)
      do i = 1, ndimgs
        write(2000,'(1I4,2F12.8)') i, gmat_1d(i,ie)
      end do ! i 
    end if 
  end do ! ie
#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  if (ierror/=MPI_SUCCESS) stop 'gathering tgmat 2'
#endif
  do ie = mpi_iebounds(1,my_rank), mpi_iebounds(2,my_rank)
    irec = ie
    !   tmat   
    write(iomain2,rec=irec) (tmat_1d(i,ie), i=1,ndimt)
  end do 
#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  if (ierror/=MPI_SUCCESS) stop 'gathering tgmat 3'
#endif

  deallocate(gmat_1d,tmat_1d,i4d_1dgs,i1d_4dgs,i4d_1dt,i1d_4dt)
  ! All done
  end subroutine gather_tgmat_mpi
