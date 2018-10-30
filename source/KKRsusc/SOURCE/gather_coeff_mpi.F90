  subroutine gather_coeff_mpi(lmaxd,ielast,ez,wz,config,my_rank,mpi_size,mpi_iebounds)

! --> modules
  use type_config
  use mod_config
  use global
  use mod_physic_params, only: cvlight

  implicit none

! --> types 
  type(config_type)             ::  config
! --> energy point label, mpi_stuff
  integer(kind=i4b), intent(in) :: lmaxd, ielast, my_rank, mpi_size
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(in) :: ez(ielast), wz(ielast)
! --> mpi_bounds
  integer(kind=i4b), intent(in) :: mpi_iebounds(2,0:mpi_size-1) 
! *********************************************************************************
  integer(kind=i4b) :: ia, ih, il, il2, im, im2, ilm, ilm2, ilmsn, ilmsn2
  integer(kind=i4b) :: nb, ib, jb, i, j, k, l, irank  ! basis
! *********************************************************************************
! Dimensions of 1d Arrays for MPI
  integer(kind=i4b)              :: ndim, ndim_all, ntot_lmsb, lmsize, ie, ierror
  integer(kind=i4b)              :: ndimgo, ndimgo_all, ndim_skip, ne_rank, ne_nusd
  complex(kind=c8b)              :: ek_susc 
! Pointers  
  integer(kind=i4b), allocatable :: i4d_1d(:,:,:,:), i1d_4d(:,:)
  integer(kind=i4b), allocatable :: i4d_1dgo(:,:,:,:), i1d_4dgo(:,:)
! 1D Arrays for MPI gathering 
! Regular solutions
  complex(kind=c8b), allocatable :: pzr_1d(:), pzr_1d_all(:)  
  complex(kind=c8b), allocatable :: pzl_1d(:), pzl_1d_all(:) 
  complex(kind=c8b), allocatable :: fzr_1d(:), fzr_1d_all(:)  
  complex(kind=c8b), allocatable :: fzl_1d(:), fzl_1d_all(:)
! Onsite green function
  complex(kind=c8b), allocatable :: gfpq_1d(:), gfpq_1d_all(:)  
! *********************************************************************************

! Define MPI constants
#ifdef MPI
  INCLUDE "mpif.h"
#endif

! Find the right size for the 1D arrray 
  ndim = 0
  do ia = 1, nasusc
    ndim = ndim + nlmsba(ia)*nlms
  end do ! nasusc

! Size gonsite
  ndimgo = 0
  do ia = 1, nasusc
    ndimgo = ndimgo + nlmsba(ia)**2
  end do ! 

! Include all energies
  ndim_all = ndim*nesusc 
  ndimgo_all = ndimgo*nesusc

! Divide on cores
  ne_rank = ceiling(nesusc*1.d0/mpi_size)
  ndim    = ne_rank*ndim
  ndimgo  = ne_rank*ndimgo

! Allocate 1d array for mpi gathering
  allocate(pzr_1d(ndim))   
  allocate(pzl_1d(ndim))
  allocate(gfpq_1d(ndimgo))   
  if (my_rank == 0) then
    allocate(pzr_1d_all(ndim*mpi_size))
    allocate(pzl_1d_all(ndim*mpi_size))
    allocate(gfpq_1d_all(ndimgo*mpi_size))
  else
    allocate(pzr_1d_all(1))
    allocate(pzl_1d_all(1))
    allocate(gfpq_1d_all(1))
  end if ! my_rank 
  if (isra == 1) then
    allocate(fzr_1d(ndim))
    allocate(fzl_1d(ndim))
    if (my_rank == 0) then
      allocate(fzr_1d_all(ndim*mpi_size))
      allocate(fzl_1d_all(ndim*mpi_size))
    else
      allocate(fzr_1d_all(1))
      allocate(fzl_1d_all(1))
    end if ! my_rank
  end if 

! Pointer for 4d <--> 1d
! Max nlmsb 
  ntot_lmsb = nlms*nbmax 
  allocate(i4d_1d(ntot_lmsb,nlms,nasusc,nesusc))
  allocate(i1d_4d(4,ndim_all))
  k = 0
  do ie = 1, nesusc
    do ia = 1, nasusc
      do i = 1, nlmsba(ia)
        do j = 1, nlms
          k  = k + 1
          i4d_1d(i,j,ia,ie) = k 
          i1d_4d(:,k) = (/i,j,ia,ie/)  
        end do ! j
      end do ! i
    end do ! ia 
  end do ! ie

! Pointers for 4d <--> 1d 
! Gonsite 
  ntot_lmsb = nlms*nbmax 
  allocate(i4d_1dgo(ntot_lmsb,ntot_lmsb,nasusc,nesusc))
  allocate(i1d_4dgo(4,ndimgo_all))
  k = 0
  do ie = 1, nesusc
    do ia = 1, nasusc
      do i = 1, nlmsba(ia)
        do j = 1, nlmsba(ia) 
          k  = k + 1
          i4d_1dgo(i,j,ia,ie) = k 
          i1d_4dgo(:,k) = (/i,j,ia,ie/)  
        end do ! j
      end do ! i
    end do ! ia 
  end do ! ie

! Loop over energies
  l = 0
  do ie = mpi_iebounds(1,my_rank), mpi_iebounds(2,my_rank)
    l = l + 1 
    do ia = 1, nasusc
      ih = iasusc(ia) 
      do i = 1, nlmsba(ia)
        do j = 1, nlms
          k = i4d_1d(i,j,ia,l) 
          pzr_1d(k) = pzr(i,j,ia,ie) 
          pzl_1d(k) = pzl(i,j,ia,ie)
          if (isra == 1) then
            fzr_1d(k) = fzr(i,j,ia,ie)
            fzl_1d(k) = fzl(i,j,ia,ie)
          end if 
        end do ! j
      end do ! i
    end do ! ia

!   Gonsite  
    do ia = 1, nasusc
      ih = iasusc(ia)
      do i = 1, nlmsba(ia)
        do j = 1, nlmsba(ia)
          k = i4d_1dgo(i,j,ia,l) 
          gfpq_1d(k) = gfpq(i,j,ia,ie)
        end do ! j
      end do ! i
    end do ! ia
  end do ! ie

! ******************************************
! Mpi gather for all energies on my_rank = 0
! ******************************************
#ifdef MPI
  call MPI_Gather(pzr_1d,ndim,MPI_DOUBLE_COMPLEX,pzr_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
  call MPI_Gather(pzl_1d,ndim,MPI_DOUBLE_COMPLEX,pzl_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
  call MPI_Gather(gfpq_1d,ndimgo,MPI_DOUBLE_COMPLEX,gfpq_1d_all,ndimgo,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
  if (isra == 1) then
    call MPI_Gather(fzr_1d,ndim,MPI_DOUBLE_COMPLEX,fzr_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
    call MPI_Gather(fzl_1d,ndim,MPI_DOUBLE_COMPLEX,fzl_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror) 
  end if
#endif

! Write out --> main thread
  if (my_rank == 0) then

    deallocate(pzr,pzl,gfpq)
    if (isra == 1) then
      deallocate(fzl,fzr,gfps,gffq,gffs)
    end if ! isra
    allocate(pzl(nlmsb,nlms,nasusc,nesusc))       
    allocate(pzr(nlmsb,nlms,nasusc,nesusc))       
    allocate(gfpq(nlmsb,nlmsb,nasusc,nesusc))    
    if (isra == 1) then
      allocate(fzl(nlmsb,nlms,nasusc,nesusc))     
      allocate(fzr(nlmsb,nlms,nasusc,nesusc))    
      allocate(gfps(nlmsb,nlmsb,nasusc,nesusc))   
      allocate(gffq(nlmsb,nlmsb,nasusc,nesusc))   
      allocate(gffs(nlmsb,nlmsb,nasusc,nesusc))   
     end if ! isra 

!   Initialize   
    ndim_skip = 0
!   Regular coefficients
    do irank = 0, mpi_size-1
!     How many non used spots 
      ne_nusd = ne_rank-(mpi_iebounds(2,irank)-mpi_iebounds(1,irank)+1)       
      do ie = mpi_iebounds(1,irank), mpi_iebounds(2,irank)
        do ia = 1, nasusc
          do i = 1, nlmsba(ia)
            do j = 1, nlms 
              k = i4d_1d(i,j,ia,ie) + ndim_skip
              ! Put back to the coefficients         
              pzr(i,j,ia,ie) = pzr_1d_all(k)
              pzl(i,j,ia,ie) = pzl_1d_all(k)
              if (isra == 1) then
                fzr(i,j,ia,ie) = fzr_1d_all(k)
                fzl(i,j,ia,ie) = fzl_1d_all(k)
              end if
            end do ! j
          end do ! i
        end do ! ia
      end do ! ie
!     Skip the garbage
      ndim_skip = ndim_skip + (ndim/ne_rank)*ne_nusd
    end do ! irank 

!   Initialize   
    ndim_skip = 0
!   Gonsite
    do irank = 0, mpi_size-1
!     How many non used spots 
      ne_nusd = ne_rank-(mpi_iebounds(2,irank)-mpi_iebounds(1,irank)+1)       
      do ie = mpi_iebounds(1,irank), mpi_iebounds(2,irank)
        do ia = 1, nasusc
          do i = 1, nlmsba(ia)
            do j = 1, nlmsba(ia)
              k = i4d_1dgo(i,j,ia,ie) + ndim_skip
              ! Put back onsite gf         
              gfpq(i,j,ia,ie) = gfpq_1d_all(k) 
              if (isra == 1) then
                !!!!!!!!!!!! add the small components !!!!!!!!!!!!!!! 
              end if
            end do ! j
          end do ! i
        end do ! ia
      end do ! ie
!     Skip the garbage
      ndim_skip = ndim_skip + (ndimgo/ne_rank)*ne_nusd
    end do ! irank 

    ! lms blocks
    lmsize = 2*(lmaxd+1)**2 
    do ie = 1, ielast 
      if (config%nsra == 1) ek_susc = sqrt(ez(ie))                                          
      if (config%nsra == 2) ek_susc = sqrt(ez(ie) + (ez(ie)/cvlight)*(ez(ie)/cvlight))      
      call out_coeffs_soc(lmsize,ie,ez(ie),ek_susc,wz(ie))
    end do

  end if ! my_rank 

! Mpi Barrier ()
#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD,ierror)
#endif 
! if(ierror/=MPI_SUCCESS) stop 'Gathering_coeff'
  deallocate(pzr_1d,pzl_1d,gfpq_1d,pzr_1d_all,pzl_1d_all,gfpq_1d_all)

! All done
  end subroutine gather_coeff_mpi
