  subroutine read_proj()

  use global

  implicit none

! ---------------------------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, il, il2, im, im2, ilm, ilm2, ilmsn, ilmsn2
  integer(kind=i4b) :: nb, ib, jb, i, j, k, l  ! basis
! ---------------------------------------------------------------------------------
! Dimensions of 1d Arrays for MPI
  integer(kind=i4b)              :: ndim, ndim_all, ntot_lmsb, lmsize, ie, ierror
  integer(kind=i4b)              :: ndimgo, ndimgo_all
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
! ---------------------------------------------------------------------------------


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
  ndim   = ceiling(nesusc*1.d0/mpi_size)*ndim
  ndimgo = ceiling(nesusc*1.d0/mpi_size)*ndimgo

! Allocate 1d array for mpi gathering
  allocate(pzr_1d(ndim))   
  allocate(pzr_1d_all(ndim_all))
  allocate(pzl_1d(ndim))
  allocate(pzl_1d_all(ndim_all))
  allocate(gfpq_1d(ndimgo))   
  allocate(gfpq_1d_all(ndimgo_all))
  if (isra == 1) then
    allocate(fzr_1d(ndim))
    allocate(fzr_1d_all(ndim_all))
    allocate(fzl_1d(ndim))
    allocate(fzl_1d_all(ndim_all))
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
          pzr_1d(k) = pzr(i,j,ih,ie) 
          pzl_1d(k) = pzl(i,j,ih,ie)
          if (isra == 1) then
            fzr_1d(k) = fzr(i,j,ih,ie)
            fzl_1d(k) = fzl(i,j,ih,ie)
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
          gfpq_1d(k) = gfpq(i,j,ih,ie)
        end do ! j
      end do ! i
    end do ! ia
  end do ! ie

! ******************************************
! Mpi gather for all energies on my_rank = 0
! ******************************************
  call MPI_Gather(pzr_1d,ndim,MPI_DOUBLE_COMPLEX,pzr_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierror)
  call MPI_Gather(pzl_1d,ndim,MPI_DOUBLE_COMPLEX,pzl_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierror)
  call MPI_Gather(gfpq_1d,ndimgo,MPI_DOUBLE_COMPLEX,gfpq_1d_all,ndimgo,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierror)
  if (isra == 1) then
    call MPI_Gather(fzr_1d,ndim,MPI_DOUBLE_COMPLEX,fzr_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierror)
    call MPI_Gather(fzl_1d,ndim,MPI_DOUBLE_COMPLEX,fzl_1d_all,ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierror) 
  end if

  if (my_rank == 0) then
 
!   Regular coefficients
    do k  = 1, ndim_all
      i  = i1d_4d(1,k)
      j  = i1d_4d(2,k)    
      ia = i1d_4d(3,k)
      ih = iasusc(ia)
      ie = i1d_4d(4,k)   
      ! Put back to the coefficients         
      pzr(i,j,ih,ie) = pzr_1d_all(k) 
      pzl(i,j,ih,ie) = pzl_1d_all(k)
      if (isra == 1) then
        pzr(i,j,ih,ie) = pzr_1d_all(k)
        pzl(i,j,ih,ie) = pzl_1d_all(k)
      end if 
    end do ! k

!   Gonsite
    do k  = 1, ndimgo_all
      i  = i1d_4dgo(1,k)
      j  = i1d_4dgo(2,k)    
      ia = i1d_4dgo(3,k)
      ie = i1d_4dgo(4,k)
      ih = iasusc(ia)
      ! Put back onsite gf         
      gfpq(i,j,ih,ie) = gfpq_1d_all(k) 
      if (isra == 1) then
        !!!!!!!!!!!! THERE WILL BE 3 MORE BUT I MUSR PROJECT THEM FIRST !!!!!!!!!!!!!!! 
      end if
    end do ! k

    ! lms blocks
    lmsize = 2*(lmaxd+1)**2 
    do ie = 1, ielast 
      if (config%nsra == 1) ek_susc = sqrt(ez(ie))                                          
      if (config%nsra == 2) ek_susc = sqrt(ez(ie) + (ez(ie)/cvlight)*(ez(ie)/cvlight))      
      call out_coeffs_soc(lmsize,ie,ez(ie),ek_susc,wz(ie))
    end do

  end if ! my_rank 

! All done
  end subroutine read_proj
