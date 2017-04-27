  subroutine structural_gf(ie,natoms,nlms,nalms,alms2i,dtmat,gstruct)
! new structural GF: (1 - Gold.dtmat).Gnew = Gold
! thread-safe if zgesv is thread-safe
  use global, only: i4b, r8b, c8b!, iodb, i2lms

  implicit none

! Energy point
  integer(kind=i4b), intent(in)    :: ie  
! Number of atoms
  integer(kind=i4b), intent(in)    :: natoms
! Maximum combined dimension of angular momentum and spin
  integer(kind=i4b), intent(in)    :: nlms
! Size of structural GF
  integer(kind=i4b), intent(in)    :: nalms
! Pointers to storage
  integer(kind=i4b), intent(in)    :: alms2i(nlms,natoms)
! Change in t-matrix from change in potential
  complex(kind=c8b), intent(in)    :: dtmat(nlms,nlms,natoms)
! Old structural GF (overwritten)
  complex(kind=c8b), intent(inout) :: gstruct(nalms,nalms)
! ----------------------------------------------------------------------
! complex parameters
  complex(kind=c8b), parameter :: cone   = ( 1.d0, 0.d0)
  complex(kind=c8b), parameter :: cminus = (-1.d0, 0.d0)
  complex(kind=c8b), parameter :: czero  = ( 0.d0, 0.d0)
! permutation array
  integer(kind=i4b), allocatable :: ipiv(:)
! work array
  complex(kind=c8b), allocatable :: work(:,:), work2(:,:)
! misc
  integer(kind=i4b) :: i, j, k, ia, ja, is, js, ilms, jlms, klms, info


!  allocate(work(nalms,nalms),work2(nalms,nalms),ipiv(nalms))
!! --> first iteration                                       
!  work = czero                                              
!  do ia=1,natoms                                            
!    do jlms=1,nlms                                          
!      j = alms2i(jlms,ia)                                   
!      do ilms=1,nlms                                        
!        i = alms2i(ilms,ia)                                 
!        work(i,j) = dtmat(ilms,jlms,ia)                     
!      end do                                                
!    end do                                                  
!  end do                                                    
!! Y = dt.G0                                                 
!  work2 = matmul(work,gstruct)                              
!! X = Y + Y.Y = dt.G0 + dt.G0.dt.G0                         
!  work = work2 + matmul(work2,work2)                        
!! Y = G0.Y                                                  
!  work2 = matmul(gstruct,work)                              
!! G = G0 + Y                                                
!  gstruct = gstruct + work2                                 
!  deallocate(work,work2,ipiv)                               
!  return                                                    
  allocate(work(nalms,nalms),ipiv(nalms))
! --> construct matrix to be inverted
! columns
  do ja=1,natoms
  do jlms=1,nlms
    j = alms2i(jlms,ja)
!   rows
    do ia=1,natoms
    do ilms=1,nlms
      i = alms2i(ilms,ia)
!     ------------------------------------------------------
      work(i,j) = czero
      if (i == j) work(i,j) = cone
!     matrix multiplication
      do klms=1,nlms
        k = alms2i(klms,ja)
        work(i,j) = work(i,j) - gstruct(i,k)*dtmat(klms,jlms,ja)
      end do
!     ------------------------------------------------------
    end do
    end do
  end do
  end do
! --> solve Dyson equation
  call zgesv(nalms,nalms,work,nalms,ipiv,gstruct,nalms,info)
  if (info /= 0) then
    write(*,*) "struct Dyson failed at ie=", ie
    stop
  end if
!  work = gstruct - transpose(gstruct)
!  do ja=1,natoms
!  do jlms=1,nlms
!    j = alms2i(jlms,ja)
!    js = i2lms(2,jlms)
!    do ia=1,natoms
!    do ilms=1,nlms
!      i = alms2i(ilms,ia)
!      is = i2lms(2,ilms)
!      if (is == js .and. abs(work(i,j)) > 1.d-6) write(iodb,'("structural_gf: asym ia, ja, ilms, jlms, diff=",4i8,2es16.8)') ia, ja, ilms, jlms, work(i,j)
!    end do
!    end do
!  end do
!  end do
  deallocate(work,ipiv)
! All done! 
  end subroutine structural_gf
