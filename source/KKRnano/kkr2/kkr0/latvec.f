      logical function latvec(n, qlat, vec) 
!- Checks if a set of vectors are lattice vectors                       
! ----------------------------------------------------------------------
!i Inputs:                                                              
!i   n     :number of vectors                                           
!i   qlat  :primitive translation vectors in reciprocal space           
!i   vec   :double-precision vector                                     
!o Outputs:                                                             
!o   latvec:.true. if all vectors are lattice vectors                   
!r Remarks:                                                             
! ----------------------------------------------------------------------
      implicit none 
      integer, intent(in) :: n 
      double precision, intent(in) :: qlat(3,*), vec(3,*) 

      integer :: i, m 
      double precision :: vdiff 
      double precision, parameter :: tol=1.d-3 
                                                                        
      latvec = .false. 
      do i = 1, n 
        do m = 1, 3 
           vdiff = vec(1,i)*qlat(1,m) + vec(2,i)*qlat(2,m)
     &           + vec(3,i)*qlat(3,m) 
           vdiff = dabs(vdiff - nint(vdiff)) 
           if (vdiff > tol) return ! false
        enddo ! m
      enddo ! i
      latvec = .true. 
      endfunction latvec

      