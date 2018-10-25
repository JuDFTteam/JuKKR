      logical function latvec(n,qlat,vec) 
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
! Passed parameters:                                                    
      integer n 
      double precision qlat(3,*),vec(3,*) 
! Local parameters:                                                     
      integer i,m 
      double precision tol,vdiff 
      parameter(tol=1.d-3) 
! Common block:                                                         
! Intrinsic functions:                                                  
      intrinsic  dabs,dnint 
                                                                        
      latvec=.false. 
      do i=1,n 
       do m=1,3 
         vdiff=vec(1,i)*qlat(1,m)+vec(2,i)*qlat(2,m)+vec(3,i)*qlat(3,m) 
         vdiff=dabs(vdiff-dnint(vdiff)) 
         if (vdiff.gt.tol) return 
       enddo 
      enddo 
      latvec=.true. 
      END                                           
