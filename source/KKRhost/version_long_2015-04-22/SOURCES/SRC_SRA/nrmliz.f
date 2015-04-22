      subroutine nrmliz(n,r,rn) 
!-  normalizes a vector                                                 
! ----------------------------------------------------------------------
!i Inputs                                                               
!i   n     :number of vectors                                           
!i   r     :vector                                                      
!o Outputs:                                                             
!o   rn    :normalized vector                                           
! ----------------------------------------------------------------------
      implicit none 
! Passed parameters:                                                    
      integer n 
      double precision r(3,*),rn(3,*) 
! Local parameters                                                      
      integer i 
      double precision d,d2 
! External calls                                                        
      external dcopy,dscal 
                                                                        
      call dcopy(3*n,r,1,rn,1) 
      do i=1,n 
        d2=r(1,i)*r(1,i)+r(2,i)*r(2,i)+r(3,i)*r(3,i) 
        d=dsqrt(d2) 
        if (d.ne.0.d0) call dscal(3,1.d0/d,rn(1,i),1) 
      enddo 
                                                                        
      END                                           
