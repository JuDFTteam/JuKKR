      subroutine dscal1(n,da,dx,incx) 
!- Scales a vector by a constant  dx(i) -> a * dx(i)                    
! ----------------------------------------------------------------------
!i Inputs:                                                              
!i   n     :lenght of dx and dy                                         
!i   da    :constant                                                    
!i   dx    :vector                                                      
!i   incx  :incrementation for x                                        
!o Outputs:                                                             
!o   dx    :vector                                                      
!r Remarks:                                                             
!r   Adapted from: jack dongarra, linpack, 3/11/78.                     
! ----------------------------------------------------------------------
      implicit none 
! Passed parameters:                                                    
      double precision da,dx(*) 
      integer incx,n 
! Local parameters:                                                     
      integer i,m,mp1,nincx 
!                                                                       
      if( n.le.0 .or. incx.le.0 )return 
      if(incx.ne.1) then 
! ----- code for increment not equal to 1                               
        nincx = n*incx 
        do i = 1,nincx,incx 
          dx(i) = da*dx(i) 
        enddo 
      else 
! ----- code for increment equal to 1                                   
        m = mod(n,5) 
        if( m .ne. 0 ) then 
          do i = 1,m 
            dx(i) = da*dx(i) 
          enddo 
          if( n .lt. 5 ) return 
        endif 
        mp1 = m + 1 
        do i = mp1,n,5 
          dx(i) = da*dx(i) 
          dx(i + 1) = da*dx(i + 1) 
          dx(i + 2) = da*dx(i + 2) 
          dx(i + 3) = da*dx(i + 3) 
          dx(i + 4) = da*dx(i + 4) 
        enddo 
      endif 
      END                                           
