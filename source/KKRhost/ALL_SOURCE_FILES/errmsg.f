      subroutine errmsg(messg,isev) 
!- Write error message to message error device                          
! ----------------------------------------------------------------------
!i Inputs:                                                              
!i   messg :error message                                               
!i   isev  :severity level                                              
!r Remarks                                                              
!r   if severity level greater or equal than error tolerance            
!r   program will stop.                                                 
! ----------------------------------------------------------------------
      implicit none 
! Passed parameters:                                                    
      integer isev 
      character*(*) messg 
! Local parameters:                                                     
      integer iline,iisev,ipos(0:20),l,nline,nunit                                                     
      character*14 c(1:4) 
! External calls:                                                       
! Intrinsic functions                                                   
      intrinsic iabs,max0,min0 
                                                                        
      data c/'Information:','Warning:','Error:','Fatal error:'/ 
                                                                        
      iisev=max0(min0(isev,4),1) 
                                                                        
      ipos(0)=1 
      nline=0 
      l = 1 
      do while (messg(l:l).ne.'$'.and.l.lt.500) 
        if (messg(l:l).eq.'|') then 
          nline=nline+1 
          ipos(nline)=l 
        endif 
        l=l+1 
      enddo 
      nline=nline+1 
      ipos(nline)=l 
                                                                        
      nunit=1 
      if (nunit.eq.1) write(*,*) 
      do iline=1,nline 
         write(*,300)c(iisev),messg(ipos(iline-1)+1:ipos(iline)-1) 
      enddo 
                                                                        
      if (iabs(isev).ge.3) stop
                                                                        
  300 format(a13,a) 
                                                                        
      END                                           
