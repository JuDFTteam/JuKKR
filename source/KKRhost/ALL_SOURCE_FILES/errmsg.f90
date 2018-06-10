SUBROUTINE errmsg(messg,isev)
!- Write error message to message error device
! ----------------------------------------------------------------------
!i Inputs:
!i   messg :error message
!i   isev  :severity level
!r Remarks
!r   if severity level greater or equal than error tolerance
!r   program will stop.
! ----------------------------------------------------------------------
      use mod_types, only: t_inc
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
                                                                        
iisev=MAX0(MIN0(isev,4),1)

ipos(0)=1
nline=0
l = 1
DO WHILE (messg(l:l) /= '$'.AND.l < 500)
  IF (messg(l:l) == '|') THEN
    nline=nline+1
    ipos(nline)=l
  END IF
  l=l+1
END DO
nline=nline+1
ipos(nline)=l

nunit=1
IF ((nunit == 1).AND.(t_inc%i_write>0)) WRITE(1337,*)
IF(t_inc%i_write>0) THEN
  DO iline=1,nline
    WRITE(1337,300)c(iisev),messg(ipos(iline-1)+1:ipos(iline)-1)
  END DO
END IF

IF (IABS(isev) >= 3) STOP

300 FORMAT(a13,a)

END SUBROUTINE errmsg
