C ***********************************************************************
      LOGICAL FUNCTION OPT(STRING)                                    
C ***********************************************************************
C                                                                      
C     OPT = 'STRING  ' IS CONTAINED IN /OPTC/.                        
C                                                                       
C ------------------------------------------------------------------------
      use mod_wunfiles, only: t_params
      IMPLICIT NONE
C                                                                    
      integer I
      character*8      STRING, OPTC(32)

      OPTC = t_params%OPTC
C                                                                  
      OPT=.FALSE.                                                     
      DO I=1,32
        IF(STRING.EQ.OPTC(I)) OPT=.TRUE.
      END DO
      RETURN                                                       
      END                                                         
