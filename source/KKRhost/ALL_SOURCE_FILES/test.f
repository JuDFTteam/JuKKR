C *********************************************************** 17.05.91 **
      LOGICAL FUNCTION TEST(STRING)                                    
C ***********************************************************************
C                                                                      
C     TEST = 'STRING  ' IS CONTAINED IN /TESTC/.                      
C                                                                    
C ------------------------------------------------------------------------
      use mod_wunfiles, only: t_params
      IMPLICIT NONE
C                                                                    
      integer I
      character*8      STRING, TESTC(32)

      TESTC = T_PARAMS%TESTC
C
      TEST=.FALSE.                                                    
      DO I=1,32
        IF(STRING.EQ.TESTC(I)) TEST=.TRUE.                         
      END DO
      RETURN                                                       
      END                                                         
