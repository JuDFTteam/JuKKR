C ***********************************************************************
      LOGICAL FUNCTION OPT(STRING)                                    
C ***********************************************************************
C                                                                      
C     OPT = 'STRING  ' IS CONTAINED IN /OPTC/.                        
C                                                                       
C ------------------------------------------------------------------------
C                                                                      
      use common_optc
C                                                                    
      character*8      STRING
C                                                                  
C                                                                      
      OPT=.FALSE.                                                     
      DO 1 I=1,8                                                     
        IF(STRING.EQ.OPTC(I)) OPT=.TRUE.
 1    END DO
      RETURN                                                       
      END                                                         
