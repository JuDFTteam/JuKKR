C ***********************************************************************
      LOGICAL FUNCTION OPT(STRING)                                    
C ***********************************************************************
C                                                                      
C     OPT = 'STRING  ' IS CONTAINED IN /OPTC/.                        
C                                                                       
C ------------------------------------------------------------------------
C                                                                      
      COMMON/OPTC/  OPTC(32)                                           
      save  /optc/
C                                                                    
      character*8      STRING   ,OPTC                               
C                                                                  
C                                                                      
      OPT=.FALSE.                                                     
      DO 1 I=1,32                                                     
        IF(STRING.EQ.OPTC(I)) OPT=.TRUE.
 1    END DO
      RETURN                                                       
      END                                                         
