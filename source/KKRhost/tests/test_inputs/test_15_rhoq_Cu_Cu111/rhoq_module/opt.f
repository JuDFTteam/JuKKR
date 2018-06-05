C ***********************************************************************
      LOGICAL FUNCTION OPT(STRING)   
      implicit none
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
      integer i
C                                                                  
      OPT=.FALSE.                                                     
      DO I=1,32
        IF(STRING.EQ.OPTC(I)) OPT=.TRUE.
      END DO
      RETURN                                                       
      END                                                         
