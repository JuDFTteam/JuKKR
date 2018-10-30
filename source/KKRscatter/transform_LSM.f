C ************************************************************************
      SUBROUTINE TRANSFORM_LSM(I1,LMAX,LMMAXD,LSM,LSM_REL,
     +                                        EZ,LSM_CHECK)

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: LMAX,I1,LMMAXD
      LOGICAL,INTENT(IN) :: LSM_CHECK

      COMPLEX,INTENT(IN) :: !TMATLL(2*LMMAXD,2*LMMAXD),
     +                      LSM(2*LMMAXD,2*LMMAXD),EZ!,
c     +                      PNS(2*LMMAXD,2*LMMAXD)

c local variables
      INTEGER            :: LM2,LM1,NLMAX,NKMMAX,NMUEMAX,
     +                      NKMPMAX,NKMAX,LINMAX 

      COMPLEX            :: LSM_REL(2*LMMAXD,2*LMMAXD),
     +                      RC(2*LMMAXD,2*LMMAXD),
     +                      CREL(2*LMMAXD,2*LMMAXD),
     +                      RREL(2*LMMAXD,2*LMMAXD),
     +                      LSM_TEST(2*LMMAXD,2*LMMAXD),
     +                      LSM_CLM(2*LMMAXD,2*LMMAXD)


      include 'inc.p'


      OPEN (UNIT=100,file="LSM_REL",form="formatted")
      OPEN (UNIT=101,file="LSM_CBAS",form="formatted")
      OPEN (UNIT=102,file="LSM_back",form="formatted")


c  transformation of the basis of real spherical harmonics to the basis
c  of complex spherical harmonics. 
      NLMAX     = LMAXD+1 
      NKMMAX    = 2*NLMAX**2
      NMUEMAX   = 2*NLMAX
      NKMPMAX   = (NKMMAX+2*NLMAX)
      NKMAX     = 2*NLMAX-1
      LINMAX    = (2*NLMAX*(2*NLMAX-1))

c      WRITE(6,*) "before DRVBASTRANS"
      CALL DRVBASTRANS(RC,CREL,RREL,
     &                       NLMAX,NKMMAX,NMUEMAX,NKMPMAX,NKMAX,LINMAX)
c      WRITE(6,*) "after DRVBASTRANS"

c  transformation of the basis of real spherical harmonics to the basis
c  of complex spherical harmonics. 

      LSM_REL=0d0

      CALL CHANGEREP(
     &         LSM(1,1),'RLM>REL',LSM_REL(1,1),
     +         NSPD*LMMAXD,NSPD*LMMAXD,
     &         RC,CREL,RREL,'L_REL ',6)

      WRITE(6,*) "after CHANGEREP LSM to rel basis" 

      LSM_TEST=0d0

      CALL CHANGEREP(LSM_REL(1,1),'REL>RLM',LSM_TEST(1,1),
     +         NSPD*LMMAXD,NSPD*LMMAXD,RC,CREL,RREL,'L_RLM ',6)

      WRITE(6,*) "after CHANGEREP LSM back"

      LSM_CLM=0d0

      CALL CHANGEREP(
     &         LSM(1,1),'RLM>CLM',LSM_CLM(1,1),
     +         NSPD*LMMAXD,NSPD*LMMAXD,
     &         RC,CREL,RREL,'L_CLM ',6)

      WRITE(6,*) "after CHANGEREP LSM to complex basis" 


      DO LM1 = 1,NSPD*LMMAXD
         DO LM2 = 1,NSPD*LMMAXD
            WRITE(100,"((e17.9),(2I5),(2e17.9))")  DREAL(EZ),
     +                                LM1,LM2,LSM_REL(LM1,LM2)
            WRITE(101,"((e17.9),(2I5),(2e17.9))")  DREAL(EZ),
     +                                LM1,LM2,LSM_CLM(LM1,LM2)
            WRITE(102,"((e17.9),(2I5),(2e17.9))") DREAL(EZ),
     +                                LM1,LM2,LSM_TEST(LM1,LM2)
         ENDDO
      ENDDO

      CLOSE(100)
      CLOSE(101)
      CLOSE(102)


 9876 FORMAT(2E16.8,3I4,2E16.8)
 9877 FORMAT(20E16.8)


      END SUBROUTINE
