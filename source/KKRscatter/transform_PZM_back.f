C ************************************************************************
      SUBROUTINE TRANSFORM_PZM(I1,LMAX,LMMAXD,PZM,PZMTRANS,
     +                           TMATREL,TMAT,EZ,LSM_CHECK,IRMAX)

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: LMAX,I1,LMMAXD,IRMAX
      LOGICAL,INTENT(IN) :: LSM_CHECK

      COMPLEX,INTENT(IN) :: TMATREL(2*LMMAXD,2*LMMAXD),
     +                      PZM(IRMAX,2*LMMAXD,2*LMMAXD),EZ!,
c     +                      PNS(2*LMMAXD,2*LMMAXD)

c local variables
      INTEGER            :: LM2,LM1,NLMAX,NKMMAX,NMUEMAX,
     +                      NKMPMAX,NKMAX,LINMAX,IR

      COMPLEX            :: PZMTRANS(IRMAX,2*LMMAXD,2*LMMAXD),
     +                      RC(2*LMMAXD,2*LMMAXD),
     +                      CREL(2*LMMAXD,2*LMMAXD),
     +                      RREL(2*LMMAXD,2*LMMAXD),
     +                      LSM_TEST(2*LMMAXD,2*LMMAXD),
     +                      LSM_CLM(2*LMMAXD,2*LMMAXD),
     +                      TMAT(2*LMMAXD,2*LMMAXD)


      include 'inc.p'


      OPEN (UNIT=100,file="PZM_DIAG",form="formatted")
      OPEN (UNIT=101,file="PZM_TRANS",form="formatted")
c      OPEN (UNIT=102,file="LSM_back",form="formatted")


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

      PZMTRANS=0d0

      DO IR=1,IRMAX

        CALL CHANGEREP(
     &         PZM(IR,1:2*LMMAXD,1:2*LMMAXD),'REL>RLM',
     +         PZMTRANS(IR,1:2*LMMAXD,1:2*LMMAXD),
     +         NSPD*LMMAXD,NSPD*LMMAXD,
     &         RC,CREL,RREL,'P_RLM ',6)

c        WRITE(6,*) "after CHANGEREP LSM to rel basis" 

c     LSM_TEST=0d0

c     CALL CHANGEREP(LSM_REL(1,1),'REL>RLM',LSM_TEST(1,1),
c    +         NSPD*LMMAXD,NSPD*LMMAXD,RC,CREL,RREL,'L_RLM ',6)

c     WRITE(6,*) "after CHANGEREP LSM back"

c     LSM_CLM=0d0

c     CALL CHANGEREP(
c    &         LSM(1,1),'RLM>CLM',LSM_CLM(1,1),
c    +         NSPD*LMMAXD,NSPD*LMMAXD,
c    &         RC,CREL,RREL,'L_CLM ',6)

c     WRITE(6,*) "after CHANGEREP LSM to complex basis" 


      ENDDO !IR=1,IRMD

      DO LM1 = 1,NSPD*LMMAXD
        DO LM2 = 1,NSPD*LMMAXD
          DO IR=1,IRMAX
            WRITE(100,"((3I5),(2e17.9))")  LM1,LM2,IR, 
     +                                     PZM(IR,LM2,LM1)
            WRITE(101,"((3I5),(2e17.9))")  LM1,LM2,IR, 
     +                                PZMTRANS(IR,LM2,LM1)
          ENDDO
        ENDDO
      ENDDO

      CLOSE(100)
      CLOSE(101)

      TMAT=0d0

      CALL CHANGEREP(TMATREL,'REL>RLM',TMAT,
     +                 NSPD*LMMAXD,NSPD*LMMAXD,
     &                 RC,CREL,RREL,'P_RLM ',6)

 9876 FORMAT(2E16.8,3I4,2E16.8)
 9877 FORMAT(20E16.8)


      END SUBROUTINE
