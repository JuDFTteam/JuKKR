C ************************************************************************
      SUBROUTINE TRANSFORM_TMAT(I1,LMAX,LMMAXD,TMATLL,LSM,PNS,
     +                                        IE,IELAST,EZ,LSM_CHECK)

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: LMAX,IE,I1,LMMAXD,IELAST
      LOGICAL,INTENT(IN) :: LSM_CHECK

      COMPLEX,INTENT(IN) :: !TMATLL(2*LMMAXD,2*LMMAXD),
     +                      LSM(2*LMMAXD,2*LMMAXD),EZ,
     +                      PNS(2*LMMAXD,2*LMMAXD)

c local variables
      INTEGER            :: LM2,LM1,NLMAX,NKMMAX,NMUEMAX,
     +                      NKMPMAX,NKMAX,LINMAX 

      COMPLEX            :: LSM_REL(2*LMMAXD,2*LMMAXD),
     +                      RC(2*LMMAXD,2*LMMAXD),
     +                      CREL(2*LMMAXD,2*LMMAXD),
     +                      RREL(2*LMMAXD,2*LMMAXD),
     +                      TMAT_REL(2*LMMAXD,2*LMMAXD),
     +                      TMAT_TEST(2*LMMAXD,2*LMMAXD),
     +                      LSM_TEST(2*LMMAXD,2*LMMAXD),
     +                      PNS_REL(2*LMMAXD,2*LMMAXD),
     +                      PNS_TEST(2*LMMAXD,2*LMMAXD),
     +                      TMATLL(2*LMMAXD,2*LMMAXD)


      include 'inc.p'


      IF (I1==1 .AND. IE==1) THEN 

        OPEN (UNIT=93,file="TMATLL_Kappa_mu",form="formatted")
        OPEN (UNIT=95,file="TMATLL_rel_abs",form="formatted")
        OPEN (UNIT=96,file="TMATLL_rel_s",form="formatted")
        OPEN (UNIT=97,file="TMATLL_rel_p",form="formatted")
        OPEN (UNIT=98,file="TMATLL_rel_d",form="formatted")

        IF (LSM_CHECK == .TRUE.) THEN
          OPEN (UNIT=92,file="TMATLL_ORIG",form="formatted")
          OPEN (UNIT=94,file="TMATLL_ruecktrafo",form="formatted")
c          OPEN (UNIT=99,file="LSM_orig",form="formatted")
          OPEN (UNIT=100,file="LSM_trans",form="formatted")
          OPEN (UNIT=101,file="LSM_back",form="formatted")
          OPEN (UNIT=102,file="PNS_trans",form="formatted")
          OPEN (UNIT=103,file="PNS_back",form="formatted")
        END IF

      END IF

      IF (LSM_CHECK == .TRUE.) THEN
        DO LM1=1,NSPD*LMMAXD
          DO LM2=1,NSPD*LMMAXD
            WRITE(92,"((e17.9),(2I5),(4e17.9))") REAL(EZ),LM2,LM1,
     +                            TMATLL(LM2,LM1),TMATLL(LM1,LM2)
          END DO
        END DO
      END IF


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

      IF (LSM_CHECK == .TRUE.) THEN

        LSM_REL=0d0

        CALL CHANGEREP(
     &           LSM(1,1),'RLM>REL',LSM_REL(1,1),
     +           NSPD*LMMAXD,NSPD*LMMAXD,
     &           RC,CREL,RREL,'L_REL ',6)

        WRITE(6,*) "after CHANGEREP LSM"

        LSM_TEST=0d0

        CALL CHANGEREP(LSM_REL(1,1),'REL>RLM',LSM_TEST(1,1),
     +           NSPD*LMMAXD,NSPD*LMMAXD,RC,CREL,RREL,'L_RLM ',6)

        WRITE(6,*) "after CHANGEREP LSM back"

        PNS_REL=0d0

        CALL CHANGEREP(
     &           PNS(1,1),'RLM>REL',PNS_REL(1,1),
     +           NSPD*LMMAXD,NSPD*LMMAXD,
     &           RC,CREL,RREL,'P_REL ',6)

        WRITE(6,*) "after CHANGEREP PNS to P_REL"

        PNS_TEST=0d0

        CALL CHANGEREP(PNS_REL(1,1),'REL>RLM',PNS_TEST(1,1),
     +           NSPD*LMMAXD,NSPD*LMMAXD,RC,CREL,RREL,'P_RLM ',6)

        WRITE(6,*) "after CHANGEREP PNS back"
      END IF

      TMAT_REL=0d0
      TMAT_TEST=0d0

      CALL CHANGEREP(TMATLL(1,1),'RLM>REL',TMAT_REL(1,1),
     +           NSPD*LMMAXD,NSPD*LMMAXD,RC,CREL,RREL,'T_REL ',6)

c        WRITE(6,*) "after CHANGEREP RLM>REL"

c Relativistic kappa-mu to real l-m basis
c (should be equal to TMATLL1, so we do it for test reasons)

      IF (LSM_CHECK == .TRUE.) THEN

        CALL CHANGEREP(TMAT_REL(1,1),'REL>RLM',TMAT_TEST(1,1),
     +           NSPD*LMMAXD,NSPD*LMMAXD,RC,CREL,RREL,'T_RLM ',6)


        DO LM1 = 1,NSPD*LMMAXD
           DO LM2 = 1,NSPD*LMMAXD
              WRITE(94,"((e17.9),(3I5),(4e17.9))")  DREAL(EZ),
     +                I1,LM1,LM2,TMAT_TEST(LM1,LM2),
     +                TMATLL(LM1,LM2)-TMAT_TEST(LM1,LM2)
              WRITE(100,"((e17.9),(2I5),(2e17.9))")  DREAL(EZ),
     +                                  LM1,LM2,LSM_REL(LM1,LM2)
              WRITE(101,"((e17.9),(2I5),(2e17.9))") DREAL(EZ),
     +                                  LM1,LM2,LSM_TEST(LM1,LM2)
              WRITE(102,"((e17.9),(2I5),(2e17.9))")  DREAL(EZ),
     +                                  LM1,LM2,PNS_REL(LM1,LM2)
              WRITE(103,"((e17.9),(2I5),(4e17.9))") DREAL(EZ),
     +                  LM1,LM2,PNS_TEST(LM1,LM2),PNS_TEST(LM2,LM1)
           ENDDO
        ENDDO

      END IF

c     DO LM1 = 1,NSPD*LMMAXD
c        DO LM2 = 1,NSPD*LMMAXD
c          TMATLL(LM2,LM1)=TMAT_TEST(LM2,LM1)
c        END DO
c     END DO

      DO LM1 = 1,NSPD*LMMAXD
        DO LM2 = 1,NSPD*LMMAXD
          WRITE(93,"((e17.9),(3I5),(2e17.9))") DREAL(EZ),
     +                          I1,LM1,LM2,TMAT_REL(LM1,LM2)
        ENDDO
      ENDDO

      WRITE(95,9877) DREAL(EZ),
     &     ABS(TMAT_REL(1,1)),
     &     ABS(TMAT_REL(2,2)),
     &     ABS(TMAT_REL(3,3)),
     &     ABS(TMAT_REL(4,4)),
     &     ABS(TMAT_REL(5,5)),
     &     ABS(TMAT_REL(6,6)),
     &     ABS(TMAT_REL(7,7)),
     &     ABS(TMAT_REL(8,8)),
     &     ABS(TMAT_REL(9,9)),
     &     ABS(TMAT_REL(10,10)),
     &     ABS(TMAT_REL(11,11)),
     &     ABS(TMAT_REL(12,12)),
     &     ABS(TMAT_REL(13,13)),
     &     ABS(TMAT_REL(14,14)),
     &     ABS(TMAT_REL(15,15)),
     &     ABS(TMAT_REL(16,16)),
     &     ABS(TMAT_REL(17,17)),
     &     ABS(TMAT_REL(18,18))
      WRITE(96,9877) DREAL(EZ),
     &     (TMAT_REL(1,1))
c     &     (TMAT_REL(2,2))
      WRITE(97,9877) DREAL(EZ),
     &     (TMAT_REL(3,3)),
c     &     (TMAT_REL(4,4)),
     &     (TMAT_REL(5,5))
c     &     (TMAT_REL(6,6)),
c     &     (TMAT_REL(7,7)),
c     &     (TMAT_REL(8,8))
      WRITE(98,9877) DREAL(EZ),
     &     (TMAT_REL(9,9)),
c     &     (TMAT_REL(10,10)),
c     &     (TMAT_REL(11,11)),
c     &     (TMAT_REL(12,12)),
     &     (TMAT_REL(13,13))
c     &     (TMAT_REL(14,14)),
c     &     (TMAT_REL(15,15)),
c     &     (TMAT_REL(16,16)),
c     &     (TMAT_REL(17,17)),
c     &     (TMAT_REL(18,18))

      IF (I1==NATYPD .AND. IE==IELAST) THEN 

        CLOSE(93)
        CLOSE(95)
        CLOSE(96)
        CLOSE(97)
        CLOSE(98)

        IF (LSM_CHECK == .TRUE.) THEN
          CLOSE(92)
          CLOSE(94)
c         CLOSE(99)
          CLOSE(100)
          CLOSE(101)
          CLOSE(102)
          CLOSE(103)
        END IF

      END IF

 9876 FORMAT(2E16.8,3I4,2E16.8)
 9877 FORMAT(20E16.8)


      END SUBROUTINE
