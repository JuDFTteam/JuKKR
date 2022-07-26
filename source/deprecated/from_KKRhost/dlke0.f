c 04.10.95 *************************************************************
      SUBROUTINE DLKE0(GLLKE,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,
     +                 BZKP,RCLS,GINP)
c **********************************************************************
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAX+1)**2)
      INTEGER ALMGF0
      PARAMETER (ALMGF0=LMGF0D*NAEZD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER NAEZ,NACLSMAX
C     ..
C     .. Array Arguments ..
c
c
      DOUBLE COMPLEX GINP(LMGF0D*NACLSMAX,LMGF0D,*),GLLKE(ALMGF0,*)
      DOUBLE PRECISION BZKP(*),RCLS(3,NACLSD,*),RR(3,0:NRD)
      INTEGER ATOM(NACLSD,*),CLS(*),EZOA(NACLSD,*),NACLS(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IC,IM,J,JN,M,N
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX GLLKE1(ALMGF0,LMGF0D)
      DOUBLE PRECISION KP(6)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,DLKE1
C     ..
C     .. Save statement ..
      SAVE
C     ..
c      write(6,*) '>>> DLKE0 : Fourier-transforms the ',
c     +           'GF of reference system'
c ----------------------------------------------------------------------

C     .. External Functions ..
      LOGICAL OPT
      EXTERNAL OPT
C     ..
      CALL CINIT(ALMGF0*ALMGF0,GLLKE(1,1))

      DO 30 I = 1,NAEZ


        KP(1) = BZKP(1)
        KP(2) = BZKP(2)
        KP(3) = BZKP(3)
        IF (OPT('COMPLEX ')) THEN
          KP(4) = BZKP(4)
          KP(5) = BZKP(5)
          KP(6) = BZKP(6)
        END IF

        IC = CLS(I)
        CALL DLKE1(GLLKE1,ALAT,NACLS,NACLSMAX,RR,EZOA(1,I),ATOM(1,I),KP,
     +             IC,GINP(1,1,IC),RCLS(1,1,IC))

        DO 20 M = 1,LMGF0D
          IM = (I-1)*LMGF0D + M
          DO 10 JN = 1,LMGF0D*NAEZ
            GLLKE(JN,IM) = GLLKE(JN,IM) + GLLKE1(JN,M)
   10     CONTINUE
   20   CONTINUE



   30 CONTINUE

c ----------------------------------------------------------------------
      IF (OPT('symG(k) ')) THEN
c
c -->   symmetrization
c
        DO 70 I = 1,NAEZ

          KP(1) = -BZKP(1)
          KP(2) = -BZKP(2)
          KP(3) = -BZKP(3)
          IF (OPT('COMPLEX ')) THEN
            KP(4) = -BZKP(4)
            KP(5) = -BZKP(5)
            KP(6) = -BZKP(6)
          END IF

          IC = CLS(I)
          CALL DLKE1(GLLKE1,ALAT,NACLS,NACLSMAX,RR,EZOA(1,I),ATOM(1,I),
     +               KP,IC,GINP(1,1,IC),RCLS(1,1,IC))

          DO 60 J = 1,NAEZ
            DO 50 M = 1,LMGF0D
              IM = (I-1)*LMGF0D + M
              DO 40 N = 1,LMGF0D
                JN = (J-1)*LMGF0D + N
                GLLKE(IM,JN) = (GLLKE(IM,JN)+GLLKE1(JN,M))/2.0D0
   40         CONTINUE
   50       CONTINUE
   60     CONTINUE

   70   CONTINUE

      END IF
c ----------------------------------------------------------------------

      RETURN

      END
