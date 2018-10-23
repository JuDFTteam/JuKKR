C*==drvbastrans.f    processed by SPAG 6.05Rc at 14:05 on 19 Sep 2001
      SUBROUTINE DRVBASTRANS(RC,CREL,RREL,SRREL,NRREL,IRREL,
     &                       NLMAX,NKMMAX,NMUEMAX,NKMPMAX,NKMAX,LINMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LINMAX,NKMAX,NKMMAX,NKMPMAX,NLMAX,NMUEMAX
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),
     &           RREL(NKMMAX,NKMMAX),SRREL(2,2,NKMMAX)
      INTEGER IRREL(2,2,NKMMAX),NRREL(2,NKMMAX)
C
C Local variables
C
      REAL*8 CGC(NKMPMAX,2)
      INTEGER I,IKM1LIN(LINMAX),IKM2LIN(LINMAX),IL,IMUE,IPRINT,
     &        KAPTAB(NMUEMAX),LTAB(NMUEMAX),MMAX,NMUETAB(NMUEMAX),
     &        NSOLLM(NLMAX,NMUEMAX)
C
C*** End of declarations rewritten by SPAG
C
      IF (NKMMAX.NE.2*NLMAX**2) 
     &     STOP ' Check NLMAX,NKMMAX in < DRVBASTRANS > '
      IF (NMUEMAX.NE.2*NLMAX) 
     &     STOP ' Check NLMAX,NMUEMAX in < DRVBASTRANS > '
      IF (NKMPMAX.NE.(NKMMAX+2*NLMAX)) 
     &     STOP ' Check NLMAX,NKMMAX,NKMPMAX in < DRVBASTRANS > '
      IF (NKMAX.NE.2*NLMAX-1) 
     &     STOP ' Check NLMAX,NKMAX in < DRVBASTRANS > '
      IF (LINMAX.NE.(2*NLMAX*(2*NLMAX-1)))
     &     STOP ' Check NLMAX,LINMAX in < DRVBASTRANS > '
C
      IPRINT = 0
C
      DO I = 1,NMUEMAX
         LTAB(I) = I/2
         IF ( 2*LTAB(I).EQ.I ) THEN
            KAPTAB(I) = LTAB(I)
         ELSE
            KAPTAB(I) = -LTAB(I) - 1
         END IF
         NMUETAB(I) = 2*ABS(KAPTAB(I))
      END DO
C
      DO IL = 1,NLMAX
         MMAX = 2*IL
         DO IMUE = 1,MMAX
            IF ( (IMUE.EQ.1) .OR. (IMUE.EQ.MMAX) ) THEN
               NSOLLM(IL,IMUE) = 1
            ELSE
               NSOLLM(IL,IMUE) = 2
            END IF
         END DO
      END DO
C
      CALL IKMLIN(IPRINT,NSOLLM,IKM1LIN,IKM2LIN,NLMAX,NMUEMAX,LINMAX,
     &            NLMAX)
C
      CALL CALCCGC(LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX)
C
C ---------------------------- now calculate the transformation matrices
C
      CALL STRSMAT(NLMAX-1,CGC,SRREL,NRREL,IRREL,NKMMAX,NKMPMAX)
C
      CALL BASTRMAT(NLMAX-1,CGC,RC,CREL,RREL,NKMMAX,NKMPMAX)
C
      RETURN
      END
C
