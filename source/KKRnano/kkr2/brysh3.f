C>    @param Y single vector containing all lm,spin components
C>    @param X spherical potential
C>    @param Z non-spherical potential
c ************************************************************************
      SUBROUTINE BRYSH3(Y,X,Z,IRMIN,IRC,IATOM,
     &                  NSPIN,IMAP,LMPOT,
C                       new parameters after inc.p removal
     &                  IRMD, IRNSD)
c*********************************************************************
c     shifts the density or potential of all mt-cell into one single
c     vector and projects out the coulomb part only.
c
c                                    s. bluegel , kfa , 1987
c
c     modified for parallelization
c                                    a. thiess , jun 2008   
c ------------------------------------------------------------------------
      IMPLICIT NONE

C      INTEGER LMPOTD
C      PARAMETER (LMPOTD= (LPOTD+1)**2)
C      INTEGER IRMIND
C      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..

      INTEGER IRMD
      INTEGER IRNSD

C     .. Scalar Arguments ..
      INTEGER IMAP,LMPOT,IATOM,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(IRMD,*),Y(*),Z((IRMD-IRNSD):IRMD,LMPOT,*)
      INTEGER IRC(*),IRMIN(*)
C     ..
C     .. Local Scalars ..
      INTEGER IP,IR,IRC1,IRMIN1,IS,LM
C     ..
      IMAP = 0
      DO 50 IS = 1,NSPIN
C        DO 40 IA = NATPS,NAEZ
          IP = IS
          IRC1 = IRC(IATOM)
          DO 10 IR = 1,IRC1
            IMAP = IMAP + 1
            Y(IMAP) = X(IR,IP)
   10     CONTINUE
c
          IF (LMPOT.GT.1) THEN
            IRMIN1 = IRMIN(IATOM)
            DO 30 LM = 2,LMPOT
              DO 20 IR = IRMIN1,IRC1
                IMAP = IMAP + 1
                Y(IMAP) = Z(IR,LM,IP)
   20         CONTINUE
   30       CONTINUE
          END IF
c
C   40   CONTINUE
   50 CONTINUE
c
      END
