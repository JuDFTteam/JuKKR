c ************************************************************************
      SUBROUTINE BRYSH2(Y,X,IRMIN,IRC,IATOM,
     &                  NAEZ,NSPIN,IMAP,LMPOT,
c                       new parameter after inc.p removal
     &                  IRMD)
c*********************************************************************
c     maps the density or potential back from one single vector into
c     the proper bins of each single mt-cell . the magnetization
c     density is also added in.
c                                    s. bluegel , kfa , 1987
c
c ------------------------------------------------------------------------
      IMPLICIT NONE
C     .. Parameters ..
c      INTEGER IRMD,LPOTD
c      PARAMETER (IRMD=424,LPOTD=8)
C      INTEGER LMPOTD
C      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..

      INTEGER IRMD

C     .. Scalar Arguments ..
      INTEGER IMAP,LMPOT,NATPS,NAEZ,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(IRMD,LMPOT,*),Y(*)
      INTEGER IRC(*),IRMIN(*)
C     ..
C     .. Local Scalars ..
      INTEGER IP,IR,IRC1,IRMIN1,IS,LM
      INTEGER IATOM
C     ..
      IMAP = 0

      DO 50 IS = 1,NSPIN
C        DO 40 IA = NATPS,NAEZ
          IP = IS
          IRC1 = IRC(IATOM)
          DO 10 IR = 1,IRC1
            IMAP = IMAP + 1
            X(IR,1,IP) = Y(IMAP)
   10     CONTINUE
c
          IF (LMPOT.GT.1) THEN
            IRMIN1 = IRMIN(IATOM)
            DO 30 LM = 2,LMPOT
              DO 20 IR = IRMIN1,IRC1
                IMAP = IMAP + 1
                X(IR,LM,IP) = Y(IMAP)
   20         CONTINUE
   30       CONTINUE
          END IF
c
C   40   CONTINUE
   50 CONTINUE
c
      END
