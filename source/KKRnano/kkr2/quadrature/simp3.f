      SUBROUTINE SIMP3(F,FINT,ISTART,IEND,DRDI)
c-----------------------------------------------------------------------
c     this subroutine does an integration from istart to iend of
c     the real function f with an extended 3-point-simpson :
c
c                          r(istart)
c
c                       fint = { f(r') dr'
c
c                           r(iend)
c
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION FINT
      INTEGER IEND,ISTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(*),F(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      A1 = 4.0D0/3.0D0
      A2 = 2.0D0/3.0D0
c
c---> initialize fint
c
      IF (MOD(IEND-ISTART,2).EQ.0) THEN
        FINT = F(ISTART)*DRDI(ISTART)/3.0D0
        IST = ISTART + 1

      ELSE
        FINT = (F(ISTART+3)*DRDI(ISTART+3)-
     +         5.0D0*F(ISTART+2)*DRDI(ISTART+2)+
     +         19.0D0*F(ISTART+1)*DRDI(ISTART+1)+
     +         9.0D0*F(ISTART)*DRDI(ISTART))/24.0D0 +
     +         F(ISTART+1)*DRDI(ISTART+1)/3.0D0
        IST = ISTART + 2
      END IF
c
c---> calculate with an extended 3-point-simpson
c
      DO 10 I = IST,IEND - 1,2
        FINT = FINT + A1*F(I)*DRDI(I) + A2*F(I+1)*DRDI(I+1)
   10 CONTINUE
      FINT = FINT - F(IEND)*DRDI(IEND)/3.0D0
c
      END
