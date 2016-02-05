C ************************************************************************
      DOUBLE COMPLEX FUNCTION EXPIDL(T,E,L)
C ************************************************************************
C     
C              1
C     t = - ------- * sin(delta) * exp(i*delta)
C           sqrt(e)
c
c     EXPIDL = exp(i*delta)
C
c ------------------------------------------------------------------------
c     .. arguments
      DOUBLE COMPLEX T,E
      INTEGER L
C     .. External Functions ..
      INTRINSIC DABS,DACOS,DIMAG,DREAL,SQRT,EXP
C     .. locals
      DOUBLE PRECISION D,DD,FAC,SINSQ
      DOUBLE COMPLEX CONE,CI
      PARAMETER(CONE=(1.0D0,0.0D0),CI=(0.0D0,1.0D0))
      DOUBLE PRECISION PI
      PARAMETER (PI= 3.14159265358979312d0)
      DATA SMALL /1.0D-10/
      SAVE
c ------------------------------------------------------------------------
      IF (DABS(DIMAG(E)).GT.SMALL .OR. DREAL(E).LE.0.0D0) THEN 
        EXPIDL = CONE
        EXPIDL = CI**L
        RETURN
      END IF

      SINSQ =-DIMAG(T*SQRT(E))
      IF (SINSQ.LE.SMALL) THEN
        D=DSQRT(SINSQ)
      ELSE
        D=DACOS(DSQRT(1-SINSQ))
      END IF

      DD = D
      IF (DREAL(-T).LE.0.0D0) DD = PI - D 
c      FAC = 1.0D0
c      IF (DIMAG(T).LE.0.0D0) FAC = -1.0D0

c      EXPIDL = ZEXP(CI*(DD+FAC*D))
      EXPIDL = EXP(CI*DD)

      RETURN
      END












