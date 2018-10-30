      MODULE MOD_BESSEL1
      CONTAINS
      SUBROUTINE BESSEL1(BJ,Y,H,ARG,LMX,LMAX,LJ,LY,LH,LCALL)
c   *****************************************************************
c   *this subroutine computes the spherical bessel functions of
c   * first ,second and third kind using a  chebychev expansion
c   * given by y.l.luke ,algorithms for the computation of
c   * mathematical functions, academic press,london 1977
c   *
c   *
c   * using subroutine cnwf01
c   *
c   *
c   *description of variables
c   *
c   *arg   -input  - argument of the bessel functions
c   *
c   *lmax  -input  - max. order of the bessel functions
c   *                (limited up to 25 in that version)
c   *
c   *lj    -input  - logical : if lj is true the spherical bessel
c   *                functions of the first kind are calculated
c   *                up to lmax
c   *
c   *ly    -input  - logical : if ly is true the spherical bessel
c   *                functions of the second kind are calculated
c   *                up to lmax
c   *
c   *lh    -input  - logical : if lh is true the spherical bessel
c   *                functions of the third kind are calculated
c   *                up to lmax
c   *
c   *lcall -input  - logical : if lh is false the chebychev coefficients
c   *                are calculated - this part has to be called once
c   *
c   *
c   *bj    -output - an array containing the bessel functions of the
c   *                first kind up to lmax if lj is true . remember ,
c   *                that bj(1) contains the function of l=0 and so on.
c   *
c   *y     -output - an array containing the bessel functions of the
c   *                second kind up to lmax if ly is true . remember ,
c   *                that y(1) contains the function of l=0 and so on.
c   *
c   *h     -output - an array containing the bessel functions of the
c   *                third kind up to lmax if lh is true . remember ,
c   *                that h(1) contains the function of l=0 and so on.
c   *
c   *        attention : contrary to abramowitz and stegun the bessel
c   *                    functions of third kind ( hankel functions)
c   *                    are definied as:
c   *                            h(l) = y(l) - i * bj(l)
c   *
c   *
c   *all other variables are for internal use
c**********************************************************************
C     .. Parameters ..
      USE MOD_CNWF011
      IMPLICIT NONE
      INTEGER NDIM
      PARAMETER (NDIM=24)
      INTEGER NDIMP1
      PARAMETER (NDIMP1=NDIM+1)
      INTEGER NDIMP2
      PARAMETER (NDIMP2=NDIM+2)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 ARG
      INTEGER LMAX,LMX
      LOGICAL LCALL,LH,LJ,LY
C     ..
C     .. Array Arguments ..
      COMPLEX*16 BJ(LMX+1),H(LMX+1),Y(LMX+1)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CI,CONE,CTWO,CZERO,RES,T0,T1,T2,TARG,TT1
      REAL*8 CPJ,CPY,FJ,FY,ONE,SUM,W,W2
      INTEGER L,LMSAVE,N,NCNW,NMAX
C     ..
C     .. Local Arrays ..
      COMPLEX*16 TN(20),ZN(NDIMP2)
      REAL*8 CNWJ(20,NDIMP1),CNWY(20,NDIMP1)
C     ..
C     .. External Subroutines ..
      EXTERNAL CNWF01,RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO,CONE,CTWO/ (0.D0,0.D0), (1.D0,0.D0), (2.D0,0.D0)/
      DATA CI/ (0.D0,1.D0)/
      DATA W2,ONE,NMAX/10.D0,1.D0,20/
C     ..

      IF (.NOT.LCALL) THEN
        IF (LMAX.GT.25) THEN
          WRITE (6,FMT=9000)
          STOP '27      '

        ELSE
          LMSAVE = LMAX + 1
          NCNW = NMAX - 2
          W = -W2*W2*.25D0
          FJ = ONE
          FY = -ONE
          DO 20 L = 1,LMSAVE
            CPJ = 0.5D0 + L
            CPY = 1.5D0 - L
            CALL CNWF011(CPJ,W,NCNW,CNWJ(1,L),SUM)
            CALL CNWF011(CPY,W,NCNW,CNWY(1,L),SUM)
            DO 10 N = 1,NMAX
              CNWJ(N,L) = CNWJ(N,L)/FJ
              CNWY(N,L) = CNWY(N,L)*FY
   10       CONTINUE
            FJ = FJ* (L+L+1)
            FY = FY* (L+L-1)
   20     CONTINUE
          LCALL = .true.
        END IF

      END IF

      IF (LMAX.GT. (LMSAVE-1) .OR. ABS(ARG).GT.W2) THEN
        WRITE (6,FMT=9010) LMAX,ABS(ARG)
        STOP '28      '

      ELSE
c
c-----calculate arg**n and tn*(-arg**2/4)
c
        ZN(1) = CONE
        DO 30 L = 2,LMSAVE
          ZN(L) = ZN(L-1)*ARG
   30   CONTINUE
        ZN(LMSAVE+1) = ZN(LMSAVE)*ARG
        T0 = CONE
        TARG = -ZN(3)*0.25D0/W
        T1 = CTWO*TARG - CONE
        TN(1) = T0
        TN(2) = T1
        TT1 = T1 + T1
        DO 40 N = 3,NMAX
          T2 = TT1*T1 - T0
          TN(N) = T2
          T0 = T1
          T1 = T2
   40   CONTINUE
        IF (LJ .OR. LH) THEN
          DO 60 L = 1,LMSAVE
            RES = CZERO
            DO 50 N = 1,NMAX
              RES = RES + TN(N)*CNWJ(N,L)
   50       CONTINUE
            BJ(L) = RES*ZN(L)
   60     CONTINUE
          IF (.NOT. (LY.OR.LH)) RETURN
        END IF

        IF (LY .OR. LH) THEN
          DO 80 L = 1,LMSAVE
            RES = CZERO
            DO 70 N = 1,NMAX
              RES = RES + TN(N)*CNWY(N,L)
   70       CONTINUE
            Y(L) = RES/ZN(L+1)
   80     CONTINUE
          IF (LH) THEN
            DO 90 L = 1,LMSAVE
c
c---> changed !   h(l) = bj(l) + ci*y(l)
c
              H(L) = Y(L) - CI*BJ(L)
   90       CONTINUE
          END IF

        ELSE
          WRITE (6,FMT=9020)
        END IF

      END IF


 9000 FORMAT (' lmax is too high, error stop in bessel')
 9010 FORMAT (' lmax is higher than previously given',
     +       '   *********   or argument too high, error stop in bessel'
     +       ,/,13x,' lmax : ',i5,' abs(arg) : ',f12.6)
 9020 FORMAT (' *******  error warning from bessel: no output ',
     +       ' required   ********',/,/)
      END SUBROUTINE
      END MODULE MOD_BESSEL1
