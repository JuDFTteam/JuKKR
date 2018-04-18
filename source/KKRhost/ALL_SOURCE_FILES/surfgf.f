      SUBROUTINE SURFGF(ML,M0,MR,X,ITERMAX,ERRMAX,ICHCK)
c============================================================
c
c solve surface green's function: f(x)=ml*(m0-x)**(-1)*mr
c method: decimation technique
c NEW VERSION (speeded up) by V.Bellini (march,1999)
c
c input:  ml,m0,mr - complex rectangular matrices
c output: x        - result, matrix of same type as before
c=============================================================
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAXD+1)**2)
      INTEGER NDIM
      PARAMETER (NDIM=LMMAXD*NPRINCD)
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER (CONE= (1.d0,0.d0),CZERO= (0.d0,0.d0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ERRMAX
      INTEGER ICHCK,ITERMAX
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX M0(NDIM,NDIM),ML(NDIM,NDIM),MR(NDIM,NDIM),
     +               X(NDIM,NDIM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ERR,SUM,XIM,XRE
      INTEGER I,INFO,ITER,J,N
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX AA(NDIM,NDIM),ALFA(NDIM,NDIM),BB(NDIM,NDIM),
     +               BETA(NDIM,NDIM),CC(NDIM,NDIM),CUNIT(NDIM,NDIM),
     +               EPS(NDIM,NDIM),TEMPIN(NDIM,NDIM),
     +               TEMPOUT(NDIM,NDIM),Y1(NDIM,NDIM),Y2(NDIM,NDIM)
      INTEGER IPVT(NDIM)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,ZAXPY,ZCOPY,ZGEMM,ZGETRF,ZGETRS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,DIMAG,DSQRT
C     ..

      CALL CINIT(NDIM*NDIM,CUNIT)
      DO N = 1,NDIM
        CUNIT(N,N) = CONE
      END DO

      CALL ZCOPY(NDIM*NDIM,M0,1,EPS,1)
      CALL ZCOPY(NDIM*NDIM,ML,1,ALFA,1)
      CALL ZCOPY(NDIM*NDIM,MR,1,BETA,1)
      CALL ZCOPY(NDIM*NDIM,M0,1,X,1)

      ITER = 1
   10 CONTINUE

      CALL ZCOPY(NDIM*NDIM,EPS,1,Y1,1)
      CALL ZCOPY(NDIM*NDIM,Y1,1,TEMPIN,1)
      CALL ZGETRF(NDIM,NDIM,TEMPIN,NDIM,IPVT,INFO)

c     aa = eps^-1 * alfa
      CALL ZCOPY(NDIM*NDIM,ALFA,1,TEMPOUT,1)
      CALL ZGETRS('N',NDIM,NDIM,TEMPIN,NDIM,IPVT,TEMPOUT,NDIM,INFO)
      CALL ZCOPY(NDIM*NDIM,TEMPOUT,1,AA,1)

c     bb = eps^-1 * beta

      CALL ZCOPY(NDIM*NDIM,BETA,1,TEMPOUT,1)
      CALL ZGETRS('N',NDIM,NDIM,TEMPIN,NDIM,IPVT,TEMPOUT,NDIM,INFO)
      CALL ZCOPY(NDIM*NDIM,TEMPOUT,1,BB,1)

c     alfa_new = alfa * aa

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,ALFA,NDIM,AA,NDIM,CZERO,Y1,
     +           NDIM)

c     beta_new = beta * bb

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,BETA,NDIM,BB,NDIM,CZERO,Y2,
     +           NDIM)

c     cc = - alfa * bb

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,ALFA,NDIM,BB,NDIM,CZERO,
     +           CC,NDIM)

c     x_new = x + cc

      CALL ZAXPY(NDIM*NDIM,CONE,CC,1,X,1)

c     cc = eps + cc

      CALL ZAXPY(NDIM*NDIM,CONE,CC,1,EPS,1)

c     eps_new = cc - beta * aa

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,BETA,NDIM,AA,NDIM,CONE,
     +           EPS,NDIM)

      CALL ZCOPY(NDIM*NDIM,Y1,1,ALFA,1)
      CALL ZCOPY(NDIM*NDIM,Y2,1,BETA,1)

      SUM = 0.d0
      DO I = 1,NDIM
        DO J = 1,NDIM
          XRE = DBLE(ALFA(I,J))
          XIM = DIMAG(ALFA(I,J))
          SUM = SUM + XRE*XRE + XIM*XIM
        END DO
      END DO

      ERR = DSQRT(SUM)
      IF (ERR.LT.ERRMAX .OR. ITER.GT.ITERMAX) GO TO 20
      ITER = ITER + 1
      GO TO 10

   20 CONTINUE

      CALL ZCOPY(NDIM*NDIM,X,1,TEMPIN,1)
      CALL ZCOPY(NDIM*NDIM,CUNIT,1,TEMPOUT,1)
      CALL ZGETRF(NDIM,NDIM,TEMPIN,NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,TEMPIN,NDIM,IPVT,TEMPOUT,NDIM,INFO)
      CALL ZCOPY(NDIM*NDIM,TEMPOUT,1,X,1)

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,X,NDIM,MR,NDIM,CZERO,
     +           TEMPIN,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,ML,NDIM,TEMPIN,NDIM,CZERO,
     +           X,NDIM)

      IF (ITER.GT.ITERMAX) THEN
        WRITE (6,FMT='('' itermax too small.  iter='',i3)') ITER
        write(6,'('' Surfgf:  iter='',i4,''  error='',d14.7)') iter,err
      END IF
      IF (ICHCK.EQ.0) RETURN
c       write(6,'('' Surfgf:  iter='',i4,''  error='',d12.7)') iter,err
c      write(6,'(/'' X matrix'')')
c      write(6,*)
c
      RETURN
      END
