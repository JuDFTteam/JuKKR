C **********************************************************************
      DOUBLE PRECISION FUNCTION DCLOCK()
C **********************************************************************
C     .. External Functions ..
CIBM  INTEGER MCLOCK
CIBM  EXTERNAL MCLOCK
      REAL ETIME,TARRY(2)
      EXTERNAL ETIME
CT3E  REAL ETIME,TARRY(2)
CT3E  EXTERNAL ETIME
CT90  REAL*8 SECOND
CT90  EXTERNAL SECOND
c
c     Next is for no timing
       DCLOCK = 0.
c
C     ..
c     next is for IBM - AIX
CIBM  DCLOCK = MCLOCK()/100.
c
c     next is for Cray UNICOS
CT90  DCLOCK = SECOND()
c
c     next is for DEC Alphas
      DCLOCK = DBLE(ETIME(TARRY))
CT3E  DCLOCK = DBLE(ETIME(TARRY))
      RETURN
      END
