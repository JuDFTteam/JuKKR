c 14.10.95 ***************************************************************
      SUBROUTINE CSIMPK(CF,CFINT,IPAN,IRCUT,DRDI)
c ************************************************************************
c     this subroutine does an integration up to rcut of an
c     complex function cf with an extended 3-point-simpson :
c
c                             rcut
c                      cfint = { cf(r') dr'
c                              0
c
c     modified for functions with kinks - at each kink the
c     integration is restarted .
c
c     attention : input cf is destroyed !
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c      INTEGER IPAND
c      PARAMETER (IPAND=4)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX CFINT
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CF(*)
      DOUBLE PRECISION DRDI(*)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST,N
C     ..
C     .. External Functions ..
      DOUBLE COMPLEX CSUM
      EXTERNAL CSUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      A1 = 4.0D0/3.0D0
      A2 = 2.0D0/3.0D0
      CFINT = 0.0D0
c
      DO 20 IP = 1,IPAN
c
c---> loop over kinks
c
        IST = IRCUT(IP-1) + 1
        IEN = IRCUT(IP)
c
        DO 10 I = IST,IEN
          CF(I) = CF(I)*DRDI(I)
   10   CONTINUE
c
        IF (MOD(IEN-IST,2).EQ.0) THEN
          CFINT = CFINT + (CF(IST)-CF(IEN))/3.0D0
          IST = IST + 1
          N = (IEN-IST+1)/2

        ELSE
c---> four point lagrange integration for the first step
          CFINT = CFINT + (9.0D0*CF(IST)+19.0D0*CF(IST+1)-
     +            5.0D0*CF(IST+2)+CF(IST+3))/24.0D0 +
     +            (CF(IST+1)-CF(IEN))/3.0D0
          IST = IST + 2
          N = (IEN-IST+1)/2
        END IF
c
c---> calculate with an extended 3-point-simpson
c
        CFINT = CFINT + A1*CSUM(N,CF(IST),2) + A2*CSUM(N,CF(IST+1),2)
   20 CONTINUE
c
      END
c 19.10.95 ***************************************************************
      COMPLEX*16 FUNCTION CSUM(N,V,IV)
c ************************************************************************
c        sum up the first N elements of the double complex
c        array V(*) with a stepwidth of IV
c ------------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IV,N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX VSUM
      INTEGER I,IBOT,ITOP
C     ..
      IF (IV.GE.0) THEN
        IBOT = 1
        ITOP = 1 + (N-1)*IV

      ELSE
        IBOT = 1 - (N-1)*IV
        ITOP = 1
      END IF

      VSUM = (0D0,0D0)
      DO 10 I = IBOT,ITOP,IV
        VSUM = VSUM + V(I)
   10 CONTINUE
      CSUM = VSUM
      RETURN
      END
c 14.10.95 ***************************************************************
      SUBROUTINE SIMP3(F,FINT,ISTART,IEND,DRDI)
c ************************************************************************
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
c 14.10.95 ***************************************************************
      SUBROUTINE SIMPK(F,FINT,IPAN,IRCUT,DRDI)
c ************************************************************************
c     this subroutine does an integration up to rcut of an
c     real function f with an extended 3-point-simpson :
c
c                            rcut
c                      fint = { f(r') dr'
c                             0
c
c     modified for functions with kinks - at each kink the
c     integration is restarted .
c
c     attention : input f is destroyed !
c
c-----------------------------------------------------------------------

C     .. Parameters ..
      include 'inc.p'
c     INTEGER IPAND
c     PARAMETER (IPAND=4)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION FINT
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(*),F(*)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST,N
C     ..
C     .. External Functions ..
      DOUBLE PRECISION SSUM
      EXTERNAL SSUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      A1 = 4.0D0/3.0D0
      A2 = 2.0D0/3.0D0
      FINT = 0.0D0
c
      DO 20 IP = 1,IPAN
c
c---> loop over kinks
c
        IST = IRCUT(IP-1) + 1
        IEN = IRCUT(IP)
c
        DO 10 I = IST,IEN
          F(I) = F(I)*DRDI(I)
   10   CONTINUE
        IF (MOD(IEN-IST,2).EQ.0) THEN
          FINT = FINT + (F(IST)-F(IEN))/3.0D0
          IST = IST + 1
          N = (IEN-IST+1)/2

        ELSE
c---> four point lagrange integration for the first step
          FINT = FINT + (9.0D0*F(IST)+19.0D0*F(IST+1)-5.0D0*F(IST+2)+
     +           F(IST+3))/24.0D0 + (F(IST+1)-F(IEN))/3.0D0
          IST = IST + 2
          N = (IEN-IST+1)/2
        END IF
c
c---> calculate with an extended 3-point-simpson
c
        FINT = FINT + A1*SSUM(N,F(IST),2) + A2*SSUM(N,F(IST+1),2)
   20 CONTINUE
c
      END
c 14.10.95 ***************************************************************
      SUBROUTINE SINWK(F,FINT,IPAN,IRCUT)
c ************************************************************************
c    this subroutine does an inwards integration of a function
c    with kinks
c
c
c                             rc
c                   fint(r) = s f(r') dr'
c                             r
c
c    at each kink the integration is restarted
c    the starting value for this integration is determined by
c    a 4 point lagrangian integration  , coefficients given by
c    m. abramowitz and i.a. stegun, handbook of mathematical functions,
c    nbs applied mathematics series 55 (1968)
c
c    the weights drdi have to be multiplied before calling this
c    subroutine .
c
c                                     b. drittler oct. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c     INTEGER IPAND
c     PARAMETER (IPAND=4)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION F(*),FINT(*)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST
C     ..
      A1 = 1.0D0/3.0D0
      A2 = 4.0D0/3.0D0
c
c---> loop over kinks
c
      DO 20 IP = IPAN,1,-1
        IST = IRCUT(IP)
        IEN = IRCUT(IP-1) + 1
c
        IF (IP.EQ.IPAN) THEN
          FINT(IST) = 0.0D0
c---> integrate fint(ist-1) with a 4 point lagrangian
          FINT(IST-1) = (F(IST-3)-5.0D0*F(IST-2)+19.0D0*F(IST-1)+
     +                  9.0D0*F(IST))/24.0D0

        ELSE
          FINT(IST) = FINT(IST+1)
c---> integrate fint(ist-1) with a 4 point lagrangian
          FINT(IST-1) = FINT(IST+1) + (F(IST-3)-5.0D0*F(IST-2)+
     +                  19.0D0*F(IST-1)+9.0D0*F(IST))/24.0D0
        END IF
c
c---> calculate fint with an extended 3-point-simpson
c
        DO 10 I = IST - 2,IEN,-1
          FINT(I) = ((FINT(I+2)+F(I+2)*A1)+F(I+1)*A2) + F(I)*A1
   10   CONTINUE
   20 CONTINUE
c
      END
c 14.10.95 *************************************************************
      SUBROUTINE SOUTK(F,FINT,IPAN,IRCUT)
c ************************************************************************
c    this subroutine does an outwards integration of a function
c    with kinks
c
c
c                             r
c                   fint(r) = s f(r') dr'
c                             0
c
c    at each kink the integration is restarted
c    the starting value for this integration is determined by
c    a 4 point lagrangian integration  , coefficients given by
c    m. abramowitz and i.a. stegun, handbook of mathematical functions,
c    nbs applied mathematics series 55 (1968)
c
c    the weights drdi have to be multiplied before calling this
c    subroutine .
c
c                                     b. drittler oct. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c      INTEGER IPAND
c      PARAMETER (IPAND=4)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION F(*),FINT(*)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST
C     ..
      A1 = 1.0D0/3.0D0
      A2 = 4.0D0/3.0D0
c
c---> loop over kinks
c
      DO 20 IP = 1,IPAN
        IEN = IRCUT(IP)
        IST = IRCUT(IP-1) + 1
c
        IF (IP.EQ.1) THEN
          FINT(IST) = 0.0D0
c---> integrate fint(ist+1) with a 4 point lagrangian
          FINT(IST+1) = (F(IST+3)-5.0D0*F(IST+2)+19.0D0*F(IST+1)+
     +                  9.0D0*F(IST))/24.0D0

        ELSE
          FINT(IST) = FINT(IST-1)
c---> integrate fint(ist+1) with a 4 point lagrangian
          FINT(IST+1) = FINT(IST-1) + (F(IST+3)-5.0D0*F(IST+2)+
     +                  19.0D0*F(IST+1)+9.0D0*F(IST))/24.0D0
        END IF

c
c---> calculate fint with an extended 3-point-simpson
c
        DO 10 I = IST + 2,IEN
          FINT(I) = ((FINT(I-2)+F(I-2)*A1)+F(I-1)*A2) + F(I)*A1
   10   CONTINUE
   20 CONTINUE
c
      END
c ************************************************************************
      DOUBLE PRECISION FUNCTION SSUM(N,V,IV)
c ************************************************************************
c        sum up the first N elements of the double precision 
c        array V(*) with a stepwidth of IV
c ------------------------------------------------------------------------
c     .. Scalar Arguments ..
      INTEGER IV,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VSUM
      INTEGER I,IBOT,ITOP
C     ..
      IF (IV.GE.0) THEN
        IBOT = 1
        ITOP = 1 + (N-1)*IV

      ELSE
        IBOT = 1 - (N-1)*IV
        ITOP = 1
      END IF

      VSUM = 0.0D0
      DO 10 I = IBOT,ITOP,IV
        VSUM = VSUM + V(I)
   10 CONTINUE
      SSUM = VSUM
      END
c ************************************************************************
c EOF 'integ.f'
c ************************************************************************
