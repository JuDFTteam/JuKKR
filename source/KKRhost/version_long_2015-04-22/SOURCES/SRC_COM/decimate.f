      SUBROUTINE DECIMATE(GLLKE,NAEZ,TINVBUP,TINVBDOWN,VACFLAG,
     +                    FACTL,NLBASIS,NRBASIS)
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
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAXD+1)**2)  ! ruess: for decimate with SOC
      INTEGER ALM,NDIM
      PARAMETER (ALM=NAEZD*LMMAXD,NDIM=NPRINCD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      INTEGER NAEZ,NLBASIS,NRBASIS
C     ..
C     .. Array Arguments ..

      DOUBLE COMPLEX GLLKE(ALM,ALM),TINVBDOWN(LMMAXD,LMMAXD,*),
     +               TINVBUP(LMMAXD,LMMAXD,*)
      DOUBLE COMPLEX FACTL(LMMAXD,LMMAXD)
      LOGICAL VACFLAG(2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ERRMAX
      INTEGER ICHCK,IHOST,II1,II2,IL1,IL2,IP1,IP1T,IP2,IP2T,ITERMAX,
     +        LDI1,LDI1T,LDI2,LDI2T,LM1,LM2,NLAYER,ICALL
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX A1(NDIM,NDIM),AN(NDIM,NDIM),B1(NDIM,NDIM),
     +               BN(NDIM,NDIM),C1(NDIM,NDIM),CN(NDIM,NDIM),
     +               X1(NDIM,NDIM),XN(NDIM,NDIM)
C     ..
      DATA ICALL /0/
C     ..
C     .. External Subroutines ..
      EXTERNAL BOFM,SURFGF
C     ..
C     .. Save statement ..
      SAVE ICALL,NLAYER,ITERMAX,ERRMAX,ICHCK

C     ..
c ----------------------------------------------------------------------
C     .. External Functions ..
      LOGICAL OPT
      EXTERNAL OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
c

      ICALL = ICALL + 1
c ----------------------------------------------------------------------
      IF (ICALL.EQ.1) THEN
         NLAYER = NAEZ/NPRINCD
c     Parameters for the "decimation" technique.
         ITERMAX = 300
         ERRMAX = 1.0D-180
         ICHCK = 1
      END IF
c ----------------------------------------------------------------------
      IF ( .NOT.VACFLAG(1) ) THEN
c
C Get the matrix B1
c
        CALL BOFM(1,1,B1,NDIM,GLLKE,ALM)

c Now Subtract t-mat of left host
        DO IP1 = 1,NPRINCD
          IHOST = MOD(IP1-1,NLBASIS) + 1
          DO LM1 = 1,LMMAXD
            DO LM2 = 1,LMMAXD
              IL1 = LMMAXD* (IP1-1) + LM1
              IL2 = LMMAXD* (IP1-1) + LM2
              B1(IL1,IL2) = (B1(IL1,IL2)-TINVBUP(LM1,LM2,IHOST))
            END DO
          END DO
        END DO

        CALL BOFM(1,2,C1,NDIM,GLLKE,ALM)
        CALL BOFM(2,1,A1,NDIM,GLLKE,ALM)

c     it performs the 'space decimation' iterative procedure.
        CALL SURFGF(A1,B1,C1,X1,ITERMAX,ERRMAX,ICHCK)
c     adds to the matrix GLLKE the elements that couples the
c     interface to the two half-spaces.
        DO IP1 = 1,NPRINCD
          DO IP2 = 1,NPRINCD
            II1 = IP1
            II2 = IP2
            DO LM1 = 1,LMMAXD
              DO LM2 = 1,LMMAXD
                LDI1 = LMMAXD* (IP1-1) + LM1
                IL1 = LMMAXD* (II1-1) + LM1
                LDI2 = LMMAXD* (IP2-1) + LM2
                IL2 = LMMAXD* (II2-1) + LM2
                GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - X1(LDI1,LDI2)
              END DO
            END DO
          END DO
        END DO
                
      END IF

      IF ( .NOT.VACFLAG(2) ) THEN

c  If 'ONEBULK' is activated then it calculates the xn decimated element
c  from the x1 element: this is just in the case of equal bulks on the

        IF ( .NOT.OPT('ONEBULK ') ) THEN

c
C     Get the matrix BN
c
          CALL BOFM(NLAYER,NLAYER,BN,NDIM,GLLKE,ALM)

c Now Substract t-mat right host
c Notes : the indexing is easier like that
          DO IP1 = 1,NPRINCD
            IHOST = NRBASIS - MOD(IP1,NRBASIS)
            IHOST = MOD(IP1-1,NRBASIS) + 1
            DO LM1 = 1,LMMAXD
              DO LM2 = 1,LMMAXD
                IL1 = LMMAXD* (IP1-1) + LM1
                IL2 = LMMAXD* (IP1-1) + LM2
                BN(IL1,IL2) = (BN(IL1,IL2)-
     +                        TINVBDOWN(LM1,LM2,IHOST))
              END DO
            END DO
          END DO

          CALL BOFM(NLAYER,NLAYER-1,AN,NDIM,GLLKE,ALM)
          CALL BOFM(NLAYER-1,NLAYER,CN,NDIM,GLLKE,ALM)

c     it performs the 'space decimation' iterative procedure.
          CALL SURFGF(CN,BN,AN,XN,ITERMAX,ERRMAX,ICHCK)
c     
c
c
        ELSE
ccccccccc
          DO IP1 = 1,NPRINCD
            DO IP2 = 1,NPRINCD
              IP1T = (NPRINCD+1) - IP2
              IP2T = (NPRINCD+1) - IP1
              DO LM1 = 1,LMMAXD
                DO LM2 = 1,LMMAXD
                  LDI1 = LMMAXD* (IP1-1) + LM1
                  LDI2 = LMMAXD* (IP2-1) + LM2
                  LDI1T = LMMAXD* (IP1T-1) + LM2
                  LDI2T = LMMAXD* (IP2T-1) + LM1
                  XN(LDI1T,LDI2T) = FACTL(LM1,LM2)*X1(LDI1,LDI2)
                END DO
              END DO
            END DO
          END DO

        END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             Added on 1.02.2000
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     adds to the matrix GLLKE the elements that couples the
c     interface to the two half-spaces.
        DO IP1 = 1,NPRINCD
          DO IP2 = 1,NPRINCD
            II1 = (NLAYER-1)*NPRINCD + IP1
            II2 = (NLAYER-1)*NPRINCD + IP2
            DO LM1 = 1,LMMAXD
              DO LM2 = 1,LMMAXD
                LDI1 = LMMAXD* (IP1-1) + LM1
                IL1 = LMMAXD* (II1-1) + LM1
                LDI2 = LMMAXD* (IP2-1) + LM2
                IL2 = LMMAXD* (II2-1) + LM2
                GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - XN(LDI1,LDI2)
              END DO
            END DO
          END DO
        END DO
        
      END IF


      RETURN

      END
