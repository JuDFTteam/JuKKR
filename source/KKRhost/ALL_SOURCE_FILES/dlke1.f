c 01.06.99 *************************************************************
      SUBROUTINE DLKE1(GLLKE,ALAT,NACLS,NACLSMAX,RR,EZOA,
     +                 ATOM,BZKP,IC,GINP,RCLS)
     
      use mod_types, only: t_inc
      implicit none
c **********************************************************************
c
c     Fourier transformation of the cluster Greens function GINP
c
c ----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
      INTEGER ALMGF0
      PARAMETER (ALMGF0=LMGF0D*NAEZD)
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IC,NACLSMAX
C     ..
C     .. Array Arguments ..
      INTEGER ATOM(*),
     +        EZOA(*),
     +        NACLS(*)
      DOUBLE COMPLEX GLLKE(ALMGF0,*),
     +               GINP(LMGF0D*NACLSMAX,*)
      DOUBLE PRECISION BZKP(*),
     +                 RR(3,0:NRD),
     +                 RCLS(3,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONVPU,TPI
      INTEGER AM,I,II,IM,LM2,M,N,N1,NL
      DOUBLE COMPLEX  EIKR,TT
      LOGICAL OPT,TEST
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX ARG(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,TEST,OPT,ZAXPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,EXP
C     ..
C     .. Save statement ..
      SAVE
C     ..
C
      II = 3
      IF (OPT('COMPLEX ')) II = 6
      IF (TEST('BZKP    ').and.(t_inc%i_write>0)) 
     &      write(1337,FMT='(6f12.6)') (bzkp(i),i=1,ii)
c
      TPI = 8.0D0*ATAN(1.0D0)                 ! = 2*PI
      CONVPU = ALAT/TPI

      CALL CINIT(LMGF0D*NAEZD*LMGF0D,GLLKE)

      DO 90 M = 1,NACLS(IC)

c
c --->  for option 'WIRE': avoid artifical couplings in the structural
c       Greens Function in in-plane-direction (perp. to c-axis)
c
        IF (ATOM(M).LT.0) GOTO 90
c
        IF (OPT('ONEBULK ')) THEN     ! added 1.02.2000
c                                       corrected on 25.02.2000
c     if the phase factor exp(ik(r-r')) is included      ~
c     in the G...so if we resolve the dyson part for the G
c     and not for the G (see Peter Lang Ph.D thesis)
c     
c     Here we do   --                           nn'
c                  \                            ii'          ii'
c                  /  exp( ik(X  -X + R  -R  ))G   (E)  =   G   (k,E)
c                  --          n'  n   i'  i    LL'          LL'
c                  n'
c                   
c In this case rcls is always (by constraction symmetric around each
c atom this means that a minus sign will not affect the result of the
c summation

           ARG(1) = -CI*TPI*RCLS(1,M)
           ARG(2) = -CI*TPI*RCLS(2,M)
           ARG(3) = -CI*TPI*RCLS(3,M)
        ELSE
c     
c     Here we do   --                  nn'
c                  \                   ii'          ii'
c                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
c                  --          n'  n   LL'          LL'
c                  n'
c  Be carefull a minus sign must be included here. RR is not
c  symmetric around each atom. The minus comes from the fact that
c  the repulsive potential GF is calculated for 0n and not n0!                   
c  and that is why we nead a minus sign extra!
c  

           ARG(1) = -CI*TPI*RR(1,EZOA(M))
           ARG(2) = -CI*TPI*RR(2,EZOA(M))
           ARG(3) = -CI*TPI*RR(3,EZOA(M))
        END IF
c
        TT = BZKP(1)*ARG(1)+BZKP(2)*ARG(2)+BZKP(3)*ARG(3)
c
        IF (OPT('COMPLEX ')) THEN
          TT = TT + CI*(BZKP(4)*ARG(1)+BZKP(5)*ARG(2)+BZKP(6)*ARG(3))
        END IF
c
c        write(6,*) BZKP(1),ARG(1),BZKP(2),ARG(2),BZKP(3),ARG(3)
c        write(6,*) BZKP(4),BZKP(5),BZKP(6)
c        write(6,*) 'm,atom(m),tt',m,atom(m),tt
c
!         EIKR = (1.0d0, 0.0d0)    ! convert to p.u.
        EIKR = EXP(TT) * CONVPU    ! convert to p.u.

c        write(6,*) 'eikr',eikr

        IM = 1 + (M-1)      *LMGF0D
        AM = 1 + (ATOM(M)-1)*LMGF0D
        DO 80 LM2 = 1,LMGF0D
          CALL ZAXPY(LMGF0D,EIKR,GINP(IM,LM2),1,GLLKE(AM,LM2),1)  
   80   CONTINUE

 90   CONTINUE                      ! M = 1,NACLS(IC)

      RETURN
 9000 format(3f12.4)
 9010 format(2f18.10)
      END
