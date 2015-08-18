c 13.10.95 ***************************************************************
      SUBROUTINE VINTRAS(CMOM,CMINST,LMAX,NSPIN,NSTART,NEND,RHO2NS,V,R,
     +                   DRDI,IRWS,IRCUT,IPAN,KSHAPE,NTCELL,ILM,IFUNM,
     +                   IMAXSH,GSH,THETAS,LMSP)
c ************************************************************************
c     calculate the electron-intracell-potentials and the charge-
c     moments of given charge densities . ( for each spin-direc-
c     tion the potential is the same in the polarized case . )
c     initialize the potential v with the electron-intracell-potentials
c     the intracell-potential is expanded into spherical harmonics .
c     the lm-term of the intracell-potential of the representive atom i
c     is given by
c                    8pi        r      r'** l
c      v(r,lm,i) =  ----- *  (  s dr' --------   rho2ns(r',lm,i,1)
c                   2*l+1       0     r **(l+1)
c
c                                 rcut    r ** l
c                               +  s dr' ---------   rho2ns(r',lm,i,1) )
c                                  r     r' **(l+1)
c
c     the lm contribution of the charge moment of the representive
c     atom i is given by
c
c                             rcut
c              cmom(lm,i) =    s dr' r'** l rho2ns(r',lm,i,1)
c                              0
c
c             (see notes by b.drittler and u.klemradt)
c
c              rcut is muffin tin or wigner seitz sphere radius,
c              depending on kshape turned on or off
c
c     attention : rho2ns(...,1) is the real charge density times r**2
c                 developed into spherical harmonics . (see deck rholm)
c
c                               b.drittler   may 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
C
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,LMAX,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CMINST(LMPOTD,*),CMOM(LMPOTD,*),DRDI(IRMD,*),
     +                 GSH(*),R(IRMD,*),RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 THETAS(IRID,NFUND,*),V(IRMD,LMPOTD,*)
      INTEGER IFUNM(NATYPD,*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
     +        IRCUT(0:IPAND,*),IRWS(*),LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FAC,PI,RL
      INTEGER I,IATYP,ICELL,IEND,IFUN,IPOT,IRC1,IRS1,ISTART,J,L,LM,LM2,
     +        LM3,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION V1(IRMD),V2(IRMD),VINT1(IRMD),VINT2(IRMD)
      INTEGER IRCUTM(0:IPAND)
C     ..
C     .. External Subroutines ..
      EXTERNAL SINWK,SOUTK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,REAL
C     ..
      PI = 4.D0*ATAN(1.D0)
c
      DO 110 IATYP = NSTART,NEND
        IF (KSHAPE.NE.0) THEN
          IRS1 = IRCUT(1,IATYP)
          IRC1 = IRCUT(IPAN(IATYP),IATYP)
          ICELL = NTCELL(IATYP)
          DO 10 I = 0,IPAN(IATYP)
            IRCUTM(I) = IRCUT(I,IATYP)
   10     CONTINUE

        ELSE

          IRS1 = IRWS(IATYP)
          IRC1 = IRS1
          IRCUTM(0) = IRCUT(0,IATYP)
          IRCUTM(1) = IRC1
        END IF
c---> determine the right potential numbers
        IPOT = NSPIN*IATYP

        DO 100 L = 0,LMAX
          FAC = 8.0D0*PI/REAL(2*L+1)
          DO 90 M = -L,L
            LM = L*L + L + M + 1
c
c---> set up of the integrands v1 and v2
c
            V1(1) = 0.0D0
            V2(1) = 0.0D0
            DO 20 I = 2,IRS1
              RL = R(I,IATYP)**L
              V1(I) = RHO2NS(I,LM,IATYP,1)*RL*DRDI(I,IATYP)
              V2(I) = RHO2NS(I,LM,IATYP,1)/R(I,IATYP)/RL*DRDI(I,IATYP)
   20       CONTINUE
c
c---> convolute charge density of interstial with shape function
c        if kshape.gt.0
c
            IF (KSHAPE.NE.0) THEN
              DO 30 I = IRS1 + 1,IRC1
                V1(I) = 0.0D0
   30         CONTINUE
              ISTART = IMAXSH(LM-1) + 1
              IEND = IMAXSH(LM)
              DO 50 J = ISTART,IEND
                LM2 = ILM(J,2)
                LM3 = ILM(J,3)
                IF (LMSP(ICELL,LM3).GT.0) THEN
                  IFUN = IFUNM(ICELL,LM3)
                  DO 40 I = IRS1 + 1,IRC1
                    V1(I) = V1(I) + GSH(J)*RHO2NS(I,LM2,IATYP,1)*
     +                      THETAS(I-IRS1,IFUN,ICELL)
   40             CONTINUE
                END IF
   50         CONTINUE

              DO 60 I = IRS1 + 1,IRC1
                RL = R(I,IATYP)**L
                V2(I) = V1(I)/R(I,IATYP)/RL*DRDI(I,IATYP)
                V1(I) = V1(I)*RL*DRDI(I,IATYP)
   60         CONTINUE
            END IF
c
c---> now integrate v1 and v2
c
            CALL SOUTK(V1,VINT1,IPAN(IATYP),IRCUTM)
            CALL SINWK(V2,VINT2,IPAN(IATYP),IRCUTM)
c
c---> gather all parts
c
            IF (LM.EQ.1) THEN
              V(1,LM,IPOT) = FAC*VINT2(1)

            ELSE

              V(1,LM,IPOT) = 0.0D0
            END IF

            DO 70 I = 2,IRC1
              RL = R(I,IATYP)**L
              V(I,LM,IPOT) = FAC* (VINT1(I)/R(I,IATYP)/RL+VINT2(I)*RL)
   70       CONTINUE
c
c---> store charge moment - in case of kshape.gt.0 this is the moment
c      of the charge in the muffin tin sphere
c
            CMOM(LM,IATYP) = VINT1(IRS1)
c
c---> store charge moment of interstial in case of kshape.gt.0
c
            IF (KSHAPE.NE.0) CMINST(LM,IATYP) = VINT1(IRC1) -
     +          VINT1(IRS1)
c
            IF (NSPIN.EQ.2) THEN
              DO 80 I = 1,IRC1
                V(I,LM,IPOT-1) = V(I,LM,IPOT)
   80         CONTINUE
            END IF

   90     CONTINUE

  100   CONTINUE

  110 CONTINUE

      RETURN

      END
