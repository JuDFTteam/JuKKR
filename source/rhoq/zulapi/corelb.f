c 20.09.95 ***************************************************************
      SUBROUTINE CORELB(IPF,NITMAX,KHYP,NONSRA,IPR,IRNUMX,IP,IRM,NSPIN,
     +                  RHOC,IMT,KSHAPE,V,RWS,ECORE,DRDI,R,Z,A,B,IRWS,
     +                  LMXC,KFG,RHYPF,HYPSUM)
c ************************************************************************
c
c     this subroutine is responsible for relaxation of core-charge 
c     density,
c     it calls intcor=integration of scalarrelativistic eqation 
c     and find the appropriate eigenvalues
c     lmxc = lmaxcore = (0,1,2,...), argon core : lmxc = 1, 
c                                    krypton    : lmxc = 2,
c                                    xenon      : lmxc = 2,
c     kfg = configuration of core states f.e. argon core: 3300=3s,3p,0d
c                                           krypton core: 4430=4s,4p,3d
c                                             xenon core: 5540=5s,5p,4d
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c
c      INTEGER NATYPD,NSPIND
c      PARAMETER (NATYPD=1,NSPIND=2)
c      INTEGER IRMD
c      PARAMETER (IRMD=424)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A(*),B(*),Z(*)
      INTEGER IP,IPF,IPR,IRM,IRNUMX,KHYP,KSHAPE,NITMAX,NONSRA,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),ECORE(20,*),R(IRMD,*),RHOC(IRMD,*),
     +                 RWS(*),V(IRMD,*)
      INTEGER IMT(*),IRWS(*),KFG(4,*),LMXC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DIFF,E,E1,E2,EDIFF,EI,FOURPI,GAUSS,QC,R0,R2RHO1,
     +                 R2RHO2,RMAX,RR,SLOPE,SUM,TOL,VALUE,ZERO
      INTEGER I,IARRAY,ID,IDD,IIR,IN,INUC,INUCP1,IPOT,IR,IRIPE,IRIPST,
     +        IRM2,IRR,IS,ISY,IU,IUU,K,L,LMP1,LP1,NC,NMAX,NN,NR,NREM
      LOGICAL VLNC
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRH(4),F(IRMD,2),G(IRMD,2),HYPSUM(10,NPOTD),
     +                 RHO(IRMD,2),RHYPF(150,NPOTD),WGT(2)
      INTEGER NRE(2)
      CHARACTER*4 SPN(2),TEXT(7)
      CHARACTER*8 TEXTAT(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL BREIT,INTCOR,RCSTOP,SIMP3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,REAL
C     ..
C     .. Save statement ..
      SAVE SPN,TEXT,TEXTAT,ZERO,FOURPI,GAUSS
C     ..
C     .. Data statements ..
      DATA SPN/'down','up  '/
      DATA TEXT/'s   ','p   ','d   ','f   ','g   ','h   ','i   '/
      DATA TEXTAT/'    host','impurity'/
      DATA ZERO,FOURPI,GAUSS/0.0D0,12.56637D0,524.D0/
C     ..
c
c      write(6,*) '>>> CORELB: core charge density, IP = ',ip
      ISY = 1
      VLNC = .false.
      VALUE = 1.D-8
      SLOPE = -1.D-8
      E2 = 50.0D0

c
      ID = NSPIN* (IP-1) + 1
      IU = ID + NSPIN - 1
      IDD = 1
      IUU = NSPIN
c
c
c---> attention renormalization confined to muffin tin sphere in case
c---> of shape-corrected calculation
      IF (KSHAPE.NE.0) THEN
        NR = IMT(IP)
        RMAX = R(NR,IP)
        IRM2 = NR

      ELSE

        NR = IRWS(IP)
        RMAX = RWS(IP)
        IRM2 = IRM
      END IF

      TOL = 1.0D-12* (Z(IP)*Z(IP)+1.D0)
      LMP1 = LMXC(IP) + 1
      NC = 0
      INUC = -IRNUMX
c
      DO 10 IR = 1,IRM2
        RHOC(IR,ID) = ZERO
        RHOC(IR,IU) = ZERO
        RHO(IR,IDD) = ZERO
        RHO(IR,IUU) = ZERO
   10 CONTINUE
      DO 20 IR = 1,IRNUMX
        HYPSUM(IR,ID) = ZERO
        HYPSUM(IR,IU) = ZERO
        RHYPF(IR,ID) = ZERO
        RHYPF(IR,IU) = ZERO
   20 CONTINUE
c
      IF (IPR.NE.0) WRITE (IPF,FMT=9000) IP,TEXTAT(ISY)
c
c    begin loop over l
c
      DO 130 LP1 = 1,LMP1
        L = LP1 - 1
        E1 = (-5.D0- ((Z(IP)+1.D0)/REAL(LP1))**2)*1.5D0 - 50.D0
        NMAX = KFG(LP1,IP)
        IF (NMAX.NE.0) THEN
c
c    core states
c
          DO 120 IN = LP1,NMAX
            NN = IN - LP1
            NC = NC + 1
            INUC = INUC + IRNUMX
            DO 30 IS = 1,NSPIN
              I = NSPIN* (IP-1) + IS
              E = ECORE(NC,I)
              EI = ECORE(NC,I)
              IF (IPR.NE.0) WRITE (IPF,FMT=9010) IN,TEXT(LP1),NN,
     +            SPN(IS),IP,TEXTAT(ISY),E
              CALL INTCOR(E1,E2,RHO(1,IS),G(1,IS),F(1,IS),V(1,I),VALUE,
     +                    SLOPE,L,NN,E,SUM,NRE(IS),VLNC,
     +                    A(IP),B(IP),Z(IP),RMAX,NR,
     +                    TOL,IRM2,IPR,NITMAX,NONSRA)
              EDIFF = E - EI
              ECORE(NC,I) = E
              WGT(IS) = REAL(L+L+1)/SUM
              IF (IPR.NE.0) WRITE (IPF,FMT=9020) EI,EDIFF,E
   30       CONTINUE
c
            IF (NSPIN.EQ.1) THEN
              NREM = NRE(IUU)

            ELSE
              NREM = MAX(NRE(IUU),NRE(IDD))
            END IF
c
            IF (KHYP.NE.0) THEN
c
              INUCP1 = INUC + 1
              RHYPF(INUCP1,ID) = ZERO
              RHYPF(INUCP1,IU) = ZERO
c
              DO 40 IR = 2,IRNUMX
                IRR = IR + INUC
                RR = R(IR,IP)
                RR = FOURPI*RR*RR
                R2RHO1 = WGT(IDD)*RHO(IR,IDD)
                R2RHO2 = WGT(IUU)*RHO(IR,IUU)
                SUM = (R2RHO1+R2RHO2)/RR
c
c                    diff=gauss*(r2rho2-r2rho1)/rr
c
                CALL BREIT(DIFF,NONSRA,IR,IUU,NREM,L,FOURPI,GAUSS,
     +                     R2RHO1,R2RHO2,RR,WGT,R(1,IP),DRDI(1,IP),G,F)
                RHYPF(IRR,IU) = DIFF
                RHYPF(IRR,ID) = SUM
   40         CONTINUE
c
              IRIPST = 2
              IRIPE = 4
              DO 90 IS = 1,NSPIN
                I = NSPIN* (IP-1) + IS
                R0 = R(1,IP)
                DO 50 IR = 1,IRIPE
                  IRR = INUC + IR
                  DRH(IR) = RHYPF(IRR,I)
   50           CONTINUE
                DO 70 K = 1,2
                  DO 60 IR = IRIPST,IRIPE - K
                    IIR = IRIPST + IRIPE - IR
                    DRH(IIR) = (DRH(IIR)-DRH(IIR-1))/
     +                         (R(IIR,IP)-R(IIR-K,IP))
   60             CONTINUE
   70           CONTINUE
                IR = 1
                IRR = INUC + IR
                RHYPF(IRR,I) = DRH(IRIPE)
                DO 80 IR = IRIPST,IRIPE - 1
                  IIR = IRIPST + IRIPE - IR - 1
                  RHYPF(IRR,I) = RHYPF(IRR,I)* (R0-R(IIR,IP)) + DRH(IIR)
   80           CONTINUE
   90         CONTINUE
              IF (L.EQ.0) THEN
                DO 100 IR = 1,IRNUMX
                  IRR = IR + INUC
                  HYPSUM(IR,ID) = HYPSUM(IR,ID) + RHYPF(IRR,ID)
                  IF (NSPIN.EQ.2) HYPSUM(IR,IU) = HYPSUM(IR,IU) +
     +                RHYPF(IRR,IU)
  100           CONTINUE
              END IF

            END IF                  ! (KHYP.NE.0)
c
c---> sum up contributions to total core charge
c
            DO 110 IR = 2,NREM
              RHOC(IR,ID) = RHOC(IR,ID) + RHO(IR,IDD)*WGT(IDD)
              RHOC(IR,IU) = RHOC(IR,IU) + RHO(IR,IUU)*WGT(IUU)
              RHO(IR,IDD) = ZERO
              RHO(IR,IUU) = ZERO
 110        CONTINUE                ! IR = 2,NREM
 120      CONTINUE                  ! IN = LP1,NMAX
        END IF                      ! (NMAX.NE.0)

 130  CONTINUE                      ! LP1 = 1,LMP1

      IARRAY = NC*IRNUMX
      IF (IARRAY.GT.150 .OR. IRNUMX.GT.10) THEN
        WRITE (IPF,FMT=9030) IARRAY,IRNUMX
        CALL RCSTOP('21      ')

      ELSE


        DO 140 IS = 1,NSPIN
          IPOT = NSPIN* (IP-1) + IS
c
c---> integrate core density to get core charge
c
          CALL SIMP3(RHOC(1,IPOT),QC,1,NR,DRDI(1,IP))

          WRITE (IPF,FMT=9040) IP,Z(IP),QC
  140   CONTINUE

      END IF
c

 9000 FORMAT (1x,5 ('*'),' core-relaxation for ',i3,'th ',a8,'-cell',
     +       '  spin=',a4,'   l = ',a4,' was done ',5 ('*'))
 9010 FORMAT (1x,77 ('*'),/,'  n = ',i1,'  l = ',a4,'   nnode = ',i1,
     +       '  spin=',a4,i5,'th ',a8,'-cell','    einput = ',1p,d16.8)
 9020 FORMAT (1x,'  einput =',1p,d16.8,'   eout - ein =',1p,d16.8,
     +       '   eoutput = ',1p,d16.8)
 9030 FORMAT (1x,
     +' space of arrays rhypf (',I6,') or hypsum (',I6,
     +     ') to small stop in subroutine corel')
 9040 FORMAT (I6,'  NUCLEAR CHARGE',F12.6,9X,'CORE CHARGE =',f13.6)
      END
