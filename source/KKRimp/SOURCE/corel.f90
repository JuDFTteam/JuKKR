SUBROUTINE COREL(NSRA,IPR,IP,RHOC,V,ECORE,LCORE,NCORE,DRDI,Z,QC,A,B,IS,NSPIN,NR,RMAX,IRMD)
  !-----------------------------------------------------------------------
  !     subroutine for core states
  !-----------------------------------------------------------------------
  !     lmxc = lmaxcore = (0,1,2,...), .e.g, argon core : lmxc = 1
  !                                        krypton core : lmxc = 2
  !     kfg = configuration of core, e.g., argon core: 3300=3s,3p,0d
  !                                      krypton core: 4430=4s,4p,3d
  !                                      xenon core: 5540=5s,5p,4d
  !-----------------------------------------------------------------------
  !     .. Parameters ..
  INTEGER NITMAX,IRNUMX
  PARAMETER (NITMAX=40,IRNUMX=10)
  DOUBLE PRECISION ZERO
  PARAMETER (ZERO=0.0D0)
  !     .. Scalar Arguments ..
  DOUBLE PRECISION A,B,QC,RMAX,Z
  INTEGER IP,IPR,IRMD,IS,NCORE,NR,NSPIN,NSRA
  !     .. Array Arguments ..
  DOUBLE PRECISION DRDI(*),ECORE(*),RHOC(*),V(*)
  INTEGER LCORE(*)
  !     .. Local Scalars ..
  DOUBLE PRECISION E,E1,E2,EDIFF,EI,SLOPE,SUM,TOL,VALUE,WGT
  INTEGER IC,IN,INUC,IR,L,LMP1,LMXC,LP1,NC,NMAX,NN,NRE
  LOGICAL VLNC
  !     .. Local Arrays ..
  DOUBLE PRECISION F(IRMD),G(IRMD),RHO(IRMD)
  INTEGER KFG(4)
  CHARACTER (len=4) SPN(2),TEXT(5)
  !     .. External Subroutines ..
  EXTERNAL INTCOR,SIMP3
  !     .. Intrinsic Functions ..
  INTRINSIC DBLE,REAL
  !     .. Save statement ..
  SAVE SPN,TEXT
  !     .. Data statements ..
  DATA SPN,TEXT/'down','up  ','s   ','p   ','d   ','f   ','g   '/

  VLNC = .false.
  VALUE = 1.D-8
  SLOPE = -1.D-8
  E2 = 50.0D0

  DO 10 IC = 1,4
    KFG(IC) = 0
  10 CONTINUE
  DO 20 IC = 1,NCORE
    IF (LCORE(IC).EQ.0) KFG(1) = KFG(1) + 1
    IF (LCORE(IC).EQ.1) KFG(2) = KFG(2) + 1
    IF (LCORE(IC).EQ.2) KFG(3) = KFG(3) + 1
    IF (LCORE(IC).EQ.3) KFG(4) = KFG(4) + 1
  20 CONTINUE
  IF (KFG(2).NE.0) KFG(2) = KFG(2) + 1
  IF (KFG(3).NE.0) KFG(3) = KFG(3) + 2
  IF (KFG(4).NE.0) KFG(4) = KFG(4) + 3
  LMXC = 0
  IF (KFG(2).NE.0) LMXC = 1
  IF (KFG(3).NE.0) LMXC = 2
  IF (KFG(4).NE.0) LMXC = 3

  TOL = 1.0D-12* (Z*Z+1.D0)
  LMP1 = LMXC + 1
  NC = 0
  INUC = -IRNUMX

  DO 30 IR = 1,IRMD
    RHOC(IR) = ZERO
    RHO(IR) = ZERO
  30 CONTINUE

  DO 70 LP1 = 1,LMP1
    L = LP1 - 1
    E1 = (-5.D0- ((Z+1.D0)/DBLE(LP1))**2)*1.5D0 - 50.D0
    NMAX = KFG(LP1)
    IF (NMAX.NE.0) THEN
      DO 60 IN = LP1,NMAX
        NN = IN - LP1
        NC = NC + 1
        INUC = INUC + IRNUMX
        E = ECORE(NC)
        EI = ECORE(NC)
        IF (IPR.NE.0) WRITE (6,FMT=9000) IN,TEXT(LP1),NN,SPN(IS), IP,E
        CALL INTCOR(E1,E2,RHO,G,F,V,VALUE,SLOPE,L,NN,E,SUM,NRE, VLNC,A,B,Z,RMAX,NR,TOL,IRMD,IPR,NITMAX,NSRA)
        EDIFF = E - EI
        ECORE(NC) = E
        WGT = REAL(L+L+1)/SUM*2.D0/REAL(NSPIN)
        IF (IPR.NE.0) WRITE (6,FMT=9010) EI,EDIFF,E
  40       CONTINUE

  !---> sum up contributions to total core charge
        DO 50 IR = 2,NRE
          RHOC(IR) = RHOC(IR) + RHO(IR)*WGT
          RHO(IR) = ZERO
  50       CONTINUE
  60     CONTINUE
    END IF

  70 CONTINUE
  IF (NC*IRNUMX.GT.150 .OR. IRNUMX.GT.10) STOP 'corel'

  !---> integrate core density to get core charge
  CALL SIMP3(RHOC,QC,1,NR,DRDI)

  9000 FORMAT (1x,90 ('*'),/,'  n = ',i1,'  l = ',a4,'   nnode = ',i1, '  spin=',a4,i5,'th cell','    einput = ',1p,d16.8)
  9010 FORMAT (1x,'  einput =',1p,d16.8,'   eout - ein =',1p,d16.8, '   eoutput = ',1p,d16.8)

END SUBROUTINE COREL
