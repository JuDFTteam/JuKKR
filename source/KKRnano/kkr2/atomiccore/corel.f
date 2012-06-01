      SUBROUTINE COREL(NSRA,IPR,IP,RHOC,V,ECORE,LCORE,NCORE,DRDI,Z,QC,
     +                   A,B,IS,NSPIN,NR,RMAX,IRMD,EBOT)
c-----------------------------------------------------------------------
c     subroutine for core states
c-----------------------------------------------------------------------
c     lmxc = lmaxcore = (0,1,2,...), .e.g, argon core : lmxc = 1
c                                        krypton core : lmxc = 2
c     kfg = configuration of core, e.g., argon core: 3300=3s,3p,0d
c                                      krypton core: 4430=4s,4p,3d
c                                      xenon core: 5540=5s,5p,4d
c-----------------------------------------------------------------------
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NITMAX,IRNUMX
      PARAMETER (NITMAX=40,IRNUMX=10)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,QC,RMAX,Z,EBOT
      INTEGER IP,IPR,IRMD,IS,NCORE,NR,NSPIN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(*),ECORE(*),RHOC(*),V(*)
      INTEGER LCORE(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION E,E1,E2,EDIFF,EI,SLOPE,SUM,TOL,VALUE,WGT
      INTEGER IC,IN,INUC,IR,L,LMP1,LMXC,LP1,NC,NMAX,NN,NRE
      LOGICAL VLNC
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION F(IRMD),G(IRMD),RHO(IRMD)
      INTEGER KFG(4)
      CHARACTER*4 SPN(2),TEXT(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL INTCOR,SIMP3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,REAL
C     ..
C     .. Save statement ..
      SAVE SPN,TEXT
C     ..
C     .. Data statements ..
      DATA SPN,TEXT/'down','up  ','s   ','p   ','d   ','f   ','g   '/
C     ..
      VLNC = .false.
      VALUE = 1.D-8
      SLOPE = -1.D-8
      E2 = 50.0D0
c
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
c
      TOL = 1.0D-12* (Z*Z+1.D0)
      LMP1 = LMXC + 1
      NC = 0
      INUC = -IRNUMX
c
      DO 30 IR = 1,IRMD
        RHOC(IR) = ZERO
        RHO(IR) = ZERO
   30 CONTINUE
c
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
            IF (IPR.NE.0) WRITE (6,FMT=9000) IN,TEXT(LP1),NN,SPN(IS),
     +          IP,E

            CALL INTCOR(E1,E2,RHO,G,F,V,VALUE,SLOPE,L,NN,E,SUM,NRE,
     +                    VLNC,A,B,Z,RMAX,NR,TOL,IRMD,IPR,NITMAX,NSRA)

        IF (E.GT.EBOT) THEN
        WRITE(6,'(''Error for L='',I1)') L
        WRITE(6,*) 'E,EBOT',E,EBOT
        WRITE(6,*) 'The program found a core state above the bottom'
        WRITE(6,*) 'of the valence-band energy contour.'
        WRITE(6,*) 'The results are very probably wrong.'
        WRITE(6,*) 'The number of core states in the input potential'
        WRITE(6,*) 'should perhaps be decreased.'
        WRITE(6,*) 'The program stops in subroutine corel.f.'
        STOP 'Error 1 in corel.f'
        END IF
            EDIFF = E - EI
            ECORE(NC) = E
            WGT = REAL(L+L+1)/SUM*2.D0/REAL(NSPIN)
            IF (IPR.NE.0) WRITE (6,FMT=9010) EI,EDIFF,E
   40       CONTINUE
c
c---> sum up contributions to total core charge
c
            DO 50 IR = 2,NRE
              RHOC(IR) = RHOC(IR) + RHO(IR)*WGT
              RHO(IR) = ZERO
   50       CONTINUE
   60     CONTINUE
        END IF
        IF (NMAX.NE.0) THEN
            IN = NMAX+1
            NN = IN - LP1
            E = ECORE(NC)/10.
            EI = ECORE(NC)/10.
            IF (IPR.NE.0) WRITE (6,FMT=9000) IN,TEXT(LP1),NN,SPN(IS),
     +          IP,E

            CALL INTCOR(E1,E2,RHO,G,F,V,VALUE,SLOPE,L,NN,E,SUM,NRE,
     +                    VLNC,A,B,Z,RMAX,NR,TOL,IRMD,IPR,NITMAX,NSRA)

        IF (E.LT.EBOT) THEN
        WRITE(6,'(''Error for L='',I1)') L
        WRITE(6,*) 'E,EBOT',E,EBOT
        WRITE(6,*) 'The program found a core state below the bottom'
        WRITE(6,*) 'of the valence-band energy contour.'
        WRITE(6,*) 'This state was not given in the input potential.'
        WRITE(6,*) 'The results are very probably wrong.'
        WRITE(6,*) 'Lower the bottom of the contour or increase the'
        WRITE(6,*) 'number of core states in the input potential.'
        WRITE(6,*) 'The program stops in subroutine corel.f.'
        STOP 'Error 2 in corel.f'
        END IF
        END IF

   70 CONTINUE
      IF (NC*IRNUMX.GT.150 .OR. IRNUMX.GT.10) STOP 'corel'
c
c---> integrate core density to get core charge
c
      CALL SIMP3(RHOC,QC,1,NR,DRDI)

 9000 FORMAT (1x,90 ('*'),/,'  n = ',i1,'  l = ',a4,'   nnode = ',i1,
     +       '  spin=',a4,i5,'th cell','    einput = ',1p,d16.8)
 9010 FORMAT (1x,'  einput =',1p,d16.8,'   eout - ein =',1p,d16.8,
     +       '   eoutput = ',1p,d16.8)
      END
