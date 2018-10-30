      SUBROUTINE VOSKO(EXC,FPIRHO,VXC,IJEND,IJD)
c-----------------------------------------------------------------------
c     calculate the spin-polarized exchange-correlation potential
c     and the spin-polarized exchange-correlation energy from
c     ceperley-alder ( parametrization of vosko, wilk and nusair )
c                                            ( m. manninen )
c     use as input the density generated on an angular mesh (see
c     subroutine vxclm) . fpirho(.,1) contains the charge density
c     times 4 pi and fpirho(.,2) the spin density times 4 pi .
c     then the ex.-cor. potential and the ex.-cor. energy on those
c     mesh points is calculated .
c     the spin-down potential is stored in vxc(.,1) .
c
c                                  b.drittler    june 1987
c-----------------------------------------------------------------------
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,IJD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EXC(*),FPIRHO(IJD,2),VXC(IJD,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AF,AP,ATNF,ATNP,BETA,BF,BP,CBRT1,CBRT2,CF,CF1,
     +                 CF2,CF3,CP,CP1,CP2,CP3,DBETA,DFS,DUC,DUC1,DUC2,
     +                 EC,ECF,ECP,FS,ONTHRD,QF,QP,RS,S,S4,SMAG,TF1,TP1,
     +                 UC0,UC1,UC10,UC2,UC20,UCF,UCP,X,XF0,XFX,XP0,XPX
      INTEGER IJ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,LOG,MAX,MIN,SIGN,SQRT
C     ..
C     .. Save statement ..
      SAVE AP,XP0,BP,CP,QP,CP1,CP2,CP3,AF,XF0,BF,CF,QF,CF1,CF2,CF3
C     ..
C     .. Data statements ..
      DATA AP,XP0,BP,CP,QP,CP1,CP2,CP3/0.0621814D0,-0.10498D0,3.72744D0,
     +     12.9352D0,6.1519908D0,1.2117833D0,1.1435257D0,-0.031167608D0/
      DATA AF,XF0,BF,CF,QF,CF1,CF2,CF3/0.0310907D0,-0.32500D0,7.06042D0,
     +     18.0578D0,4.7309269D0,2.9847935D0,2.7100059D0,-0.1446006D0/
C     ..
c
      ONTHRD = 1.0D0/3.0D0
c
c---> loop over the angular mesh points
c
      DO 10 IJ = 1,IJEND
        FPIRHO(IJ,1) = MAX(1.0D-10,FPIRHO(IJ,1))
        SMAG = SIGN(1.0D0,FPIRHO(IJ,2))
        FPIRHO(IJ,2) = SMAG*MIN(FPIRHO(IJ,1)-1.0D-10,ABS(FPIRHO(IJ,2)))
        RS = (3.D0/FPIRHO(IJ,1))**ONTHRD
        S = FPIRHO(IJ,2)/FPIRHO(IJ,1)
        X = SQRT(RS)
        XPX = X*X + BP*X + CP
        XFX = X*X + BF*X + CF
        S4 = S**4 - 1.D0
        CBRT1 = (1.D0+S)** (1.D0/3.D0)
        CBRT2 = (1.D0-S)** (1.D0/3.D0)
        FS = ((1.D0+S)** (4.D0/3.D0)+ (1.D0-S)** (4.D0/3.D0)-2.D0)/
     +       (2.D0** (4.D0/3.D0)-2.D0)
        BETA = 1.D0/ (2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
        DFS = 4.D0/3.D0* (CBRT1-CBRT2)/ (2.D0** (4.D0/3.D0)-2.D0)
        DBETA = - (0.27402D0*X+0.09873D0+1.591D0/X)*BETA**2
        ATNP = ATAN(QP/ (2.D0*X+BP))
        ATNF = ATAN(QF/ (2.D0*X+BF))
        ECP = AP* (LOG(X*X/XPX)+CP1*ATNP-
     +        CP3* (LOG((X-XP0)**2/XPX)+CP2*ATNP))
        ECF = AF* (LOG(X*X/XFX)+CF1*ATNF-
     +        CF3* (LOG((X-XF0)**2/XFX)+CF2*ATNF))
        EC = ECP + FS* (ECF-ECP)* (1.D0+S4*BETA)
c
c---> calculate ex.-cor. energy
c
        EXC(IJ) = EC - 0.9163306D0/RS - 0.2381735D0/RS*FS
        TP1 = (X*X+BP*X)/XPX
        TF1 = (X*X+BF*X)/XFX
        UCP = ECP - AP/3.D0* (1.D0-TP1-CP3* (X/ (X-XP0)-TP1-XP0*X/XPX))
        UCF = ECF - AF/3.D0* (1.D0-TF1-CF3* (X/ (X-XF0)-TF1-XF0*X/XFX))
        UC0 = UCP + (UCF-UCP)*FS
        UC10 = UC0 - (ECF-ECP)* (S-1.D0)*DFS
        UC20 = UC0 - (ECF-ECP)* (S+1.D0)*DFS
        DUC = (UCF-UCP)*BETA*S4*FS + (ECF-ECP)* (-RS/3.D0)*DBETA*S4*FS
        DUC1 = DUC - (ECF-ECP)*BETA* (S-1.D0)* (4.D0*S**3*FS+S4*DFS)
        DUC2 = DUC - (ECF-ECP)*BETA* (S+1.D0)* (4.D0*S**3*FS+S4*DFS)
        UC1 = UC10 + DUC1
        UC2 = UC20 + DUC2
c
c---> calculate exc.-cor. potential
c
        VXC(IJ,2) = UC1 - 1.221774D0/RS*CBRT1
        VXC(IJ,1) = UC2 - 1.221774D0/RS*CBRT2
        IF (ABS(FPIRHO(IJ,1)).LE.1.0D-10) THEN
          VXC(IJ,1) = 0.0D0
          VXC(IJ,2) = 0.0D0
        END IF

   10 CONTINUE
      END
