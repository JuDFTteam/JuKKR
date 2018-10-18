  !-------------------------------------------------------------------------------
  !> Summary: Calculate the spin-polarized exchange-correlation potential and the spin-polarized exchange-correlation energy .
  !> Author: B. Drittler
  !> Date: June 1987
  !> Category: xc-potential, KKRimp
  !> Deprecated: False 
  !> Calculate the spin-polarized exchange-correlation potential and the 
  !> spin-polarized exchange-correlation energy kxc=0 means : spin-polarized 
  !> exchange-correlation potential U. Von Barth and l.hedin, J. Phys. C5,1629 (1972)
  !> with parametrization of Moruzzi, Janak, Williams 
  !>
  !> kxc=1 means : spin-polarized exchange-correlation potential 
  !> U. Von Barth and L. Hedin, J. Phys. C5,1629 (1972)
  !> with parametrization of Von Barth, Hedin
  !> 
  !> use as input the density generated on an angular mesh (see subroutine `vxclm`).
  !> `fpirho(.,1)` contains the charge density times \(4\pi\) and `fpirho(.,2)` 
  !> the spin density times \(4 \pi\).
  !> Then the ex.-cor. potential and the ex.-cor. energy on those mesh points is calculated .
  !> the spin-down potential is stored in `vxc(.,1)`.
  !-------------------------------------------------------------------------------
      SUBROUTINE VXCSPO(EXC,FPIRHO,VXC,KXC,IJEND,IJD)
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,KXC,IJD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EXC(*),FPIRHO(IJD,2),VXC(IJD,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CEX,CF,CFLN,CFMJW,CFVBH,CP,CPLN,CPMJW,CPVBH,D1,
     +                 D2,DCFX,EXCFRS,EXCPRS,EXFRS,EXPRS,FAC,FF,ONTHRD,
     +                 RF,RFMJW,RFVBH,RP,RPMJW,RPVBH,RS,SMAG,TE1B3,VXCC,
     +                 X,XFAC
      INTEGER IJ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DSIGN,LOG,MAX,MIN
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION F
C     ..
C     .. Save statement ..
      SAVE CPMJW,CFMJW,RPMJW,RFMJW,CPVBH,CFVBH,RPVBH,RFVBH,FF,CEX,
     +     ONTHRD,TE1B3
C     ..
C     .. Data statements ..
c
c---> ff=1/(2**(1/3)-1) , cex=2*(3/(2*pi))**(2/3) , te1b3=2**(1/3)
c
      DATA CPMJW,CFMJW,RPMJW,RFMJW/0.045D0,0.0225D0,21.D0,
     +     52.916684096D0/
      DATA CPVBH,CFVBH,RPVBH,RFVBH/0.0504D0,0.0254D0,30.D0,75.D0/
      DATA FF,CEX/3.847322101863D0,1.221774115422D0/
      DATA ONTHRD,TE1B3/0.333333333333D0,1.259921049899D0/
C     ..
C     .. Statement Function definitions ..
c
      F(X) = (1.D0+X*X*X)*LOG(1.D0+1.D0/X) + 0.5D0*X - X*X - 1.0D0/3.0D0
C     ..
c
c---> get key dependent the right parameters
c
      IF (KXC.EQ.1) THEN
        CP = CPVBH
        CF = CFVBH
        RP = RPVBH
        RF = RFVBH

      ELSE
        CP = CPMJW
        CF = CFMJW
        RP = RPMJW
        RF = RFMJW
      END IF
c
c---> loop over the angular mesh points
c
      DO 10 IJ = 1,IJEND
        FPIRHO(IJ,1) = MAX(1.0D-10,FPIRHO(IJ,1))
        SMAG = DSIGN(1.0D0,FPIRHO(IJ,2))
        FPIRHO(IJ,2) = SMAG*MIN(FPIRHO(IJ,1)-1.0D-10,ABS(FPIRHO(IJ,2)))
        RS = (3.D0/FPIRHO(IJ,1))**ONTHRD
        CPLN = CP*LOG(1.D0+RP/RS)
        CFLN = CF*LOG(1.D0+RF/RS)
        DCFX = (CF*F(RS/RF)-CP*F(RS/RP))*4.D0*ONTHRD
        D1 = (1.D0+FPIRHO(IJ,2)/FPIRHO(IJ,1))**ONTHRD
        D2 = (1.D0-FPIRHO(IJ,2)/FPIRHO(IJ,1))**ONTHRD
        FAC = (D1**4+D2**4-2.D0)*0.5D0
c
c---> calculate ex.-cor. energy
c
        EXPRS = -0.75D0*CEX/RS
        EXFRS = EXPRS*TE1B3
        EXCPRS = EXPRS - CP*F(RS/RP)
        EXCFRS = EXFRS - CF*F(RS/RF)
        EXC(IJ) = EXCPRS + (EXCFRS-EXCPRS)*FAC*FF
c
c---> calculate ex.-cor. potential
c
        VXCC = -CPLN + (FAC* (CPLN-CFLN+DCFX)+DCFX)*FF
        XFAC = -CEX/RS - DCFX*FF
        VXC(IJ,2) = VXCC + D1*XFAC
        VXC(IJ,1) = VXCC + D2*XFAC
   10 CONTINUE
      END SUBROUTINE
