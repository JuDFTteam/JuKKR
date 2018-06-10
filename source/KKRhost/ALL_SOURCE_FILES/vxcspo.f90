SUBROUTINE vxcspo(exc,fpirho,vxc,kxc,ijend,ijd)
!-----------------------------------------------------------------------
!     calculate the spin-polarized exchange-correlation potential
!     and the spin-polarized exchange-correlation energy .
!     kxc=0 means : spin-polarized exchange-correlation potential
!                   u. von barth and l.hedin, j.phys.c5,1629 (1972)
!                   with parametrization of moruzzi,janak,williams
!     kxc=1 means : spin-polarized exchange-correlation potential
!                   u. von barth and l.hedin, j.phys.c5,1629 (1972)
!                   with parametrization of von barth,hedin

!     use as input the density generated on an angular mesh (see
!     subroutine vxclm) . fpirho(.,1) contains the charge density
!     times 4 pi and fpirho(.,2) the spin density times 4 pi .
!     then the ex.-cor. potential and the ex.-cor. energy on those
!     mesh points is calculated .
!     the spin-down potential is stored in vxc(.,1) .

!                                  b.drittler    june 1987
!-----------------------------------------------------------------------
!..
!.. Scalar Arguments ..
      INTEGER IJEND,KXC,IJD
!..
!.. Array Arguments ..
      DOUBLE PRECISION EXC(*),FPIRHO(IJD,2),VXC(IJD,2)
!..
!.. Local Scalars ..
DOUBLE PRECISION CEX,CF,CFLN,CFMJW,CFVBH,CP,CPLN,CPMJW,CPVBH,D1, &
                 D2,DCFX,EXCFRS,EXCPRS,EXFRS,EXPRS,FAC,FF,ONTHRD, &
                 RF,RFMJW,RFVBH,RP,RPMJW,RPVBH,RS,SMAG,TE1B3,VXCC, &
                 X,XFAC
INTEGER IJ
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DSIGN,LOG,MAX,MIN
!..
!.. Statement Functions ..
      DOUBLE PRECISION F
!..
!.. Save statement ..
SAVE CPMJW,CFMJW,RPMJW,RFMJW,CPVBH,CFVBH,RPVBH,RFVBH,FF,CEX, &
     ONTHRD,TE1B3
!..
!.. Data statements ..
!
!---> ff=1/(2**(1/3)-1) , cex=2*(3/(2*pi))**(2/3) , te1b3=2**(1/3)
!
DATA CPMJW,CFMJW,RPMJW,RFMJW/0.045D0,0.0225D0,21.D0, &
     52.916684096D0/
DATA CPVBH,CFVBH,RPVBH,RFVBH/0.0504D0,0.0254D0,30.D0,75.D0/
DATA FF,CEX/3.847322101863D0,1.221774115422D0/
DATA ONTHRD,TE1B3/0.333333333333D0,1.259921049899D0/
!     ..
!     .. Statement Function definitions ..

f(x) = (1.d0+x*x*x)*LOG(1.d0+1.d0/x) + 0.5D0*x - x*x - 1.0D0/3.0D0
!     ..

!---> get key dependent the right parameters

IF (kxc == 1) THEN
  cp = cpvbh
  cf = cfvbh
  rp = rpvbh
  rf = rfvbh
  
ELSE
  cp = cpmjw
  cf = cfmjw
  rp = rpmjw
  rf = rfmjw
END IF

!---> loop over the angular mesh points

DO  ij = 1,ijend
  fpirho(ij,1) = MAX(1.0D-10,fpirho(ij,1))
  smag = DSIGN(1.0D0,fpirho(ij,2))
  fpirho(ij,2) = smag*MIN(fpirho(ij,1)-1.0D-10,ABS(fpirho(ij,2)))
  rs = (3.d0/fpirho(ij,1))**onthrd
  cpln = cp*LOG(1.d0+rp/rs)
  cfln = cf*LOG(1.d0+rf/rs)
  dcfx = (cf*f(rs/rf)-cp*f(rs/rp))*4.d0*onthrd
  d1 = (1.d0+fpirho(ij,2)/fpirho(ij,1))**onthrd
  d2 = (1.d0-fpirho(ij,2)/fpirho(ij,1))**onthrd
  fac = (d1**4+d2**4-2.d0)*0.5D0
  
!---> calculate ex.-cor. energy
  
  exprs = -0.75D0*cex/rs
  exfrs = exprs*te1b3
  excprs = exprs - cp*f(rs/rp)
  excfrs = exfrs - cf*f(rs/rf)
  exc(ij) = excprs + (excfrs-excprs)*fac*ff
  
!---> calculate ex.-cor. potential
  
  vxcc = -cpln + (fac* (cpln-cfln+dcfx)+dcfx)*ff
  xfac = -cex/rs - dcfx*ff
  vxc(ij,2) = vxcc + d1*xfac
  vxc(ij,1) = vxcc + d2*xfac
END DO
END SUBROUTINE vxcspo
