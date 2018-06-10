SUBROUTINE vosko(exc,fpirho,vxc,ijend,ijd)
!-----------------------------------------------------------------------
!     calculate the spin-polarized exchange-correlation potential
!     and the spin-polarized exchange-correlation energy from
!     ceperley-alder ( parametrization of vosko, wilk and nusair )
!                                            ( m. manninen )
!     use as input the density generated on an angular mesh (see
!     subroutine vxclm) . fpirho(.,1) contains the charge density
!     times 4 pi and fpirho(.,2) the spin density times 4 pi .
!     then the ex.-cor. potential and the ex.-cor. energy on those
!     mesh points is calculated .
!     the spin-down potential is stored in vxc(.,1) .

!                                  b.drittler    june 1987
!-----------------------------------------------------------------------
implicit none
!..
!.. Scalar Arguments ..
      INTEGER IJEND,IJD
!..
!.. Array Arguments ..
      DOUBLE PRECISION EXC(*),FPIRHO(IJD,2),VXC(IJD,2)
!..
!.. Local Scalars ..
DOUBLE PRECISION AF,AP,ATNF,ATNP,BETA,BF,BP,CBRT1,CBRT2,CF,CF1, &
                 CF2,CF3,CP,CP1,CP2,CP3,DBETA,DFS,DUC,DUC1,DUC2, &
                 EC,ECF,ECP,FS,ONTHRD,QF,QP,RS,S,S4,SMAG,TF1,TP1, &
                 UC0,UC1,UC10,UC2,UC20,UCF,UCP,X,XF0,XFX,XP0,XPX
INTEGER IJ
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,LOG,MAX,MIN,SIGN,SQRT
!..
!.. Save statement ..
      SAVE AP,XP0,BP,CP,QP,CP1,CP2,CP3,AF,XF0,BF,CF,QF,CF1,CF2,CF3
!..
!.. Data statements ..
DATA ap,xp0,bp,cp,qp,cp1,cp2,cp3/0.0621814D0,-0.10498D0,3.72744D0,  &
    12.9352D0,6.1519908D0,1.2117833D0,1.1435257D0,-0.031167608D0/
DATA af,xf0,bf,cf,qf,cf1,cf2,cf3/0.0310907D0,-0.32500D0,7.06042D0,  &
    18.0578D0,4.7309269D0,2.9847935D0,2.7100059D0,-0.1446006D0/
!     ..

onthrd = 1.0D0/3.0D0

!---> loop over the angular mesh points

DO  ij = 1,ijend
  fpirho(ij,1) = MAX(1.0D-10,fpirho(ij,1))
  smag = SIGN(1.0D0,fpirho(ij,2))
  fpirho(ij,2) = smag*MIN(fpirho(ij,1)-1.0D-10,ABS(fpirho(ij,2)))
  rs = (3.d0/fpirho(ij,1))**onthrd
  s = fpirho(ij,2)/fpirho(ij,1)
  x = SQRT(rs)
  xpx = x*x + bp*x + cp
  xfx = x*x + bf*x + cf
  s4 = s**4 - 1.d0
  cbrt1 = (1.d0+s)** (1.d0/3.d0)
  cbrt2 = (1.d0-s)** (1.d0/3.d0)
  fs = ((1.d0+s)** (4.d0/3.d0)+ (1.d0-s)** (4.d0/3.d0)-2.d0)/  &
      (2.d0** (4.d0/3.d0)-2.d0)
  beta = 1.d0/ (2.74208D0+3.182D0*x+0.09873D0*x*x+0.18268D0*x**3)
  dfs = 4.d0/3.d0* (cbrt1-cbrt2)/ (2.d0** (4.d0/3.d0)-2.d0)
  dbeta = - (0.27402D0*x+0.09873D0+1.591D0/x)*beta**2
  atnp = ATAN(qp/ (2.d0*x+bp))
  atnf = ATAN(qf/ (2.d0*x+bf))
  ecp = ap* (LOG(x*x/xpx)+cp1*atnp- cp3* (LOG((x-xp0)**2/xpx)+cp2*atnp))
  ecf = af* (LOG(x*x/xfx)+cf1*atnf- cf3* (LOG((x-xf0)**2/xfx)+cf2*atnf))
  ec = ecp + fs* (ecf-ecp)* (1.d0+s4*beta)
  
!---> calculate ex.-cor. energy
  
  exc(ij) = ec - 0.9163306D0/rs - 0.2381735D0/rs*fs
  tp1 = (x*x+bp*x)/xpx
  tf1 = (x*x+bf*x)/xfx
  ucp = ecp - ap/3.d0* (1.d0-tp1-cp3* (x/ (x-xp0)-tp1-xp0*x/xpx))
  ucf = ecf - af/3.d0* (1.d0-tf1-cf3* (x/ (x-xf0)-tf1-xf0*x/xfx))
  uc0 = ucp + (ucf-ucp)*fs
  uc10 = uc0 - (ecf-ecp)* (s-1.d0)*dfs
  uc20 = uc0 - (ecf-ecp)* (s+1.d0)*dfs
  duc = (ucf-ucp)*beta*s4*fs + (ecf-ecp)* (-rs/3.d0)*dbeta*s4*fs
  duc1 = duc - (ecf-ecp)*beta* (s-1.d0)* (4.d0*s**3*fs+s4*dfs)
  duc2 = duc - (ecf-ecp)*beta* (s+1.d0)* (4.d0*s**3*fs+s4*dfs)
  uc1 = uc10 + duc1
  uc2 = uc20 + duc2
  
!---> calculate exc.-cor. potential
  
  vxc(ij,2) = uc1 - 1.221774D0/rs*cbrt1
  vxc(ij,1) = uc2 - 1.221774D0/rs*cbrt2
  IF (ABS(fpirho(ij,1)) <= 1.0D-10) THEN
    vxc(ij,1) = 0.0D0
    vxc(ij,2) = 0.0D0
  END IF
  
END DO
END SUBROUTINE vosko
