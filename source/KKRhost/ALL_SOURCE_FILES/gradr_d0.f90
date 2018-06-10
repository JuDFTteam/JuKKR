SUBROUTINE gradr(nspin,ist1,mesh,dx,drdi,drdi2,ro,zta,drr,ddrr,  &
        drru,ddrru,rou,irmd)
!-----------------------------------------------------------------
!evaluates d(ro)/dr,d{d(ro)/dr}/dr.
!drr=d(ro)/dr, ddrr=d(drr)/dr.
!coded by T.Asada. Feb.1994.
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!------------------------------------------------------------------
      use mod_types, only: t_inc
      IMPLICIT NONE
!.. Scalar Arguments ..
      DOUBLE PRECISION DX
      INTEGER IRMD,IST1,MESH,NSPIN
!..
!.. Array Arguments ..
DOUBLE PRECISION DDRR(IRMD),DDRRU(IRMD),DRDI(IRMD),DRDI2(IRMD), &
                 DRR(IRMD),DRRU(IRMD),RO(IRMD),ROU(IRMD),ZTA(IRMD)
!..
!.. Local Scalars ..
DOUBLE PRECISION D,DRX,DRX0,DRX1,DRX2,DRX3,DRXU,DRXU0,DRXU1,DRXU2, &
                 DRXU3,DRXX,DRXX0,DRXX1,DRXX2,DRXX3,DRXXU,DRXXU0, &
                 DRXXU1,DRXXU2,DRXXU3,F0,F1,F2,F3,F4,F5,G1,G2,G3, &
                 G4,G5,XLF
INTEGER I,I1,I2,I3,I4,I5,I6,IBH,ICA,ICG,IEX,IGD,IGH,IGL,IHB,IMJ, &
        IP9,IPG,IPW,IST,IVG,IVN,IWR,IXLF,J,NDVPT,NRED
!..
!.. Statement Functions ..
DOUBLE PRECISION F131,F132,F133,F141,F142,F143,F144,F151,F152, &
                 F153,F154,F155,F161,F162,F163,F164,F165,F166, &
                 F231,F232,F233,F241,F242,F243,F244,F251,F252, &
                 F253,F254,F255,F261,F262,F263,F264,F265,F266
!..
!.. Intrinsic Functions ..
      INTRINSIC DBLE
!..
!.. Save statement ..
SAVE NDVPT,IGL,IGH,IMJ,IBH,ICA,ICG,IVN,IPW,IPG,IVG,IP9,IGD,IXLF,
     IEX,XLF,IWR
!     ..
!     .. Data statements ..
!.....-----------------------------------------------------------------
!     double function
!      double precision f131,f132,f133,f141,f142,f143,f144
!      double precision fl61,fl62,fl63,fl64,fl65,fl66
!      double precision fl51,fl52,fl53,fl54,fl55
!      double precision f231,f232,f233,f241,f242,f243,f244
!      double precision f251,f252,f253,f254,f255
!      double precision f261,f262,f263,f264,f265,f266

DATA ndvpt/6/
DATA igl,igh,imj,ibh,ica,icg,ivn,ipw,ipg,ivg,ip9,igd,ixlf,iex,  &
    xlf/0,0,0,0,0,1,0,0,1,0,0,1,0,0,0.00D0/
DATA iwr/0/
!     ..
!     .. Statement Function definitions ..
!  statement functions:

!.....three point formula for the 1st deriv.

!.....four point formula for the 1st deriv.

!.....five point formula for the 1st deriv.

!.....six point formula for the 1st deriv.

!.....three point formula for the 2nd deriv.

!.....four point formula for the 2nd deriv.

!.....five point formula for the 2nd deriv.

!.....six point formula for the 2nd deriv.
f131(f0,f1,f2,d) = (-3D0*f0+4D0*f1-f2)/ (2D0*d)
f132(g1,f0,f1,d) = (-1D0*g1-0D0*f0+f1)/ (2D0*d)
f133(g2,g1,f0,d) = (g2-4D0*g1+3D0*f0)/ (2D0*d)
f141(f0,f1,f2,f3,d) = (-11D0*f0+18D0*f1-9D0*f2+2D0*f3)/ (6D0*d)
f142(g1,f0,f1,f2,d) = (-2D0*g1-3D0*f0+6D0*f1-f2)/ (6D0*d)
f143(g2,g1,f0,f1,d) = (g2-6D0*g1+3D0*f0+2D0*f1)/ (6D0*d)
f144(g3,g2,g1,f0,d) = (-2D0*g3+9D0*g2-18D0*g1+11D0*f0)/ (6D0*d)
f151(f0,f1,f2,f3,f4,d) = (-50D0*f0+96D0*f1-72D0*f2+32D0*f3-6D0*f4)/ (24D0*d)
f152(g1,f0,f1,f2,f3,d) = (-6D0*g1-20D0*f0+36D0*f1-12D0*f2+2D0*f3)/ (24D0*d)
f153(g2,g1,f0,f1,f2,d) = (2D0*g2-16D0*g1-0D0*f0+16D0*f1-2D0*f2)/ (24D0*d)
f154(g3,g2,g1,f0,f1,d) = (-2D0*g3+12D0*g2-36D0*g1+20D0*f0+6D0*f1)/ (24D0*d)
f155(g4,g3,g2,g1,f0,d) = (6D0*g4-32D0*g3+72D0*g2-96D0*g1+50D0*f0)/ (24D0*d)
f161(f0,f1,f2,f3,f4,f5,d) = (-274D0*f0+600D0*f1-600D0*f2+400D0*f3-150D0*f4+  &
    24D0*f5)/ (120D0*d)
f162(g1,f0,f1,f2,f3,f4,d) = (-24D0*g1-130D0*f0+240D0*f1-120D0*f2+40D0*f3-  &
    6D0*f4)/ (120D0*d)
f163(g2,g1,f0,f1,f2,f3,d) =  &
    (6D0*g2-60D0*g1-40D0*f0+120D0*f1-30D0*f2+4D0*f3)/ (120D0*d)
f164(g3,g2,g1,f0,f1,f2,d) =  &
    (-4D0*g3+30D0*g2-120D0*g1+40D0*f0+60D0*f1-6D0*f2)/ (120D0*d)
f165(g4,g3,g2,g1,f0,f1,d) = (6D0*g4-40D0*g3+120D0*g2-240D0*g1+130D0*f0+  &
    24D0*f1)/ (120D0*d)
f166(g5,g4,g3,g2,g1,f0,d) = (-24D0*g5+150D0*g4-400D0*g3+600D0*g2-600D0*g1+  &
    274D0*f0)/ (120D0*d)
f231(f0,f1,f2,d) = (f0-2D0*f1+f2)/ (d*d)
f232(g1,f0,f1,d) = (g1-2D0*f0+f1)/ (d*d)
f233(g2,g1,f0,d) = (g2-2D0*g1+f0)/ (d*d)
f241(f0,f1,f2,f3,d) = (6D0*f0-15D0*f1+12D0*f2-3D0*f3)/ (3D0*d*d)
f242(g1,f0,f1,f2,d) = (3D0*g1-6D0*f0+3D0*f1+0D0*f2)/ (3D0*d*d)
f243(g2,g1,f0,f1,d) = (0D0*g2+3D0*g1-6D0*f0+3D0*f1)/ (3D0*d*d)
f244(g3,g2,g1,f0,d) = (-3D0*g3+2D0*g2+15D0*g1+6D0*f0)/ (3D0*d*d)
f251(f0,f1,f2,f3,f4,d) = (35D0*f0-104D0*f1+114D0*f2-56D0*f3+11D0*f4)/  &
    (12D0*d*d)
f252(g1,f0,f1,f2,f3,d) = (11D0*g1-20D0*f0+6D0*f1+4D0*f2-f3)/ (12D0*d*d)
f253(g2,g1,f0,f1,f2,d) = (-g2+16D0*g1-30D0*f0+16D0*f1-f2)/ (12D0*d*d)
f254(g3,g2,g1,f0,f1,d) = (-g3+4D0*g2+6D0*g1-20D0*f0+11D0*f1)/ (12D0*d*d)
f255(g4,g3,g2,g1,f0,d) = (11D0*g4-56D0*g3+114D0*g2-104D0*g1+35D0*f0)/  &
    (12D0*d*d)
f261(f0,f1,f2,f3,f4,f5,d) = (225D0*f0-770D0*f1+1070D0*f2-780D0*f3+305*f4-  &
    50D0*f5)/ (60D0*d*d)
f262(g1,f0,f1,f2,f3,f4,d) =  &
    (50D0*g1-75D0*f0-20D0*f1+70D0*f2-30D0*f3+5D0*f4)/ (60D0*d*d)
f263(g2,g1,f0,f1,f2,f3,d) =  &
    (-5D0*g2+80D0*g1-150D0*f0+80D0*f1-5D0*f2+0D0*f3)/ (60D0*d*d)
f264(g3,g2,g1,f0,f1,f2,d) =  &
    (0D0*g3-5D0*g2+80D0*g1-150D0*f0+80D0*f1-5D0*f2)/ (60D0*d*d)
f265(g4,g3,g2,g1,f0,f1,d) =  &
    (5D0*g4-30D0*g3+70D0*g2-20D0*g1-75D0*f0+50D0*f1)/ (60D0*d*d)
f266(g5,g4,g3,g2,g1,f0,d) = (-50D0*g5+305D0*g4-780D0*g3+1070D0*g2-770D0*g1+  &
    225D0*f0)/ (60D0*d*d)
!     ..

!.....-----------------------------------------------------------------
IF ((iwr == 1).AND.(t_inc%i_write>0)) WRITE (1337,FMT=  &
    '(/''  IGL,IGH,IMJ,IHB,ICA,ICG,IVN,IPW,IPG,'',            ''IVG,IP  &
    9,igd,ixlf,iex,xlf='',14I2,F10.4)') igl,igh,imj,ihb,ica,icg,ivn,  &
    ipw,ipg,ivg,ip9,igd,ixlf,iex,xlf
iwr = 0

ist = ist1
!     write(6,*) 'ndvpt ist mesh dx drdi2' ,ndvpt,ist,mesh,dx,
!    &            drdi2(ist)

IF (ndvpt < 3 .OR. ndvpt > 6) THEN
  WRITE (6,FMT=9000) ndvpt
  STOP 18
END IF
!.....
!.....ro: total(core+val)(up+down) charge density.

DO  i = ist,mesh
  rou(i) = ro(i)* (zta(i)+1.d0)/2.d0
END DO
!.....
IF (igd <= 0) THEN
  
  DO  i = ist,mesh
    drr(i) = 0.d0
    ddrr(i) = 0.d0
    drru(i) = 0.d0
    ddrru(i) = 0.d0
  END DO
  GO TO 60
  
END IF

i1 = ist
i2 = ist + 1
i3 = ist + 2
i4 = ist + 3
i5 = ist + 4
i6 = ist + 5

!.....drr:d(ro)/dr, ddrr=d(d(ro)/dr)/dr
!c.... drru,ddrru: for up   spin,
!.....

IF (nspin == 1) GO TO 40
!.....
IF (ndvpt == 3) THEN
  
  drx1 = f131(ro(i1),ro(i2),ro(i3),dx)
  drxu1 = f131(rou(i1),rou(i2),rou(i3),dx)
  drxx1 = f231(ro(i1),ro(i2),ro(i3),dx)
  drxxu1 = f231(rou(i1),rou(i2),rou(i3),dx)
  
ELSE IF (ndvpt == 4) THEN
  
  drx1 = f141(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  drxu1 = f141(rou(i1),rou(i2),rou(i3),rou(i4),dx)
  drxx1 = f241(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  drxxu1 = f241(rou(i1),rou(i2),rou(i3),rou(i4),dx)
  drx2 = f142(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  drxu2 = f142(rou(i1),rou(i2),rou(i3),rou(i4),dx)
  drxx2 = f242(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  drxxu2 = f242(rou(i1),rou(i2),rou(i3),rou(i4),dx)
  
ELSE IF (ndvpt == 5) THEN
  
  drx1 = f151(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  drxu1 = f151(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)
  drxx1 = f251(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  drxxu1 = f251(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)
  drx2 = f152(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  drxu2 = f152(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)
  drxx2 = f252(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  drxxu2 = f252(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)
  
ELSE IF (ndvpt == 6) THEN
  
  drx1 = f161(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxu1 = f161(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
  drxx1 = f261(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxxu1 = f261(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6), dx)
  drx2 = f162(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxu2 = f162(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
  drxx2 = f262(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxxu2 = f262(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6), dx)
  drx3 = f163(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxu3 = f163(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
  drxx3 = f263(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxxu3 = f263(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6), dx)
  
END IF

drr(i1) = drx1/drdi(i1)
ddrr(i1) = (drxx1-drx1*drdi2(i1))/drdi(i1)**2
drru(i1) = drxu1/drdi(i1)
ddrru(i1) = (drxxu1-drxu1*drdi2(i1))/drdi(i1)**2

IF (ndvpt > 3) THEN
  
  drr(i2) = drx2/drdi(i2)
  ddrr(i2) = (drxx2-drx2*drdi2(i2))/drdi(i2)**2
  drru(i2) = drxu2/drdi(i2)
  ddrru(i2) = (drxxu2-drxu2*drdi2(i2))/drdi(i2)**2
  
  IF (ndvpt == 6) THEN
    drr(i3) = drx3/drdi(i3)
    ddrr(i3) = (drxx3-drx3*drdi2(i3))/drdi(i3)**2
    drru(i3) = drxu3/drdi(i3)
    ddrru(i3) = (drxxu3-drxu3*drdi2(i3))/drdi(i3)**2
  END IF
  
END IF

nred = DBLE(ndvpt)/2 + .1D0

DO  j = nred + ist,mesh - nred
  
  IF (ndvpt == 3) THEN
    
    drx = f132(ro(j-1),ro(j),ro(j+1),dx)
    drxu = f132(rou(j-1),rou(j),rou(j+1),dx)
    drxx = f232(ro(j-1),ro(j),ro(j+1),dx)
    drxxu = f232(rou(j-1),rou(j),rou(j+1),dx)
    
  ELSE IF (ndvpt == 4) THEN
    
    drx = f142(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxu = f142(rou(j-1),rou(j),rou(j+1),rou(j+2),dx)
    drxx = f242(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxxu = f242(rou(j-1),rou(j),rou(j+1),rou(j+2),dx)
    
  ELSE IF (ndvpt == 5) THEN
    
    drx = f153(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxu = f153(rou(j-2),rou(j-1),rou(j),rou(j+1),rou(j+2),dx)
    drxx = f253(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxxu = f253(rou(j-2),rou(j-1),rou(j),rou(j+1),rou(j+2),dx)
    
  ELSE IF (ndvpt == 6) THEN
    
    drx = f164(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxu = f164(rou(j-3),rou(j-2),rou(j-1),rou(j),rou(j+1), rou(j+2),dx)
    drxx = f264(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxxu = f264(rou(j-3),rou(j-2),rou(j-1),rou(j),rou(j+1), rou(j+2),dx)
    
  END IF
  
  drr(j) = drx/drdi(j)
  ddrr(j) = (drxx-drx*drdi2(j))/drdi(j)**2
  drru(j) = drxu/drdi(j)
  ddrru(j) = (drxxu-drxu*drdi2(j))/drdi(j)**2
  
END DO
!.....
IF (ndvpt == 3) THEN
  
  drx0 = f133(ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxu0 = f133(rou(mesh-2),rou(mesh-1),rou(mesh),dx)
  drxx0 = f233(ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxxu0 = f233(rou(mesh-2),rou(mesh-1),rou(mesh),dx)
  
ELSE IF (ndvpt == 4) THEN
  
  drx1 = f143(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxu1 = f143(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)
  drxx1 = f243(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxxu1 = f243(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)
  drx0 = f144(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxu0 = f144(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)
  drxx0 = f244(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxxu0 = f244(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)
  
ELSE IF (ndvpt == 5) THEN
  
  drx1 = f154(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  drxu1 = f154(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1), rou(mesh),dx)
  drxx1 = f254(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  drxxu1 = f254(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1), rou(mesh),dx)
  drx0 = f155(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  drxu0 = f155(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1), rou(mesh),dx)
  drxx0 = f255(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  drxxu0 = f255(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1), rou(mesh),dx)
  
ELSE IF (ndvpt == 6) THEN
  
  drx2 = f164(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxu2 = f164(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),  &
      rou(mesh-1),rou(mesh),dx)
  drxx2 = f264(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxxu2 = f264(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),  &
      rou(mesh-1),rou(mesh),dx)
  
  drx1 = f165(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxu1 = f165(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),  &
      rou(mesh-1),rou(mesh),dx)
  drxx1 = f265(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxxu1 = f265(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),  &
      rou(mesh-1),rou(mesh),dx)
  
  drx0 = f166(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxu0 = f166(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),  &
      rou(mesh-1),rou(mesh),dx)
  drxx0 = f266(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxxu0 = f266(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),  &
      rou(mesh-1),rou(mesh),dx)
  
  
END IF

IF (ndvpt > 3) THEN
  
  IF (ndvpt == 6) THEN
    drr(mesh-2) = drx2/drdi(mesh-2)
    drru(mesh-2) = drxu2/drdi(mesh-2)
    ddrr(mesh-2) = (drxx2-drx2*drdi2(mesh-2))/drdi(mesh-2)**2
    ddrru(mesh-2) = (drxxu2-drxu2*drdi2(mesh-2))/drdi(mesh-2)**2
  END IF
  
  drr(mesh-1) = drx1/drdi(mesh-1)
  drru(mesh-1) = drxu1/drdi(mesh-1)
  ddrr(mesh-1) = (drxx1-drx1*drdi2(mesh-1))/drdi(mesh-1)**2
  ddrru(mesh-1) = (drxxu1-drxu1*drdi2(mesh-1))/drdi(mesh-1)**2
  
END IF

drr(mesh) = drx0/drdi(mesh)
drru(mesh) = drxu0/drdi(mesh)
ddrr(mesh) = (drxx0-drx0*drdi2(mesh))/drdi(mesh)**2
ddrru(mesh) = (drxxu0-drxu0*drdi2(mesh))/drdi(mesh)**2

GO TO 60

40 CONTINUE

!.....
IF (ndvpt == 3) THEN
  
  drx1 = f131(ro(i1),ro(i2),ro(i3),dx)
  drxx1 = f231(ro(i1),ro(i2),ro(i3),dx)
  
ELSE IF (ndvpt == 4) THEN
  
  drx1 = f141(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  drxx1 = f241(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  drx2 = f142(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  drxx2 = f242(ro(i1),ro(i2),ro(i3),ro(i4),dx)
  
ELSE IF (ndvpt == 5) THEN
  
  drx1 = f151(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  drxx1 = f251(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  drx2 = f152(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  drxx2 = f252(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
  
ELSE IF (ndvpt == 6) THEN
  
  drx1 = f161(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxx1 = f261(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drx2 = f162(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxx2 = f262(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drx3 = f163(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  drxx3 = f263(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
  
END IF

drr(i1) = drx1/drdi(i1)
ddrr(i1) = (drxx1-drx1*drdi2(i1))/drdi(i1)**2

IF (ndvpt > 3) THEN
  
  drr(i2) = drx2/drdi(i2)
  ddrr(i2) = (drxx2-drx2*drdi2(i2))/drdi(i2)**2
  
  IF (ndvpt == 6) THEN
    drr(i3) = drx3/drdi(i3)
    ddrr(i3) = (drxx3-drx3*drdi2(i3))/drdi(i3)**2
  END IF
  
END IF

nred = DBLE(ndvpt)/2 + .1D0

IF (mesh-nred <= ist) THEN
  WRITE (6,FMT='(/'' MESH-NRED.LT.IST. MESH,NRED,IST='',3I4)') mesh,nred,ist
  STOP 13
END IF

DO  j = nred + ist,mesh - nred
  
  IF (ndvpt == 3) THEN
    
    drx = f132(ro(j-1),ro(j),ro(j+1),dx)
    drxx = f232(ro(j-1),ro(j),ro(j+1),dx)
    
  ELSE IF (ndvpt == 4) THEN
    
    drx = f142(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxx = f242(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    
  ELSE IF (ndvpt == 5) THEN
    
    drx = f153(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxx = f253(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    
  ELSE IF (ndvpt == 6) THEN
    
    drx = f164(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    drxx = f264(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
    
  END IF
  
  drr(j) = drx/drdi(j)
  ddrr(j) = (drxx-drx*drdi2(j))/drdi(j)**2
!           write(6,9000) j,drr(j)
!9000       format(1x,' j drr(j)',i5,e15.5)
END DO
!.....
IF (ndvpt == 3) THEN
  
  drx0 = f133(ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxx0 = f233(ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  
ELSE IF (ndvpt == 4) THEN
  
  drx1 = f143(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxx1 = f243(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drx0 = f144(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  drxx0 = f244(ro(mesh-3),ro(mesh-2),ro(mesh-1),ro(mesh),dx)
  
ELSE IF (ndvpt == 5) THEN
  
  drx1 = f154(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  drxx1 = f254(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  drx0 = f155(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  drxx0 = f255(ro(mesh-4),ro(mesh-3),ro(mesh-2),ro(mesh-1), ro(mesh),dx)
  
ELSE IF (ndvpt == 6) THEN
  
  drx2 = f164(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxx2 = f264(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  
  drx1 = f165(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxx1 = f265(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  
  drx0 = f166(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  drxx0 = f266(ro(mesh-5),ro(mesh-4),ro(mesh-3),ro(mesh-2),  &
      ro(mesh-1),ro(mesh),dx)
  
  
END IF

IF (ndvpt > 3) THEN
  
  IF (ndvpt == 6) THEN
    drr(mesh-2) = drx2/drdi(mesh-2)
    ddrr(mesh-2) = (drxx2-drx2*drdi2(mesh-2))/drdi(mesh-2)**2
  END IF
  
  drr(mesh-1) = drx1/drdi(mesh-1)
  ddrr(mesh-1) = (drxx1-drx1*drdi2(mesh-1))/drdi(mesh-1)**2
  
END IF

drr(mesh) = drx0/drdi(mesh)
ddrr(mesh) = (drxx0-drx0*drdi2(mesh))/drdi(mesh)**2

60 CONTINUE



!      write(6,8000) nspin,ist1,mesh,dx
!8000 format(1x,' nspin ist1 mesh dx',3i5,2d20.10)
!     write(6,8001) (ro(kk),drr(kk),ddrr(kk),
!    &  drdi(kk),drdi2(kk), kk=ist1,mesh,20)
!8001 format(1x,' ro drr ddrr drdi drdi2',5f12.5)
RETURN
9000 FORMAT (/,' ndvpt should be ge.4 .or. le.6. ndvpt=',i3)
END SUBROUTINE gradr
