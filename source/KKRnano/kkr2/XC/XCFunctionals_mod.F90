!> @author Modularisation: Marcel Bornemann

module XCFunctionals_mod
  implicit none
  private

  public :: mkxcpe_pw91
  public :: mkxcpe_pbe
  public :: fpexcpbe
! public :: excpbex
! public :: excpbec
  public :: excgcor2
 
  contains

SUBROUTINE mkxcpe_pw91(nspin,ir,np,l1max,rv,rholm,vxcp,excp,thet,ylm,  &
        dylmt1,dylmt2,dylmf1,dylmf2,dylmtf,drrl,ddrrl,  &
        drrul,ddrrul,irmd,lmpotd)

! Code converted using TO_F90 by Alan Miller
! Date: 2015-07-02  Time: 16:01:44

!     .. Parameters ..

INTEGER, PARAMETER :: ijd=434

!     ..
!     .. Scalar Arguments ..

INTEGER, INTENT(IN)                      :: nspin
INTEGER, INTENT(IN)                      :: ir
INTEGER, INTENT(IN)                      :: np
INTEGER, INTENT(IN)                      :: l1max
DOUBLE PRECISION, INTENT(IN)             :: rv
INTEGER, INTENT(IN OUT)                  :: irmd
INTEGER, INTENT(IN OUT)                  :: lmpotd

!     ..
!     .. Array Arguments ..

DOUBLE PRECISION, INTENT(IN)             :: rholm(lmpotd,2)
DOUBLE PRECISION, INTENT(OUT)            :: vxcp(ijd,2)
DOUBLE PRECISION, INTENT(OUT)            :: excp(ijd)
DOUBLE PRECISION, INTENT(IN OUT)         :: thet(ijd)
DOUBLE PRECISION, INTENT(IN)             :: ylm(ijd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: dylmt1(ijd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: dylmt2(ijd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: dylmf1(ijd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: dylmf2(ijd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: dylmtf(ijd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: drrl(irmd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: ddrrl(irmd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: drrul(irmd,lmpotd)
DOUBLE PRECISION, INTENT(IN)             :: ddrrul(irmd,lmpotd)

!     ..
!     .. Local Scalars ..
DOUBLE PRECISION :: agr,agrd,agru,cedg,cedl,chg,cosx,dagrf,dagrfd,  &
    dagrfu,dagrr,dagrrd,dagrru,dagrt,dagrtd,dagrtu,  &
    ddrr,ddrrd,ddrru,df1,df2,drr,drrd,drru,dt1,dt2,  &
    dtf,dzdfs,dzdr,dzdtr,etot0,etota0,g2r,g2rd,g2ru,  &
    gggr,gggrd,gggru,grf,grfd,grfu,grgrd,grgru,grr,  &
    grrd,grru,grt,grtd,grtu,gzgr,ro,rod,rou,  &
    rv2,rv3,rvsin1,rvsin2,ry2,rylm,sint1,sint2,smag,  &
    spi,tant1,vcg1,vcg2,vcl1,vcl2,vtot1,vtot2,vtota1,  &
    vtota2,vxg1,vxg2,vxl1,vxl2,xedg,xedl,zta
DOUBLE PRECISION, parameter :: rdspr = 9.d0
INTEGER :: idspr,im,ip,l1,ll,llmax,lm,lmax,nn,nn1
!     ..
!     .. Local Arrays ..
DOUBLE PRECISION :: ddry(ijd),ddryd(ijd),ddryu(ijd),drdf(ijd),  &
    drdfd(ijd),drdfu(ijd),drdt(ijd),drdtd(ijd),  &
    drdtu(ijd),dry(ijd),dryd(ijd),dryu(ijd),  &
    rdf1(ijd),rdf1d(ijd),rdf1u(ijd),rdf2(ijd),  &
    rdf2d(ijd),rdf2u(ijd),rdt1(ijd),rdt1d(ijd),  &
    rdt1u(ijd),rdt2(ijd),rdt2d(ijd),rdt2u(ijd),  &
    rdtf(ijd),rdtfd(ijd),rdtfu(ijd),ry(ijd),ryd(ijd), ryu(ijd)
!     ..
!     .. External Subroutines ..
EXTERNAL gxcpt
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ABS,COS,MAX,MIN,SIGN,SIN,SQRT,TAN

llmax = l1max*l1max
lmax = l1max - 1


DO  ip = 1,np
  
  ry(ip) = 0.d0
  dry(ip) = 0.d0
  ddry(ip) = 0.d0
  rdt1(ip) = 0.d0
  rdt2(ip) = 0.d0
  rdf1(ip) = 0.d0
  rdf2(ip) = 0.d0
  rdtf(ip) = 0.d0
  drdt(ip) = 0.d0
  drdf(ip) = 0.d0
  
  ryu(ip) = 0.d0
  dryu(ip) = 0.d0
  ddryu(ip) = 0.d0
  rdt1u(ip) = 0.d0
  rdt2u(ip) = 0.d0
  rdf1u(ip) = 0.d0
  rdf2u(ip) = 0.d0
  rdtfu(ip) = 0.d0
  drdtu(ip) = 0.d0
  drdfu(ip) = 0.d0
  
  ryd(ip) = 0.d0
  dryd(ip) = 0.d0
  ddryd(ip) = 0.d0
  rdt1d(ip) = 0.d0
  rdt2d(ip) = 0.d0
  rdf1d(ip) = 0.d0
  rdf2d(ip) = 0.d0
  rdtfd(ip) = 0.d0
  drdtd(ip) = 0.d0
  drdfd(ip) = 0.d0
  
END DO

!     write(6,'(/'' nspin,mesh,np,l1max='',4i5)') nspin,mesh,np,l1max

lm = 0

DO  l1 = 1,l1max
  
  ll = l1 - 1
  
  
  
  DO  im = -ll,ll
    
    lm = lm + 1
    
    
    ro = rholm(lm,1)*2.d0
    rou = ro/2.d0
    rod = rou
    
    IF (nspin /= 1) THEN
      
      ro = rholm(lm,1) + rholm(lm,2)
      rou = rholm(lm,2)
      rod = rholm(lm,1)
!        write(6,9001) ro,rou,rod
    END IF
    drr = drrl(ir,lm)
    ddrr = ddrrl(ir,lm)
    drru = drrul(ir,lm)
    ddrru = ddrrul(ir,lm)
    drrd = drr - drru
    ddrrd = ddrr - ddrru
    
    
    DO  ip = 1,np
      
      rylm = ylm(ip,lm)
      dt1 = dylmt1(ip,lm)
      dt2 = dylmt2(ip,lm)
      df1 = dylmf1(ip,lm)
      df2 = dylmf2(ip,lm)
      dtf = dylmtf(ip,lm)
      
      ry(ip) = ry(ip) + ro*rylm
      dry(ip) = dry(ip) + drr*rylm
      ddry(ip) = ddry(ip) + ddrr*rylm
      
      ryu(ip) = ryu(ip) + rou*rylm
      dryu(ip) = dryu(ip) + drru*rylm
      ddryu(ip) = ddryu(ip) + ddrru*rylm
      
      ryd(ip) = ryd(ip) + rod*rylm
      dryd(ip) = dryd(ip) + drrd*rylm
      ddryd(ip) = ddryd(ip) + ddrrd*rylm
      
      rdt1(ip) = rdt1(ip) + ro*dt1
      rdt2(ip) = rdt2(ip) + ro*dt2
      rdf1(ip) = rdf1(ip) + ro*df1
      rdf2(ip) = rdf2(ip) + ro*df2
      rdtf(ip) = rdtf(ip) + ro*dtf
      drdt(ip) = drdt(ip) + drr*dt1
      drdf(ip) = drdf(ip) + drr*df1
      
      rdt1u(ip) = rdt1u(ip) + rou*dt1
      rdt2u(ip) = rdt2u(ip) + rou*dt2
      rdf1u(ip) = rdf1u(ip) + rou*df1
      rdf2u(ip) = rdf2u(ip) + rou*df2
      rdtfu(ip) = rdtfu(ip) + rou*dtf
      drdtu(ip) = drdtu(ip) + drru*dt1
      drdfu(ip) = drdfu(ip) + drru*df1
      
      rdt1d(ip) = rdt1d(ip) + rod*dt1
      rdt2d(ip) = rdt2d(ip) + rod*dt2
      rdf1d(ip) = rdf1d(ip) + rod*df1
      rdf2d(ip) = rdf2d(ip) + rod*df2
      rdtfd(ip) = rdtfd(ip) + rod*dtf
      drdtd(ip) = drdtd(ip) + drrd*dt1
      drdfd(ip) = drdfd(ip) + drrd*df1
!rc             if (ip.eq.5.or.ip.eq.6) then
!             write(6,9907) lm,rylm,dt1,dt2,df1,df2,dtf
!9907         format(1x,' lmt ',i3,6d12.6)
!             write(6,*) 'nikos',dry(ip),ddry(ip)
!             end if
      
      
    END DO
  END DO
END DO



DO  ip = 1,np
  sint1 = SIN(thet(ip))
  sint2 = sint1**2
  tant1 = TAN(thet(ip))
  IF (sint1 == 0.d0) THEN
    vxcp(ip,1) = 0.d0
    vxcp(ip,2) = 0.d0
    excp(ip) = 0.d0
!          WRITE (6,FMT=*) 'interpolate'
    
! set values later
  ELSE
    rv2 = rv**2
    rv3 = rv**3
    
    
    rvsin1 = rv*sint1
    rvsin2 = rv2*sint2
    
    grr = dry(ip)
    grt = rdt1(ip)/rv
    grf = rdf1(ip)/rvsin1
    ry2 = ry(ip)**2
    
    agr = SQRT(grr**2+grt**2+grf**2)
    
    dagrr = (dry(ip)*ddry(ip)*rv3+ rdt1(ip)* (drdt(ip)*rv-rdt1(ip))+  &
        rdf1(ip)* (drdf(ip)*rv-rdf1(ip))/sint2)/agr/rv3
    
    dagrt = (dry(ip)*drdt(ip)*rv2+rdt1(ip)*rdt2(ip)+  &
        rdf1(ip)* (-rdf1(ip)/tant1+rdtf(ip))/sint2)/ (agr*rv3)
    
    dagrf = (dry(ip)*drdf(ip)*rv2+rdt1(ip)*rdtf(ip)+  &
        rdf1(ip)*rdf2(ip)/sint2)/ (agr*rv3*sint1)
    
    
    dzdr = ((dryu(ip)-dryd(ip))*ry(ip)- (ryu(ip)-ryd(ip))*dry(ip))/ry2
    
    dzdtr = ((rdt1u(ip)-rdt1d(ip))*ry(ip)- (ryu(ip)-ryd(ip))*rdt1(ip))/ry2/rv
    
    dzdfs = ((rdf1u(ip)-rdf1d(ip))*ry(ip)-  &
        (ryu(ip)-ryd(ip))*rdf1(ip))/ry2/rvsin1
    
    g2r = ddry(ip) + 2.d0*dry(ip)/rv +  &
        (rdt2(ip)+rdt1(ip)/tant1+rdf2(ip)/sint2)/rv2
    
    gggr = grr*dagrr + grt*dagrt + grf*dagrf
    
    gzgr = dzdr*grr + dzdtr*grt + dzdfs*grf
    
    
    chg = ry(ip)
    spi = ryu(ip) - ryd(ip)
    chg = MAX(1.0D-12,chg)
    smag = SIGN(1.0D0,spi)
    spi = smag*MIN(chg-1.0D-12,ABS(spi))
    zta = spi/chg
    
    
    grru = dryu(ip)
    grtu = rdt1u(ip)/rv
    grfu = rdf1u(ip)/rvsin1
    
    agru = SQRT(grru**2+grtu**2+grfu**2)
    
    dagrru = (dryu(ip)*ddryu(ip)*rv3+ rdt1u(ip)* (drdtu(ip)*rv-rdt1u(ip))+  &
        rdf1u(ip)* (drdfu(ip)*rv-rdf1u(ip))/sint2)/agru/rv3
    
    dagrtu = (dryu(ip)*drdtu(ip)*rv2+rdt1u(ip)*rdt2u(ip)+  &
        rdf1u(ip)* (-rdf1u(ip)/tant1+rdtfu(ip))/sint2)/ (agru*rv3)
    
    dagrfu = (dryu(ip)*drdfu(ip)*rv2+rdt1u(ip)*rdtfu(ip)+  &
        rdf1u(ip)*rdf2u(ip)/sint2)/ (agru*rv3*sint1)
    
    
    
    g2ru = ddryu(ip) + 2.d0*dryu(ip)/rv +  &
        (rdt2u(ip)+rdt1u(ip)/tant1+rdf2u(ip)/sint2)/rv2
    
    gggru = grru*dagrru + grtu*dagrtu + grfu*dagrfu
    
    grgru = grr*grru + grt*grtu + grf*grfu
    
    
    grrd = dryd(ip)
    grtd = rdt1d(ip)/rv
    grfd = rdf1d(ip)/rvsin1
    
    agrd = SQRT(grrd**2+grtd**2+grfd**2)
    
    dagrrd = (dryd(ip)*ddryd(ip)*rv3+ rdt1d(ip)* (drdtd(ip)*rv-rdt1d(ip))+  &
        rdf1d(ip)* (drdfd(ip)*rv-rdf1d(ip))/sint2)/agrd/rv3
    
    dagrtd = (dryd(ip)*drdtd(ip)*rv2+rdt1d(ip)*rdt2d(ip)+  &
        rdf1d(ip)* (-rdf1d(ip)/tant1+rdtfd(ip))/sint2)/ (agrd*rv3)
    
    dagrfd = (dryd(ip)*drdfd(ip)*rv2+rdt1d(ip)*rdtfd(ip)+  &
        rdf1d(ip)*rdf2d(ip)/sint2)/ (agrd*rv3*sint1)
    
    
    
    g2rd = ddryd(ip) + 2.d0*dryd(ip)/rv +  &
        (rdt2d(ip)+rdt1d(ip)/tant1+rdf2d(ip)/sint2)/rv2
    
    gggrd = grrd*dagrrd + grtd*dagrtd + grfd*dagrfd
    
    grgrd = grr*grrd + grt*grtd + grf*grfd
    
    
    idspr = 0
    IF (rv > rdspr) idspr = 1
    
    
! for debug
    CALL gxcpt(idspr,chg,zta,agr,agru,agrd,g2r,g2ru,g2rd,gggr,  &
        gggru,gggrd,grgru,grgrd,gzgr,vxcp(ip,2),vxcp(ip,1),  &
        excp(ip),vxl1,vxl2,vcl1,vcl2,xedl,cedl,vxg1,vxg2, vcg1,vcg2,xedg,cedg)
    
    
  END IF
  
END DO

! This is expected to work only for the Lebedev points
nn = 0
nn1 = 0
vtot1 = 0.d0
vtot2 = 0.d0
etot0 = 0.d0
vtota1 = 0.d0
vtota2 = 0.d0
etota0 = 0.d0
DO  ip = 1,np
  cosx = COS(thet(ip))
  IF (cosx > 0.99D0 .AND. cosx /= 1.d0) THEN
    nn = nn + 1
    vtot1 = vtot1 + vxcp(ip,1)
    vtot2 = vtot2 + vxcp(ip,2)
    etot0 = etot0 + excp(ip)
!        write(6,*) 'more',ip,vxcp(ip,1),nn
  END IF
  IF (cosx < -0.99D0 .AND. cosx /= -1.d0) THEN
    nn1 = nn1 + 1
    vtota1 = vtota1 + vxcp(ip,1)
    vtota2 = vtota2 + vxcp(ip,2)
    etota0 = etota0 + excp(ip)
!           write(6,*) 'less',ip,vxcp(ip,1),nn1
  END IF
END DO
DO  ip = 1,np
  cosx = COS(thet(ip))
  IF (cosx == 1.d0) THEN
    vxcp(ip,1) = vtot1/nn
    vxcp(ip,2) = vtot2/nn
    excp(ip) = etot0/nn
!     write(6,*) 'averaging ',ip,vxcp(ip,1),vxcp(ip,2),excp(ip)
!     write(6,*) 'averaging1 ',vtot1,vtot2,etot0,nn
  END IF
  IF (cosx == -1.d0) THEN
    vxcp(ip,1) = vtota1/nn1
    vxcp(ip,2) = vtota2/nn1
    excp(ip) = etota0/nn1
!     write(6,*) 'averaging ',ip,cosx,vxcp(ip,1),vxcp(ip,2),excp(ip)
!     write(6,*)'averaging2 ',vtota1,vtota2,etota0,nn1
  END IF
END DO
RETURN
END SUBROUTINE mkxcpe_pw91

SUBROUTINE mkxcpe_pbe(ir,np,rv,rholm,vxcp,excp,ylm,dylmt1,dylmf1,  &
        dylmf2,dylmtf,drrl,ddrrl,drrul,ddrrul,irmd,  &
        lmpotd,lmmax,use_sol)

! Code converted using TO_F90 by Alan Miller
! Date: 2015-07-02  Time: 16:01:52

!     ------------------------------------------------------------------
!     Calculation of the exchange-correlation potential.
!     coded by M. Ogura, Apr. 2015, Munich
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER, PARAMETER :: ijd=434

INTEGER, INTENT(IN)                      :: ir
INTEGER, INTENT(IN)                      :: np
REAL*8, INTENT(IN)                       :: rv
REAL*8, INTENT(IN)                       :: rholm(lmpotd,2)
REAL*8, INTENT(IN OUT)                   :: vxcp(ijd,2)
REAL*8, INTENT(IN OUT)                   :: excp(ijd)
REAL*8, INTENT(IN)                       :: ylm(ijd,lmpotd)
REAL*8, INTENT(IN)                       :: dylmt1(ijd,lmpotd)
REAL*8, INTENT(IN)                       :: dylmf1(ijd,lmpotd)
REAL*8, INTENT(IN)                       :: dylmf2(ijd,lmpotd)
REAL*8, INTENT(IN)                       :: dylmtf(ijd,lmpotd)
REAL*8, INTENT(IN OUT)                   :: drrl(irmd,lmpotd)
REAL*8, INTENT(IN OUT)                   :: ddrrl(irmd,lmpotd)
REAL*8, INTENT(IN)                       :: drrul(irmd,lmpotd)
REAL*8, INTENT(IN)                       :: ddrrul(irmd,lmpotd)
INTEGER, INTENT(IN OUT)                  :: irmd
INTEGER, INTENT(IN OUT)                  :: lmpotd
INTEGER, INTENT(IN)                      :: lmmax
LOGICAL, INTENT(IN)                      :: use_sol



REAL*8 c,pi,s
INTEGER :: ip,ispin,l1,lm,n
REAL*8 d(2),d1(3,2),d2(5,2),dl(2)
REAL*8 :: um
REAL*8 :: bet

pi=4D0*ATAN(1D0)

! Set 'UM' for subroutine 'excpbex' and 'BET' for subroutine 'excpbec' for PBE or PBEsol
IF (use_sol) THEN
um=0.123456790123456D0
bet=0.046D0
ELSE
um=0.2195149727645171D0
bet=0.06672455060314922D0
END IF


!     --- surface integration
DO  ip=1,np
  DO  ispin=1,2
    d(ispin)=0D0
    dl(ispin)=0D0
    DO  n=1,3
      d1(n,ispin)=0D0
    END DO
    DO  n=1,5
      d2(n,ispin)=0D0
    END DO
  END DO
  DO  lm=1,lmmax
    l1=SQRT(DBLE(lm)-5D-1)
    d(1)=d(1)+rholm(lm,1)*ylm(ip,lm)
    d(2)=d(2)+rholm(lm,2)*ylm(ip,lm)
    dl(1)=dl(1)+DBLE(l1*(l1+1))*rholm(lm,1)*ylm(ip,lm)
    dl(2)=dl(2)+DBLE(l1*(l1+1))*rholm(lm,2)*ylm(ip,lm)
    d1(1,2)=d1(1,2)+drrul(ir,lm)*ylm(ip,lm)
    d1(1,1)=d1(1,1)+(drrl(ir,lm)-drrul(ir,lm))*ylm(ip,lm)
    d1(2,1)=d1(2,1)+rholm(lm,1)*dylmt1(ip,lm)
    d1(2,2)=d1(2,2)+rholm(lm,2)*dylmt1(ip,lm)
    d1(3,1)=d1(3,1)+rholm(lm,1)*dylmf1(ip,lm)
    d1(3,2)=d1(3,2)+rholm(lm,2)*dylmf1(ip,lm)
    d2(1,2)=d2(1,2)+ddrrul(ir,lm)*ylm(ip,lm)
    d2(1,1)=d2(1,1)+(ddrrl(ir,lm)-ddrrul(ir,lm))*ylm(ip,lm)
    d2(2,2)=d2(2,2)+drrul(ir,lm)*dylmt1(ip,lm)
    d2(2,1)=d2(2,1)+(drrl(ir,lm)-drrul(ir,lm))*dylmt1(ip,lm)
    d2(3,2)=d2(3,2)+drrul(ir,lm)*dylmf1(ip,lm)
    d2(3,1)=d2(3,1)+(drrl(ir,lm)-drrul(ir,lm))*dylmf1(ip,lm)
    d2(4,1)=d2(4,1)+rholm(lm,1)*dylmtf(ip,lm)
    d2(4,2)=d2(4,2)+rholm(lm,2)*dylmtf(ip,lm)
    d2(5,1)=d2(5,1)+rholm(lm,1)*dylmf2(ip,lm)
    d2(5,2)=d2(5,2)+rholm(lm,2)*dylmf2(ip,lm)
  END DO
  c=ylm(ip,3)*SQRT(4D0*pi/3D0)
  s=SQRT(1D0-c**2)
  DO  ispin=1,2
    dl(ispin)=dl(ispin)/rv**2
    d1(2,ispin)=d1(2,ispin)/rv
    d2(2,ispin)=d2(2,ispin)/rv
    d2(4,ispin)=d2(4,ispin)/rv**2
    IF(s > 1D-8)THEN
      d1(3,ispin)=d1(3,ispin)/rv/s
      d2(3,ispin)=d2(3,ispin)/rv/s
      d2(5,ispin)=d2(5,ispin)/rv**2/s
    ELSE
      d1(3,ispin)=0D0
      d2(3,ispin)=0D0
      d2(5,ispin)=0D0
    END IF
  END DO
  CALL fpexcpbe(d,dl,d1,d2,rv,s,c,vxcp(ip,1),vxcp(ip,2),excp(ip),um,bet)
END DO

RETURN
END SUBROUTINE mkxcpe_pbe

SUBROUTINE fpexcpbe(ro,rol,ro1,ro2,xr,s,c,v1,v2,exc,um,bet)
!----------------------------------------------------------------------
!     driver routine for PBE GGA subroutines.
!     based on excpbe.f in Munich code (version on 20 Dec 2009)
!     coded by M. Ogura, Jun. 2011, Munich
!----------------------------------------------------------------------

IMPLICIT NONE

REAL*8, INTENT(IN)                       :: ro(2)
REAL*8, INTENT(IN)                       :: rol(2)
REAL*8, INTENT(IN)                       :: ro1(3,2)
REAL*8, INTENT(IN)                       :: ro2(5,2)
REAL*8, INTENT(IN)                       :: xr
REAL*8, INTENT(IN)                       :: s
REAL*8, INTENT(IN)                       :: c
REAL*8, INTENT(OUT)                      :: v1
REAL*8, INTENT(OUT)                      :: v2
REAL*8, INTENT(OUT)                      :: exc
REAL*8, INTENT(IN)                       :: um
REAL*8, INTENT(IN)                       :: bet




REAL*8 conf,conrs,d,drv1,drv2,drv2s,drv3,drv4,ec,ex,fk,g,pi,rs,sk  &
    ,ss,thrd,thrd2,tt,uu,vcdn,vcup,vv,vx,vxcdn,vxcup,ww,x,xd,xu ,y,z,zet
INTEGER :: jsp,llda

pi=4D0*ATAN(1D0)
thrd=1D0/3D0
thrd2=2D0/3D0
conf=(3D0*pi**2)**thrd
conrs=(3D0/(4D0*pi))**thrd
llda=0
exc=0D0
vxcup=0D0
vxcdn=0D0
IF(ro(1) > 1D-12 .AND. ro(2) > 1D-12 )THEN
  
!     ---begin the spin loop for exchange
  IF(ro(1) <= 1D-6 .OR. ro(2) <= 1D-6)llda=1
  DO  jsp=1,2
    d=2D0*ro(jsp)
    fk=conf*d**thrd
    x=ro1(1,jsp)
    y=ro1(2,jsp)
    z=ro1(3,jsp)
    drv1=SQRT(x**2+y**2+z**2)*2D0
    IF(ABS(drv1) < 1D-8)THEN
      drv2=0D0
    ELSE
      drv2s=2D0*y*z*ro2(4,jsp)+(z**2-y**2)*ro2(5,jsp)-c/xr*y*(z**2+y**2)
      drv2=x**2*ro2(1,jsp)+2D0*x*y*ro2(2,jsp)+2D0*x*z*ro2(3,jsp)  &
          -x*y**2/xr-x*z**2/xr-y**2*rol(jsp)
      IF(ABS(drv2s) >= 1D-10)drv2=drv2+drv2s/s
      drv2=drv2/drv1*8D0
    END IF
    drv3=ro2(1,jsp)+2D0/xr*x-rol(jsp)
    drv3=drv3*2D0
    ss=drv1/(d*2D0*fk)
    uu=drv2/(d**2*(2D0*fk)**3)
    vv=drv3/(d*(2D0*fk)**2)
    CALL excpbex(d,ss,uu,vv,ex,vx,llda,um) !,bet
    exc=exc+ex*(d/2D0)/(ro(1)+ro(2))
    IF(jsp == 1)vxcup=vx
    IF(jsp == 2)vxcdn=vx
  END DO
  
!     ---correlation
  d=ro(1)+ro(2)
  zet=(ro(1)-ro(2))/d
  rs=conrs/d**thrd
  fk=1.91915829D0/rs
  sk=DSQRT(4D0*fk/pi)
  g=((1D0+zet)**thrd2+(1D0-zet)**thrd2)/2D0
  x=ro1(1,1)+ro1(1,2)
  y=ro1(2,1)+ro1(2,2)
  z=ro1(3,1)+ro1(3,2)
  drv1=SQRT(x**2+y**2+z**2)
  IF(drv1 < 1D-8)THEN
    drv2=0D0
  ELSE
    drv2s=2D0*y*z*(ro2(4,1)+ro2(4,2))+(z**2-y**2)*(ro2(5,1)+ro2(5,2))  &
        -c/xr*y*(z**2+y**2)
    drv2=x**2*(ro2(1,1)+ro2(1,2))+2D0*x*y*(ro2(2,1)+ro2(2,2))  &
        +2D0*x*z*(ro2(3,1)+ro2(3,2))-x*y**2/xr-x*z**2/xr -y**2*(rol(1)+rol(2))
    IF(ABS(drv2s) >= 1D-10)drv2=drv2+drv2s/s
    drv2=drv2/drv1
  END IF
  drv3=ro2(1,1)+ro2(1,2)+2D0/xr*x-rol(1)-rol(2)
  drv4=x*(ro1(1,1)-ro1(1,2)-zet*x)+y*(ro1(2,1)-ro1(2,2)-zet*y)  &
      +z*(ro1(3,1)-ro1(3,2)-zet*z)
  tt=drv1/(d*2D0*sk*g)
  uu=drv2/(d**2*(2D0*sk*g)**3)
  vv=drv3/(d*(2D0*sk*g)**2)
  ww=drv4/(d**2*(2D0*sk*g)**2)
  CALL excpbec(rs,zet,tt,uu,vv,ww,ec,vcup,vcdn,llda,bet) !,um
  exc=exc+ec
  vxcup=vxcup+vcup
  vxcdn=vxcdn+vcdn
  
!     ---convert from h to ry
  exc=2D0*exc
  xu=2D0*vxcup
  xd=2D0*vxcdn
  
  v1=xu
  v2=xd
END IF
END SUBROUTINE fpexcpbe

!*==excpbex.f    processed by SPAG 6.55Rc at 08:17 on 20 Dec 2009

SUBROUTINE excpbex(rho,s,u,v,ex,vx,llda,um) !,bet
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burke's modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
!   (for U,V, see PW86(24))
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!----------------------------------------------------------------------
! References:
! [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
!     Phys. Rev. Lett. 77, 3865 (1996).
! [b] J.P. Perdew and Y. Wang,
!     Phys. Rev. B33, 8800 (1986); B40, 3399 (1989) (E).
!----------------------------------------------------------------------
! Formulas:
!    e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
! e_x[PBE]=e_x[unif]*FxPBE(s)
! FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------

!  All input and output is in atomic units!

!  Modifications by: E. Engel
!  Last revision:    May 9, 2001
!engel

IMPLICIT NONE

REAL*8, INTENT(IN)                       :: rho
REAL*8, INTENT(IN)                       :: s
REAL*8, INTENT(IN OUT)                   :: u
REAL*8, INTENT(IN)                       :: v
REAL*8, INTENT(OUT)                      :: ex
REAL*8, INTENT(OUT)                      :: vx
REAL*8, INTENT(IN)                       :: um
!REAL*8, INTENT(IN)                       :: bet
INTEGER, INTENT(IN)                      :: llda

!*** Start of declarations rewritten by SPAG

! PARAMETER definitions

REAL*8, PARAMETER :: thrd=1.d0/3.d0
REAL*8, PARAMETER :: thrd4=4.d0/3.d0
REAL*8, PARAMETER :: ax=-0.738558766382022405884230032680836D0
REAL*8, PARAMETER :: uk=0.8040D0
REAL*8            :: ul


! Local variables

REAL*8 exunif,fs,fss,fxpbe,p0,s2

!*** End of declarations rewritten by SPAG

! Set value of 'ul'
ul=um/uk

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
exunif = ax*rho**thrd
IF ( llda == 1 ) THEN
  ex = exunif
  vx = ex*thrd4
  RETURN
END IF
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
s2 = s*s
p0 = 1.d0 + ul*s2
fxpbe = 1D0 + uk - uk/p0
ex = exunif*fxpbe
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
fs = 2.d0*uk*ul/(p0*p0)
fss = -4.d0*ul*s*fs/p0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24)
vx = exunif*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
END SUBROUTINE excpbex
!*==excpbec.f    processed by SPAG 6.55Rc at 08:17 on 20 Dec 2009

SUBROUTINE excpbec(rs,zeta,t,uu,vv,ww,ec,vcup,vcdn,llda,bet) ! removed um
!engel
!  This subroutine evaluates the correlation energy per particle and
!  spin-up and spin-dn correlation potentials within the Perdew-Burke-
!  Ernzerhof GGA. It is a slightly modified version of K. Burke's
!  official PBE subroutine.

!  Input:  RS   = WIGNER-SEITZ RADIUS = ( 3 / (4*PI*(DUP+DDN)) )**(1/3)
!          ZETA = RELATIVE SPIN POLARIZATION = (DUP-DDN)/(DUP+DDN)
!          T    = ABS(GRAD D) / ( (2*SK*G) * D )
!          UU   = (GRAD D)*GRAD(ABS(GRAD D)) / ( (2*SK*G)**3 * D**2 )
!          VV   = (LAPLACIAN D) / ( (2*SK*G)**2 * D )
!          WW   = (GRAD D)*(GRAD ZETA) / ( (2*SK*G)**2 * D )
!  where:  FK   = LOCAL FERMI MOMENTUM = (3*PI**2*(DUP+DDN))**(1/3)
!          SK   = LOCAL SCREENING MOMENTUM = (4*FK/PI)**(1/2)

!  Output: EC   = correlation energy per particle
!          VCUP = spin-up correlation potential
!          VCDN = spin-dn correlation potential

!  All input and output is in atomic units!

! References:
! [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
!     Phys. Rev. Lett. 77, 3865 (1996).
! [b] J. P. Perdew, K. Burke, and Y. Wang,
!     Phys. Rev. B54, 16533 (1996).
! [c] J. P. Perdew and Y. Wang,
!     Phys. Rev. B45, 13244 (1992).


!  Last revision:    May 9, 2001
!  Written by:       K. Burke, May 14, 1996.
!  Modifications by: E. Engel
!engel

IMPLICIT NONE

REAL*8, INTENT(IN)                       :: rs
REAL*8, INTENT(IN)                       :: zeta
REAL*8, INTENT(IN)                       :: t
REAL*8, INTENT(IN)                       :: uu
REAL*8, INTENT(IN)                       :: vv
REAL*8, INTENT(IN)                       :: ww
REAL*8, INTENT(OUT)                      :: ec
REAL*8, INTENT(OUT)                      :: vcup
REAL*8, INTENT(OUT)                      :: vcdn
!REAL*8, INTENT(IN)                       :: um
REAL*8, INTENT(IN)                       :: bet
INTEGER, INTENT(IN)                      :: llda

!*** Start of declarations rewritten by SPAG

! PARAMETER definitions

REAL*8, PARAMETER :: thrd=1.d0/3.d0
REAL*8, PARAMETER :: thrdm=-thrd
REAL*8, PARAMETER :: thrd2=2.d0*thrd
REAL*8, PARAMETER :: sixthm=thrdm/2.d0
REAL*8, PARAMETER :: thrd4=4.d0*thrd
REAL*8, PARAMETER :: gam=0.5198420997897463295344212145565D0
REAL*8, PARAMETER :: fzz=8.d0/(9.d0*gam)
REAL*8, PARAMETER :: gamma=0.03109069086965489503494086371273D0
REAL*8, PARAMETER :: eta=1.d-12


! Local variables

REAL*8 alfm,alfrsm,b,b2,bec,bg,comm,ecrs,eczeta,ep,eprs,eu,eurs,f,  &
    fac,fact0,fact1,fact2,fact3,fact5,fz,g,g3,g4,gz,h,hb,hbt,  &
    hrs,hrst,ht,htt,hz,hzt,pon,pref,q4,q5,q8,q9,rsthrd,rtrs,t2, &
    t4,t6,z4,delt

!*** End of declarations rewritten by SPAG

! Set value for 'delt'
delt=bet/gamma

! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
IF ( rs < 3.d5 ) THEN
  rtrs = SQRT(rs)
  CALL excgcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,  &
      0.49294D0,rtrs,eu,eurs)
  CALL excgcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,  &
      3.3662D0,0.62517D0,rtrs,ep,eprs)
  CALL excgcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,  &
      0.88026D0,0.49671D0,rtrs,alfm,alfrsm)
  z4 = zeta**4
  f = ((1.d0+zeta)**thrd4+(1.d0-zeta)**thrd4-2.d0)/gam
  ec = eu*(1.d0-f*z4) + ep*f*z4 - alfm*f*(1.d0-z4)/fzz
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZETA=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
  ecrs = eurs*(1.d0-f*z4) + eprs*f*z4 - alfrsm*f*(1.d0-z4)/fzz
  fz = thrd4*((1.d0+zeta)**thrd-(1.d0-zeta)**thrd)/gam
  eczeta = 4.d0*(zeta**3)*f*(ep-eu+alfm/fzz)  &
      + fz*(z4*ep-z4*eu-(1.d0-z4)*alfm/fzz)
  comm = ec - rs*ecrs/3.d0 - zeta*eczeta
  vcup = comm + eczeta
  vcdn = comm - eczeta
  IF ( llda == 1 ) RETURN
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
  g = ((1.d0+zeta)**thrd2+(1.d0-zeta)**thrd2)/2.d0
  g3 = g**3
  pon = -ec/(g3*gamma)
  b = delt/(EXP(pon)-1.d0)
  b2 = b*b
  t2 = t*t
  t4 = t2*t2
  q4 = 1.d0 + b*t2
  q5 = 1.d0 + b*t2 + b2*t4
  h = g3*(bet/delt)*LOG(1.d0+delt*q4*t2/q5)
  ec = ec + h
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
  g4 = g3*g
  t6 = t4*t2
  rsthrd = rs/3.d0
  gz = (((1.d0+zeta)**2+eta)**sixthm-((1.d0-zeta)**2+eta) **sixthm)/3.d0
  fac = delt/b + 1.d0
  bg = -3.d0*b2*ec*fac/(bet*g4)
  bec = b2*fac/(bet*g3)
  q8 = q5*q5 + delt*q4*q5*t2
  q9 = 1.d0 + 2.d0*b*t2
  hb = -bet*g3*b*t6*(2.d0+b*t2)/q8
  hrs = -rsthrd*hb*bec*ecrs
  fact0 = 2.d0*delt - 6.d0*b
  fact1 = q5*q9 + q4*q9*q9
  hbt = 2.d0*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
  hrst = rsthrd*t2*hbt*bec*ecrs
  hz = 3.d0*gz*h/g + hb*(bg*gz+bec*eczeta)
  ht = 2.d0*bet*g3*q9/q8
  hzt = 3.d0*gz*ht/g + hbt*(bg*gz+bec*eczeta)
  fact2 = q4*q5 + b*t2*(q4*q9+q5)
  fact3 = 2.d0*b*q5*q9 + delt*fact2
  htt = 4.d0*bet*g3*t*(2.d0*b/q8-(q9*fact3/q8)/q8)
  comm = h + hrs + hrst + t2*ht/6.d0 + 7.d0*t2*t*htt/6.d0
  pref = hz - gz*t2*ht/g
  fact5 = gz*(2.d0*ht+t*htt)/g
  comm = comm - pref*zeta - uu*htt - vv*ht - ww*(hzt-fact5)
  vcup = vcup + comm + pref
  vcdn = vcdn + comm - pref
ELSE
  vcup = 0.d0
  vcdn = 0.d0
END IF
END SUBROUTINE excpbec
!*==excgcor2.f    processed by SPAG 6.55Rc at 08:17 on 20 Dec 2009

SUBROUTINE excgcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.

IMPLICIT NONE

REAL*8, INTENT(IN)                       :: a
REAL*8, INTENT(IN)                       :: a1
REAL*8, INTENT(IN)                       :: b1
REAL*8, INTENT(IN)                       :: b2
REAL*8, INTENT(IN)                       :: b3
REAL*8, INTENT(IN)                       :: b4
REAL*8, INTENT(IN)                       :: rtrs
REAL*8, INTENT(OUT)                      :: gg
REAL*8, INTENT(OUT)                      :: ggrs

!*** Start of declarations rewritten by SPAG

! Local variables

REAL*8 q0,q1,q2,q3

!*** End of declarations rewritten by SPAG

q0 = -2.d0*a*(1.d0+a1*rtrs*rtrs)
q1 = 2.d0*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
q2 = LOG(1.d0+1.d0/q1)
gg = q0*q2
q3 = a*(b1/rtrs+2.d0*b2+rtrs*(3.d0*b3+4.d0*b4*rtrs))
ggrs = -2.d0*a*a1*q2 - q0*q3/(q1*(1.d0+q1))
END SUBROUTINE excgcor2


end module XCFunctionals_mod

