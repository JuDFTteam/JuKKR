SUBROUTINE mkxcpe2(ir,np,rv,rholm,vxcp,excp,ylm,dylmt1,dylmf1,  &
        dylmf2,dylmtf,drrl,ddrrl,drrul,ddrrul,irmd,  &
        lmpotd,lmmax,use_sol)
!     ------------------------------------------------------------------
!     Calculation of the exchange-correlation potential.
!     coded by M. Ogura, Apr. 2015, Munich
!     ------------------------------------------------------------------
implicit none
integer ijd
parameter(ijd=434)

double precision rv
integer ir,irmd,lmmax,lmpotd,np
double precision ddrrl(irmd,lmpotd),ddrrul(irmd,lmpotd), &
       drrl(irmd,lmpotd) &
      ,drrul(irmd,lmpotd),dylmf1(ijd,lmpotd),dylmf2(ijd,lmpotd) &
      ,dylmt1(ijd,lmpotd),dylmtf(ijd,lmpotd),excp(ijd) &
      ,rholm(lmpotd,2),vxcp(ijd,2),ylm(ijd,lmpotd)

double precision c,pi,s
integer ip,ispin,l1,lm,n
double precision d(2),d1(3,2),d2(5,2),dl(2)
!use_sol=0 -> PBE, use_sol=1 -> PBEsol
logical use_sol
double precision :: um,bet

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
      CALL fpexcpbe(d,dl,d1,d2,rv,s,c,vxcp(ip,1),vxcp(ip,2),excp(ip), um,bet)
    ELSE
      d1(3,ispin)=0D0
      d2(3,ispin)=0D0
      d2(5,ispin)=0D0
      vxcp(ip,1)=0D0
      vxcp(ip,2)=0D0
      excp(ip)=0D0
    END IF
  END DO
END DO

RETURN
END SUBROUTINE mkxcpe2

SUBROUTINE fpexcpbe(ro,rol,ro1,ro2,xr,s,c,v1,v2,exc,um,bet)
!----------------------------------------------------------------------
!     driver routine for PBE GGA subroutines.
!     based on excpbe.f in Munich code (version on 20 Dec 2009)
!     coded by M. Ogura, Jun. 2011, Munich
!----------------------------------------------------------------------
implicit none

double precision c,exc,s,v1,v2,xr
double precision ro(2),rol(2),ro1(3,2),ro2(5,2)
double precision um,bet

double precision conf,conrs,d,drv1,drv2,drv2s,drv3,drv4,ec,ex,fk, &
       g,pi,rs,sk &
      ,ss,thrd,thrd2,tt,uu,vcdn,vcup,vv,vx,vxcdn,vxcup,ww,x,xd,xu &
      ,y,z,zet
integer jsp,llda

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
    CALL excpbex(d,ss,uu,vv,ex,vx,llda,um)
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
  CALL excpbec(rs,zet,tt,uu,vv,ww,ec,vcup,vcdn,llda,bet)
  exc=exc+ec
  vxcup=vxcup+vcup
  vxcdn=vxcdn+vcdn
  
END IF
!     ---convert from h to ry
exc=2D0*exc
xu=2D0*vxcup
xd=2D0*vxcdn

v1=xu
v2=xd
END SUBROUTINE fpexcpbe

SUBROUTINE excpbex(rho,s,u,v,ex,vx,llda,um)
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
!
!  All input and output is in atomic units!
!
!  Modifications by: E. Engel
!  Last revision:    May 9, 2001
!engel
IMPLICIT NONE

! PARAMETER definitions
REAL*8, PARAMETER :: thrd=1.d0/3.d0
REAL*8, PARAMETER :: thrd4=4.d0/3.d0
REAL*8, PARAMETER :: ax=-0.738558766382022405884230032680836D0
REAL*8, PARAMETER :: uk=0.8040D0
DOUBLE PRECISION            :: ul

! Dummy arguments
DOUBLE PRECISION EX,RHO,S,U,V,VX,UM
INTEGER LLDA

! Local variables
DOUBLE PRECISION EXUNIF,FS,FSS,FXPBE,P0,S2

!----------------------------------------------------------------------
!    Define UL with via UM and UK
ul=um/uk

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

SUBROUTINE excpbec(rs,zeta,t,uu,vv,ww,ec,vcup,vcdn,llda,bet)
!engel
!  This subroutine evaluates the correlation energy per particle and
!  spin-up and spin-dn correlation potentials within the Perdew-Burke-
!  Ernzerhof GGA. It is a slightly modified version of K. Burke's
!  official PBE subroutine.
!
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
!
!  All input and output is in atomic units!
!
! References:
! [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
!     Phys. Rev. Lett. 77, 3865 (1996).
! [b] J. P. Perdew, K. Burke, and Y. Wang,
!     Phys. Rev. B54, 16533 (1996).
! [c] J. P. Perdew and Y. Wang,
!     Phys. Rev. B45, 13244 (1992).
!
!
!  Last revision:    May 9, 2001
!  Written by:       K. Burke, May 14, 1996.
!  Modifications by: E. Engel
!engel
IMPLICIT NONE

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

! Dummy arguments
DOUBLE PRECISION EC,RS,T,UU,VCDN,VCUP,VV,WW,ZETA,BET
INTEGER LLDA

! Local variables
DOUBLE PRECISION ALFM,ALFRSM,B,B2,BEC,BG,COMM,ECRS,ECZETA,EP, &
       EPRS,EU,EURS,F, &
       FAC,FACT0,FACT1,FACT2,FACT3,FACT5,FZ,G,G3,G4,GZ,H,HB,HBT, &
       HRS,HRST,HT,HTT,HZ,HZT,PON,PREF,Q4,Q5,Q8,Q9,RSTHRD,RTRS,T2, &
       T4,T6,Z4,DELT
EXTERNAL EXCGCOR2

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
!    Define DELT via BET and GAMMA
delt=bet/gamma

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

! Dummy arguments
DOUBLE PRECISION A,A1,B1,B2,B3,B4,GG,GGRS,RTRS

! Local variables
DOUBLE PRECISION Q0,Q1,Q2,Q3

q0 = -2.d0*a*(1.d0+a1*rtrs*rtrs)
q1 = 2.d0*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
q2 = LOG(1.d0+1.d0/q1)
gg = q0*q2
q3 = a*(b1/rtrs+2.d0*b2+rtrs*(3.d0*b3+4.d0*b4*rtrs))
ggrs = -2.d0*a*a1*q2 - q0*q3/(q1*(1.d0+q1))
END SUBROUTINE excgcor2
