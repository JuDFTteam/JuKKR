SUBROUTINE cpw91(fk,sk,gz,ec,ecrs,eczta,rs,zta,t,uu,vv,ww,h,dvcup,  &
        dvcdn)
!-----------------------------------------------------------------
!gga91 correlation
!-----------------------------------------------------------------
!input
!      rs: seitz radius
!    zta: relative spin polarization
!       t: abs(grad d)/(d*2.*ks*gz)
!      uu: (grad d)*grad(abs(grad d))/(d**2 * (2*ks*gz)**3)
!      vv: (laplacian d)/(d * (2*ks*gz)**2)
!      ww: (grad d)*(gradzta)/(d * (2*ks*gz)**2
! output
!             h: nonlocal part of correlation energy per electron
!     dvcup,-dn: nonlocal parts of correlation potentials.

! with ks=sqrt(4*kf/pai), gz=[(1+zta)**(2/3)+(1-zta)**(2/3)]/2, &
!      kf=cbrt(3*pai**2*d).
!-----------------------------------------------------------------
!.. Scalar Arguments ..
DOUBLE PRECISION DVCDN,DVCUP,EC,ECRS,ECZTA,FK,GZ,H,RS,SK,T,UU,VV, &
                 WW,ZTA
!..
!.. Local Scalars ..
DOUBLE PRECISION A4,ALF,ARGMX,B,B2,B2FAC,BEC,BET,BG,C1,C2,C3,C4, &
                 C5,C6,CC,CC0,CCRS,COEFF,COMM,CX,DELT,FACT0,FACT1, &
                 FACT2,FACT3,FACT4,FACT5,GM,GZ3,GZ4,H0,H0B,H0BT, &
                 H0RS,H0RST,H0T,H0TT,H0Z,H0ZT,H1,H1RS,H1RST,H1T, &
                 H1TT,H1Z,H1ZT,HRS,HRST,HT,HTT,HZ,HZT,PON,PREF,Q4, &
                 Q5,Q6,Q7,Q8,Q9,R0,R1,R2,R3,R4,RS2,RS3,RSTHRD,T2, &
                 T4,T6,THRD2,THRDM,XNU
!..
!.. Intrinsic Functions ..
      INTRINSIC EXP,LOG
!..
!.. Save statement ..
      SAVE XNU,CC0,CX,ALF,C1,C2,C3,C4,C5,C6,A4,THRDM,THRD2
!..
!.. Data statements ..
!-----------------------------------------------------------------
      DATA XNU,CC0,CX,ALF/15.75592d0,0.004235d0,-0.001667212d0,0.09d0/
      DATA C1,C2,C3,C4/0.002568d0,0.023266d0,7.389d-6,8.723d0/
      DATA C5,C6,A4/0.472d0,7.389d-2,100.d0/
      DATA THRDM,THRD2/-0.333333333333d0,0.666666666667d0/
!..
!-----------------------------------------------------------------
argmx = 174.0D0
bet = xnu*cc0
delt = 2.d0*alf/bet
gz3 = gz**3
gz4 = gz3*gz
pon = -delt*ec/ (gz3*bet)
IF (pon > argmx) THEN
  b = 0.d0
ELSE
  b = delt/ (EXP(pon)-1.d0)
END IF
b2 = b*b
t2 = t*t
t4 = t2*t2
t6 = t4*t2
rs2 = rs*rs
rs3 = rs2*rs
q4 = 1.d0 + b*t2
q5 = 1.d0 + b*t2 + b2*t4
q6 = c1 + c2*rs + c3*rs2
q7 = 1.d0 + c4*rs + c5*rs2 + c6*rs3
cc = -cx + q6/q7
r0 = (sk/fk)**2
r1 = a4*r0*gz4
coeff = cc - cc0 - 3.d0*cx/7.d0
r2 = xnu*coeff*gz3
r3 = EXP(-r1*t2)
h0 = gz3* (bet/delt)*LOG(1.d0+delt*q4*t2/q5)
h1 = r3*r2*t2
h = h0 + h1
!  local correlation option:
!     h = 0.0d0

!  energy done. now the potential:

ccrs = (c2+2.d0*c3*rs)/q7 - q6* (c4+2.d0*c5*rs+3.d0*c6*rs2)/q7**2
rsthrd = rs/3.d0
r4 = rsthrd*ccrs/coeff
gm = ((1.d0+zta)**thrdm- (1.d0-zta)**thrdm)/3.d0
IF (pon > argmx) THEN
  b2fac = 0.d0
ELSE
  b2fac = b2* (delt/b+1.d0)
END IF
bg = -3.d0*ec*b2fac/ (bet*gz4)
bec = b2fac/ (bet*gz3)
q8 = q5*q5 + delt*q4*q5*t2
q9 = 1.d0 + 2.d0*b*t2
h0b = -bet*gz3*b*t6* (2.d0+b*t2)/q8
h0rs = -rsthrd*h0b*bec*ecrs
fact0 = 2.d0*delt - 6.d0*b
fact1 = q5*q9 + q4*q9*q9
h0bt = 2.d0*bet*gz3*t4* ((q4*q5*fact0-delt*fact1)/q8)/q8
h0rst = rsthrd*t2*h0bt*bec*ecrs
h0z = 3.d0*gm*h0/gz + h0b* (bg*gm+bec*eczta)
h0t = 2.d0*bet*gz3*q9/q8
h0zt = 3.d0*gm*h0t/gz + h0bt* (bg*gm+bec*eczta)
fact2 = q4*q5 + b*t2* (q4*q9+q5)
fact3 = 2.d0*b*q5*q9 + delt*fact2
h0tt = 4.d0*bet*gz3*t* (2.d0*b/q8- (q9* (fact3/q8))/q8)
h1rs = r3*r2*t2* (-r4+r1*t2/3.d0)
fact4 = 2.d0 - r1*t2
h1rst = r3*r2*t2* (2.d0*r4* (1.d0-r1*t2)-thrd2*r1*t2*fact4)
h1z = gm*r3*r2*t2* (3.d0-4.d0*r1*t2)/gz
h1t = 2.d0*r3*r2* (1.d0-r1*t2)
h1zt = 2.d0*gm*r3*r2* (3.d0-11.d0*r1*t2+4.d0*r1*r1*t4)/gz
h1tt = 4.d0*r3*r2*r1*t* (-2.d0+r1*t2)
hrs = h0rs + h1rs
hrst = h0rst + h1rst
ht = h0t + h1t
htt = h0tt + h1tt
hz = h0z + h1z
hzt = h0zt + h1zt
comm = h + hrs + hrst + t2*ht/6.d0 + 7.d0*t2*t*htt/6.d0
pref = hz - gm*t2*ht/gz
fact5 = gm* (2.d0*ht+t*htt)/gz
comm = comm - pref*zta - uu*htt - vv*ht - ww* (hzt-fact5)
dvcup = comm + pref
dvcdn = comm - pref

!  local correlation option:
!     dvcup = 0.0d0
!     dvcdn = 0.0d0

RETURN
END SUBROUTINE cpw91
