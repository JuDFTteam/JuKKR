SUBROUTINE cradwf(eryd,ek,nsra,alpha,ipan,ircut,cvlight,rs,sl,  &
        pz,fz,qz,sz,tmat,vm2z,drdi,rmesh,zat,lirrsol,  &
        idoldau,lopt,wldauav,cutoff)
!-----------------------------------------------------------------------
!  subroutine for radial wave functions of spherical potentials

!             the generalized phase shifts are calculated by
!             a wronski relation :

!                 alpha(z,l) =-sqrt(z)*wronski{hl(r;z),rl(r;z)}; r->0

!             where hl is the free hankel function and rl the regular
!             solution . Using the analytical behaviour of rl at the
!             origin (rl = alphal * r**(l+1)  ; r->0),
!             the generalized phase shifts can be calculated
!             directly with the renormalization alphal .
!                                           b.drittler nov.1987

!   LDA+U added, March 2003 - Dec 2004, Munich/Juelich
!-----------------------------------------------------------------------
IMPLICIT NONE
!     .. Parameters ..
INCLUDE 'inc.p'
      INTEGER LMAXP1
      PARAMETER (LMAXP1=LMAXD+1)
      DOUBLE COMPLEX CI,CZERO
      PARAMETER (CI= (0.D0,1.D0),CZERO= (0.0D0,0.0D0))
!..
!.. Scalar Arguments ..
      DOUBLE COMPLEX ERYD,EK
      DOUBLE PRECISION CVLIGHT,ZAT
      DOUBLE PRECISION WLDAUAV
      INTEGER IPAN,NSRA,IDOLDAU,LOPT
      LOGICAL LIRRSOL
!..
!.. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),FZ(IRMD,0:LMAXD),PZ(IRMD,0:LMAXD), &
                     QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),TMAT(0:LMAXD)
      DOUBLE PRECISION DRDI(IRMD),RMESH(IRMD),RS(IRMD,0:LMAXD), &
                       SL(0:LMAXD),VM2Z(IRMD),CUTOFF(IRMD)
      INTEGER IRCUT(0:IPAND)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX ALPHAL,ARG,BL,EKLFAC,HL,PN,QF,SLOPE,TLSQEZ,VALUE
      DOUBLE PRECISION RIRC,RIRC1,RSIRC,S1
      INTEGER I,IR,IRC1,L1
!..
!.. Local Arrays ..
      DOUBLE COMPLEX BESSJW(0:LMAXP1),BESSYW(0:LMAXP1),DLOGDP(0:LMAXD), &
                     HAMF(IRMD,0:LMAXD),HANKWS(0:LMAXP1),MASS(IRMD)
      DOUBLE PRECISION DROR(IRMD)
!..
!.. External Subroutines ..
      EXTERNAL BESHAN,IRWSOL,REGSOL
!..
!.. Intrinsic Functions ..
      INTRINSIC DBLE

irc1 = ircut(ipan)
DO ir = 2,irc1
  dror(ir) = drdi(ir)/rmesh(ir)
END DO
rirc = rmesh(irc1)
rirc1 = 1D0/rirc
arg = rirc*ek
CALL beshan(hankws,bessjw,bessyw,arg,lmaxp1)

!---> calculate regular wavefunctions

CALL regsol(cvlight,eryd,nsra,dlogdp,fz,hamf,mass,pz,dror,rmesh,  &
    sl,vm2z,zat,ipan,ircut,idoldau,lopt,wldauav,cutoff, irmd,ipand,lmaxd)

eklfac = ek
! ======================================================================
DO l1 = 0,lmaxd
  
!---> determine t - matrix
  
  qf = DBLE(l1)*rirc1
  hl = hankws(l1) * dlogdp(l1)
  bl = bessjw(l1) * dlogdp(l1)
  hl = qf*hankws(l1) - ek*hankws(l1+1) - hl
  bl = bl - qf*bessjw(l1) + ek*bessjw(l1+1)
  hl = hl * ek
  tmat(l1) = ci * bl/hl
  
!---> determine the renormalization
  
  tlsqez = tmat(l1) * ek
  s1 = sl(l1)
  rsirc = rs(irc1,l1)
  eklfac = eklfac/ek*DBLE(2*l1+1)
  pn = pz(irc1,l1)*rsirc
  alphal = (bessjw(l1) - ci*hankws(l1)*tlsqez)*rirc/pn
  
!---> determine the alpha matrix
  
  alpha(l1) = alphal*eklfac
  
  DO i = 2,irc1
    pz(i,l1) = pz(i,l1)*alphal
    fz(i,l1) = fz(i,l1)*alphal
  END DO
  
  value = -ci*hankws(l1)*rirc*rsirc
  slope = DBLE(l1+1)*hankws(l1) - rirc*ek*hankws(l1+1)
  slope = (-ci*slope*rsirc+s1/rirc*value)
  qz(irc1,l1) = value
  sz(irc1,l1) = (slope*rirc - (s1+1.0D0)*value)/mass(irc1) * dror(irc1)
END DO
! ======================================================================

! -> calculate irregular wavefunctions

IF ( lirrsol ) CALL irwsol(ek,fz,hamf,mass,pz,qz,sz,dror,sl,  &
    ipan,ircut,irmd,ipand,lmaxd)
! ======================================================================
DO l1 = 0,lmaxd
  IF (nsra == 2) THEN
    DO i = 2,irc1
      pz(i,l1) = pz(i,l1)*rs(i,l1)
      qz(i,l1) = qz(i,l1)/rs(i,l1)
      fz(i,l1) = fz(i,l1)*rs(i,l1)/cvlight
      sz(i,l1) = sz(i,l1)/rs(i,l1)/cvlight
    END DO
  ELSE
    DO i = 2,irc1
      pz(i,l1) = pz(i,l1)*rs(i,l1)
      qz(i,l1) = qz(i,l1)/rs(i,l1)
      fz(i,l1) = czero
      sz(i,l1) = czero
    END DO
  END IF
END DO
! ======================================================================
END SUBROUTINE cradwf
