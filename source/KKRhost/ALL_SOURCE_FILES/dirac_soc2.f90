SUBROUTINE dirabmsoc2(getirrsol,c,socscl,it,e,l,mj,kap1,kap2,pis,  &
    cg1,cg2,cg4,cg5,cg8,v,b,z,nucleus,r,drdi,  &
    dovr,nmesh,dxp,pr,qr,pi,qi,dp,dq,nrmax)
!   ********************************************************************
!   *                                                                  *
!   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
!   *                                                                  *
!   *   with the SOC-operator manipulated                              *
!   *                                                                  *
!   *           H_soc = xi [ (sig_x*l_x+sig*l_y) + sig_z*l_z ]         *
!   *                                                                  *
!   *   SOCSCL = -1  ==  IXY = 0 : neglect 1st term of SOC-operator    *
!   *   SOCSCL = -2: ==  IXY = 1 : neglect 2nd term of SOC-operator    *
!   *                                                                  *
!   *   the outward integration is started by a power expansion        *
!   *   and continued by ADAMS-BASHFORTH-MOULTON - pred./corr.-method  *
!   *   NABM = 4(5) selects the 4(5)-point formula                     *
!   *                                                                  *
!   *   the inward integration is started analytically                 *
!   *                                                                  *
!   *   returns the wave functions up to the mesh point NMESH          *
!   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
!   *   and    R/I standing for regular/irregular solution             *
!   *                                                                  *
!   *  12/12/97  HE                                                    *
!   *  21/01/98  HE  finite nucelus                                    *
!   ********************************************************************

use mod_types, only: t_inc
IMPLICIT NONE

! PARAMETER definitions

INTEGER MPSMAX,NPEMAX,NABM
PARAMETER (MPSMAX=40,NPEMAX=4)
PARAMETER ( NABM   =      4 )
!PARAMETER ( NABM   =      5 )
DOUBLE COMPLEX CZ
PARAMETER (CZ=(0.0D0,0.0D0))
DOUBLE PRECISION TOL
PARAMETER (TOL=1.0D-9)
INTEGER ITMAX
PARAMETER (ITMAX=50)

! Dummy arguments
DOUBLE PRECISION C,CG1,CG2,CG4,CG5,CG8,MJ,SOCSCL
DOUBLE COMPLEX E
LOGICAL GETIRRSOL
INTEGER IT,KAP1,KAP2,L,NMESH,NRMAX,NUCLEUS,Z
DOUBLE COMPLEX PIS
DOUBLE PRECISION B(NRMAX),DOVR(NRMAX),DRDI(NRMAX),R(NRMAX), &
                 V(NRMAX)
DOUBLE COMPLEX DP(2,2,NRMAX),DQ(2,2,NRMAX),DXP(2,2), &
           PI(2,2,NRMAX), &
           PR(2,2,NRMAX),QI(2,2,NRMAX),QR(2,2,NRMAX)

! Local variables
DOUBLE COMPLEX AA11,AA12,AA21,AA22,ARG,BB1,BB2,BPP,BQQ,CFAC,D14, &
           DB14,DBH,DH,DIFFA,DIFFB,DV14,DVH,EMVPP,EMVQQ,S0,T0
DOUBLE PRECISION ACORR(0:NABM-1),ACORR0(0:NABM-1), &
       APRED(NABM),APRED0(NABM), &
       ASTEP,B14,BC(0:NPEMAX),BH,BHLP(NABM+4),CGD(2),CGMD(2),CGO, &
       CGOZ,CGZ(2),CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),CSQR, &
       DBDR(NRMAX),DBHLP(NABM+4),DHLP(NABM+4),DVDR(NRMAX), &
       DVHLP(NABM+4),GAM(2),GPM,HLP(NABM+4),HLP1,KAP(2),KPX(2), &
       KPY(2),R14,RH,RHLP(NABM+4),RPWGPM,RR,SK(2),SK1,SK2,SO2,SO6, &
       SRK,TZ,V14,VC(0:NPEMAX),VH,VHLP(NABM+4),X14,XH
DOUBLE COMPLEX CJLZ
DOUBLE PRECISION DABS,DBLE,DSQRT
DOUBLE COMPLEX DETD,MP1(2,2),MP2(2,2),MP3(2,2),MP4(2,2),MQ1(2,2), &
           MQ2(2,2),MQ3(2,2),MQ4(2,2),P1(2,2),P2(2,2),P3(2,2), &
           P4(2,2),PC(2,2,-NPEMAX:MPSMAX),PNEW(2,2),POLD(2,2), &
           Q1(2,2),Q2(2,2),Q3(2,2),Q4(2,2),QC(2,2,-NPEMAX:MPSMAX), &
           QNEW(2,2),QOLD(2,2),SOCPP(2),ZZ
INTEGER I,IC,IP,IRK,ISK1,ISK2,IV,IXY,J,JCORR,K,LB(2),LB1,LB2,M, &
        MPS,N,NACORR,NDIV,NHLP,NM,NPE,NSOL,NTOP
INTEGER INT,ISIGN,NINT
DOUBLE PRECISION YLAG

DATA APRED0/55.0D0, - 59.0D0, + 37.0D0, - 9.0D0/
DATA ACORR0/9.0D0, + 19.0D0, - 5.0D0, + 1.0D0/
DATA ASTEP/24.0D0/

csqr = c*c
cfac = pis*c/(e+csqr)

! find   NPE  expansion coefficients for the potential and b-field
npe = 4

tz = DBLE(2*z)

DO iv = 1,npe
  DO n = 1,npe
    cm(n,iv) = r(n)**(iv-1)
  END DO
END DO

CALL rinvgj(cmi,cm,npemax,npe)

DO iv = 1,npe
  vc(iv-1) = 0.0D0
  DO n = 1,npe
    vc(iv-1) = vc(iv-1) + cmi(iv,n)*(v(n)+tz/r(n))
  END DO
END DO

DO iv = 1,npe
  bc(iv-1) = 0.0D0
  DO n = 1,npe
    bc(iv-1) = bc(iv-1) + cmi(iv,n)*b(n)
  END DO
END DO

!*HF* IXY=0 ==> ONLY ZZ-PART, IXY=1 ==> ONLY XY-PART
IF ( nint(socscl) == -1 ) THEN
  ixy = 0
ELSE
  ixy = 1
END IF
!*HF*

!    calculate g-coefficients of b-field

isk1 = ISIGN(1,kap1)
isk2 = ISIGN(1,kap2)
sk1 = DBLE(isk1)
sk2 = DBLE(isk2)
lb1 = l - isk1
lb2 = l - isk2

cg1 = -mj/(kap1+0.5D0)
cg5 = -mj/(-kap1+0.5D0)
cgd(1) = cg1
cgmd(1) = cg5
kap(1) = DBLE(kap1)
gam(1) = DSQRT(kap(1)**2-(tz/c)**2)
IF ( ixy == 0 ) THEN
  cgz(1) = kap(1) + 1 + mj*cgd(1) - 0.5D0
ELSE
  cgz(1) = -(mj*cgd(1)-0.5D0)
END IF
lb(1) = lb1
sk(1) = sk1
IF ( DABS(mj) > l ) THEN
  cg2 = 0.0D0
  cg4 = 0.0D0
  cg8 = 0.0D0
  nsol = 1
  cgd(2) = 0.0D0
  cgo = 0.0D0
  cgoz = 0.0D0
  cgmd(2) = 0.0D0
  gam(2) = 0.0D0
  kap(2) = 0.0D0
  cgz(2) = 0.0D0
  lb(2) = 0
  sk(2) = 0.0D0
ELSE
  cg2 = -DSQRT(1.0D0-(mj/(kap1+0.5D0))**2)
  cg4 = -mj/(kap2+0.5D0)
  cg8 = -mj/(-kap2+0.5D0)
  nsol = 2
  cgd(2) = cg4
  cgo = cg2
  cgmd(2) = cg8
  kap(2) = DBLE(kap2)
  gam(2) = DSQRT(kap(2)**2-(tz/c)**2)
  IF ( ixy == 0 ) THEN
    cgz(2) = kap(2) + 1 + mj*cgd(2) - 0.5D0
    cgoz = mj*cgo
  ELSE
    cgz(2) = -(mj*cgd(2)-0.5D0)
    cgoz = -mj*cgo
  END IF
  lb(2) = lb2
  sk(2) = sk2
END IF
DO i = 1,nsol
  kpx(i) = kap(i)
  kpy(i) = kap(i)
END DO

CALL dvdrspline(v,r,dvdr,nmesh)
CALL derspl(nmesh,r,b,dbdr)

DO ip = 1,nabm
  ic = ip - 1
  apred(ip) = apred0(ip)/astep
  acorr(ic) = acorr0(ic)/astep
END DO
nacorr = nabm - 1

DO i = 1,2
  DO j = 1,2
    DO ip = -npemax,mpsmax
      pc(i,j,ip) = cz
      qc(i,j,ip) = cz
    END DO
  END DO
END DO

! ======================================================================

IF ( (tz >= 2) .AND. (nucleus == 0) ) THEN
  
  DO j = 1,nsol
    i = 3 - j
    pc(j,j,0) = DSQRT(ABS(kpy(j))-gam(j))
    qc(j,j,0) = (kpy(j)+gam(j))*(csqr/tz)*pc(j,j,0)
    pc(i,j,0) = cz
    qc(i,j,0) = cz
  END DO
  
!  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS
  
  mps = 40
  
  aa12 = -tz/csqr
  aa21 = tz
  emvqq = (e-vc(0)+csqr)/csqr
  emvpp = -e + vc(0)
  bqq = bc(0)/csqr
  DO i = 1,nsol
    socpp(i) = cz
  END DO
  
  DO j = 1,nsol
    
    DO m = 1,mps
      DO i = 1,nsol
        k = 3 - i
        bb1 = (emvqq+bqq*cgmd(i))*qc(i,j,m-1)
        bb2 = (emvpp+bc(0)*cgd(i))*pc(i,j,m-1) + bc(0) *cgo*pc(k,j,m-1)
        DO ip = 1,npe - 1
          bb1 = bb1 + (-vc(ip)+bc(ip)*cgmd(i))*qc(i,j,m-1-ip) /csqr
          bb2 = bb2 + (+vc(ip)+bc(ip)*cgd(i))*pc(i,j,m-1-ip)  &
              + bc(ip)*cgo*pc(k,j,m-1-ip)
        END DO
        
        aa11 = gam(j) + m + kpy(i)
        aa22 = gam(j) + m - kpy(i)
        detd = aa11*aa22 - aa12*aa21
        pc(i,j,m) = (bb1*aa22-aa12*bb2)/detd
        qc(i,j,m) = (aa11*bb2-bb1*aa21)/detd
      END DO
    END DO
    
  END DO
  
  
!  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
!  FOR THE FIRST   NABM   R - MESH - POINTS
  
  DO n = 1,nabm
    rr = r(n)
    
    DO j = 1,nsol
      rpwgpm = rr**gam(j)
      
      DO i = 1,nsol
        pr(i,j,n) = pc(i,j,0)*rpwgpm
        qr(i,j,n) = qc(i,j,0)*rpwgpm
        dp(i,j,n) = pc(i,j,0)*rpwgpm*gam(j)*dovr(n)
        dq(i,j,n) = qc(i,j,0)*rpwgpm*gam(j)*dovr(n)
      END DO
      
      DO m = 1,mps
        rpwgpm = rpwgpm*rr
        gpm = gam(j) + m
        
        DO i = 1,nsol
          pr(i,j,n) = pr(i,j,n) + pc(i,j,m)*rpwgpm
          qr(i,j,n) = qr(i,j,n) + qc(i,j,m)*rpwgpm
          dp(i,j,n) = dp(i,j,n) + pc(i,j,m) *rpwgpm*gpm*dovr(n)
          dq(i,j,n) = dq(i,j,n) + qc(i,j,m) *rpwgpm*gpm*dovr(n)
        END DO
        
      END DO
    END DO
  END DO
! ======================================================================
!                                  == EMPTY SPHERE  or FINITE NUCLEUS ==
ELSE
  
!        assume constant pot: V=V(1)   ignore coupling: B=0
  
  t0 = e - v(1)
  s0 = (e-v(1))/csqr + 1
  
  DO n = 1,nabm
    
    DO j = 1,nsol
      DO i = 1,nsol
        pr(i,j,n) = cz
        qr(i,j,n) = cz
        dp(i,j,n) = cz
        dq(i,j,n) = cz
      END DO
    END DO
    
    zz = CDSQRT(s0*t0)*r(n)
    
    DO j = 1,nsol
      pr(j,j,n) = cjlz(l,zz)*r(n)
      dp(j,j,n) = (DBLE(l+1)*cjlz(l,zz)-zz*cjlz(l+1,zz)) *drdi(n)
      
      qr(j,j,n) = (dp(j,j,n)/drdi(n)+pr(j,j,n)*(kap(j)/r(n))) /s0
      dq(j,j,n) = qr(j,j,n)*(kap(j)/r(n)) - pr(j,j,n)*t0
    END DO
  END DO
  
END IF
! ===================================================================

! =============================================================== N ====
!     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-BASHFORTH-MOULTON)

DO n = nabm + 1,nmesh
  
!    EVALUATE PREDICTOR
  
  DO j = 1,nsol
    DO i = 1,nsol
      pnew(i,j) = pr(i,j,n-1)
      qnew(i,j) = qr(i,j,n-1)
      
      DO ip = 1,nabm
        pnew(i,j) = pnew(i,j) + apred(ip)*dp(i,j,n-ip)
        qnew(i,j) = qnew(i,j) + apred(ip)*dq(i,j,n-ip)
      END DO
    END DO
  END DO
  
  emvqq = (e-v(n)+csqr)*drdi(n)/csqr
  emvpp = -(e-v(n))*drdi(n)
  bqq = b(n)*drdi(n)/csqr
  bpp = b(n)*drdi(n)
  DO i = 1,nsol
    socpp(i) = dovr(n)*csqr*(dvdr(n)-dbdr(n)*cgmd(i))  &
        /(e-v(n)+csqr+b(n)*cgmd(i))**2
    
  END DO
  
!    EVALUATE CORRECTOR
  
  
  DO jcorr = 1,itmax
    
    DO j = 1,nsol
      DO i = 1,nsol
        k = 3 - i
        pold(i,j) = pnew(i,j)
        qold(i,j) = qnew(i,j)
        dp(i,j,n) = -kpx(i)*pnew(i,j)*dovr(n) + (emvqq+bqq*cgmd(i))*qnew(i,j)
        dq(i,j,n) = kpx(i)*qnew(i,j)*dovr(n) + (emvpp+bpp*cgd(i))*pnew(i,j)  &
            + (bpp*cgo+socpp(k)*cgoz)*pnew(k,j) + socpp(i)*cgz(i)*pnew(i,j)
        pnew(i,j) = pr(i,j,n-1)
        qnew(i,j) = qr(i,j,n-1)
        DO ic = 0,nacorr
          pnew(i,j) = pnew(i,j) + acorr(ic)*dp(i,j,n-ic)
          qnew(i,j) = qnew(i,j) + acorr(ic)*dq(i,j,n-ic)
        END DO
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        diffa = pold(i,j) - pnew(i,j)
        IF ( ABS(dreal(diffa)) > (tol*ABS(dreal(pnew(i,j)))) ) GO TO 50
        IF ( ABS(DIMAG(diffa)) > (tol*ABS(DIMAG(pnew(i,j)))) ) GO TO 50
        
        diffb = qold(i,j) - qnew(i,j)
        IF ( ABS(dreal(diffb)) > (tol*ABS(dreal(qnew(i,j)))) ) GO TO 50
        IF ( ABS(DIMAG(diffb)) > (tol*ABS(DIMAG(qnew(i,j)))) ) GO TO 50
      END DO
    END DO
    GO TO 100
    
  50      END DO
  IF(t_inc%i_write>0) WRITE (1337,99001) kap1,n,r(n),diffa,  &
      diffb,it,l,INT(2*mj),'REG'
  
!                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
  
  100     CONTINUE
  DO j = 1,nsol
    DO i = 1,nsol
      k = 3 - i
      pr(i,j,n) = pnew(i,j)
      qr(i,j,n) = qnew(i,j)
      dp(i,j,n) = -kpx(i)*pnew(i,j)*dovr(n) + (emvqq+bqq*cgmd(i))*qnew(i,j)
      dq(i,j,n) = kpx(i)*qnew(i,j)*dovr(n) + (emvpp+bpp*cgd(i))  &
          *pnew(i,j) + (bpp*cgo+socpp(k)*cgoz)  &
          *pnew(k,j) + socpp(i)*cgz(i)*pnew(i,j)
    END DO
  END DO
  
  
END DO
! =============================================================== N ====

DO i = 1,2
  DO j = 1,2
    dxp(i,j) = dp(i,j,nmesh)
  END DO
END DO


CALL cinit(2*2*nrmax,pi)
CALL cinit(2*2*nrmax,qi)

IF ( .NOT.getirrsol ) RETURN

! #####################################################################

!             IRREGULAR SOLUTION IRREGULAR SOLUTION  IRREGULAR SOLUTION

!  CALCULATE THE INITIAL VALUES OF THE WAVEFUNCTION AT THE SPHERE
!  BOUNDARY

DO n = nmesh,nmesh + nabm
  arg = pis*r(n)
  
  DO j = 1,nsol
    i = 3 - j
    pi(j,j,n) = cjlz(l,arg)*r(n)
    qi(j,j,n) = cfac*sk(j)*cjlz(lb(j),arg)*r(n)*c
    dp(j,j,n) = (DBLE(l+1)*cjlz(l,arg)-arg*cjlz(l+1,arg)) *drdi(n)
    m = lb(j)
    dq(j,j,n) = cfac*sk(j) *(DBLE(m+1)*cjlz(m,arg)-arg*cjlz(m+1,arg))  &
        *drdi(n)*c
    
    pi(i,j,n) = cz
    qi(i,j,n) = cz
    dp(i,j,n) = cz
    dq(i,j,n) = cz
  END DO
END DO
!        ------------------------------------------------------------
!              INITIALIZE INWARD INTEGRATION WITH RUNGE - KUTTA
!        ------------------------------------------------------------
ndiv = 60
IF ( ndiv /= 0 ) THEN
  
  srk = 1.0D0/DBLE(ndiv)
  so2 = srk/2.0D0
  so6 = srk/6.0D0
  
  n = nmesh
  
  emvqq = (e-v(n)+csqr)*drdi(n)/csqr
  emvpp = -(e-v(n))*drdi(n)
  bqq = b(n)*drdi(n)/csqr
  bpp = b(n)*drdi(n)
  DO i = 1,nsol
    socpp(i) = dovr(n)*csqr*(dvdr(n)-dbdr(n)*cgmd(i))  &
        /(e-v(n)+csqr+b(n)*cgmd(i))**2
    
  END DO
  
! *** reinitialize Q using only DP and PI
  
  DO j = 1,nsol
    i = 3 - j
    qi(j,j,n) = (dp(j,j,n)+kpx(j)*pi(j,j,n)*dovr(n)) /(emvqq+bqq*cgmd(j))
    qi(i,j,n) = cz
  END DO
  
  DO j = 1,nsol
    DO i = 1,nsol
      k = 3 - i
      dq(i,j,n) = kpx(i)*qi(i,j,n)*dovr(n) + (emvpp+bpp*cgd(i))  &
          *pi(i,j,n) + (bpp*cgo+socpp(k)*cgoz)  &
          *pi(k,j,n) + socpp(i)*cgz(i)*pi(i,j,n)
    END DO
  END DO
  
  DO j = 1,nsol
    DO i = 1,nsol
      p1(i,j) = pi(i,j,n)
      q1(i,j) = qi(i,j,n)
      mp1(i,j) = dp(i,j,n)
      mq1(i,j) = dq(i,j,n)
    END DO
  END DO
  
  x14 = DBLE(n)
  nhlp = nabm + 4
  hlp1 = DBLE(nmesh-nhlp)
  DO i = 1,nhlp
    hlp(i) = DBLE(i)
    vhlp(i) = v(nmesh-nhlp+i)
    bhlp(i) = b(nmesh-nhlp+i)
    dhlp(i) = drdi(nmesh-nhlp+i)
    rhlp(i) = r(nmesh-nhlp+i)
    dvhlp(i) = dvdr(nmesh-nhlp+i)
    dbhlp(i) = dbdr(nmesh-nhlp+i)
  END DO
  
  DO irk = 1,(nabm-1)*ndiv
    
    xh = x14 - so2
    vh = ylag(xh-hlp1,hlp,vhlp,0,3,nhlp)
    bh = ylag(xh-hlp1,hlp,bhlp,0,3,nhlp)
    rh = ylag(xh-hlp1,hlp,rhlp,0,3,nhlp)
    dh = ylag(xh-hlp1,hlp,dhlp,0,3,nhlp)
    dvh = ylag(xh-hlp1,hlp,dvhlp,0,3,nhlp)
    dbh = ylag(xh-hlp1,hlp,dbhlp,0,3,nhlp)
    
    emvqq = (e-vh+csqr)*dh/csqr
    emvpp = -(e-vh)*dh
    bqq = bh*dh/csqr
    bpp = bh*dh
    DO i = 1,nsol
      socpp(i) = (dh/rh)*csqr*(dvh-dbh*cgmd(i)) /(e-vh+csqr+bh*cgmd(i))**2
      
    END DO
    n = nmesh - irk/ndiv - n + nint(xh)
    
    DO j = 1,nsol
      DO i = 1,nsol
        p2(i,j) = p1(i,j) - so2*mp1(i,j)
        q2(i,j) = q1(i,j) - so2*mq1(i,j)
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        k = 3 - i
        mp2(i,j) = -kpx(i)*p2(i,j)*dh/rh + (emvqq+bqq*cgmd(i)) *q2(i,j)
        mq2(i,j) = kpx(i)*q2(i,j)*dh/rh + (emvpp+bpp*cgd(i))  &
            *p2(i,j) + (bpp*cgo+socpp(k)*cgoz)*p2(k,j) + socpp(i)*cgz(i)*p2(i,j)
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        p3(i,j) = p1(i,j) - so2*mp2(i,j)
        q3(i,j) = q1(i,j) - so2*mq2(i,j)
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        k = 3 - i
        mp3(i,j) = -kpx(i)*p3(i,j)*dh/rh + (emvqq+bqq*cgmd(i)) *q3(i,j)
        mq3(i,j) = kpx(i)*q3(i,j)*dh/rh + (emvpp+bpp*cgd(i))  &
            *p3(i,j) + (bpp*cgo+socpp(k)*cgoz)*p3(k,j) + socpp(i)*cgz(i)*p3(i,j)
      END DO
    END DO
    
    x14 = x14 - srk
    v14 = ylag(x14-hlp1,hlp,vhlp,0,3,nhlp)
    b14 = ylag(x14-hlp1,hlp,bhlp,0,3,nhlp)
    r14 = ylag(x14-hlp1,hlp,rhlp,0,3,nhlp)
    d14 = ylag(x14-hlp1,hlp,dhlp,0,3,nhlp)
    dv14 = ylag(x14-hlp1,hlp,dvhlp,0,3,nhlp)
    db14 = ylag(x14-hlp1,hlp,dbhlp,0,3,nhlp)
    
    
    emvqq = (e-v14+csqr)*d14/csqr
    emvpp = -(e-v14)*d14
    bqq = b14*d14/csqr
    bpp = b14*d14
    DO i = 1,nsol
      socpp(i) = (d14/r14)*csqr*(dv14-db14*cgmd(i))  &
          /(e-v14+csqr+b14*cgmd(i))**2
      
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        p4(i,j) = p1(i,j) - srk*mp3(i,j)
        q4(i,j) = q1(i,j) - srk*mq3(i,j)
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        k = 3 - i
        mp4(i,j) = -kpx(i)*p4(i,j)*d14/r14 + (emvqq+bqq*cgmd(i))*q4(i,j)
        mq4(i,j) = kpx(i)*q4(i,j)*d14/r14 + (emvpp+bpp*cgd(i))  &
            *p4(i,j) + (bpp*cgo+socpp(k)*cgoz)*p4(k,j) + socpp(i)*cgz(i)*p4(i,j)
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        p1(i,j) = p1(i,j) - so6*(mp1(i,j)+2*(mp2(i,j)+mp3(i,j))  &
            +mp4(i,j))
        q1(i,j) = q1(i,j) - so6*(mq1(i,j)+2*(mq2(i,j)+mq3(i,j))  &
            +mq4(i,j))
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        k = 3 - i
        mp1(i,j) = -kpx(i)*p1(i,j)*d14/r14 + (emvqq+bqq*cgmd(i))*q1(i,j)
        mq1(i,j) = kpx(i)*q1(i,j)*d14/r14 + (emvpp+bpp*cgd(i))  &
            *p1(i,j) + (bpp*cgo+socpp(k)*cgoz)*p1(k,j) + socpp(i)*cgz(i)*p1(i,j)
      END DO
    END DO
    
    IF ( MOD(irk,ndiv) == 0 ) THEN
      n = nmesh - irk/ndiv
      IF ( ABS(x14-DBLE(n)) > 1.0D-5 ) THEN
        WRITE (*,*) ' <DIRAC> RUNGE-KUTTA: ',irk,ndiv,n,x14
        STOP
      END IF
      DO j = 1,nsol
        DO i = 1,nsol
          pi(i,j,n) = p1(i,j)
          qi(i,j,n) = q1(i,j)
          dp(i,j,n) = mp1(i,j)
          dq(i,j,n) = mq1(i,j)
        END DO
      END DO
      
    END IF
    
  END DO
  
END IF

! =============================================================== N ====

!     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-BASHFORTH-MOULTON)

IF ( ndiv /= 0 ) THEN
  ntop = nmesh - nabm
ELSE
  ntop = nmesh
END IF

DO nm = 1,ntop
  n = 1 + ntop - nm
  
!    EVALUATE PREDICTOR
  
  DO j = 1,nsol
    DO i = 1,nsol
      pnew(i,j) = pi(i,j,n+1)
      qnew(i,j) = qi(i,j,n+1)
      
      DO ip = 1,nabm
        pnew(i,j) = pnew(i,j) - apred(ip)*dp(i,j,n+ip)
        qnew(i,j) = qnew(i,j) - apred(ip)*dq(i,j,n+ip)
      END DO
    END DO
  END DO
  
  emvqq = (e-v(n)+csqr)*drdi(n)/csqr
  emvpp = -(e-v(n))*drdi(n)
  bqq = b(n)*drdi(n)/csqr
  bpp = b(n)*drdi(n)
  DO i = 1,nsol
    socpp(i) = dovr(n)*csqr*(dvdr(n)-dbdr(n)*cgmd(i))  &
        /(e-v(n)+csqr+b(n)*cgmd(i))**2
  END DO
  
!    EVALUATE CORRECTOR
  
  DO jcorr = 1,itmax
    DO j = 1,nsol
      DO i = 1,nsol
        k = 3 - i
        pold(i,j) = pnew(i,j)
        qold(i,j) = qnew(i,j)
        dp(i,j,n) = -kpx(i)*pnew(i,j)*dovr(n) + (emvqq+bqq*cgmd(i))*qnew(i,j)
        dq(i,j,n) = kpx(i)*qnew(i,j)*dovr(n) + (emvpp+bpp*cgd(i))*pnew(i,j)  &
            + (bpp*cgo+socpp(k)*cgoz)*pnew(k,j) + socpp(i)*cgz(i)*pnew(i,j)
        
        pnew(i,j) = pi(i,j,n+1)
        qnew(i,j) = qi(i,j,n+1)
        DO ic = 0,nacorr
          pnew(i,j) = pnew(i,j) - acorr(ic)*dp(i,j,n+ic)
          qnew(i,j) = qnew(i,j) - acorr(ic)*dq(i,j,n+ic)
        END DO
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        diffa = pold(i,j) - pnew(i,j)
        IF ( ABS(dreal(diffa)) > (tol*ABS(dreal(pnew(i,j)))) ) GO TO 150
        IF ( ABS(DIMAG(diffa)) > (tol*ABS(DIMAG(pnew(i,j)))) ) GO TO 150
        
        diffb = qold(i,j) - qnew(i,j)
        IF ( ABS(dreal(diffb)) > (tol*ABS(dreal(qnew(i,j)))) ) GO TO 150
        IF ( ABS(DIMAG(diffb)) > (tol*ABS(DIMAG(qnew(i,j)))) ) GO TO 150
      END DO
    END DO
    GO TO 200
    
  150     END DO
  IF(t_inc%i_write>0) WRITE (1337,99001) kap1,n,r(n),diffa,  &
      diffb,it,l,INT(2*mj),'IRR'
  
!                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
  
  200     CONTINUE
  DO j = 1,nsol
    DO i = 1,nsol
      k = 3 - i
      pi(i,j,n) = pnew(i,j)
      qi(i,j,n) = qnew(i,j)
      dp(i,j,n) = -kpx(i)*pnew(i,j)*dovr(n) + (emvqq+bqq*cgmd(i))*qnew(i,j)
      dq(i,j,n) = kpx(i)*qnew(i,j)*dovr(n) + (emvpp+bpp*cgd(i))  &
          *pnew(i,j) + (bpp*cgo+socpp(k)*cgoz)  &
          *pnew(k,j) + socpp(i)*cgz(i)*pnew(i,j)
    END DO
  END DO
  
END DO
! =============================================================== N ====

!     the minor component for the soc-manipulated wf is meaningless
!     =>  set it to zero

CALL cinit(2*2*nrmax,qr)
CALL cinit(2*2*nrmax,qi)
99001 FORMAT (' PRE/CORR NOT CONV. IN <DIRABMSOC2> ',2I4,f10.7,2X,  &
    4E12.4,3I2,'/2 ',a3)

END SUBROUTINE dirabmsoc2
!*==dvdrspline.f    processed by SPAG 6.05Rc at 13:35 on  4 Dec 2002

SUBROUTINE dvdrspline(v,r,dvdr,n)
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE DERIVATIVE OF THE POTENTIAL                          *
!   *   DV/DR= ( D(R**2*V)/DR-2.R*V ) / R**2                           *
!   *                                                                  *
!   *   calling J.REDINGER's routine <DERSPL>                          *
!   *                                                                  *
!   * 27/10/94  HE  adopted for SPRKKR-package                         *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER N
DOUBLE PRECISION DVDR(N),R(N),V(N)

! Local variables
DOUBLE PRECISION DR2VDR(N),R2V(N)
INTEGER I

DO i = 1,n
  r2v(i) = v(i)*r(i)**2
END DO

CALL derspl(n,r,r2v,dr2vdr)

DO i = 1,n
  dvdr(i) = (dr2vdr(i)-2.d0*r2v(i)/r(i))/r(i)**2
END DO
END SUBROUTINE dvdrspline
!*==derspl.f    processed by SPAG 6.05Rc at 13:35 on  4 Dec 2002

SUBROUTINE derspl(n,x,f,d)
!   ********************************************************************
!   *                                                                  *
!   *   F(I) ARE THE FUNCTION VALUES AT THE POINTS X(I) FOR I=1,N      *
!   *   AND THE SPLINE DERIVATIVES D(I) ARE FOUND.                     *
!   *   THE DIMENSION OF A MUST NOT BE LESS THAN 3*N.                  *
!   *                                                                  *
!   *   sepp redinger 1985                                             *
!   *                                                                  *
!   ********************************************************************
use mod_types, only: t_inc
IMPLICIT NONE

! Dummy arguments
INTEGER N
DOUBLE PRECISION D(N),F(N),X(N)

! Local variables
DOUBLE PRECISION A(N*3),H1,H2,P
INTEGER I,J,K

DO i = 2,n
  IF ( x(i) <= x(i-1) ) THEN
    IF(t_inc%i_write>0) WRITE (1337,99001) i
    a(1) = 1.d0
    RETURN
  END IF
END DO
DO i = 1,n
  j = 2
  IF ( i /= 1 ) THEN
    j = n - 1
    IF ( i /= n ) THEN
      h1 = 1.d0/(x(i)-x(i-1))
      h2 = 1.d0/(x(i+1)-x(i))
      a(3*i-2) = h1
      a(3*i-1) = 2.d0*(h1+h2)
      a(3*i) = h2
      d(i) = 3*(f(i+1)*h2*h2+f(i)*(h1*h1-h2*h2)-f(i-1)*h1*h1)
      GO TO 100
    END IF
  END IF
  h1 = 1.d0/(x(j)-x(j-1))
  h2 = 1.d0/(x(j+1)-x(j))
  a(3*i-2) = h1*h1
  a(3*i-1) = h1*h1 - h2*h2
  a(3*i) = -h2*h2
  d(i) = 2.d0*(f(j)*(h2*h2*h2+h1*h1*h1)-f(j+1)*h2*h2*h2-f(j-1) *h1*h1*h1)
100  END DO
p = a(4)/a(1)
a(5) = a(5) - p*a(2)
a(6) = a(6) - p*a(3)
d(2) = d(2) - p*d(1)
DO i = 3,n
  k = 3*i - 4
  p = a(k+2)/a(k)
  a(k+3) = a(k+3) - p*a(k+1)
  d(i) = d(i) - p*d(i-1)
  IF ( i == n-1 ) THEN
    p = a(k+5)/a(k)
    a(k+5) = a(k+6) - p*a(k+1)
    a(k+6) = a(k+7)
    d(n) = d(n) - p*d(n-2)
  END IF
END DO
d(n) = d(n)/a(3*n-1)
DO i = 3,n
  j = n + 2 - i
  d(j) = (d(j)-a(3*j)*d(j+1))/a(3*j-1)
END DO
d(1) = (d(1)-d(2)*a(2)-d(3)*a(3))/a(1)
a(1) = 0.d0
99001 FORMAT (' RETURN FROM DERSPL  ',i3,' OUT OF ORDER')
END SUBROUTINE derspl
