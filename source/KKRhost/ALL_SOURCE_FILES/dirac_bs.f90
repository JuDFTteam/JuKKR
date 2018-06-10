SUBROUTINE dirbs(getirrsol,c,e,l,mj,kap1,kap2,pis,cg1,cg2,cg4,  &
        cg5,cg8,v,b,z,nucleus,r,drdi,dovr,nmesh,pr,qr,pi,  &
        qi,dp,dq)
!   ********************************************************************
!   *                                                                  *
!   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
!   *                                                                  *
!   *   the outward integration is started by a power expansion        *
!   *   and the inward integration is started analytically             *
!   *   the integration itself is done by the BURLISCH-STOER method    *
!   *   see: numerical recipes chapter 15.4                            *
!   *                                                                  *
!   *   returns the wave functions up to the mesh point NMESH          *
!   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
!   *   and    R/I standing for regular/irregular solution             *
!   *                                                                  *
!   *   bug fixed 93/11/24                                             *
!   *  31/10/94  HE  arg. list changed - return P,Q instead of g,f     *
!   *  06/12/94  HE  CM real                                           *
!   *  29/04/95  MB  Adopted for finite nucleus                        *
!   ********************************************************************
IMPLICIT NONE
INCLUDE 'sprkkr_rmesh.dim'

! PARAMETER definitions
INTEGER MPSMAX,NPEMAX,NABM
PARAMETER (MPSMAX=40,NPEMAX=4,NABM=5)
DOUBLE COMPLEX C0
PARAMETER (C0=(0.0D0,0.0D0))
DOUBLE PRECISION EPSBS
PARAMETER (EPSBS=2.0D-7)

! COMMON variables
DOUBLE PRECISION CGD(2),CGMD(2),CGO,CSQR,KAP(2)
DOUBLE COMPLEX EBS
INTEGER NRADBS,NSOLBS
COMMON /COMMBS/ EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,NRADBS

! Dummy arguments
DOUBLE PRECISION C,CG1,CG2,CG4,CG5,CG8,MJ
DOUBLE COMPLEX E
LOGICAL GETIRRSOL
INTEGER KAP1,KAP2,L,NMESH,NUCLEUS,Z
DOUBLE COMPLEX PIS
DOUBLE PRECISION B(NRMAX),DOVR(NRMAX),DRDI(NRMAX),R(NRMAX), &
                 V(NRMAX)
DOUBLE COMPLEX DP(2,2,NRMAX),DQ(2,2,NRMAX),PI(2,2,NRMAX), &
               PR(2,2,NRMAX)
DOUBLE COMPLEX QI(2,2,NRMAX),QR(2,2,NRMAX)

! Local variables
DOUBLE COMPLEX A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA, &
               BB1,BB2,BETA
DOUBLE COMPLEX BQQ,CFAC,EMVPP,EMVQQ,W1,W2,W3,W4,W5,W6,W7
DOUBLE PRECISION BC(0:NPEMAX),CM(NPEMAX,NPEMAX), &
                 CMI(NPEMAX,NPEMAX),DIX
DOUBLE PRECISION GAM(2),GPM,HBS,RPWGPM,RR,SK(2),SK1,SK2,TZ, &
                 VC(0:NPEMAX),X
DOUBLE COMPLEX CJLZ
 
DOUBLE COMPLEX DETD,DY(NCFMAX),EFAC,FY(NCFMAX), &
               PC(2,2,-NPEMAX:MPSMAX)
DOUBLE COMPLEX QC(2,2,-NPEMAX:MPSMAX),ZZ

INTEGER I,IP,ISK1,ISK2,IV,J,K,LB(2),LB1,LB2,M,MPS,N,NFY,NPE,NSOL
INTEGER ISIGN
SAVE a11,a12,a21,a22,aa11,aa12,aa21,aa22,alpha,bb1,bb2,bc,beta,  &
    bqq,cfac,cm,cmi,detd,dix,dy,efac,emvpp,emvqq,fy,gam,gpm,hbs,  &
    i,ip,isk1,isk2,iv,j,k,lb,lb1,lb2,m,mps,n,nfy,npe,nsol,pc,qc,  &
    rpwgpm,rr,sk,sk1,sk2,tz,vc,w1,w2,w3,w4,w5,w6,w7,x,zz

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
    IF ( nucleus == 0 ) THEN
      vc(iv-1) = vc(iv-1) + cmi(iv,n)*(v(n)+tz/r(n))
    ELSE
      vc(iv-1) = vc(iv-1) + cmi(iv,n)*v(n)
    END IF
  END DO
END DO

DO iv = 1,npe
  bc(iv-1) = 0.0D0
  DO n = 1,npe
    bc(iv-1) = bc(iv-1) + cmi(iv,n)*b(n)
  END DO
END DO



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
! MB
IF ( nucleus == 0 ) THEN
  gam(1) = DSQRT(kap(1)**2-(tz/c)**2)
ELSE
  gam(1) = DABS(kap(1))
END IF
! MB
lb(1) = lb1
sk(1) = sk1
IF ( DABS(mj) > l ) THEN
  cg2 = 0.0D0
  cg4 = 0.0D0
  cg8 = 0.0D0
  nsol = 1
  cgd(2) = 0.0D0
  cgo = 0.0D0
  cgmd(2) = 0.0D0
  gam(2) = 0.0D0
  kap(2) = 0.0D0
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
! MB
  IF ( nucleus == 0 ) THEN
    gam(2) = DSQRT(kap(2)**2-(tz/c)**2)
  ELSE
    gam(2) = DABS(kap(2))
  END IF
! MB
  lb(2) = lb2
  sk(2) = sk2
END IF

nsolbs = nsol
ebs = e

DO i = 1,2
  DO j = 1,2
    DO ip = -npemax,mpsmax
      pc(i,j,ip) = c0
      qc(i,j,ip) = c0
    END DO
  END DO
END DO

! ======================================================================
IF ( tz >= 2 ) THEN
  
  DO j = 1,nsol
    i = 3 - j
    pc(j,j,0) = DSQRT(ABS(kap(j))-gam(j))
    qc(j,j,0) = (kap(j)+gam(j))*(csqr/tz)*pc(j,j,0)
    pc(i,j,0) = c0
    qc(i,j,0) = c0
  END DO
  
!  determine higher expansion coefficients for the wave functions
  
  mps = 40
  
  aa12 = -tz/csqr
  aa21 = tz
  emvqq = (e-vc(0)+csqr)/csqr
  emvpp = -e + vc(0)
  bqq = bc(0)/csqr
!MBA
  IF ( nucleus == 0 ) THEN
!MBE
    
    DO j = 1,nsol
      
      DO m = 1,mps
        DO i = 1,nsol
          k = 3 - i
          bb1 = (emvqq+bqq*cgmd(i))*qc(i,j,m-1)
          bb2 = (emvpp+bc(0)*cgd(i))*pc(i,j,m-1) + bc(0) *cgo*pc(k,j,m-1)
          DO ip = 1,npe - 1
            bb1 = bb1 + (-vc(ip)+bc(ip)*cgmd(i)) *qc(i,j,m-1-ip)/csqr
            bb2 = bb2 + (+vc(ip)+bc(ip)*cgd(i))  &
                *pc(i,j,m-1-ip) + (+bc(ip)*cgo) *pc(k,j,m-1-ip)
          END DO
          
          aa11 = gam(j) + m + kap(i)
          aa22 = gam(j) + m - kap(i)
          detd = aa11*aa22 - aa12*aa21
          pc(i,j,m) = (bb1*aa22-aa12*bb2)/detd
          qc(i,j,m) = (aa11*bb2-bb1*aa21)/detd
        END DO
      END DO
      
    END DO
! MBA
  ELSE
! EXPANSION ADAPTED FOR POTENTIALS WITH FINITE NUCLEUS
! EXPANSION OF POTENTIAL UP TO SECOND ORDER: V_O+V_1*R+V_2*R*R
    DO j = 1,nsol
      i = 3 - j
      IF ( kap(j) > 0 ) THEN
! ARBITRARY STARTING VALUES
        alpha = 0.0D0
        beta = 0.174D0
      ELSE
        beta = 0.0D0
        alpha = 0.174D0
      END IF
      pc(j,j,0) = alpha
      qc(j,j,0) = beta
      pc(i,j,0) = 0.0D0
      qc(i,j,0) = 0.0D0
    END DO
    w4 = bc(0)*cgo
    w2 = vc(1)/csqr
    w5 = vc(1)
    w6 = vc(2)/csqr
    w7 = vc(2)
    DO j = 1,nsol
      DO i = 1,nsol
        w1 = emvqq + bqq*cgmd(i)
        w3 = -emvpp + bc(0)*cgd(i)
        a11 = gam(j) + kap(i) + 1D0
        a12 = gam(j) - kap(i) + 1D0
        IF ( a11 /= 0 ) pc(i,j,1) = w1/a11*qc(i,j,0)
        IF ( a12 /= 0 ) qc(i,j,1) = (-w3*pc(i,j,0)+w4*pc(3-i,j,0))/a12
      END DO
    END DO
    DO j = 1,nsol
      DO i = 1,nsol
        w1 = emvqq + bqq*cgmd(i)
        w3 = -emvpp + bc(0)*cgd(i)
        a11 = gam(j) + kap(i) + 2D0
        a12 = gam(j) - kap(i) + 2D0
        IF ( a11 /= 0 ) pc(i,j,2) = (w1*qc(i,j,1)-w2*qc(i,j,0))/a11
        IF ( a12 /= 0 ) qc(i,j,2)  &
            = (-w3*pc(i,j,1)+w4*pc(3-i,j,1)+w5*pc(i,j,0))/a12
      END DO
    END DO
    DO j = 1,nsol
      DO m = 3,mps
        DO i = 1,nsol
          w1 = emvqq + bqq*cgmd(i)
          w3 = -emvpp + bc(0)*cgd(i)
          a21 = gam(j) + kap(i) + DBLE(m)
          a22 = gam(j) - kap(i) + DBLE(m)
          IF ( a21 /= 0 ) pc(i,j,m) = (w1*qc(i,j,m-1)-w2*qc(i,j,m-2)  &
              -w6*qc(i,j,m-3))/a21
          IF ( a22 /= 0 ) qc(i,j,m) = (-w3*pc(i,j,m-1)+w4*pc(3-i,j,m-1)  &
              +w5*pc(i,j,m-2)+w7*pc(i,j,m-3))/a22
        END DO
      END DO
    END DO
  END IF
!MBE
  
  
  
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
!                                                     == EMPTY SPHERE ==
ELSE
!                     assume constant pot: V=V(1)   ignore coupling: B=0
  
  DO n = 1,nabm
    zz = CDSQRT(e-v(1))*r(n)
    efac = (zz/r(n))*c/(e+csqr)
    
    DO j = 1,nsol
      i = 3 - j
      pr(j,j,n) = cjlz(l,zz)*r(n)
      qr(j,j,n) = efac*sk(j)*cjlz(lb(j),zz)*r(n)*c
      dp(j,j,n) = (DBLE(l+1)*cjlz(l,zz)-zz*cjlz(l+1,zz)) *drdi(n)
      m = lb(j)
      dq(j,j,n) = efac*sk(j) *(DBLE(m+1)*cjlz(m,zz)-zz*cjlz(m+1,zz))  &
          *drdi(n)*c
      
      pr(i,j,n) = c0
      qr(i,j,n) = c0
      dp(i,j,n) = c0
      dq(i,j,n) = c0
    END DO
  END DO
  
END IF

! =============================================================== n ====
!     DO 400 J=1,NSOL

nfy = 0
DO j = 1,nsol
  DO i = 1,nsol
    fy(nfy+1) = pr(i,j,1)
    fy(nfy+2) = qr(i,j,1)
    dy(nfy+1) = dp(i,j,1)
    dy(nfy+2) = dq(i,j,1)
    nfy = nfy + 2
  END DO
END DO
x = 1.0D0
dix = 1.0D0

DO n = 2,nmesh
  
  nradbs = n
  
  CALL dirbsstp(fy,dy,nfy,x,dix,epsbs,fy,b,v,r,drdi,nmesh)
  
  
  nfy = 0
  DO j = 1,nsol
    DO i = 1,nsol
      pr(i,j,n) = fy(nfy+1)
      qr(i,j,n) = fy(nfy+2)
      nfy = nfy + 2
    END DO
  END DO
  
END DO


IF ( .NOT.getirrsol ) RETURN


! =============================================================== n ====

! ######################################################################
!                     irregular solution
! ######################################################################

!         calculate the initial values of the wavefunction
!                     at the sphere boundary

n = nmesh

zz = pis*r(n)

DO j = 1,nsol
  pi(j,j,n) = cjlz(l,zz)*r(n)
  qi(j,j,n) = cfac*sk(j)*cjlz(lb(j),zz)*r(n)*c
  dp(j,j,n) = (DBLE(l+1)*cjlz(l,zz)-zz*cjlz(l+1,zz)) *drdi(n)
  
  m = lb(j)
  dq(j,j,n) = cfac*sk(j) *(DBLE(m+1)*cjlz(m,zz)-zz*cjlz(m+1,zz))  &
      *drdi(n)*c
  
  i = 3 - j
  pi(i,j,n) = c0
  qi(i,j,n) = c0
  dp(i,j,n) = c0
  dq(i,j,n) = c0
END DO

! =============================================================== n ====
hbs = -1.0D0

nfy = 0
DO j = 1,nsol
  DO i = 1,nsol
    fy(nfy+1) = pi(i,j,nmesh)
    fy(nfy+2) = qi(i,j,nmesh)
    dy(nfy+1) = dp(i,j,nmesh)
    dy(nfy+2) = dq(i,j,nmesh)
    nfy = nfy + 2
  END DO
END DO
x = DBLE(nmesh)

nradbs = nmesh
CALL dirbsrad(x,fy,dy,drdi,b,v,r,nmesh)

DO n = nmesh - 1,1, - 1
  nradbs = n
  
  CALL dirbsstp(fy,dy,nfy,x,hbs,epsbs,fy,b,v,r,drdi,nmesh)
  
  nfy = 0
  DO j = 1,nsol
    DO i = 1,nsol
      pi(i,j,n) = fy(nfy+1)
      qi(i,j,n) = fy(nfy+2)
!              DP(I,J,N) = DY(NFY+1)
!              DQ(I,J,N) = DY(NFY+2)
      nfy = nfy + 2
    END DO
  END DO
  
END DO

! =============================================================== n ====
END SUBROUTINE dirbs
