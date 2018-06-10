SUBROUTINE coredir(it,c,e,l,mj,way,vv,bb,rc,drdic,dovrc,nmatch,  &
        nzero,gc,fc,dp,dq,wp,wq,pow,qow,piw,qiw,cgd,  &
        cgmd,cgo,nrc,z,nucleus)
!   ********************************************************************
!   *                                                                  *
!   *   ROUTINE TO SOLVE RADIAL SPIN-POLARISED DIRAC EQUATIONS         *
!   *   FOR THE CORE WAVE FUNCTIONS                                    *
!   *                                                                  *
!   *   SIMILAR TO LOUCKS' METHOD TO SOLVE THE COUPLED SET OF          *
!   *   DIFFERENTIAL EQUATIONS                                         *
!   *                                    HE JAN. 1989                  *
!   *                                                                  *
!   *   ADOPTED FOR FINITE NUCLEUS       MB MAR. 1995                  *
!   *                                                                  *
!   *                                                                  *
!   ********************************************************************
use mod_types, only: t_inc
IMPLICIT NONE

! PARAMETER definitions
INTEGER MPSMAX,NPEMAX,INVMAX
PARAMETER (MPSMAX=20,NPEMAX=20,INVMAX=3)
DOUBLE PRECISION TOL
PARAMETER (TOL=1.0D-9)
INTEGER ITMAX
PARAMETER (ITMAX=50)

! Dummy arguments
DOUBLE PRECISION C,CGO,MJ
DOUBLE PRECISION E
INTEGER IT,L,NMATCH,NRC,NUCLEUS,NZERO,Z
CHARACTER*3 WAY
DOUBLE PRECISION BB(NRC),CGD(2),CGMD(2),DOVRC(NRC),DP(2,2,NRC), &
       DQ(2,2,NRC), &
       DRDIC(NRC),FC(2,2,NRC),GC(2,2,NRC),PIW(2,2),POW(2,2), &
       QIW(2,2),QOW(2,2),RC(NRC),VV(NRC),WP(2,2,NRC),WQ(2,2,NRC)

! Local variables
DOUBLE PRECISION A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,BB1, &
       BB2,BETA, &
       BOVA,BPP,BQQ,DIFFA,DIFFB,DMUE,EMVPP,EMVQQ,W1,W2,W3,W4,W5, &
       W6,W7
DOUBLE PRECISION BC(0:NPEMAX),CG1,CG2,CG4,CG5,CG8,CSQR,DET,DVC, &
       GAM(2),GPM, &
       H24,KAP(2),PC(2,2,0:MPSMAX),PNEW(2,2),POLD(2,2), &
       QC(2,2,0:MPSMAX),QNEW(2,2),QOLD(2,2),RPWGPM,RR,TZ, &
       VC(0:NPEMAX)
DOUBLE PRECISION DABS,DBLE,DEXP,DSQRT
INTEGER I,IV,J,JCORR,K,KAP1,KAP2,M,MPS,N,NN,NSOL
INTEGER INT,NINT
SAVE A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,BB1,BB2,BC,BETA, &
     BOVA,BPP,BQQ,CG1,CG2,CG4,CG5,CG8,CSQR,DET,DIFFA,DIFFB,DMUE, &
     DVC,EMVPP,EMVQQ,GAM,GPM,H24,I,IV,J,JCORR,K,KAP,KAP1,KAP2,M, &
     MPS,N,NN,NSOL,PC,PNEW,POLD,QC,QNEW,QOLD,RPWGPM,RR,TZ,VC,W1, &
     W2,W3,W4,W5,W6,W7

! MB
!     DOUBLE PRECISION CM(INVMAX,INVMAX),CMI(INVMAX,INVMAX)
! MB


h24 = 1.0D0/24.0D0
dvc = c
csqr = dvc*dvc

!     EXPANSION COEFFICIENTS FOR THE POTENTIAL AND B-FIELD
! MB
IF ( nucleus == 0 ) THEN
  tz = DBLE(nint(-vv(1)*rc(1)))
  vc(0) = vv(1) - (-tz)/rc(1)
ELSE
  tz = 2.0D0*DBLE(z)
  vc(0) = vv(1)
END IF
DO i = 1,2
  DO j = 1,2
    DO k = 1,npemax
      pc(i,j,k) = 0.0D0
      qc(i,j,k) = 0.0D0
    END DO
  END DO
END DO
! MB

bc(0) = bb(1)


!    CALCULATE G-COEFFICIENTS OF B-FIELD

kap1 = -l - 1
kap2 = +l

cg1 = -mj/(kap1+0.5D0)
cg5 = -mj/(-kap1+0.5D0)
cgd(1) = cg1
cgmd(1) = cg5
kap(1) = DBLE(kap1)
! MB
IF ( nucleus == 0 ) THEN
  gam(1) = DSQRT(kap(1)**2-(tz/dvc)**2)
ELSE
  gam(1) = DABS(kap(1))
END IF
! MB
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
ELSE
  cg2 = -DSQRT(1.0D0-(mj/(kap1+0.5D0))**2)
  cg4 = -mj/(kap2+0.5D0)
  cg8 = -mj/(-kap2+0.5D0)
  nsol = 2
  cgd(2) = cg4
  cgo = cg2
  cgmd(2) = cg8
  kap(2) = DBLE(kap2)
!MBA
  IF ( nucleus == 0 ) THEN
    gam(2) = DSQRT(kap(2)**2-(tz/dvc)**2)
  ELSE
    gam(2) = DABS(kap(2))
  END IF
!MBE
END IF


IF ( way == 'INW' ) THEN
  
! #####################################################################
! #####################################################################
! #####################################################################
  
!             INWARD INTEGRATION
  
  
  
  
  dmue = SQRT(-e-e*e/csqr)
  bova = -dmue/(1.0D0+e/csqr)
  
  DO n = (nzero-3),nzero
    
    rr = rc(n)
    
    DO j = 1,nsol
      i = 3 - j
      wp(j,j,n) = DEXP(-dmue*rr)
      dp(j,j,n) = -dmue*drdic(n)*wp(j,j,n)
      wq(j,j,n) = bova*wp(j,j,n)
      dq(j,j,n) = bova*dp(j,j,n)
      
      wp(i,j,n) = 0.0D0
      wq(i,j,n) = 0.0D0
      dp(i,j,n) = 0.0D0
      dq(i,j,n) = 0.0D0
    END DO
  END DO
  
! =============================================================== N ====
  
!     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)
  
  DO nn = 1,(nzero-3-nmatch)
    n = nzero - 3 - nn
    
!    EVALUATE PREDICTOR
    
    DO j = 1,nsol
      DO i = 1,nsol
        pnew(i,j) = wp(i,j,n+1) - h24*(55.0D0*dp(i,j,n+1)-59.0D0*dp(i,j,  &
            n+2)+37.0D0*dp(i,j,n+3)-9.0D0*dp(i,j,n+4))
        qnew(i,j) = wq(i,j,n+1) - h24*(55.0D0*dq(i,j,n+1)-59.0D0*dq(i,j,  &
            n+2)+37.0D0*dq(i,j,n+3)-9.0D0*dq(i,j,n+4))
      END DO
    END DO
    
    emvqq = (e-vv(n)+csqr)*drdic(n)/csqr
    emvpp = -(e-vv(n))*drdic(n)
    bqq = bb(n)*drdic(n)/csqr
    bpp = bb(n)*drdic(n)
    
!    EVALUATE CORRECTOR
    
    DO jcorr = 1,itmax
      
      DO j = 1,nsol
        DO i = 1,nsol
          pold(i,j) = pnew(i,j)
          qold(i,j) = qnew(i,j)
          dp(i,j,n) = -kap(i)*pnew(i,j)*dovrc(n)  &
              + (emvqq+bqq*cgmd(i))*qnew(i,j)
          dq(i,j,n) = kap(i)*qnew(i,j)*dovrc(n)  &
              + (emvpp+bpp*cgd(i))*pnew(i,j) + bpp*cgo*pnew(3-i,j)
          
          pnew(i,j) = wp(i,j,n+1) - h24*(9.0D0*dp(i,j,n)+19.0D0*dp(i,j,  &
              n+1)-5.0D0*dp(i,j,n+2)+dp(i,j,n+3))
          qnew(i,j) = wq(i,j,n+1) - h24*(9.0D0*dq(i,j,n)+19.0D0*dq(i,j,  &
              n+1)-5.0D0*dq(i,j,n+2)+dq(i,j,n+3))
        END DO
      END DO
      
      DO j = 1,nsol
        DO i = 1,nsol
          diffa = pold(i,j) - pnew(i,j)
          IF ( ABS(diffa) > (tol*ABS(pnew(i,j))) ) GO TO 20
          IF ( ABS(diffa) > (tol*ABS(pnew(i,j))) ) GO TO 20
          
          diffb = qold(i,j) - qnew(i,j)
          IF ( ABS(diffb) > (tol*ABS(qnew(i,j))) ) GO TO 20
          IF ( ABS(diffb) > (tol*ABS(qnew(i,j))) ) GO TO 20
        END DO
      END DO
      GO TO 40
      
    20         END DO
    IF(t_inc%i_write>0) WRITE (1337,99001) kap1,n,rc(n),diffa,  &
        diffb,it,l,INT(2*mj),' IN'
    
!                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
    
    
    
    40         CONTINUE
    DO j = 1,nsol
      DO i = 1,nsol
        wp(i,j,n) = pnew(i,j)
        wq(i,j,n) = qnew(i,j)
        dp(i,j,n) = -kap(i)*pnew(i,j)*dovrc(n) + (emvqq+bqq*cgmd(i))*qnew(i,j)
        dq(i,j,n) = kap(i)*qnew(i,j)*dovrc(n)  &
            + (emvpp+bpp*cgd(i))*pnew(i,j) + bpp*cgo*pnew(3-i,j)
      END DO
    END DO
    
  END DO
! =============================================================== N ====
  
  
!     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
  
  DO n = nmatch,nzero
    DO j = 1,nsol
      DO i = 1,nsol
        gc(i,j,n) = wp(i,j,n)/rc(n)
        fc(i,j,n) = wq(i,j,n)/(rc(n)*dvc)
      END DO
    END DO
  END DO
  
  DO j = 1,nsol
    DO i = 1,nsol
      piw(i,j) = wp(i,j,nmatch)
      qiw(i,j) = wq(i,j,nmatch)
    END DO
  END DO
  GO TO 99999
END IF

! #####################################################################
! #####################################################################
! #####################################################################

!             OUTWARD INTEGRATION





!  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS

mps = 20

aa12 = -tz/csqr
aa21 = tz
emvqq = (e-vc(0)+csqr)/csqr
emvpp = -e + vc(0)
bqq = bc(0)/csqr
!MBA
IF ( nucleus == 0 ) THEN
!MBE
  DO j = 1,nsol
    i = 3 - j
    pc(j,j,0) = DSQRT(ABS(kap(j))-gam(j))
    qc(j,j,0) = (kap(j)+gam(j))*(csqr/tz)*pc(j,j,0)
    pc(i,j,0) = 0.0D0
    qc(i,j,0) = 0.0D0
  END DO
  
  DO j = 1,nsol
    
    DO m = 1,mps
      DO i = 1,nsol
        bb1 = (emvqq+bqq*cgmd(i))*qc(i,j,m-1)
        bb2 = (emvpp+bc(0)*cgd(i))*pc(i,j,m-1) + bc(0) *cgo*pc(3-i,j,m-1)
        aa11 = gam(j) + m + kap(i)
        aa22 = gam(j) + m - kap(i)
        det = aa11*aa22 - aa12*aa21
        pc(i,j,m) = (bb1*aa22-aa12*bb2)/det
        qc(i,j,m) = (aa11*bb2-bb1*aa21)/det
      END DO
    END DO
    
  END DO
! MBA
ELSE
! EXPANSION ADAPTED FOR POTENTIALS WITH FINITE NUCLEUS
! EXPANSION OF POTENTIAL actually UP TO zeroth ORDER
  
!       DO IV=1,INVMAX
!        DO N=1,INVMAX
!         CM(N,IV)=RC(N)**(IV-1)
!        ENDDO
!       ENDDO
  
!       CALL RINVGJ(CMI,CM,INVMAX,INVMAX)
  DO iv = 1,invmax
    vc(iv-1) = 0.0D0
!        DO N=1,INVMAX
!         VC(IV-1)=VC(IV-1)+CMI(IV,N)*VV(N)
!        ENDDO
  END DO
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
      IF ( a11 /= 0 ) pc(i,j,2) = (w1*qc(i,j,1)-w2*qc(i,j,0)) /a11
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
        IF ( a21 /= 0 ) pc(i,j,m)  &
            = (w1*qc(i,j,m-1)-w2*qc(i,j,m-2)-w6*qc(i,j,m-3)) /a21
        IF ( a22 /= 0 ) qc(i,j,m) = (-w3*pc(i,j,m-1)+w4*pc(3-i,j,m-1)  &
            +w5*pc(i,j,m-2)+w7*pc(i,j,m-3))/a22
      END DO
    END DO
  END DO
END IF
!MBE

!  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
!  FOR THE FIRST 4 R - MESH - POINTS

DO n = 1,4
  rr = rc(n)
  
  DO j = 1,nsol
    rpwgpm = rr**gam(j)
    
    DO i = 1,nsol
      wp(i,j,n) = pc(i,j,0)*rpwgpm
      wq(i,j,n) = qc(i,j,0)*rpwgpm
!      print*,WP(i,j,N), WQ(I,J,N),0,i,j
      dp(i,j,n) = pc(i,j,0)*rpwgpm*gam(j)*dovrc(n)
      dq(i,j,n) = qc(i,j,0)*rpwgpm*gam(j)*dovrc(n)
    END DO
    
    DO m = 1,mps
      rpwgpm = rpwgpm*rr
      gpm = gam(j) + m
      
      DO i = 1,nsol
        wp(i,j,n) = wp(i,j,n) + pc(i,j,m)*rpwgpm
        wq(i,j,n) = wq(i,j,n) + qc(i,j,m)*rpwgpm
!       print*,WP(i,j,N),WQ(I,J,N) ,m,i,j
        dp(i,j,n) = dp(i,j,n) + pc(i,j,m)*rpwgpm*gpm*dovrc(n)
        dq(i,j,n) = dq(i,j,n) + qc(i,j,m)*rpwgpm*gpm*dovrc(n)
      END DO
      
    END DO
!      if((nsol.eq.2).and.(N.gt.1))stop
  END DO
END DO



! =============================================================== N ====
!     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)

DO n = 5,nmatch
  
!    EVALUATE PREDICTOR
  
  DO j = 1,nsol
    DO i = 1,nsol
      pnew(i,j) = wp(i,j,n-1) + h24*(55.0D0*dp(i,j,n-1)-59.0D0*dp(i,j,n-2)  &
          +37.0D0*dp(i,j,n-3)-9.0D0*dp(i,j,n-4))
      qnew(i,j) = wq(i,j,n-1) + h24*(55.0D0*dq(i,j,n-1)-59.0D0*dq(i,j,n-2)  &
          +37.0D0*dq(i,j,n-3)-9.0D0*dq(i,j,n-4))
    END DO
  END DO
  
  emvqq = (e-vv(n)+csqr)*drdic(n)/csqr
  emvpp = -(e-vv(n))*drdic(n)
  bqq = bb(n)*drdic(n)/csqr
  bpp = bb(n)*drdic(n)
  
!    EVALUATE CORRECTOR
  
  
  DO jcorr = 1,itmax
    
    DO j = 1,nsol
      DO i = 1,nsol
        pold(i,j) = pnew(i,j)
        qold(i,j) = qnew(i,j)
        dp(i,j,n) = -kap(i)*pnew(i,j)*dovrc(n) + (emvqq+bqq*cgmd(i))*qnew(i,j)
        dq(i,j,n) = kap(i)*qnew(i,j)*dovrc(n)  &
            + (emvpp+bpp*cgd(i))*pnew(i,j) + bpp*cgo*pnew(3-i,j)
        
        pnew(i,j) = wp(i,j,n-1) + h24*(9.0D0*dp(i,j,n)+19.0D0*dp(i,j,n-1)  &
            -5.0D0*dp(i,j,n-2)+dp(i,j,n-3))
        qnew(i,j) = wq(i,j,n-1) + h24*(9.0D0*dq(i,j,n)+19.0D0*dq(i,j,n-1)  &
            -5.0D0*dq(i,j,n-2)+dq(i,j,n-3))
      END DO
    END DO
    
    DO j = 1,nsol
      DO i = 1,nsol
        diffa = pold(i,j) - pnew(i,j)
        IF ( ABS(diffa) > (tol*ABS(pnew(i,j))) ) GO TO 50
        IF ( ABS(diffa) > (tol*ABS(pnew(i,j))) ) GO TO 50
        
        diffb = qold(i,j) - qnew(i,j)
        IF ( ABS(diffb) > (tol*ABS(qnew(i,j))) ) GO TO 50
        IF ( ABS(diffb) > (tol*ABS(qnew(i,j))) ) GO TO 50
      END DO
    END DO
    GO TO 100
    
  50      END DO
  IF(t_inc%i_write>0) WRITE (1337,99001) kap1,n,rc(n),diffa,  &
      diffb,it,l,INT(2*mj),'OUT'
  
!                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
  
  
  100     CONTINUE
  DO j = 1,nsol
    DO i = 1,nsol
      wp(i,j,n) = pnew(i,j)
      wq(i,j,n) = qnew(i,j)
      dp(i,j,n) = -kap(i)*pnew(i,j)*dovrc(n) + (emvqq+bqq*cgmd(i))*qnew(i,j)
      dq(i,j,n) = kap(i)*qnew(i,j)*dovrc(n) + (emvpp+bpp*cgd(i))*pnew(i,j)  &
          + bpp*cgo*pnew(3-i,j)
    END DO
  END DO
  
END DO
! =============================================================== N ====

!     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS


DO n = 1,nmatch
  DO j = 1,nsol
    DO i = 1,nsol
      gc(i,j,n) = wp(i,j,n)/rc(n)
      
      
      fc(i,j,n) = wq(i,j,n)/(rc(n)*dvc)
    END DO
  END DO
END DO

DO j = 1,nsol
  DO i = 1,nsol
    pow(i,j) = wp(i,j,nmatch)
    qow(i,j) = wq(i,j,nmatch)
  END DO
END DO

RETURN

99001 FORMAT (' P/C NOT CONV. IN <DIRAC> ',2I4,2X,f10.7,2X,2E12.4,3I2,  &
    '/2 ',a3)
99999 CONTINUE
END SUBROUTINE coredir
