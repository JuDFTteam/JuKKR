SUBROUTINE intcor(f1,f2,rho,g,f,v,value,slope,l,nn,e,sum,nre,  &
        vlnc,a,b,z,rn,nr,tol,irm,ipr,nitmax,nsra)
      use mod_types, only: t_inc
      IMPLICIT NONE
!.. Scalar Arguments ..
      DOUBLE PRECISION A,B,E,F1,F2,RN,SLOPE,SUM,TOL,VALUE,Z
      INTEGER IPR,IRM,L,NITMAX,NN,NR,NRE,NSRA
      LOGICAL VLNC
!..
!.. Array Arguments ..
      DOUBLE PRECISION F(*),G(*),RHO(*),V(*)
!..
!.. Local Scalars ..
DOUBLE COMPLEX ARG,CAPPAI,DOFE
DOUBLE PRECISION CVLIGHT,DE,DG1,DG2,DPSI1,DPSI2,DRDIKC,E1,E2,EA, &
                 GKC2,PI,PKC1,PKC2,PSI1,PSI2,Q,QKC1,QKC2,RATIO, &
                 RATIO1,RE,RKC,RPB,SLOP,TSRME,VALU,VME,XXX,ZZ
INTEGER IR,K,K2,KC,NITER,NNE,NREM1,NREM2
!..
!.. Local Arrays ..
      DOUBLE COMPLEX HL(6)
!..
!.. External Subroutines ..
      EXTERNAL HANKEL,INTIN,INTOUT
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DCMPLX,DSQRT,EXP,LOG,MAX0,MIN0,REAL,SQRT
!..
pi = 4.d0*ATAN(1.d0)
zz = z + z
cvlight = 274.0720442D0
IF (nsra == 1) cvlight = 1.0D0
ea = EXP(a)
niter = 0
e1 = f1
e2 = f2
IF (ipr == 2.AND.(t_inc%i_write>0)) WRITE (1337,FMT=9000)  &
    l,nn,nr,f1,e,f2,tol,value,slope
10 CONTINUE
niter = niter + 1
DO  ir = 1,irm
  g(ir) = 0.0D0
  f(ir) = 0.0D0
END DO
IF (niter > nitmax) THEN
  GO TO 80
  
ELSE
  IF (e <= e1 .OR. e >= e2) e = .5D0* (e1+e2)
  nre = nr
  IF (e <= -1.d-8) THEN
    tsrme = 2.d0*SQRT(-e)
    re = (LOG(-tsrme*e/1.d-8)/tsrme-zz/e)*2.d0
    nre = LOG(re/b+1.d0)/a + 1.d0
    nre = (nre/2)*2 + 1
    nre = MIN0(nre,nr)
    nre = MAX0(nre,35)
  END IF
  xxx = 1.d0
  valu = 1.d-1
  slop = -1.d-1
  IF (nre < nr .AND. niter == 1 .AND. ipr /= 0 .AND. (t_inc%i_write>0))  &
      WRITE (1337,FMT=9010)
  IF (nre >= nr) THEN
    valu = value
    slop = slope
    IF (.NOT.vlnc) THEN
!--->   single site  boundary condition
      vme = -e
      IF (nsra == 1) THEN
        cappai = DCMPLX(0.d0,DSQRT(vme))
      ELSE
        cappai = DCMPLX(0.d0,DSQRT((1.d0- vme/cvlight/cvlight)*vme))
      END IF
      arg = cappai*rn
      CALL hankel(hl,l+2,arg)
      dofe = REAL(l+1)/rn - cappai*hl(l+2)/hl(l+1)
      valu = 1.d-10
      slop = valu*dofe
    END IF
    
  END IF
  k2 = 30
  IF (nn == 0) k2 = nre/3
  nne = 0
  
  CALL intin(g,f,v,e,l,nne,valu,slop,nre,k2,kc,dg2,a,b,z,nsra)
  
  rkc = b*EXP(a*kc-a) - b
  drdikc = a* (rkc+b)
  gkc2 = g(kc)
  psi2 = g(kc)
  dpsi2 = dg2/drdikc
  qkc2 = psi2*psi2 + dpsi2*dpsi2*rkc*rkc
  pkc2 = .5D0 - ATAN(rkc*dpsi2/psi2)/pi
  
  CALL intout(g,f,v,e,l,nne,kc,dg1,a,b,z,nsra)
  
  psi1 = g(kc)
  dpsi1 = dg1/drdikc
  qkc1 = psi1*psi1 + dpsi1*dpsi1*rkc*rkc
  pkc1 = .5D0 - ATAN(rkc*dpsi1/psi1)/pi
  IF (nne == 9) nne = 0
  IF (nne == nn) THEN
    ratio1 = gkc2/g(kc)
    ratio = SQRT(qkc2/qkc1)
    IF (ratio1 < 0.d0) ratio = -ratio
    DO  k = 1,kc
      g(k) = g(k)*ratio
      f(k) = f(k)*ratio
    END DO
    sum = 0.d0
    IF (nsra == 1) THEN
      DO  k = 1,nre
        f(k) = 0.0D0
      END DO
    END IF
    rpb = b/ea
    q = ea*ea
    nrem1 = nre - 1
    DO  k = 2,nrem1,2
      rpb = rpb*q
      sum = sum + rpb* (g(k)*g(k)+f(k)*f(k))
    END DO
    rpb = b
    sum = sum + sum
    nrem2 = nre - 2
    DO  k = 3,nrem2,2
      rpb = rpb*q
      sum = sum + rpb* (g(k)*g(k)+f(k)*f(k))
    END DO
    sum = sum + sum + rpb*q* (g(nre)*g(nre)+f(nre)*f(nre))
    sum = a*sum/3.d0
    de = pi*qkc2* (pkc2-pkc1)/sum/rkc
    IF (niter >= nitmax-10 .OR. ipr == 2  &
        .AND. (t_inc%i_write>0)) WRITE (1337,  &
        FMT=9020) niter,nne,nre,kc,e1,e,e2,de
    IF (de > 0.d0) e1 = e
    IF (de < 0.d0) e2 = e
    e = e + de
    IF (ABS(de) > tol .AND. niter < nitmax) GO TO 10
    
  ELSE
    IF (niter >= nitmax-10 .OR. ipr == 2  &
        .AND. (t_inc%i_write>0)) WRITE (1337, FMT=9020) niter,nne,nre,kc,e1,e,e2
    IF (nne > nn) e2 = e
    IF (nne < nn) e1 = e
    e = .5D0* (e1+e2)
    GO TO 10
    
  END IF
  
END IF

e = e - de
DO  k = 1,nre
  rho(k) = g(k)*g(k) + f(k)*f(k)
END DO
IF (xxx <= 0.d0 .AND. (t_inc%i_write>0)) WRITE (1337,FMT=9030)
IF (niter >= nitmax-10 .OR. ipr >= 1 .OR. xxx <= 0.d0  &
    .AND. (t_inc%i_write>0))  &
    WRITE (1337,FMT=9040) l,nn,niter,kc,nre,valu,slop,e,de,sum
RETURN

80 WRITE (6,FMT=9050) nitmax
STOP 'INTCOR'


9000 FORMAT (' l=',i3,'  nn=',i2,'  nr=',i4,'  f1/e/f2=',3F10.3,/,  &
    ' tol=',1P,d12.3,'  value/slope=',2D12.3)
9010 FORMAT (13X,'  no boundary condition had to be used')
9020 FORMAT (2I3,2I4,1P,3D16.8,1X,2D9.2)
9030 FORMAT (/,' **** int: 0-pressure bcs not real')
9040 FORMAT (' state',i2,',',i1,':',i4,'x,',i5,'/',i3,',  bc=',1P,  &
    2D12.3,/,14X,'e=',d14.6,'   de=',d11.2,'   sum=',d12.4)
9050 FORMAT (' *** int: stop after',i4,' iterations')
END SUBROUTINE intcor
