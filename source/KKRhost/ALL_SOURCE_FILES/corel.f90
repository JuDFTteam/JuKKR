SUBROUTINE corel(nsra,ipr,ip,rhoc,v,ecore,lcore,ncore,drdi,z,qc,  &
        a,b,is,nspin,nr,rmax,irmd)
!-----------------------------------------------------------------------
!     subroutine for core states
!-----------------------------------------------------------------------
!     lmxc = lmaxcore = (0,1,2,...), .e.g, argon core : lmxc = 1
!                                        krypton core : lmxc = 2
!     kfg = configuration of core, e.g., argon core: 3300=3s,3p,0d
!                                      krypton core: 4430=4s,4p,3d
!                                      xenon core: 5540=5s,5p,4d
!-----------------------------------------------------------------------
      use mod_types, only: t_inc
      IMPLICIT NONE
!.. Parameters ..
      INTEGER NITMAX,IRNUMX
      PARAMETER (NITMAX=40,IRNUMX=10)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION A,B,QC,RMAX,Z
      INTEGER IP,IPR,IRMD,IS,NCORE,NR,NSPIN,NSRA
!..
!.. Array Arguments ..
      DOUBLE PRECISION DRDI(*),ECORE(*),RHOC(*),V(*)
      INTEGER LCORE(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION E,E1,E2,EDIFF,EI,SLOPE,SUM,TOL,VALUE,WGT
      INTEGER IC,IN,INUC,IR,L,LMP1,LMXC,LP1,NC,NMAX,NN,NRE
      LOGICAL VLNC
!..
!.. Local Arrays ..
      DOUBLE PRECISION F(IRMD),G(IRMD),RHO(IRMD)
      INTEGER KFG(4)
      CHARACTER*4 SPN(2),TEXT(5)
!..
!.. External Subroutines ..
      EXTERNAL INTCOR,SIMP3
!..
!.. Intrinsic Functions ..
      INTRINSIC DBLE,REAL
!..
!.. Save statement ..
      SAVE SPN,TEXT
!..
!.. Data statements ..
      DATA SPN,TEXT/'down','up  ','s   ','p   ','d   ','f   ','g   '/
!..
vlnc = .false.
value = 1.d-8
slope = -1.d-8
e2 = 50.0D0

DO  ic = 1,4
  kfg(ic) = 0
END DO
DO  ic = 1,ncore
  IF (lcore(ic) == 0) kfg(1) = kfg(1) + 1
  IF (lcore(ic) == 1) kfg(2) = kfg(2) + 1
  IF (lcore(ic) == 2) kfg(3) = kfg(3) + 1
  IF (lcore(ic) == 3) kfg(4) = kfg(4) + 1
END DO
IF (kfg(2) /= 0) kfg(2) = kfg(2) + 1
IF (kfg(3) /= 0) kfg(3) = kfg(3) + 2
IF (kfg(4) /= 0) kfg(4) = kfg(4) + 3
lmxc = 0
IF (kfg(2) /= 0) lmxc = 1
IF (kfg(3) /= 0) lmxc = 2
IF (kfg(4) /= 0) lmxc = 3

tol = 1.0D-12* (z*z+1.d0)
lmp1 = lmxc + 1
nc = 0
inuc = -irnumx

DO  ir = 1,irmd
  rhoc(ir) = zero
  rho(ir) = zero
END DO

DO  lp1 = 1,lmp1
  l = lp1 - 1
  e1 = (-5.d0- ((z+1.d0)/DBLE(lp1))**2)*1.5D0 - 50.d0
  nmax = kfg(lp1)
  IF (nmax /= 0) THEN
    DO  in = lp1,nmax
      nn = in - lp1
      nc = nc + 1
      inuc = inuc + irnumx
      e = ecore(nc)
      ei = ecore(nc)
      IF((t_inc%i_write>0).AND.(ipr /= 0))  &
          WRITE (1337,FMT=9000) in,text(lp1),nn,spn(is),ip,e
      CALL intcor(e1,e2,rho,g,f,v,value,slope,l,nn,e,sum,nre,  &
          vlnc,a,b,z,rmax,nr,tol,irmd,ipr,nitmax,nsra)
      ediff = e - ei
      ecore(nc) = e
      wgt = REAL(l+l+1)/sum*2.d0/REAL(nspin)
      IF((t_inc%i_write>0).AND.(ipr /= 0)) WRITE (1337,FMT=9010) ei,ediff,e
      40       CONTINUE
      
!---> sum up contributions to total core charge
      
      DO  ir = 2,nre
        rhoc(ir) = rhoc(ir) + rho(ir)*wgt
        rho(ir) = zero
      END DO
    END DO
  END IF
  
END DO
IF (nc*irnumx > 150 .OR. irnumx > 10) STOP 'corel'

!---> integrate core density to get core charge

CALL simp3(rhoc,qc,1,nr,drdi)

9000 FORMAT (1X,90 ('*'),/,'  n = ',i1,'  l = ',a4,'   nnode = ',i1,  &
    '  spin=',a4,i5,'th cell','    einput = ',1P,d16.8)
9010 FORMAT (1X,'  einput =',1P,d16.8,'   eout - ein =',1P,d16.8,  &
    '   eoutput = ',1P,d16.8)
END SUBROUTINE corel
