SUBROUTINE decitmat(eryd,zat,ipan,rr,dror,visp,ircut,rirc,  &
        krel,nsra,ins,tmatll,loflm,  &
        idoldau,lopt,wldauav,  &
        solver,soctl,ctl,zrel,vtrel,btrel,drdi,r2drdi,  &
        ipand,irmd,lmaxd,lmaxdp1,lm2d,lmmaxd)
! **********************************************************************
! *                                                                    *
! * A modified form of the CALCTMAT routine to deal with the host      *
! * t-matrices in case of decimation                                   *
! *                                                                    *
! * Non-spherical potential not implemented yet, neither LDA+U         *
! *                                                                    *
! **********************************************************************
IMPLICIT NONE

!Parameters ..
DOUBLE PRECISION CVLIGHT
PARAMETER ( CVLIGHT=274.0720442D0 )
DOUBLE COMPLEX CI
PARAMETER ( CI = (0D0,1D0) )

!Scalar arguments ..
INTEGER IDOLDAU,IPAN,KREL,LOPT,NSRA,INS,ZREL
INTEGER IPAND,IRMD,LM2D,LMAXD,LMAXDP1,LMMAXD
DOUBLE PRECISION ZAT,RIRC,WLDAUAV
DOUBLE COMPLEX ERYD
CHARACTER*10 SOLVER

!Array arguments ..
INTEGER IRCUT(0:IPAND),LOFLM(LM2D)
DOUBLE PRECISION RR(IRMD),DROR(IRMD),VISP(IRMD)
DOUBLE COMPLEX TMATLL(LMMAXD,LMMAXD)
DOUBLE PRECISION SOCTL(KREL*LMAXD+1)
DOUBLE PRECISION CTL(KREL*LMAXD+1)
DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
DOUBLE PRECISION DRDI(IRMD),R2DRDI(IRMD*KREL+(1-KREL))

!Local scalars ..
INTEGER LL,LM1
DOUBLE PRECISION RIRC1
DOUBLE COMPLEX EK,CARG,QF,HLW,BLW

!Local arrays ..
DOUBLE PRECISION CUTOFF(IRMD)
DOUBLE PRECISION RS(:,:),S(:)
DOUBLE COMPLEX BESSJW(:),BESSYW(:),HANKWS(:),DLOGDP(:)
DOUBLE COMPLEX TMAT(:),MASS(:),HAMF(:,:),FZ(:,:),PZ(:,:)
ALLOCATABLE RS,S
ALLOCATABLE BESSJW,BESSYW,HANKWS,DLOGDP
ALLOCATABLE TMAT,MASS,HAMF,FZ,PZ

!External subroutines ..
EXTERNAL BESHAN,CINIT,REGSOL,WFMESH


CALL cinit(lmmaxd*lmmaxd,tmatll)
! ================================================================= KREL
IF ( krel == 0 ) THEN
  allocate (bessjw(0:lmaxdp1),bessyw(0:lmaxdp1),stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Allocate BESSJW/BESSYW'
  allocate (hankws(0:lmaxdp1),dlogdp(0:lmaxd),stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Allocate HANKWS/DLOGFP'
  allocate (tmat(0:lmaxd),mass(irmd),stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Allocate TMAT/MASS'
  allocate (hamf(irmd,0:lmaxd),fz(irmd,0:lmaxd),stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Allocate HAMF/FZ'
  allocate (pz(irmd,0:lmaxd),stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Allocate PZ'
  allocate (rs(irmd,0:lmaxd),s(0:lmaxd),stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Allocate RS/S'
  rirc1 = 1D0/rirc
  CALL wfmesh(eryd,ek,cvlight,nsra,zat,rr,s,rs,ircut(ipan), irmd,lmaxd)
  
  carg = rirc * ek
  CALL beshan(hankws,bessjw,bessyw,carg,lmaxdp1)
  DO ll = 0,lmaxdp1
    hankws(ll) = bessyw(ll) - ci*bessjw(ll)
  END DO
  
  CALL regsol(cvlight,eryd,nsra,dlogdp,fz,hamf,mass,pz,  &
      dror,rr,s,visp,zat,ipan,ircut,idoldau,lopt,wldauav, cutoff,irmd,ipand,lmaxd)
  
! ----------------------------------------------------------------------
! --> determine KREL=0 t - matrix
  
  DO ll = 0,lmaxd
    qf = DBLE(ll)*rirc1
    hlw = hankws(ll) * dlogdp(ll)
    blw = bessjw(ll) * dlogdp(ll)
    
    hlw = qf*hankws(ll) - ek*hankws(ll+1) - hlw
    blw = blw - qf*bessjw(ll) + ek*bessjw(ll+1)
    hlw = hlw * ek
    tmat(ll) = blw/hlw
  END DO
  
! --> spherical/non-spherical
  
  IF ( ins == 0 ) THEN
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    DO lm1 = 1,lmmaxd
      tmatll(lm1,lm1) = tmat(loflm(lm1))
    END DO
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    STOP ' not implemented'
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  END IF
  deallocate (bessjw,bessyw,hankws,dlogdp,stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Deallocate'
  deallocate (tmat,mass,hamf,fz,pz,stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Deallocate'
  deallocate (rs,s,stat=lm1)
  IF ( lm1 /= 0 ) STOP '    Deallocate'
! ----------------------------------------------------------------------
ELSE                      ! KREL
  CALL drvreltmat(eryd,tmatll,vtrel,btrel,rr,  &
      drdi,r2drdi,zrel,ircut(ipan),solver,soctl, ctl,lmmaxd,lmaxd,irmd)
END IF
! ================================================================= KREL
END SUBROUTINE decitmat
