!     This subroutine computes a matrix that is the basis for constructing
!     the KKR representation of the torque operator. It is adapted from the
!     CALC_RHO_LL_SS subroutine, but the spin dependent part, i.e., the exhange
!     field, replaces the shape function in the integration.
!
!                                     Guillaume Geranton, September 2014
SUBROUTINE calc_torq_ll_ss(lmmax,rll,ircut,ipan,icell,  &
    cleb,icleb,iend,ifunm,lmsp,irws,drdi,dens, visp,nspin,iatom,vins,irmin)

IMPLICIT NONE

INCLUDE 'inc.p'
INTEGER :: lmmaxd
PARAMETER        (lmmaxd= (lmaxd+1)**2)
INTEGER :: lmpotd
PARAMETER        (lmpotd= (lpotd+1)**2)
INTEGER :: irmind
PARAMETER        (irmind=irmd-irnsd)
!     ..
!     .. Scalar Arguments ..
INTEGER :: iend,lmmax,irws,nspin,iatom,irmin!,LM1,LM2
!     ..
!     .. Array Arguments ..
DOUBLE COMPLEX   rll(irmd,lmmaxd,lmmaxd), &  ! non-sph. eigen states of single pot  &
    dens
DOUBLE PRECISION :: cleb(*),  &
    drdi(irmd),  &                          ! derivative dr/di  &
    visp(irmd,*),& !              spherical part of the potential  &
    vins(irmind:irmd,lmpotd,*) ! non-sph. part of the potential
INTEGER :: icleb(ncleb,4),ifunm(natypd,lmpotd),  &
    lmsp(natypd,*),ircut(0:ipand),ipan, icell,ifun


! local variables
DOUBLE PRECISION ::   c0ll
DOUBLE COMPLEX               ::   clt
DOUBLE COMPLEX, allocatable  ::   rsp(:),rges(:)
INTEGER ::   lm1p,lm2p,lm3p,ir,j,i
INTEGER ::   ircutm(0:ipand)

EXTERNAL test
LOGICAL :: test

!     ..
!  ---> first calculate only the spherically symmetric contribution
!       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
!       multiplied with the shape functions...

!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)

allocate(rges(irmd))
allocate(rsp(irmd))

c0ll = 1.0D0/DSQRT(16.0D0*DATAN(1.0D0))
rsp=0D0
rges=0D0

!     Compute spherical contribution to the torque (LM=1)
!     Sph. potential has to be multiplied by sqrt(4 PI) !
DO lm1p = 1,lmmax
  DO ir=1,irmd
    rsp(ir)=rsp(ir)+rll(ir,lm1p,lm1p)*c0ll*(-1)*  &
        (visp(ir,nspin*(iatom-1)+2)-visp(ir,nspin*(iatom-1)+1))*0.5*  &
        DSQRT(16.0D0*DATAN(1.0D0))
    
  END DO
END DO

DO  ir = 1,irmd
! cut contributions from outside the MT if recquired
  IF (test('ONLYMT  ') .AND. (ir > ircut(1))) THEN
    rges(ir) = 0
  ELSE
    rges(ir) = rsp(ir)
  END IF
END DO

IF (.NOT.test('ONLYSPH ')) THEN
  DO  j = 1,iend
    lm1p = icleb(j,1)
    lm2p = icleb(j,2)
    lm3p = icleb(j,3)   ! always >= 2 here
    clt = cleb(j)
    
!--->   calculate the non spherically symmetric contribution
    IF (ipan > 1 .AND. lmsp(icell,lm3p) > 0) THEN
      ifun = ifunm(icell,lm3p)
      IF (lm1p == lm2p ) THEN
!             DO 150 IR = IRCUT(1)+1,IRCUT(IPAN)
        DO  ir = irmin,ircut(ipan)
          rges(ir) = rges(ir)+rll(ir,lm2p,lm1p)*cleb(j)*(-1)*  &
              (vins(ir,lm3p,nspin* (iatom-1) + 2) -  &
              vins(ir,lm3p,nspin* (iatom-1) + 1))*0.5
        END DO
      ELSE
!             DO IR = IRCUT(1)+1,IRCUT(IPAN)
        DO ir = irmin,ircut(ipan)
          rges(ir) = rges(ir)+  &
              cleb(j)*(-1)*(vins(ir,lm3p,nspin* (iatom-1) + 2) -  &
              vins(ir,lm3p,nspin* (iatom-1) + 1))*0.5*  &
              (rll(ir,lm2p,lm1p)+rll(ir,lm1p,lm2p))
        END DO
      END IF
    END IF
    
  END DO
END IF

IF (ipan == 1) THEN
  ircutm(0) = 0
  ircutm(1) = irws
ELSE
  DO  i = 0,ipan
    ircutm(i) = ircut(i)
  END DO
END IF

CALL csimpk(rges(:),dens,ipan,ircutm,drdi)

deallocate(rges)
deallocate(rsp)

END SUBROUTINE
