SUBROUTINE calc_rho_ll_ss_lmdos(rll,ircut,ipan,icell,  &
        thetas,cleb,icleb,iend,ifunm,lmsp,irws,drdi,dens,  &
        lmdos)
IMPLICIT NONE

include 'inc.p'
INTEGER          LMMAXD
PARAMETER        (LMMAXD= (LMAXD+1)**2)
INTEGER          LMPOTD
PARAMETER        (LMPOTD= (LPOTD+1)**2)
INTEGER          IRMIND
PARAMETER        (IRMIND=IRMD-IRNSD)
!..
!.. Scalar Arguments ..
INTEGER          IEND,IRWS,LMDOS
!..
!.. Array Arguments ..
DOUBLE COMPLEX   RLL(IRMD,LMMAXD,LMMAXD), &  ! non-sph. eigen states of single pot 
                 DENS
DOUBLE PRECISION CLEB(*), &
                 THETAS(IRID,NFUND,*), &
                 DRDI(IRMD)                            ! derivative dr/di
INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD), &
                 LMSP(NATYPD,*),IRCUT(0:IPAND),IPAN, &
                 ICELL,IFUN


! local variables

DOUBLE PRECISION             ::   C0LL
DOUBLE COMPLEX               ::   CLT 
DOUBLE COMPLEX, ALLOCATABLE  ::   RSP(:),RGES(:)
integer                      ::   LM1P,LM2P,LM3P,IR,J,I
INTEGER                      ::   IRCUTM(0:IPAND)

!     ..
!  ---> first calculate only the spherically symmetric contribution
!       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
!       multiplied with the shape functions...

!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)

allocate(rges(irmd))
allocate(rsp(irmd))

!      WRITE(6,*) "In rho ll"

c0ll = 1.0D0/DSQRT(16.0D0*DATAN(1.0D0))
rsp=0D0
rges=0D0

DO ir=1,irmd
  rsp(ir)=rsp(ir)+rll(ir,lmdos,lmdos)
END DO

DO  ir = 1,ircut(ipan)
  rges(ir) = rsp(ir)
END DO

IF (ipan > 1) THEN
  DO  ir = ircut(1)+1,ircut(ipan)
    rges(ir) = rsp(ir)*c0ll*thetas(ir-ircut(1),1,icell)
  END DO
END IF

!      STOP " "
!      WRITE(6,*) "IRCUT(1)",IRCUT(1)
!      WRITE(6,*) "IRCUT(IPAN)",IRCUT(IPAN)
!      WRITE(6,*) "IRCUT(IPAN)-IRMIND",IRCUT(IPAN)-IRMIND
!      WRITE(6,*) "IRMIND",IRMIND

DO  j = 1,iend
  lm1p = icleb(j,1)
  lm2p = icleb(j,2)
  lm3p = icleb(j,3)
  clt = cleb(j)
  
!---> calculate the non spherically symmetric contribution
  
  IF (ipan > 1 .AND. lmsp(icell,lm3p) > 0) THEN
    ifun = ifunm(icell,lm3p)
!          WRITE(156,*) "IFUN",IFUN
    IF (lm1p == lm2p.AND.lm1p == lmdos ) THEN
      DO  ir = ircut(1)+1,ircut(ipan)
        rges(ir) = rges(ir)+rll(ir,lm2p,lm1p)*  &
            cleb(j)*thetas(ir-ircut(1),ifun,icell)
      END DO
!          ELSE
!            DO IR = IRCUT(1)+1,IRCUT(IPAN)
!              RGES(IR) = RGES(IR)+
!     +             CLEB(J)*THETAS(IR-IRCUT(1),IFUN,ICELL)*
!     +          (RLL(IR,LM2P,LM1P)+RLL(IR,LM1P,LM2P))
!            END DO
    END IF
  END IF
  
END DO

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
