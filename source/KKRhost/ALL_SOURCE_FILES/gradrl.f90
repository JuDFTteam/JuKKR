SUBROUTINE gradrl(nspin,mesh,l1max,dx,rhol,rv,drdi,ipan,ipand,  &
        ircut,drrl,ddrrl,drrul,ddrrul,irmd,lmpotd)
!------------------------------------------------------------------
!gradient of rl with rl defined by charge density=sum(rl*ylm).
!mesh,l1max: max of mesh and l+1.
!IRMD,LMPOTD: maxima of corresponding dimension parameters.
!drrl=d(rl)/dr, ddrrl=d(drrl)/dr, drrul=d(rl-up)/dr,
!ztal: zeta for each l-component necessary to get down-components.
!------------------------------------------------------------------
!------------------------------------------------------------------
use mod_types, only: t_inc
IMPLICIT NONE
!.. Parameters ..
DOUBLE PRECISION ZERO,ZERO1
PARAMETER (ZERO=0.D0,ZERO1=1.D-12)
!..
!.. Scalar Arguments ..
DOUBLE PRECISION DX
INTEGER IPAN,IPAND,IRMD,L1MAX,LMPOTD,MESH,NSPIN
!..
!.. Array Arguments ..
DOUBLE PRECISION DDRRL(IRMD,LMPOTD),DDRRUL(IRMD,LMPOTD), &
                 DRDI(IRMD),DRRL(IRMD,LMPOTD),DRRUL(IRMD,LMPOTD), &
                 RHOL(IRMD,2,LMPOTD),RV(IRMD)
INTEGER IRCUT(0:IPAND)
!..
!.. Local Scalars ..
DOUBLE PRECISION CHGDEN,PI,R2,S4,SPIDEN
INTEGER I1,IEN,IP,IR,IST,LLMAX
!..
!.. Local Arrays ..
DOUBLE PRECISION DRDI2(IRMD),RL1(IRMD),RL1UDM(IRMD),ZTAL(IRMD)
!..
!.. External Subroutines ..
EXTERNAL GRADR
!..
!.. Intrinsic Functions ..
INTRINSIC ABS,ACOS,SQRT
!..
!------------------------------------------------------------------
pi = COS(-1.d0)
s4 = SQRT(4.d0*pi)
llmax = l1max*l1max

DO  ip = 1,ipan
  ist = ircut(ip-1) + 1
  ien = ircut(ip)
  IF(t_inc%i_write>0) WRITE (1337,FMT=9010) ip,ist,ien
  IF (ip == 1) THEN
    DO  ir = ist,ien
      drdi2(ir) = dx
    END DO
  ELSE
    DO  ir = ist,ien
      drdi2(ir) = zero
    END DO
  END IF
END DO


DO  i1 = 1,llmax
  
  IF (nspin == 1) GO TO 50
  
  DO  ir = 2,mesh
    r2 = rv(ir)*rv(ir)
    chgden = rhol(ir,1,i1) + rhol(ir,2,i1)
    spiden = rhol(ir,2,i1) - rhol(ir,1,i1)
    IF (ABS(chgden) >= zero1) THEN
      rl1(ir) = chgden
      ztal(ir) = spiden/chgden
    ELSE
      rl1(ir) = zero
      ztal(ir) = zero
    END IF
  END DO
  
  GO TO 70
  
  
  50   CONTINUE
  DO  ir = 2,mesh
    r2 = rv(ir)*rv(ir)
    rl1(ir) = rhol(ir,1,i1) + rhol(ir,2,i1)
    ztal(ir) = zero
  END DO
  
  70   CONTINUE
  
  rl1(1) = rl1(2)
  ztal(1) = ztal(2)
  
  
  DO  ip = 1,ipan
    ist = ircut(ip-1) + 1
    ien = ircut(ip)
    
    CALL gradr(nspin,ist,ien,1.d0,drdi,drdi2,rl1,ztal,drrl(1,i1),  &
        ddrrl(1,i1),drrul(1,i1),ddrrul(1,i1),rl1udm,irmd)
    
    IF (ip == 1) THEN
      DO  ir = 1,4
        drrl(ir,i1) = drrl(5,i1)
        ddrrl(ir,i1) = ddrrl(5,i1)
        drrul(ir,i1) = drrul(5,i1)
        ddrrul(ir,i1) = ddrrul(5,i1)
      END DO
    END IF
    
    IF (nspin == 1) THEN
      DO  ir = ist,ien
        drrul(ir,i1) = drrl(ir,i1)/2.d0
        ddrrul(ir,i1) = ddrrl(ir,i1)/2.d0
      END DO
    END IF
    
  END DO
  
END DO

RETURN
9000 FORMAT (1X,' l1max=',i5,' mesh=',i5,'nspi=',i5,' ipan=',i5)
9010 FORMAT (1X,'  ip ist ien',3I5)
END SUBROUTINE gradrl
