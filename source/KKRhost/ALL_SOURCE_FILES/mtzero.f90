! ************************************************************************
SUBROUTINE mtzero(lmpot,natyp,conc,nspin,v,vbc,z,r,drdi,imt,ircut,  &
    ipan,ntcell,lmsp,ifunm,thetas,irws,eshift, ishift,nshell,lsurf)
! ************************************************************************

!     determine muffin tin zero and shift potential to muffin tin zero

!     for spin polarized calculations muffin tin zero is related to
!         the average of the 2 spins

!                                            may,2000 (new version)

!-----------------------------------------------------------------------
implicit none
!.. Parameters ..
include 'inc.p'
INTEGER LMPOTD
PARAMETER (LMPOTD= (LPOTD+1)**2)
!..
!.. Scalar Arguments ..
DOUBLE PRECISION ESHIFT,VBC(*)
INTEGER ISHIFT,LMPOT,NATYP,NSPIN
!..
!.. Array Arguments ..
DOUBLE PRECISION  &
     DRDI(IRMD,*),CONC(NATYPD), &
     R(IRMD,*), &
     THETAS(IRID,NFUND,*), &
     V(IRMD,LMPOTD,*), &
     Z(*)
INTEGER &
     IFUNM(NATYPD,*),IMT(*),IPAN(*),IRCUT(0:IPAND,*),IRWS(*), &
     LMSP(NATYPD,*), &
     NTCELL(*),NSHELL(0:NSHELD)
LOGICAl LSURF
!..
!.. Local Scalars ..
DOUBLE PRECISION FPI,RFPI,VAV0,VOL0,ZZOR
INTEGER ICELL,IFUN,IH,IMT1,IPAN1,IPOT,IR,IRC1,IRH,IS,LM
!..
!.. Local Arrays ..
DOUBLE PRECISION V1(IRMD),V2(IRMD),VAV1(2),VOL1(2)
!..
!.. External Subroutines ..
LOGICAL TEST,OPT
EXTERNAL SIMP3,SIMPK,TEST
!..
!.. Intrinsic Functions ..
INTRINSIC ATAN,SQRT
!..
fpi = 16.0D0*ATAN(1.0D0)
rfpi = SQRT(fpi)

vav0 = 0.0D0
vol0 = 0.0D0
vav1(1) = 0.d0
vav1(2) = 0.d0
vol1(1) = 0.d0
vol1(2) = 0.d0
DO  ih = 1,natyp
  
  DO  ir = 1,irmd
    v1(ir) = 0.0D0
    v2(ir) = 0.0D0
  END DO
  DO is = 1,nspin
    ipot = nspin* (ih-1) + is
    ipan1 = ipan(ih)
    imt1 = imt(ih)
    
    IF (ipan1 == 1) THEN
      
!---  >     muffin tin or atomic sphere calculation
      
      irc1 = irws(ih)
      DO  ir = imt1,irc1
        v2(ir) = fpi*r(ir,ih)**2
        zzor = 2.0D0*z(ih)/r(ir,ih)
        v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
      END DO
      
      CALL simp3(v1,vav1(is),imt1,irc1,drdi(1,ih))
      CALL simp3(v2,vol1(is),imt1,irc1,drdi(1,ih))
      
    ELSE                 ! (IPAN1.EQ.1)
      
!---  >     full potential calculation
      
      irc1 = ircut(ipan1,ih)
      icell = ntcell(ih)
      imt1 = imt(ih)
      DO  ir = imt1 + 1,irc1
        v2(ir) = r(ir,ih)**2*thetas(ir-imt1,1,icell)*rfpi
        zzor = 2.0D0*z(ih)/r(ir,ih)
        v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
      END DO
      DO  lm = 2,lmpot
        IF (lmsp(icell,lm) > 0) THEN
          ifun = ifunm(icell,lm)
          
          DO  ir = imt1 + 1,irc1
            irh = ir - imt1
            v1(ir) = v1(ir) + r(ir,ih)**2*v(ir,lm,ipot)*  &
                thetas(irh,ifun,icell)
          END DO
          
        END IF
        
      END DO
      
      CALL simpk(v1,vav1(is),ipan1,ircut(0,ih),drdi(1,ih))
      CALL simpk(v2,vol1(is),ipan1,ircut(0,ih),drdi(1,ih))
      
    END IF               ! (IPAN1.EQ.1)
    
  END DO                  ! SPIN LOOP
  IF (nspin == 1) THEN
    vav1(2) = vav1(1)
    vol1(2) = vol1(1)
  END IF
  
  
!     19.5.99   Nikos
!     This way it is compatible with old kkr and tb-kkr
  IF(lsurf.AND.(ih == 1)) WRITE(1337,*) 'Vacancies are ignored for VBC'
  
  IF (lsurf.AND.(z(ih) < 1.d0)) CYCLE
  vav0 = vav0 + conc(ih)*nshell(ih)*(vav1(1)+vav1(2))/2.d0
  vol0 = vol0 + conc(ih)*nshell(ih)*(vol1(1)+vol1(2))/2.d0
END DO
IF (.NOT.(opt('DECIMATE'))) THEN ! added 10.11.99 to fix vbc
  vbc(1) = 0.0D0
  IF (ABS(vav0) > 1D-10) vbc(1) = -vav0/vol0
  IF (ishift > 0) vbc(1) = vbc(1) + eshift
END IF

WRITE (1337,FMT=9000) vol0,vav0,vbc(1)
vbc(2) = vbc(1)

!---  > shift potential to muffin tin zero

DO  is = 1,nspin
  DO  ih = 1,natyp
    ipot = nspin* (ih-1) + is
    DO  ir = 1,ircut(ipan(ih),ih)
      v(ir,1,ipot) = v(ir,1,ipot) + rfpi*vbc(is)
    END DO
    
  END DO
END DO

RETURN

9000 FORMAT ('  VOL INT.',f16.9,'  VAV INT.',f16.9,'  VMT ZERO',f16.9)
9010 FORMAT ('  ATOM ',i4,' VMT ZERO :',f16.9)
END SUBROUTINE mtzero
