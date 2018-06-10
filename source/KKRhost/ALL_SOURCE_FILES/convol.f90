! ************************************************************************
SUBROUTINE convol(imt1,irc1,icell,imaxsh,ilm_map,ifunm,lmpot,gsh,  &
    thetas,thesme,z,rfpi,r,vons,vspsmo,lmsp)
! ************************************************************************
!.. Parameters ..
include 'inc.p'
INTEGER LMPOTD
PARAMETER (LMPOTD= (LPOTD+1)**2)
!..
!.. Scalar Arguments ..
DOUBLE PRECISION RFPI,Z
INTEGER ICELL,IMAXSH,IMT1,IRC1,LMPOT
!..
!.. Array Arguments ..
DOUBLE PRECISION GSH(*),R(*),THETAS(IRID,NFUND,*),VONS(IRMD,*), &
                 THESME(IRID,NFUND,*),VSPSMO(IRMD)
INTEGER IFUNM(NATYPD,*),ILM_MAP(NGSHD,3),LMSP(NATYPD,*)
!..
!.. Local Scalars ..
DOUBLE PRECISION ZZOR
INTEGER I,IFUN,IR,IRH,LM,LM1,LM2,LM3
!..
!.. Local Arrays ..
DOUBLE PRECISION VSTORE(IRID,LMPOTD),VSTSME(IRID,LMPOTD)
!..

DO  lm = 1,lmpot
  DO  ir = 1,irc1 - imt1
    vstore(ir,lm) = 0.0D0
    vstsme(ir,lm) = 0.0D0
  END DO
END DO

DO  ir = imt1 + 1,irc1
  zzor = 2.0D0*z/r(ir)*rfpi
  vons(ir,1) = vons(ir,1) - zzor
END DO

DO  i = 1,imaxsh
  lm1 = ilm_map(i,1)
  lm2 = ilm_map(i,2)
  lm3 = ilm_map(i,3)
  IF (lmsp(icell,lm3) > 0) THEN
    ifun = ifunm(icell,lm3)
    DO  ir = imt1 + 1,irc1
      irh = ir - imt1
      vstore(irh,lm1) = vstore(irh,lm1) +  &
          gsh(i)*vons(ir,lm2)*thetas(irh,ifun,icell)
      vstsme(irh,lm1) = vstsme(irh,lm1) +  &
          gsh(i)*vons(ir,lm2)*thesme(irh,ifun,icell)
    END DO
  END IF
END DO

DO  ir = imt1 + 1,irc1
  irh = ir - imt1
  zzor = 2.0D0*z/r(ir)*rfpi
  vons(ir,1) = vstore(irh,1) + zzor
  vspsmo(ir) = (vstsme(irh,1) + zzor) /rfpi
END DO

!     COPY THE PART INSIDE THE MT SPHERE
DO ir = 1,imt1
  vspsmo(ir) = vons(ir,1)/rfpi
END DO

DO  lm = 2,lmpot
  DO  ir = imt1 + 1,irc1
    irh = ir - imt1
    vons(ir,lm) = vstore(irh,lm)
  END DO
END DO

RETURN

END SUBROUTINE convol
