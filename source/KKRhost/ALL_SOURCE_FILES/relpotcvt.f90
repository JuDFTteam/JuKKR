SUBROUTINE relpotcvt(icall,vm2z,zin,rin,drdiin,ircut,vtrel,btrel,  &
    zrel,rmrel,jwsrel,drdirel,r2drdirel,irshift, ipand,irmd,npotd,natypd)
!   ********************************************************************
!   *                                                                  *
!   * driving routine to convert the TB-KKR potential from the non-    *
!   * relativistic representation VM2Z(IRMD,NPOTD), with IPOTD the     *
!   * combined index for ATOM and SPIN to the relativistic one         *
!   *      VTREL =  (VUP+VDN)/2.0D0                                    *
!   *      BTREL =  (VUP-VDN)/2.0D0                                    *
!   *                                                                  *
!   * IMPORTANT 1:  because this routine is called only IF KREL.EQ.0,  *
!   *               the number of spins in VM2Z is always 2!           *
!   * IMPORTANT 2:  so far, only SPHERICAL part implemented            *
!   *                                                                  *
!   * Additionally, for compatibility with the relativistic routines   *
!   * included in the package, VTREL includes the Coulomb term, and    *
!   * the auxiliary arrays                                             *
!   *      ZREL, RMREL, JWSREL, DRDI, R2DRDI and IRSHIFT               *
!   * are created. IRSHIFT(NATYPD) accounts for the shift in the       *
!   * radial mesh, since the first point (sometimes first two points)  *
!   * of VM2Z ( = 0D0 ) are skipped. The relativistic routines require *
!   * an ODD number of radial points (Simpson integration routine)     *
!   *                                                                  *
!   * v.popescu, munich, may 2004                                      *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

!PARAMETER definitions
INTEGER NSPINPOT
PARAMETER (NSPINPOT=2)

! Scalar arguments
INTEGER ICALL,IPAND,IRMD,NPOTD,NATYPD

! Array arguments
DOUBLE PRECISION VM2Z(IRMD,NPOTD)
DOUBLE PRECISION ZIN(NATYPD),RIN(IRMD,NATYPD)
DOUBLE PRECISION DRDIIN(IRMD,NATYPD)
INTEGER IRCUT(0:IPAND,NATYPD)

DOUBLE PRECISION VTREL(IRMD,NATYPD),BTREL(IRMD,NATYPD)
DOUBLE PRECISION DRDIREL(IRMD,NATYPD),R2DRDIREL(IRMD,NATYPD)
DOUBLE PRECISION RMREL(IRMD,NATYPD)
INTEGER IRSHIFT(NATYPD),JWSREL(NATYPD),ZREL(NATYPD)

! Local scalars
DOUBLE PRECISION VDN, VUP
INTEGER IT,IR,IP,IPOT,ISHIFT,JR

! Intrinsic Functions
INTRINSIC NINT

! External Subroutines
EXTERNAL RINIT

! ------------------------------------------------------- INITIALISATION
IF ( icall == 1 ) THEN
  CALL rinit(irmd*natypd,rmrel)
  CALL rinit(irmd*natypd,drdirel)
  CALL rinit(irmd*natypd,r2drdirel)
  DO it = 1,natypd
    jwsrel(it) = 0
    irshift(it) = 0
    zrel(it) = 0
  END DO
END IF
CALL rinit(irmd*natypd,vtrel)
CALL rinit(irmd*natypd,btrel)
! ------------------------------------------------------- INITIALISATION

! *************************************************************** NATYPD
DO it = 1,natypd
! ================================================================ ICALL
!                                       variables require init only once
  IF ( icall == 1 ) THEN
    
! skip first mesh point and also the second if IRCUT(1,IT) = WS-rad odd,
! since JWSREL(IT) must be odd
    
    ishift = 1
    IF ( MOD(ircut(1,it),2) == 1 ) ishift = 2
    ir = 0
! ----------------------------------------------------------------------
    DO jr = 1 + ishift,ircut(1,it)
      ir = ir + 1
      rmrel(ir,it) = rin(jr,it)
      drdirel(ir,it) = drdiin(jr,it)
      r2drdirel(ir,it) = rmrel(ir,it)*rmrel(ir,it) *drdirel(ir,it)
    END DO
! ----------------------------------------------------------------------
    jwsrel(it) = ir
    irshift(it) = ishift
    zrel(it) = nint(zin(it))
  END IF
! ================================================================ ICALL
  ipot = (it-1)*nspinpot + 1
  ishift = irshift(it)
! ----------------------------------------------------------------------
  DO ir = 1,jwsrel(it)
    ip = ir + ishift
    vdn = -2D0*zin(it)/rin(ip,it) + vm2z(ip,ipot)
    vup = -2D0*zin(it)/rin(ip,it) + vm2z(ip,ipot+1)
    vtrel(ir,it) = (vup+vdn)/2.0D0
    btrel(ir,it) = (vup-vdn)/2.0D0
  END DO
! ----------------------------------------------------------------------
END DO
! *************************************************************** NATYPD
END SUBROUTINE relpotcvt
