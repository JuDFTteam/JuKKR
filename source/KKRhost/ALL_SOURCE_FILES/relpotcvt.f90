subroutine relpotcvt(icall, vm2z, zin, rin, drdiin, ircut, vtrel, btrel, zrel, &
  rmrel, jwsrel, drdirel, r2drdirel, irshift, ipand, irmd, npotd, natypd)
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

  implicit none

!PARAMETER definitions
  integer :: nspinpot
  parameter (nspinpot=2)

! Scalar arguments
  integer :: icall, ipand, irmd, npotd, natypd

! Array arguments
  double precision :: vm2z(irmd, npotd)
  double precision :: zin(natypd), rin(irmd, natypd)
  double precision :: drdiin(irmd, natypd)
  integer :: ircut(0:ipand, natypd)

  double precision :: vtrel(irmd, natypd), btrel(irmd, natypd)
  double precision :: drdirel(irmd, natypd), r2drdirel(irmd, natypd)
  double precision :: rmrel(irmd, natypd)
  integer :: irshift(natypd), jwsrel(natypd), zrel(natypd)

! Local scalars
  double precision :: vdn, vup
  integer :: it, ir, ip, ipot, ishift, jr

! Intrinsic Functions
  intrinsic :: nint

! External Subroutines
  external :: rinit

! ------------------------------------------------------- INITIALISATION
  if (icall==1) then
    call rinit(irmd*natypd, rmrel)
    call rinit(irmd*natypd, drdirel)
    call rinit(irmd*natypd, r2drdirel)
    do it = 1, natypd
      jwsrel(it) = 0
      irshift(it) = 0
      zrel(it) = 0
    end do
  end if
  call rinit(irmd*natypd, vtrel)
  call rinit(irmd*natypd, btrel)
! ------------------------------------------------------- INITIALISATION

! *************************************************************** NATYPD
  do it = 1, natypd
! ================================================================ ICALL
!                                       variables require init only once
    if (icall==1) then

! skip first mesh point and also the second if IRCUT(1,IT) = WS-rad odd,
! since JWSREL(IT) must be odd

      ishift = 1
      if (mod(ircut(1,it),2)==1) ishift = 2
      ir = 0
! ----------------------------------------------------------------------
      do jr = 1 + ishift, ircut(1, it)
        ir = ir + 1
        rmrel(ir, it) = rin(jr, it)
        drdirel(ir, it) = drdiin(jr, it)
        r2drdirel(ir, it) = rmrel(ir, it)*rmrel(ir, it)*drdirel(ir, it)
      end do
! ----------------------------------------------------------------------
      jwsrel(it) = ir
      irshift(it) = ishift
      zrel(it) = nint(zin(it))
    end if
! ================================================================ ICALL
    ipot = (it-1)*nspinpot + 1
    ishift = irshift(it)
! ----------------------------------------------------------------------
    do ir = 1, jwsrel(it)
      ip = ir + ishift
      vdn = -2d0*zin(it)/rin(ip, it) + vm2z(ip, ipot)
      vup = -2d0*zin(it)/rin(ip, it) + vm2z(ip, ipot+1)
      vtrel(ir, it) = (vup+vdn)/2.0d0
      btrel(ir, it) = (vup-vdn)/2.0d0
    end do
! ----------------------------------------------------------------------
  end do
! *************************************************************** NATYPD
end subroutine
