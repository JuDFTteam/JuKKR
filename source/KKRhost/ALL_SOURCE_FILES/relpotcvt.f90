    Subroutine relpotcvt(icall, vm2z, zin, rin, drdiin, ircut, vtrel, btrel, &
      zrel, rmrel, jwsrel, drdirel, r2drdirel, irshift, ipand, irmd, npotd, &
      natypd)
      Use mod_datatypes, Only: dp
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

      Implicit None

!PARAMETER definitions
      Integer :: nspinpot
      Parameter (nspinpot=2)

! Scalar arguments
      Integer :: icall, ipand, irmd, npotd, natypd

! Array arguments
      Real (Kind=dp) :: vm2z(irmd, npotd)
      Real (Kind=dp) :: zin(natypd), rin(irmd, natypd)
      Real (Kind=dp) :: drdiin(irmd, natypd)
      Integer :: ircut(0:ipand, natypd)

      Real (Kind=dp) :: vtrel(irmd, natypd), btrel(irmd, natypd)
      Real (Kind=dp) :: drdirel(irmd, natypd), r2drdirel(irmd, natypd)
      Real (Kind=dp) :: rmrel(irmd, natypd)
      Integer :: irshift(natypd), jwsrel(natypd), zrel(natypd)

! Local scalars
      Real (Kind=dp) :: vdn, vup
      Integer :: it, ir, ip, ipot, ishift, jr

! Intrinsic Functions
      Intrinsic :: nint

! External Subroutines
      External :: rinit

! ------------------------------------------------------- INITIALISATION
      If (icall==1) Then
        Call rinit(irmd*natypd, rmrel)
        Call rinit(irmd*natypd, drdirel)
        Call rinit(irmd*natypd, r2drdirel)
        Do it = 1, natypd
          jwsrel(it) = 0
          irshift(it) = 0
          zrel(it) = 0
        End Do
      End If
      Call rinit(irmd*natypd, vtrel)
      Call rinit(irmd*natypd, btrel)
! ------------------------------------------------------- INITIALISATION

! *************************************************************** NATYPD
      Do it = 1, natypd
! ================================================================ ICALL
!                                       variables require init only once
        If (icall==1) Then

! skip first mesh point and also the second if IRCUT(1,IT) = WS-rad odd,
! since JWSREL(IT) must be odd

          ishift = 1
          If (mod(ircut(1,it),2)==1) ishift = 2
          ir = 0
! ----------------------------------------------------------------------
          Do jr = 1 + ishift, ircut(1, it)
            ir = ir + 1
            rmrel(ir, it) = rin(jr, it)
            drdirel(ir, it) = drdiin(jr, it)
            r2drdirel(ir, it) = rmrel(ir, it)*rmrel(ir, it)*drdirel(ir, it)
          End Do
! ----------------------------------------------------------------------
          jwsrel(it) = ir
          irshift(it) = ishift
          zrel(it) = nint(zin(it))
        End If
! ================================================================ ICALL
        ipot = (it-1)*nspinpot + 1
        ishift = irshift(it)
! ----------------------------------------------------------------------
        Do ir = 1, jwsrel(it)
          ip = ir + ishift
          vdn = -2E0_dp*zin(it)/rin(ip, it) + vm2z(ip, ipot)
          vup = -2E0_dp*zin(it)/rin(ip, it) + vm2z(ip, ipot+1)
          vtrel(ir, it) = (vup+vdn)/2.0E0_dp
          btrel(ir, it) = (vup-vdn)/2.0E0_dp
        End Do
! ----------------------------------------------------------------------
      End Do
! *************************************************************** NATYPD
    End Subroutine
