! ************************************************************************
    Subroutine mtzero(lmpot, natyp, conc, nspin, v, vbc, z, r, drdi, imt, &
      ircut, ipan, ntcell, lmsp, ifunm, thetas, irws, eshift, ishift, nshell, &
      lsurf)
      Use mod_datatypes, Only: dp
! ************************************************************************

!     determine muffin tin zero and shift potential to muffin tin zero

!     for spin polarized calculations muffin tin zero is related to
!         the average of the 2 spins

!                                            may,2000 (new version)

!-----------------------------------------------------------------------
      Implicit None
!.. Parameters ..
      Include 'inc.p'
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: eshift, vbc(*)
      Integer :: ishift, lmpot, natyp, nspin
!..
!.. Local Arrays ..
      Real (Kind=dp) :: drdi(irmd, *), conc(natypd), r(irmd, *), &
        thetas(irid, nfund, *), v(irmd, lmpotd, *), z(*)
      Integer :: ifunm(natypd, *), imt(*), ipan(*), ircut(0:ipand, *), &
        irws(*), lmsp(natypd, *), ntcell(*), nshell(0:nsheld)
      Logical :: lsurf
!..
!.. External Subroutines ..
      Real (Kind=dp) :: fpi, rfpi, vav0, vol0, zzor
      Integer :: icell, ifun, ih, imt1, ipan1, ipot, ir, irc1, irh, is, lm
!..
!.. Intrinsic Functions ..
      Real (Kind=dp) :: v1(irmd), v2(irmd), vav1(2), vol1(2)
!..

      Logical :: test, opt
      External :: simp3, simpk, test


      Intrinsic :: atan, sqrt

      fpi = 16.0E0_dp*atan(1.0E0_dp)
      rfpi = sqrt(fpi)
!---  >     muffin tin or atomic sphere calculation
      vav0 = 0.0E0_dp
      vol0 = 0.0E0_dp
      vav1(1) = 0.E0_dp
      vav1(2) = 0.E0_dp
      vol1(1) = 0.E0_dp
      vol1(2) = 0.E0_dp
      Do ih = 1, natyp

        Do ir = 1, irmd
          v1(ir) = 0.0E0_dp
          v2(ir) = 0.0E0_dp
        End Do
        Do is = 1, nspin
          ipot = nspin*(ih-1) + is
          ipan1 = ipan(ih)
          imt1 = imt(ih)

          If (ipan1==1) Then

! (IPAN1.EQ.1)

            irc1 = irws(ih)
            Do ir = imt1, irc1
              v2(ir) = fpi*r(ir, ih)**2
              zzor = 2.0E0_dp*z(ih)/r(ir, ih)
              v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
            End Do
!---  >     full potential calculation
            Call simp3(v1, vav1(is), imt1, irc1, drdi(1,ih))
            Call simp3(v2, vol1(is), imt1, irc1, drdi(1,ih))

          Else



            irc1 = ircut(ipan1, ih)
            icell = ntcell(ih)
            imt1 = imt(ih)
            Do ir = imt1 + 1, irc1
              v2(ir) = r(ir, ih)**2*thetas(ir-imt1, 1, icell)*rfpi
              zzor = 2.0E0_dp*z(ih)/r(ir, ih)
              v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
            End Do
            Do lm = 2, lmpot
              If (lmsp(icell,lm)>0) Then
                ifun = ifunm(icell, lm)

                Do ir = imt1 + 1, irc1
                  irh = ir - imt1
                  v1(ir) = v1(ir) + r(ir, ih)**2*v(ir, lm, ipot)*thetas(irh, &
                    ifun, icell)
                End Do
! (IPAN1.EQ.1)
              End If

            End Do
! SPIN LOOP
            Call simpk(v1, vav1(is), ipan1, ircut(0,ih), drdi(1,ih))
            Call simpk(v2, vol1(is), ipan1, ircut(0,ih), drdi(1,ih))

          End If
!     19.5.99   Nikos
        End Do !     This way it is compatible with old kkr and tb-kkr
        If (nspin==1) Then
          vav1(2) = vav1(1)
          vol1(2) = vol1(1)
        End If

! added 10.11.99 to fix vbc


        If (lsurf .And. (ih==1)) Write (1337, *) &
          'Vacancies are ignored for VBC'
!---  > shift potential to muffin tin zero
        If (lsurf .And. (z(ih)<1.E0_dp)) Cycle
        vav0 = vav0 + conc(ih)*nshell(ih)*(vav1(1)+vav1(2))/2.E0_dp
        vol0 = vol0 + conc(ih)*nshell(ih)*(vol1(1)+vol1(2))/2.E0_dp
      End Do
      If (.Not. (opt('DECIMATE'))) Then
        vbc(1) = 0.0E0_dp
        If (abs(vav0)>1E-10_dp) vbc(1) = -vav0/vol0
        If (ishift>0) vbc(1) = vbc(1) + eshift
      End If

      Write (1337, Fmt=100) vol0, vav0, vbc(1)
      vbc(2) = vbc(1)


! ************************************************************************
      Do is = 1, nspin
        Do ih = 1, natyp
          ipot = nspin*(ih-1) + is
          Do ir = 1, ircut(ipan(ih), ih)
            v(ir, 1, ipot) = v(ir, 1, ipot) + rfpi*vbc(is)
          End Do
! ************************************************************************
        End Do
      End Do

      Return
!     determine muffin tin zero and shift potential to muffin tin zero
100   Format ('  VOL INT.', F16.9, '  VAV INT.', F16.9, '  VMT ZERO', F16.9)
110   Format ('  ATOM ', I4, ' VMT ZERO :', F16.9)
    End Subroutine
