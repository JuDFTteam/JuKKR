! ************************************************************************
    Subroutine convol(imt1, irc1, icell, imaxsh, ilm_map, ifunm, lmpot, gsh, &
      thetas, thesme, z, rfpi, r, vons, vspsmo, lmsp)
      Use mod_datatypes, Only: dp
! ************************************************************************
!.. Parameters ..
      Include 'inc.p'
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: rfpi, z
      Integer :: icell, imaxsh, imt1, irc1, lmpot
!..
!.. Local Arrays ..
      Real (Kind=dp) :: gsh(*), r(*), thetas(irid, nfund, *), vons(irmd, *), &
        thesme(irid, nfund, *), vspsmo(irmd)
      Integer :: ifunm(natypd, *), ilm_map(ngshd, 3), lmsp(natypd, *)
!..

      Real (Kind=dp) :: zzor
      Integer :: i, ifun, ir, irh, lm, lm1, lm2, lm3


      Real (Kind=dp) :: vstore(irid, lmpotd), vstsme(irid, lmpotd)


      Do lm = 1, lmpot
        Do ir = 1, irc1 - imt1
          vstore(ir, lm) = 0.0E0_dp
          vstsme(ir, lm) = 0.0E0_dp
        End Do
      End Do
!     COPY THE PART INSIDE THE MT SPHERE
      Do ir = imt1 + 1, irc1
        zzor = 2.0E0_dp*z/r(ir)*rfpi
        vons(ir, 1) = vons(ir, 1) - zzor
      End Do

      Do i = 1, imaxsh
        lm1 = ilm_map(i, 1)
        lm2 = ilm_map(i, 2)
        lm3 = ilm_map(i, 3)
        If (lmsp(icell,lm3)>0) Then
          ifun = ifunm(icell, lm3)
          Do ir = imt1 + 1, irc1
            irh = ir - imt1
            vstore(irh, lm1) = vstore(irh, lm1) + gsh(i)*vons(ir, lm2)*thetas( &
              irh, ifun, icell)
            vstsme(irh, lm1) = vstsme(irh, lm1) + gsh(i)*vons(ir, lm2)*thesme( &
              irh, ifun, icell)
          End Do
        End If
      End Do

      Do ir = imt1 + 1, irc1
        irh = ir - imt1
        zzor = 2.0E0_dp*z/r(ir)*rfpi
        vons(ir, 1) = vstore(irh, 1) + zzor
        vspsmo(ir) = (vstsme(irh,1)+zzor)/rfpi
      End Do

! ************************************************************************
      Do ir = 1, imt1
        vspsmo(ir) = vons(ir, 1)/rfpi
      End Do
! ************************************************************************
      Do lm = 2, lmpot
        Do ir = imt1 + 1, irc1
          irh = ir - imt1
          vons(ir, lm) = vstore(irh, lm)
        End Do
      End Do
!.. Parameters ..
      Return
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
    End Subroutine
