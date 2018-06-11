! LLY Lloyd  &
    Subroutine renorm_lly(cdos_lly, ielast, nspin, natyp, cden, lmaxp1, conc, &
      iestart, ieend, wez, ircut, ipan, ez, zat, rho2ns, r2nef, denef, &
      denefat, espv)
      Use mod_datatypes, Only: dp
! Renormalize the valence charge according to Lloyd's formula.
! Find renormalization constant per energy, then renormalize
! charge/atom/energy, then integrate over energies to find
! the renormalized charge/atom. Use it to renormalize the density.
! Phivos Mavropoulos, July 2014
      Implicit None
      Include 'inc.p'
      Integer :: lmaxd1
      Parameter (lmaxd1=lmaxd+1)
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
      Integer :: npotd
      Parameter (npotd=(2*(krel+korbit)+(1-(krel+korbit))*nspind)*natypd)
! Concentration (for cpa)
      Integer :: lmaxp1, natyp, nspin
      Integer :: iestart, ieend, ielast ! Non-renormalized density per atom (density=-cden/pi)
      Integer :: ircut(0:ipand, natypd), ipan(natypd) ! DOS according to Lloyd's formula
      Real (Kind=dp) :: conc(natypd) ! Input/Output:
      Complex (Kind=dp) :: cden(0:lmaxd1, ielast, npotd) ! Internal:
      Complex (Kind=dp) :: cdos_lly(iemxd, nspind) ! 1: charge renormalization per atom (energy-integrated)
      Complex (Kind=dp) :: wez(iemxd), ez(iemxd)
      Real (Kind=dp) :: zat(natypd)
! 2: same for spin moment
      Real (Kind=dp) :: rho2ns(irmd, lmpotd, natypd, 2)
      Real (Kind=dp) :: r2nef(irmd, lmpotd, natypd, 2)
      Real (Kind=dp) :: denef, denefat(natypd)
      Real (Kind=dp) :: espv(0:lmaxd1, npotd)
!  Density from local summation
      Integer :: ll, ie, i1, ispin, ipot, spindegen, irc1, signsp, idim
      Real (Kind=dp) :: renorm_at(natypd, 2) !      and from Lloyd's formula
! Renormalization constant for charge and spin density
      Complex (Kind=dp) :: cdos_loc(iemxd, (1+krel)*nspind) ! Atomic charge per spin (local summation and renormalized)
      Complex (Kind=dp) :: cdos_locvc(iemxd, (1+krel)*nspind)
      Real (Kind=dp) :: cren(iemxd, 2)
      Real (Kind=dp) :: charge(natypd, 2), charge_lly(natypd, 2)
      Complex (Kind=dp) :: chadd(iemxd, natypd, nspind), cdos_add ! Integration step for charge/atom/spin
      Complex (Kind=dp) :: qlly(2), qstar(2)
      Real (Kind=dp) :: sum0(2), sum1(2)
      Complex (Kind=dp) :: czero
      Real (Kind=dp) :: pi
      Logical :: opt, test
      External :: opt, test
! Spin degeneracy, 2 if nspin=1, 1 if nspin=2

      czero = (0.E0_dp, 0.E0_dp)
      pi = 4.E0_dp*atan(1.E0_dp)

      spindegen = 3 - nspin ! First find renormalization factor per energy and atomic charges
! Factor 1/pi included in Wez
      cren(:, :) = 0E0_dp
      renorm_at(:, :) = 1.E0_dp
      charge_lly(:, :) = 0.E0_dp
      charge(:, :) = 0.E0_dp
      qlly(:) = czero
      qstar(:) = czero
! Complex charge
! I1=1,NATYP
      cdos_loc = czero
      cdos_locvc = czero
      chadd = czero
      Do ie = iestart, ieend
        Do ispin = 1, nspin
          Do i1 = 1, natyp
            ipot = (i1-1)*nspin + ispin
            cdos_add = czero
            Do ll = 0, lmaxp1
              cdos_add = cdos_add + conc(i1)*cden(ll, ie, ipot) ! ISPIN = 1,NSPIN
            End Do
            If (zat(i1)>1E-06_dp) Then
              cdos_loc(ie, ispin) = cdos_loc(ie, ispin) + cdos_add
            Else
              cdos_locvc(ie, ispin) = cdos_locvc(ie, ispin) + cdos_add
            End If
            chadd(ie, i1, ispin) = wez(ie)*cdos_add ! IE = IESTART,IEEND
            charge(i1, ispin) = charge(i1, ispin) + aimag(chadd(ie,i1,ispin))/ &
              real(nspin, kind=dp)
          End Do ! Now the locally-summed charge/energy is in cdos_loc, charge/energy/atom in chadd
          cdos_loc(ie, ispin) = -cdos_loc(ie, ispin)/pi
          cdos_locvc(ie, ispin) = -cdos_locvc(ie, ispin)/pi
        End Do ! Renormalization factor per energy:
      End Do ! Apply to DOS of each atom:
! ISPIN = 1,NSPIN
      If (.Not. opt('NEWSOSOL')) Then
        Do ie = iestart, ieend
          Do ispin = 1, nspin
! IE = IESTART,IEEND
            cren(ie, ispin) = aimag((cdos_lly(ie,ispin)-cdos_locvc(ie, &
              ispin))*wez(ie))/aimag(cdos_loc(ie,ispin)*wez(ie))
! Renormalization factor per energy:
            Do i1 = 1, natypd
              If (zat(i1)>1E-06_dp) Then
                charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
                  cren(ie, ispin)*aimag(chadd(ie,i1,ispin))/ &
                  real(nspin, kind=dp)
              Else
                charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
                  aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
              End If
            End Do
          End Do ! Apply to DOS of each atom:
        End Do
      Else
        Do ie = iestart, ieend
! IE = IESTART,IEEND
          cren(ie, 1) = aimag((cdos_lly(ie,1)-cdos_locvc(ie,1)-cdos_locvc(ie, &
            2))*wez(ie))/aimag((cdos_loc(ie,1)+cdos_loc(ie,2))*wez(ie))
! add term from sum from l>lmax to infinity
          Do ispin = 1, nspin
            Do i1 = 1, natypd
              If (zat(i1)>1E-06_dp) Then
                charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
                  cren(ie, 1)*aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
              Else
                charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
                  aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
              End If
            End Do
          End Do
        End Do !            DO I1=1,NATYPD
      End If
!             DO ISPIN=1,NSPIN
!             CHARGE_LLY(I1,ISPIN)=CHARGE_LLY(I1,ISPIN)-DIMAG(CDOS2(I1))
!             ENDDO
!            ENDDO


! Now apply renormalization to energy-integrated density
! If spins are coupled, then only charge density
      If (nspin==1 .Or. opt('NEWSOSOL')) cren(:, 2) = cren(:, 1)
! Index of outmost radial point


      If (nspin==1) Then
        Do i1 = 1, natyp
          If (charge(i1,1)>0) Then
            renorm_at(i1, 1) = charge_lly(i1, 1)/charge(i1, 1)
          Else
            renorm_at(i1, 1) = 1.0E0_dp
          End If
          renorm_at(i1, 2) = renorm_at(i1, 1)
          irc1 = ircut(ipan(i1), i1)
          rho2ns(1:irc1, 1:lmpotd, i1, 1) = rho2ns(1:irc1, 1:lmpotd, i1, 1)* &
            renorm_at(i1, 1)
          r2nef(1:irc1, 1:lmpotd, i1, 1) = r2nef(1:irc1, 1:lmpotd, i1, 1)* &
            renorm_at(i1, 1)
        End Do
! First decouple charge and spin density to the density of each channels
      Else
! Index of outmost radial point
! Second merge density of each channels to charge and spin density
        idim = irmd*lmpotd*natypd
        Call daxpy(idim, 1.0E0_dp, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
        Call dscal(idim, 0.5E0_dp, rho2ns(1,1,1,2), 1)
        Call daxpy(idim, -1.0E0_dp, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)

        Call daxpy(idim, 1.0E0_dp, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
        Call dscal(idim, 0.5E0_dp, r2nef(1,1,1,2), 1)
        Call daxpy(idim, -1.0E0_dp, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
        Do i1 = 1, natyp
          irc1 = ircut(ipan(i1), i1) 
          Do ispin = 1, nspin
            renorm_at(i1, ispin) = charge_lly(i1, ispin)/charge(i1, ispin)
            rho2ns(1:irc1, 1:lmpotd, i1, ispin) = rho2ns(1:irc1, 1:lmpotd, i1, &
              ispin)*renorm_at(i1, ispin)
            r2nef(1:irc1, 1:lmpotd, i1, ispin) = r2nef(1:irc1, 1:lmpotd, i1, &
              ispin)*renorm_at(i1, ispin)
          End Do
        End Do
! calculate density at Fermi level
        Call dscal(idim, 2.0E0_dp, rho2ns(1,1,1,1), 1)
        Call daxpy(idim, -0.5E0_dp, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
        Call daxpy(idim, 1.0E0_dp, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)
! LL
        Call dscal(idim, 2.0E0_dp, r2nef(1,1,1,1), 1)
        Call daxpy(idim, -0.5E0_dp, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
        Call daxpy(idim, 1.0E0_dp, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
      End If
! ISPIN
! I1
      denef = 0E0_dp
      Do i1 = 1, natyp
        denefat(i1) = 0E0_dp
        Do ispin = 1, nspin
          ipot = (i1-1)*nspin + ispin
          Do ll = 0, lmaxp1
            If (zat(i1)>1E-06_dp) Then
              denefat(i1) = denefat(i1) - 2.0E0_dp*conc(i1)*cren(ielast, ispin &
                )*aimag(cden(ll,ielast,ipot))/pi/real(nspin, kind=dp)
            Else
              denefat(i1) = denefat(i1) - 2.0E0_dp*conc(i1)*aimag(cden(ll, &
                ielast,ipot))/pi/real(nspin, kind=dp)
            End If
            espv(ll, ipot) = 0E0_dp
            If (zat(i1)>1E-06_dp) Then
              Do ie = 1, ielast
                espv(ll, ipot) = espv(ll, ipot) + cren(ie, ispin)*aimag(ez(ie) &
                  *cden(ll,ie,ipot)*wez(ie)/real(nspin,kind=dp))
              End Do
            Else
              Do ie = 1, ielast
                espv(ll, ipot) = espv(ll, ipot) + aimag(ez(ie)*cden(ll,ie,ipot &
                  )*wez(ie)/real(nspin,kind=dp))
              End Do
            End If
          End Do ! Write out renormalization factors
        End Do
        denef = denef + denefat(i1)
      End Do
! -1,+1 for spin down,up (ispin=1,2)
      Write (1337, *) 'Information on renormalization by Lloyds formula'
      Write (1337, *) 'RENORM_LLY: Complex renormalization factor per energy:'
      Write (1337, Fmt='(A5,2A32)') 'IE', 'Spin 1 (down)           ', &
        'Spin 2 (up)           '
      Do ie = iestart, ieend
        Write (1337, Fmt='(I5,4F16.12)') ie, (cren(ie,ispin), ispin=1, nspin)
      End Do
      Write (1337, *) 'RENORM_LLY: renormalization factor per atom:'
      Write (1337, Fmt='(A5,2A16)') 'IAT', 'Spin down', 'Spin up'
      Do i1 = 1, natyp
        Write (1337, Fmt='(I5,2E17.9)') i1, (renorm_at(i1,ispin), ispin=1, &
          nspin)
      End Do
      Write (1337, *) 'RENORM_LLY: Renormalized charge per atom:'
      Write (1337, Fmt='(A5,2A16)') 'IAT', 'Spin down', 'Spin up'
      Do i1 = 1, natyp
        Write (1337, Fmt='(I5,2F16.12)') i1, (charge_lly(i1,ispin), ispin=1, &
          nspin)
      End Do
      sum0(:) = 0.E0_dp
      sum1(:) = 0.E0_dp
      Do ispin = 1, nspin
        signsp = 2*ispin - 3
        If (nspin==1) signsp = 1
        Do i1 = 1, natyp
          sum0(ispin) = sum0(ispin) + signsp*conc(i1)*charge(i1, ispin)
          sum1(ispin) = sum1(ispin) + signsp*conc(i1)*charge_lly(i1, ispin)

        End Do
      End Do
      Write (1337, Fmt='(A45,2E17.9)') &
        'RENORM_LLY: Locally summed charge and moment:', &
        (sum0(ispin), ispin=1, nspin)
      Write (1337, Fmt='(A45,2E17.9)') &
        'RENORM_LLY: Renormalized charge and moment:  ', &
        (sum1(ispin), ispin=1, nspin)
      Write (1337, Fmt='(A50,2E17.9)') &
        'RENORM_LLY: Renormalization factor of total charge:', sum1(1)/sum0(1)
! LLY Lloyd  &
! Renormalize the valence charge according to Lloyd's formula.
! Find renormalization constant per energy, then renormalize
    End Subroutine
