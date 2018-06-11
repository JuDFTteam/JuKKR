! 13.10.95 ***************************************************************
    Subroutine mixstr(rmsavq, rmsavm, ins, lpot, lmpot, natref, nshell, &
      nstart, nend, conc, nspin, itc, rfpi, fpi, ipf, mixing, fcm, irc, irmin, &
      r, drdi, vons, visp, vins, vspsmo, vspsme, lsmear)
      Use mod_datatypes, Only: dp
! ************************************************************************
implicit none
!.. Parameters ..
      Include 'inc.p'
      Integer :: lmpotd, irmind
      Parameter (lmpotd=(lpotd+1)**2, irmind=irmd-irnsd)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: fcm, fpi, mixing, rfpi, rmsavm, rmsavq
      Integer :: ins, ipf, itc, lmpot, lpot, natref, nend, nspin, nstart
      Integer :: lsmear
!..
!.. Intrinsic Functions ..
      Real (Kind=dp) :: drdi(irmd, natypd), r(irmd, natypd), &
        vins(irmind:irmd, lmpotd, *), visp(irmd, *), vons(irmd, lmpotd, *), &
        conc(natypd), vspsmo(irmd, nspotd), vspsme(irmd, nspotd)
      Integer :: irc(natypd), irmin(natypd), nshell(0:nsheld)
!     ..

      Real (Kind=dp) :: fac, rmserm, rmserq, vmn, vnm, vnp, voldm, voldp, vpn
      Real (Kind=dp) :: natom
      Integer :: i, ih, ihp1, irc1, irmin1, j, lm, np
!---> final construction of the potentials
!     attention : the spherical averaged potential is the lm=1
      Intrinsic :: mod, real, sqrt
!                     component of vons times sqrt(4 pi).
      rmsavq = 0.0E0_dp
      rmsavm = 0.0E0_dp

!     first mixing scheme : straight mixing
!---> determination of the root mean sqare error





      natom = 0.0E0_dp

      Do np = nstart, nend

        i = np - natref
        natom = natom + real(nshell(i), kind=dp)*conc(i)
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
        If (nspin==2) Then
          ih = 2*np - 1
          ihp1 = ih + 1
        Else
          ih = np
          ihp1 = ih
        End If

        irc1 = irc(np)
        rmserq = 0.0E0_dp
        rmserm = 0.0E0_dp
        fac = 0.5E0_dp/rfpi
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
        Do j = 1, irc1
          vnp = fac*(vons(j,1,ih)+vons(j,1,ihp1))
          vnm = fac*(vons(j,1,ih)-vons(j,1,ihp1))
          voldp = 0.5E0_dp*(visp(j,ih)+visp(j,ihp1))
          voldm = 0.5E0_dp*(visp(j,ih)-visp(j,ihp1))
          rmserq = rmserq + 2.0E0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, &
            np)*drdi(j, np)*(vnp-voldp)*(vnp-voldp)
          rmserm = rmserm + 2.0E0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, &
            np)*drdi(j, np)*(vnm-voldm)*(vnm-voldm)
          vpn = voldp + mixing*(vnp-voldp)
          vmn = voldm + fcm*mixing*(vnm-voldm)
          vons(j, 1, ihp1) = vpn - vmn
          vons(j, 1, ih) = vpn + vmn
        End Do


        If (lsmear>=3) Then
          Do j = 1, irc1
            vnp = 0.5E0_dp*(vspsmo(j,ih)+vspsmo(j,ihp1))
            vnm = 0.5E0_dp*(vspsmo(j,ih)-vspsmo(j,ihp1))
            voldp = 0.5E0_dp*(vspsme(j,ih)+vspsme(j,ihp1))
            voldm = 0.5E0_dp*(vspsme(j,ih)-vspsme(j,ihp1))
            vpn = voldp + mixing*(vnp-voldp)
            vmn = voldm + fcm*mixing*(vnm-voldm)
            vspsmo(j, ihp1) = vpn - vmn
            vspsmo(j, ih) = vpn + vmn
          End Do
        End If

        If ((lsmear==1) .Or. (lsmear==2)) Then
          Do j = 1, irc1
            vspsme(j, ihp1) = vspsmo(j, ihp1)
            vspsme(j, ih) = vspsmo(j, ih)
          End Do
        End If


        rmserq = rmserq/(r(irc1,np)**3)
        rmserm = rmserm/(r(irc1,np)**3)
        rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
        rmsavm = rmsavm + rmserm*nshell(i)*conc(i)

        If (nspin==2) Then
          Write (ipf, Fmt=100) i, sqrt(rmserq), sqrt(rmserm)
        Else
          Write (ipf, Fmt=120) i, sqrt(rmserq)
        End If

        If (ins/=0 .And. lpot>0) Then

          rmserq = 0.0E0_dp
          rmserm = 0.0E0_dp
          irmin1 = irmin(np)
          Do lm = 2, lmpot
            Do j = irmin1, irc1
              vnp = 0.5E0_dp*(vons(j,lm,ih)+vons(j,lm,ihp1))
              vnm = 0.5E0_dp*(vons(j,lm,ih)-vons(j,lm,ihp1))
              voldp = 0.5E0_dp*(vins(j,lm,ih)+vins(j,lm,ihp1))
              voldm = 0.5E0_dp*(vins(j,lm,ih)-vins(j,lm,ihp1))
              rmserq = rmserq + 2.0E0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r( &
                j, np)*drdi(j, np)*(vnp-voldp)*(vnp-voldp)
              rmserm = rmserm + 2.0E0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r( &
                j, np)*drdi(j, np)*(vnm-voldm)*(vnm-voldm)
              vpn = voldp + mixing*(vnp-voldp)
              vmn = voldm + fcm*mixing*(vnm-voldm)
              vons(j, lm, ihp1) = vpn - vmn
              vons(j, lm, ih) = vpn + vmn
            End Do
          End Do
          rmserq = rmserq/(r(irc1,np)**3)/fpi
          rmserm = rmserm/(r(irc1,np)**3)/fpi
          rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
          rmsavm = rmsavm + rmserm*nshell(i)*conc(i)

          If (nspin==2) Then
            Write (ipf, Fmt=110) i, sqrt(rmserq), sqrt(rmserm)
          Else
            Write (ipf, Fmt=130) i, sqrt(rmserq)
          End If

        End If

      End Do
! 13.10.95 ***************************************************************
! ************************************************************************
      rmsavq = sqrt(rmsavq/natom)
      rmsavm = sqrt(rmsavm/natom)
!.. Parameters ..
      Write (1337, '(79("-"),/)')
      If (nspin==2) Then
        Write (ipf, Fmt=140) itc, rmsavq, rmsavm
        Write (6, Fmt=140) itc, rmsavq, rmsavm
      Else
        Write (ipf, Fmt=150) itc, rmsavq
        Write (6, Fmt=150) itc, rmsavq
      End If
      Write (1337, '(79("-"))')
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
100   Format (5X, ' rms-error for atom', I3, 1X, ':', 'v+ + v- = ', 1P, D11.4, &
        2X, ',', 2X, 'v+ - v- = ', 1P, D11.4)
110   Format (5X, ' rms-error non spherical contribution for atom ', I3, 1X, &
        ':', 'v+ + v- = ', 1P, D11.4, 02X, ',', 2X, 'v+ - v- = ', 1P, D11.4)
120   Format (5X, ' rms-error for atom', I3, 1X, ':', 'v+ + v- = ', 1P, D11.4)
130   Format (5X, ' rms-error non spherical contribution for atom ', I3, 1X, &
        ':', 'v+ + v- = ', 1P, D11.4)
140   Format ('      ITERATION', I4, ' average rms-error : v+ + v- = ', 1P, &
        D11.4, /, 39X, ' v+ - v- = ', 1P, D11.4)
150   Format ('      ITERATION', I4, ' average rms-error : v+ + v- = ', 1P, &
        D11.4)
    End Subroutine
