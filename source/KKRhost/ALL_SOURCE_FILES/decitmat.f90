    Subroutine decitmat(eryd, zat, ipan, rr, dror, visp, ircut, rirc, krel, &
      nsra, ins, tmatll, loflm, idoldau, lopt, wldauav, solver, soctl, ctl, &
      zrel, vtrel, btrel, drdi, r2drdi, ipand, irmd, lmaxd, lmaxdp1, lm2d, &
      lmmaxd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * A modified form of the CALCTMAT routine to deal with the host      *
! * t-matrices in case of decimation                                   *
! *                                                                    *
! * Non-spherical potential not implemented yet, neither LDA+U         *
! *                                                                    *
! **********************************************************************
      Implicit None

!Parameters ..
      Real (Kind=dp) :: cvlight
      Parameter (cvlight=274.0720442E0_dp)
      Complex (Kind=dp) :: ci
      Parameter (ci=(0E0_dp,1E0_dp))

!Scalar arguments ..
      Integer :: idoldau, ipan, krel, lopt, nsra, ins, zrel
      Integer :: ipand, irmd, lm2d, lmaxd, lmaxdp1, lmmaxd
      Real (Kind=dp) :: zat, rirc, wldauav
      Complex (Kind=dp) :: eryd
      Character (Len=10) :: solver

!Array arguments ..
      Integer :: ircut(0:ipand), loflm(lm2d)
      Real (Kind=dp) :: rr(irmd), dror(irmd), visp(irmd)
      Complex (Kind=dp) :: tmatll(lmmaxd, lmmaxd)
      Real (Kind=dp) :: soctl(krel*lmaxd+1)
      Real (Kind=dp) :: ctl(krel*lmaxd+1)
      Real (Kind=dp) :: vtrel(irmd*krel+(1-krel))
      Real (Kind=dp) :: btrel(irmd*krel+(1-krel))
      Real (Kind=dp) :: drdi(irmd), r2drdi(irmd*krel+(1-krel))

!Local scalars ..
      Integer :: ll, lm1
      Real (Kind=dp) :: rirc1
      Complex (Kind=dp) :: ek, carg, qf, hlw, blw

!Local arrays ..
      Real (Kind=dp) :: cutoff(irmd)
      Real (Kind=dp) :: rs(:, :), s(:)
      Complex (Kind=dp) :: bessjw(:), bessyw(:), hankws(:), dlogdp(:)
      Complex (Kind=dp) :: tmat(:), mass(:), hamf(:, :), fz(:, :), pz(:, :)
      Allocatable :: rs, s
      Allocatable :: bessjw, bessyw, hankws, dlogdp
      Allocatable :: tmat, mass, hamf, fz, pz

!External subroutines ..
      External :: beshan, cinit, regsol, wfmesh


      Call cinit(lmmaxd*lmmaxd, tmatll)
! ================================================================= KREL
      If (krel==0) Then
        Allocate (bessjw(0:lmaxdp1), bessyw(0:lmaxdp1), Stat=lm1)
        If (lm1/=0) Stop '    Allocate BESSJW/BESSYW'
        Allocate (hankws(0:lmaxdp1), dlogdp(0:lmaxd), Stat=lm1)
        If (lm1/=0) Stop '    Allocate HANKWS/DLOGFP'
        Allocate (tmat(0:lmaxd), mass(irmd), Stat=lm1)
        If (lm1/=0) Stop '    Allocate TMAT/MASS'
        Allocate (hamf(irmd,0:lmaxd), fz(irmd,0:lmaxd), Stat=lm1)
        If (lm1/=0) Stop '    Allocate HAMF/FZ'
        Allocate (pz(irmd,0:lmaxd), Stat=lm1)
        If (lm1/=0) Stop '    Allocate PZ'
        Allocate (rs(irmd,0:lmaxd), s(0:lmaxd), Stat=lm1)
        If (lm1/=0) Stop '    Allocate RS/S'
        rirc1 = 1E0_dp/rirc
        Call wfmesh(eryd, ek, cvlight, nsra, zat, rr, s, rs, ircut(ipan), &
          irmd, lmaxd)

        carg = rirc*ek
        Call beshan(hankws, bessjw, bessyw, carg, lmaxdp1)
        Do ll = 0, lmaxdp1
          hankws(ll) = bessyw(ll) - ci*bessjw(ll)
        End Do

        Call regsol(cvlight, eryd, nsra, dlogdp, fz, hamf, mass, pz, dror, rr, &
          s, visp, zat, ipan, ircut, idoldau, lopt, wldauav, cutoff, irmd, &
          ipand, lmaxd)

! ----------------------------------------------------------------------
! --> determine KREL=0 t - matrix

        Do ll = 0, lmaxd
          qf = real(ll, kind=dp)*rirc1
          hlw = hankws(ll)*dlogdp(ll)
          blw = bessjw(ll)*dlogdp(ll)

          hlw = qf*hankws(ll) - ek*hankws(ll+1) - hlw
          blw = blw - qf*bessjw(ll) + ek*bessjw(ll+1)
          hlw = hlw*ek
          tmat(ll) = blw/hlw
        End Do

! --> spherical/non-spherical

        If (ins==0) Then
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          Do lm1 = 1, lmmaxd
            tmatll(lm1, lm1) = tmat(loflm(lm1))
          End Do
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Else
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          Stop ' not implemented'
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        End If
        Deallocate (bessjw, bessyw, hankws, dlogdp, Stat=lm1)
        If (lm1/=0) Stop '    Deallocate'
        Deallocate (tmat, mass, hamf, fz, pz, Stat=lm1)
        If (lm1/=0) Stop '    Deallocate'
        Deallocate (rs, s, Stat=lm1)
        If (lm1/=0) Stop '    Deallocate'
! ----------------------------------------------------------------------
      Else ! KREL
        Call drvreltmat(eryd, tmatll, vtrel, btrel, rr, drdi, r2drdi, zrel, &
          ircut(ipan), solver, soctl, ctl, lmmaxd, lmaxd, irmd)
      End If
! ================================================================= KREL
    End Subroutine
