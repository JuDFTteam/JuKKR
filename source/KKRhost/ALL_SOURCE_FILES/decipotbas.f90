    Subroutine decipotbas(ihost, iqoff, itoff, nq, nt, rbasis, qmtet, qmphi, &
      noq, kaoez, zat, iqat, conc, irws, ipan, ircut, rr, drdi, visp, nspin, &
      krel, solver, socscl, cscl, vtrel, btrel, irmd, ipand, nembd1, ntmax, &
      nspind, lmaxd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * reads in the potential data for the host atoms from the potential  *
! * file 'decimate.pot'                                                *
! *                                        v.popescu - munich, Dec 04  *
! *                                                                    *
! * Note: so far, only  SPHERICAL case implemented                     *
! *                                                                    *
! **********************************************************************
      Implicit None
!..
!.. Arguments
      Integer :: ihost, ipand, iqoff, irmd, itoff, krel, lmaxd, nembd1, nq, &
        nspin, nspind, nt, ntmax
      Character (Len=10) :: solver
      Real (Kind=dp) :: btrel(irmd*krel+(1-krel), ntmax), conc(ntmax), &
        cscl(krel*lmaxd+1, krel*ntmax+(1-krel)), drdi(irmd, ntmax), &
        qmphi(nembd1), qmtet(nembd1), rbasis(3, nembd1), rr(irmd, ntmax), &
        socscl(krel*lmaxd+1, krel*ntmax+(1-krel)), visp(irmd, ntmax*nspind), &
        vtrel(irmd*krel+(1-krel), ntmax), zat(ntmax)
      Integer :: ipan(ntmax), iqat(nembd1, ntmax), ircut(0:ipand, ntmax), &
        irws(ntmax), kaoez(nembd1, nembd1), noq(nembd1)
!..
!.. Locals
      Integer :: i, ih, ihf, il, ipot1, ipot2
      Integer :: nint
      Real (Kind=dp) :: rmt(ntmax), rws(ntmax)
      Character (Len=3) :: txtt(nt)
! ......................................................................
! --> read basis

      Do ih = 1, nq
        ihf = ih + iqoff
        Read (36+ihost, 100) il, (rbasis(i,ihf), i=1, 3)
        If (ih/=il) Stop ' Inconsistent data '
        Write (1337, 100) il, (rbasis(i,ihf), i=1, 3)
      End Do
      Read (36+ihost, *)
      Write (1337, 120)
      Do ih = 1, nq
        ihf = ih + iqoff
        Read (36+ihost, Fmt=110) il, qmtet(ihf), qmphi(ihf), noq(ihf), &
          (kaoez(i,ihf), i=1, noq(ihf))
        If (ih/=il) Stop ' Inconsistent data '
        Write (1337, 130) ih, qmtet(ihf), qmphi(ihf), noq(ihf), &
          (kaoez(i,ihf), i=1, noq(ihf))
      End Do
      Write (1337, 140)
      If (krel==1) Read (36+ihost, '(7X,A10)') solver
      Read (36+ihost, *)

! --> read atoms

      Do ih = 1, nt
        ihf = ih + itoff
        Read (36+ihost, 150) il, zat(ihf), iqat(1, ihf), conc(ihf), irws(ihf), &
          ipan(ihf), (ircut(i,ihf), i=0, ipan(ihf))
        If (ih/=il) Stop ' Inconsistent data '

        If (krel==1) Then
          Read (36+ihost, 160) socscl(1, ihf), cscl(1, ihf)
          Do il = 2, lmaxd + 1
            socscl(il, ihf) = socscl(1, ihf)
            cscl(il, ihf) = cscl(1, ihf)
          End Do
        End If

      End Do
      Do ih = 1, nt
        ihf = ih + itoff
        ipot1 = (ihf-1)*nspin + 1
        ipot2 = ipot1 + 1
        Read (36+ihost, *)
        Read (36+ihost, 170) il, txtt(ih), rmt(ihf), rws(ihf)
        If (ih/=il) Stop ' Inconsistent data '
        Write (1337, 180) ih, txtt(ih), nint(zat(ihf)), conc(ihf), irws(ihf), &
          rws(ihf)
        If (krel==0) Then
          Read (36+ihost, *)
          Do i = 1, irws(ihf)
            Read (36+ihost, 190) rr(i, ihf), drdi(i, ihf), &
              (visp(i,il), il=ipot1, ipot2)
          End Do
        Else
          Read (36+ihost, '(7X,I3)') il
          irws(ihf) = irws(ihf) - il
          ircut(ipan(ihf), ihf) = irws(ihf)
          Read (36+ihost, *)
          Do i = 1, irws(ihf)
            Read (36+ihost, 190) rr(i, ihf), drdi(i, ihf), vtrel(i, ihf), &
              btrel(i, ihf)
          End Do
        End If
      End Do

100   Format (9X, I3, 3F12.8)
110   Format (7X, I3, 2(7X,F9.4), 7X, I3, 7X, 8I3)
120   Format (9X, 39('-'), /, 9X, '   THETA   ', '   PHI   ', 'OCC', ' IT')
130   Format (9X, I3, 2(F9.4), I3, 8I3)
140   Format (9X, 39('-'), /, 10X, 'ATOMS', /, 15X, 'Z   CONC  IWS    RWS')
150   Format (7X, I3, 7X, F4.0, 7X, I3, 7X, F7.4, /, 17X, I4, 7X, I3, 7X, 6I4)
160   Format (17X, F10.6, 7X, D13.6)
170   Format (7X, I3, 1X, A3, 2(/,7X,F12.8))
180   Format (9X, I3, 1X, A3, I3, F7.4, I4, F10.6)
190   Format (1P, 4D20.12)
    End Subroutine
