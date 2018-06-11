    Subroutine decitset(alat, bravsys, ez, ielast, nlbasis, nrbasis, fileleft, &
      fileright, ins, kvrel, krel, nspin, kmrot, vref, rmtref, nref, refpot, &
      lefttinv, righttinv, vacflag, nembd1, iemxd, irmd, ipand, lmaxd, lmgf0d, &
      lmmaxd, lm2d, nspind)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * This subroutine is thought as an alternative to the < decimaread > *
! * which requires an a priori calculated set of single-site matrices  *
! * over a fixed energy mesh. It is using the potential as written out *
! * in < outpothost > routine and determines the matrix                *
! *                                                                    *
! *            /         \-1   /             \-1                       *
! *            | Delta t |   = |  t    - t   |                         *
! *            \         /     \  sys    ref /                         *
! *                                                                    *
! * for the left and the right host, using the energy mesh as read in  *
! * from the input-file                                                *
! *                                                                    *
! *                                        v.popescu - munich, Dec 04  *
! *                                                                    *
! * Notes: - no charge moments are calculated -- thus this option CAN  *
! *          NOT be used in SCF calculations                           *
! *        - non-spherical case not implemented, neither LDA+U (al-    *
! *          though the interface to regsol in decitmat is supplied)   *
! *        - CPA case not implemented - requires BZ integration        *
! *                                                                    *
! **********************************************************************
      Implicit None
!..
!.. Scalars arguments ..
      Integer :: iemxd, nembd1, lmmaxd, ipand, nspind, irmd, lmaxd, lm2d, &
        lmgf0d
      Integer :: ielast, kmrot, nlbasis, nrbasis, nref, ins, kvrel, krel, &
        nspin
      Real (Kind=dp) :: alat
      Character (Len=40) :: fileleft, fileright
!..
!.. Array arguments ..
      Integer :: refpot(nembd1)
      Real (Kind=dp) :: vref(*), rmtref(*), bravsys(3, 3)
      Complex (Kind=dp) :: ez(iemxd)
      Complex (Kind=dp) :: lefttinv(lmmaxd, lmmaxd, nembd1, nspind, iemxd), &
        righttinv(lmmaxd, lmmaxd, nembd1, nspind, iemxd)
      Logical :: vacflag(2)
!..
!.. Local scalars ..
      Integer :: ihost, i, ll, mm, lngstring, nqhost, ilhost
      Integer :: nhost
      Integer :: nq, nt, iqoff, itoff, ie, ih, iqh, ioq, info
      Integer :: ipot, i1, ispin, nsra, lm1, lm2, irc1, iref
      Integer :: ntleft, ntright, nthost
!.. LDA+U
      Integer :: idoldau, lopt
      Real (Kind=dp) :: wldauav
!..
      Real (Kind=dp) :: efermi, rirc
      Complex (Kind=dp) :: eryd, carg, cfctor
      Logical :: test
      Character (Len=40) :: filehost
      Character (Len=10) :: solver
!..
!.. Local arrays
      Integer :: krelh(2), nspinh(2), insh(2), ipvt(lmmaxd)
      Integer :: noq(nembd1), kaoez(nembd1, nembd1), inhost(2)
      Real (Kind=dp) :: bravais(3, 3, 2), rbasis(3, nembd1), qmtet(nembd1), &
        qmphi(nembd1)
      Character (Len=5) :: chhost(2)
      Character (Len=9) :: txts(2)
!..
!.. Allocatable local arrays 
      Integer :: ntmax
      Real (Kind=dp) :: zat(:), rws(:), rmt(:), conc(:)
      Real (Kind=dp) :: rr(:, :), drdi(:, :), visp(:, :), dror(:, :)
      Real (Kind=dp) :: socscl(:, :), cscl(:, :)
      Integer :: irws(:), ipan(:), iqat(:, :), ircut(:, :), loflm(:)
      Complex (Kind=dp) :: trefll(:, :, :), tmatll(:, :), dhmat(:, :, :)
      Complex (Kind=dp) :: dtrefll(:, :, :) ! LLY Lloyd
      Complex (Kind=dp) :: alpharef(:, :), dalpharef(:, :) ! LLY Lloyd Alpha matrix and deriv.
      Real (Kind=dp) :: vtrel(:, :), btrel(:, :), r2drdirel(:, :)
      Integer :: zrel(:)
      Allocatable :: zat, rws, rmt, conc, rr, drdi, visp, dror, socscl, cscl
      Allocatable :: vtrel, btrel, r2drdirel
      Allocatable :: irws, ipan, iqat, ircut, loflm, zrel
      Allocatable :: trefll, tmatll, dhmat
      Allocatable :: dtrefll, alpharef, dalpharef ! LLY
!.. 
!.. External subroutines
      External :: calctref13, changerep, cinit, cmatstr, decipotbas, &
        decipothead, decitmat, zaxpy, zcopy, zgetrf, zgetri
!..
!.. External Functions ..
      External :: lngstring, test
!..
!.. Data statements
      Data chhost/'LEFT ', 'RIGHT'/
      Data txts/'spin   UP', 'spin DOWN'/
! ......................................................................

      cfctor = alat/(8.E0_dp*atan(1.0E0_dp)) ! = ALAT/(2*PI)

      idoldau = 0
      lopt = -1
      wldauav = 0E0_dp
      Allocate (loflm(lm2d), Stat=i1)
      If (i1/=0) Stop '    Allocate LOFLM'
      Write (6, '(5X,A,/,8X,65("-"))') 'Reading in host potentials'
      vacflag(1) = .False.
      vacflag(2) = .False.
      nsra = 1
      If (kvrel>=1) nsra = 2
      i = 1
      Do ll = 0, 2*lmaxd
        Do mm = -ll, ll
          loflm(i) = ll
          i = i + 1
        End Do
      End Do
      ntleft = 0
      ntright = 0
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      nhost = 0
      Do ihost = 1, 2
        filehost = fileleft
        nqhost = nlbasis
        If (ihost==2) Then
          filehost = fileright
          nqhost = nrbasis
        End If
        ilhost = lngstring(filehost, 40)
        Call decipothead(ihost, filehost, ilhost, nqhost, vacflag, alat, &
          bravsys, nq, nt, bravais(1,1,ihost), efermi, insh(ihost), &
          krelh(ihost), nspinh(ihost), ins, krel, nspin, kmrot)

        If (.Not. vacflag(ihost)) Then
          nhost = nhost + 1
          inhost(nhost) = ihost
          If (ihost==1) Then
            ntleft = nt
          Else
            ntright = nt
          End If
        End If


      End Do

      If (ntleft+ntright<=0) Then
        Write (6, &
          '(8X,"Vacuum will be considered on both sides",/, 8X,65("-"))')
        Return
      End If

      ntmax = ntleft + ntright
      Allocate (zat(ntmax), rws(ntmax), rmt(ntmax), conc(ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate ZAT/RWS/RMT/CONC'
      Allocate (rr(irmd,ntmax), drdi(irmd,ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate RR/DRDI'
      Allocate (visp(irmd,ntmax*nspind), Stat=i1)
      If (i1/=0) Stop '    Allocate VISP'
      Allocate (irws(ntmax), ipan(ntmax), iqat(nembd1,ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate IRWS/IPAN/IQAT'
      Allocate (ircut(0:ipand,ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate IRCUT'
      Allocate (socscl(krel*lmaxd+1,krel*ntmax+(1-krel)), Stat=i1)
      If (i1/=0) Stop '    Allocate SOCSCL'
      Allocate (cscl(krel*lmaxd+1,krel*ntmax+(1-krel)), Stat=i1)
      If (i1/=0) Stop '    Allocate CSCL'
      Allocate (vtrel(irmd*krel+(1-krel),ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate VTREL'
      Allocate (btrel(irmd*krel+(1-krel),ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate BTREL'

! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      Do ihost = 1, 2
        Write (6, '(8X,A5," side host: ")', Advance='no') chhost(ihost)
        iqoff = 0
        itoff = 0
        nqhost = nlbasis
        nthost = ntleft
        filehost = fileleft
        If (ihost==2) Then
          nqhost = nrbasis
          nthost = ntright
          iqoff = nlbasis
          itoff = ntleft
          filehost = fileright
        End If
        ilhost = lngstring(filehost, 40)

        If (filehost(1:7)=='vacuum') Then
          Write (6, '(A,/,8X,65("-"))') 'VACUUM will be used'
        Else
          Write (6, '(A,/)') filehost(1:ilhost)
          Write (6, 130) krelh(ihost), nspinh(ihost), insh(ihost), kmrot, &
            nqhost, alat, efermi
          Write (6, 140)((bravais(ll,mm,ihost),mm=1,3), ll=1, 3)
          Call decipotbas(ihost, iqoff, itoff, nqhost, nthost, rbasis, qmtet, &
            qmphi, noq, kaoez, zat, iqat, conc, irws, ipan, ircut, rr, drdi, &
            visp, nspinh(ihost), krelh(ihost), solver, socscl, cscl, vtrel, &
            btrel, irmd, ipand, nembd1, ntmax, nspind, lmaxd)
          Write (6, '(8X,65("-"))')
        End If
      End Do
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      Allocate (dror(irmd,ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate DROR'
      Allocate (r2drdirel(irmd*krel+(1-krel),ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate R2DRDIREL'
      Allocate (zrel(ntmax), Stat=i1)
      If (i1/=0) Stop '    Allocate ZREL'

      Write (6, '(/,5X,A,/)') 'Calculating host (Delta_t)^(-1) matrices'
      If (krel==0) Then
        Do i = 1, ntleft + ntright
          irc1 = ircut(ipan(i), i)
          Do i1 = 2, irc1
            dror(i1, i) = drdi(i1, i)/rr(i1, i)
          End Do
        End Do
      Else
        Do i = 1, ntleft + ntright
          irc1 = ircut(ipan(i), i)
          Do i1 = 1, irc1
            r2drdirel(i1, i) = rr(i1, i)*rr(i1, i)*drdi(i1, i)
          End Do
          zrel(i) = nint(zat(i))
        End Do
      End If

! ******************************************************* energy loop IE
      Allocate (trefll(lmmaxd,lmmaxd,nref), Stat=i1)
      If (i1/=0) Stop '    Allocate TREFLL'
      Allocate (dtrefll(lmmaxd,lmmaxd,nref), Stat=i1) ! LLY
      If (i1/=0) Stop '    Allocate DTREFLL' ! LLY
      Allocate (tmatll(lmmaxd,lmmaxd), dhmat(lmmaxd,lmmaxd,2), Stat=i1)
      If (i1/=0) Stop '    Allocate TMATLL/DHMAT'
      Allocate (alpharef(0:lmaxd,nref), dalpharef(0:lmaxd,nref), Stat=i1) ! LLY Lloyd Alpha matrix AND deriv.
      If (i1/=0) Stop '    Allocate ALPHAREF/DALPHAREF'

      Do ie = 1, ielast
        eryd = ez(ie)

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! -> set up t matrices for the reference system

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        If (krel==0) Then
          Do i1 = 1, nref
            Call calctref13(eryd, vref(i1), rmtref(i1), lmaxd, ih, &
              trefll(1,1,i1), dtrefll(1,1,i1), alpharef(0,i1), &
              dalpharef(0,i1), lmaxd+1, lmgf0d)
          End Do
        End If
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
        Do ilhost = 1, nhost
          ihost = inhost(ilhost)
          iqoff = 0
          itoff = 0
          nqhost = nlbasis
          If (ihost==2) Then
            nqhost = nrbasis
            iqoff = nlbasis
            itoff = ntleft
          End If
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ sites in host
          Do ih = 1, nqhost

! -> assign Delta_t = -t_ref

!  Note: REFPOT(1) = REFPOT(NAEZ) (i.e., of the 2D system)

            iqh = iqoff + ih
            iref = refpot(iqh+1)
            Do lm2 = 1, lmmaxd
              Do lm1 = 1, lmmaxd
                dhmat(lm1, lm2, 1) = -trefll(lm1, lm2, iref)
              End Do
            End Do

            If (nspinh(ihost)>1) Then
              Do lm2 = 1, lmmaxd
                Call zcopy(lmmaxd, dhmat(1,lm2,1), 1, dhmat(1,lm2,2), 1)
              End Do
            End If
! ====================================================== spins and atoms
            Do ispin = 1, nspinh(ihost)
! ----------------------------------------------------------------------
              Do ioq = 1, noq(iqh)
                i1 = kaoez(ioq, iqh) + itoff
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> calculate t_sys for the atom I1 located on site IH

                ipot = (i1-1)*nspinh(ihost) + ispin
                irc1 = ircut(ipan(i1), i1)
                rirc = rr(irc1, i1)

                Call decitmat(eryd, zat(i1), ipan(i1), rr(1,i1), dror(1,i1), &
                  visp(1,ipot), ircut(0,i1), rirc, krel, nsra, ins, tmatll, &
                  loflm, idoldau, lopt, wldauav, solver, socscl(1,krel*i1+(1- &
                  krel)), cscl(1,krel*i1+(1-krel)), zrel(i1), vtrel(1,i1), &
                  btrel(1,i1), drdi(1,i1), r2drdirel(1,i1), ipand, irmd, &
                  lmaxd, lmaxd+1, lm2d, lmmaxd)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tmat calculated

! -> Delta_t = Delta_t + CONC(I1)*t_mat(I1)

                carg = conc(i1)
                Do lm2 = 1, lmmaxd
                  Call zaxpy(lmmaxd, carg, tmatll(1,lm2), 1, &
                    dhmat(1,lm2,ispin), 1)
                End Do
              End Do
! ----------------------------------------------------------------------
! tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
              If (test('tmat    ')) Then
                Write (1337, *)
                Write (1337, 100, Advance='no') &
                  '      ---> Delta_t  matrix for site: ', iqh
                If (krel==0) Write (1337, 110, Advance='no') txts(ispin)
                Write (1337, 120) ', energy: ', eryd
                Call cmatstr(' ', 1, dhmat(1,1,ispin), lmmaxd, lmmaxd, &
                  2*krel+1, 2*krel+1, 0, 1E-8_dp, 6)
                Write (1337, *)
              End If
! tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt

! --> inversion

              Call zgetrf(lmmaxd, lmmaxd, dhmat(1,1,ispin), lmmaxd, ipvt, &
                info)
              Call zgetri(lmmaxd, dhmat(1,1,ispin), lmmaxd, ipvt, tmatll, &
                lmmaxd*lmmaxd, info)
            End Do
! ======================================================================

! --> scaling the host t-matrices to p.u.

            If (ihost==1) Then
              Do ispin = 1, nspinh(ihost)
                Do lm2 = 1, lmmaxd
                  Do lm1 = 1, lmmaxd
                    lefttinv(lm1, lm2, ih, ispin, ie) = cfctor* &
                      dhmat(lm1, lm2, ispin)
                  End Do
                End Do
              End Do
            Else
              Do ispin = 1, nspinh(ihost)
                Do lm2 = 1, lmmaxd
                  Do lm1 = 1, lmmaxd
                    righttinv(lm1, lm2, ih, ispin, ie) = cfctor* &
                      dhmat(lm1, lm2, ispin)
                  End Do
                End Do
              End Do
            End If
          End Do
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        End Do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      End Do
! **********************************************************************
      Deallocate (zat, rws, rmt, conc, rr, drdi, visp, Stat=i1)
      If (i1/=0) Stop '   Deallocate ZAT/RWS/RMT/.../VISP'
      Deallocate (irws, ipan, iqat, ircut, loflm, Stat=i1)
      If (i1/=0) Stop '   Deallocate IRWS/IPAN/IQAT/IRCUT/LOFLM'
      Deallocate (trefll, tmatll, dhmat, Stat=i1)
      If (i1/=0) Stop '   Deallocate TREFLL/TMATLL/DHMAT'
      Deallocate (socscl, cscl, vtrel, btrel, Stat=i1)
      If (i1/=0) Stop '   Deallocate SOCSCL/CSCL/VTREL/BTREL'
      Deallocate (alpharef, dalpharef, Stat=i1)
      If (i1/=0) Stop '   Deallocate ALPHAREF/DALPHAREF'
      If (krel==0) Then
        Deallocate (dror, Stat=i1)
        If (i1/=0) Stop '   Deallocate DROR'
      Else
        Deallocate (r2drdirel, zrel, Stat=i1)
        If (i1/=0) Stop '   Deallocate R2DRDIREL/ZREL'
      End If
100   Format (A, I3)
110   Format (', ', A)
120   Format (A, 2F10.6)
130   Format (10X, 'KREL= ', I1, ' NSPIN= ', I1, ' INS= ', I1, ' KMROT= ', I1, &
        /, 10X, 'NAEZ=', I3, ' ALAT= ', F9.6, ' EFERMI= ', F9.6)
140   Format (10X, 'BRAVAIS ', /, 10X, 3F8.4, /, 10X, 3F8.4, /, 10X, 3F8.4, /, &
        10X, 'RBASIS')
    End Subroutine
