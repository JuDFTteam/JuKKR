!-------------------------------------------------------------------------------
! SUBROUTINE: RHOOUTNEW
!> @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine rhooutnew(nsra, lmax, gmatll, ek, lmpot, df, &
      npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, &
      lmsp, rll, rllleft, sllleft, cden, cdenlm, cdenns, rho2nsc, corbital, &
      gflle_part, rpan_intervall, ipan_intervall)

      Use constants
      Use profiling
      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

      Integer, Intent (In) :: nsra
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: iend !< Number of nonzero gaunt coefficients
      Integer, Intent (In) :: imt1
      Integer, Intent (In) :: ncheb !< Number of Chebychev pannels for the new solver
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: irmdnew
      Integer, Intent (In) :: corbital
      Integer, Intent (In) :: npan_tot
      Complex (Kind=dp), Intent (In) :: ek
      Complex (Kind=dp), Intent (In) :: df
      Integer, Dimension (*), Intent (In) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      Integer, Dimension (*), Intent (In) :: ifunm
      Integer, Dimension (0:ntotd), Intent (In) :: ipan_intervall
      Integer, Dimension (ncleb, 4), Intent (In) :: icleb !< Pointer array
      Real (Kind=dp), Dimension (*), Intent (In) :: cleb !< GAUNT coefficients (GAUNT)
      Real (Kind=dp), Dimension (0:ntotd), Intent (In) :: rpan_intervall
      Real (Kind=dp), Dimension (ntotd*(ncheb+1), nfund), &
        Intent (In) :: thetasnew
      Complex (Kind=dp), Dimension (lmmaxso, lmmaxso), Intent (In) :: gmatll !< GMATLL = diagonal elements of the G matrix (system)
! Note that SLL is not needed for calculation of density, only needed for calculation of Green function
      Complex (Kind=dp), Dimension (nsra*lmmaxso, lmmaxso, irmdnew), &
        Intent (In) :: rll
      Complex (Kind=dp), Dimension (nsra*lmmaxso, lmmaxso, irmdnew), &
        Intent (In) :: rllleft
      Complex (Kind=dp), Dimension (nsra*lmmaxso, lmmaxso, irmdnew), &
        Intent (In) :: sllleft

! .. Output variables
      Complex (Kind=dp), Dimension (irmdnew, 4), Intent (Out) :: cdenns
      Complex (Kind=dp), Dimension (lmmaxso, lmmaxso), &
        Intent (Out) :: gflle_part ! lmlm-dos
      Complex (Kind=dp), Dimension (irmdnew, 0:lmax, 4), Intent (Out) :: cden
      Complex (Kind=dp), Dimension (irmdnew, lmmaxd, 4), &
        Intent (Out) :: cdenlm

! .. In/Out variables
      Complex (Kind=dp), Dimension (irmdnew, lmpot, 4), &
        Intent (Inout) :: rho2nsc

! .. Local variables
      Integer :: ir, jspin, lm1, lm2, lm3, m1, l1, j, ifun
      Integer :: i_stat, i_all
      Real (Kind=dp) :: c0ll
      Complex (Kind=dp) :: cltdf
      Integer, Dimension (4) :: lmshift1
      Integer, Dimension (4) :: lmshift2
      Complex (Kind=dp), Dimension (lmmaxso, lmmaxso, 3) :: loperator

! .. Local allocatable arrays
      Complex (Kind=dp), Dimension (:), Allocatable :: cwr ! lmlm-dos
      Complex (Kind=dp), Dimension (:, :), Allocatable :: qnsi
      Complex (Kind=dp), Dimension (:, :), Allocatable :: pnsi
      Complex (Kind=dp), Dimension (:, :, :), Allocatable :: wr
      Complex (Kind=dp), Dimension (:, :, :), Allocatable :: wr1 ! LDAU
! .. External routines
      Logical :: test, opt
      External :: test, opt

      Allocate (wr(lmmaxso,lmmaxso,irmdnew), Stat=i_stat)
      Call memocc(i_stat, product(shape(wr))*kind(wr), 'WR', 'RHOOUTNEW')
      wr = czero
      Allocate (cwr(irmdnew), Stat=i_stat)
      Call memocc(i_stat, product(shape(cwr))*kind(cwr), 'CWR', 'RHOOUTNEW')
      cwr = czero
      Allocate (wr1(lmmaxso,lmmaxso,irmdnew), Stat=i_stat)
      Call memocc(i_stat, product(shape(wr1))*kind(wr1), 'WR1', 'RHOOUTNEW')
      wr1 = czero
      Allocate (qnsi(lmmaxso,lmmaxso), Stat=i_stat)
      Call memocc(i_stat, product(shape(qnsi))*kind(qnsi), 'QNSI', &
        'RHOOUTNEW')
      qnsi = czero
      Allocate (pnsi(lmmaxso,lmmaxso), Stat=i_stat)
      Call memocc(i_stat, product(shape(pnsi))*kind(pnsi), 'PNSI', &
        'RHOOUTNEW')
      pnsi = czero

! set LMSHIFT value which is need to construct CDEN
      lmshift1(1) = 0
      lmshift1(2) = lmmaxd
      lmshift1(3) = 0
      lmshift1(4) = lmmaxd
      lmshift2(1) = 0
      lmshift2(2) = lmmaxd
      lmshift2(3) = lmmaxd
      lmshift2(4) = 0

! for orbital moment
      If (corbital/=0) Then
        Call calc_orbitalmoment(lmax, lmmaxso, loperator)
      End If

      c0ll = 1E0_dp/sqrt(16E0_dp*atan(1E0_dp))
      cden = czero
      cdenlm = czero

      Do ir = 1, irmdnew
        Do lm1 = 1, lmmaxso
          Do lm2 = 1, lmmaxso
            qnsi(lm1, lm2) = sllleft(lm1, lm2, ir)
!          PNSI(LM1,LM2)=RLL(LM1,LM2,IR)
            pnsi(lm1, lm2) = rllleft(lm1, lm2, ir)
          End Do
        End Do
!        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
!     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
        Call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, &
          gmatll, lmmaxso, ek, qnsi, lmmaxso)
        Do lm1 = 1, lmmaxso
          Do lm2 = 1, lmmaxso
            pnsi(lm1, lm2) = rll(lm1, lm2, ir)
          End Do
        End Do
        Call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, &
          qnsi, lmmaxso, czero, wr(1,1,ir), lmmaxso)
!
        If (nsra==2) Then
          Do lm1 = 1, lmmaxso
            Do lm2 = 1, lmmaxso
!          QNSI(LM1,LM2)=SLLLEFT(LM1+LMMAXSO,LM2,IR)
              qnsi(lm1, lm2) = -sllleft(lm1+lmmaxso, lm2, ir)
!          PNSI(LM1,LM2)=RLLLEFT(LM1+LMMAXSO,LM2,IR)
              pnsi(lm1, lm2) = -rllleft(lm1+lmmaxso, lm2, ir)
            End Do
          End Do
!        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
!     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
          Call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, &
            gmatll, lmmaxso, ek, qnsi, lmmaxso)
          Do lm1 = 1, lmmaxso
            Do lm2 = 1, lmmaxso
              pnsi(lm1, lm2) = rll(lm1+lmmaxso, lm2, ir)
            End Do
          End Do
          Call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, &
            qnsi, lmmaxso, cone, wr(1,1,ir), lmmaxso)
        End If
!
! For orbital moment
        If (corbital/=0) Then
          Call zgemm('N', 'N', lmmaxso, lmmaxso, lmmaxso, cone, &
            loperator(1,1,corbital), lmmaxso, wr(1,1,ir), lmmaxso, czero, &
            pnsi, lmmaxso)
          Do lm1 = 1, lmmaxso
            Do lm2 = 1, lmmaxso
              wr(lm1, lm2, ir) = pnsi(lm1, lm2)
            End Do
          End Do
        End If
        Do lm1 = 1, lmmaxso
          Do lm2 = 1, lmmaxso
            wr1(lm1, lm2, ir) = wr(lm1, lm2, ir)
          End Do
        End Do
        Do lm1 = 1, lmmaxso
          Do lm2 = 1, lm1 - 1
            wr1(lm1, lm2, ir) = wr1(lm1, lm2, ir) + wr1(lm2, lm1, ir)
          End Do
        End Do
!
        Do jspin = 1, 4
          Do lm1 = 1, lmmaxd
            Do lm2 = 1, lm1 - 1
              wr(lm1+lmshift1(jspin), lm2+lmshift2(jspin), ir) &
                = wr(lm1+lmshift1(jspin), lm2+lmshift2(jspin), ir) + &
                wr(lm2+lmshift1(jspin), lm1+lmshift2(jspin), ir)
            End Do
          End Do
        End Do ! JSPIN
      End Do !IR

! IF lmdos or LDAU
      If (opt('lmlm-dos') .Or. opt('LDA+U   ')) Then ! lmlm-dos
! Integrate only up to muffin-tin radius.                                  ! lmlm-dos
        gflle_part = czero ! lmlm-dos
        Do lm2 = 1, lmmaxso ! lmlm-dos
          Do lm1 = 1, lmmaxso ! lmlm-dos
! For integration up to MT radius do this:                           ! lmlm-dos
! CWR(1:IMT1) = WR(LM1,LM2,1:IMT1)                                   ! lmlm-dos
! CWR(IMT1+1:IRMDNEW) = CZERO                                        ! lmlm-dos
! CALL INTCHEB_CELL(CWR,GFLLE_PART(LM1,LM2),RPAN_INTERVALL,&         ! lmlm-dos
!     IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)                         ! lmlm-dos
! For full cell integration replace loop content with this:          ! lmlm-dos
            cwr(1:irmdnew) = wr1(lm1, lm2, 1:irmdnew) ! lmlm-dos
! If LDAU, integrate only up to MT
            Do ir = imt1 + 1, irmdnew
              If (opt('LDA+U   ')) Then
                cwr(ir) = czero ! LDAU
              Else
                cwr(ir) = cwr(ir)*thetasnew(ir, 1)*c0ll ! lmlm-dos
              End If
            End Do
            Call intcheb_cell(cwr, gflle_part(lm1,lm2), rpan_intervall, &
              ipan_intervall, npan_tot, ncheb, irmdnew)
          End Do
        End Do
      End If ! OPT('lmlm-dos').OR.OPT('LDA+U   ')
!
!      DO IR = 1,IRMDNEW
!       DO JSPIN = 1,4
!        DO LM1 = 1,LMMAXD
!         DO LM2 = 1,LM1-1
!          WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)=
!    +           WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)+
!    +           WR(LM2+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
!         ENDDO
!        ENDDO
!       ENDDO ! JSPIN
!      ENDDO !IR
!
! First calculate the spherical symmetric contribution
!
      Do l1 = 0, lmax
        Do m1 = -l1, l1
          lm1 = l1*(l1+1) + m1 + 1
          Do ir = 1, irmdnew
            Do jspin = 1, 4
              cden(ir, l1, jspin) = cden(ir, l1, jspin) + &
                wr(lm1+lmshift1(jspin), lm1+lmshift2(jspin), ir)
              cdenlm(ir, lm1, jspin) = wr(lm1+lmshift1(jspin), &
                lm1+lmshift2(jspin), ir)
            End Do ! JPSIN
          End Do ! IR
        End Do ! M1
!
        Do jspin = 1, 4
          Do ir = 1, irmdnew
            rho2nsc(ir, 1, jspin) = rho2nsc(ir, 1, jspin) + &
              c0ll*(cden(ir,l1,jspin)*df)
          End Do ! IR
!
          Do ir = imt1 + 1, irmdnew
            cden(ir, l1, jspin) = cden(ir, l1, jspin)*thetasnew(ir, 1)*c0ll
            Do m1 = -l1, l1
              lm1 = l1*(l1+1) + m1 + 1
              cdenlm(ir, lm1, jspin) = cdenlm(ir, lm1, jspin)* &
                thetasnew(ir, 1)*c0ll
            End Do ! M1
          End Do ! IR
        End Do ! JSPIN
      End Do ! L1
!
      cdenns = czero
!
      Do j = 1, iend
        lm1 = icleb(j, 1)
        lm2 = icleb(j, 2)
        lm3 = icleb(j, 3)
        cltdf = df*cleb(j)
        Do jspin = 1, 4
          Do ir = 1, irmdnew
            rho2nsc(ir, lm3, jspin) = rho2nsc(ir, lm3, jspin) + &
              (cltdf*wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir))
          End Do
!
          If (lmsp(lm3)>0) Then
            ifun = ifunm(lm3)
            Do ir = imt1 + 1, irmdnew
              cdenns(ir, jspin) = cdenns(ir, jspin) + &
                cleb(j)*wr(lm1+lmshift1(jspin), lm2+lmshift2(jspin), ir)* &
                thetasnew(ir, ifun)
            End Do
          End If
        End Do ! JSPIN
      End Do ! J

      i_all = -product(shape(wr))*kind(wr)
      Deallocate (wr, Stat=i_stat)
      Call memocc(i_stat, i_all, 'WR', 'RHOOUTNEW')
      i_all = -product(shape(wr1))*kind(wr1)
      Deallocate (wr1, Stat=i_stat)
      Call memocc(i_stat, i_all, 'WR1', 'RHOOUTNEW')
      i_all = -product(shape(cwr))*kind(cwr)
      Deallocate (cwr, Stat=i_stat)
      Call memocc(i_stat, i_all, 'CWR', 'RHOOUTNEW')
      i_all = -product(shape(qnsi))*kind(qnsi)
      Deallocate (qnsi, Stat=i_stat)
      Call memocc(i_stat, i_all, 'QNSI', 'RHOOUTNEW')
      i_all = -product(shape(pnsi))*kind(pnsi)
      Deallocate (pnsi, Stat=i_stat)
      Call memocc(i_stat, i_all, 'PNSI', 'RHOOUTNEW')

    End Subroutine
