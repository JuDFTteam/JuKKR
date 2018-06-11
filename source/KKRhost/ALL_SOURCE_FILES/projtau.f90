    Subroutine projtau(icpaflag, cpachng, kmrot, wrtau, wrtaumq, ifiltau, &
      eryd, nt, nq, nkmq, msst, mssq, nlinq, iqat, conc, tauq, taut, tautlin, &
      ikm1lin, ikm2lin, drotq, ntmax, nqmax, nkmmax, linmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   calculate the component projected TAU - matrices               *
!   *                                                                  *
!   *      TAU(IT) =  TAU(IQ) * ( 1 + (m(t)-m(c))*TAU(IQ) )**(-1)      *
!   *                                                                  *
!   *   NOTE: it is assumed that all equivalent sites  IQ  have the    *
!   *   same TAU-matrix  TAUQ(IQ). To get  TAU(IT)  the first site IQ  *
!   *   occupied by type IT is taken to be representative for          *
!   *   all other (NAT(IT)-1) sites occupied by IT                     *
!   *                                                                  *
!   *   allows an atom type IT to have different orientation of        *
!   *   its moment on different but equivalent sites  IQ               *
!   *                                                                  *
!   * 01/11/00                                                         *
!   ********************************************************************
      Implicit None

! PARAMETER definitions
      Real (Kind=dp) :: tol
      Parameter (tol=1.0E-6_dp)
      Complex (Kind=dp) :: c0, c1
      Parameter (c0=(0.0E0_dp,0.0E0_dp), c1=(1.0E0_dp,0.0E0_dp))

! Dummy arguments
      Real (Kind=dp) :: cpachng
      Complex (Kind=dp) :: eryd
      Integer :: icpaflag, ifiltau, kmrot, linmax, nkmmax, nq, nqmax, nt, &
        ntmax
      Logical :: wrtau, wrtaumq
      Real (Kind=dp) :: conc(ntmax)
      Complex (Kind=dp) :: drotq(nkmmax, nkmmax, nqmax), &
        mssq(nkmmax, nkmmax, nqmax), msst(nkmmax, nkmmax, ntmax), &
        tauq(nkmmax, nkmmax, nqmax), taut(nkmmax, nkmmax, ntmax), &
        tautlin(linmax, ntmax)
      Integer :: ikm1lin(linmax), ikm2lin(linmax), iqat(ntmax), nkmq(nqmax), &
        nlinq(nqmax)

! Local variables
      Real (Kind=dp) :: cpac
      Complex (Kind=dp) :: dmamc(nkmmax, nkmmax), dmattg(nkmmax, nkmmax), &
        dtiltg(nkmmax, nkmmax), rmss, rtau, w1(nkmmax, nkmmax)
      Integer :: i, icpaf, iq, it, j, lin, m, n

      Do it = 1, nt

! ---------- pick first site IQ occupied by type IT to be representative
! ----------- all other (NAT(IT)-1) occupied sites have to be equivalent

        iq = iqat(it)
        m = nkmmax
        n = nkmq(iq)

        If (conc(it)<0.995_dp) Then

! ------------------------- rotate the single site m-matrix if necessary
          If (kmrot/=0) Then

            Call rotate(msst(1,1,it), 'L->G', w1, n, drotq(1,1,iq), m)

            Call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), &
              w1, m)

          Else

            Call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), &
              msst(1,1,it), m)

          End If

!     -------------------------------------------
!              TAU(t) = TAU * D~(t)
!     ----------------------------------------
          Call zgemm('N', 'N', n, n, n, c1, tauq(1,1,iq), m, dtiltg, m, c0, &
            taut(1,1,it), m)

          icpaf = icpaflag
          cpac = cpachng

        Else

!     CONC > 0.995:  COPY TAU TO TAUTLIN

          Do j = 1, n
            Call zcopy(n, tauq(1,j,iq), 1, taut(1,j,it), 1)
          End Do

          icpaf = 0
          cpac = 0E0_dp

        End If

!     -------------------------------------------
!            rotate  TAU(t)  if required
!     -------------------------------------------

        If (kmrot/=0) Then

          Do j = 1, n
            Call zcopy(n, taut(1,j,it), 1, w1(1,j), 1)
          End Do

          Call rotate(w1, 'G->L', taut(1,1,it), n, drotq(1,1,iq), m)

        End If

!     -------------------------------------------
!        STORE TAU(t) IN LINEAR ARRAY TAUTLIN
!     -------------------------------------------

        Do lin = 1, nlinq(iq)
          tautlin(lin, it) = taut(ikm1lin(lin), ikm2lin(lin), it)
        End Do

        If (wrtau) Then
          Write (ifiltau, 100) eryd, it, iq, icpaf, cpac
          Do i = 1, n
            Do j = 1, n
              If (i==j) Then
                Write (ifiltau, 120) i, j, taut(i, j, it)
              Else
                If (abs(taut(i,j,it)/taut(i,i,it))>tol) Write (ifiltau, 120) i &
                  , j, taut(i, j, it)
              End If
            End Do
          End Do
        End If

      End Do
!================================================================= IT ==

      If (wrtaumq) Then
        Do iq = 1, nq
          Write (ifiltau, 110) eryd, iq, icpaflag, cpachng
          Do i = 1, n
            Do j = 1, n
              If (i==j) Then
                Write (ifiltau, 120) i, j, tauq(i, j, iq), mssq(i, j, iq)
              Else
                rtau = tauq(i, j, iq)/tauq(i, i, iq)
                rmss = mssq(i, j, iq)/mssq(i, i, iq)
                If ((abs(rtau)>tol) .Or. (abs(rmss)>tol)) Write (ifiltau, 120) &
                  i, j, tauq(i, j, iq), mssq(i, j, iq)
              End If

            End Do
          End Do

        End Do
      End If
!--------------------------------------------------------- FORMAT IFMT=2
100   Format (/, 80('*'), /, 2F21.15, ' RYD   TAU FOR IT=', I2, '  IQ=', I2, &
        :, '  CPA:', I2, F15.6)
110   Format (/, 80('*'), /, 2F21.15, ' RYD   TAU-C M-C  FOR IQ=', I2, :, &
        '  CPA:', I2, F15.6)
120   Format (2I5, 1P, 4E22.14)

    End Subroutine
