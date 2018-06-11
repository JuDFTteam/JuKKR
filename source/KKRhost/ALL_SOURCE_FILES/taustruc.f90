    Subroutine taustruct(drot, nsym, symunitary, nkm, nq, nqmax, nkmmax, &
      iprint, irel)
!   ********************************************************************
!   *                                                                  *
!   *   find the structure of the site-diagonal TAU - matrices  TAUQ   *
!   *                                                                  *
!   ********************************************************************
      Use mod_datatypes
      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: c0
      Parameter (c0=(0.0E0_dp,0.0E0_dp))

! Dummy arguments
      Integer :: iprint, irel, nkm, nkmmax, nq, nqmax, nsym
      Complex (Kind=dp) :: drot(nkmmax, nkmmax, 48)
      Logical :: symunitary(48)

! Local variables
      Real (Kind=dp) :: abst, x
      Integer :: i, i0, imweight, iq, isym, iw, iwr, j, k, l, lin, nelmt, &
        nkmq(nqmax), nkmtop, nlin, non0(nqmax)
      Complex (Kind=dp) :: st(nkmmax, nkmmax), tauk(nkmmax, nkmmax, nqmax)

      Do iq = 1, nq
        nkmq(iq) = nkm
      End Do

      imweight = 0
      nelmt = 0
      nlin = 0
      iw = 6

      Do iq = 1, nq
        non0(iq) = 0
        nkmtop = nkmq(iq)

        If (iprint>0) Write (1337, 120) iq
        Do i = 1, nkmtop
          Do j = 1, nkmtop
            st(i, j) = 0.0E0_dp

            Call cinit(nkmmax*nkmmax*nqmax, tauk)

            Do isym = 1, nsym
              i0 = iq

              If (symunitary(isym)) Then
                Do l = 1, nkmtop
                  Do k = 1, nkmtop
                    tauk(k, l, i0) = tauk(k, l, i0) + &
                      drot(i, k, isym)*conjg(drot(j,l,isym))
                  End Do
                End Do
              Else
                Do l = 1, nkmtop
                  Do k = 1, nkmtop
                    tauk(l, k, i0) = tauk(l, k, i0) + &
                      drot(i, k, isym)*conjg(drot(j,l,isym))
                  End Do
                End Do
              End If
            End Do

            lin = 0
            iwr = 0
            Do k = 1, nkmq(iq)
              Do l = 1, nkmq(iq)
                abst = abs(tauk(k,l,iq))
                st(i, j) = st(i, j) + abst
                If (abst>1E-8_dp) Then
                  If (aimag(tauk(k,l,iq))>1E-5_dp) Then
                    If (iprint>0) Write (1337, *) ' Im(Weight) > 1D-5 ', i, j, &
                      k, l
                    imweight = 1
                  End If
                  x = real(tauk(k,l,iq))/real(nsym, kind=dp)

                  If (iprint>1) Then
                    If (iwr==0) Then
                      iwr = 1
                      Write (iw, 100) i, j, iq, x, k + (iq-1)*nkm, &
                        l + (iq-1)*nkm
                    Else
                      Write (iw, 110) x, k + (iq-1)*nkm, l + (iq-1)*nkm
                    End If
                  End If
                  lin = lin + 1
                End If

              End Do
            End Do

            If (lin>0) Then
              nlin = nlin + lin
              nelmt = nelmt + 1
              non0(iq) = non0(iq) + 1
            End If

            If (abs(st(i,j))>1E-5_dp) st(i, j) = 2

          End Do
        End Do

        If (iprint>1) Call cmatstr('TAU-MAT', 7, st, nkmtop, nkmmax, irel, &
          irel, 0, 1E-8_dp, 6)
      End Do

      Write (1337, 130) nelmt, (non0(iq), iq=1, nq)
      Write (1337, 140) nlin

      If (imweight/=0) Write (1337, 150)

!-----------------------------------------------------------------------

      Return
100   Format ('     TAUQ(', I2, ',', I2, ',', I2, ') =  ', '   ', F10.4, &
        ' * <', I3, '|T(K)|', I3, '>')
110   Format (23X, ' + ', F10.4, ' * <', I3, '|T(K)|', I3, '>')
120   Format (/, /, &
        ' ===========================================================', /, &
        '   structure of  TAU-matrix   INT <i|t(k)|j>     IQ=', I3, /, &
        ' ===========================================================', /)
130   Format (/, 5X, 'non-0 TAU-elements          ', I5, '   Q:', 80I4)
140   Format (5X, 'terms to sum up             ', I5, /)
150   Format (/, 5X, 50('#'), /, 5X, 'WARNING: complex TAU weights found', /, &
        5X, 'this may occur for rotated magnetic moments', /, 5X, &
        'relevant only for tetrahedron BZ-integration', /, 5X, 50('#'), /)
    End Subroutine
