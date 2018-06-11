    Subroutine wrmoms(krel, natyp, nspin, texts, textl, textns, charge, muorb, &
      lmaxd, lmaxd1)
      Use mod_datatypes, Only: dp

      Implicit None

! Dummy arguments
      Integer :: krel, lmaxd, lmaxd1, natyp, nspin
      Character (Len=5) :: textns
      Real (Kind=dp) :: charge(0:lmaxd1, natyp, 2)
      Character (Len=4) :: textl(0:6)
      Character (Len=7) :: texts(3)

! Local variables
      Real (Kind=dp) :: chtot(natyp), chval
      Real (Kind=dp) :: muorb(0:lmaxd1+1, 3, natyp)
      Real (Kind=dp) :: muspin(natyp, 0:lmaxd1+1), mutot(natyp), &
        sumch(natyp, 2)
      Character (Len=80) :: fmt1, fmt2, fmt31, fmt32
      Integer :: is, ispin, it, l, lf1, lf2

      Write (1337, *)

      If ((krel==1) .Or. (nspin==2)) Then
        Write (1337, '(78("#"))')
        Write (1337, 100)
        Write (1337, '(78("#"))')
      Else
        Write (1337, '(44("#"))')
        Write (1337, 110)
        Write (1337, '(44("#"))')
      End If

      Write (1337, *)
      Write (1337, 120, Advance='no')
      Write (6, 130, Advance='no')
      Do it = 1, natyp
        muspin(it, lmaxd1+1) = 0E0_dp
        sumch(it, 1) = 0E0_dp
        sumch(it, 2) = 0E0_dp
        Do l = 0, lmaxd1
          Do ispin = 1, nspin
            sumch(it, ispin) = sumch(it, ispin) + charge(l, it, ispin)
          End Do
          muspin(it, l) = charge(l, it, 2) - charge(l, it, 1)
          muspin(it, lmaxd1+1) = muspin(it, lmaxd1+1) + muspin(it, l)
        End Do
        chtot(it) = sumch(it, 1) + sumch(it, 2)
      End Do

      If (krel==1) Then
        Do it = 1, natyp
          mutot(it) = muspin(it, lmaxd1+1) + muorb(lmaxd1+1, 3, it)
        End Do
      End If

      is = 0
      If (nspin==1) is = is + 2
      Do ispin = 1, nspin
        is = is + 1
        Write (1337, 140, Advance='no') texts(is)
        Write (6, 140, Advance='no') texts(is)
      End Do

      If (krel==1) Then
        Write (1337, 150, Advance='no')
        Write (6, 150, Advance='no')
        Write (1337, 160)
        Write (6, 160)
      Else
        If (nspin==2) Write (1337, 150, Advance='no')
        If (nspin==2) Write (6, 150, Advance='no')
        Write (1337, *)
        Write (6, *)
      End If

      Write (1337, '(3X,26("="))', Advance='no')
      If (krel==1) Then
        Write (1337, '(46("="))')
      Else
        If (nspin==2) Write (1337, '(23("="))', Advance='no')
        Write (1337, *)
      End If

      fmt1 = '(4X,I3,2X,A4,2(F12.8),2X,F8.4'
      fmt2 = '(9X,A4,2(F12.8),2X,F8.4'
      fmt31 = '(4X,I3,2X,A4,F12.8)'
      fmt32 = '(9X,A4,F12.8)'
      lf1 = 30
      lf2 = 24

      If (krel==1) Then
        fmt1 = fmt1(1:lf1) // ',2X,3F8.4)'
        fmt2 = fmt2(1:lf2) // ',2X,3F8.4)'
      Else
        If (nspin==2) Then
          fmt1 = fmt1(1:lf1) // ')'
          fmt2 = fmt2(1:lf2) // ')'
        Else
          fmt1 = fmt31
          fmt2 = fmt32
        End If
      End If

      Do it = 1, natyp
        If (krel==1) Then
          Write (1337, Fmt=fmt1) it, textl(0), (charge(0,it,ispin), ispin=1, &
            nspin), muspin(it, 0), muorb(0, 3, it), (muorb(0,ispin,it), ispin= &
            1, nspin)
        Else
          If (nspin==2) Then
            Write (1337, Fmt=fmt1) it, textl(0), (charge(0,it,ispin), ispin=1, &
              nspin), muspin(it, 0)
          Else
            Write (1337, Fmt=fmt1) it, textl(0), charge(0, it, 1)
          End If
        End If

        Do l = 1, lmaxd
          If (krel==1) Then
            Write (1337, Fmt=fmt2) textl(l), (charge(l,it,ispin), ispin=1, &
              nspin), muspin(it, l), muorb(l, 3, it), &
              (muorb(l,ispin,it), ispin=1, nspin)
          Else
            If (nspin==2) Then
              Write (1337, Fmt=fmt2) textl(l), (charge(l,it,ispin), ispin=1, &
                nspin), muspin(it, l)
            Else
              Write (1337, Fmt=fmt2) textl(l), charge(l, it, 1)
            End If
          End If

        End Do

        If (krel==1) Then
          Write (1337, Fmt=fmt2) textns, (charge(lmaxd1,it,ispin), ispin=1, &
            nspin), muspin(it, lmaxd1), muorb(lmaxd1, 3, it), &
            (muorb(lmaxd1,ispin,it), ispin=1, nspin)
        Else
          If (nspin==2) Then
            Write (1337, Fmt=fmt2) textns, (charge(lmaxd1,it,ispin), ispin=1, &
              nspin), muspin(it, lmaxd1)
          Else
            Write (1337, Fmt=fmt2) textns, charge(lmaxd1, it, 1)
          End If
        End If

        Write (1337, '(10x,19("-"))', Advance='no')
        If (krel==1) Then
          Write (1337, '(44("-"))')
          Write (1337, Fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), &
            muspin(it, lmaxd1+1), muorb(lmaxd1+1, 3, it), &
            (muorb(lmaxd1+1,ispin,it), ispin=1, nspin)
          Write (6, Fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), &
            muspin(it, lmaxd1+1), muorb(lmaxd1+1, 3, it), &
            (muorb(lmaxd1+1,ispin,it), ispin=1, nspin)
          Write (1337, '(25X,F12.8,12X,F8.4)') chtot(it), mutot(it)
        Else
          If (nspin==2) Then
            Write (1337, '(17("-"))')
            Write (1337, Fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), &
              muspin(it, lmaxd1+1)
            Write (6, Fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), &
              muspin(it, lmaxd1+1)
            Write (1337, '(25X,F12.8)') chtot(it)
          Else
            Write (1337, *)
            Write (1337, Fmt=fmt2) ' TOT', sumch(it, 1)
            Write (6, Fmt=fmt2) ' TOT', sumch(it, 1)
          End If
        End If

        If (it/=natyp) Then
          Write (1337, '(3X,26("="))', Advance='no')
          If (krel==1) Then
            Write (1337, '(40("="))')
          Else
            If (nspin==2) Write (1337, '(17("="))', Advance='no')
            Write (1337, *)
          End If
        End If
      End Do

      Write (1337, *)
      If ((krel==1) .Or. (nspin==2)) Then
        Write (1337, '(78("#"))')
      Else
        Write (1337, '(44("#"))')
      End If
      Write (1337, *)

      chval = 0.E0_dp
      Do it = 1, natyp
        chval = chval + chtot(it)
      End Do
      Write (1337, *) 'Sum of valence charges of atoms (local summation)', &
        chval

      Return

100   Format (15X, 'l-decomposed valence charges and magnetic moments')
110   Format (8X, 'l-decomposed valence charges')
120   Format (3X, 'ATOM      ')
130   Format (3X, '          ')
140   Format (2X, 'Ne ', A7)
150   Format ('    m_spin')
160   Format ('    m_orb   spin dn  spin up')
    End Subroutine
