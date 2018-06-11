    Subroutine rmatstr(str, lstr, a, n, m, mlin, mcol, tolp, nfil)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   writes structure of REAL      NxN   matrix   A                 *
!   *                                                                  *
!   *   M           is the actual array - size used for   A            *
!   *   MLIN/COL    MODE for line and comlun indexing                  *
!   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
!   *   TOL         tolerance for difference                           *
!   *                                                                  *
!   *                                                         03/01/96 *
!   ********************************************************************
      Implicit None

      Integer, Parameter :: ndifmax = 250
      Integer, Intent (In) :: lstr, m, n, mlin, mcol, nfil
      Real (Kind=dp), Intent (In) :: tolp

      Character (Len=lstr) :: str
      Character (Len=150) :: fmt1, fmt2, fmt3, fmt4
      Character (Len=1) :: ctab(0:ndifmax), vz(-1:+1)
      Integer :: iw(m), ilsep(20)
      Integer :: icauc, iczuc, natoz, icalc, icilc, k, nk, nm, nm2, nm1, nm3, &
        ic0, n3, n2, n1, lf, l3, mm, nsl, il, i, nnon0, nd, nalf, j, id, icnd, &
        isl
      Real (Kind=dp) :: tol
      Real (Kind=dp) :: a(m, m), cnum, ca, cb, dtab(0:ndifmax)
      Logical :: csmall, csame

      Save :: vz
      Data vz/'-', ' ', ' '/

      csmall(cnum) = abs(cnum*tol) < 1.0E0_dp

      csame(ca, cb) = csmall(1.0E0_dp-ca/cb)

      tol = 1.0E0_dp/tolp

      icauc = ichar('A')
      iczuc = ichar('Z')
      natoz = iczuc - icauc + 1
      icalc = ichar('a')
      icilc = ichar('i')

!---------------------------------------------------------------- header
      ic0 = ichar('0')
      n3 = n/100
      n2 = n/10 - n3*10
      n1 = n - n2*10 - n3*100

      fmt1 = '(8X,I3,''|'','
      fmt2 = '( 9X,''--|'','
      fmt3 = '( 9X,'' #|'','
      fmt4 = '( 9X,''  |'','

      lf = 11
      l3 = 11
      If (mcol==0) Then
        fmt1 = fmt1(1:lf) // char(ic0+n3) // char(ic0+n2) // char(ic0+n1) // &
          '( 2A1),''|'',I3)'
        fmt2 = fmt2(1:lf) // char(ic0+n3) // char(ic0+n2) // char(ic0+n1) // &
          '(''--''),''|'',I3)'
        fmt3 = fmt3(1:lf) // '60(2X,I2))'
        fmt4 = fmt4(1:lf) // '60(I2,2X))'
        lf = 21
      Else
        If (mcol==1) Then
          nk = nint(sqrt(real(n,kind=dp)))
        Else If (mcol==2) Then
          nk = nint(sqrt(real(n/2,kind=dp)))
        Else If (mcol==3) Then
          nk = 2*nint(sqrt(real(n/2,kind=dp))) - 1
        End If
        Do k = 1, nk
          If (mcol<=2) Then
            nm = 2*k - 1
          Else
            nm = 2*((k+1)/2)
          End If
          nm2 = nm/10
          nm1 = nm - nm2*10
          nm3 = nm/2
          fmt1 = fmt1(1:lf) // char(ic0+nm2) // char(ic0+nm1) // &
            '( 2A1),''|'','
          fmt2 = fmt2(1:lf) // char(ic0+nm2) // char(ic0+nm1) // &
            '(''--''),''|'','

          If (mcol<=2) Then
            Do mm = 1, nm
              If (mod(mm,2)==mod(k,2)) Then
                fmt3 = fmt3(1:l3) // '2X,'
                fmt4 = fmt4(1:l3) // 'I2,'
              Else
                fmt3 = fmt3(1:l3) // 'I2,'
                fmt4 = fmt4(1:l3) // '2X,'
              End If
              l3 = l3 + 3
            End Do
            fmt3 = fmt3(1:l3) // '''|'','
            fmt4 = fmt4(1:l3) // '''|'','
            l3 = l3 + 4
          Else
            fmt3 = fmt3(1:lf) // char(ic0+nm3) // '(2X,I2),''|'','
            fmt4 = fmt4(1:lf) // char(ic0+nm3) // '(I2,2X),''|'','
            l3 = l3 + 13
          End If
          lf = lf + 13
        End Do
        If (mcol==2) Then
          fmt1 = fmt1(1:lf) // fmt1(12:lf)
          fmt2 = fmt2(1:lf) // fmt2(12:lf)
          fmt3 = fmt3(1:l3) // fmt4(12:l3)
          fmt4 = fmt4(1:l3) // fmt3(12:l3)
          lf = 2*lf - 11
        End If
        fmt1 = fmt1(1:lf) // 'I3)'
        fmt2 = fmt2(1:lf) // 'I3)'
        fmt3 = fmt3(1:l3) // 'I3)'
        fmt4 = fmt4(1:l3) // 'I3)'
      End If
      If (mlin==0) Then
        nsl = 1
        ilsep(1) = n
      Else If (mlin==1) Then
        nsl = nint(sqrt(real(n,kind=dp)))
        Do il = 1, nsl
          ilsep(il) = il**2
        End Do
      Else If (mlin==2) Then
        nsl = nint(sqrt(real(n/2,kind=dp)))
        Do il = 1, nsl
          ilsep(il) = il**2
        End Do
        Do il = 1, nsl
          ilsep(nsl+il) = ilsep(nsl) + il**2
        End Do
        nsl = 2*nsl
      Else If (mlin==3) Then
        nsl = 2*nint(sqrt(real(n/2,kind=dp))) - 1
        ilsep(1) = 2
        Do k = 2, nsl
          ilsep(k) = ilsep(k-1) + 2*((k+1)/2)
        End Do
      End If


      Write (nfil, 110) str(1:lstr)
      Write (nfil, fmt3)(i, i=2, n, 2)
      Write (nfil, fmt4)(i, i=1, n, 2)
      Write (nfil, Fmt=fmt2)
!------------------------------------------------------------ header end
      nnon0 = 0
      nd = 0
      nalf = 0
      ctab(0) = ' '
      dtab(0) = 9999E0_dp

      Do i = 1, n
        Do j = 1, n
          If (.Not. csmall(a(i,j))) Then
            nnon0 = nnon0 + 1
            Do id = 1, nd
              If (csame(a(i,j),+dtab(id))) Then
                iw(j) = +id
                Go To 100
              End If
              If (csame(a(i,j),-dtab(id))) Then
                iw(j) = -id
                Go To 100
              End If
            End Do
!----------------------------------------------------------- new element
            nd = nd + 1
            If (nd>ndifmax) Then
              Write (nfil, '(''nd>array size ndifmax='', i3)') ndifmax
              Stop
            End If
            iw(j) = nd
            dtab(nd) = a(i, j)
            If (abs(dtab(nd)-1.0E0_dp)*tol<1.0E0_dp) Then
              ctab(nd) = '1'
            Else If (abs(dtab(nd)+1.0E0_dp)*tol<1.0E0_dp) Then
              dtab(nd) = +1.0E0_dp
              ctab(nd) = '1'
              iw(j) = -nd
            Else
              nalf = nalf + 1
              If (nalf<=natoz) Then
                ctab(nd) = char(icauc+nalf-1)
              Else
                icnd = icalc + nalf - natoz - 1
                If (icnd<icilc) Then
                  ctab(nd) = char(icnd)
                Else
                  ctab(nd) = char(icnd+1)
                End If
              End If
            End If
100         Continue
          Else
            iw(j) = 0
          End If
        End Do
!------------------------------------------------------------ write line
        Write (nfil, Fmt=fmt1) i, (vz(sign(1,iw(j))), ctab(abs(iw(j))), j=1, n &
          ), i


        Do isl = 1, nsl
          If (i==ilsep(isl)) Write (nfil, Fmt=fmt2)
        End Do
      End Do

!------------------------------------------------------------------ foot

      Write (nfil, fmt4)(i, i=1, n, 2)
      Write (nfil, fmt3)(i, i=2, n, 2)

      Write (nfil, 120)(id, ctab(id), dtab(id), id=1, nd)
      Write (nfil, 130) nnon0, n*n - nnon0

110   Format (/, 8X, A, /)
120   Format (/, 8X, 'symbols used:', /, (8X,I3,3X,A1,2X,F20.12))
130   Format (/, 8X, 'elements <> 0:', I4, /, 8X, 'elements  = 0:', I4)
      Return
    End Subroutine
