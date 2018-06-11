    Subroutine cmatstr(str, lstr, a, n, m, mlin, mcol, ijq, tolp, k_fmt_fil)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   writes structure of COMPLEX   NxN   matrix   A                 *
!   *                                                                  *
!   *   M           is the actual array - size used for   A            *
!   *   MLIN/COL    MODE for line and column indexing                  *
!   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
!   *   TOL         tolerance for difference                           *
!   *   IJQ         if IJQ > 1000    pick  IQ-JQ-block matrix          *
!   *               assuming  IJQ = IQ*1000 + JQ                       *
!   *               else: no IQ-JQ-indexing                            *
!   *   K_FMT_FIL   output channel                                     *
!   *               a negative sign suppresses table at the end        *
!   *                                                                  *
!   *   any changes should be done in RMATSTR as well !!!!!!!!!!!!!!!  *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: ci
      Parameter (ci=(0.0E0_dp,1.0E0_dp))

! Dummy arguments
      Integer :: ijq, k_fmt_fil, lstr, m, mcol, mlin, n
      Character (Len=lstr) :: str
      Real (Kind=dp) :: tolp
      Complex (Kind=dp) :: a(m, m)

! Local variables
      Complex (Kind=dp) :: b(n, n), ca, cb, arg, dtab(0:n*n)
      Character :: char
      Logical :: same, small
      Character (Len=1) :: ctab(0:n*n), vz(-1:+1)
      Real (Kind=dp) :: dble
      Character (Len=150) :: fmt1, fmt2, fmt3, fmt4
      Integer :: i, i1, ic0, id, il, ilsep(20), ipt(218), iq, isl, iw(m), j, &
        j0, jp, jq, k, l3, lf, mm, n1, n2, n3, nc, nd, nfil, nk, nm, nm1, nm2, &
        nm3, nnon0, nsl
      Integer :: ichar, isign, nint
      Real (Kind=dp) :: tol

      Data vz/'-', ' ', ' '/

      small(arg) = abs(arg*tol) < 1.0E0_dp

      same(ca, cb) = small(1.0E0_dp-ca/cb)

      nfil = abs(k_fmt_fil)

      tol = 1.0E0_dp/tolp

!----------------------------------------------- set block indices IQ JQ

      If (ijq>1000) Then
        iq = ijq/1000
        jq = ijq - iq*1000
        If (iq*n>m .Or. iq*n>m) Then
          Write (1337, 120) ijq, iq, jq, iq*n, jq*n, n, m
          Return
        End If
      Else
        iq = 1
        jq = 1
      End If

!----------------------------------------------------- copy matrix block

      j0 = n*(jq-1)
      Do j = 1, n
        i1 = n*(iq-1) + 1
        jp = j0 + j
        Call zcopy(n, a(i1,jp), 1, b(1,j), 1)
      End Do

!------------------------------------------------ set up character table

      nc = 0
      Do i = 1, 26
        nc = nc + 1
        ipt(nc) = 62 + i
      End Do
      Do i = 1, 8
        nc = nc + 1
        ipt(nc) = 96 + i
      End Do
      Do i = 10, 26
        nc = nc + 1
        ipt(nc) = 96 + i
      End Do
      Do i = 191, 218
        nc = nc + 1
        ipt(nc) = i
      End Do
      Do i = 35, 38
        nc = nc + 1
        ipt(nc) = i
      End Do
      Do i = 40, 42
        nc = nc + 1
        ipt(nc) = i
      End Do
      Do i = 91, 93
        nc = nc + 1
        ipt(nc) = i
      End Do

!---------------------------------------------------------------- header
      ic0 = ichar('0')
      n3 = n/100
      n2 = n/10 - n3*10
      n1 = n - n2*10 - n3*100

      If (n<=18) Then
        fmt1 = '(8X,I3,''|'','
        fmt2 = '( 9X,''--|'','
        fmt3 = '( 9X,'' #|'','
        fmt4 = '( 9X,''  |'','
      Else
        fmt1 = '(   I4,''|'','
        fmt2 = '( 2X,''--|'','
        fmt3 = '( 2X,'' #|'','
        fmt4 = '( 2X,''  |'','
      End If

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
          nk = nint(sqrt(dble(n)))
        Else If (mcol==2) Then
          nk = nint(sqrt(dble(n/2)))
        Else If (mcol==3) Then
          nk = 2*nint(sqrt(dble(n/2))) - 1
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
          fmt3 = fmt3(1:l3) // fmt3(12:l3)
          fmt4 = fmt4(1:l3) // fmt4(12:l3)
          lf = 2*lf - 11
          l3 = 2*l3 - 11
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
        nsl = nint(sqrt(dble(n)))
        Do il = 1, nsl
          ilsep(il) = il**2
        End Do
      Else If (mlin==2) Then
        nsl = nint(sqrt(dble(n/2)))
        Do il = 1, nsl
          ilsep(il) = il**2
        End Do
        Do il = 1, nsl
          ilsep(nsl+il) = ilsep(nsl) + il**2
        End Do
        nsl = 2*nsl
      Else If (mlin==3) Then
        nsl = 2*nint(sqrt(dble(n/2))) - 1
        ilsep(1) = 2
        Do k = 2, nsl
          ilsep(k) = ilsep(k-1) + 2*((k+1)/2)
        End Do
      End If


      Write (nfil, 110) str(1:lstr)
      If (ijq>1000) Write (nfil, 130) iq, jq
      Write (nfil, fmt3)(i, i=2, n, 2)
      Write (nfil, fmt4)(i, i=1, n, 2)
      Write (nfil, Fmt=fmt2)
!------------------------------------------------------------ header end
      nnon0 = 0
      nd = 0
      ctab(0) = ' '
      dtab(0) = 9999E0_dp

      Do i = 1, n
        Do j = 1, n
          If (.Not. small(b(i,j))) Then
            nnon0 = nnon0 + 1
            Do id = 1, nd
              If (same(b(i,j),+dtab(id))) Then
                iw(j) = +id
                Go To 100
              End If
              If (same(b(i,j),-dtab(id))) Then
                iw(j) = -id
                Go To 100
              End If
            End Do
!----------------------------------------------------------- new element
            nd = nd + 1
            iw(j) = nd
            dtab(nd) = b(i, j)
            If (abs(dtab(nd)-1.0E0_dp)*tol<1.0E0_dp) Then
              ctab(nd) = '1'
            Else If (abs(dtab(nd)+1.0E0_dp)*tol<1.0E0_dp) Then
              dtab(nd) = +1.0E0_dp
              ctab(nd) = '1'
              iw(j) = -nd
            Else If (abs(dtab(nd)-ci)*tol<1.0E0_dp) Then
              ctab(nd) = 'i'
            Else If (abs(dtab(nd)+ci)*tol<1.0E0_dp) Then
              dtab(nd) = +ci
              ctab(nd) = 'i'
              iw(j) = -nd
            Else
              ctab(nd) = char(ipt(1+mod((nd+1),nc)))
            End If
          Else
            iw(j) = 0
          End If
100     End Do
!------------------------------------------------------------ write line
        Write (nfil, Fmt=fmt1) i, (vz(isign(1,iw(j))), ctab(abs(iw(j))), j=1, &
          n), i

        Do isl = 1, nsl
          If (i==ilsep(isl)) Write (nfil, Fmt=fmt2)
        End Do
      End Do

!------------------------------------------------------------------ foot

      Write (nfil, fmt4)(i, i=1, n, 2)
      Write (nfil, fmt3)(i, i=2, n, 2)

      If (k_fmt_fil>0) Then
        Write (nfil, 140)(id, ctab(id), dtab(id), id=1, nd)
        Write (nfil, 150) nnon0, tolp, n*n - nnon0, tolp
      Else
        Write (nfil, *) ' '
      End If

110   Format (/, 8X, A, /)
120   Format (/, 1X, 79('*'), /, 10X, 'inconsistent call of <CMATSTR>', /, &
        10X, 'argument IJQ =', I8, '  implies IQ=', I3, '   JQ=', I3, /, 10X, &
        'IQ*N=', I6, ' > M   or   JQ*N=', I6, ' > M   for N =', I4, ' M=', I4, &
        /, 1X, 79('*'), /)
130   Format (8X, 'IQ-JQ-block  for  IQ = ', I3, '   JQ = ', I3, /)
140   Format (/, 8X, 'symbols used:', /, (8X,I3,3X,A1,2X,2F20.12))
150   Format (/, 8X, I5, ' elements   >', 1P, E9.1, /, 8X, I5, &
        ' elements   <', 1P, E9.1, /)
    End Subroutine
