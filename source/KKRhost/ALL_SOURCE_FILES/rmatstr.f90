module mod_rmatstr

contains

subroutine rmatstr(str, lstr, a, n, m, mlin, mcol, tolp, nfil)
  use :: mod_datatypes, only: dp
  ! ********************************************************************
  ! *                                                                  *
  ! *   writes structure of REAL      NxN   matrix   A                 *
  ! *                                                                  *
  ! *   M           is the actual array - size used for   A            *
  ! *   MLIN/COL    MODE for line and comlun indexing                  *
  ! *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
  ! *   TOL         tolerance for difference                           *
  ! *                                                                  *
  ! *                                                         03/01/96 *
  ! ********************************************************************
  implicit none

  integer, parameter :: ndifmax = 250
  integer, intent (in) :: lstr, m, n, mlin, mcol, nfil
  real (kind=dp), intent (in) :: tolp

  character (len=lstr) :: str
  character (len=150) :: fmt1, fmt2, fmt3, fmt4
  character (len=1) :: ctab(0:ndifmax), vz(-1:+1)
  integer :: iw(m), ilsep(20)
  integer :: icauc, iczuc, natoz, icalc, icilc, k, nk, nm, nm2, nm1, nm3, ic0, &
    n3, n2, n1, lf, l3, mm, nsl, il, i, nnon0, nd, nalf, j, id, icnd, isl
  real (kind=dp) :: tol
  real (kind=dp) :: a(m, m), cnum, ca, cb, dtab(0:ndifmax)
  logical :: csmall, csame

  save :: vz
  data vz/'-', ' ', ' '/

  csmall(cnum) = abs(cnum*tol) < 1.0e0_dp

  csame(ca, cb) = csmall(1.0e0_dp-ca/cb)

  tol = 1.0e0_dp/tolp

  icauc = ichar('A')
  iczuc = ichar('Z')
  natoz = iczuc - icauc + 1
  icalc = ichar('a')
  icilc = ichar('i')

  ! ---------------------------------------------------------------- header
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
  if (mcol==0) then
    fmt1 = fmt1(1:lf) // char(ic0+n3) // char(ic0+n2) // char(ic0+n1) // &
      '( 2A1),''|'',I3)'
    fmt2 = fmt2(1:lf) // char(ic0+n3) // char(ic0+n2) // char(ic0+n1) // &
      '(''--''),''|'',I3)'
    fmt3 = fmt3(1:lf) // '60(2X,I2))'
    fmt4 = fmt4(1:lf) // '60(I2,2X))'
    lf = 21
  else
    if (mcol==1) then
      nk = nint(sqrt(real(n,kind=dp)))
    else if (mcol==2) then
      nk = nint(sqrt(real(n/2,kind=dp)))
    else if (mcol==3) then
      nk = 2*nint(sqrt(real(n/2,kind=dp))) - 1
    end if
    do k = 1, nk
      if (mcol<=2) then
        nm = 2*k - 1
      else
        nm = 2*((k+1)/2)
      end if
      nm2 = nm/10
      nm1 = nm - nm2*10
      nm3 = nm/2
      fmt1 = fmt1(1:lf) // char(ic0+nm2) // char(ic0+nm1) // '( 2A1),''|'','
      fmt2 = fmt2(1:lf) // char(ic0+nm2) // char(ic0+nm1) // '(''--''),''|'','

      if (mcol<=2) then
        do mm = 1, nm
          if (mod(mm,2)==mod(k,2)) then
            fmt3 = fmt3(1:l3) // '2X,'
            fmt4 = fmt4(1:l3) // 'I2,'
          else
            fmt3 = fmt3(1:l3) // 'I2,'
            fmt4 = fmt4(1:l3) // '2X,'
          end if
          l3 = l3 + 3
        end do
        fmt3 = fmt3(1:l3) // '''|'','
        fmt4 = fmt4(1:l3) // '''|'','
        l3 = l3 + 4
      else
        fmt3 = fmt3(1:lf) // char(ic0+nm3) // '(2X,I2),''|'','
        fmt4 = fmt4(1:lf) // char(ic0+nm3) // '(I2,2X),''|'','
        l3 = l3 + 13
      end if
      lf = lf + 13
    end do
    if (mcol==2) then
      fmt1 = fmt1(1:lf) // fmt1(12:lf)
      fmt2 = fmt2(1:lf) // fmt2(12:lf)
      fmt3 = fmt3(1:l3) // fmt4(12:l3)
      fmt4 = fmt4(1:l3) // fmt3(12:l3)
      lf = 2*lf - 11
    end if
    fmt1 = fmt1(1:lf) // 'I3)'
    fmt2 = fmt2(1:lf) // 'I3)'
    fmt3 = fmt3(1:l3) // 'I3)'
    fmt4 = fmt4(1:l3) // 'I3)'
  end if
  if (mlin==0) then
    nsl = 1
    ilsep(1) = n
  else if (mlin==1) then
    nsl = nint(sqrt(real(n,kind=dp)))
    do il = 1, nsl
      ilsep(il) = il**2
    end do
  else if (mlin==2) then
    nsl = nint(sqrt(real(n/2,kind=dp)))
    do il = 1, nsl
      ilsep(il) = il**2
    end do
    do il = 1, nsl
      ilsep(nsl+il) = ilsep(nsl) + il**2
    end do
    nsl = 2*nsl
  else if (mlin==3) then
    nsl = 2*nint(sqrt(real(n/2,kind=dp))) - 1
    ilsep(1) = 2
    do k = 2, nsl
      ilsep(k) = ilsep(k-1) + 2*((k+1)/2)
    end do
  end if


  write (nfil, 110) str(1:lstr)
  write (nfil, fmt3)(i, i=2, n, 2)
  write (nfil, fmt4)(i, i=1, n, 2)
  write (nfil, fmt=fmt2)
  ! ------------------------------------------------------------ header end
  nnon0 = 0
  nd = 0
  nalf = 0
  ctab(0) = ' '
  dtab(0) = 9999e0_dp

  do i = 1, n
    do j = 1, n
      if (.not. csmall(a(i,j))) then
        nnon0 = nnon0 + 1
        do id = 1, nd
          if (csame(a(i,j),+dtab(id))) then
            iw(j) = +id
            go to 100
          end if
          if (csame(a(i,j),-dtab(id))) then
            iw(j) = -id
            go to 100
          end if
        end do
        ! ----------------------------------------------------------- new
        ! element
        nd = nd + 1
        if (nd>ndifmax) then
          write (nfil, '(''nd>array size ndifmax='', i3)') ndifmax
          stop
        end if
        iw(j) = nd
        dtab(nd) = a(i, j)
        if (abs(dtab(nd)-1.0e0_dp)*tol<1.0e0_dp) then
          ctab(nd) = '1'
        else if (abs(dtab(nd)+1.0e0_dp)*tol<1.0e0_dp) then
          dtab(nd) = +1.0e0_dp
          ctab(nd) = '1'
          iw(j) = -nd
        else
          nalf = nalf + 1
          if (nalf<=natoz) then
            ctab(nd) = char(icauc+nalf-1)
          else
            icnd = icalc + nalf - natoz - 1
            if (icnd<icilc) then
              ctab(nd) = char(icnd)
            else
              ctab(nd) = char(icnd+1)
            end if
          end if
        end if
100     continue
      else
        iw(j) = 0
      end if
    end do
    ! ------------------------------------------------------------ write line
    write (nfil, fmt=fmt1) i, (vz(sign(1,iw(j))), ctab(abs(iw(j))), j=1, n), i


    do isl = 1, nsl
      if (i==ilsep(isl)) write (nfil, fmt=fmt2)
    end do
  end do

  ! ------------------------------------------------------------------ foot

  write (nfil, fmt4)(i, i=1, n, 2)
  write (nfil, fmt3)(i, i=2, n, 2)

  write (nfil, 120)(id, ctab(id), dtab(id), id=1, nd)
  write (nfil, 130) nnon0, n*n - nnon0

110 format (/, 8x, a, /)
120 format (/, 8x, 'symbols used:', /, (8x,i3,3x,a1,2x,f20.12))
130 format (/, 8x, 'elements <> 0:', i4, /, 8x, 'elements  = 0:', i4)
  return
end subroutine rmatstr

end module mod_rmatstr
