!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_cmatstr

  private
  public :: cmatstr

contains

  !-------------------------------------------------------------------------------
  !> Summary: Write complex matrix to file
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Writes structure of COMPLEX   NxN   matrix   A             
  !>                                                            
  !> M           is the actual array - size used for   A        
  !> MLIN/COL    MODE for line and column indexing              
  !>             0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)
  !> TOL         tolerance for difference                       
  !> IJQ         if IJQ > 1000    pick  IQ-JQ-block matrix      
  !>             assuming  IJQ = IQ*1000 + JQ                   
  !>             else: no IQ-JQ-indexing                        
  !> K_FMT_FIL   output channel                                 
  !>             a negative sign suppresses table at the end    
  !>           
  !> @note                                                 
  !> any changes should be done in RMATSTR as well !!!!!!!!!!!!!
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine cmatstr(str, lstr, a, n, m, mlin, mcol, ijq, tolp, k_fmt_fil)

    use :: mod_constants, only: cone
    use :: mod_datatypes, only: dp
    implicit none

    ! Dummy arguments
    integer :: ijq, k_fmt_fil, lstr, m, mcol, mlin, n
    character (len=lstr) :: str
    real (kind=dp) :: tolp
    complex (kind=dp) :: a(m, m)

    ! Local variables
    complex (kind=dp) :: b(n, n), ca, cb, arg, dtab(0:n*n)
    logical :: same, small
    character (len=1) :: ctab(0:n*n), vz(-1:+1)
    character (len=150) :: fmt1, fmt2, fmt3, fmt4
    integer :: i, i1, ic0, id, il, ilsep(20), ipt(218), iq, isl, iw(m), j, j0, jp, jq, k, l3, lf, mm, n1, n2, n3, nc, nd, nfil, nk, nm, nm1, nm2, nm3, nnon0, nsl
    real (kind=dp) :: tol

    data vz/'-', ' ', ' '/

    small(arg) = abs(arg*tol) < 1.0e0_dp

    same(ca, cb) = small(1.0e0_dp-ca/cb)

    nfil = abs(k_fmt_fil)

    tol = 1.0e0_dp/tolp

    ! ----------------------------------------------- set block indices IQ JQ

    if (ijq>1000) then
      iq = ijq/1000
      jq = ijq - iq*1000
      if (iq*n>m .or. iq*n>m) then
        write (1337, 120) ijq, iq, jq, iq*n, jq*n, n, m
        return
      end if
    else
      iq = 1
      jq = 1
    end if

    ! ----------------------------------------------------- copy matrix block

    j0 = n*(jq-1)
    do j = 1, n
      i1 = n*(iq-1) + 1
      jp = j0 + j
      call zcopy(n, a(i1,jp), 1, b(1,j), 1)
    end do

    ! ------------------------------------------------ set up character table

    nc = 0
    do i = 1, 26
      nc = nc + 1
      ipt(nc) = 62 + i
    end do
    do i = 1, 8
      nc = nc + 1
      ipt(nc) = 96 + i
    end do
    do i = 10, 26
      nc = nc + 1
      ipt(nc) = 96 + i
    end do
    do i = 191, 218
      nc = nc + 1
      ipt(nc) = i
    end do
    do i = 35, 38
      nc = nc + 1
      ipt(nc) = i
    end do
    do i = 40, 42
      nc = nc + 1
      ipt(nc) = i
    end do
    do i = 91, 93
      nc = nc + 1
      ipt(nc) = i
    end do

    ! ---------------------------------------------------------------- header
    ic0 = ichar('0')
    n3 = n/100
    n2 = n/10 - n3*10
    n1 = n - n2*10 - n3*100

    if (n<=18) then
      fmt1 = '(8X,I3,''|'','
      fmt2 = '( 9X,''--|'','
      fmt3 = '( 9X,'' #|'','
      fmt4 = '( 9X,''  |'','
    else
      fmt1 = '(   I4,''|'','
      fmt2 = '( 2X,''--|'','
      fmt3 = '( 2X,'' #|'','
      fmt4 = '( 2X,''  |'','
    end if

    lf = 11
    l3 = 11
    if (mcol==0) then
      fmt1 = fmt1(1:lf) // char(ic0+n3) // char(ic0+n2) // char(ic0+n1) // '( 2A1),''|'',I3)'
      fmt2 = fmt2(1:lf) // char(ic0+n3) // char(ic0+n2) // char(ic0+n1) // '(''--''),''|'',I3)'
      fmt3 = fmt3(1:lf) // '60(2X,I2))'
      fmt4 = fmt4(1:lf) // '60(I2,2X))'
      lf = 21
    else
      if (mcol==1) then
        nk = nint(sqrt(real(n, kind=dp)))
      else if (mcol==2) then
        nk = nint(sqrt(real(n/2, kind=dp)))
      else if (mcol==3) then
        nk = 2*nint(sqrt(real(n/2, kind=dp))) - 1
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
        fmt3 = fmt3(1:l3) // fmt3(12:l3)
        fmt4 = fmt4(1:l3) // fmt4(12:l3)
        lf = 2*lf - 11
        l3 = 2*l3 - 11
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
      nsl = nint(sqrt(real(n, kind=dp)))
      do il = 1, nsl
        ilsep(il) = il**2
      end do
    else if (mlin==2) then
      nsl = nint(sqrt(real(n/2, kind=dp)))
      do il = 1, nsl
        ilsep(il) = il**2
      end do
      do il = 1, nsl
        ilsep(nsl+il) = ilsep(nsl) + il**2
      end do
      nsl = 2*nsl
    else if (mlin==3) then
      nsl = 2*nint(sqrt(real(n/2, kind=dp))) - 1
      ilsep(1) = 2
      do k = 2, nsl
        ilsep(k) = ilsep(k-1) + 2*((k+1)/2)
      end do
    end if


    write (nfil, 110) str(1:lstr)
    if (ijq>1000) write (nfil, 130) iq, jq
    write (nfil, fmt3)(i, i=2, n, 2)
    write (nfil, fmt4)(i, i=1, n, 2)
    write (nfil, fmt=fmt2)
    ! ------------------------------------------------------------ header end
    nnon0 = 0
    nd = 0
    ctab(0) = ' '
    dtab(0) = 9999e0_dp

    do i = 1, n
      do j = 1, n
        if (.not. small(b(i,j))) then
          nnon0 = nnon0 + 1
          do id = 1, nd
            if (same(b(i,j),+dtab(id))) then
              iw(j) = +id
              go to 100
            end if
            if (same(b(i,j),-dtab(id))) then
              iw(j) = -id
              go to 100
            end if
          end do
          ! ----------------------------------------------------------- new
          ! element
          nd = nd + 1
          iw(j) = nd
          dtab(nd) = b(i, j)
          if (abs(dtab(nd)-1.0e0_dp)*tol<1.0e0_dp) then
            ctab(nd) = '1'
          else if (abs(dtab(nd)+1.0e0_dp)*tol<1.0e0_dp) then
            dtab(nd) = +1.0e0_dp
            ctab(nd) = '1'
            iw(j) = -nd
          else if (abs(dtab(nd)-cone)*tol<1.0e0_dp) then
            ctab(nd) = 'i'
          else if (abs(dtab(nd)+cone)*tol<1.0e0_dp) then
            dtab(nd) = +cone
            ctab(nd) = 'i'
            iw(j) = -nd
          else
            ctab(nd) = char(ipt(1+mod((nd+1),nc)))
          end if
        else
          iw(j) = 0
        end if
100   end do
      ! ------------------------------------------------------------ write line
      write (nfil, fmt=fmt1) i, (vz(isign(1,iw(j))), ctab(abs(iw(j))), j=1, n), i

      do isl = 1, nsl
        if (i==ilsep(isl)) write (nfil, fmt=fmt2)
      end do
    end do

    ! ------------------------------------------------------------------ foot

    write (nfil, fmt4)(i, i=1, n, 2)
    write (nfil, fmt3)(i, i=2, n, 2)

    if (k_fmt_fil>0) then
      write (nfil, 140)(id, ctab(id), dtab(id), id=1, nd)
      write (nfil, 150) nnon0, tolp, n*n - nnon0, tolp
    else
      write (nfil, *) ' '
    end if

110 format (/, 8x, a, /)
120 format (/, 1x, 79('*'), /, 10x, 'inconsistent call of <CMATSTR>', /, 10x, 'argument IJQ =', i8, '  implies IQ=', i3, '   JQ=', i3, /, 10x, 'IQ*N=', i6, ' > M   or   JQ*N=', i6, &
      ' > M   for N =', i4, ' M=', i4, /, 1x, 79('*'), /)
130 format (8x, 'IQ-JQ-block  for  IQ = ', i3, '   JQ = ', i3, /)
140 format (/, 8x, 'symbols used:', /, (8x,i3,3x,a1,2x,2f20.12))
150 format (/, 8x, i5, ' elements   >', 1p, e9.1, /, 8x, i5, ' elements   <', 1p, e9.1, /)
  end subroutine cmatstr

end module mod_cmatstr
