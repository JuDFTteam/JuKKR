module mod_wrmoms
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine wrmoms(krel, natyp, nspin, texts, textl, textns, charge, muorb, lmaxd, lmaxd1)

    implicit none

    ! Dummy arguments
    integer :: krel, lmaxd, lmaxd1, natyp, nspin
    character (len=5) :: textns
    real (kind=dp) :: charge(0:lmaxd1, natyp, 2)
    character (len=4) :: textl(0:6)
    character (len=7) :: texts(3)

    ! Local variables
    real (kind=dp) :: chtot(natyp), chval
    real (kind=dp) :: muorb(0:lmaxd1+1, 3, natyp)
    real (kind=dp) :: muspin(natyp, 0:lmaxd1+1), mutot(natyp), sumch(natyp, 2)
    character (len=80) :: fmt1, fmt2, fmt31, fmt32
    integer :: is, ispin, it, l, lf1, lf2

    write (1337, *)

    if ((krel==1) .or. (nspin==2)) then
      write (1337, '(78("#"))')
      write (1337, 100)
      write (1337, '(78("#"))')
    else
      write (1337, '(44("#"))')
      write (1337, 110)
      write (1337, '(44("#"))')
    end if

    write (1337, *)
    write (1337, 120, advance='no')
    write (6, 130, advance='no')
    do it = 1, natyp
      muspin(it, lmaxd1+1) = 0e0_dp
      sumch(it, 1) = 0e0_dp
      sumch(it, 2) = 0e0_dp
      do l = 0, lmaxd1
        do ispin = 1, nspin
          sumch(it, ispin) = sumch(it, ispin) + charge(l, it, ispin)
        end do
        muspin(it, l) = charge(l, it, 2) - charge(l, it, 1)
        muspin(it, lmaxd1+1) = muspin(it, lmaxd1+1) + muspin(it, l)
      end do
      chtot(it) = sumch(it, 1) + sumch(it, 2)
    end do

    if (krel==1) then
      do it = 1, natyp
        mutot(it) = muspin(it, lmaxd1+1) + muorb(lmaxd1+1, 3, it)
      end do
    end if

    is = 0
    if (nspin==1) is = is + 2
    do ispin = 1, nspin
      is = is + 1
      write (1337, 140, advance='no') texts(is)
      write (6, 140, advance='no') texts(is)
    end do

    if (krel==1) then
      write (1337, 150, advance='no')
      write (6, 150, advance='no')
      write (1337, 160)
      write (6, 160)
    else
      if (nspin==2) write (1337, 150, advance='no')
      if (nspin==2) write (6, 150, advance='no')
      write (1337, *)
      write (6, *)
    end if

    write (1337, '(3X,26("="))', advance='no')
    if (krel==1) then
      write (1337, '(46("="))')
    else
      if (nspin==2) write (1337, '(23("="))', advance='no')
      write (1337, *)
    end if

    fmt1 = '(4X,I3,2X,A4,2(F12.8),2X,F8.4'
    fmt2 = '(9X,A4,2(F12.8),2X,F8.4'
    fmt31 = '(4X,I3,2X,A4,F12.8)'
    fmt32 = '(9X,A4,F12.8)'
    lf1 = 30
    lf2 = 24

    if (krel==1) then
      fmt1 = fmt1(1:lf1) // ',2X,3F8.4)'
      fmt2 = fmt2(1:lf2) // ',2X,3F8.4)'
    else
      if (nspin==2) then
        fmt1 = fmt1(1:lf1) // ')'
        fmt2 = fmt2(1:lf2) // ')'
      else
        fmt1 = fmt31
        fmt2 = fmt32
      end if
    end if

    do it = 1, natyp
      if (krel==1) then
        write (1337, fmt=fmt1) it, textl(0), (charge(0,it,ispin), ispin=1, nspin), muspin(it, 0), muorb(0, 3, it), (muorb(0,ispin,it), ispin=1, nspin)
      else
        if (nspin==2) then
          write (1337, fmt=fmt1) it, textl(0), (charge(0,it,ispin), ispin=1, nspin), muspin(it, 0)
        else
          write (1337, fmt=fmt1) it, textl(0), charge(0, it, 1)
        end if
      end if

      do l = 1, lmaxd
        if (krel==1) then
          write (1337, fmt=fmt2) textl(l), (charge(l,it,ispin), ispin=1, nspin), muspin(it, l), muorb(l, 3, it), (muorb(l,ispin,it), ispin=1, nspin)
        else
          if (nspin==2) then
            write (1337, fmt=fmt2) textl(l), (charge(l,it,ispin), ispin=1, nspin), muspin(it, l)
          else
            write (1337, fmt=fmt2) textl(l), charge(l, it, 1)
          end if
        end if

      end do

      if (krel==1) then
        write (1337, fmt=fmt2) textns, (charge(lmaxd1,it,ispin), ispin=1, nspin), muspin(it, lmaxd1), muorb(lmaxd1, 3, it), (muorb(lmaxd1,ispin,it), ispin=1, nspin)
      else
        if (nspin==2) then
          write (1337, fmt=fmt2) textns, (charge(lmaxd1,it,ispin), ispin=1, nspin), muspin(it, lmaxd1)
        else
          write (1337, fmt=fmt2) textns, charge(lmaxd1, it, 1)
        end if
      end if

      write (1337, '(10x,19("-"))', advance='no')
      if (krel==1) then
        write (1337, '(44("-"))')
        write (1337, fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), muspin(it, lmaxd1+1), muorb(lmaxd1+1, 3, it), (muorb(lmaxd1+1,ispin,it), ispin=1, nspin)
        write (6, fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), muspin(it, lmaxd1+1), muorb(lmaxd1+1, 3, it), (muorb(lmaxd1+1,ispin,it), ispin=1, nspin)
        write (1337, '(25X,F12.8,12X,F8.4)') chtot(it), mutot(it)
      else
        if (nspin==2) then
          write (1337, '(17("-"))')
          write (1337, fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), muspin(it, lmaxd1+1)
          write (6, fmt=fmt2) ' TOT', (sumch(it,ispin), ispin=1, nspin), muspin(it, lmaxd1+1)
          write (1337, '(25X,F12.8)') chtot(it)
        else
          write (1337, *)
          write (1337, fmt=fmt2) ' TOT', sumch(it, 1)
          write (6, fmt=fmt2) ' TOT', sumch(it, 1)
        end if
      end if

      if (it/=natyp) then
        write (1337, '(3X,26("="))', advance='no')
        if (krel==1) then
          write (1337, '(40("="))')
        else
          if (nspin==2) write (1337, '(17("="))', advance='no')
          write (1337, *)
        end if
      end if
    end do

    write (1337, *)
    if ((krel==1) .or. (nspin==2)) then
      write (1337, '(78("#"))')
    else
      write (1337, '(44("#"))')
    end if
    write (1337, *)

    chval = 0.e0_dp
    do it = 1, natyp
      chval = chval + chtot(it)
    end do
    write (1337, *) 'Sum of valence charges of atoms (local summation)', chval

    return

100 format (15x, 'l-decomposed valence charges and magnetic moments')
110 format (8x, 'l-decomposed valence charges')
120 format (3x, 'ATOM      ')
130 format (3x, '          ')
140 format (2x, 'Ne ', a7)
150 format ('    m_spin')
160 format ('    m_orb   spin dn  spin up')
  end subroutine wrmoms

end module mod_wrmoms
