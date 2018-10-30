!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Write charges and magnetic and orbital moments to file
!> Author: 
!> Write charges and magnetic and orbital moments to file. The output is l-decomposed
!------------------------------------------------------------------------------------
module mod_wrmoms
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Write charges and magnetic and orbital moments to file
  !> Author: 
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Write charges and magnetic and orbital moments to file. The output is l-decomposed
  !-------------------------------------------------------------------------------
  subroutine wrmoms(krel,natyp,nspin,texts,textl,textns,charge,muorb,lmaxd,lmaxd1)

    implicit none

    ! Dummy arguments
    integer, intent(in) :: krel !! Switch for non- (or scalar-) relativistic/relativistic (Dirac) program (0/1). Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent(in) :: lmaxd  !! Maximum l component in wave function expansion
    integer, intent(in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent(in) :: lmaxd1 !! lmax+1
    character (len=5), intent(in) :: textns
    real (kind=dp), dimension(0:lmaxd1, natyp, 2), intent(in) :: charge
    character (len=4), dimension(0:6), intent(in) :: textl
    character (len=7), dimension(3), intent(in) :: texts

    ! Local variables
    real (kind=dp) :: chval
    real (kind=dp), dimension(natyp) :: mutot
    real (kind=dp), dimension(natyp) :: chtot
    real (kind=dp), dimension(natyp, 2) :: sumch
    real (kind=dp), dimension(natyp, 0:lmaxd1+1) :: muspin
    real (kind=dp), dimension(0:lmaxd1+1, 3, natyp) :: muorb
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
