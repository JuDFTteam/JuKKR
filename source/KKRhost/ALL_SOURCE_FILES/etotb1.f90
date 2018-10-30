!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_etotb1

contains

  !-------------------------------------------------------------------------------
  !> Summary: Collects total energy of cluster
  !> Author: B. Drittler, P. Zahn, V. Popescu
  !> Category: KKRhost, total-energy
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculate the total energy of the cluster.
  !> gather all energy-parts which are calculated in different
  !> subroutines .
  !> since the program uses group theory only shell-indices
  !> are used instead of atom-indices .
  !>
  !> B. Drittler   May 1987
  !>
  !> modified for supercells with nshell(i) atoms of type i in the
  !> unit cell
  !> P. Zahn       Oct. 95
  !>
  !> adopted for more atoms per site (CPA) V. Popescu Feb. 02
  !-------------------------------------------------------------------------------
  subroutine etotb1(ecou, epotin, espc, espv, exc, kpre, lmax, lpot, lcoremax, nspin, natyp, nshell, conc, idoldau, lopt, eu, edcldau)

    use :: mod_types, only: t_inc
    use :: mod_datatypes, only: dp
    use :: global_variables, only: natypd, npotd, lmaxd
    implicit none

    integer :: kpre, lmax, lpot, natyp, nspin, idoldau
    real (kind=dp) :: conc(natypd)
    real (kind=dp) :: ecou(0:lpot, *), epotin(*), espc(0:3, npotd), espv(0:(lmaxd+1), npotd), exc(0:lpot, *), eu(*), edcldau(*)
    integer :: lcoremax(*), nshell(*), lopt(*)
    real (kind=dp) :: bandesum, bandet, ecous, edc, efctor, et, etot, excs
    real (kind=dp) :: etotldau
    real (kind=dp) :: dble
    integer :: iatyp, ipot, is, ispin, l
    character (len=4) :: textl(0:6)
    character (len=5) :: textns
    character (len=13) :: texts(3)

    logical :: test
    external :: test

    data textl/' s =', ' p =', ' d =', ' f =', ' g =', ' h =', ' i ='/
    data texts/' spin down   ', ' spin  up    ', ' paramagnetic'/
    data textns/' ns ='/
    ! ---> loop over host atoms
    efctor = 1.0d0/13.6058d0

    etot = 0.0d0
    bandesum = 0.0d0
    etotldau = 0.0d0

    if (kpre==1 .and. (t_inc%i_write>0)) write (1337, fmt=100)



    do iatyp = 1, natyp

      if (kpre==1 .and. (t_inc%i_write>0)) write (1337, fmt=110) iatyp

      edc = 0.0d0
      et = 0.0d0
      bandet = 0.0d0
      ecous = 0.0d0
      excs = 0.0d0

      is = 0
      if (nspin==1) is = is + 2
      do ispin = 1, nspin
        is = is + 1
        ipot = (iatyp-1)*nspin + ispin
        ! -> LDA+U
        if (kpre==1 .and. (t_inc%i_write>0)) then
          write (1337, fmt=120) texts(is)
          write (1337, fmt=130)(textl(l), espc(l,ipot), l=0, lcoremax(iatyp))
          write (1337, fmt=140)(textl(l), espv(l,ipot), l=0, lmax)
          write (1337, fmt=150) textns, espv((lmaxd+1), ipot)
        end if

        do l = 0, lcoremax(iatyp)
          et = et + espc(l, ipot)
        end do

        do l = 0, lmax
          bandet = bandet + espv(l, ipot)
          et = et + espv(l, ipot)
        end do
        bandet = bandet + espv((lmaxd+1), ipot)
        et = et + espv((lmaxd+1), ipot)
      end do
      ! --->  sum up Coulomb and Ex.-Corel. contribution


      et = et + eu(iatyp)
      bandet = bandet + eu(iatyp)
      if (kpre==1 .and. idoldau==1 .and. lopt(iatyp)>=0 .and. (t_inc%i_write>0)) write (1337, 280) eu(iatyp)



      do l = 0, lpot
        ecous = ecous + ecou(l, iatyp)
        excs = excs + exc(l, iatyp)
      end do

      if (kpre==1 .and. (t_inc%i_write>0)) then
        write (1337, fmt=160) et
        write (1337, fmt=170) bandet
        write (1337, fmt=180)(l, ecou(l,iatyp), l=0, lpot)
        write (1337, fmt=190)
        write (1337, fmt=270) ecous
        write (1337, fmt=200)(l, exc(l,iatyp), l=0, lpot)
        write (1337, fmt=190)
        write (1337, fmt=260) excs
        write (1337, fmt=240) epotin(iatyp)
      end if

      if (.not. (test('NoMadel '))) then

        et = et + ecous + excs
        edc = edc + ecous + excs

        et = et + epotin(iatyp) - edcldau(iatyp)
        edc = edc + epotin(iatyp) - edcldau(iatyp)

        if (kpre==1 .and. (t_inc%i_write>0)) then
          if (idoldau==1 .and. lopt(iatyp)>=0) write (1337, 290) - edcldau(iatyp)
          write (1337, fmt=250) edc
        end if
        ! IATYP = 1,NATYP
      end if

      if (natyp>1 .or. nshell(iatyp)>1) then
        if (t_inc%i_write>0) write (1337, fmt=210) iatyp, et
        if (kpre==1 .and. idoldau==1 .and. lopt(iatyp)>=0 .and. (t_inc%i_write>0)) write (1337, 300) eu(iatyp) - edcldau(iatyp)
        write (1337, fmt=310)
      end if

      etot = etot + et*dble(nshell(iatyp))*conc(iatyp)
      bandesum = bandesum + bandet*dble(nshell(iatyp))*conc(iatyp)

    end do

    if (t_inc%i_write>0) write (1337, fmt=220) bandesum
    if (t_inc%i_write>0) write (1337, fmt=230) etot, etot/efctor
    write (*, fmt=320) etot


    return
    ! 17.10.95 ***************************************************************
100 format (32('='), ' TOTAL ENERGIES ', 31('='), /)
110 format (3x, 'Total energies atom ', i3, /, 3x, 23('-'))
120 format (5x, 'single particle energies ', a13)
130 format (7x, '  core   contribution : ', 2(a4,f15.8), /, (31x,2(a4,f15.8)))
140 format (7x, 'valence  contribution : ', 2(a4,f15.8), /, (31x,2(a4,f15.8)))
150 format (7x, '                        ', a4, f15.8)
160 format (5x, 68('-'), /, 5x, 'total contribution of the single particle energies :', 1x, f15.8)
170 format (5x, '                              band energy per atom :', 1x, f15.10, /)
180 format (5x, 'coulomb  contribution : ', 2(i3,1x,f15.8), /, (29x,2(i3,1x,f15.8)))
190 format (5x, 68('-'))
    ! ************************************************************************
200 format (5x, 'ex.-cor. contribution : ', 2(i3,1x,f15.8), /, (29x,2(i3,1x,f15.8)))
210 format (/, 3x, 'Total contribution of atom', i3, ' =', f15.8)
220 format (5x, '                              sum of band energies :', 1x, f15.10, /, 3x, 70('-'))
230 format (/, 3x, 70('+'), /, 15x, 'TOTAL ENERGY in ryd. : ', f25.8, /, 15x, '                 eV  : ', f25.8, /, 3x, 70('+'))
240 format (5x, 'eff. pot. contribution     : ', f15.8)
250 format (5x, 'total double counting contribution                 :', 1x, f15.8)
260 format (5x, 'tot. ex.-cor. contribution : ', f15.8, /)
270 format (5x, 'tot. coulomb contribution : ', f15.8, /)
280 format (/, 5x, 'LDA+U correction to the single particle energy     :', f16.8)
290 format (/, 5x, 'LDA+U double counting contribution                 :', f16.8)
    ! calculate the total energy of the cluster .
300 format (3x, '   including LDA+U correction :', f15.8)
310 format (3x, 70('-'))
    ! gather all energy-parts which are calculated in different
320 format ('TOTAL ENERGY in ryd. : ', f25.8, /, 15x)
  end subroutine etotb1

end module mod_etotb1
