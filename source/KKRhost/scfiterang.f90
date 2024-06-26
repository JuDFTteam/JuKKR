!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Self-consistent mixing of the spin angle specifying the **local** frame of reference
!> Author: 
!> Self-consistent mixing of the spin angle specifying the **local** frame of reference
!------------------------------------------------------------------------------------
module mod_scfiterang

contains

  !-------------------------------------------------------------------------------
  !> Summary: Self-consistent mixing of the spin angle specifying the **local** frame of reference
  !> Author:
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Self-consistent mixing of the spin angle specifying the **local** frame of reference.
  !> As the self-consistent process evolves the local spin directions are varied and a
  !> new local frame of reference is found.
  !-------------------------------------------------------------------------------
  subroutine scfiterang(itrscf,itoq,fact,mvphi,mvtet,mvgam,qmphi,qmtet,qmgam,nq,nk, &
    erravang,nqmax,ntmax,nmvecmax,nkmmax)

    use :: mod_wunfiles, only: t_params
    use :: mod_datatypes, only: dp
    use :: mod_calcrotmat
    implicit none

    ! PARAMETER definitions

    integer :: ixtrmax
    parameter (ixtrmax=4)
    real (kind=dp) :: dangmax
    parameter (dangmax=3.0_dp)

    ! Dummy arguments

    real (kind=dp) :: erravang
    integer :: itrscf, nk, nkmmax, nmvecmax, nq, nqmax, ntmax
    real (kind=dp) :: fact(0:100), mvgam(ntmax, nmvecmax), mvphi(ntmax, nmvecmax), mvtet(ntmax, nmvecmax), qmgam(nqmax), qmphi(nqmax), qmtet(nqmax)
    integer :: itoq(ntmax, nqmax)

    ! Local variables

    real (kind=dp) :: a, b, c, phixtr, tetxtr, d12, d23, d3x
    real (kind=dp) :: delphi, deltet, lasterr, mixing, qmphimix, qmtetmix, wn, wo
    integer :: i, imv, iprev, iprint, iq, it, itab, ixtr
    real (kind=dp) :: qmgamtab(nqmax, 3), qmphitab(nqmax, 3), qmtettab(nqmax, 3)
    complex (kind=dp) :: drotq(nkmmax, nkmmax, nqmax)

    mixing = 1.0_dp
    wo = 1.0_dp - mixing
    wn = mixing
    iprint = 0
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! ITERMDIR I/O

    qmtet = t_params%qmtet
    qmphi = t_params%qmphi
    qmphitab = t_params%qmphitab
    qmtettab = t_params%qmtettab
    qmgamtab = t_params%qmgamtab
    itab = t_params%itab
    lasterr = t_params%lasterr
    ! READ (67) QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ITAB,LASTERR
    ! REWIND (67)

    ! ITERMDIR I/O
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    write (1337, fmt=110)

    if (itrscf==0) then

      ! ----------------------------------------------------- store old angles
      ! --------------------------------------- dummy for tbkkr, done in main0
      do iq = 1, nq
        qmphitab(iq, 1) = qmphi(iq)
        qmtettab(iq, 1) = qmtet(iq)
        qmgamtab(iq, 1) = qmgam(iq)
        do i = 2, 3
          qmphitab(iq, i) = 0.0_dp
          qmtettab(iq, i) = 0.0_dp
          qmgamtab(iq, i) = 0.0_dp
        end do
      end do

      return

    else

      ! ------------------------- set new local frame of reference according to
      ! -------------------------------------------- orientation of spin moment

      imv = 1

      do iq = 1, nq

        it = itoq(1, iq)

        qmphi(iq) = mvphi(it, imv)
        qmtet(iq) = mvtet(it, imv)
        qmgam(iq) = mvgam(it, imv)

      end do

      ! =======================================================================

      if (itrscf<=2) then
        itab = itrscf
        do iq = 1, nq
          qmphi(iq) = wo*qmphitab(iq, itab) + wn*qmphi(iq)
          qmtet(iq) = wo*qmtettab(iq, itab) + wn*qmtet(iq)
          qmgam(iq) = wo*qmgamtab(iq, itab) + wn*qmgam(iq)
        end do

        itab = itrscf + 1
        do iq = 1, nq

          qmphitab(iq, itab) = qmphi(iq)
          qmtettab(iq, itab) = qmtet(iq)
          qmgamtab(iq, itab) = qmgam(iq)
        end do

        iprev = itab - 1

      else
        if (lasterr>2.0_dp) then
          ixtr = 2
        else
          ixtr = ixtrmax
        end if

        do iq = 1, nq

          qmphimix = wo*qmphitab(iq, 3) + wn*qmphi(iq)
          qmtetmix = wo*qmtettab(iq, 3) + wn*qmtet(iq)

          a = qmphitab(iq, 2)
          b = (qmphitab(iq,3)-qmphitab(iq,1))*0.5_dp
          c = (qmphitab(iq,3)+qmphitab(iq,1)-2.0_dp*a)*0.5_dp
          phixtr = a + b*real(ixtr, kind=dp) + c*real(ixtr, kind=dp)**2
          if (phixtr>=0.0_dp .and. phixtr<=360.0_dp .and. abs(phixtr-qmphitab(iq,3))<dangmax) then
            qmphi(iq) = phixtr
          else
            qmphi(iq) = qmphimix
          end if

          if (qmphi(iq)<0.0_dp) qmphi(iq) = qmphi(iq) + 360.0_dp
          if (qmphi(iq)>360.0_dp) qmphi(iq) = qmphi(iq) - 360.0_dp

          a = qmtettab(iq, 2)
          b = (qmtettab(iq,3)-qmtettab(iq,1))*0.5_dp
          c = (qmtettab(iq,3)+qmtettab(iq,1)-2.0_dp*a)*0.5_dp
          tetxtr = a + b*real(ixtr, kind=dp) + c*real(ixtr, kind=dp)**2
          d12 = qmtettab(iq, 1) - qmtettab(iq, 2)
          d23 = qmtettab(iq, 2) - qmtettab(iq, 3)
          d3x = qmtettab(iq, 3) - tetxtr
          if (tetxtr>=0.0_dp .and. tetxtr<=180.0_dp .and. abs(tetxtr-qmtettab(iq,3))<dangmax .and. (d12*d23>0.0_dp) .and. (d23*d3x>0.0_dp)) then
            qmtet(iq) = tetxtr
          else
            qmtet(iq) = qmtetmix
          end if

          a = qmgamtab(iq, 2)
          b = (qmgamtab(iq,3)-qmgamtab(iq,1))*0.5_dp
          c = (qmgamtab(iq,3)+qmgamtab(iq,1)-2.0_dp*a)*0.5_dp
          qmgam(iq) = a + b*real(ixtr, kind=dp) + c*real(ixtr, kind=dp)**2

          if (iprint>0) write (1337, 100) iq, ixtr, lasterr, ('PHI', i, qmphitab(iq,i), i=1, 3), 'PHIMIX', qmphimix, 'PHIXTR', phixtr, 'PHINEW', qmphi(iq), &
            ('TET', i, qmtettab(iq,i), i=1, 3), 'TETMIX', qmtetmix, 'TETXTR', tetxtr, 'TETNEW', qmtet(iq)

          do i = 1, 2
            qmphitab(iq, i) = qmphitab(iq, i+1)
            qmtettab(iq, i) = qmtettab(iq, i+1)
            qmgamtab(iq, i) = qmgamtab(iq, i+1)
          end do
          qmphitab(iq, 3) = qmphi(iq)
          qmtettab(iq, 3) = qmtet(iq)
          qmgamtab(iq, 3) = qmgam(iq)

        end do

        iprev = 2
      end if

      write (1337, 130)

      erravang = 0.0_dp
      do iq = 1, nq

        delphi = abs(qmphi(iq)-qmphitab(iq,iprev))
        deltet = abs(qmtet(iq)-qmtettab(iq,iprev))
        erravang = max(erravang, delphi, deltet)

        write (1337, 140) iq, qmphitab(iq, iprev), qmtettab(iq, iprev), qmphi(iq), qmtet(iq), delphi, deltet

        ! --> update the rotation matrices DROTQ for the new angles

        call calcrotmat(nk, 3, qmphi(iq), qmtet(iq), 0.0_dp, drotq(1,1,iq), fact, nkmmax)

      end do

      write (1337, fmt=120) itrscf, erravang
      write (1337, '(I5,4F10.3,'' #  ANGLES'',/,79("+"),/)') itrscf, (qmphi(iq), qmtet(iq), iq=1, min(2,nq))

      lasterr = erravang
      ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ! ITERMDIR I/O

      ! WRITE (67) QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ITAB,LASTERR
      ! WRITE (67) DROTQ
      t_params%qmtet = qmtet
      t_params%qmphi = qmphi
      t_params%qmphitab = qmphitab
      t_params%qmtettab = qmtettab
      t_params%qmgamtab = qmgamtab
      t_params%itab = itab
      t_params%lasterr = lasterr
      t_params%drotq = drotq

      ! ITERMDIR I/O
      ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    end if
100 format (/, 5x, 'IQ', i3, '  IXTR', i3, '  LAST ERROR', f12.4, /, 3(3(5x,a,5x,'TAB',i2,2x,f12.6,/),3(5x,a,9x,f12.6,/)))
    ! ====================================================================
110 format (79('+'), /, 28x, 'ANGLE - mixing scheme')
120 format (5x, 'iter.', i4, '     max. CHANGE = ', f12.6)
130 format (/, 5x, ' setting new LOCAL frame of reference ', /, /, 5x, 'IQ  old   phi      tet    --> new   phi      tet', '     del   phi      tet')
140 format (5x, i2, 4x, 2f9.4, 8x, 2f9.4, 5x, 2f9.4)
  end subroutine scfiterang

end module mod_scfiterang
