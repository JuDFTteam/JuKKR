!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Write density of states to file
!> Author: People who wrote it
!> Write density of states to file. Both complex DOS and the real part of the DOS
!> are printed to file in an l-decomposed fashion.
!------------------------------------------------------------------------------------
module mod_wrldos
  use :: mod_datatypes, only: dp
  private :: dp

contains
 
  !-------------------------------------------------------------------------------
  !> Summary: Write density of states to file
  !> Author: 
  !> Category: physical-observables, KKRhost 
  !> Deprecated: False 
  !> Write density of states to file. Both complex DOS and the real part of the DOS
  !> are printed to file in an l-decomposed fashion.
  !-------------------------------------------------------------------------------
  !> @note Jonathan Chico: It might be a good idea to improve the headers and to set
  !> printing flags so that the types of files that one wants to plot are actually plotted,
  !> e.g. the complex DOS is printed if one passes a certain flag in the `inputcard`
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine wrldos(den,ez,wez,lmaxd1,iemxd,npotd,ititle,efermi,e1,e2,alatc,tk,     &
    nacls1,nspinpot,natyp,conc,ielast,intervx,intervy,intervz,dostot)
    use :: mod_version_info
    use :: mod_datatypes
    use :: mod_constants, only: pi, kb, ryd, czero
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer, intent(in) :: natyp    !! Number of kinds of atoms in unit cell
    integer, intent(in) :: npotd    !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent(in) :: iemxd    !! Dimension for energy-dependent arrays
    integer, intent(in) :: lmaxd1   !! lmax+1
    integer, intent(in) :: nacls1
    integer, intent(in) :: ielast
    integer, intent(in) :: intervx  !! Number of intervals in x-direction for k-net in IB of the BZ
    integer, intent(in) :: intervy  !! Number of intervals in y-direction for k-net in IB of the BZ
    integer, intent(in) :: intervz  !! Number of intervals in z-direction for k-net in IB of the BZ
    integer, intent(in) :: nspinpot !! krel*2 + (1-krel)*nspin
    real (kind=dp), intent(in) :: e1      !! Lower value (in Ryd) for the energy contour
    real (kind=dp), intent(in) :: e2      !! Maximum value (in Ryd) for the DOS calculation Controls also [NPT2] in some cases
    real (kind=dp), intent(in) :: tk      !! Temperature
    real (kind=dp), intent(in) :: alatc   !! Lattice constant in a.u.
    real (kind=dp), intent(in) :: efermi  !! Fermi energy
    integer, dimension(20, npotd), intent(in) :: ititle
    real (kind=dp), dimension(*), intent(in)  :: conc !! Concentration of a given atom
    complex (kind=dp), dimension(iemxd), intent(in) :: ez
    complex (kind=dp), dimension(iemxd), intent(in) :: wez
    complex (kind=dp), dimension(0:lmaxd1, ielast, npotd), intent(in) :: den
    ! .. In/Out variables
    real (kind=dp), dimension(0:lmaxd1, 2), intent(inout) :: dostot
    ! ..
    ! .. Local Scalars ..
    integer :: i1, ia, ie, ipot, ispin, l
    real (kind=dp) :: dos, dossgn, efctor
    complex (kind=dp) :: doscmplx
    character (len=8) :: dosfl0
    character (len=11) :: dosfl
    ! ..
    ! .. External Functions ..
    logical :: test
    external :: test
    ! ..
    dosfl0 = 'dos.atom'
    efctor = 1.0e0_dp
    do ispin = 1, nspinpot
      do l = 0, lmaxd1
        dostot(l, ispin) = 0.0e0_dp
      end do
    end do
    do i1 = 1, natyp
      if (i1<10) write (dosfl, fmt='(A8,I1)') dosfl0, i1
      if (i1>=10 .and. i1<100) write (dosfl, fmt='(A8,I2)') dosfl0, i1
      if (i1>=100) write (dosfl, fmt='(A8,I3)') dosfl0, i1
      open (48, file=trim(dosfl), form='formatted')
      call version_print_header(48)
      do ispin = 1, nspinpot
        ipot = nspinpot*(i1-1) + ispin
        dossgn = 1.0e0_dp
        if (ispin/=nspinpot) dossgn = -1.0e0_dp

        write (48, fmt=110)(ititle(ia,ipot), ia=1, 19)
        write (48, fmt=120) i1
        write (48, fmt=130) ispin, ielast, e1, e2, efermi, efctor
        write (48, fmt=140) efermi
        write (48, fmt=150) tk, pi*kb*tk, alatc, intervx, intervy, intervz, nacls1
        do ie = 1, ielast
          dos = 0.0e0_dp
          do l = 0, lmaxd1
            dos = dos - 2.0e0_dp*aimag(den(l,ie,ipot))/pi/real(nspinpot, kind=dp)
            dostot(l, ispin) = dostot(l, ispin) + aimag(wez(ie)*den(l,ie,ipot))
          end do
          write (48, fmt=160) real(ez(ie))*efctor, dos*dossgn/efctor, (-2.0e0_dp*aimag(den(l,ie,ipot))*dossgn/efctor/pi/real(nspinpot,kind=dp), l=0, lmaxd1)
        end do
        write (48, fmt=180)(dostot(l,ispin)/efctor/real(nspinpot,kind=dp), l=0, lmaxd1)
        if (ispin/=nspinpot) write (48, fmt=100)
      end do
      close (48)
    end do

    !--------------------------------------------------------------------------------
    ! Write complex DOS in unit 49:
    !--------------------------------------------------------------------------------
    open (49, file='complex.dos', form='formatted')
    call version_print_header(49)
    write (49, *) natyp*nspinpot
    write (49, *) ielast
    write (49, *) lmaxd1
    do i1 = 1, natyp

      do ispin = 1, nspinpot
        ipot = nspinpot*(i1-1) + ispin
        dossgn = 1.0e0_dp
        if (ispin/=nspinpot) dossgn = -1.0e0_dp

        write (49, fmt=110)(ititle(ia,ipot), ia=1, 19)
        write (49, fmt=120) i1
        write (49, fmt=130) ispin, ielast, e1, e2, efermi, efctor
        write (49, fmt=140) efermi
        write (49, fmt=150) tk, pi*kb*tk, alatc, intervx, intervy, intervz, nacls1
        do ie = 1, ielast
          doscmplx = czero
          do l = 0, lmaxd1
            doscmplx = doscmplx - 2.0e0_dp*den(l, ie, ipot)/pi/real(nspinpot, kind=dp)
          end do
          write (49, fmt=170) ez(ie)*efctor, (-2.0e0_dp*den(l,ie,ipot)*dossgn/efctor/pi/real(nspinpot,kind=dp), l=0, lmaxd1), doscmplx*dossgn/efctor
        end do
        if (ispin/=nspinpot .or. i1/=natyp) write (49, fmt=100)
      end do
    end do
    close (49)

    !--------------------------------------------------------------------------------
    ! Write total DOS summed over atoms and spins(complex)
    !--------------------------------------------------------------------------------
    open (49, file='total_cmplx.dos', form='formatted')
    call version_print_header(49)
    write (49, fmt='(4A16)') '# Real(E)', '  Im(E)', ' Re(DEN)', ' Im(DEN)'
    do ie = 1, ielast
      doscmplx = czero 
      do i1 = 1, natyp
        do ispin = 1, nspinpot
          ipot = nspinpot*(i1-1) + ispin
          do l = 0, lmaxd1
            doscmplx = doscmplx - conc(i1)*2.0e0_dp*den(l, ie, ipot)/pi/real(nspinpot, kind=dp)
          end do
        end do
      end do
      write (49, fmt='(10E16.8)') ez(ie), doscmplx
    end do
    close (49)


    return
    ! ccc 9000 FORMAT ('&')
100 format (' ')
110 format ('#', 19a4)
120 format ('# I1    :', i8)
130 format ('# ISPIN :', i8, '   IELAST :', i5, /, '# E1,E2 :', 2f12.5, ' EFERMI :', f12.5, '   EFCTR', f10.6)
140 format ('# FERMI :', f12.5)
150 format ('# TK    =', f8.1, '   Kelvin =', 3p, f8.3, ' mRyd', 0p, /, '# ALAT   :', f12.5, /, '# INTERV X,Y,Z  :', 3i5, /, '# NACLS :', i8)
160 format (1p, 8e15.7)
170 format (16('(',e16.8,',',e16.8,')'))
180 format ('# Integrated DOS ', 1p, d10.3, 7d11.3)
190 format ('&')
  end subroutine wrldos

end module mod_wrldos
