module mod_wrldos
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine wrldos(den, ez, wez, lmaxd1, iemxd, npotd, ititle, efermi, e1, e2, alatc, tk, nacls1, nspinpot, natyp, conc, ielast, intervx, intervy, intervz, dostot)
    use :: mod_version_info
    use :: mod_datatypes
    implicit none
    ! .. Parameters ..
    real (kind=dp) :: kb
    parameter (kb=0.6333659e-5_dp)
    ! ..
    ! .. Scalar Arguments ..
    real (kind=dp) :: alatc, e1, e2, efermi, tk
    integer :: ielast, iemxd, intervx, intervy, intervz, lmaxd1, nacls1, natyp, npotd
    ! ===== uses spin-up and down also in the REL mode (KREL=1)
    integer :: nspinpot
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: den(0:lmaxd1, ielast, npotd), ez(iemxd), wez(iemxd)
    real (kind=dp) :: dostot(0:lmaxd1, 2), conc(*) ! CONC(NATYPD)
    integer :: ititle(20, npotd)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: doscmplx
    real (kind=dp) :: dos, dossgn, efctor, pi
    integer :: i1, ia, ie, ipot, ispin, l
    character (len=8) :: dosfl0
    character (len=11) :: dosfl
    ! ..
    ! .. External Functions ..
    logical :: test
    external :: test
    ! ..
    pi = 4.0e0_dp*atan(1.0e0_dp)
    dosfl0 = 'dos.atom'
    efctor = 1.0e0_dp
    if (test('EV      ')) efctor = 13.6058e0_dp
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

    ! Write complex DOS in unit 49:
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
          doscmplx = cmplx(0.0e0_dp, 0.e0_dp, kind=dp)
          do l = 0, lmaxd1
            doscmplx = doscmplx - 2.0e0_dp*den(l, ie, ipot)/pi/real(nspinpot, kind=dp)
          end do
          write (49, fmt=170) ez(ie)*efctor, (-2.0e0_dp*den(l,ie,ipot)*dossgn/efctor/pi/real(nspinpot,kind=dp), l=0, lmaxd1), doscmplx*dossgn/efctor
        end do
        if (ispin/=nspinpot .or. i1/=natyp) write (49, fmt=100)
      end do
    end do
    close (49)

    ! Write total DOS summed over atoms and spins(complex)
    open (49, file='total_cmplx.dos', form='formatted')
    call version_print_header(49)
    write (49, fmt='(4A16)') '# Real(E)', '  Im(E)', ' Re(DEN)', ' Im(DEN)'
    do ie = 1, ielast
      doscmplx = cmplx(0.0e0_dp, 0.e0_dp, kind=dp)
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
