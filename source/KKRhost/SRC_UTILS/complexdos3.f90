module mod_complexdos3

contains

program complexdos
  use :: mod_datatypes, only: dp
  implicit none
  ! Principle of DOS here: Two-point contour integration
  ! for DOS in the middle of the two points. The input DOS
  ! and energy must be complex. Parameter deltae should be
  ! of the order of magnitude of eim.


  ! <-2*deltae->   _
  ! /\        |     DOS=(n(1)+n(2))/2 + (n(1)-n(2))*eim/deltae
  ! /  \       |
  ! (1)  (2)   2*i*eim=2*i*pi*Kb*Tk
  ! /      \     |
  ! /        \    |
  ! ------------------------ (Real E axis)
  integer *4 :: iemaxd, lmaxd
  parameter (iemaxd=1000, lmaxd=10)
  integer *4 :: npot, iemax, lmax
  real (kind=dp) :: eim, deltae, tk, kb, pi
  ! Total dos stored in DOS(LMAX+1,IE)
  complex (kind=dp) :: dos(0:lmaxd+1, iemaxd), ez(iemaxd)
  complex (kind=dp) :: dosnew(0:lmaxd+1)
  real (kind=dp) :: temp(2*lmaxd+6)
  real (kind=dp) :: enew, ef, ev
  integer *4 :: ie, ii, ll, iheader
  character (len=80) :: text
  character (len=45) :: text1
  character (len=20) :: text2

  ev = 13.6058_dp

  ! If only Tk is known, use Bolzmann constant, pi to find eim:
  ! Kb=0.6333659D-5
  ! pi=3.14159265358979312d0
  ! eim=pi*Kb*Tk
  open (49, file='complex.dos', form='formatted')
  open (50, file='new3.dos', form='formatted')
  open (51, file='new3_eV_EF.dos', form='formatted')
  read (49, *) text2               ! dummy readin of header, may be replaced
                                   ! later
  read (49, *) npot
  read (49, *) iemax
  read (49, *) lmax
  if (iemax>iemaxd) stop 'IEMAX.gt.IEMAXD'
  if (lmax>lmaxd) stop 'LMAX.gt.LMAXD'


  do ii = 1, npot
    write (*, *) 'Reading potential', ii
    ! Read header:
    do iheader = 1, 3
      read (49, 120) text
      write (50, 120) text
      write (51, 120) text
      if (iheader==1) write (*, 120) text
    end do
    read (49, fmt='(A45,F10.6,A20)') text1, ef, text2
    write (50, fmt='(A45,F10.6,A20)') text1, ef, text2
    write (51, fmt='(A45,F10.6,A20)') text1, ef, text2
    do iheader = 5, 9
      read (49, 120) text
      write (50, 120) text
      write (51, 120) text
    end do
    ! Read dos: (total dos stored at DOS(LMAX+1,IE))
    do ie = 1, iemax
      read (49, *) ez(ie), (dos(ll,ie), ll=0, lmax+1)
    end do
    ! Compute and write out corrected dos at new (middle) energy points:
    do ie = 2, iemax - 1
      deltae = real(ez(ie+1)-ez(ie))
      eim = aimag(ez(ie))
      enew = real(ez(ie))          ! Real quantity

      do ll = 0, lmax + 1
        dosnew(ll) = dos(ll, ie) + 0.5e0_dp*(dos(ll,ie-1)-dos(ll,ie+1))*cmplx( &
          0.e0_dp, eim, kind=dp)/deltae
      end do

      write (50, 110) enew, aimag(dosnew(lmax+1)), (aimag(dosnew(ll)), ll=0, &
        lmax)
      write (51, 110)(enew-ef)*ev, aimag(dosnew(lmax+1))/ev, &
        (aimag(dosnew(ll))/ev, ll=0, lmax)

    end do                         ! IE=2,IEMAX-1

    if (ii/=npot) then
      read (49, *)
      write (50, 120) ' '
      write (51, 120) ' '
    end if

  end do                           ! II=1,NPOT

  close (49)
  close (50)
  close (51)

  stop 'telos'
110 format (10e14.6)
120 format (a80)
130 format (8(2e12.4))
end program complexdos

end module mod_complexdos3
