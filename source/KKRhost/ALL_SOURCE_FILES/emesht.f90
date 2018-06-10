!-------------------------------------------------------------------------------
! SUBROUTINE: EMESHT
!> @brief This subroutine provides the energy mesh in array EZ and the
!> appropriate integration weights in array DF.
!> @details Poles of the Fermi function C (Matsubara frequencies) and a contour in
!> the complex energy are used as described in (????).
!> The contour consists of three straight lines with NPNT1, NPNT2, and NPNT3
!> integration points and is determined by the input arguments: EBOT, EMU, TK, and NPOL.
!> The three lines are defined by:
!> 1. The line from EBOT to \f$ EBOT+2*NPOL*\pi*i*k*TK \f$ with NPNT1 integration points (Gauss-Legendre rule)
!> 2. The line from \f$ EBOT+2*NPOL*\pi*i*k*TK\f$ to \f$ EMU+(2*NPOL*\pi*i-30)*k*TK\f$ with NPNT2 integration points (Gauss-Legendre rule)
!> 3. The line from \f$ EMU+(2*NPOL*\pi*i-30)*k*TK\f$ to \f$ \infty \f$
!>
!> The total number of integration points is given by: \f$ NPNT=NPNT1+NPNT2+NPNT3+NPOL\f$
!> The integration points and weights on three lines are chosen according to Gauss integration rules. Only in third interval
!> the Fermi function matters since \f$ e^x < 10^{-10} \f$ for \f$ x < -25\f$.
!> There are two special cases determined by NPOL = 0 and NPOL < 0.
!> - NPOL = 0 leads to density-of-states calculations with constant integration weights and equally distributed points
!> between \f$ EBOT - \pi*i*k*TK\f$ and \f$ EMU - \pi*i*k*TK\f$.
!> The total number of integration points is given by: NPNT=NPNT2
!> - NPOL < 0 is meant for calculations where the Fermi-Dirac function is replaced by a step function with step at EMU. When
!> this option is used no poles of the Fermi-Dirac function are used and the contour consists of the three straight lines:
!> 1. The line from \f$EBOT\f$ to \f$ EBOT-2*NPOL*\pi*i*k*TK\f$ with NPNT1 integration points (Gauss-Legendre rule)
!> 2. The line from \f$EBOT-2*NPOL*\pi*i*k*TK\f$ to \f$EMU-2*NPOL*\pi*i*k*TK\f$ with NPNT2 integration points (Gauss-Legendre rule)
!> 3. The line from \f$ EMU-2*NPOL*\pi*i*k*TK\f$ to \f$EMU\f$ with NPNT3 integration points (Gauss-Legendre rule)
!>
!> The total number of integration points is given by: \f$NPNT=NPNT1+NPNT2+NPNT3\f$
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine emesht(ez, df, npnt, ebot, emu, efermi, tk, npol, npnt1, npnt2, &
  npnt3, iemxd)

  use :: mod_types, only: t_inc
  use :: constants

  implicit none
! ..
! .. Input variables
  integer, intent (in) :: npol !< Number of Matsubara Poles (EMESHT)
  integer, intent (in) :: npnt1 !< number of E points (EMESHT) for the contour integration
  integer, intent (in) :: npnt2 !< number of E points (EMESHT) for the contour integration
  integer, intent (in) :: npnt3 !< number of E points (EMESHT) for the contour integration
  integer, intent (in) :: iemxd !< Dimension for energy-dependent arrays
  double precision, intent (in) :: tk !< Temperature
  double precision, intent (in) :: emu !< Top of the contour
  double precision, intent (in) :: ebot !< Bottom of the contour
  double precision, intent (in) :: efermi !< Fermi energy
! .. Input/Output variables
  integer, intent (inout) :: npnt
  double complex, dimension (iemxd), intent (inout) :: df
  double complex, dimension (iemxd), intent (inout) :: ez
! .. Local Scalars ..
  integer :: i
  double complex :: de
  double precision :: er, etk
! .. Local Arrays ..
  double precision, dimension (128) :: wi, xi
! .. External Subroutines ..
  logical :: opt
  external :: gaufd, gauleg, opt
! ..
! .. Intrinsic Functions ..
  intrinsic :: dcmplx
! ..
!----------------------------------------------------------------------------
! OUTPUT
!----------------------------------------------------------------------------
  if (t_inc%i_write>0) then
    write (1337, '(5X,A,F12.6," (Ry)",8X,A,F12.6," (Ry)")') 'E min = ', ebot, &
      'Fermi energy = ', efermi
    write (1337, '(5X,A,F12.6," (Ry)",8X,A,F12.6," (K )",/,5X,62(1H-))') &
      'E max = ', emu, 'Temperature  = ', tk
  end if
!----------------------------------------------------------------------------
! OUTPUT
!----------------------------------------------------------------------------
  etk = pi*kb*tk
!----------------------------------------------------------------------------
  if (npol==0) then
    de = (emu-ebot)
    if (npnt2>1) then
      de = de/(npnt2-1)
    else
      de = dcmplx(1.0d0, 0.0d0)
    end if
    npnt = 0
    do i = 1, npnt2
      npnt = npnt + 1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      er = ebot + (i-1)*de
      ez(npnt) = dcmplx(er, etk)
      df(npnt) = de
    end do ! I
    if (t_inc%i_write>0) write (1337, fmt=100) npnt, etk, etk*ryd
!----------------------------------------------------------------------------
! NPOL > 0
!----------------------------------------------------------------------------
  else if (npol>0) then
    call gauleg(xi, wi, npnt1)
    de = npol*dcmplx(0.0d0, etk)
    npnt = 0
    do i = 1, npnt1
      npnt = npnt + 1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      ez(npnt) = xi(i)*de + de + ebot
      df(npnt) = wi(i)*de
    end do ! I -> NPNT1
    call gauleg(xi, wi, npnt2)
    de = (emu-30*kb*tk-ebot)*0.5d0
    do i = 1, npnt2
      npnt = npnt + 1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      ez(npnt) = xi(i)*de + de + ebot + 2*npol*dcmplx(0.0d0, etk)
      df(npnt) = wi(i)*de
    end do ! I -> NPTN2
    call gaufd(xi, wi, npnt3)
    de = 30*kb*tk
    do i = 1, npnt3
      npnt = npnt + 1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      ez(npnt) = xi(i)*de + emu + 2*npol*dcmplx(0.0d0, etk)
      df(npnt) = wi(i)*de
    end do ! I - >NPTN3
    do i = npol, 1, -1
      npnt = npnt + 1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      ez(npnt) = emu + (2*i-1)*dcmplx(0.0d0, etk)
      df(npnt) = -2*dcmplx(0.0d0, etk)
    end do
    if (t_inc%i_write>0) write (1337, 110) npnt, npol, npnt1, npnt2, npnt3
!-------------------------------------------------------------------------
! NPOL < 0
!-------------------------------------------------------------------------
  else
    if (npnt1>0) call gauleg(xi, wi, npnt1)
    de = -npol*dcmplx(0.0d0, etk)
    npnt = 0
    do i = 1, npnt1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      npnt = npnt + 1
      ez(npnt) = xi(i)*de + de + ebot
      df(npnt) = wi(i)*de
    end do ! I -> NPNT1
    call gauleg(xi, wi, npnt2)
    de = (emu-ebot)*0.5d0
    do i = 1, npnt2
      npnt = npnt + 1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      ez(npnt) = xi(i)*de + de + ebot - 2*npol*dcmplx(0.0d0, etk)
      if (opt('GF-EF   ')) ez(npnt) = emu + npol*dcmplx(0.0d0, etk)
      df(npnt) = wi(i)*de
    end do ! I -> NPNT2
    if (npnt3>0) call gauleg(xi, wi, npnt3)
    de = -npol*dcmplx(0.0d0, etk)
    do i = npnt3, 1, -1
      npnt = npnt + 1
      if (npnt>iemxd) then
        write (6, '(/,5X,2A,I4)') &
          'Dimension ERROR: Increase IEMXD in the inputcard to ', 'at least ', &
          npnt
        stop '     < EMESHT >'
      end if
      ez(npnt) = xi(i)*de + de + emu
      df(npnt) = -wi(i)*de
    end do ! I -> NPNT3
    if (t_inc%i_write>0) write (1337, 120) npnt, -npol, npnt1, npnt2, npnt3
  end if
!----------------------------------------------------------------------------
  if (t_inc%i_write>0) write (1337, *)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correction Factor for the weight in the integration according to Phivos Idea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1, npnt
    df(i) = df(i)
  end do
!*(8.5D0/8.48686D0)*(8.75D0/8.74083D0)
! GaCrN*(8.5D0/8.49286D0)
!*(8.5D0/8.48969D0)
!*(8.5D0/8.48823D0)
!*(8.75D0/8.73983D0)
!*(8.5D0/8.48686D0)
!*(8.75D0/8.75659D0)
!*(6.5D0/6.55253D0)*(7.5D0/7.47798D0)
!     *(8.75D0/8.75659D0)
!*(8.5D0/8.54963D0)
!*(6.5D0/6.41299D0)
!*(8.5D0/8.47767D0)
!*(6.5D0/6.45787D0)
!*(4.D0/4.01579D0)
!*(8.8D0/8.80272D0)*(8.8D0/8.78691D0)
!*(4.0D0/4.0419213D0)*(17.5D0/17.508D0)
!     &  *(8.0D0/7.9885D0)*(8.75D0/8.74682D0)*(8.0D0/7.9246D0)
!     &  *(8.25D0/8.24085)
! 90        write(*,*)'DF=',I,DF(I)
!*************************************************************
!**********************************************************

  return
100 format (5x, 'Density-of-States calculation', /, 5x, &
    'Number of energy points :', i4, 4x, 'broadening =', 3p, f9.3, ' ( mRy )', &
    /, 48x, ' =', 3p, f9.3, ' ( meV )')
110 format (5x, 'GF integration rectangular contour ( ImE > 0 )', /, 5x, &
    'Number of energy points :', i4, 13x, 'poles =', i2, /, 23x, &
    'contour: N1 =', i2, ', N2 =', i4, ', N3 =', i2)
120 format (5x, 'GF integration rectangular contour ( ImE < 0 )', /, 5x, &
    'Number of energy points :', i4, 13x, 'poles =', i2, /, 23x, &
    'contour: N1 =', i2, ', N2 =', i4, ', N3 =', i2)
end subroutine
