module mod_outtmathost

contains

subroutine outtmathost(alat, ins, krel, kmrot, nspin, naez, lmmax, bravais, &
  rbasis, qmtet, qmphi, e2in, tk, npol, npnt1, npnt2, npnt3)
  ! **********************************************************************
  ! *                                                                    *
  ! *  Writes out the header of the t-matrices decimation file           *
  ! *                                                                    *
  ! **********************************************************************
  use :: mod_version_info
  use :: mod_datatypes, only: dp
  implicit none
  ! ..
  ! .. Arguments ..
  integer :: ins, krel, kmrot, nspin, naez, lmmax, npol, npnt1, npnt2, npnt3
  real (kind=dp) :: alat, e2in, tk
  real (kind=dp) :: bravais(3, 3), rbasis(3, *), qmtet(*), qmphi(*)
  ! ..
  ! .. Locals ..
  integer :: i, ih
  ! ----------------------------------------------------------------------
  write (1337, '(5X,A,/)') '< DECIOPT > : writing header of decimation file'

  open (37, file='decifile', status='unknown')
  call version_print_header(37)
  write (37, fmt=*) 'INVERSE T-MATRIX AND CMOMS'
  write (37, fmt=100)
  write (37, fmt=110) alat, nspin, naez, lmmax, ins, krel, kmrot
  write (37, fmt=120) bravais
  if (krel==0) then
    write (37, fmt=130)
    do ih = 1, naez
      write (37, fmt=150)(rbasis(i,ih), i=1, 3)
    end do
  else
    write (37, fmt=140)
    do ih = 1, naez
      write (37, fmt=160)(rbasis(i,ih), i=1, 3), qmtet(ih), qmphi(ih)
    end do
  end if
  write (37, fmt=170) e2in, tk
  write (37, fmt=180) npnt1, npnt2, npnt3, npol
  close (37)
  ! ----------------------------------------------------------------------
100 format (' Vectors in lattice constant units', /, &
    '                                 ')
110 format ('ALAT=', f9.6, ' NSPIN=', i2, '  NAEZ=', i3, ' LMMAX=', i3, &
    ' INS=', i1, ' KREL=', i1, ' KMROT=', i1)
120 format ('BRAVAIS ', /, 3f8.4, /, 3f8.4, /, 3f8.4)
130 format ('RBASIS')
140 format ('RBASIS', 20x, 'MAGNETISATION ANGLES THETA/PHI')
150 format (3f8.4)
160 format (3f8.4, 2f9.4)
170 format ('EF=', f10.6, ' TEMP=', f10.4, ' Kelvin')
180 format ('N1=', i3, ' N2=', i3, ' N3=', i3, ' NPOL=', i3)
end subroutine outtmathost

end module mod_outtmathost
