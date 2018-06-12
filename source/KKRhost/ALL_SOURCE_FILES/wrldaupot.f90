subroutine wrldaupot(itrunldau, lopt, ueff, jeff, erefldau, natyp, wldau, &
  uldau, phildau, irmd, natypd, nspind, mmaxd, irws)
  ! **********************************************************************
  ! *                                                                    *
  ! * Writes out LDA+U arrays into formatted file 'ldaupot'              *
  ! *                                                                    *
  ! **********************************************************************
  use :: mod_version_info
  use :: mod_datatypes, only: dp
  implicit none
  ! ..
  integer :: irmd, mmaxd, natypd, nspind, irws(natypd)
  ! ..
  ! .. Arguments ..
  integer :: itrunldau, natyp
  integer :: lopt(natypd)
  real (kind=dp) :: ueff(natypd), jeff(natypd), erefldau(natypd)
  real (kind=dp) :: wldau(mmaxd, mmaxd, nspind, natypd)
  real (kind=dp) :: uldau(mmaxd, mmaxd, mmaxd, mmaxd, natypd)
  complex (kind=dp) :: phildau(irmd, natypd)
  ! ..
  ! ..  Locals
  integer :: ir, m1, m2, m3, m4, it, is
  ! ======================================================================

  open (67, file='ldaupot_new', form='FORMATTED')
  call version_print_header(67)
  write (1337, 100)
  write (67, 110) itrunldau, '    ITRUNLDAU'
  write (67, 110) natyp, '    NATYP'
  write (67, 120) natyp
  write (67, 130)(lopt(it), it=1, natyp)
  write (67, 140)
  do it = 1, natyp
    if (lopt(it)+1/=0) write (67, 150) it, ueff(it), jeff(it), erefldau(it)
  end do
  do it = 1, natyp
    if (lopt(it)+1/=0) then
      write (67, 110) it, '    WLDAU'
      do is = 1, nspind
        do m1 = 1, mmaxd
          write (67, 160)(wldau(m1,m2,is,it), m2=1, mmaxd)
        end do
      end do
      write (67, 110) it, '    ULDAU'
      write (67, 160)((((uldau(m1,m2,m3,m4,it),m4=1,mmaxd),m3=1, &
        mmaxd),m2=1,mmaxd), m1=1, mmaxd)
      write (67, 110) it, '    PHILDAU'
      write (67, 160)(phildau(ir,it), ir=1, irws(it))
    end if
  end do
  close (67)

100 format (/, 5x, '< WRLDAUPOT > : ', &
    'Writing out LDA+U potential (file ldaupot_new)', /)
110 format (i6, a)
120 format ('LOPT 1..', i3)
130 format (16i3)
140 format ('IAT', 6x, 'UEFF', 12x, 'JEFF', 12x, 'EREF')
150 format (i3, 3(1x,e15.8))
160 format (5e16.8)
end subroutine wrldaupot
