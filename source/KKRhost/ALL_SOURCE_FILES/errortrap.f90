subroutine errortrap(routine, k, istop)

  implicit none
  integer, parameter :: kmax = 14
  integer, intent (in) :: k, istop
  character (len=*), intent (in) :: routine
  character (len=60) :: t(kmax), text

  data t/'program KKRSCF called for TASK <> SCF    check input file   ', &
    'TASK = SCF in input   but not program KKRSCF called         ', &
    'ITEST should be <=4                                         ', &
    'DATASET not initialized           >>>> check input file     ', &
    'POTFIL  not initialized           >>>> check input file     ', &
    'SCFSTART for   PROGRAM <> KKRSCF  >>>>  call  kkrscf instead', &
    'SUM OF CONC <> 1                                            ', &
    'all SOCTL should be >= 0          >>>> check input file     ', &
    'all SOCTL should be < 0           >>>> check input file     ', &
    'use BZINT = POINTS    for SP-SREL case and for spin spirals ', &
    'I < > NLM   for IREL = 2                                    ', &
    'routine called for IREL = 2   deal with that case outside   ', &
    'routine called for IREL = 3   and   even  NK                ', &
    'anti-unitary symmetry matrices created for IREL = 2         '/

  if (k>=1 .and. k<=kmax) then
    text = t(k)
  else
    text = 'unknown reason   key out of bounds'
  end if

  if (istop==1) then
    write (6, 100) routine, text
    stop
  else
    write (6, 110) routine, text
  end if

100 format (/, 1x, 79('*'), /, 1x, 79('*'), /, 5x, 'STOP in subroutine <', a, &
    '>', /, 5x, a, /, 1x, 79('*'), /, 1x, 79('*'), /)
110 format (/, 1x, 79('<'), /, 5x, 'warning from subroutine <', a, '>', /, 5x, &
    a, /, 1x, 79('>'), /)
  return
end subroutine
