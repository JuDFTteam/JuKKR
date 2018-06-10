SUBROUTINE errortrap(routine,k,istop)
 
IMPLICIT NONE
INTEGER, PARAMETER           :: KMAX=14
INTEGER, intent(in)          :: K, ISTOP
CHARACTER(len=*), intent(in) :: ROUTINE
CHARACTER(len=60)            :: T(KMAX), TEXT

DATA t / 'program KKRSCF called for TASK <> SCF    check input file   ',  &
    'TASK = SCF in input   but not program KKRSCF called         ',  &
    'ITEST should be <=4                                         ',  &
    'DATASET not initialized           >>>> check input file     ',  &
    'POTFIL  not initialized           >>>> check input file     ',  &
    'SCFSTART for   PROGRAM <> KKRSCF  >>>>  call  kkrscf instead',  &
    'SUM OF CONC <> 1                                            ',  &
    'all SOCTL should be >= 0          >>>> check input file     ',  &
    'all SOCTL should be < 0           >>>> check input file     ',  &
    'use BZINT = POINTS    for SP-SREL case and for spin spirals ',  &
    'I < > NLM   for IREL = 2                                    ',  &
    'routine called for IREL = 2   deal with that case outside   ',  &
    'routine called for IREL = 3   and   even  NK                ',  &
    'anti-unitary symmetry matrices created for IREL = 2         '/

IF( k >= 1 .AND. k <= kmax) THEN
  text = t(k)
ELSE
  text = 'unknown reason   key out of bounds'
END IF

IF( istop == 1 ) THEN
  WRITE(6,9000) routine, text
  STOP
ELSE
  WRITE(6,9001) routine, text
END IF

9000 FORMAT(/,1X,79('*'),/,1X,79('*'),/,5X,'STOP in subroutine <',  &
    a,'>',/,5X,a,/,1X,79('*'),/,1X,79('*'),/)
9001 FORMAT(/,1X,79('<'),/,5X,'warning from subroutine <',a,'>',/,  &
    5X,a,/,1X,79('>'),/)
RETURN
END SUBROUTINE errortrap
