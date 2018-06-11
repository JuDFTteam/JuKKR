    Subroutine errortrap(routine, k, istop)

      Use mod_datatypes, Only: dp
      Implicit None
      Integer, Parameter :: kmax = 14
      Integer, Intent (In) :: k, istop
      Character (Len=*), Intent (In) :: routine
      Character (Len=60) :: t(kmax), text

      Data t/'program KKRSCF called for TASK <> SCF    check input file   ', &
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

      If (k>=1 .And. k<=kmax) Then
        text = t(k)
      Else
        text = 'unknown reason   key out of bounds'
      End If

      If (istop==1) Then
        Write (6, 100) routine, text
        Stop
      Else
        Write (6, 110) routine, text
      End If

100   Format (/, 1X, 79('*'), /, 1X, 79('*'), /, 5X, 'STOP in subroutine <', &
        A, '>', /, 5X, A, /, 1X, 79('*'), /, 1X, 79('*'), /)
110   Format (/, 1X, 79('<'), /, 5X, 'warning from subroutine <', A, '>', /, &
        5X, A, /, 1X, 79('>'), /)
      Return
    End Subroutine
