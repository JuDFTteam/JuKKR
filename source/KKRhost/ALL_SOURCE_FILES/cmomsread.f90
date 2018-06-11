    Subroutine cmomsread(nlbasis, nrbasis, naez, cmomhost, vacflag, kaoez, &
      natypd, nembd1, lmpotd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * This subroutine reads in the CMOMHOST from the decimation files    *
! * Note that they are needed only in case of an SCF decimation cal-   *
! * culation (SCFSTEPS > 1 )                                           *
! *                                                                    *
! * The t-matrices are writen out in kloopz1  (option 'deci-out')      *
! *                                                                    *
! * This subroutine must be called after the t-matrices for all the    *
! * energies are read in (see < decimaread > )                         *
! * It returns the CMOMHOST array. First the left bulk (unit 37) then  *
! * the right bulk (unit 38) are indexed.                              *
! * CMOMHOST(*,NEMBD1) =                                               *
! *                CMOMHOST(*,1..NLBASIS,NLBASIS+1..NLBASIS+NRBASIS)   *
! * Condider this mapping for further use.                             *
! *                                                                    *
! *                                                  29.10.99          *
! *                                                  05.06.04          *
! **********************************************************************
      Implicit None
!..
!.. Arguments
      Integer :: lmpotd, natypd, nembd1
      Integer :: naez, nlbasis, nrbasis
      Real (Kind=dp) :: cmomhost(lmpotd, nembd1)
      Integer :: kaoez(natypd, *)
      Logical :: vacflag(2)
!..
!.. Local variables ..
      Real (Kind=dp) :: c00(lmpotd)
      Character (Len=5) :: chhost(2)
      Integer :: ih, ih1, ihl, ihost, lm, lmpotl, naezl, nathost
!..
!.. Data statements
      Data chhost/'LEFT ', 'RIGHT'/
!..
      Write (1337, '(5X,A,/,8X,30("-"),/,8X,3A6,A10,/,8X,30("-"))') &
        'Reading in host charge moments ( SCFSTEPS > 1 )', ' HOST ', '  IBAS', &
        '  ATOM', '   CMOM(1)'
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      Do ihost = 1, 2
        nathost = nlbasis
        If (ihost==2) nathost = nrbasis
        Write (1337, '(8X,A5,1X)', Advance='no') chhost(ihost)
! ----------------------------------------------------------------------
        If (vacflag(ihost)) Then
          Do ih = 1, nlbasis
            Do lm = 1, lmpotd
              cmomhost(lm, (ihost-1)*nlbasis+ih) = 0.E0_dp
! mapping the CMOMHOST array, ordering is important
            End Do
          End Do
          Write (1337, '(A)') ' Vacuum setting    0.000'
          If (ihost==1) Then
            Write (1337, '(14X,24("-"))')
          Else
            Write (1337, '(8X,30("-"))')
          End If
! ----------------------------------------------------------------------
        Else
          Read (36+ihost, 110) naezl, lmpotl
! ......................................................................
          If (naezl/=nathost) Then
            Write (6, '(/,5X,2A)') 'ERROR: ', &
              'host not compatible with your input.'
            Write (6, '(/,12X,A,I3,A,I3)') 'Charge moments tabulated for', &
              naezl, ' host atoms, input NBASIS =', nathost
            Stop '       < CMOMSREAD > '
          End If
! ......................................................................
          Do ih = 1, naezl
            Read (36+ihost, *) ihl
            If (ihl/=ih) Then
              Write (6, '(/,5X,2A,/)') 'ERROR reading host file', &
                ' basis indexing wrong'
              Stop '       < CMOMSREAD > '
            End If
            Read (36+ihost, 100)(c00(lm), lm=1, lmpotl)
            ih1 = kaoez(1, naez+(ihost-1)*nlbasis+ih)

            Do lm = 1, lmpotl
              cmomhost(lm, (ihost-1)*nlbasis+ih) = c00(lm)
! mapping the CMOMHOST array, ordering is important
            End Do

            If (ih==1) Then
              Write (1337, '(1X,2I6,D12.4)') ih, ih1, c00(1)
            Else
              Write (1337, '(14X,2I6,D12.4)') ih, ih1, c00(1)
            End If
          End Do
! ......................................................................
          If (ihost==1) Then
            Write (1337, '(14X,24("-"))')
          Else
            Write (1337, '(8X,30("-"))')
          End If
        End If
! ----------------------------------------------------------------------
      End Do
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      Write (1337, *)

100   Format (4D22.14)
110   Format (5X, 2I6)
    End Subroutine
