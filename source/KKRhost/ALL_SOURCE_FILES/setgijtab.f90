    Subroutine setgijtab(linterface, icc, naez, iqat, rbasis, bravais, &
      natomimp, atomimp, rclsimp, nofgij, ijtabcalc, iofgij, jofgij, nqcalc, &
      iqcalc, natomimpd, ijtabcalc_i)
      Use mod_datatypes, Only: dp
! **********************************************************************
! * Task-specific settings of Gij elements that need to be calculated  *
! * Subroutine (called for ICC=-1) sets up the arrays                  *
! * NATOMIMP    : number of different sites i,j = 1,NATOMIMP           *
! * RCLSIMP     : site coordinates                                     *
! * ATOMIMP     : index of the corresponding site in the unit cell     *
! * IJTABCALC   : flag specifying wehter pair (I,J) needs to be        *
! *               calculated - linear pointer (I-1)*NATOMIMP + J = 1/0 *
! *               for YES/NO                                           *
! * NOFGIJ      : number of all (I,J) pairs - sum of all non-zero I,J  *
! * IOFGIJ      : I index in the list 1..NATOMIMP for pair I,J         *
! * JOFGIJ      : J index                                              *
! **********************************************************************
      Implicit None

!Scalar arguments
      Integer :: icc, naez, natomimp, natomimpd, nofgij, nqcalc
      Logical :: linterface

!Array arguments
      Integer :: atomimp(*), ijtabcalc(*), ijtabcalc_i(*), iofgij(*), iqat(*), &
        iqcalc(*), jofgij(*)
      Real (Kind=dp) :: bravais(3, 3), rbasis(3, *), rclsimp(3, *)

!Local scalars
      Integer :: i, ido, ii, j, jj, nn
      Logical :: opt

!External subroutines
      External :: gijcond, gijxcpl, opt


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(79("="),/,15X,A)') &
        'SETGIJTAB: setting task-specific Gij pairs'
      Write (1337, '(79("="),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      ido = 0
! ======================================================================
      If (opt('CONDUCT ')) Call gijcond(ido, naez, rbasis, iqat, natomimp, &
        rclsimp, atomimp, ijtabcalc, natomimpd)
! ======================================================================
      If (opt('XCPL    ')) Call gijxcpl(ido, naez, rbasis, bravais, &
        linterface, nqcalc, iqcalc, natomimp, rclsimp, atomimp, ijtabcalc, &
        ijtabcalc_i, natomimpd)
! ======================================================================
      If (ido==0) Then
        icc = 0
        Write (6, 110)
        Return
      End If
! ======================================================================
      nofgij = 0
      Do i = 1, natomimp
        nn = (i-1)*natomimp
        Do j = 1, natomimp
          If (ijtabcalc(nn+j)>0) Then
            nofgij = nofgij + 1
            If (nofgij>natomimpd*natomimpd) Then
              Write (6, 100) 'NATOMIMPD', nofgij/natomimp
              Stop
            End If
            iofgij(nofgij) = i
            jofgij(nofgij) = j
          End If
        End Do
      End Do
      If (nofgij==0) Then
        icc = 0
        Write (6, 110)
        Return
      End If

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 120) natomimp, nofgij
      Write (1337, 130)
      Write (1337, 140)
      Write (1337, 130)
      Do i = 1, nofgij
        ii = iofgij(i)
        jj = jofgij(i)
        Write (1337, 150) i, ii, atomimp(ii), (rclsimp(j,ii), j=1, 3), jj, &
          atomimp(jj), (rclsimp(j,jj), j=1, 3)
      End Do
      Write (1337, 130)
      Write (1337, *)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

100   Format (6X, 'brahim ERROR: please increase the global parameter', /, 6X, &
        A, ' to a value >=', I5, /)
110   Format (6X, 'WARNING: Subroutine entered with invalid task ', &
        'specification', /, 6X, &
        '         ICC will be set to 0 - no Gij calculated - ', &
        'input check? ', /)
120   Format (6X, 'Number of different sites (NATOMIMP) :', I4, /, 6X, &
        'Number of pairs set       (NOFGIJ)   :', I4)
130   Format (8X, 71('-'))
140   Format (9X, 'pair|', ' I  IQ           position', 9X, &
        'J  JQ           position')
150   Format (9X, I3, ' |', 2(I3,1X), 3F8.4, 1X, 2(I3,1X), 3F8.4)
160   Format (I5, 2(I5,1X), 3F10.6, 1X, 2(I5,1X), 3F10.6)
    End Subroutine
