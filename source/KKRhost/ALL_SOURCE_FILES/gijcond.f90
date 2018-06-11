    Subroutine gijcond(ido, naez, rbasis, iqat, natomimp, rclsimp, atomimp, &
      ijtabcalc, natomimpd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * In case of tasks requiring Gij blocks calculation, set variables:  *
! *                                                                    *
! * NATOMIMP, RCLSIMP(3,1..NATOMIMP), ATOMIMP(1..NATOMIMP)             *
! * IJTABCALC flag to which pair is needed: I,J --> (I-1)*NATOMIMP + J *
! * IDO takes on the value 1 or 0 if setting up process was OK or not  *
! *                                                                    *
! * CONDUCTANCE calculation case                                       *
! *             still to be implemente the correct read in             *
! **********************************************************************
      Implicit None

!Parameters
      Integer :: ncpaird
      Parameter (ncpaird=10)

!Arguments
      Integer :: ido, naez, natomimp, natomimpd
      Integer :: atomimp(*), ijtabcalc(*), iqat(*)
      Real (Kind=dp) :: rbasis(3, *), rclsimp(3, *)

!Locals
      Integer :: i, iat, iatcondl(ncpaird), iatcondr(ncpaird), j, jat, &
        ncondpair, nn

      ido = 0
      Do i = 1, ncpaird
        iatcondl(i) = 0
        iatcondr(i) = 0
      End Do

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 100)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

!     ---------------------------------------------------- dummy
!     settings so far, need to be replaced by conductance input
!     and some output
      ncondpair = 4
      If (ncondpair>ncpaird) Then
        Write (6, 110) 'local', 'NCPAIRD', ncondpair
        Stop
      End If
      iatcondl(1) = 1
      iatcondr(1) = 2
      iatcondl(2) = 1
      iatcondr(2) = 2
      iatcondl(3) = 2
      iatcondr(3) = 1
      iatcondl(4) = 2
      iatcondr(4) = 1

!     ---------------------------------------------------- dummy
      If (ncondpair==0) Return
      Do i = 1, ncondpair
        If ((iatcondl(i)<=0) .Or. (iatcondl(i)>naez)) Return
        If ((iatcondr(i)<=0) .Or. (iatcondr(i)>naez)) Return
      End Do

      natomimp = 2*ncondpair
      If (natomimp>natomimpd) Then
        Write (6, 110) 'global', 'NATOMIMPD', natomimp
        Stop
      End If

      Do i = 1, natomimp
        nn = (i-1)*natomimp
        Do j = 1, natomimp
          ijtabcalc(nn+j) = 0
        End Do
      End Do

      nn = 0
      Do i = 1, ncondpair
        iat = iqat(iatcondl(i)) ! left lead
        nn = nn + 1
        Do j = 1, 3
          rclsimp(j, nn) = rbasis(j, iat)
        End Do
        atomimp(nn) = iat
        iat = nn

        jat = iqat(iatcondr(i)) ! right lead
        nn = nn + 1
        Do j = 1, 3
          rclsimp(j, nn) = rbasis(j, jat)
        End Do
        atomimp(nn) = jat
        jat = nn
        ijtabcalc((iat-1)*natomimp+jat) = 1
      End Do
      If (natomimp/=nn) Then
        Write (6, '(6X,A,/,6X,A,/)') &
          'ERROR: Found some inconsistencies in IATCOND arrays', &
          '       Please check your CONDUCTANCE input'
        Stop
      End If
      ido = 1
100   Format (5X, '< GIJCOND > : Conductance/conductivity calculation', /)
110   Format (6X, 'Dimension ERROR: please increase the ', A, ' parameter', /, &
        6X, A, ' to a value >=', I5, /)
    End Subroutine
