    Subroutine gfmask(linterface, icheck, icc, invmod, nsh1, nsh2, naez, &
      nshell, naezd, nprincd)
! **********************************************************************
! *                                                                    *
! * This subroutine prepares the ICHECK matrix that is used for        *
! * calculating the proper off-diagonal GF matrix elements ( e.g.      *
! * impurity) in case of no full inversion algorithm                   *
! *                                                                    *
! * ICHECK(I,J) points to the block (I,J) of the GF matrix having the  *
! * size NPRINCD                                                       *
! *                                                                    *
! *                                            29.02.2000              *
! *                                                                    *
! **********************************************************************
      Use mod_datatypes, Only: dp
      Implicit None
!     ..
!     .. Scalar arguments
      Integer :: naezd, nprincd
      Integer :: icc, invmod, nlayer, naez, nshell
      Logical :: linterface
!     ..
!     .. Array arguments
      Integer :: icheck(naezd/nprincd, naezd/nprincd)
      Integer :: nsh1(*), nsh2(*)
!     .. Local variables
      Integer :: icouple(naezd, naezd)
      Integer :: i, j, k, ii, istep1, ilt1, istep2, ilt2, il2, il1, lfchk
      Character (Len=80) :: fmtchk
      Character (Len=35) :: invalg(0:2)
!     ..
!     .. External functions
      Logical :: opt, test
      External :: opt, test
!     ..
!     .. Data statements
      Data invalg/'FULL MATRIX                        ', &
        'BANDED MATRIX (slab)               ', &
        'BANDED + CORNERS MATRIX (supercell)'/

      Write (1337, 100)

! --> set default inversion to SUPERCELL mode = banded matrix + corners

      invmod = 2

! --> LINTERFACE = use band diagonal mode

      If (linterface) invmod = 1

! --> full inversion is performed ONLY BY EXPLICIT request

      If (opt('full inv')) invmod = 0

      If ((invmod/=0) .And. (mod(naez,nprincd)/=0)) Then
        Write (6, 110) naez, nprincd
        Stop
      End If

      Write (1337, 120) invalg(invmod)

      nlayer = naez/nprincd
! ----------------------------------------------------------- INVMOD = 1
!                                                   band-diagonal matrix
      If (invmod==1) Then
        Do i = 1, nlayer
          Do j = 1, nlayer
            If (i==j) Then
              icheck(i, j) = 1
            Else
              icheck(i, j) = 0
            End If
          End Do
        End Do
      End If
! ----------------------------------------------------------- INVMOD = 2
!              band-diagonal matrix with corners (slab periodic along z)
      If (invmod==2) Then
        Do i = 1, nlayer
          Do j = 1, nlayer
            If ((i==j) .Or. (j==nlayer) .Or. (i==nlayer)) Then
              icheck(i, j) = 1
            Else
              icheck(i, j) = 0
            End If
          End Do
        End Do
      End If
! ================================================= INVMOD = 1, ICC <> 0
!                   band-diagonal matrix, off-diagonal G elements needed

! --> prepare the matrix ICOUPLE which has 1 in all nn' blocks
!     (atomic sites) that are needed

! ======================================================================
      If ((icc/=0) .And. (invmod==1)) Then
        Do i = 1, naez
          Do j = 1, naez
            icouple(i, j) = 0

            Do ii = 1, nshell
              If (((nsh1(ii)==i) .And. (nsh2(ii)==j)) .Or. ((nsh1(ii)== &
                j) .And. (nsh2(ii)==i))) icouple(i, j) = 1
            End Do
          End Do
        End Do
!cccC ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!cccC                                               conductivity calculation
!ccc         IF (OPT('CONDUCT ')) THEN
!ccc            DO I=1,NLAYER
!ccc               DO J=1,NLAYER
!ccc                  ICHECK(I,J)=0
!ccc               ENDDO
!ccc            ENDDO
!ccc            DO I=1,NAEZ
!ccc               DO J=1,NAEZ
!ccc                  ICOUPLE(I,J) = 0
!ccc               END DO
!ccc            END DO
!ccc            DO II=1,NCONDPAIR
!ccc               I = IATCONDL(II)
!ccc               J = IATCONDR(II)
!ccc               ICOUPLE(I,J) = 1
!ccc               ICOUPLE(J,I) = 1
!ccc            ENDDO
!ccc         END IF                 ! Conductivity calculation
!cccC ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! ----------------------------------------------------------------------
! Now given the matrix ICOUPLE prepare the matrix ICHECK which has 1 in
! all principal-layer blocks that we need -- this will be used in the
! matrix inversion
! ----------------------------------------------------------------------
        istep1 = 0
        ilt1 = 1
! ----------------------------------------------------------------------
        Do il1 = 1, naez
          istep1 = istep1 + 1

          If (istep1>nprincd) Then
            ilt1 = ilt1 + 1
            istep1 = 1
          End If

          ilt2 = 1
          istep2 = 0
! ......................................................................
          Do il2 = 1, naez
            istep2 = istep2 + 1

            If (istep2>nprincd) Then
              ilt2 = ilt2 + 1
              istep2 = 1
            End If

            If (icouple(il1,il2)==1) icheck(ilt1, ilt2) = 1
          End Do
! ......................................................................
        End Do
! ----------------------------------------------------------------------
! in the case of calculation of single blocks it has to put the correct
! value to ICHECK in order to calculate all the elements also necessary
! to calculate that single block          ?????
! ----------------------------------------------------------------------
        Do j = 1, nlayer

! --> loop over the element ICHECK(I,J) with fixed J and I < J

          If (j/=1) Then
            Do i = 1, j - 1
              If (icheck(i,j)==1) Then
                Do k = i + 1, j
                  icheck(k, j) = 1
                End Do
                Do k = j, nlayer
                  icheck(k, k) = 1
                End Do
              End If
            End Do
          End If

          If (.Not. opt('CONDUCT ')) Then

! --> loop over the element ICHECK(I,J) with fixed J and I > J

            If (j/=nlayer) Then
              Do i = nlayer, j + 1, -1
                If (icheck(i,j)==1) Then
                  Do k = i - 1, j, -1
                    icheck(k, j) = 1
                  End Do
                End If
              End Do
            End If
          End If
        End Do
! ----------------------------------------------------------------------
      End If
! ======================================================================

      If (test('ICHECK  ')) Then

        fmtchk = ' '
        lfchk = 1
        Do i = 1, min(35, nlayer)
          fmtchk = fmtchk(1:lfchk) // '--'
          lfchk = lfchk + 2
        End Do

        Write (1337, '(8X,A,/,8X,A)') 'ICHECK matrix :', fmtchk(1:lfchk)
        Do i = 1, nlayer
          Write (1337, '(9X,35I2)')(icheck(i,j), j=1, min(35,nlayer))
        End Do
        Write (1337, '(8X,A,/)') fmtchk(1:lfchk)
      End If

100   Format (5X, '< GFMASK > : set KKR matrix inversion algorithm', /)
110   Format (6X, 'ERROR: Number of sites (NAEZ) =', I3, &
        ' not an integer multiplier', /, 6X, 'of principal layers (NPRINCD) =' &
        , I3, /, 6X, 'Use ONLY  full inversion in this case')
120   Format (8X, 'INVERSION algorithm used : ', A35, /)
    End Subroutine
