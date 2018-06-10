FUNCTION ikapmue(kappa,muem05)
!   ********************************************************************
!   *                                                                  *
!   *  INDEXING OF MATRIX-ELEMENTS:                                    *
!   *                                                                  *
!   *  I = 2*L*(J+1/2) + J + MUE + 1                                   *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER KAPPA,MUEM05
INTEGER IKAPMUE

! Local variables
INTEGER IABS
INTEGER JP05,L

jp05 = IABS(kappa)

IF ( kappa < 0 ) THEN
  l = -kappa - 1
ELSE
  l = kappa
END IF

ikapmue = 2*l*jp05 + jp05 + muem05 + 1

END FUNCTION ikapmue
