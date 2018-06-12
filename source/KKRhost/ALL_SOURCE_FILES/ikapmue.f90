function ikapmue(kappa, muem05)
  ! ********************************************************************
  ! *                                                                  *
  ! *  INDEXING OF MATRIX-ELEMENTS:                                    *
  ! *                                                                  *
  ! *  I = 2*L*(J+1/2) + J + MUE + 1                                   *
  ! *                                                                  *
  ! ********************************************************************
  implicit none

  ! Dummy arguments
  integer :: kappa, muem05
  integer :: ikapmue

  ! Local variables
  integer :: iabs
  integer :: jp05, l

  jp05 = iabs(kappa)

  if (kappa<0) then
    l = -kappa - 1
  else
    l = kappa
  end if

  ikapmue = 2*l*jp05 + jp05 + muem05 + 1

end function ikapmue
