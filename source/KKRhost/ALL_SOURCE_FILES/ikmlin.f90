subroutine ikmlin(iprint, nsollm, ikm1lin, ikm2lin, nlmax, nmuemax, linmax, &
  nl)
!   ********************************************************************
!   *                                                                  *
!   * SETUP TABLE OF INDICES    IKM(INT)                               *
!   *                                                                  *
!   *  IKM IS STANDARD INDEX IN  (KAPPA,MUE)-REPRESENTATION            *
!   *  IKM = 2*L*(J+1/2) + J + MUE + 1                                 *
!   *                                                                  *
!   *  INT NUMBERS LINEARLY ONLY NON-VANISHING ELEMENTS OF M-SS        *
!   *  USED TO CALCULATE DOS ...                                       *
!   *                                                                  *
!   ********************************************************************
  use :: mod_types, only: t_inc
  implicit none


! Dummy arguments
  integer :: iprint, linmax, nl, nlmax, nmuemax
  integer :: ikm1lin(linmax), ikm2lin(linmax), nsollm(nlmax, nmuemax)

! Local variables
  integer :: i, il, imue, k1, k2, kap(2), l, lin, muem05, nsol
  integer :: ikapmue

  lin = 0

  do il = 1, nl
    l = il - 1
    muem05 = -il - 1
    kap(1) = -l - 1
    kap(2) = +l

    do imue = 1, 2*il
      muem05 = muem05 + 1
      nsol = nsollm(il, imue)

      do k2 = 1, nsol
        do k1 = 1, nsol
          lin = lin + 1
          ikm1lin(lin) = ikapmue(kap(k1), muem05)
          ikm2lin(lin) = ikapmue(kap(k2), muem05)
        end do
      end do

    end do
  end do

  if (iprint<2) return
  if (t_inc%i_write>0) then
    write (1337, fmt='('' INT='',I3,''  IKM=('',I3,'','',I3,'')'')')(i, &
      ikm1lin(i), ikm2lin(i), i=1, lin)
  end if
end subroutine
