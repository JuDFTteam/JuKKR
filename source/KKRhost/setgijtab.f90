!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper for the calculation of the elements of the Greens function \(G_{ij}\)
!> Author: 
!> rapper for the calculation of the elements of the Greens function \(G_{ij}\)
!> will calculate the pairs depending on wheter one wishes to calculate the conductivity
!> of the pair-wise Heisenberg exchange interactions.
!------------------------------------------------------------------------------------
module mod_setgijtab
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for the calculation of the elements of the Greens function \(G_{ij}\)
  !> Author: 
  !> Category: input-output, reference-system, structural-greensfunction, KKRhost
  !> Deprecated: False 
  !> Wrapper for the calculation of the elements of the Greens function \(G_{ij}\)
  !> will calculate the pairs depending on wheter one wishes to calculate the conductivity
  !> of the pair-wise Heisenberg exchange interactions.
  !-------------------------------------------------------------------------------
  subroutine setgijtab(linterface,icc,naez,iqat,rbasis,bravais,natomimp,atomimp,    &
    rclsimp,ijtabcalc,iofgij,jofgij,nqcalc,iqcalc,natomimpd,ijtabcalc_i)
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
    use :: mod_gijcond
    use :: mod_gijxcpl
    use :: mod_runoptions, only: calc_exchange_couplings, use_cond_LB
    implicit none

    ! Scalar arguments
    integer :: icc, naez, natomimp, natomimpd, nqcalc
    logical :: linterface

    ! Array arguments
    integer :: atomimp(*), ijtabcalc(*), ijtabcalc_i(*), iofgij(*), iqat(*), iqcalc(*), jofgij(*)
    real (kind=dp) :: bravais(3, 3), rbasis(3, *), rclsimp(3, *)

    ! Local scalars
    integer :: i, ido, ii, j, jj, nn, nofgij
    ! external funcitons
    logical, external :: opt

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, '(79("="),/,15X,A)') 'SETGIJTAB: setting task-specific Gij pairs'
    write (1337, '(79("="),/)')
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    ido = 0
    ! ======================================================================
    if (use_cond_LB) call gijcond(ido, naez, rbasis, iqat, natomimp, rclsimp, atomimp, ijtabcalc, natomimpd)
    ! ======================================================================
    if (calc_exchange_couplings) call gijxcpl(ido, naez, rbasis, bravais, linterface, nqcalc, iqcalc, natomimp, rclsimp, atomimp, ijtabcalc, ijtabcalc_i, natomimpd)
    ! ======================================================================
    if (ido==0) then
      icc = 0
      write (6, 110)
      return
    end if
    ! ======================================================================
    nofgij = 0
    do i = 1, natomimp
      nn = (i-1)*natomimp
      do j = 1, natomimp
        if (ijtabcalc(nn+j)>0) then
          nofgij = nofgij + 1
          if (nofgij>natomimpd*natomimpd) then
            write (6, 100) 'NATOMIMPD', nofgij/natomimp
            stop
          end if
          iofgij(nofgij) = i
          jofgij(nofgij) = j
        end if
      end do
    end do
    if (nofgij==0) then
      icc = 0
      write (6, 110)
      return
    end if

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 120) natomimp, nofgij
    write (1337, 130)
    write (1337, 140)
    write (1337, 130)
    do i = 1, nofgij
      ii = iofgij(i)
      jj = jofgij(i)
      write (1337, 150) i, ii, atomimp(ii), (rclsimp(j,ii), j=1, 3), jj, atomimp(jj), (rclsimp(j,jj), j=1, 3)
    end do
    write (1337, 130)
    write (1337, *)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

100 format (6x, 'brahim ERROR: please increase the global parameter', /, 6x, a, ' to a value >=', i5, /)
110 format (6x, 'WARNING: Subroutine entered with invalid task ', 'specification', /, 6x, '         ICC will be set to 0 - no Gij calculated - ', 'input check? ', /)
120 format (6x, 'Number of different sites (NATOMIMP) :', i4, /, 6x, 'Number of pairs set       (NOFGIJ)   :', i4)
130 format (8x, 71('-'))
140 format (9x, 'pair|', ' I  IQ           position', 9x, 'J  JQ           position')
150 format (9x, i3, ' |', 2(i3,1x), 3f8.4, 1x, 2(i3,1x), 3f8.4)
160 format (i5, 2(i5,1x), 3f10.6, 1x, 2(i5,1x), 3f10.6)
  end subroutine setgijtab

end module mod_setgijtab
