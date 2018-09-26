module mod_gijcond
  
  private
  public :: gijcond

contains

  !-------------------------------------------------------------------------------
  !> Summary: Construct structure of G_ij structure needed for conductivity calculation
  !> Author: 
  !> Category: KKRhost, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> In case of tasks requiring Gij blocks calculation, set variables: 
  !>
  !> NATOMIMP, RCLSIMP(3,1..NATOMIMP), ATOMIMP(1..NATOMIMP)
  !> IJTABCALC flag to which pair is needed: I,J --> (I-1)*NATOMIMP + J
  !> IDO takes on the value 1 or 0 if setting up process was OK or not 
  !>
  !> CONDUCTANCE calculation case
  !>             still to be implemente the correct read in
  !-------------------------------------------------------------------------------
  subroutine gijcond(ido, naez, rbasis, iqat, natomimp, rclsimp, atomimp, ijtabcalc, natomimpd)
    use :: mod_datatypes, only: dp
    implicit none

    ! Parameters
    integer :: ncpaird
    parameter (ncpaird=10)

    ! Arguments
    integer :: ido, naez, natomimp, natomimpd
    integer :: atomimp(*), ijtabcalc(*), iqat(*)
    real (kind=dp) :: rbasis(3, *), rclsimp(3, *)

    ! Locals
    integer :: i, iat, iatcondl(ncpaird), iatcondr(ncpaird), j, jat, ncondpair, nn

    ido = 0
    do i = 1, ncpaird
      iatcondl(i) = 0
      iatcondr(i) = 0
    end do

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 100)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    ! ---------------------------------------------------- dummy
    ! settings so far, need to be replaced by conductance input
    ! and some output
    ncondpair = 4
    if (ncondpair>ncpaird) then
      write (6, 110) 'local', 'NCPAIRD', ncondpair
      stop
    end if
    iatcondl(1) = 1
    iatcondr(1) = 2
    iatcondl(2) = 1
    iatcondr(2) = 2
    iatcondl(3) = 2
    iatcondr(3) = 1
    iatcondl(4) = 2
    iatcondr(4) = 1

    ! ---------------------------------------------------- dummy
    if (ncondpair==0) return
    do i = 1, ncondpair
      if ((iatcondl(i)<=0) .or. (iatcondl(i)>naez)) return
      if ((iatcondr(i)<=0) .or. (iatcondr(i)>naez)) return
    end do

    natomimp = 2*ncondpair
    if (natomimp>natomimpd) then
      write (6, 110) 'global', 'NATOMIMPD', natomimp
      stop
    end if

    do i = 1, natomimp
      nn = (i-1)*natomimp
      do j = 1, natomimp
        ijtabcalc(nn+j) = 0
      end do
    end do

    nn = 0
    do i = 1, ncondpair
      iat = iqat(iatcondl(i))      ! left lead
      nn = nn + 1
      do j = 1, 3
        rclsimp(j, nn) = rbasis(j, iat)
      end do
      atomimp(nn) = iat
      iat = nn

      jat = iqat(iatcondr(i))      ! right lead
      nn = nn + 1
      do j = 1, 3
        rclsimp(j, nn) = rbasis(j, jat)
      end do
      atomimp(nn) = jat
      jat = nn
      ijtabcalc((iat-1)*natomimp+jat) = 1
    end do
    if (natomimp/=nn) then
      write (6, '(6X,A,/,6X,A,/)') 'ERROR: Found some inconsistencies in IATCOND arrays', '       Please check your CONDUCTANCE input'
      stop
    end if
    ido = 1
100 format (5x, '< GIJCOND > : Conductance/conductivity calculation', /)
110 format (6x, 'Dimension ERROR: please increase the ', a, ' parameter', /, 6x, a, ' to a value >=', i5, /)
  end subroutine gijcond

end module mod_gijcond
