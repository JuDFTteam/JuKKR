module mod_cdnlzdz
  
  private
  public :: cdnlzdz

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates derivative of Neumann function
  !> Author: 
  !> Category: KKRhost, special-functions, dirac
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the derivative of the Neumann function 
  !> d n(L,Z) / dz
  !> analytically
  !-------------------------------------------------------------------------------
  function cdnlzdz(l, z, mode)
    use :: mod_datatypes, only: dp
    use :: mod_cnlz
    implicit none

    ! Dummy arguments
    integer :: l, mode
    complex (kind=dp) :: z
    complex (kind=dp) :: cdnlzdz

    if (mode==1) then

      if (l==0) then

        cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
      else
        cdnlzdz = (l*cnlz(l-1,z)-(l+1)*cnlz(l+1,z))/real(2*l+1, kind=dp)
        return
      end if
    else if (mode==2) then

      if (l==0) then
        cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
      else
        cdnlzdz = cnlz(l-1, z) - (l+1)*cnlz(l, z)/z
        return
      end if
    else
      cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
    end if
  end function cdnlzdz

end module mod_cdnlzdz
