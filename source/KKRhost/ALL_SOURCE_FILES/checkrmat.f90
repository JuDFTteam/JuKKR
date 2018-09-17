module mod_checkrmat
  use :: mod_datatypes, only: dp
  private :: dp

contains

  function checkrmat(rmat, co1, si1, co2, si2, co3, si3, i, j)
    ! ********************************************************************
    ! *                                                                  *
    ! *  check whether the values of the cosinus and sinus found for the *
    ! *  Euler angles TET1, TET2, TET3 are consistent with the           *
    ! *  rotation matrix   RMAT                                          *
    ! *                                                                  *
    ! ********************************************************************

    implicit none

    ! Dummy arguments
    real (kind=dp) :: co1, co2, co3, si1, si2, si3
    integer :: i, j
    logical :: checkrmat
    real (kind=dp) :: rmat(3, 3)

    ! Local variables
    real (kind=dp) :: a, b
    logical :: equal
    logical :: result

    equal(a, b) = (abs(a-b)<1e-7_dp)

    result = .false.

    if (i==1) then
      if (j==1) then
        result = equal(rmat(1,1), co3*co2*co1-si3*si1)
      else if (j==2) then
        result = equal(rmat(1,2), co3*co2*si1+si3*co1)
      else if (j==3) then
        result = equal(rmat(1,3), -co3*si2)
      end if
    else if (i==2) then
      if (j==1) then
        result = equal(rmat(2,1), -si3*co2*co1-co3*si1)
      else if (j==2) then
        result = equal(rmat(2,2), -si3*co2*si1+co3*co1)
      else if (j==3) then
        result = equal(rmat(2,3), si3*si2)
      end if
    else if (j==1) then
      result = equal(rmat(3,1), si2*co1)
    else if (j==2) then
      result = equal(rmat(3,2), si2*si1)
    else if (j==3) then
      result = equal(rmat(3,3), co2)
    end if
    checkrmat = result
  end function checkrmat

end module mod_checkrmat
