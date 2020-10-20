!------------------------------------------------------------------------------------
!> Summary: Broyden spin mixing for non-collinear magnetic moments
!> Author: Philipp Ruessmann
!> similar to mixbroydenspin of impurity code
!------------------------------------------------------------------------------------
module mod_mixbroydenspin

  private
  public :: spinmix_broyden

contains

  !-------------------------------------------------------------------------------
  !> Summary: Do Broyden spin mixing for noncollinear magnetic moment directions
  !> Author: Philipp Ruessmann
  !> Category: KKRhost 
  !> Deprecated: False 
  !> Perform spin mixing of nonco angles by first creating magentization vectors
  !> that are then mixed (should be more stable than mixing the angles directly)
  !-------------------------------------------------------------------------------
  subroutine spinmix_broyden(iter, natyp, new_angles, iounit)
    use mod_datatypes, only: dp
    use mod_constants, only: pi
    use mod_wunfiles, only: t_params
    use global_variables, only: mixfac_broydenspin, ninit_broydenspin, memlen_broydenspin
    use mod_runoptions, only: write_angles_alliter
    use mod_broyden, only: broyden
    implicit none
    ! interface
    integer, intent(in) :: iter !! iteration counter
    integer, intent(in) :: natyp !! number of atoms (array dimensions)
    real (kind=dp), intent(in) :: new_angles(2, natyp) !! new directions of the nonco angles
    integer, intent(in) :: iounit !! output unit where the nonco_angles output file is written
    ! local
    integer :: nfixed !! number of fixed angles
    integer :: ipos !! for loop indices etc.
    integer :: vlen !! length of vector
    real (kind=dp), allocatable :: vector(:, :) !! magnetization directions (first index=1 is old values, =2 is new values)
    real (kind=dp) :: rms, moment(3), totmoment, totxymoment, theta, phi
    integer :: i1

    ! get number of fixed angles (we need to use the rest only for the mixing)
    nfixed = 0
    do i1 = 1, natyp
        if (t_params%fixdir(i1)) nfixed = nfixed + 1
    end do

    ! allocate working arrays
    vlen = 3*(natyp-nfixed)
    allocate(vector(vlen, 2), stat=ipos)
    if (ipos>0) stop 'Error allocating vector in spinmix_broyden'
    vector(:,:) = 0.0_dp

    ! convert angles to magnetization direction vectors 
    ipos = 0
    do i1 = 1, natyp
      if (.not. t_params%fixdir(i1)) then
        ! set old vector
        theta = t_params%theta(i1)
        phi = t_params%phi(i1)
        vector(1+ipos*3, 1) = cos(phi)*sin(theta)
        vector(2+ipos*3, 1) = sin(phi)*sin(theta)
        vector(3+ipos*3, 1) = cos(theta)
        ! new vector
        theta = new_angles(1, i1)
        phi = new_angles(2, i1)
        vector(1+ipos*3, 2) = cos(phi)*sin(theta)
        vector(2+ipos*3, 2) = sin(phi)*sin(theta)
        vector(3+ipos*3, 2) = cos(theta)
        ! update loop counter
        ipos = ipos + 1
      end if ! fixdir
    end do

    ! find rms value
    rms = 0.0_dp
    do ipos = 1, vlen
      rms = rms + (vector(ipos,2) - vector(ipos,1))**2
    end do
    rms = sqrt(rms)

    ! now do Broyden mixing
    call broyden (vector, vlen, mixfac_broydenspin, rms, iter, &
                  ninit_broydenspin, memlen_broydenspin, vlen)

    write(*,*) 'spinmix_broyden', iter, natyp-nfixed, rms
    
    ! output 
    ipos = 0
    do i1 = 1, natyp
      if (.not. t_params%fixdir(i1)) then
        ! transform back to angles (vector(:,2) now is the output direction vector)
        moment(1) = vector(1+3*ipos, 2)
        moment(2) = vector(2+3*ipos, 2)
        moment(3) = vector(3+3*ipos, 2)

        totmoment = sqrt(moment(1)**2+moment(2)**2+moment(3)**2)
        totxymoment = sqrt(moment(1)**2+moment(2)**2)

        ! theta not 0 or pi
        if (abs(totxymoment)>1d-05) then
          theta = acos(moment(3)/totmoment)
          phi = atan2(moment(2), moment(1))
        ! theta is 0 or pi
        else
          if (moment(3) < 0.0_dp) then
            theta = pi
          else
            theta = 0.0_dp
          end if
          phi = 0.0_dp
        end if

        ! now set mixed angles as output
        t_params%theta(i1) = theta
        t_params%phi(i1) = phi

        ! update loop counter
        ipos = ipos + 1
      end if ! fixdir
      ! write output nonco angles file
      write (iounit, *) t_params%theta(i1)/pi*180.0_dp, t_params%phi(i1)/pi*180.0_dp, t_params%fixdir(i1)
      if (write_angles_alliter) write (iounit+1, *) t_params%theta(i1)/pi*180.0_dp, t_params%phi(i1)/pi*180.0_dp, t_params%fixdir(i1)
    end do

    ! cleanup allocations of working arrays
    deallocate(vector, stat=ipos)

  end subroutine spinmix_broyden

end module mod_mixbroydenspin
