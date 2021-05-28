!------------------------------------------------------------------------------------
!> Summary:Spin mixing for non-collinear magnetic moments
!> Author: Philipp Ruessmann
!> See also mixbroydenspin of impurity code for reference to broyden spin mixing
!------------------------------------------------------------------------------------
module mod_mixnocospin

  private
  public :: spinmix_noco

contains


  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for spin miging of nonco angles
  !> Author: Philipp Ruessmann
  !> Category: KKRhost 
  !> Deprecated: False 
  !> Perform spin mixing of nonco angles either with setting output angles to new input
  !> (this is the default behavior) or do broyden mixing for the nonco moments to get
  !> the new directions.
  !-------------------------------------------------------------------------------
  subroutine spinmix_noco(iter, natyp, theta, phi, fixdir, angles_new, totmoment, iounit)
    use mod_datatypes, only: dp
    use mod_constants, only: pi
    use mod_runoptions, only: use_broyden_spinmix, write_angles_alliter, disable_print_serialnumber, fix_nonco_angles, decouple_spins_cheby
    use global_variables, only: qbound_spin, angles_cutoff  !! MdSD: criterion to fix angles (see also rhovalnew)
    use mod_version_info, only: version_print_header
    use mod_wunfiles, only: t_params
    implicit none
    ! interface
    integer, intent(in) :: iter !! iteration counter
    integer, intent(in) :: natyp !! number of atoms (array dimensions)
    real (kind=dp), dimension (natyp), intent(in) :: theta !! theta nonco angles in rad
    real (kind=dp), dimension (natyp), intent(in) :: phi !! phi nonco angles in rad
    logical, dimension(natyp), intent(in) :: fixdir  !! booleans to fix nocon angles
    real (kind=dp), intent(in) :: totmoment(natyp) !! length of the magentization vectors
    real (kind=dp), intent(in) :: angles_new(2, natyp) !! new directions of the nonco angles
    integer, intent(in) :: iounit !! output unit where the nonco_angles output file is written
    ! local
    integer :: i1 !! atom index
    real (kind=dp) :: diff_phi !! for difference in phi nonco angles
    real (kind=dp), dimension (natyp) :: diff_angles_sq !! squared difference in nonco angles in rad
    logical :: all_fixed !! True if all nonco angles are fixed


    ! MdSD,PR: write information on new angles to output file
    if (.not.decouple_spins_cheby) then

      write (1337,*)
      write (1337, '("      I1    In/Out THETA[deg]       In/Out PHI[deg]        FIXDIR[boolean]   RMS(angles)[deg]")')
      do i1 = 1, natyp
        if (.not.fixdir(i1)) then
          ! renormalize to difference in phi in range [0,180]
          diff_phi = phi(i1) - angles_new(2, i1)
          if (diff_phi>pi) diff_phi = diff_phi - pi
          ! now compute squared difference of the angles
          diff_angles_sq(i1) = (theta(i1)-angles_new(1, i1))**2 + diff_phi**2
        else
          diff_angles_sq(i1) = 0.0_dp
        end if
        write (1337, '(I8,4F12.6,3x,1L,17x,1F12.6)') i1, theta(i1)*180.0_dp/pi, angles_new(1, i1)/pi*180.0_dp, &
            phi(i1)*180.0_dp/pi, angles_new(2, i1)/pi*180.0_dp, fixdir(i1), sqrt(diff_angles_sq(i1))/pi*180.0_dp
        if (i1==1) all_fixed = .true. ! initialize
        if(.not.fixdir(i1)) all_fixed = .false. ! if at least one is not fixed we set this to False
      end do                     ! i1
      write (1337, '(A, 1F12.6)') "Total RMS(angles) [deg]:", sqrt(sum(diff_angles_sq))/pi*180.0_dp
      ! write also to std out
      write (*, '(A, 1F12.6)') "Total RMS(angles) [deg]:", sqrt(sum(diff_angles_sq))/pi*180.0_dp

      ! check if all angles should be fixed (see SPINMIXQBOUND input parameter)
      if ((.not.all_fixed) .and. sqrt(sum(diff_angles_sq))/pi*180.0_dp < qbound_spin) then
        write(*,*) 'TOTAL RMS(angles) < qbound_spin:', sqrt(sum(diff_angles_sq))/pi*180.0_dp, qbound_spin
        write(*,*) 'Fix all angles from now on'
        do i1 = 1, natyp
          t_params%fixdir(i1) = .true.
        end do
      end if

    end if ! .not.decouple_spins_cheby
  
    ! rewrite new theta and phi to nonco_angle_out.dat, nonco_angle.dat is the input
    if (.not. fix_nonco_angles) then

      ! open nocno_angles files
      open (unit=iounit, file='nonco_angle_out.dat', form='formatted')
      call version_print_header(iounit, disable_print=disable_print_serialnumber)
      if (write_angles_alliter) open(unit=iounit+1, file='nonco_angle_out_all_iter.dat', form='formatted', position='append')

      if (.not.use_broyden_spinmix) then

        ! conventional scheme: use output angles as input
        do i1 = 1, natyp
          ! save to file in converted units (degrees)
          if (t_params%fixdir(i1)) then
            ! keep the old angles
            write (iounit, *) t_params%theta(i1)/pi*180.0_dp, t_params%phi(i1)/pi*180.0_dp, t_params%fixdir(i1)
          else
            ! MdSD: don't change the angles if the spin moment is too small
            if (totmoment(i1) > angles_cutoff) then
              ! update angles
              write (iounit, *) angles_new(1, i1)/pi*180.0_dp, angles_new(2, i1)/pi*180.0_dp, t_params%fixdir(i1)
              ! use internal units here
              t_params%theta(i1) = angles_new(1, i1)
              t_params%phi(i1) = angles_new(2, i1)
            else
              write (iounit, *) t_params%theta(i1)/pi*180.0_dp, t_params%phi(i1)/pi*180.0_dp, t_params%fixdir(i1)
            end if
          end if
          if (write_angles_alliter) write (iounit+1, *) t_params%theta(i1)/pi*180.0_dp, t_params%phi(i1)/pi*180.0_dp, t_params%fixdir(i1)
        end do

      else ! use_broyden_spinmix

        ! use spin moxing as in impurity code
        call spinmix_broyden(iter, natyp, angles_new, totmoment, iounit) 

      end if ! use_broyden_spinmix

      ! close output noco_angles files
      close (iounit)
      if (write_angles_alliter) close(iounit+1)

    end if ! .not.fix_nonco_angles

  end subroutine spinmix_noco



  !-------------------------------------------------------------------------------
  !> Summary: Do Broyden spin mixing for noncollinear magnetic moment directions
  !> Author: Philipp Ruessmann
  !> Category: KKRhost 
  !> Deprecated: False 
  !> Perform spin mixing of nonco angles by first creating magentization vectors
  !> that are then mixed (should be more stable than mixing the angles directly)
  !-------------------------------------------------------------------------------
  subroutine spinmix_broyden(iter, natyp, angles_new, totmoment_atoms, iounit)
    use mod_datatypes, only: dp
    use mod_constants, only: pi
    use mod_wunfiles, only: t_params
    use global_variables, only: mixfac_broydenspin, ninit_broydenspin, memlen_broydenspin, angles_cutoff  !! MdSD: criterion to fix angles (see also rhovalnew)
    use mod_runoptions, only: write_angles_alliter
    use mod_broyden, only: broyden
    implicit none
    ! interface
    integer, intent(in) :: iter !! iteration counter
    integer, intent(in) :: natyp !! number of atoms (array dimensions)
    real (kind=dp), intent(in) :: angles_new(2, natyp) !! new directions of the nonco angles
    real (kind=dp), intent(in) :: totmoment_atoms(natyp) !! length of the magentization vectors
    integer, intent(in) :: iounit !! output unit where the nonco_angles output file is written
    ! local
    real (kind=dp), parameter :: tol = 1.0e-12_dp !! MdSD: shifting problems with Broyden to the future
    integer :: nfixed !! number of fixed angles
    integer :: ipos !! for loop indices etc.
    integer :: vlen !! length of vector
    real (kind=dp), allocatable :: vector(:, :) !! magnetization directions (first index=1 is old values, =2 is new values)
    real (kind=dp) :: rms, moment(3), totmoment, totxymoment, theta, phi, alpha
    integer :: i1

    ! get number of fixed angles (we need to use the rest only for the mixing)
    nfixed = 0
    do i1 = 1, natyp
        if (t_params%fixdir(i1) .or. totmoment_atoms(i1) < angles_cutoff) nfixed = nfixed + 1
    end do

    ! allocate working arrays
    vlen = 3*(natyp-nfixed)
    allocate(vector(vlen, 2), stat=ipos)
    if (ipos>0) stop 'Error allocating vector in spinmix_broyden'
    vector(:,:) = 0.0_dp

    ! convert angles to magnetization direction vectors 
    ipos = 0
    do i1 = 1, natyp
      if (.not. (t_params%fixdir(i1) .or. totmoment_atoms(i1) < angles_cutoff)) then
        ! we use the same vector length (i.e. the output length since we can assume that it will not change much)
        ! set old vector
        theta = t_params%theta(i1)
        phi = t_params%phi(i1)
        vector(1+ipos*3, 1) = totmoment_atoms(i1)*cos(phi)*sin(theta)
        vector(2+ipos*3, 1) = totmoment_atoms(i1)*sin(phi)*sin(theta)
        vector(3+ipos*3, 1) = totmoment_atoms(i1)*cos(theta)
        ! new vector
        theta = angles_new(1, i1)
        phi = angles_new(2, i1)
        vector(1+ipos*3, 2) = totmoment_atoms(i1)*cos(phi)*sin(theta)
        vector(2+ipos*3, 2) = totmoment_atoms(i1)*sin(phi)*sin(theta)
        vector(3+ipos*3, 2) = totmoment_atoms(i1)*cos(theta)
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
    write (1337,'("spinmix_broyden: iter=",i8,"  ninit=",i8,"  rms=",es16.8)') iter, ninit_broydenspin, rms

    ! different alpha for simple mixing and broyden mixing steps
    alpha = mixfac_broydenspin
    if (iter<=ninit_broydenspin) alpha = 1.0_dp ! always use alpha=1 for simple mixing steps


    ! MdSD: there are special high-symmetry situations where the rms for the angles is zero
    ! MdSD: linear mixing is fine with that, but if broyden is called it will cause a NaN
    ! MdSD: this line delays using Broyden
    if (rms < tol .and. iter == ninit_broydenspin) ninit_broydenspin = ninit_broydenspin + 1

    ! now do Broyden mixing
    call broyden (vector, vlen, alpha, rms, iter, &
                  ninit_broydenspin, memlen_broydenspin, vlen)

    ! output 
    ipos = 0
    do i1 = 1, natyp
      if (.not. (t_params%fixdir(i1) .or. totmoment_atoms(i1) < angles_cutoff)) then
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
          if (moment(3) < 0.0_dp .and. abs(moment(3)) > 1e-14_dp) then
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

end module mod_mixnocospin
