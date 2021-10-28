!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Module storing the run options and the paramters for bfields and constraining fields
!> 
!> Author: Sascha Brinker, Eduardo Mendive, Nicolas Essing
!>   Written by Sascha Brinker, ported from KKRhost to KKRnano by the other two.
!------------------------------------------------------------------------------------
module mod_bfield
  
  use :: NonCollinearMagnetism_Helpers_mod, only : rotatematrix

  implicit none

  private
  public :: bfield_data, load_bfields_from_disk, init_bfield, add_bfield

  ! In KKRhost, there is a type named 'type_bfield',
  ! which contains the input parameters and the information from 'bfield_data'.
  ! Here, one instance of 'bfield_data' is stored for each atom, while
  ! in KKRhost the components of that type were arrays (respectively arrays with one
  ! dimension more than here) containing the information for all atoms.
  
  !-------------------------------------------------------------------------------
  !> Summary: A type storing information on magnetic fields for a single atom
  !-------------------------------------------------------------------------------
  type :: bfield_data
    double precision, dimension(3) :: bfield_ext    !! external magnetic field in cartesian coordinates
    double precision, dimension(3) :: bfield_constr !! constraining field in cartesian coordinates
    double precision, dimension(3) :: mag_torque    !! Magnetic torque
    double precision, dimension(3) :: mag_mom       !! Magnetic moment

    double precision, dimension(:,:,:), allocatable :: thetallmat !! shapefun in the ll' expansion
  end type

contains

  !> Load the external noncollinear magnetic field (if used) and the initial
  !> guess for the constraining fields from disk.
  !> This subroutine loads information for the locally treated atoms into a
  !> given array of structs of the right size.
  subroutine load_bfields_from_disk(bfields, external_bfield, verbosity, naez, &
                                    atom_ids, fix_angle_modes)
    type(bfield_data), intent(inout) :: bfields(:)        !! Array for the read fields
    logical,           intent(in)    :: external_bfield   !! Whether to apply external nonco fields
    integer,           intent(in)    :: verbosity         !! Output verbosity
    integer,                       intent(in) :: naez     !! Number of atoms in the unit cell
    integer(kind=4), dimension(:), intent(in) :: atom_ids !! Locally treated atoms
    integer(kind=1), dimension(:), intent(in) :: fix_angle_modes
    
    integer :: num_local_atoms ! number of atoms treated in this MPI rank
    integer :: ila             ! atom index
    type(bfield_data), dimension(:), allocatable :: all_bfields ! array of all bfields to read

    ! Get the number of locally treated atoms form the size of the output array
    num_local_atoms = size(bfields)

    ! Allocate an array for all atoms in the unit cell to read the files (which of course
    ! contain the fields of every atom)
    allocate(all_bfields(naez))

    ! Read the files
    call read_external_bfield(all_bfields, external_bfield, verbosity)
    call read_constr_bfield_initial_guess(all_bfields, verbosity, fix_angle_modes)

    ! Copy the locally treated atoms to the output
    do ila = 1, num_local_atoms
      bfields(ila) = all_bfields(atom_ids(ila))
    end do

    ! Forget the other atoms
    deallocate(all_bfields)
  end subroutine

  !> Initialize a bfield_data struct. Allocates and calculates the shapefunction
  !> used in other subroutines.
  subroutine init_bfield(bfield, lmax, npan_log, npan_eq, ipan_intervall, &
                         thetasnew, iend, icleb, cleb, ifunm)
    type(bfield_data), intent(inout) :: bfield !! The bfield data
    integer, intent(in) :: lmax     !! Angular momentum cutoff
    integer, intent(in) :: npan_log !! Chebyshev mesh resolution
    integer, intent(in) :: npan_eq  !! Chebyshev mesh resolution
    integer, dimension(:), intent(in) :: ipan_intervall !! Indices for important radial points in the mesh
    double precision, dimension(:,:), intent(in)  :: thetasnew !! Interpolated shape function in Chebychev mesh
    integer,                          intent (in) :: iend      !! Number of nonzero gaunt coefficients
    integer,          dimension(:,:), intent (in) :: icleb     !! Mapping from array index to angular momentum indices
    double precision, dimension(:),   intent (in) :: cleb      !! gaunt coefficients
    integer, dimension(:),            intent (in) :: ifunm     !! Switch for use of gaunt coefficients
    
    integer :: ncleb, irmdnew
    integer :: i_stat

    ncleb = size(cleb)
    irmdnew = size(thetasnew, 1) ! npan_tot*(ncheb+1)

    ! Calculate the LL' expansion of the shape function in the new mesh which
    ! is needed to convolute the magnetic field (only done once and stored to
    ! safe computing time)
    allocate (bfield%thetallmat((lmax+1)**2,(lmax+1)**2,irmdnew), stat=i_stat)
    call calc_thetallmat(bfield%thetallmat, lmax, ipan_intervall(npan_log+npan_eq) + 1, &
                          iend, irmdnew, thetasnew, ifunm, icleb, cleb)

  end subroutine init_bfield

  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise initial guess for the constraining magnetic
  !> field from bconstr.dat
  !> Author: MdSD, Nicolas Essing
  !>
  !> The file should contain three floating point values in each line and a line
  !> for each atom. The first line is treated as a header and skipped.
  !> Intepreted as initial constraining bfield in cartesian coordinates, in Ry.
  !-------------------------------------------------------------------------------
  subroutine read_constr_bfield_initial_guess(bfields, verbosity, fix_angle_modes)
    type(bfield_data), dimension(:), intent(inout) :: bfields
    integer, intent(in) :: verbosity
    integer(kind=1), dimension(:), intent(in) :: fix_angle_modes
    
    integer :: number_of_atoms
    integer :: iatom
    integer :: iostat
    logical :: file_exists

    number_of_atoms = size(bfields)

    inquire(file='bconstr_in.dat', exist=file_exists)

    if (file_exists) then
      open(unit=57493215, file='bconstr_in.dat', iostat=iostat)
      read(57493215, *)  ! skip header
      do iatom = 1, number_of_atoms
        read(57493215, *, iostat=iostat) bfields(iatom)%bfield_constr(:)
        if (iostat /= 0) then
          write(*,*) "Error reading bconstr_in.dat"
          stop
        end if

        ! In case constraint magnetism is not used for this atom, set the initial
        ! guess to zero. Will not be updated in that case and can be added to the
        ! potential without further checking on input parameters.
        ! If an initial guess was provided that is not zero, give a warning.
        if (.not. (fix_angle_modes(iatom) == 2 .or. fix_angle_modes(iatom) == 3)) then
          if (any(bfields(iatom)%bfield_constr /= 0) .and. verbosity >= 0) then
            write(*,'(2A,I3,2A)') 'Warning: Initial guess for constraint magnetic field ', &
                    'for atom ', iatom, ' was not zero, but no constraint magnetism is ', &
                    'used for this atom. Will be set to zero.'
          end if
          bfields(iatom)%bfield_constr(:) = 0
        end if
      end do
      close(57493215)
    else
      ! No 'bconstr_in.dat' given, use default 0
      do iatom = 1, number_of_atoms
        bfields(iatom)%bfield_constr(:) = 0.
      end do
    end if

    if (verbosity >= 3) then
      ! Write detailed information
      write(*,'(79("#"))')
      write(*,'(2X,A)') 'Initial guess for constraining fields'
      write(*,'(79("#"))')
      write(*,'(2X,A4,3(7X,A7,3X))') 'atom','Bx [Ry]','By [Ry]','Bz [Ry]'
      do iatom = 1, number_of_atoms
        write(*,'(2X,I4,3(2X,E15.8))') iatom, bfields(iatom)%bfield_constr(:)
      end do
      write(*,'(79("#"))')
    else if (verbosity >= 2) then
      ! Give an overview
      if (file_exists) then
        write(*,'(A)') 'Initial constraining magnetic fields read from file.'
      else
        write(*,'(A)') 'Initial constraining magnetic fields set to zero.'
      end if
    else
      ! No output
    end if
  end subroutine


  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise external magnetic field from bfield.dat,
  !> if it should be used according to the input parameters.
  !> Author: Nicolas Essing
  !> 
  !> If external fields shall not be used, the external field is set to zero and
  !> an eventually present file is ignored.
  !> If external fields shall be used and no file is present, the program stops
  !> with an error.
  !> The file should contain three numbers per line and one line per atom.
  !> Read as 'theta  phi  bfield_strength' with the angles in degrees and the
  !> strength in Ry.
  !> Lines can be commented out with a # as first character.
  !>------------------------------------------------------------------------------
  subroutine read_external_bfield(bfields, external_bfield, verbosity)
    type(bfield_data), dimension(:), intent(inout) :: bfields
    logical, intent(in) :: external_bfield
    integer, intent(in) :: verbosity

    integer        :: number_of_atoms
    integer        :: iatom, iostat
    character(256) :: linebuffer
    double precision, dimension(:), allocatable :: phi, theta, strength

    number_of_atoms = size(bfields)

    allocate(phi(number_of_atoms), theta(number_of_atoms), strength(number_of_atoms))

    if (external_bfield) then
      open(unit=57493215, file='bfield.dat', iostat=iostat)
      if (iostat /= 0) then
        write(*,*) "I/O Error opening bfield.dat"
        write(*,*) "If you do not want external nonco bfields, turn it of in the input.conf"
        stop
      end if

      ! manual loop over iatom because comments are allowed
      iatom = 1
      do while (iatom <= number_of_atoms)
        read(57493215, '(A)', iostat=iostat) linebuffer
        if (iostat /= 0) then
          write(*,*) "I/O Error while reading bfield.dat"
          stop
        end if
        if (linebuffer(1:1) == '#') cycle ! input line commented out
        read(linebuffer, *, iostat=iostat) theta(iatom), phi(iatom), strength(iatom)
        if (iostat /= 0) then
          write(*,*) "Error parsing a line in bfield.dat"
          stop
        end if
        iatom = iatom + 1
      end do
      close(57493215)
    else
      ! No 'bfield.dat' given, use default 0
      theta(:)    = 0.
      phi(:)      = 0.
      strength(:) = 0.
    end if

    ! Output before converting to radians and to carthesian coordinates.
    if (verbosity >= 3) then
      ! Write detailed information
      write(*,'(79("#"))')
      write(*,'(16X,A)') 'external non-collinear magnetic fields'
      write(*,'(79("#"))')
      write(*,'(2X,A4,4X,A10,6X,A8,6X,A11)') 'atom', 'theta [°]', 'phi [°]', 'bfield [Ry]'
      do iatom = 1, number_of_atoms
        write(*,'(2X,I4,2(2X,F12.8),2X,E15.8)') iatom, theta(iatom), phi(iatom), strength(iatom)
      end do
      write(*,'(79("#"))')
    else if (verbosity >= 2) then
      ! Give an overview
      if (external_bfield) then
        write(*,'(A)') 'External non-collinear magnetic fields loaded from file.'
        write(*,'(2X,A,E16.8)') 'Mean magnitude [Ry]:', sum(abs(strength)) / number_of_atoms
      else
        write(*,'(A)') 'Not applying external non-collinear magnetic fields.'
      end if
    else
      ! No output
    end if

    theta(:) = theta(:) / 360.0d0 * 8.d0 * datan(1.d0)
    phi(:)   = phi(:)   / 360.0d0 * 8.d0 * datan(1.d0)
    do iatom = 1, number_of_atoms
      bfields(iatom)%bfield_ext(1) = strength(iatom) * sin(theta(iatom)) * cos(phi(iatom))
      bfields(iatom)%bfield_ext(2) = strength(iatom) * sin(theta(iatom)) * sin(phi(iatom))
      bfields(iatom)%bfield_ext(3) = strength(iatom) * cos(theta(iatom))
    end do
  end subroutine


  !>------------------------------------------------------------------------------
  !> Summary: Adds the magnetic field to the potential
  !> Author: Sascha Brinker, Nicolas Essing
  !> 
  !> The field is added to the potential in LL' expansion. Both the external and
  !> the constraint field are added, if they are activated.
  !> The potential is updated as H = H - sigma * B with sigma the vector of
  !> pauli matrices and B the combined bfield.
  !>------------------------------------------------------------------------------
  subroutine add_bfield(bfield, vnspll, theta, phi, imt, iteration_number, &
                        itscf0, itscf1, lbfield_trans, &
                        lbfield_mt, transpose_bfield)
    type(bfield_data), intent(in) :: bfield
    double complex, dimension(:,:,:), intent(inout) :: vnspll ! The potential to add to
    double precision,  intent(in) :: theta, phi ! angles of the magnetic moment, not to be confused with theta and phi in bfield
    integer, intent(in) :: imt ! MT radius (index in cheb mesh)
    integer, intent(in) :: iteration_number !TODO this, or just a logical and do the check outside?
    integer, intent(in) :: itscf0, itscf1   !TODO ^
    logical, intent(in) :: lbfield_trans ! Apply only transveral
    logical, intent(in) :: lbfield_mt ! Apply only up do MT radius
    logical, intent(in) :: transpose_bfield ! Transpose the bfield (for left solutions)

    double complex, parameter :: cplx_i = (0.d0, 1.d0)
    integer :: lmmax, irmd, iend, ir ! loop boundaries and indices
    double precision, dimension(3) :: combined_bfields, dir ! vector of combined bfields and unit vector of magnetic moment direction
    double complex, dimension(2,2) :: bfield_mat ! bfield times pauli matrices
    double complex :: temp ! used to transpose the matrix

    ! If the current iteration is not in the window the magnetic fields should be
    ! applied, return without changing the potential
    if (iteration_number < itscf0 .or. iteration_number > itscf1) then
      return
    end if

    lmmax = size(bfield%thetallmat, 1) ! size(vnspll, 1) is 2*lmmax
    irmd = size(vnspll, 3)

    ! Add external and constraint field. If one of them is tured off by input
    ! parameters or mode, it is zero, so this distinction does not have to be
    ! done here.
    combined_bfields(:) = bfield%bfield_ext(:) + bfield%bfield_constr(:)

    if (lbfield_trans) then
      dir(1) = sin(theta) * cos(phi)
      dir(3) = sin(theta) * sin(phi)
      dir(2) = cos(theta)
      combined_bfields = combined_bfields - dir * dot_product(dir, combined_bfields)
    end if

    ! Fill potential matrix from bfield vector
    bfield_mat(1,1) = - combined_bfields(3)
    bfield_mat(1,2) = combined_bfields(1) + cplx_i * combined_bfields(2)
    bfield_mat(2,1) = combined_bfields(1) - cplx_i * combined_bfields(2)
    bfield_mat(2,2) = combined_bfields(3)

    ! Rotate to local frame
    call rotatematrix(bfield_mat, theta, phi, 1, 1)

    ! For the left solutions, transpose the bfield
    if (transpose_bfield) then
      temp = bfield_mat(1,2)
      bfield_mat(1,2) = bfield_mat(2,1)
      bfield_mat(2,1) = temp
    end if

    ! Define the loop boundary for the next loop. If the fields are only applied
    ! in the muffin tin region, integrate to that index, otherwise all the way out
    if (lbfield_mt) then
      iend = imt
    else
      iend = irmd
    end if

    ! Add the bfield to the potential
    do ir = 1, iend
      vnspll(1:lmmax,1:lmmax,ir) = vnspll(1:lmmax,1:lmmax,ir) - bfield_mat(1,1) * bfield%thetallmat(:,:,ir)
      vnspll(1:lmmax,lmmax+1:2*lmmax,ir) = vnspll(1:lmmax,lmmax+1:2*lmmax,ir) - bfield_mat(1,2) * bfield%thetallmat(:,:,ir)
      vnspll(lmmax+1:2*lmmax,1:lmmax,ir) = vnspll(lmmax+1:2*lmmax,1:lmmax,ir) - bfield_mat(2,1) * bfield%thetallmat(:,:,ir)
      vnspll(lmmax+1:2*lmmax,lmmax+1:2*lmmax,ir) = vnspll(lmmax+1:2*lmmax,lmmax+1:2*lmmax,ir) - bfield_mat(2,2) * bfield%thetallmat(:,:,ir)
    end do

  end subroutine

  !------------------------------------------------------------------------------------
  !> Summary: Shape function LL' expansion
  !> Author: Sascha Brinker
  !> Calculates the LL' expansion of the shape function similarly to vllmat_new.f90
  !> @note The input shapefunction (single L) uses pointer arrays for the lm index.
  !> The output does not need pointers!
  !> @endnote
  !------------------------------------------------------------------------------------
  subroutine calc_thetallmat(thetansll, lmax, imt1, iend, irmdnew, thetasnew, ifunm, icleb, cleb)
    double precision, dimension(:,:,:), intent (out) :: thetansll !! LL' expansion of the shapefunction
    integer,                            intent (in)  :: lmax      !! Angular momentum cut-off
    integer,                            intent (in)  :: imt1      !! index muffin-tin radius
    integer,                            intent (in)  :: iend      !! Number of nonzero gaunt coefficients
    integer,          dimension(:,:),   intent (in)  :: icleb     !! l and m values for the Gaunt coefficients
    double precision, dimension(:),     intent (in)  :: cleb      !! Gaunt coefficients
    double precision, dimension(:, :),  intent (in)  :: thetasnew !! shapefun on the Cheby mesh
    integer         , dimension(:),     intent (in)  :: ifunm     !! pointer array for shapefun

    double precision, parameter :: rfpi = 3.5449077018110318 ! sqrt(4*pi)
    double precision, parameter :: c0ll = 1.d0 / rfpi
    integer :: irmdnew       ! number of radials point on the Chebyshev mesh
    integer :: lmmax         ! number of angular momentum entries
    integer :: lmmax2
    integer :: ifun
    integer :: lm1, lm2, lm3 ! angular momentum indices
    integer :: ir            ! radial index
    integer :: j             ! index for loop over Gaunt coeffs
    double precision, dimension(:,:), allocatable :: shapefun_mod

    lmmax   = (lmax+1)**2
    lmmax2  = (2*lmax+1)**2
    irmdnew = size(thetasnew, 1) ! npan_tot*(ncheb+1)

    allocate(shapefun_mod(irmdnew, lmmax2))

    ! Build the shapefun array. Start with all zero. Inside muffin tin only l=0
    ! component is /= 0 (and constant), copy l=0 component for r outside of mt
    ! from thetasnew.
    shapefun_mod(:,:)        = 0.d0
    shapefun_mod(1:imt1,1)   = rfpi ! is multipled by C_LL^0
    shapefun_mod(imt1+1:, 1) = thetasnew(imt1+1:irmdnew,1)

    ! Copy other components from thetasnew. Convert from pointer indices to
    ! normal (l,m)-index
    do lm1 = 2, lmmax2
      ifun = ifunm(lm1)
      if(ifun /= 0) then !shapefun%lmused(lm1)==1) then
        shapefun_mod(imt1+1:, lm1) = thetasnew(imt1+1:, ifun)
      end if
    end do

    ! Initialize result
    thetansll(:,:,:) = 0.d0

    ! Diagonal part (not contained in gaunt-coeff)
    do lm1 = 1, lmmax
      thetansll(lm1,lm1,:) = shapefun_mod(:,1) * c0ll
    end do

    ! Offdiagonal part. This is effectively a loop over angular momentum indices
    ! lm1,lm2,lm3.  Iterate instead over the flattened array cleb
    ! containing the Gaunt coefficients and extract the angular momentum
    ! indices for each j.
    do j = 1, iend
      lm1 = icleb(j, 1)
      lm2 = icleb(j, 2)
      lm3 = icleb(j, 3)
      if (lm1 > lmmax .or. lm2 > lmmax .or. lm3 > lmmax2) then
        ! I think this should not happen, as icleb should only map to valid values
        ! for lm1,lm2,lm3.
        !TODO check that and remove this warning
        write(*,*) "Warning from calc_thetallmat in bfield.F90: Maybe invalid values in icleb:"
        write(*,*) j, lm1, lm2, lm3, lmmax, lmmax2
        cycle ! Gaunt coeffs zero
      end if
      if (ifunm(lm3) == 0) then
        cycle ! Some symmetry related to shapefunction, I think
      end if
      do ir = 1, irmdnew
        thetansll(lm1,lm2,ir) = thetansll(lm1,lm2,ir) + cleb(j) * shapefun_mod(ir,lm3)
        thetansll(lm2,lm1,ir) = thetansll(lm2,lm1,ir) + cleb(j) * shapefun_mod(ir,lm3)
      end do
    end do

  end subroutine calc_thetallmat

end module mod_bfield
