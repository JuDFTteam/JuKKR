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
  
  use :: NonCollinearMagnetism_mod, only : rotatematrix

  implicit none

  ! In KKRhost, there is a type named 'type_bfield',
  ! which contains the input parameters and the information from 'bfield_data'.
  ! Here, one instance of 'bfield_data' is stored for each atom, while
  ! in KKRhost the components of that type were arrays (respectively arrays with one
  ! dimension more than here) containing the information for all atoms.
  
  !-------------------------------------------------------------------------------
  !> Summary: A type storing information on magnetic fields for a single atom
  !-------------------------------------------------------------------------------
  type :: bfield_data
    double precision, dimension(3) :: bfield !! external magnetic field in cartesian coordinates
    double precision               :: bfield_strength !! absolute value of the external magnetic field, dimensions 
    double precision, dimension(3) :: bfield_constr !! constraining field in cartesian coordinates
    double precision               :: theta !! polar angle of the magnetic field
    double precision               :: phi   !! azimuthal angle of the magnetic field
    double precision, dimension(:,:,:), allocatable :: thetallmat !! shapefun in the ll' expansion
    double precision, dimension(3) :: mag_torque !! Magnetic torque 
  end type

contains

  subroutine load_bfields_from_disk(bfields, lbfield, lbfield_constr)
    type(bfield_data), allocatable :: bfields(:)
    logical, intent(in) :: lbfield
    logical, intent(in) :: lbfield_constr
    
    integer :: number_of_atoms

    number_of_atoms = size(bfields)
    
    if (lbfield) then
      call read_bfield(bfields, number_of_atoms)
    end if
    if (lbfield_constr) then
      call read_bconstr(bfields, number_of_atoms)
    end if
    
  end subroutine

  subroutine init_bfield(lbfield, lbfield_constr, bfield, lmax, &
                         npan_log, npan_eq, ipan_intervall, thetasnew, &
                         iend, icleb, cleb, ifunm)
    logical, intent(in) :: lbfield        !! Whether to use external nonco bfields
    logical, intent(in) :: lbfield_constr !! Whether to use constraint bfields
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
    ! TODO: the ".or. lbfield_constr" is not present in host. Is it only needed for externals?
    if(lbfield .or. lbfield_constr) then
      allocate (bfield%thetallmat((lmax+1)**2,(lmax+1)**2,irmdnew), stat=i_stat)
      call calc_thetallmat(bfield%thetallmat, lmax, ipan_intervall(npan_log+npan_eq) + 1, &
                           iend, irmdnew, thetasnew, ifunm, icleb, cleb)
    end if

  end subroutine init_bfield

  !TODO save_bconst()

  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise constraining field from bconstr.dat
  !> Author: MdSD, Nicolas Essing
  !>
  !> The file should contain three floating point values in each line and a line
  !> for each atom. The first line is treated as a header and skipped.
  !> Intepreted as initial constraining bfield in cartesian coordinates, in Ry.
  !-------------------------------------------------------------------------------
  subroutine read_bconstr(bfields, number_of_atoms)
    type(bfield_data), dimension(:), intent(inout) :: bfields
    integer                        , intent(in)  :: number_of_atoms
    
    integer :: iatom
    integer :: iostat
    logical :: file_exists

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
      end do
      close(57493215)
    else
      ! No 'bconstr_in.dat' given, use default 0
      do iatom = 1, number_of_atoms
        bfields(iatom)%bfield_constr(:) = 0.
      end do
    end if

    write(*,*) '  ############################################################'
    write(*,*) '  input constraining fields'
    write(*,*) '  ############################################################'
    write(*,*) '  iatom      Bx              By              Bz        (in Ry)'
    do iatom = 1, number_of_atoms
      write(*,'(2X,I4,3(E16.8))') iatom, bfields(iatom)%bfield_constr(:)
    end do
  end subroutine read_bconstr


  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise magnetic field from bfield.dat
  !> Author: Nicolas Essing
  !> 
  !> The file should contain three numbers per line and one line per atom.
  !> Read as 'theta  phi  bfield_strength' with the angles in degrees and the
  !> strength in Ry ( ! might get changed to Tesla ! )
  !> Lines can be commented out with a # as first character.
  !>------------------------------------------------------------------------------
  subroutine read_bfield(bfields,number_of_atoms)
    type(bfield_data), dimension(:), intent(inout) :: bfields
    integer                        , intent(in)   :: number_of_atoms
    
    integer        :: iatom, iostat
    character(256) :: linebuffer
   
    open(unit=57493215, file='bfield.dat', iostat=iostat)
    if (iostat /= 0) then
      write(*,*) "Error reading bfield.dat."
      write(*,*) "Provide correct file or turn off magnetic fields."
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
      read(linebuffer, *, iostat=iostat) bfields(iatom)%theta, bfields(iatom)%phi, bfields(iatom)%bfield_strength
      if (iostat /= 0) then
        write(*,*) "Error parsing a line in bfield.dat"
        stop
      end if
      iatom = iatom + 1
    end do
    close(57493215)

    write(*,*) '  ###############################################'
    write(*,*) '  external non-collinear magnetic fields'
    write(*,*) '  ###############################################'
    write(*,*) '  iatom      theta       phi         bfield (in Ry)'
    do iatom = 1, number_of_atoms
      bfields(iatom)%theta           = bfields(iatom)%theta / 360.0d0 * 8.d0 * datan(1.d0)
      bfields(iatom)%phi             = bfields(iatom)%phi   / 360.0d0 * 8.d0 * datan(1.d0)
      bfields(iatom)%bfield_strength = bfields(iatom)%bfield_strength ! / 235051.787 ! conversion from Tesla to Ry
      bfields(iatom)%bfield(1)       = bfields(iatom)%bfield_strength*cos(bfields(iatom)%phi)*sin(bfields(iatom)%theta)
      bfields(iatom)%bfield(2)       = bfields(iatom)%bfield_strength*sin(bfields(iatom)%phi)*sin(bfields(iatom)%theta)
      bfields(iatom)%bfield(3)       = bfields(iatom)%bfield_strength*cos(bfields(iatom)%theta)
      write(*,'(2X,I4,3(E16.8))') iatom, bfields(iatom)%theta, bfields(iatom)%phi, bfields(iatom)%bfield_strength
   end do
  end subroutine read_bfield

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
                        itscf0, itscf1, lbfield, lbfield_constr, lbfield_trans, &
                        lbfield_mt, transpose_bfield)
    type(bfield_data), intent(in) :: bfield
    double complex, dimension(:,:,:), intent(inout) :: vnspll ! The potential to add to
    double precision,  intent(in) :: theta, phi ! angles or the magnetic moments, not to be confused with theta and phi in bfield
    integer, intent(in) :: imt ! MT radius (index in cheb mesh)
    integer, intent(in) :: iteration_number !TODO this, or just a logical and do the check outside?
    integer, intent(in) :: itscf0, itscf1   !TODO ^
    logical, intent(in) :: lbfield, lbfield_constr !! Wheter to apply the fields
    logical, intent(in) :: lbfield_trans ! Apply only transveral
    logical, intent(in) :: lbfield_mt ! Apply only up do MT radius
    logical, intent(in) :: transpose_bfield ! Transpose the bfield (for left solutions)

    double complex, parameter :: cplx_i = (0.d0, 1.d0)
    integer :: lmmax, irmd, iend, ir ! loop boundaries and indices
    double complex, dimension(3) :: combined_bfields, dir ! vector of combined bfields and unit vector of magnetic moment direction
    double complex, dimension(2,2) :: bfield_mat ! bfield times pauli matrices
    double complex :: temp ! used to transpose the matrix

    lmmax = size(bfield%thetallmat, 1) ! size(vnspll, 1) is 2*lmmax
    irmd = size(vnspll, 3)

    combined_bfields(:) = bfield%bfield(:) ! start with external, is zero if unused
    if (lbfield_constr .and. itscf0 <= iteration_number .and. iteration_number <= itscf1) then
      ! Add constraint field
      combined_bfields(:) = combined_bfields(:) + bfield%bfield_constr(:)
    end if

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
  !> Category: KKRhost, geometry, new-mesh, shapefun
  !> Deprecated: False 
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

    double precision, parameter :: rfpi = 3.5449077018110318
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

    shapefun_mod(:,:)        = 0.d0
    shapefun_mod(1:imt1,1)   = rfpi ! is multipled by C_LL^0
    shapefun_mod(imt1+1:, 1) = thetasnew(imt1+1:irmdnew,1)

    ! convert from pointer to real l,m
    do lm1 = 2, lmmax2
      ifun = ifunm(lm1)
      if(.not. ifun == 0) then !shapefun%lmused(lm1)==1) then
        shapefun_mod(imt1+1:, lm1) = thetasnew(imt1+1:, ifun)
      end if
    end do
    
    thetansll(:,:,:) = 0.d0

    ! diagonal part (not contained in gaunt-coeff)
    do lm1 = 1, lmmax
      thetansll(lm1,lm1,:) = shapefun_mod(:,1) * c0ll
    end do

    ! This is effectivly a summation over one angular momentum index lm3 for each
    ! combination of lm1, lm2. Iterate instead over the flattened array cleb
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
