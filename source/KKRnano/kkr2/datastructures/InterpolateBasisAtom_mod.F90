!> Module to interpolate potential of a basis atom to different mesh
!> using cubic spline interpolation.
!> @author Elias Rabel
!
module InterpolateBasisAtom_mod
  implicit none
  private
  public :: interpolateBasisAtom

  contains

  !----------------------------------------------------------------------------
  ! Interpolate 'old_atom' to 'new_mesh' and return result in 'new_atom' using
  ! lpot = lpot_new (can be different from that in old_atom).
  ! new mesh must be persistent
  ! new_atom gets allocated
  ! old_atom unmodified
  subroutine interpolateBasisAtom(new_atom, old_atom, new_mesh, lpot_new)
    use BasisAtom_mod, only: BasisAtom, associateBasisAtomMesh
    use BasisAtom_mod, only: createBasisAtom ! deprecated
    use RadialMeshData_mod, only: RadialMeshData
    type(BasisAtom), intent(inout) :: new_atom
    type(BasisAtom), intent(in) :: old_atom
    type(RadialMeshData), intent(in) :: new_mesh
    integer, intent(in) :: lpot_new

    integer :: nspin
    integer :: lmpot_min
    integer :: irmin_old, irmin_new
    integer :: irws_old, irws_new
    integer :: ii, lm
    type (RadialMeshData), pointer :: old_mesh
    logical :: do_interpolation
    double precision, parameter :: TOL = 1.d-8

    do_interpolation = .true.

    nspin = old_atom%nspin
    old_mesh => old_atom%mesh_ptr

    call createBasisAtom(new_atom, old_atom%atom_index, lpot_new, nspin, &
                         new_mesh%irmin, new_mesh%irmd)

    lmpot_min = min(old_atom%potential%lmpot, new_atom%potential%lmpot)

    new_atom%atom_index = old_atom%atom_index
    new_atom%cell_index = old_atom%cell_index
    new_atom%Z_nuclear = old_atom%Z_nuclear
    new_atom%radius_muffin_tin = old_atom%radius_muffin_tin
    new_atom%core%LCORE = old_atom%core%LCORE
    new_atom%core%NCORE = old_atom%core%NCORE
    new_atom%core%ECORE = old_atom%core%ECORE
    new_atom%core%ITITLE = old_atom%core%ITITLE

    new_atom%cell_ptr => old_atom%cell_ptr

    call associateBasisAtomMesh(new_atom, new_mesh)

    ! get old and new non-spherical mesh bounds
    irmin_old = old_mesh%irmin
    irmin_new = new_mesh%irmin
    irws_old = old_mesh%irws
    irws_new = new_mesh%irws

    ! check if interpolation is really necessary
    ! Check if mesh has changed!
    ! criterion: change in number of points OR
    !            sum(abs(mesh%r - old_mesh%r)) > 1.d-8

    do_interpolation = .true.
    if (size(old_mesh%r) == size(new_mesh%r) .and. &
        old_atom%potential%lmpot == new_atom%potential%lmpot) then
      if (sum(abs(old_mesh%r - new_mesh%r)) < TOL .and. &
          irmin_old == irmin_new .and. irws_old == irws_new) then
        ! mesh has not changed - don't interpolate
        ! this avoids numerical inaccuracies - arising from
        ! spline interpolation routine
        do_interpolation = .false.
      endif
    endif

    if (do_interpolation .eqv. .true.) then

      new_atom%potential%VINS = 0.0d0
      new_atom%potential%VISP = 0.0d0

      ! TODO: scale?
      ! interpolate non-spherical potential - use only non-spherical part of mesh
      do ii = 1, nspin
        do lm = 1, lmpot_min
          call interpolate(old_mesh%r(irmin_old:irws_old), &
                           old_atom%potential%VINS(irmin_old:irws_old,lm,ii), &
                           new_mesh%r(irmin_new:irws_new), &
                           new_atom%potential%VINS(irmin_new:irws_new,lm,ii))
        enddo ! lm
      enddo ! ii

      ! be careful when mesh partition has changed
      if (irmin_new < irmin_old) then
        new_atom%potential%VINS(irmin_new:min(irmin_old, irws_new),:,:) = 0.0d0
      endif

      ! interpolate spherical potential
      do ii = 1, nspin
        call interpolate(old_mesh%r, old_atom%potential%VISP(:,ii), &
                         new_mesh%r, new_atom%potential%VISP(:,ii))
      enddo ! ii

    else
      ! no interpolation - just copy
      new_atom%potential = old_atom%potential

    endif
  endsubroutine interpolateBasisAtom


  !----------------------------------------------------------------------------
  !> Interpolate function values with cubic splines
  !> @param[out] ynew  interpolated function values
  subroutine interpolate(xval, yval, xnew, ynew)
    ! use Splines_mod, only: spline
    double precision, intent(in) :: xval(:)
    double precision, intent(in) :: yval(:)
    double precision, intent(in) :: xnew(:)
    double precision, intent(out) :: ynew(:)

    !-----------------
    double precision, allocatable :: xarray(:)
    double precision, allocatable :: yarray(:)
    double precision, allocatable :: y2ndder(:)
    double precision :: yderiv
    integer :: num
    integer :: counter, ii
    double precision, parameter :: TOL = 1.d-10

    num = size(xval)
    allocate(xarray(num))
    allocate(yarray(num))
    allocate(y2ndder(num))

    counter = 1
    xarray(1) = xval(1)
    yarray(1) = yval(1)
    ! remove duplicate x - values (panels!)
    ! assumption: no discontinuities (jumps)!
    do ii = 2, num
      if (abs(xval(ii) - xval(ii-1)) > TOL) then
        counter = counter + 1
        xarray(counter) = xval(ii)
        yarray(counter) = yval(ii)
      endif
    enddo ! ii

    ! note: 1st derivative at upper boundary forced to 0.0d0                                                  
    call spline(num, xarray, yarray, counter, 1.d35, 0.0d0, y2ndder)

    do ii = 1, size(ynew)

      if (xnew(ii) <= xarray(counter)) then
        call splint(xarray, yarray, y2ndder, counter, xnew(ii), ynew(ii), yderiv)
      else
        ! use constant value for x > x_max_old
        ynew(ii) = yarray(counter)
      endif

    enddo ! ii

  endsubroutine interpolate
  
endmodule InterpolateBasisAtom_mod
