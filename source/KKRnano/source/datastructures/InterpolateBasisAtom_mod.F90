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
  ! Interpolate 'olda' to 'new_mesh' and return result in 'newa' using
  ! lpot = lpot_new (can be different from that in olda).
  ! new mesh must be persistent
  ! newa gets allocated
  ! olda unmodified
  subroutine interpolateBasisAtom(newa, olda, new_mesh, lpot_new)
    use BasisAtom_mod, only: BasisAtom, create, associateBasisAtomMesh
    use RadialMeshData_mod, only: RadialMeshData
    type(BasisAtom), intent(inout) :: newa
    type(BasisAtom), intent(in) :: olda
    type(RadialMeshData), intent(in) :: new_mesh
    integer, intent(in) :: lpot_new

    type(RadialMeshData), pointer :: old_mesh
    double precision, parameter :: TOL = 1.d-8
    integer :: nspin
    integer :: lmpot_min
    integer :: irmin_old, irmin_new
    integer :: irws_old, irws_new
    integer :: ii, lm
    logical :: do_interpolation

    do_interpolation = .true.

    nspin = olda%nspin
    old_mesh => olda%mesh_ptr

    call create(newa, olda%atom_index, lpot_new, nspin, new_mesh%irmin, new_mesh%irmd) ! createBasisAtom

    lmpot_min = min(olda%potential%lmpot, newa%potential%lmpot)

    newa%atom_index  = olda%atom_index
    newa%cell_index  = olda%cell_index
    newa%Z_nuclear   = olda%Z_nuclear
    newa%radius_muffin_tin = olda%radius_muffin_tin
    newa%core%LCORE  = olda%core%LCORE
    newa%core%NCORE  = olda%core%NCORE
    newa%core%ECORE  = olda%core%ECORE
    newa%core%ITITLE = olda%core%ITITLE

    newa%cell_ptr => olda%cell_ptr

    call associateBasisAtomMesh(newa, new_mesh)

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
        olda%potential%lmpot == newa%potential%lmpot) then
      if (sum(abs(old_mesh%r - new_mesh%r)) < TOL .and. &
          irmin_old == irmin_new .and. irws_old == irws_new) then
        ! mesh has not changed - do not interpolate
        ! this avoids numerical inaccuracies - arising from
        ! spline interpolation routine
        do_interpolation = .false.
      endif
    endif

    if (do_interpolation) then

      newa%potential%VINS = 0.d0
      newa%potential%VISP = 0.d0

      ! TODO: scale?
      ! interpolate non-spherical potential - use only non-spherical part of mesh
      do ii = 1, nspin
        do lm = 1, lmpot_min
          call interpolate(old_mesh%r(irmin_old:irws_old), &
                           olda%potential%VINS(irmin_old:irws_old,lm,ii), &
                           new_mesh%r(irmin_new:irws_new), &
                           newa%potential%VINS(irmin_new:irws_new,lm,ii))
        enddo ! lm
      enddo ! ii

      ! be careful when mesh partition has changed
      if (irmin_new < irmin_old) then
        newa%potential%VINS(irmin_new:min(irmin_old, irws_new),:,:) = 0.d0
      endif

      ! interpolate spherical potential
      do ii = 1, nspin
        call interpolate(old_mesh%r, olda%potential%VISP(:,ii), &
                         new_mesh%r, newa%potential%VISP(:,ii))
      enddo ! ii

    else
      ! no interpolation - just copy
      newa%potential = olda%potential

    endif
  endsubroutine ! interpolate


  !----------------------------------------------------------------------------
  !> Interpolate function values with cubic splines
  !> @param[out] ynew  interpolated function values
  subroutine interpolate(xval, yval, xnew, ynew)
    ! use Splines_mod, only: spline, splint
    double precision, intent(in) :: xval(:)
    double precision, intent(in) :: yval(:)
    double precision, intent(in) :: xnew(:)
    double precision, intent(out) :: ynew(:)

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

    ! note: 1st derivative at upper boundary forced to 0.d0                                                  
    call spline(num, xarray, yarray, counter, 1.d35, 0.d0, y2ndder)

    do ii = 1, size(ynew)

      if (xnew(ii) <= xarray(counter)) then
        call splint(xarray, yarray, y2ndder, counter, xnew(ii), ynew(ii), yderiv)
      else
        ! use constant value for x > x_max_old
        ynew(ii) = yarray(counter)
      endif

    enddo ! ii

  endsubroutine ! interpolate
  
endmodule ! InterpolateBasisAtom_mod
