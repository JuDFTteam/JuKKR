module InterpolateBasisAtom_mod
  use BasisAtom_mod
  use RadialMeshData_mod
  use PotentialData_mod
  use AtomicCoreData_mod

  implicit none

  CONTAINS

  ! new mesh must be persistent
  ! new_atom gets allocated
  ! old_atom unmodified
  subroutine interpolateBasisAtom(new_atom, old_atom, new_mesh)
    implicit none
    type(BasisAtom), intent(inout) :: new_atom
    type(BasisAtom), intent(in) :: old_atom
    type(RadialMeshData), intent(in) :: new_mesh

    integer :: nspin
    integer :: lpot
    integer :: irmin_old, irmin_new
    integer :: irws_old, irws_new
    integer :: ii, lm
    type (RadialMeshData), pointer :: old_mesh

    nspin = old_atom%nspin
    lpot = old_atom%potential%lpot
    old_mesh => old_atom%mesh_ptr

    call createBasisAtom(new_atom, old_atom%atom_index, lpot, nspin, &
                         new_mesh%irmin, new_mesh%irmd)


    new_atom%atom_index = old_atom%atom_index
    new_atom%cell_index = old_atom%cell_index
    new_atom%cluster_index = old_atom%cluster_index
    new_atom%Z_nuclear = old_atom%Z_nuclear
    new_atom%RMTref = old_atom%RMTref
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

    ! TODO: scale?
    ! interpolate non-spherical potential
    do lm = 1, old_atom%potential%lmpot
      do ii = 1, nspin
        call interpolate(old_mesh%r(irmin_old:irws_old), &
                         old_atom%potential%VINS(irmin_old:irws_old,lm,ii), &
                         new_mesh%r(irmin_new:irws_new), &
                         new_atom%potential%VINS(irmin_new:irws_new,lm,ii))
      end do
    end do

    ! interpolate spherical potential
    do ii = 1, nspin
      call interpolate(old_mesh%r, old_atom%potential%VISP(:,ii), &
                       new_mesh%r, new_atom%potential%VISP(:,ii))
    end do
  end subroutine


  !----------------------------------------------------------------------------
  !> Interpolate function values with cubic splines
  subroutine interpolate(xval, yval, xnew, ynew)
    implicit none
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
    double precision, parameter :: TOL = 1.d-14

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
      end if
    end do

    call spline(num, xarray, yarray, counter, 1.d35, 1.d35, y2ndder)

    do ii = 1, size(ynew)
      call splint(xarray, yarray, y2ndder, counter, xnew(ii), ynew(ii), yderiv)
    end do

  end subroutine
end module InterpolateBasisAtom_mod
