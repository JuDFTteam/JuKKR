!> Module to read formatted potential.
!>
!> @author Elias Rabel, Marcel Bornemann
!> 2015

module read_formatted_mod
  implicit none
  private

  ! use the following 2 routines to read one potential entry from a file.
  public :: PotentialEntry, create, destroy
  public :: create_read_PotentialEntry, destroy_PotentialEntry

  type PotentialHeader
    integer :: ititle(20)
    double precision :: rmt
    double precision :: alat
    double precision :: rmtnew
    double precision :: z_nuclear
    double precision :: rws
    double precision :: efermi
    double precision :: vbc
    integer :: irws
    double precision :: a_log_mesh
    double precision :: b_log_mesh
  endtype

  type CoreStatesBlock
     integer :: ncore
     integer :: inew
     integer :: lcore(20)
     double precision :: ecore(20)
  endtype

  type SphericalBlock
     integer :: irt1p
     integer :: irns
     integer :: lmpot
     integer :: isave
     double precision, allocatable :: visp(:)
  endtype

  type NonSphericalBlocks
     double precision, allocatable :: vins(:,:)
  endtype

  type PotentialEntry
    type(PotentialHeader) :: header
    type(CoreStatesBlock) :: csblock
    type(SphericalBlock) :: sblock
    type(NonSphericalBlocks) :: nsblocks
  endtype
  
  interface create
    module procedure create_read_PotentialEntry
  endinterface
  
  interface destroy
    module procedure destroy_PotentialEntry
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  !> Read header of potential entry.
  subroutine read_PotentialHeader(header, unit)
    type(PotentialHeader), intent(out) :: header
    integer, intent(in) :: unit

    read (unit,fmt="(20a4)") header%ititle(1:20)

!---  >read muffin-tin radius , lattice constant and new muffin radius
!      (not used)
    read (unit,fmt="(3f12.8)") header%rmt, header%alat, header%rmtnew

!---> read nuclear charge
!     wigner seitz radius (not used), fermi energy and energy difference
!     between electrostatic zero and muffin tin zero (not used)

    read (unit,fmt="(f10.5,/,f10.5,2f15.10)") header%z_nuclear, header%rws, header%efermi, header%vbc

!---> read : number of radial mesh points
!     (in case of ws input-potential: last mesh point corresponds
!     to ws-radius, in case of shape-corrected input-potential
!     last mesh point of the exponential mesh corresponds to
!     mt-radius/nevertheless this point is always in the array
!     irws(ih)),number of points for the radial non-muffin-tin
!     mesh  needed for shape nctions, the constants a and b
!     for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
!     the no. of different core states and some other stuff

      read (unit,fmt="(i3,/,2d15.8)") header%irws, header%a_log_mesh, header%b_log_mesh

  endsubroutine

  !----------------------------------------------------------------------------
  !> Read core state block.
  subroutine read_CoreStateBlock(block, unit)
    type(CoreStatesBlock), intent(out) :: block
    integer, intent(in) :: unit

    integer :: icore

    read (unit,fmt="(2i2)") block%ncore, block%inew

! read the different core states : l and energy
    if (block%ncore  >  20) stop "Error: More than 20 core states."

    block%lcore(:) = -1
    block%ecore(:) = 9999.0d0

    do icore = 1, block%ncore
      read(unit,fmt="(i5,1p,d20.11)") block%lcore(icore), block%ecore(icore)
    enddo ! icore

  endsubroutine

  !----------------------------------------------------------------------------
  !> Read SphericalBlock, do not forget to run destroy_SphericalBlock.
  subroutine create_read_SphericalBlock(block, unit)
    type(SphericalBlock), intent(out) :: block
    integer, intent(in) :: unit

    read (unit,fmt="(10i5)") block%irt1p, block%irns, block%lmpot, block%isave
    allocate(block%visp(block%irt1p))
    read (unit,fmt="(1p,4d20.13)") block%visp(1:block%irt1p)

  endsubroutine ! create

    !----------------------------------------------------------------------------
  !> Deallocate array for VISP data
  subroutine destroy_SphericalBlock(block)
    type(SphericalBlock), intent(inout) :: block
    deallocate(block%VISP)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Read NonSphericalBlocks, do not forget to run destroy_NonSphericalBlocks.
  subroutine create_read_NonSphericalBlocks(blocks, sb, unit)
    type(NonSphericalBlocks), intent(out) :: blocks
    type(SphericalBlock), intent(in) :: sb ! only integers members from the spherical block are read-accessed
    integer, intent(in) :: unit

    integer :: irmin, lm, lm1

    irmin = sb%irt1p - sb%irns

    allocate(blocks%vins(irmin:sb%irt1p,1:sb%lmpot))
    blocks%vins = 0.d0

    lm1 = 2
    do lm = 2, sb%lmpot
      if (lm1 /= 1) then

        if (sb%isave == 1) then
          read (unit,fmt="(10i5)") lm1
        else
          lm1 = lm
        endif

        if (lm1 > 1) then

          if (lm1 < 1 .or. lm1 > sb%lmpot) then
            write(*,*) "ERROR: potential file is not correctly formatted."
            write(*,*) "Error when trying to read entry: ", LM1
            stop
          endif

          read (unit,fmt="(1p,4d20.13)") blocks%vins(irmin:sb%irt1p,lm1)
        endif
      endif
    enddo ! lm

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Deallocate array for VINS data.
  subroutine destroy_NonSphericalBlocks(blocks)
    type(NonSphericalBlocks), intent(inout) :: blocks
    deallocate(blocks%VINS)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Create a PotentialEntry by reading from file 'unit'.
  subroutine create_read_PotentialEntry(pe, unit)
    type(PotentialEntry), intent(out) :: pe
    integer, intent(in) :: unit

    call read_PotentialHeader(pe%header, unit)
    call read_CoreStateBlock(pe%csblock, unit)
    call create_read_SphericalBlock(pe%sblock, unit)
    call create_read_NonSphericalBlocks(pe%nsblocks, pe%sblock, unit)
  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Destroy a PotentialEntry.
  subroutine destroy_PotentialEntry(pe)
    type(PotentialEntry), intent(inout) :: pe

    call destroy_SphericalBlock(pe%sblock)
    call destroy_NonSphericalBlocks(pe%nsblocks)
  endsubroutine ! destroy

endmodule read_formatted_mod

#ifdef TEST_READ_FORMATTED_MOD
program test_read_formatted
  use read_formatted_mod, only: PotentialEntry, create_read_PotentialEntry, destroy_PotentialEntry
  implicit none

  type(PotentialEntry) :: pe
  integer, parameter :: fu=42
  integer :: lm

  open(fu, form='formatted', file='potential')
  call create_read_PotentialEntry(pe, fu)

  write(*,fmt="(' <#',20a4)") pe%header%ITITLE
  write(*,*) pe%sblock%VISP
  write(*,*) "---------------------------------------------------------------"
  write(*,*) "Number of non-spherical components: ", pe%sblock%LMPOT

  do lm = 1, pe%sblock%LMPOT
    write(*,*) "---------------------------------------------------------------"
    write(*,*) "LM = ", lm
    write(*,*) "---------------------------------------------------------------"
    write(*,*) pe%nsblocks%VINS(:,lm)
  enddo ! lm

  call destroy_PotentialEntry(pe)
  close(fu)
endprogram
#endif

