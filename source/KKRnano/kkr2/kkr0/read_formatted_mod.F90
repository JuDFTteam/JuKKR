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
    integer :: ITITLE(20)
    double precision :: RMT
    double precision :: ALAT
    double precision :: RMTNEW
    double precision :: Z_nuclear
    double precision :: RWS
    double precision :: EFERMI
    double precision :: VBC
    integer :: IRWS
    double precision :: A_log_mesh
    double precision :: B_log_mesh
  endtype

  type CoreStatesBlock
     integer :: NCORE
     integer :: INEW
     integer :: LCORE (20)
     double precision :: ECORE (20)
  endtype

  type SphericalBlock
     integer :: IRT1P
     integer :: IRNS
     integer :: LMPOT
     integer :: ISAVE
     double precision, allocatable :: VISP(:)
  endtype

  type NonSphericalBlocks
     double precision, allocatable :: VINS(:,:)
  endtype

  type PotentialEntry
    type (PotentialHeader) :: header
    type (CoreStatesBlock) :: csblock
    type (SphericalBlock)  :: sblock
    type (NonSphericalBlocks)  :: nsblocks
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
    type (PotentialHeader), intent(out) :: header
    integer, intent(in) :: unit

    READ (unit,FMT="(20a4)") header%ITITLE(1:20)

!---  >read muffin-tin radius , lattice constant and new muffin radius
!      (not used)
    READ (unit,FMT="(3f12.8)") header%RMT, header%ALAT, header%RMTNEW

!---> read nuclear charge
!     wigner seitz radius (not used), fermi energy and energy difference
!     between electrostatic zero and muffin tin zero (not used)

    READ (unit,FMT="(f10.5,/,f10.5,2f15.10)") header%Z_nuclear, header%RWS, header%EFERMI, header%VBC

!---> read : number of radial mesh points
!     (in case of ws input-potential: last mesh point corresponds
!     to ws-radius, in case of shape-corrected input-potential
!     last mesh point of the exponential mesh corresponds to
!     mt-radius/nevertheless this point is always in the array
!     irws(ih)),number of points for the radial non-muffin-tin
!     mesh  needed for shape functions, the constants a and b
!     for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
!     the no. of different core states and some other stuff

      READ (unit,FMT="(i3,/,2d15.8)") header%IRWS, header%A_log_mesh, header%B_log_mesh

  endsubroutine

  !----------------------------------------------------------------------------
  !> Read core state block.
  subroutine read_CoreStateBlock(block, unit)
    type (CoreStatesBlock), intent(out) :: block
    integer, intent(in) :: unit

    integer :: ICORE

    READ (unit,FMT="(2i2)") block%NCORE, block%inew

! read the different core states : l and energy

!          check: e.r.
    if (block%ncore .GT. 20) THEN
       write(*,*) "Error: More than 20 core states."
       STOP
    endif

    block%LCORE = -1
    block%ECORE = 9999.0d0

    IF (block%NCORE.GE.1) THEN
       DO ICORE=1,block%NCORE
          READ (unit,FMT="(i5,1p,d20.11)") block%LCORE(ICORE), block%ECORE(ICORE)
       ENDDO ! ICORE
    ENDIF

  endsubroutine

  !----------------------------------------------------------------------------
  !> Read SphericalBlock, do not forget to run destroy_SphericalBlock.
  subroutine create_read_SphericalBlock(block, unit)
    type (SphericalBlock), intent(out) :: block
    integer, intent(in) :: unit

    integer IR, NR

    READ (unit,FMT="(10i5)") block%IRT1P,block%IRNS,block%LMPOT,block%ISAVE

    NR = block%IRT1P
    allocate(block%VISP(NR))

    READ (unit,FMT="(1p,4d20.13)") (block%VISP(IR), IR=1,NR)

  endsubroutine ! create

    !----------------------------------------------------------------------------
  !> Deallocate array for VISP data
  subroutine destroy_SphericalBlock(block)
    type (SphericalBlock), intent(inout) :: block
    deallocate (block%VISP)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Read NonSphericalBlocks, do not forget to run destroy_NonSphericalBlocks.
  subroutine create_read_NonSphericalBlocks(blocks, spherical_block, unit)
    type (NonSphericalBlocks), intent(out) :: blocks
    type (SphericalBlock), intent(in) :: spherical_block
    integer, intent(in) :: unit

    integer LMPOT, ISAVE, IRMIN, IRMD
    integer LM, LM1, IR

    LMPOT = spherical_block%LMPOT
    ISAVE = spherical_block%ISAVE
    IRMD = spherical_block%IRT1P
    IRMIN = spherical_block%IRT1P - spherical_block%IRNS

    allocate(blocks%VINS(IRMIN:IRMD, LMPOT))
    blocks%VINS = 0.0d0

    IF (LMPOT.GT.1) THEN
      LM1 = 2
      DO LM = 2, LMPOT
        IF (LM1.NE.1) THEN

          IF (ISAVE.EQ.1) THEN
            READ (unit,FMT="(10i5)") LM1
          ELSE
            LM1 = LM
          ENDIF

          IF (LM1.GT.1) THEN

            if (LM1 < 1 .or. LM1 > LMPOT) then
              write(*,*) "ERROR: potential file is not correctly formatted."
              write(*,*) "Error when trying to read entry: ", LM1
              STOP
            endif

            READ (unit,FMT="(1p,4d20.13)") (blocks%VINS(IR, LM1),IR=IRMIN,IRMD)
          ENDIF
        ENDIF
      ENDDO ! LM
    ENDIF

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Deallocate array for VINS data.
  subroutine destroy_NonSphericalBlocks(blocks)
    type (NonSphericalBlocks), intent(inout) :: blocks
    deallocate(blocks%VINS)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Create a PotentialEntry by reading from file 'unit'.
  subroutine create_read_PotentialEntry(potential_entry, unit)
    type (PotentialEntry), intent(out) :: potential_entry
    integer, intent(in) :: unit

    call read_PotentialHeader(potential_entry%header, unit)
    call read_CoreStateBlock(potential_entry%csblock, unit)
    call create_read_SphericalBlock(potential_entry%sblock, unit)
    call create_read_NonSphericalBlocks(potential_entry%nsblocks, potential_entry%sblock, unit)
  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Destroy a PotentialEntry.
  subroutine destroy_PotentialEntry(potential_entry)
    type (PotentialEntry), intent(inout) :: potential_entry

    call destroy_SphericalBlock(potential_entry%sblock)
    call destroy_NonSphericalBlocks(potential_entry%nsblocks)
  endsubroutine ! destroy

endmodule read_formatted_mod

#ifdef TEST_READ_FORMATTED_MOD
program test_read_formatted
  use read_formatted_mod, only: PotentialEntry, create_read_PotentialEntry, destroy_PotentialEntry
  implicit none

  type(PotentialEntry) :: pe
  integer, parameter :: UNIT=42
  integer :: LM

  open(UNIT, form='formatted', file='potential')
  call create_read_PotentialEntry(pe, UNIT)

  write(*,fmt="(' <#',20a4)") pe%header%ITITLE
  write(*,*) pe%sblock%VISP
  write(*,*) "---------------------------------------------------------------"
  write(*,*) "Number of non-spherical components: ", pe%sblock%LMPOT

  do LM = 1, pe%sblock%LMPOT
    write(*,*) "---------------------------------------------------------------"
    write(*,*) "LM = ", LM
    write(*,*) "---------------------------------------------------------------"
    write(*,*) pe%nsblocks%VINS(:,LM)
  enddo ! LM

  call destroy_PotentialEntry(pe)
  close(UNIT)
endprogram
#endif

