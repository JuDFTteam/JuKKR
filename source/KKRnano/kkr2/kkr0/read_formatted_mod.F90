module read_formatted_mod
  implicit none

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
  end type

  type CoreStatesBlock
     integer :: NCORE
     integer :: INEW
     integer :: LCORE (20)
     double precision :: ECORE (20)
  end type

  type SphericalBlock
     integer :: IRT1P
     integer :: IRNS
     integer :: LMPOT
     integer :: ISAVE
     double precision, allocatable :: VISP(:)
  end type

  type NonSphericalBlocks
     double precision, allocatable :: VINS(:,:)
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  !> Read header of potential entry.
  subroutine read_PotentialHeader(header, unit)
    type (PotentialHeader), intent(out) :: header
    integer, intent(in) :: unit

    integer IA

    READ (unit,FMT=9020) (header%ITITLE(IA),IA=1,20)

!---  >read muffin-tin radius , lattice constant and new muffin radius
!      (not used)
    READ (unit,FMT=9030) header%RMT, header%ALAT, header%RMTNEW

!---> read nuclear charge
!     wigner seitz radius (not used), fermi energy and energy difference
!     between electrostatic zero and muffin tin zero (not used)

    READ (unit,FMT=9040) header%Z_nuclear, header%RWS, header%EFERMI, header%VBC

!---> read : number of radial mesh points
!     (in case of ws input-potential: last mesh point corresponds
!     to ws-radius, in case of shape-corrected input-potential
!     last mesh point of the exponential mesh corresponds to
!     mt-radius/nevertheless this point is always in the array
!     irws(ih)),number of points for the radial non-muffin-tin
!     mesh  needed for shape functions, the constants a and b
!     for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
!     the no. of different core states and some other stuff

      READ (unit,FMT=9050) header%IRWS, header%A_log_mesh, header%B_log_mesh

      9020 FORMAT (20a4)
      9030 FORMAT (3f12.8)
      9040 FORMAT (f10.5,/,f10.5,2f15.10)
      9050 FORMAT (i3,/,2d15.8)
  end subroutine

  !----------------------------------------------------------------------------
  !> Read core state block.
  subroutine read_CoreStateBlock(block, unit)
    type (CoreStatesBlock), intent(out) :: block
    integer, intent(in) :: unit

    integer :: ICORE

    READ (unit,FMT=9055) block%NCORE, block%inew

! read the different core states : l and energy

!          check: e.r.
    if (block%ncore .GT. 20) THEN
       write(*,*) "Error: More than 20 core states."
       STOP
    endif

    IF (block%NCORE.GE.1) THEN
       DO ICORE=1,block%NCORE
          READ (unit,FMT=9070) block%LCORE(ICORE), block%ECORE(ICORE)
       END DO
    END IF

    9055 FORMAT (2i2)
    9070 FORMAT (i5,1p,d20.11)

  end subroutine

  !----------------------------------------------------------------------------
  !> Read SphericalBlock, do not forget to run destroy_SphericalBlock.
  subroutine create_read_SphericalBlock(block, unit)
    type (SphericalBlock), intent(out) :: block
    integer, intent(in) :: unit

    integer IR, NR

    READ (unit,FMT=9090) block%IRT1P,block%IRNS,block%LMPOT,block%ISAVE

    NR = block%IRT1P
    allocate(block%VISP(NR))

    READ (unit,FMT=9100) (block%VISP(IR), IR=1,NR)

9090 FORMAT (10i5)
9100 FORMAT (1p,4d20.13)

  end subroutine

    !----------------------------------------------------------------------------
  !> Deallocate array for VISP data
  subroutine destroy_SphericalBlock(block)
    type (SphericalBlock), intent(inout) :: block
    deallocate (block%VISP)
  end subroutine

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
      DO 50 LM = 2, LMPOT
        IF (LM1.NE.1) THEN

          IF (ISAVE.EQ.1) THEN
            READ (unit,FMT=9090) LM1
          ELSE
            LM1 = LM
          END IF

          IF (LM1.GT.1) THEN

            if (LM1 < 1 .or. LM1 > LMPOT) then
              write(*,*) "ERROR: potential file is not correctly formatted."
              write(*,*) "Error when trying to read entry: ", LM1
              STOP
            endif

            READ (unit,FMT=9100) (blocks%VINS(IR, LM1),IR=IRMIN,IRMD)
          END IF

        END IF

50    CONTINUE

    END IF

9090 FORMAT (10i5)
9100 FORMAT (1p,4d20.13)

  end subroutine

end module read_formatted_mod
