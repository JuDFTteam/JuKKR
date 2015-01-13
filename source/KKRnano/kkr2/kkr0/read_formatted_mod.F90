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

  CONTAINS

  !----------------------------------------------------------------------------
  !> Read header of potential entry.
  subroutine read_potential_header(header, unit)
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

  !> Read header of potential entry.
  subroutine read_CoreStateBlock(block, unit)
    type (CoreStatesBlock), intent(out) :: block
    integer, intent(in) :: unit

    integer :: ICORE

    READ (IFILE,FMT=9055) block%NCORE, block%inew

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

end module read_formatted_mod
