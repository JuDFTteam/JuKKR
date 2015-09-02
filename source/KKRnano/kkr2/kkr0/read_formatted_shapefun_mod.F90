#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!> Module to read shape function file.
!>
!> @author Elias Rabel, Marcel Bornemann
!> 2015
module read_formatted_shapefun_mod
  implicit none
  private
  
  ! use the following 2 routines to read a shape function file.
  public :: ShapefunFile, create, destroy
  public :: create_read_ShapefunFile, destroy_ShapefunFile ! deprecated

  type Intermesh
    integer :: NPAN
    integer :: MESHN
    integer, allocatable :: NM(:)
    double precision, allocatable :: XRN(:)
    double precision, allocatable :: DRN(:)
  endtype

  type Shapefunction
    integer :: NFU
    integer, allocatable :: LLMSP(:)
    double precision, allocatable :: THETAS(:,:)
  endtype

  type ShapefunFile
    integer NCELL
    type(Intermesh), allocatable :: mesh(:)
    type(Shapefunction), allocatable :: shapes(:)
  endtype

  interface create
    module procedure create_read_ShapefunFile
  endinterface
  
  interface destroy
    module procedure destroy_ShapefunFile
  endinterface
  
  contains

  !--------------------------------------------------------------------------
  subroutine create_read_intermesh (inter, unit)
    type(Intermesh), intent(out) :: inter
    integer, intent(in) :: unit

    integer :: IPAN, IR

    READ (unit,FMT="(16i5)") inter%NPAN,inter%MESHN

    CHECKASSERT(inter%NPAN>=0)
    CHECKASSERT(inter%MESHN>=0)

    allocate (inter%NM(inter%NPAN))
    allocate (inter%XRN(inter%MESHN))
    allocate (inter%DRN(inter%MESHN))

    READ (unit,FMT="(16i5)") (inter%NM(IPAN), IPAN=1,inter%NPAN)
    READ (unit,FMT="(4d20.12)") (inter%XRN(IR), inter%DRN(IR), IR=1,inter%MESHN)

  endsubroutine ! create

  !--------------------------------------------------------------------------
  subroutine destroy_intermesh (inter)
    type(Intermesh), intent(inout) :: inter

    deallocate (inter%NM)
    deallocate (inter%XRN)
    deallocate (inter%DRN)

  endsubroutine ! destroy

  !--------------------------------------------------------------------------
  subroutine create_read_shapefunction (shapef, inter, unit)
    type(Shapefunction), intent(inout) :: shapef
    type(Intermesh), intent(in) :: inter
    integer, intent(in) :: unit

    integer :: IFUN, N, LM

    READ (unit,FMT="(16i5)") shapef%NFU

    allocate (shapef%LLMSP(shapef%NFU))
    allocate (shapef%THETAS(inter%MESHN,shapef%NFU))

    DO IFUN = 1, shapef%NFU
          READ (unit,FMT="(16i5)") LM
          CHECKASSERT(LM > 0)
          shapef%LLMSP(IFUN) = LM
          READ (unit,FMT="(4d20.12)") (shapef%THETAS(N,IFUN),N=1,inter%MESHN)
    ENDDO ! IFUN
    
  endsubroutine ! create

  !--------------------------------------------------------------------------
  subroutine destroy_shapefunction (shapef)
    type(Shapefunction), intent(inout) :: shapef

    deallocate (shapef%LLMSP)
    deallocate (shapef%THETAS)
  endsubroutine ! destroy

  subroutine create_read_shapefunFile (sfile, unit)
    type(ShapefunFile), intent (out) :: sfile
    integer, intent(in) :: unit

    integer :: ICELL
    double precision :: dummy

    READ (unit,FMT="(16i5)") sfile%NCELL

    READ (unit, FMT="(4d20.12)") (dummy, ICELL=1,sfile%NCELL)

    allocate (sfile%mesh(sfile%NCELL))
    allocate (sfile%shapes(sfile%NCELL))

    DO ICELL = 1, sfile%NCELL
      CALL create_read_intermesh(sfile%mesh(ICELL), unit)
      CALL create_read_shapefunction(sfile%shapes(ICELL), sfile%mesh(ICELL), unit)
    ENDDO ! ICELL

  endsubroutine ! create

  subroutine destroy_shapefunFile (sfile)
    type(ShapefunFile), intent (inout) :: sfile

    integer :: ICELL

    DO ICELL = 1, sfile%NCELL
      CALL destroy_intermesh(sfile%mesh(ICELL))
      CALL destroy_shapefunction(sfile%shapes(ICELL))
    ENDDO ! ICELL

    deallocate (sfile%mesh)
    deallocate (sfile%shapes)

  endsubroutine ! destroy

endmodule read_formatted_shapefun_mod

#ifdef TEST_READ_FORMATTED_SHAPEFUN_MOD
program test_read_formatted
  use read_formatted_shapefun_mod, only: create_read_ShapefunFile, destroy_ShapefunFile
  implicit none

  type(ShapefunFile) :: sfile
  integer, parameter :: UNIT = 42
  integer :: IFUN, ICELL

  open(UNIT, form='formatted', file='shapefun')
  call create_read_ShapefunFile(sfile, UNIT)

  do ICELL = 1, sfile%NCELL

    write(*,*) "NM :", sfile%mesh(ICELL)%NM

    do IFUN = 1, sfile%shapes(ICELL)%NFU
      write(*,*) "LM = ", sfile%shapes(ICELL)%LLMSP(IFUN)
      write(*,*) "-----------------------------------------------"
      write(*,*) sfile%shapes(ICELL)%THETAS(:, IFUN)
    enddo

  enddo

  call destroy_ShapefunFile(sfile)
  close(UNIT)

endprogram
#endif
