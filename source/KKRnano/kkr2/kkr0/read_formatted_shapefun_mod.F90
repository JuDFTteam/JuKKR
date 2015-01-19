#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!> Module to read shape function file.
!>
!> @author Elias Rabel, Marcel Bornemann
!> 2015
module read_formatted_shapefun_mod
    implicit none

    ! use the following 2 routines to read a shape function file.
    public :: create_read_ShapefunFile
    public :: destroy_ShapefunFile

    type Intermesh
      integer :: NPAN
      integer :: MESHN
      integer, allocatable :: NM (:)
      double precision, allocatable :: XRN (:)
      double precision, allocatable :: DRN (:)
    end type

    type Shapefunction
      integer :: NFU
      integer, allocatable :: LLMSP (:)
      double precision, allocatable :: THETAS (:,:)
    end type

    type ShapefunFile
      integer NCELL
      type (Intermesh), allocatable :: mesh (:)
      type (Shapefunction), allocatable :: shapes (:)
    end type

    contains

    !--------------------------------------------------------------------------
    subroutine create_read_intermesh (inter, unit)

      type (Intermesh), intent(out) :: inter
      integer, intent(in) :: unit

      integer :: IPAN, IR

      READ (unit,FMT=9000) inter%NPAN,inter%MESHN

      CHECKASSERT(inter%NPAN>=0)
      CHECKASSERT(inter%MESHN>=0)

      allocate (inter%NM(inter%NPAN))
      allocate (inter%XRN(inter%MESHN))
      allocate (inter%DRN(inter%MESHN))

      READ (unit,FMT=9000) (inter%NM(IPAN), IPAN=1,inter%NPAN)
      READ (unit,FMT=9010) (inter%XRN(IR), inter%DRN(IR), IR=1,inter%MESHN)

      9000 FORMAT (16i5)
      9010 FORMAT (4d20.12)
    end subroutine

    !--------------------------------------------------------------------------
    subroutine destroy_intermesh (inter)
      type (Intermesh), intent(inout) :: inter

      deallocate (inter%NM)
      deallocate (inter%XRN)
      deallocate (inter%DRN)

    end subroutine

    !--------------------------------------------------------------------------
    subroutine create_read_shapefunction (shapef, inter, unit)
      type (Shapefunction), intent(inout) :: shapef
      type (Intermesh), intent(in) :: inter
      integer, intent(in) :: unit

      integer :: IFUN, N, LM

      READ (unit,FMT=9000) shapef%NFU

      allocate (shapef%LLMSP(shapef%NFU))
      allocate (shapef%THETAS(inter%MESHN,shapef%NFU))

      DO IFUN = 1, shapef%NFU
            READ (unit,FMT=9000) LM
            CHECKASSERT(LM > 0)
            shapef%LLMSP(IFUN) = LM
            READ (unit,FMT=9010) (shapef%THETAS(N,IFUN),N=1,inter%MESHN)
      ENDDO

      9000 FORMAT (16i5)
      9010 FORMAT (4d20.12)
    end subroutine

    !--------------------------------------------------------------------------
    subroutine destroy_shapefunction (shapef)
      type (Shapefunction), intent(inout) :: shapef

      deallocate (shapef%LLMSP)
      deallocate (shapef%THETAS)
    end subroutine

    subroutine create_read_shapefunfile (sfile, unit)
      type(ShapefunFile), intent (out) :: sfile
      integer, intent(in) :: unit

      integer :: ICELL
      double precision :: dummy

      READ (unit,FMT=9000) sfile%NCELL

      READ (unit, FMT=9010) (dummy, ICELL=1,sfile%NCELL)

      allocate (sfile%mesh(sfile%NCELL))
      allocate (sfile%shapes(sfile%NCELL))

      DO ICELL = 1, sfile%NCELL
        CALL create_read_intermesh(sfile%mesh(ICELL), unit)
        CALL create_read_shapefunction(sfile%shapes(ICELL), sfile%mesh(ICELL), unit)
      ENDDO

      9000 FORMAT (16i5)
      9010 FORMAT (4d20.12)
    end subroutine

    subroutine destroy_shapefunfile (sfile)
      type(ShapefunFile), intent (inout) :: sfile

      integer :: ICELL

      DO ICELL = 1, sfile%NCELL
        CALL destroy_intermesh(sfile%mesh(ICELL))
        CALL destroy_shapefunction(sfile%shapes(ICELL))
      ENDDO

      deallocate (sfile%mesh)
      deallocate (sfile%shapes)

    end subroutine

end module read_formatted_shapefun_mod

#ifdef TEST_READ_FORMATTED_SHAPEFUN_MOD
program test_read_formatted
  use read_formatted_shapefun_mod
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

    call destroy_Shapefunfile(sfile)
  close(UNIT)

end program
#endif
