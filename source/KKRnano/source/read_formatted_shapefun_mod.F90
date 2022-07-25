!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Module to read shape function file.
!>
!> @author Elias Rabel, Marcel Bornemann
!> 2015
module read_formatted_shapefun_mod
#include "macros.h"
  implicit none
  private
  
  ! use the following 2 routines to read a shape function file.
  public :: ShapefunFile, create, destroy

  type Intermesh
    integer :: npan
    integer :: meshn
    integer, allocatable :: nm(:)
    double precision, allocatable :: xrn(:)
    double precision, allocatable :: drn(:)
  endtype

  type Shapefunction
    integer :: nfu
    integer, allocatable :: llmsp(:)
    double precision, allocatable :: thetas(:,:)
  endtype

  type ShapefunFile
    integer :: ncell
    type(Intermesh), allocatable :: mesh(:)
    type(Shapefunction), allocatable :: shapes(:)
  endtype

  interface create
    module procedure create_read_ShapefunFile
  endinterface
  
  interface destroy
    module procedure destroy_ShapefunFile, destroy_intermesh, destroy_shapefunction
  endinterface
  
  contains

  subroutine create_read_ShapefunFile(sfile, unit)
    type(ShapefunFile), intent(out) :: sfile
    integer, intent(in) :: unit

    integer :: icell
    double precision :: dummy

    read(unit, fmt="(16i5)") sfile%ncell

    read(unit, fmt="(4d20.12)") (dummy, icell=1,sfile%ncell)

    allocate(sfile%mesh(sfile%ncell), sfile%shapes(sfile%ncell))

    do icell = 1, sfile%ncell
      call create_read_intermesh(sfile%mesh(icell), unit)
      call create_read_shapefunction(sfile%shapes(icell), sfile%mesh(icell), unit)
    enddo ! icell

  endsubroutine ! create
  
  subroutine create_read_intermesh(inter, unit)
    type(Intermesh), intent(out) :: inter
    integer, intent(in) :: unit

    integer :: ir

    read(unit, fmt="(16i5)") inter%npan, inter%meshn

    CHECKASSERT(inter%npan >= 0)
    CHECKASSERT(inter%meshn >= 0)

    allocate (inter%nm(inter%npan))
    allocate (inter%xrn(inter%meshn))
    allocate (inter%drn(inter%meshn))

    read(unit, fmt="(16i5)") inter%nm(1:inter%npan)
    read(unit, fmt="(4d20.12)") (inter%xrn(ir), inter%drn(ir), ir=1,inter%meshn)

  endsubroutine ! create

  subroutine create_read_shapefunction(shapef, inter, unit)
    type(Shapefunction), intent(inout) :: shapef
    type(Intermesh), intent(in) :: inter
    integer, intent(in) :: unit

    integer :: ifun, lm

    read(unit, fmt="(16i5)") shapef%nfu

    allocate(shapef%llmsp(shapef%nfu), shapef%thetas(inter%meshn,shapef%nfu))

    do ifun = 1, shapef%nfu
      read(unit, fmt="(16i5)") lm
      CHECKASSERT(lm > 0)
      shapef%llmsp(ifun) = lm
      read(unit, fmt="(4d20.12)") shapef%thetas(1:inter%meshn,ifun)
    enddo ! ifun
    
  endsubroutine ! create

  elemental subroutine destroy_intermesh(inter)
    type(Intermesh), intent(inout) :: inter
    
    integer :: ist
    deallocate (inter%nm, inter%xrn, inter%drn, stat=ist)
  endsubroutine ! destroy
  
  elemental subroutine destroy_shapefunction(shapef)
    type(Shapefunction), intent(inout) :: shapef

    integer :: ist
    deallocate(shapef%llmsp, shapef%thetas, stat=ist)
  endsubroutine ! destroy
  
  elemental subroutine destroy_ShapefunFile(sfile)
    type(ShapefunFile), intent (inout) :: sfile

    integer :: icell, ist

    do icell = 1, sfile%ncell
      call destroy_intermesh(sfile%mesh(icell))
      call destroy_shapefunction(sfile%shapes(icell))
    enddo ! icell

    deallocate(sfile%mesh, sfile%shapes, stat=ist)
  endsubroutine ! destroy
  
endmodule read_formatted_shapefun_mod

#ifdef TEST_READ_FORMATTED_SHAPEFUN_MOD
program test_read_formatted
  use read_formatted_shapefun_mod, only: create_read_ShapefunFile, destroy_ShapefunFile
  implicit none

  type(ShapefunFile) :: sfile
  integer, parameter :: fu = 42
  integer :: ifun, icell

  open(fu, form='formatted', file='shapefun', action='read', status='old')
  call create_read_ShapefunFile(sfile, fu)

  do icell = 1, sfile%ncell

    write(*,*) "nm :", sfile%mesh(icell)%nm

    do ifun = 1, sfile%shapes(icell)%nfu
      write(*,*) "lm = ", sfile%shapes(icell)%llmsp(ifun)
      write(*,*) "-----------------------------------------------"
      write(*,*) sfile%shapes(icell)%thetas(:, ifun)
    enddo ! ifun

  enddo ! icell

  call destroy_ShapefunFile(sfile)
  close(fu)

endprogram
#endif
