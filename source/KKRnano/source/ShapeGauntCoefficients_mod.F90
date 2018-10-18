!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module ShapeGauntCoefficients_mod
!-------------------------------------------------------------------------------
!> Summary: Calculates the Gaunt coefficients for shape functions
!> Author: Elias Rabel, Paul F Baumeister
!> Category: KKRnano, shape-functions, special-functions
!>
!> This separated datastructure was necessary since it seems that a different
!> convention for the Gaunt coefficients was used in shape function related code.
!> Some macros for checked allocation/deallocation they need an integer variable 
!> named memory_stat declared in each routine they are used.
!-------------------------------------------------------------------------------
#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)
  implicit none
  private
  public :: ShapeGauntCoefficients, create, destroy

  type ShapeGauntCoefficients
    double precision, allocatable :: GSH(:)
    integer, allocatable :: ILM(:,:)
    integer, allocatable :: IMAXSH(:)
    integer :: NGSHD
    integer :: lmax
    integer :: lmpotd
  endtype
  
  interface create
    module procedure createShapeGauntCoefficients
  endinterface
  
  interface destroy
    module procedure destroyShapeGauntCoefficients
  endinterface

  contains

  !----------------------------------------------------------------------------
  subroutine createShapeGauntCoefficients(self, lmax)
    use Harmonics_mod, only: Gaunt2 ! initialization of wg and yrg
    use Harmonics_mod, only: SHAPEG_count
    use Harmonics_mod, only: SHAPEG
    type(ShapeGauntCoefficients), intent(inout) :: self
    integer, intent(in) :: lmax
    !---------------------------
    integer :: memory_stat
    integer :: LASSLD
    integer :: LPOT
    integer :: LMPOTD
    integer :: NGSHD
    double precision, allocatable :: WG(:) 
    double precision, allocatable :: YRG(:,:,:)

    LPOT = 2*lmax
    LASSLD = 4*lmax

    LMPOTD = (LPOT+1) ** 2

    self%lmax = lmax
    self%lmpotd = LMPOTD

    ALLOCATECHECK(WG(LASSLD))
    ALLOCATECHECK(YRG(LASSLD,0:LASSLD,0:LASSLD))

    call GAUNT2(WG, YRG, lmax)
    call SHAPEG_count(LPOT, WG, YRG, LMAX, NGSHD) ! determine number of coefficients

    ALLOCATECHECK(self%GSH(NGSHD))
    ALLOCATECHECK(self%ILM(NGSHD,3))
    ALLOCATECHECK(self%IMAXSH(0:LMPOTD))

    self%GSH = 0.0d0
    self%ILM = -1
    self%IMAXSH = -1

    call SHAPEG(LPOT, self%GSH, self%ILM, self%IMAXSH, WG, YRG, LMAX, NGSHD)

    self%NGSHD = NGSHD

    DEALLOCATECHECK(WG)
    DEALLOCATECHECK(YRG)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyShapeGauntCoefficients(self)
    type(ShapeGauntCoefficients), intent(inout) :: self
    integer :: ist ! ignore status
    deallocate(self%GSH, self%ILM, self%IMAXSH, stat=ist)
  endsubroutine ! destroy

endmodule ! ShapeGauntCoefficients_mod
