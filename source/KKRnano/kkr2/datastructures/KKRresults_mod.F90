!> This module defines a struct that contains the results of the scattering calculations
!> for ONE SPECIFIC atom.
!>
!> @author Elias Rabel

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module KKRresults_mod
  implicit none
  private
  public :: KKRresults, create, destroy
  public :: createKKRresults, destroyKKRresults ! deprecated

  type KKRresults
    double complex, allocatable :: TMATN(:,:,:)  
    double complex, allocatable :: DTDE(:,:,:)  
    double complex, allocatable :: TREFLL(:,:,:)  
    double complex, allocatable :: DTREFLL(:,:,:)  
    double complex, allocatable :: DGREFN(:,:,:,:)
    double complex, allocatable :: GMATN(:,:,:,:)
    double complex, allocatable :: LLY_G0TR(:)   
    double complex, allocatable :: LLY_GRDT(:,:)  
    double complex, allocatable :: TR_ALPH(:)  
    integer  :: NOITER

    integer :: LMMAXD
    integer :: NSPIND
    integer :: NACLSD
    integer :: IEMXD
    integer :: ekmd
    integer :: nguessd
    integer :: smpid
  endtype ! KKRresults

  interface create
    module procedure createKKRresults
  endinterface
  
  interface destroy
    module procedure destroyKKRresults
  endinterface
  
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a KKRresults object.
  !> @param[inout] self    The KKRresults object to construct.
  !> @param[in]    dims
  !> @param[in]    naclsd  maximal number of cluster atoms
  subroutine createKKRresults(self, dims, naclsd)
    use DimParams_mod, only: DimParams
    type (KKRresults), intent(inout) :: self
    type (DimParams),  intent(in)    :: dims
    integer, intent(in)              :: naclsd

    call createKKRresultsImpl(self, dims%LMMAXD, dims%NSPIND, NACLSD, dims%IEMXD, dims%nguessd, dims%ekmd, dims%smpid)
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Constructs a KKRresults object.
  !> @param[inout] self    The KKRresults object to construct.
  !> @param[in]    LMMAXD
  !> @param[in]    NSPIND
  !> @param[in]    NACLSD
  !> @param[in]    IEMXD
  subroutine createKKRresultsImpl(self, LMMAXD,NSPIND,NACLSD,IEMXD, nguessd, ekmd, smpid)
    type (KKRresults), intent(inout) :: self
    integer, intent(in) ::  LMMAXD
    integer, intent(in) ::  NSPIND
    integer, intent(in) ::  NACLSD
    integer, intent(in) ::  IEMXD
    integer, intent(in) ::  nguessd
    integer, intent(in) ::  ekmd
    integer, intent(in) ::  smpid

    integer :: memory_stat
    double complex, parameter :: CZERO=(0.d0, 0.d0)

    self%LMMAXD = LMMAXD
    self%NSPIND = NSPIND
    self%NACLSD = NACLSD
    self%IEMXD = IEMXD
    self%nguessd = nguessd
    self%ekmd = ekmd
    self%smpid = smpid

    if (naclsd < 1) then
      write(*,*) "ERROR: Number of atoms in cluster <= 0", __FILE__, __LINE__
      stop
    endif

    ALLOCATECHECK(self%TMATN(LMMAXD,LMMAXD,NSPIND))
    ALLOCATECHECK(self%DTDE(LMMAXD,LMMAXD,NSPIND))
    ALLOCATECHECK(self%TREFLL(LMMAXD,LMMAXD,naclsd))
    ALLOCATECHECK(self%DTREFLL(LMMAXD,LMMAXD,naclsd))
    ALLOCATECHECK(self%DGREFN(LMMAXD,LMMAXD,NACLSD,1))
    ALLOCATECHECK(self%GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND))
    ALLOCATECHECK(self%LLY_G0TR(IEMXD))
    ALLOCATECHECK(self%LLY_GRDT(IEMXD,NSPIND))
    ALLOCATECHECK(self%TR_ALPH(NSPIND))

    self%noiter = 0

    self%TREFLL = CZERO
    self%DTREFLL = CZERO

    ! use garbage values for initialisation
    self%GMATN = dcmplx(99999.0d0, 99999.0d0)
    self%LLY_GRDT = dcmplx(99999.0d0, 99999.0d0)
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a KKRresults object.
  !> @param[inout] self    The KKRresults object to destroy.
  subroutine destroyKKRresults(self)
    type (KKRresults), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%TMATN)
    DEALLOCATECHECK(self%DTDE)
    DEALLOCATECHECK(self%TREFLL)
    DEALLOCATECHECK(self%DTREFLL)
    DEALLOCATECHECK(self%DGREFN)
    DEALLOCATECHECK(self%GMATN)
    DEALLOCATECHECK(self%LLY_G0TR)
    DEALLOCATECHECK(self%LLY_GRDT)
    DEALLOCATECHECK(self%TR_ALPH)

  endsubroutine ! destroy

endmodule
