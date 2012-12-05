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

  type KKRresults
    double complex , allocatable, dimension(:,:,:)  :: TMATN
    double complex , allocatable, dimension(:,:,:)  :: DTDE
    double complex , allocatable, dimension(:,:,:)  :: TREFLL
    double complex , allocatable, dimension(:,:,:)  :: DTREFLL
    double complex , allocatable, dimension(:,:,:,:)  :: DGREFN
    double complex , allocatable, dimension(:,:,:,:)  :: GREFN
    double complex , allocatable, dimension(:,:,:,:)  :: GMATN
    double complex , allocatable, dimension(:,:)  :: LLY_G0TR
    double complex , allocatable, dimension(:,:)  :: LLY_GRDT
    double complex , allocatable, dimension(:)  :: TR_ALPH
    integer  :: NOITER

    integer :: LMMAXD
    integer :: NSPIND
    integer :: NREFD
    integer :: NACLSD
    integer :: NCLSD
    integer :: IEMXD
  end type KKRresults

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a KKRresults object.
  !> @param[inout] self    The KKRresults object to construct.
  !> @param[in]    dims
  subroutine createKKRresults(self, dims)
    use DimParams_mod
    implicit none
    type (KKRresults), intent(inout) :: self
    type (DimParams),  intent(in)    :: dims

    call createKKRresultsImpl(self, dims%LMMAXD, dims%NSPIND, &
                              dims%NREFD, dims%NACLSD, dims%NCLSD, dims%IEMXD)
  end subroutine

  !-----------------------------------------------------------------------------
  !> Constructs a KKRresults object.
  !> @param[inout] self    The KKRresults object to construct.
  !> @param[in]    LMMAXD
  !> @param[in]    NSPIND
  !> @param[in]    NREFD
  !> @param[in]    NACLSD
  !> @param[in]    NCLSD
  !> @param[in]    IEMXD
  subroutine createKKRresultsImpl(self, LMMAXD,NSPIND,NREFD,NACLSD,NCLSD,IEMXD)
    implicit none
    type (KKRresults), intent(inout) :: self
    integer, intent(in) ::  LMMAXD
    integer, intent(in) ::  NSPIND
    integer, intent(in) ::  NREFD
    integer, intent(in) ::  NACLSD
    integer, intent(in) ::  NCLSD
    integer, intent(in) ::  IEMXD

    integer :: memory_stat
    double complex, parameter :: CZERO = (0.0d0, 0.0d0)

    self%LMMAXD = LMMAXD
    self%NSPIND = NSPIND
    self%NREFD = NREFD
    self%NACLSD = NACLSD
    self%NCLSD = NCLSD
    self%IEMXD = IEMXD

    ALLOCATECHECK(self%TMATN(LMMAXD,LMMAXD,NSPIND))
    ALLOCATECHECK(self%DTDE(LMMAXD,LMMAXD,NSPIND))
    ALLOCATECHECK(self%TREFLL(LMMAXD,LMMAXD,NREFD))
    ALLOCATECHECK(self%DTREFLL(LMMAXD,LMMAXD,NREFD))
    ALLOCATECHECK(self%DGREFN(LMMAXD,LMMAXD,NACLSD,NCLSD))
    ALLOCATECHECK(self%GREFN(LMMAXD,LMMAXD,NACLSD,NCLSD))
    ALLOCATECHECK(self%GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND))
    ALLOCATECHECK(self%LLY_G0TR(NCLSD,IEMXD))
    ALLOCATECHECK(self%LLY_GRDT(IEMXD,NSPIND))
    ALLOCATECHECK(self%TR_ALPH(NSPIND))

    self%noiter = 0

    self%TREFLL = CZERO
    self%DTREFLL = CZERO

    ! use garbage values for initialisation
    self%GMATN = dcmplx(99999.0d0, 99999.0d0)
    self%LLY_GRDT = dcmplx(99999.0d0, 99999.0d0)
  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a KKRresults object.
  !> @param[inout] self    The KKRresults object to destroy.
  subroutine destroyKKRresults(self)
    implicit none
    type (KKRresults), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%TMATN)
    DEALLOCATECHECK(self%DTDE)
    DEALLOCATECHECK(self%TREFLL)
    DEALLOCATECHECK(self%DTREFLL)
    DEALLOCATECHECK(self%DGREFN)
    DEALLOCATECHECK(self%GREFN)
    DEALLOCATECHECK(self%GMATN)
    DEALLOCATECHECK(self%LLY_G0TR)
    DEALLOCATECHECK(self%LLY_GRDT)
    DEALLOCATECHECK(self%TR_ALPH)

  end subroutine

end module
