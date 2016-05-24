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

  type KKRresults
    double precision, allocatable :: rMTref(:)
    double complex, allocatable :: tmatN(:,:,:)  
    double complex, allocatable :: dtdE(:,:,:)  
    double complex, allocatable :: trefLL(:,:,:)  
    double complex, allocatable :: dtrefLL(:,:,:)  
    double complex, allocatable :: dGrefN(:,:,:,:)
    double complex, allocatable :: GmatN(:,:,:,:)
    double complex, allocatable :: Lly_G0Tr(:)   
    double complex, allocatable :: Lly_Grdt(:,:)  
    double complex, allocatable :: Tr_alph(:)  
    integer :: noiter

    integer :: lmmaxd
    integer :: nspind
    integer :: naclsd
    integer :: iemxd
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
    type(KKRresults), intent(inout) :: self
    type(DimParams),  intent(in)    :: dims
    integer, intent(in)              :: naclsd

    call createKKRresultsImpl(self, dims%lmmaxd, dims%nspind, naclsd, dims%iemxd, dims%nguessd, dims%ekmd, dims%smpid)
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Constructs a KKRresults object.
  !> @param[inout] self    The KKRresults object to construct.
  !> @param[in]    LMMAXD
  !> @param[in]    nspind
  !> @param[in]    naclsd
  !> @param[in]    iemxd
  subroutine createKKRresultsImpl(self, lmmaxd, nspind, naclsd, iemxd, nguessd, ekmd, smpid)
    type(KKRresults), intent(inout) :: self
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: nspind
    integer, intent(in) :: naclsd
    integer, intent(in) :: iemxd
    integer, intent(in) :: nguessd
    integer, intent(in) :: ekmd
    integer, intent(in) :: smpid

    integer :: memory_stat
    double complex, parameter :: CZERO=(0.d0, 0.d0)

    self%lmmaxd = LMMAXD
    self%nspind = nspind
    self%naclsd = naclsd
    self%iemxd = iemxd
    self%nguessd = nguessd
    self%ekmd = ekmd
    self%smpid = smpid

    if (naclsd < 1) then
      write(*,*) "ERROR: Number of atoms in cluster <= 0", __FILE__, __LINE__
      stop
    endif

    ALLOCATECHECK(self%rMTref(naclsd))
    ALLOCATECHECK(self%tmatN(LMMAXD,LMMAXD,nspind))
    ALLOCATECHECK(self%dtdE(LMMAXD,LMMAXD,nspind))
    ALLOCATECHECK(self%trefLL(LMMAXD,LMMAXD,naclsd))
    ALLOCATECHECK(self%dtrefLL(LMMAXD,LMMAXD,naclsd))
    ALLOCATECHECK(self%dGrefN(LMMAXD,LMMAXD,naclsd,1))
    ALLOCATECHECK(self%GmatN(LMMAXD,LMMAXD,iemxd,nspind))
    ALLOCATECHECK(self%Lly_G0Tr(iemxd))
    ALLOCATECHECK(self%Lly_Grdt(iemxd,nspind))
    ALLOCATECHECK(self%Tr_alph(nspind))

    self%noiter = 0

    self%trefLL  = CZERO
    self%dtrefLL = CZERO

    ! use garbage values for initialisation
    self%GmatN    = dcmplx(9e9, 9e9)
    self%Lly_Grdt = dcmplx(9e9, 9e9)
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a KKRresults object.
  !> @param[inout] self    The KKRresults object to destroy.
  subroutine destroyKKRresults(self)
    type(KKRresults), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%rMTref)
    DEALLOCATECHECK(self%tmatN)
    DEALLOCATECHECK(self%dtdE)
    DEALLOCATECHECK(self%trefLL)
    DEALLOCATECHECK(self%dtrefLL)
    DEALLOCATECHECK(self%dGrefN)
    DEALLOCATECHECK(self%GmatN)
    DEALLOCATECHECK(self%Lly_G0Tr)
    DEALLOCATECHECK(self%Lly_Grdt)
    DEALLOCATECHECK(self%Tr_alph)

  endsubroutine ! destroy

endmodule
