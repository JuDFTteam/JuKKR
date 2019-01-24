!> This module defines a struct that contains the results of the scattering calculations
!> for ONE SPECIFIC atom.
!>
!> @author Elias Rabel

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)

module KKRresults_mod
  implicit none
  private
  public :: KKRresults, create, destroy

  type KKRresults
    double precision, allocatable :: rMTref(:)
    double complex, allocatable :: TmatN(:,:,:)  
    double complex, allocatable :: dTmatN(:,:,:) !< energy derivative of TmatN
    double complex, allocatable :: tref_ell(:,:) !< dim(0:lmax,1:)
    double complex, allocatable :: dtref_ell(:,:)  !< dim(0:lmax,1:)
    double complex, allocatable :: dGrefN(:,:,:,:)
    double complex, allocatable :: GmatN(:,:,:,:)
    double complex, allocatable :: Lly_G0Tr(:)   
    double complex, allocatable :: Lly_Grdt(:,:)  
    double complex, allocatable :: Tr_alph(:)  

    integer :: lmmaxd
    integer :: lmmaxd_noco ! NOCO
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
    integer, intent(in)             :: naclsd

    call createKKRresultsImpl(self, dims%lmaxd, dims%lmmaxd, dims%korbit, dims%nspind, naclsd, dims%iemxd, dims%nguessd, dims%ekmd, dims%smpid)
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Constructs a KKRresults object.
  !> @param[inout] self    The KKRresults object to construct.
  !> @param[in]    LMMAXD
  !> @param[in]    nspind
  !> @param[in]    naclsd
  !> @param[in]    iemxd
  subroutine createKKRresultsImpl(self, lmaxd, lmmaxd, korbit, nspind, naclsd, iemxd, nguessd, ekmd, smpid)
    type(KKRresults), intent(inout) :: self
    integer, intent(in) :: lmaxd
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: korbit ! NOCO
    integer, intent(in) :: nspind
    integer, intent(in) :: naclsd
    integer, intent(in) :: iemxd
    integer, intent(in) :: nguessd
    integer, intent(in) :: ekmd
    integer, intent(in) :: smpid

    integer :: memory_stat
    integer :: lmaxd_noco, lmmaxd_noco ! NOCO
    double complex, parameter :: CZERO=(0.d0, 0.d0)

    lmmaxd_noco = lmmaxd*(korbit+1) ! NOCO, matrix size must be doubled if korbit==1

    self%lmmaxd = lmmaxd
    self%lmmaxd_noco = lmmaxd_noco ! NOCO
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
    ALLOCATECHECK(self%TmatN(lmmaxd_noco,lmmaxd_noco,nspind))
    ALLOCATECHECK(self%dTmatN(lmmaxd_noco,lmmaxd_noco,nspind))
    ALLOCATECHECK(self%tref_ell(0:lmaxd,naclsd))
    ALLOCATECHECK(self%dtref_ell(0:lmaxd,naclsd))
    ALLOCATECHECK(self%dGrefN(lmmaxd,lmmaxd,naclsd,1))
    ALLOCATECHECK(self%GmatN(lmmaxd_noco,lmmaxd_noco,iemxd,nspind))
    ALLOCATECHECK(self%Lly_G0Tr(iemxd))
    ALLOCATECHECK(self%Lly_Grdt(iemxd,nspind))
    ALLOCATECHECK(self%Tr_alph(nspind))

    ! use garbage values for initialisation
    self%GmatN    = dcmplx(9e9, 9e9)
    self%Lly_Grdt = dcmplx(9e9, 9e9)
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a KKRresults object.
  !> @param[inout] self    The KKRresults object to destroy.
  elemental subroutine destroyKKRresults(self)
    type(KKRresults), intent(inout) :: self

    integer :: ist

    deallocate(self%rMTref, stat=ist)
    deallocate(self%TmatN, stat=ist)
    deallocate(self%dTmatN, stat=ist)
    deallocate(self%tref_ell, stat=ist)
    deallocate(self%dtref_ell, stat=ist)
    deallocate(self%dGrefN, stat=ist)
    deallocate(self%GmatN, stat=ist)
    deallocate(self%Lly_G0Tr, stat=ist)
    deallocate(self%Lly_Grdt, stat=ist)
    deallocate(self%Tr_alph, stat=ist)

  endsubroutine ! destroy

endmodule ! KKRresults_mod
