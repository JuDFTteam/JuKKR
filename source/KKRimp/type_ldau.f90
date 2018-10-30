!-------------------------------------------------------------------------------
!> Summary: Type holding LDA+U information
!> Author: 
!> Category: KKRimp, xc-potential
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_ldau

  use nrtype, only: dp, dpc

  type :: ldau_type
    integer                       :: idoldau               !! switch turning lda+u on/off (1,0)
    integer                       :: lopt                  !! angular momentum value to apply lda+u
    integer                       :: ieldaustart,ieldauend !! energies between which lda+u should be applied
    real(kind=dp)                 :: ueff,jeff             !! effective u and j
    real(kind=dp)                 :: erefldau              !! energy for calculation of basis set phi
    real(kind=dp)                 :: eu,edc                !! total-energy terms (single-particle and double-counting)
    real(kind=dp),allocatable     :: cutoff(:)             !! cutoff for applying wldau in radial mesh (nrmax)
    real(kind=dp),allocatable     :: uldau(:,:,:,:)        !! coulomb matrix (mmax,mmax,mmax,mmax)
    real(kind=dp),allocatable     :: wldau(:,:,:)          !! interaction potential (mmax,mmax,nspin) 
    real(kind=dp)                 :: wldauav(2)            !! averaged interaction potential (nspin)
    complex(kind=dp),allocatable  :: phi(:)                !! basis function for lda+u (nrmax)         
    complex(kind=dpc),allocatable :: denmatc(:,:,:)        !! density matrix (mmax,mmax,nspin)
  end type ldau_type

end module type_ldau
