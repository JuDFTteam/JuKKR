module type_ldau
use nrtype
 TYPE                               ::  LDAU_TYPE
   INTEGER                          ::  IDOLDAU
   INTEGER                          ::  LOPT                      ! Angular momentum value to apply lda+u
   INTEGER                          ::  IELDAUSTART,IELDAUEND     ! Energies between which lda+u should be applied
   REAL(KIND=DP)                    ::  UEFF,JEFF                 ! Effective U and J
   REAL(KIND=DP)                    ::  EREFLDAU                  ! Energy for calculation of basis set phi
   REAL(KIND=DP)                    ::  EU,EDC                    ! Total-energy terms (single-particle and double-counting)
   REAL(KIND=DP),ALLOCATABLE        ::  CUTOFF(:)                 ! Cutoff for applying WLDAU in radial mesh (nrmax)
   REAL(KIND=DP),ALLOCATABLE        ::  ULDAU(:,:,:,:)            ! Coulomb matrix (mmax,mmax,mmax,mmax)
   REAL(KIND=DP),ALLOCATABLE        ::  WLDAU(:,:,:)              ! Interaction potential (mmax,mmax,nspin) 
   REAL(KIND=DP)                    ::  WLDAUAV(2)                ! Averaged interaction potential (nspin)
   COMPLEX(KIND=DP),ALLOCATABLE     ::  PHI(:)                    ! Basis function for LDA+U (nrmax)         
   COMPLEX(KIND=DPC),ALLOCATABLE    ::  DENMATC(:,:,:)            ! Density matrix  (mmax,mmax,nspin)
 END TYPE LDAU_TYPE

end module type_ldau
