!-------------------------------------------------------------------------------
!> Summary: Type holding wavefunctions and corresponding array dimensions and save-flags
!> Author: 
!> Category: KKRimp, 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_wavefunction

  type :: wavefunction_type

    integer                     :: lmsize,lmsize2,nrmaxnew
    integer                     :: nvec
    ! MdSD: sll never seems to be used, better skip their allocation
    double complex,allocatable  :: sll(:,:,:,:), rll(:,:,:,:), ull(:,:,:,:)
    double complex,allocatable  :: sllleft(:,:,:,:), rllleft(:,:,:,:)

    integer                     :: deallocate = 0
    integer                     :: rll_saved
    integer                     :: sll_saved
    integer                     :: rllleft_saved
    integer                     :: sllleft_saved

  end type wavefunction_type

end module type_wavefunction
