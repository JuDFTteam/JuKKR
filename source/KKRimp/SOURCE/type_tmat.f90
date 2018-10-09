!-------------------------------------------------------------------------------
!> Summary: Data type holding single-site tmatix and Delta-t-matrix (for Jijs)
!> Author: 
!> Category: KKRimp, single-site
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_tmat

  type ::  tmat_type
    double complex, allocatable :: tmat(:,:)          !! single-site t-matrix
    double complex, allocatable :: deltat_jij(:,:,:)  !! \[\Delta t\] matrix for Jij calculation
  end type tmat_type

end module type_tmat