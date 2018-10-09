!-------------------------------------------------------------------------------
!> Summary: Type holding screening (tight-binding) cluster information
!> Author: 
!> Category: KKRimp, reference-system
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_tbcluster

  use nrtype, only: dp

  type                      :: tbcluster_type

  integer                   :: naclsd           !! maximum number of atoms in a cluster
  integer                   :: nclsd            !! maximum number of clusters
  integer                   :: ncls             !! number of clusters
  real(kind=dp),allocatable :: rcls(:,:,:)      !! real space position of atom in cluster (:,:,:natom)
  integer,allocatable       :: cls(:)           !! sort of cluster around atom (natom) 
  integer,allocatable       :: nacls(:)         !! number of atoms in cluster(nclsd) 
  integer,allocatable       :: naclsimpmax(:)   !! number of impurity atoms in cluster
  integer,allocatable       :: atom(:,:)        !! (naclsd,nr)
  integer,allocatable       :: naclsimp(:)      !!
  integer,allocatable       :: iatom2index(:,:) !!(natom,nclsd)
  integer,allocatable       :: index2iatom(:,:) !!(naclsd*lmaxd,nclsd)
  integer,allocatable       :: nclsgf(:)        !!(nclsd) 

  end type tbcluster_type

end module type_tbcluster