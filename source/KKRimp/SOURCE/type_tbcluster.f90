!-------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> Category: KKRimp, 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_tbcluster
use nrtype, only: dp
 TYPE                      :: TBCLUSTER_TYPE

 INTEGER                   :: NACLSD !! maximum number of atoms in a cluster
 INTEGER                   :: NCLSD !! maximum number of clusters
 INTEGER                   :: NCLS !! number of clusters
 REAL(kind=DP),allocatable :: RCLS(:,:,:) !! real space position of atom in cluster (:,:,:NATOM)
 INTEGER,allocatable       :: CLS(:) !! sort of cluster around atom (NATOM) 
 INTEGER,allocatable       :: NACLS(:) !! number of atoms in cluster(NCLSD) 
 INTEGER,allocatable       :: NACLSIMPMAX(:) !! number of impurity atoms in cluster
 INTEGER,allocatable       :: ATOM(:,:) !! (NACLSD,NR)
 INTEGER,allocatable       :: NACLSIMP(:)
 INTEGER,allocatable       :: IATOM2INDEX(:,:) !!(NATOM,NCLSD)
 INTEGER,allocatable       :: INDEX2IATOM(:,:) !!(NACLSD*LMAXD,NCLSD)
 INTEGER,allocatable       :: NCLSGF(:) !!(NCLSD) 

 END TYPE TBCLUSTER_TYPE

end module type_tbcluster