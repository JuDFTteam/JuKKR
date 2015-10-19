module type_tbcluster
use nrtype
 TYPE                      ::      TBCLUSTER_TYPE

 INTEGER                   ::         NACLSD,              &    ! maximum number of atoms in a cluster
                                      NCLSD                     ! maximum number of clusters
 INTEGER                   ::         NCLS                     ! number of clusters
 REAL(kind=DP),allocatable ::         RCLS(:,:,:) !(:,:,:NATOM)       ! real space position of atom in cluster
 INTEGER,allocatable       ::         CLS(:),  &!(NATOM),             &    ! sort of cluster around atom
                                      NACLS(:), &!(NCLSD),        &    ! number of atoms in cluster
                                      NACLSIMPMAX(:),&                   ! number of impurity atoms in cluster
                                      ATOM(:,:), &!(NACLSD,NR)           ! ??????????
                                      NACLSIMP(:)
 INTEGER,allocatable       ::         IATOM2INDEX(:,:), & !(NATOM,NCLSD), &
                                      INDEX2IATOM(:,:), & !(NACLSD*LMAXD,NCLSD), &
                                      NCLSGF(:) !(NCLSD) 

 END TYPE TBCLUSTER_TYPE

end module type_tbcluster