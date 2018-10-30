module type_tmat
use nrtype
 TYPE                              ::  TMAT_TYPE
   DOUBLE COMPLEX,ALLOCATABLE              ::  TMAT(:,:)
   DOUBLE COMPLEX,ALLOCATABLE              ::  deltaT_Jij(:,:,:)

 END TYPE TMAT_TYPE

end module type_tmat