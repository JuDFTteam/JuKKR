module type_wavefunction
use nrtype
 TYPE                              ::  wavefunction_TYPE
       INTEGER                     :: lmsize,lmsize2,nrmaxnew
       INTEGER                     :: NVEC
       DOUBLE COMPLEX,allocatable  ::  SLL(:,:,:,:), RLL(:,:,:,:)
       DOUBLE COMPLEX,allocatable  ::  SLLleft(:,:,:,:), RLLleft(:,:,:,:)

       INTEGER                     :: deallocate
       INTEGER                     :: rll_saved
       INTEGER                     :: sll_saved
       INTEGER                     :: rllleft_saved
       INTEGER                     :: sllleft_saved



 END TYPE wavefunction_TYPE

end module type_wavefunction