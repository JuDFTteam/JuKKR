module type_ENERGYPARTS
 TYPE                              ::  ENERGYPARTS_TYPE
! C ----------------------------------------------------------------------
! C   ECOU(0:LPOTD,NATYPD)    ! Coulomb energy                    
! C   EPOTIN(NATYPD),         ! energy of input potential (EPOTINB
! C   ESPC(0:3,NPOTD),        ! energy single particle core       
! C   ESPV(0:LMAXD1,NPOTD)    ! energy single particle valence    
! C                           ! both changed for the relativistic 
! C                           ! case
! C   EXC(0:LPOTD,NATYPD),    ! E_xc
! C ----------------------------------------------------------------------
      REAL(kind=8),allocatable  ::    EXC(:,:)
      REAL(kind=8),allocatable  ::    ECOU(:,:) !,EPOTIN(NATYPD),      
      REAL(kind=8),allocatable  ::    EPOTIN(:)
      REAL(kind=8),allocatable  ::    ESPC(:,:,:)
      REAL(kind=8),allocatable  ::    ESPV(:,:,:)
 END TYPE ENERGYPARTS_TYPE

end module type_ENERGYPARTS