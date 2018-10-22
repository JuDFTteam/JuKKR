!------------------------------------------------------------------------------------
!> Summary: Energy type
!> Author: 
!> Category: KKRimp, total-energy
!> Deprecated: False 
!> Stores the different parts of the energy
!------------------------------------------------------------------------------------
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
      REAL(kind=8),allocatable  ::    EXC(:,:)     !! exchange-correlation energy
      REAL(kind=8),allocatable  ::    ECOU(:,:)    !! Coulomb energy
      REAL(kind=8),allocatable  ::    EPOTIN(:)    !! energy of the input potential (epotinb)
      REAL(kind=8),allocatable  ::    ESPC(:,:,:)  !! single particle core energy 
      REAL(kind=8),allocatable  ::    ESPV(:,:,:)  !! single particle valence energy (changed for the relativistic case)
 END TYPE ENERGYPARTS_TYPE

end module type_ENERGYPARTS
