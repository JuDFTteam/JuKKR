!------------------------------------------------------------------------------------
!> Summary: NOT USED
!> Author: 
!> Category: KKRimp
!> Deprecated: True
!> Non-sense type, which is luckily not used
!------------------------------------------------------------------------------------
module type_cell
 TYPE                              ::  CELL_ENERGYPARTS
C ----------------------------------------------------------------------
C   ECOU(0:LPOTD,NATYPD)    ! Coulomb energy                    
C   EPOTIN(NATYPD),         ! energy of input potential (EPOTINB
C   ESPC(0:3,NPOTD),        ! energy single particle core       
C   ESPV(0:LMAXD1,NPOTD)    ! energy single particle valence    
C                           ! both changed for the relativistic 
C                           ! case
C   EXC(0:LPOTD,NATYPD),    ! E_xc
C ----------------------------------------------------------------------
      REAL(kind=8),allocatable  ::    EXC(:,:)
!       REAL(kind=8)  :: ECOU(0:LPOTD,NATYPD),EPOTIN(NATYPD),      
!                        ESPC(0:3,NPOTD),ESPV(0:LMAXD1,NPOTD),
!                        EXC(0:LPOTD,NATYPD)

 END TYPE CELL_TYPE

end module type_cell
