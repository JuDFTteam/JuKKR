module mod_calctmatfull
  contains

SUBROUTINE CALCTMATFULL(ERYD,VPOT1,VPOT2,CELL,ZATOM,LMAXATOM,TMATLL,config,nspin,lmaxd,iatom,inpsusc)! & ! susc |||| 'iatom,inpsusc' added by Benedikt (2014/06)
  use mod_calctmat
  use type_cell
  use type_config
  use type_inpsusc ! susc |||| added by Benedikt (2014/06)
  implicit none
!interface variables
  complex(kind=dpc),intent(in)              ::  eryd
  type(cell_type),intent(in)                ::  cell
  real(kind=dp),intent(in)                  ::  vpot1(1:cell%nrmaxd,(2*lmaxatom+1)**2)
  real(kind=dp),intent(in)                  ::  vpot2(1:cell%nrmaxd,(2*lmaxatom+1)**2)
  real(kind=dp),intent(in)                  ::  zatom
  integer,intent(in)                        ::  lmaxatom
  complex(kind=dpc)                         ::  tmatll(2*(lmaxatom+1)**2,2*(lmaxatom+1)**2)
  type(config_type),intent(in)              ::  config
  integer,intent(in)                        ::  nspin
  integer,intent(in)                        ::  iatom ! susc |||| 'iatom' added by Benedikt (2014/06)
!   integer,intent(in)                        ::  natom
  integer,intent(in)                        ::  lmaxd


  type(inpsusc_type),intent(inout)          :: inpsusc   ! susc |||| added by Benedikt (2014/06)

  integer                       ::  ispin,lmmax,bound1,bound2
  complex(kind=dpc)                         ::  tmatll_temp((lmaxatom+1)**2,(lmaxatom+1)**2)

lmmax=(lmaxatom+1)**2
! write(*,*) 'lmmax',lmmax
do ispin=1,nspin
bound1=lmmax*(ispin-1)+1
! write(*,*) 'bound1',bound1

bound2=lmmax*ispin
! write(*,*) 'bound2',bound2

! write(*,*) 'CALCTMAT'
 if (ispin==1) then
        CALL  CALCTMAT(ERYD,VPOT1,CELL,ZATOM,lmaxatom, &
                       tmatll_temp,config,ispin,nspin,iatom,inpsusc ) !  & ! susc |||| 'iatom,inpsusc' added by Benedikt (2014/06)
! write(*,*) 'end CALCTMAT'
elseif (ispin==2) then
        CALL  CALCTMAT(ERYD,VPOT2,CELL,ZATOM,lmaxatom, &
                       tmatll_temp,config,ispin,nspin,iatom,inpsusc ) !  & ! susc |||| 'iatom,inpsusc' added by Benedikt (2014/06)
else 
stop '[CALCTMATFULL] error'
end if

         tmatll(bound1:bound2,bound1:bound2)=tmatll_temp

! write(*,*) 'Ctt'


end do

! write(*,*) 'Ctt'

end subroutine calctmatfull

end module mod_calctmatfull
