MODULE MOD_RHOVALfull
  CONTAINS
      SUBROUTINE RHOVALFULL(IE,IELAST, ERYD ,WEZ, GMATLL, NSPIN, &
                        IATOM,CELL,VPOT_UP,VPOT_DN,SHAPEFUN,GAUNTCOEFF, ZATOM, DENSITY, &
                         LMAXATOM,LMMAXATOM,config,lmaxd,energyparts,efermi)
      use mod_rhoval
      use type_cell
      use type_gauntcoeff
      use type_shapefun
      use type_density
      use type_config
      use type_energyparts
      implicit none
!interface variables
      integer,intent(in)                    ::   lmaxd
      integer,intent(in)                    ::   ie
      integer                               ::   ielast
      complex(kind=dpc),intent(in)          ::   eryd
      complex(kind=dpc),intent(in)          ::   wez
      complex(kind=dpc),intent(in)          ::   gmatll(2*lmmaxatom,2*lmmaxatom)
      integer,intent(in)                    ::   nspin
      integer,intent(in)                    ::   iatom
      type(cell_type),intent(in)            ::   cell
      real(kind=dp),intent(in)              ::   vpot_up(1:cell%nrmaxd,(2*lmaxd+1)**2)
      real(kind=dp),intent(in)              ::   vpot_dn(1:cell%nrmaxd,(2*lmaxd+1)**2)
      type(shapefun_type),intent(in)        ::   shapefun
      type(gauntcoeff_type),intent(in)      ::   gauntcoeff
      real(kind=dp),intent(in)              ::   zatom
      type(density_type)                    ::   density
!       complex(kind=dpc)                     ::   den(0:LMAXatom+1,300) !(0:lmaxatom+1) !,ez(iemxd), &
!       real(kind=dp)                         ::   rho2ns(cell%nrmax,(2*lmaxatom+1)**2,2)
!       real(kind=dp)                         ::   espv(0:lmaxatom+1,2)
      integer,intent(in)                    ::   lmaxatom
      integer,intent(in)                    ::   lmmaxatom
      type(config_type)                    ::   config
      type(energyparts_type)                    ::   energyparts
      double precision              ::  efermi


       complex(kind=dpc)                     ::   gmatll_temp(lmmaxatom,lmmaxatom)
      integer                               ::   ispin,lmmax,bound1,bound2
      
lmmax=lmmaxatom
! write(*,*) 'lmmax',lmmax
do ispin=1,nspin
  bound1=lmmax*(ispin-1)+1 
!   write(*,*) 'bound1',bound1
  
  bound2=lmmax*ispin
!   write(*,*) 'bound2',bound2

  gmatll_temp=gmatll(bound1:bound2,bound1:bound2)

if (ispin==1) then
  call RHOVAL(ie,ielast,eryd ,WEZ, gmatll_temp,ISPIN, NSPIN, &
                  IATOM,CELL,VPOT_UP(:,:),SHAPEFUN, GAUNTCOEFF, ZATOM,&
                  density, LMAXATOM, LMMAXATOM,config ,lmaxd,energyparts,efermi)
elseif (ispin==2) then
  call RHOVAL(ie,ielast,eryd ,WEZ, gmatll_temp,ISPIN, NSPIN, &
                  IATOM,CELL,VPOT_DN(:,:),SHAPEFUN, GAUNTCOEFF, ZATOM,&
                  density, LMAXATOM, LMMAXATOM,config ,lmaxd,energyparts,efermi)
else
  stop '[rhovalfull] ispin error'
end if
end do !ispin


end subroutine rhovalfull

end module mod_rhovalfull
