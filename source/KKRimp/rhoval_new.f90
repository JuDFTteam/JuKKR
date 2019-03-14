module mod_rhoval_new

contains
!-------------------------------------------------------------------------
!> Summary: Main driver for valence charge density, new solver
!> Category: physical-observables, KKRimp
!>
!-------------------------------------------------------------------------
subroutine rhoval_new(eryd,ie,wez,cellnew,wavefunction,cell,gmatll,iatom,ispin,nspin,SHAPEFUN,GAUNTCOEFF, ZATOM, DENSITY, &
                         LMAXATOM,LMMAXATOM,config,lmaxd,energyparts,kspinorbit,use_fullgmat,nspinden,efermi, &
                         ldau) ! lda+u
use nrtype, only: dp, dpc
use type_cell, only: cell_type
use type_cellnew, only: cell_typenew
use type_shapefun, only: shapefun_type
use type_gauntcoeff, only: gauntcoeff_type
use type_density, only: density_type
use type_config, only: config_type
use type_wavefunction, only: wavefunction_type
use type_energyparts, only: energyparts_type
use mod_cheb2oldgrid, only: cheb2oldgrid
use mod_checknan, only: checknan
use mod_rllsll, only: rllsll
use mod_rhooutnew, only: rhooutnew
use mod_config, only: config_testflag
use mod_physic_params, only: cvlight
use mod_config, only: config_testflag
use type_ldau, only: ldau_type                          ! lda+u
use mod_intcheb_cell, only: intcheb_cell   ! lda+u

implicit none

!interface
complex(kind=dpc)                         :: eryd
integer                                   :: ie
complex(kind=dpc)                         :: wez
type(cell_typenew)                        :: cellnew
type(wavefunction_type)                   :: wavefunction
type(cell_type),intent(in)                :: cell
complex(kind=dpc),intent(in)              :: gmatll(:,:)
integer                                   :: iatom
integer                                   :: ispin
integer                                   :: nspin
type(shapefun_type),intent(in)            :: shapefun
type(gauntcoeff_type),intent(in)          :: gauntcoeff
real(kind=dp),intent(in)                  :: zatom
type(density_type)                        :: density
integer                                   :: lmaxatom
integer,intent(in)                        :: lmmaxatom
type(config_type)                         :: config
integer                                   :: lmaxd
type(energyparts_type)                    :: energyparts
integer                                   :: kspinorbit
integer                                   :: use_fullgmat
integer                                   :: nspinden
double precision                          :: efermi
type(ldau_type)                           :: ldau                       ! lda+u variables

!local
double complex,allocatable                :: rho2ns_complex(:,:,:)
double complex,allocatable                :: rho2ns_complex_temp(:)
!real(kind=dp)                             ::   espv(0:lmaxatom+1,2)
complex(kind=dpc)                         ::   df,ek
double complex                            :: rho2ns_integrated(4),temp1
integer                                   :: lval,imt1

double complex                            :: cden(cellnew%nrmaxnew,0:lmaxd,nspinden),cdenns(cellnew%nrmaxnew,nspinden), &
                                             cdenlm(cellnew%nrmaxnew,lmmaxatom,nspinden)  ! lm-dos
double complex,allocatable                :: gflle_part(:,:) ! lda+u

integer                                   :: lm1,ialpha, mval
integer                                   :: jspin
integer                                   :: lmslo,lmshi            ! lda+u
integer                                   :: nspinstart, nspinstop
! double precision              :: C0LL,fac

!       C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))

if (use_fullgmat==1) then
!   lmsize=2*lmmaxatom
!   nspinden=4
  nspinstart=1
  nspinstop=nspinden
else
!   nspinden=nspin
!   lmsize=lmmaxatom
  nspinstart=ispin
  nspinstop=ispin
end if

allocate(rho2ns_complex(cellnew%nrmaxnew,(2*lmaxatom+1)**2,nspinden))
allocate(rho2ns_complex_temp(cellnew%nrmaxnew))
allocate(gflle_part(wavefunction%lmsize,wavefunction%lmsize)) ! lda+u

if (ubound(gmatll,1)/=wavefunction%lmsize ) then
  stop 'error in rhovalnew'
end if

rho2ns_complex=(0.0d0,0.0d0)

df = wez/dble(nspin)


!#######################################################
! Define the prefactor for the greens function
! THE DEFINITION USED HERE DIFFERS BY A FACTOR OF 
! (1.0D0+eryd/ (cvlight*cvlight))
! from the other Juelich KKR codes
!#######################################################
if (config%nsra.eq.1) ek = sqrt(eryd)
if (config%nsra.eq.2) ek = sqrt(eryd+eryd*eryd/ (cvlight*cvlight))*(1.0D0+eryd/ (cvlight*cvlight))


!#######################################################
! if the wavefunctions are stored in memory, then there
! is no need to recalculate them
!#######################################################

if ( .not. allocated(wavefunction%rll)) then
  stop '[rhoval_new] rll not allocated'
end if





    imt1=cellnew%ipan_intervall(cellnew%npan_log+cellnew%npan_eq)+1




call rhooutnew(gauntcoeff,df,gmatll,ek,cellnew,wavefunction,rho2ns_complex(:,:,:), &
               config%nsra, &
               lmaxd,lmaxatom,lmmaxatom,wavefunction%lmsize,wavefunction%lmsize2,(2*lmaxd+1)**2,cellnew%nrmaxnew, &
               ispin,nspinden,imt1,cden,cdenlm,cdenns,shapefun,0, &
               gflle_part)  ! lda+u

! Copy to correct position in array gflle                                                ! lda+u
if (use_fullgmat==1) then                                                                ! lda+u
   lmslo = 1                                                                             ! lda+u
   lmshi = wavefunction%lmsize                                                           ! lda+u
else                                                                                     ! lda+u
   lmslo = (ispin - 1) * lmmaxatom + 1                                                   ! lda+u
   lmshi = lmslo - 1 + lmmaxatom                                                         ! lda+u
endif                                                                                    ! lda+u
density%gflle(lmslo:lmshi,lmslo:lmshi,ie) = &                                            ! lda+u
                     gflle_part(1:wavefunction%lmsize,1:wavefunction%lmsize)             ! lda+u
! Integration step for gfint.                                                            ! lda+u
density%gfint(lmslo:lmshi,lmslo:lmshi) = density%gfint(lmslo:lmshi,lmslo:lmshi) + &      ! lda+u
                     gflle_part(1:wavefunction%lmsize,1:wavefunction%lmsize) * df        ! lda+u


do jspin=nspinstart,nspinstop

    do lval=0,lmaxd
      call intcheb_cell(cden(:,lval,jspin),rho2ns_integrated(jspin), cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
      density%rho2ns_integrated(jspin)=density%rho2ns_integrated(jspin)+rho2ns_integrated(jspin)*df !*C0LL
      if (jspin<=2) then 
        density%den(lval,jspin,ie)=density%den(lval,jspin,ie)+rho2ns_integrated(jspin)
      end if
    end do

    if (jspin<=2) then 
      do lm1=1,(lmaxd+1)**2
        call intcheb_cell(cdenlm(:,lm1,jspin),density%denlm(lm1,jspin,ie), cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
      end do 
      call intcheb_cell(cdenns(:,jspin),density%den(lmaxd+1,jspin,ie), cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
      density%rho2ns_integrated(jspin)=density%rho2ns_integrated(jspin)+density%den(lmaxd+1,jspin,ie)*df !*C0LL
    end if


    call cheb2oldgrid(cell%nrmax, cellnew%nrmaxnew, (2*lmaxatom+1)**2, cell%rmesh, cellnew%ncheb, cellnew%npan_tot, cellnew%rpan_intervall, cellnew%ipan_intervall, rho2ns_complex(:,:,jspin), density%rho2ns_complex(:,:,jspin), cell%nrmax)

    if (config_testflag('write_rho2complex')) then
      write(4420,'(50000E25.14)') cellnew%rmeshnew
      write(4421,'(50000E25.14)') cell%rmesh
      do lm1=1,(2*lmaxatom+1)**2
        write(4424,'(50000E25.14)') rho2ns_complex(:,lm1,jspin)
        write(4425,'(50000E25.14)') density%rho2ns_complex(:,lm1,jspin)
      end do
    end if


end do

 do jspin=nspinstart,nspinstop
   if (jspin<=2) then
      do lval = 0,lmaxatom+1
        density%ncharge(lval,jspin) = density%ncharge(lval,jspin) + DIMAG(density%den(lval,jspin,ie)*df)
      end do
      do lval = 0,lmaxatom+1
        energyparts%espv(lval,jspin,iatom) = energyparts%espv(lval,jspin,iatom) + dimag( (eryd-efermi)*density%den(lval,jspin,ie)*df)
      end do
   end if
 end do


    if (config_testflag('write_rho2nscompnew')) then
     if (.not. allocated ( density%rho2ns_complexnew ) ) then 
       allocate (density%rho2ns_complexnew(cellnew%nrmaxnew,(2*lmaxatom+1)**2,nspinden))
       density%rho2ns_complexnew=(0.0D0,0.0D0)
     end if
     density%rho2ns_complexnew=density%rho2ns_complexnew+rho2ns_complex
    end if




if (.not. config_testflag('noscatteringmoment')) then
!                                    EK=0.0 here !!!!!!
  call rhooutnew(gauntcoeff,df,gmatll,(0.0D0,0.0D0),cellnew,wavefunction,rho2ns_complex(:,:,:), &
                 config%nsra, &
                 lmaxd,lmaxatom,lmmaxatom,wavefunction%lmsize,wavefunction%lmsize2,(2*lmaxd+1)**2,cellnew%nrmaxnew, &
                 ispin,nspinden,imt1,cden,cdenlm,cdenns,shapefun,0, &
                 gflle_part)  ! lda+u

  do jspin=nspinstart,nspinstop

    do lval=0,lmaxd
      call intcheb_cell(cden(:,lval,jspin),rho2ns_integrated(jspin), cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
      density%rho2ns_integrated_scattering(jspin)=density%rho2ns_integrated_scattering(jspin)+rho2ns_integrated(jspin)*df !*C0LL
    end do

    if (jspin<=2) then 
      call intcheb_cell(cdenns(:,jspin),temp1, cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
      density%rho2ns_integrated_scattering(jspin)=density%rho2ns_integrated_scattering(jspin)+temp1*df !*C0LL
    end if

  end do

end if

if (config%calcorbitalmoment==1) then
!                                    EK=0.0 here !!!!!!
  do ialpha=1,3
    call rhooutnew(gauntcoeff,df,gmatll,ek,cellnew,wavefunction,rho2ns_complex(:,:,:), &
                  config%nsra, &
                  lmaxd,lmaxatom,lmmaxatom,wavefunction%lmsize,wavefunction%lmsize2,(2*lmaxd+1)**2,cellnew%nrmaxnew, &
                  ispin,nspinden,imt1,cden,cdenlm,cdenns,shapefun,ialpha, &
                  gflle_part)  ! lda+u

    do jspin=nspinstart,nspinstop
      if (jspin<=2) then 

        do lval=0,lmaxd
          call intcheb_cell(cden(:,lval,jspin),rho2ns_integrated(jspin), cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
          density%orbitalmom(ialpha)=density%orbitalmom(ialpha) + rho2ns_integrated(jspin)*df!*C0LL
          density%orbitalmom_sp(jspin,ialpha)=density%orbitalmom_sp(jspin,ialpha) + rho2ns_integrated(jspin)*df!*C0LL
          density%orbitalmom_lm(lval,ialpha)=density%orbitalmom_lm(lval,ialpha) + rho2ns_integrated(jspin)*df!*C0LL
        end do
        if (jspin<=2) then 
          call intcheb_cell(cdenns(:,jspin),temp1, cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
          density%orbitalmom_ns(ialpha)=density%orbitalmom_ns(ialpha) + temp1*df!*C0LL
        end if
      end if
    end do
  end do !alpha

end if

end subroutine rhoval_new

end module mod_rhoval_new
