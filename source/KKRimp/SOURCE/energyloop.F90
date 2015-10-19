module mod_energyloop
  contains
subroutine energyloop(my_rank,mpi_size,ITSCF,cell, vpot, shapefun,zatom,natom,nspin,lmaxatom, &
                      lmaxd,density,ielast,ez,wez,config,gmat,gmatonsite,tmat,energyparts, &
                      ldau)                                                                         ! lda+u
      use nrtype
      use mod_rhoval
      use mod_rhoval_new
      use mod_rhovalfull
      use type_cell
      use type_cellnew
      use type_density
      use type_tmat
      use type_gmat
      use type_config
      use type_gmatonsite
      use type_shapefun
      use type_energyparts
      use type_wavefunction
      use type_ldau    ! lda+u
      use mod_read_spinorbit
      use mod_calctmat
      use mod_calctmatfull
      use mod_calctmat_bauernew
      use mod_gdyson
      use mod_interpolatecell
      use mod_rotatespin
      use mod_rotatespinframe
      use mod_checknan
!       use mod_config, only: config_testflag
      use mod_config
      use mod_timing
      use mod_complexdos3
      use mod_log, only: log_write
      use mod_read_angle
      use mod_mixbroydenspin
      use mod_wavefunctodisc
      use mod_gauntharmonics, only: gauntcoeff
      use mod_mpienergy, only: mpienergy_distribute
      use mod_calccouplingconstants, only: calcJijmatrix,calccouplingconstants_writeoutJij
      use mod_checkinterpolation
      implicit none
      integer                         :: ITSCF
      type(cell_type)                 ::  cell(natom)
      type(cell_typenew)                 ::  cellnew(natom)
      type(wavefunction_type),allocatable   ::  wavefunction(:,:)
      real(kind=dp)                   ::  vpot(:,:,:,:)
      type(shapefun_type)             ::  shapefun(natom)
      type(density_type),allocatable  ::  density(:)
      type(ldau_type)                 ::  ldau(:)                                  ! lda+u variables
      real(kind=dp)                   ::  zatom(natom)
      integer                         :: natom,nspin
      integer                         ::  lmaxatom(natom) !=3
      integer                         ::  lmaxd !=(lmaxatom+1)**2
      complex(kind=dpc)               ::  ez(ielast),wez(ielast)
      integer                         :: ielast
      type(config_type)             ::  config
      complex(kind=dpc)               ::  gmatll((lmaxd+1)**2,(lmaxd+1)**2)
!       real(kind=dp)                   ::  rho2ns(cell%nrmax,(2*lmaxatom+1)**2,nspin,natom), espv(0:lmaxatom+1,2)
!       double complex                  ::  den(0:lmaxatom+1) !,ez(iemxd), &
      integer                         ::  ie,ispin,iatom,jatom,ialpha
      integer                         ::  lm1,lm2
      type(tmat_type),allocatable    :: tmat(:,:) !(lmmaxd,lmmaxd)
      type(energyparts_type)          :: energyparts !(lmmaxd,lmmaxd)
!       integer                         :: itscf



      type(gmatonsite_type),allocatable   :: gmatonsite(:,:)
      integer                         ::  use_fullgmat
      type(gmat_type)               :: gmat
      integer                        :: ilval,ilm,ilm2,ilmatom,lmsmax !lda+u
      integer                        :: nlmhost
!       integer                        :: kgrefsoc

      integer                        :: natomtemp,nlmhosttemp
      integer                        :: ispinfac,ierror

      integer                        :: writeout_ie
      integer                        :: writeout_step
      integer                        :: nspinden
! write(*,*) 'test'
      double precision,allocatable    :: Jijmatrix(:,:,:,:)
      double precision,allocatable    :: Jijmatrix_temp(:,:,:,:)
      double precision,allocatable    :: Aimatrix(:,:)
      double precision,allocatable    :: Aimatrix_temp(:,:)

      double precision              ::  efermi
      integer                        :: saveGmat

!mpi
      integer,allocatable                      :: mpi_iebounds(:,:)
      integer                                  :: my_rank
      integer                                  :: mpi_size
      double complex,allocatable       ::  rho2ns_complex_temp(:,:,:)
      double precision,allocatable       ::  rho2ns_temp(:,:,:)
      double complex                   ::  rho2ns_integrated_temp(4)
      real(kind=dp),allocatable       ::  espv_temp(:,:,:)
      real(kind=dp)                   ::  ncharge_temp(0:lmaxd+1,nspin)
      double complex,allocatable      ::   den_temp(:,:,:) !,ez(iemxd), &
      double complex,allocatable      ::   denlm_temp(:,:,:) !,ez(iemxd), &
      double complex,allocatable      ::  gflle_temp(:,:,:),gfint_temp(:,:)              ! lda+u
      character(len=20)                        :: ctemp
     double complex      ::   orbitalmom_temp(3) !,ez(iemxd), &
     double complex      ::   orbitalmom_temp2(10,3) !,ez(iemxd), &
     double complex      ::   orbitalmom_temp3(2,3) !,ez(iemxd), &
     integer            :: idotime 
     double precision,allocatable   :: testpotdummy(:,:)
#ifdef MPI
       INCLUDE "mpif.h"
#endif

call timing_start('energyloop')

use_fullgmat=0

if(      config_runflag('force_fullgmat') &
    .or. config%kspinorbit==1 & 
    .or. config%ncoll==1 & 
    .or. config%nsra==3                     ) then

     use_fullgmat=1

end if

if (config%calcJijmat==1) then
  saveGmat=1
else 
  saveGmat=0
end if

if (use_fullgmat==1) then
  ispinfac=2
  nspinden=4
else
  ispinfac=1
  nspinden=nspin
end if

efermi=dreal( ez(ielast) )

! ###################################
! # prepare the iatom2nlmindex array
! # array returns for a given atom the first and last
! # index of the greens function
! ###################################
if (ITSCF==1) then
  ilm=1
  ilm2=1
  nlmhost=0

  call gdyson_read_kgrefsoc(gmat%kgrefsoc)

  allocate(gmat%iatom2nlmindex(3,natom))
  allocate(gmat%iatom2nlmindexhost(2,natom))
  do iatom=1,natom
    ilmatom=(lmaxatom(iatom)+1)**2 
    nlmhost=nlmhost + ilmatom * (gmat%kgrefsoc+1)
    gmat%iatom2nlmindex(1,iatom)=ilm
    gmat%iatom2nlmindex(2,iatom)=ilm+ilmatom-1
    gmat%iatom2nlmindex(3,iatom)=ilm+ispinfac*ilmatom-1
    ilm=ilm+ispinfac*ilmatom
    gmat%iatom2nlmindexhost(1,iatom)=ilm2
    gmat%iatom2nlmindexhost(2,iatom)=ilm2+ilmatom-1
    ilm2=ilm2+ilmatom
  end do !iatom=1,natom
  gmat%gmathostdim=nlmhost
  write(1337,*) '****************************************'
  write(1337,*) 'GHOST dim is', gmat%gmathostdim
  write(1337,*) '****************************************'
  gmat%gmatdim=gmat%iatom2nlmindex(3,natom)
  write(1337,*) 'GMAT matrix dimension'
  write(1337,*) 'GMAT dim is', gmat%gmatdim
  do iatom=1,natom
    write(1337,*) 'iatom=',iatom,'matrix index', gmat%iatom2nlmindex(1,iatom), gmat%iatom2nlmindex(3,iatom)
  end do
  write(1337,*) '****************************************'

end if !ITSCF==1

if (itscf<=3) then
  idotime=1
else
  idotime=0
end if

! nlmhost=0
! do iatom=1,natom
!   nlmhost=nlmhost+(lmaxatom(iatom)+1)**2 * (gmat%kgrefsoc+1)
! end do !iatom=1,natom


if (ITSCF==1) then


  open(unit=34536254,file='out_Jijmatrix')
  open(unit=34536256,file='out_JijDij')
  open(unit=34536258,file='out_Jijmatrix_local')
  open(unit=34536259,file='out_JijDij_local')
  open(unit=34536268,file='out_Aimatrix')

  open(unit=22349375,file='out_energytotal_eV')
  open(unit=22349376,file='out_energysp_eV')
  open(unit=22349378,file='out_energytotal_per_atom_eV')
  open(unit=22349379,file='out_energysp_per_atom_eV')

  if ( config_testflag('write_density') ) then
  open(unit=2342345,file='test_rho2ns')
  end if !( config_testflag('write_density') ) then
  if ( config_testflag('write_gmatonsite') ) then
    write(ctemp,'(I03.3)') my_rank
    open(unit=73467345,file='test_gmatonsite'//trim(ctemp)//'.txt')
  end if

  open(unit=2342348,file='test_rmesh')
  open(unit=2342349,file='test_rmeshnew')
  if (config%ncoll==1) then
    open(unit=23452324,file='out_magneticmoments_angle_nomix')
    open(unit=23452325,file='out_magneticmoments_nomix')
    open(unit=23452326,file='out_magneticmoments_angle')
    open(unit=23452327,file='out_magneticmoments')
  end if

  if (config%calcorbitalmoment==1) then
    open(unit=88943362,file='out_orbitalmoments')
    open(unit=88943363,file='out_orbitalmoments_ns')
    open(unit=88943364,file='out_orbitalmoments_lm')
    open(unit=88943365,file='out_orbitalmoments_sp')
  end if

end if

if (config%calcJijmat==1) then
  allocate(Jijmatrix(3,3,natom,natom))
  allocate(Aimatrix(3,natom))
  Jijmatrix=0.0D0
  Aimatrix=0.0D0
end if


! ###################################
! # gmatonsite allocation
! ###################################
if (ITSCF==1) then
  allocate( gmatonsite(natom,nspin), stat=ierror)
  do iatom=1,natom
    do ispin=1,nspin
!       write(*,*) iatom,ispin
      allocate( gmatonsite(iatom,ispin)%gmat( ispinfac*(1+lmaxatom(iatom))**2, ispinfac*(1+lmaxatom(iatom))**2 ) )
    end do !ispin
  end do !iatom=1,natom
end if

! ###################################
! # tmatll allocation
! ###################################
if (ITSCF==1) then
  allocate ( tmat(natom,nspin) )
  do iatom=1,natom
    do ispin=1,nspin
      allocate( tmat(iatom,ispin)%tmat( ispinfac*(1+lmaxatom(iatom))**2, ispinfac*(1+lmaxatom(iatom))**2 ) )
    end do !ispin
  end do !iatom=1,natom
end if
! ###################################
! # density allocation
! ###################################
if (ITSCF==1) then
  allocate(density(natom))
end if

if (.not. allocated(wavefunction)) then
  allocate(wavefunction(natom,nspin-use_fullgmat))
end if

if (config%ncoll==1) then
  if (ITSCF==1 .and. config_runflag('force_angles') .and. my_rank==0  ) then
    write(*,*) '##############################################################'
    write(*,*) '# Force angles in each iteration to satisfy angles listed in #'
    write(*,*) '# file kkrflex_angles                                        #'
    write(*,*) '##############################################################'
  end if

  if (ITSCF==1 .or. config_runflag('force_angles')  ) then
    call read_angle(natom,my_rank,density)
  end if

#ifdef MPI
      write(1337,*) 'Distributed new angles from my_rank 0 to others'
    do iatom=1,natom
      call mpi_bcast(density(iatom)%theta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
      call mpi_bcast(density(iatom)%phi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
      write(1337,*) 'Atom ',iatom
      write(1337,*) 'new theta',density(iatom)%theta/2.0D0/pi*360
      write(1337,*) 'new phi',density(iatom)%phi/2.0D0/pi*360
    end do
#endif

 call check_angle(density)

end if ! (config%ncoll==1)

!     density(1)%theta= 3.1415926D0/4.0D0 !/4.0D0
!     density(1)%phi=0.0D0! 3.1415926D0/4.0D0


do iatom=1,natom
  if (ITSCF==1) then
    allocate ( density(iatom)%rho2ns(cell(iatom)%nrmax,(2*lmaxatom(iatom)+1)**2,2), &
               density(iatom)%rho2ns_complex(cell(iatom)%nrmax,(2*lmaxatom(iatom)+1)**2,nspinden), &
!                density(iatom)%espv(0:lmaxatom(iatom)+1,2), &
               density(iatom)%den(0:lmaxatom(iatom)+2,nspin,ielast), &
               density(iatom)%denlm((lmaxatom(iatom)+1)**2,nspin,ielast), &
               density(iatom)%ncharge(0:lmaxatom(iatom)+1,nspin ), &
               density(iatom)%gfint(2*(lmaxatom(iatom)+1)**2,2*(lmaxatom(iatom)+1)**2),&        ! lda+u big matrix in lm*spin space
               density(iatom)%gflle(2*(lmaxatom(iatom)+1)**2,2*(lmaxatom(iatom)+1)**2,ielast) )  ! lda+u big matrix in lm*spin space
  end if !itscf
  density(iatom)%den=(0.0D0,0.0D0)
  density(iatom)%denlm=(0.0D0,0.0D0)
  density(iatom)%rho2ns=0.0D0
  density(iatom)%rho2ns_complex=(0.0D0,0.0D0)
  energyparts%ESPV=0.0D0
  density(iatom)%ncharge=0.0D0
  density(iatom)%rho2ns_integrated=(0.0D0,0.0D0)
  density(iatom)%rho2ns_integrated_scattering=(0.0D0,0.0D0)

  density(iatom)%orbitalmom=(0.0D0,0.0D0)
  density(iatom)%orbitalmom_lm=(0.0D0,0.0D0)
  density(iatom)%orbitalmom_ns=(0.0D0,0.0D0)
  density(iatom)%orbitalmom_sp=(0.0D0,0.0D0)

end do !iatom=1,natom

!----------------------------------------------------------------------------------------------------

! Initialize lda+u density matrix
do iatom = 1,natom                                                                  ! lda+u
   if (ldau(iatom)%lopt.ge.0) ldau(iatom)%denmatc(:,:,:) = (0.d0,0.d0)              ! lda+u
enddo                                                                               ! lda+u
!-----------------------------------------------------------------                  ! lda+u

! ###################################
! # open GHOST file / TMAT file
! ###################################
if (ITSCF==1) then
  open (3434560,access='direct',recl=wlength*4*gmat%gmathostdim*gmat%gmathostdim,file='kkrflex_greennew',form='unformatted')
  if ( config_testflag('write_tmat') ) then
    write(ctemp,'(I03.3)') my_rank
    open(43892349,file='out_tmat'//trim(ctemp)//'.txt')
  end if
end if


! ###################################
! # Distribute energy points IE to processors
! ###################################
call mpienergy_distribute(my_rank,mpi_size,ielast,mpi_iebounds)


!############################################### 
!# calculate the t-matrix out of the potential #
!###############################################
if (my_rank==0) then
   write(*,*) ' ###############################################'
   write(*,*) ' ###    Starting energy loop'
   write(*,*) ' ###############################################'
end if
writeout_ie=1
!####################################################################
!#################################################################### 
!#                      Start the energy loop
!####################################################################
!####################################################################
#ifdef MPI
call mpi_barrier(MPI_COMM_WORLD,ierror)
if (ITSCF==1) then 
  write(*,'(A,I3,A,I4,A,I4,A,I2,A)') '[energyloop] proc= ',my_rank,' Energy ie=', &
                    mpi_iebounds(1,my_rank),' to',mpi_iebounds(2,my_rank),' (',mpi_iebounds(2,my_rank)-mpi_iebounds(1,my_rank)+1,' energies)'
  call mpi_barrier(MPI_COMM_WORLD,ierror)
end if
#endif

if (ITSCF==1 .and. config%kspinorbit==1) then 
  call read_spinorbit(natom,cellnew,my_rank)
end if

!####################################################################
!####################################################################
if ( .not. config_testflag('nodisplayeloop') ) then
  if (my_rank==0) then
    write(*,'(A$)') '0%'
    do ie=1,ielast
      write(*,'(A$)') '*'
    end do
    write(*,'(A)') '100%'
    write(*,'(A$)') '  '
  end if
#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD,ierror)
#endif
end if
!####################################################################
!####################################################################


do ie=mpi_iebounds(1,my_rank),   mpi_iebounds(2,my_rank)

  if ( .not. config_testflag('nodisplayeloop') ) then
    write(*,'(A$)') '|'
  end if

  write(1337,*) '[energyloop] Energy ie=',ie

  if ( config_testflag('Jij(E)') .and. config%calcJijmat==1) then
    write(ctemp,'(I03.3)') ie
    open(unit=234932875,file='test_Jijmatrix_Eres_IE_int'//trim(ctemp)//'.dat')
    open(unit=234932876,file='test_Jijmatrix_Eres_IE'//trim(ctemp)//'.dat')
  end if 



! stop
! #ifndef MPI
!   writeout_step = 1+(mpi_iebounds(2,my_rank)-mpi_iebounds(1,my_rank))/4
!   if (mod(ie,writeout_step)==0 .or. ie==mpi_iebounds(2,my_rank)) then 
!     write(*,'(A,I4,A,I4,A)') '[energyloop] Energy ie=',writeout_ie,' to',ie,' done'
!     writeout_ie=ie
!   end if
! #endif


  ! ###################################
  ! # Start the spin loop
  ! ###################################
  do ispin=1,nspin-use_fullgmat
    ! ###################################
    ! # Start the atom loop
    ! ###################################
    call timing_start('vpot->tmat')
    do iatom=1,natom
      ! ###################################
      ! # calculate the T-matrix out of the potential
      ! ###################################

        if ( config_testflag('VPOT=0') ) then
          if (iatom==1 .and. ispin==1 .and.  my_rank==0 .and. ie==1 .and. itscf==1) then
            write(*,*) 'setting VPOT=0'
          end if
          VPOT(:,:,:,:)=0.0D0
        end if

        if ( config_testflag('Z=0') ) then
          if (iatom==1 .and. ispin==1 .and. my_rank==0 .and. ie==1 .and. itscf==1) then
            write(*,*) 'setting ZATOM=0'
          end if
          ZATOM(:)=0.0D0
        end if


        if ( config_testflag('VPOT=Vsph') ) then
          if (iatom==1 .and. ispin==1 .and. my_rank==0 .and. ie==1 .and. itscf==1) then
            write(*,*) 'setting VPOT=Vsph'
          end if
          VPOT(:,2:,:,:)=0.0D0
        end if


        if ( config_testflag('checkinterpol') ) then
         call  checkinterpolation(cell(iatom),cellnew(iatom),vpot(:,1,ispin,iatom),config)
        end if




      if (use_fullgmat==1 .and. .not. config_testflag('tmatnew')) then
        tmat(iatom,ispin)%tmat=(0.0D0,0.0D0)
        write(1337,*) 'call CALCTMATFULL'
        CALL  CALCTMATFULL(EZ(IE),VPOT(:,:,1,IATOM),VPOT(:,:,2,IATOM),CELL(IATOM),ZATOM(IATOM),lmaxatom(iatom), &
                            tmat(iatom,ispin)%tmat,config,nspin,lmaxd ) !  &
        write(1337,*) 'call CALCTMATFULL'

        if ( config_testflag('write_tmat') ) then
          write(43892349,'(2I4,5000e24.16)') ie,iatom,tmat(iatom,ispin)%tmat
        end if

      else
        tmat(iatom,ispin)%tmat=(0.0D0,0.0D0)
        write(1337,*) 'call CALCTMAT'

        if ( .not. config_testflag('tmatnew') ) then

           CALL  CALCTMAT(EZ(IE),VPOT(:,:,ISPIN,IATOM),CELL(IATOM),ZATOM(IATOM),lmaxatom(iatom), &
                tmat(iatom,ispin)%tmat,config,ispin,nspin ) !  &
           write(1337,*) 'call CALCTMAT'

           if ( config_testflag('write_tmat') ) then
              write(43892349,'(3I4,5000g24.16)') ie,ispin,iatom,tmat(iatom,ispin)%tmat
           end if

        else
           if (ie==mpi_iebounds(1,my_rank)) then


              call interpolatecell(cell(iatom)%nrcut(1)+1,shapefun(iatom)%thetas(:,:),(4*lmaxatom(iatom)+1)**2,CELL(IATOM), lmaxatom(iatom),cellnew(iatom),ispin,nspin,config,'shapefun',testpotdummy)





              call interpolatecell(1,VPOT(:,:,ISPIN,IATOM),(2*lmaxatom(iatom)+1)**2,CELL(IATOM),lmaxatom(iatom),cellnew(iatom),ispin,nspin,config,'potential',testpotdummy)
              if (use_fullgmat==1) then
                 call interpolatecell(1,VPOT(:,:,2,IATOM),(2*lmaxatom(iatom)+1)**2,CELL(IATOM),lmaxatom(iatom),cellnew(iatom),2,nspin,config,'potential',testpotdummy)
              end if
           end if

           if (ie==1) then
              write(2342348,'(50000E)') cell(iatom)%rmesh
              write(2342349,'(50000E)') cellnew(iatom)%rmeshnew
           end if


           if ( .not. config_testflag('calctmatfirstIter') .or.  itscf==1 ) then
! Pass  array wldaumat into subroutine, to be added to vpotll in the subroutine ! lda+u
              call calctmat_bauernew (cell(IATOM),tmat(iatom,ispin),lmaxatom(iatom),ez(IE),ZATOM(iatom), &
                   cellnew(iatom),wavefunction(iatom,ispin),ispin,nspin,config%kspinorbit, &
                   use_fullgmat,density(IATOM)%theta,density(IATOM)%phi,config%ncoll,config%nsra,config,idotime, &
                   ie,ldau(iatom))        ! lda+u


              if ( iatom > config%wavefunc_recalc_threshhold ) then
                 if ( config_testflag('rlltodisc') ) then
                    call wavefunctodisc_write(wavefunction(iatom,ispin),cellnew(iatom),iatom,ispin,my_rank)
                 else
                    deallocate(wavefunction(iatom,ispin)%rll,     wavefunction(iatom,ispin)%sll)
                    if (config%kspinorbit==1) then
                       deallocate(wavefunction(iatom,ispin)%rllleft, wavefunction(iatom,ispin)%sllleft)
                    end if
                 end if
                 wavefunction(iatom,ispin)%deallocate=1
              end if

           end if

           if (config%ncoll==1) then
              call rotatematrix(tmat(iatom,ispin)%tmat,density(iatom)%theta, density(iatom)%phi,(lmaxatom(iatom)+1)**2,'loc->glob')
           end if

        if ( config_testflag('write_tmat') ) then
          write(43892349,'(3I4,5000g24.16)') ie,ispin,iatom,tmat(iatom,ispin)%tmat
        end if


        if ( config_testflag('calctmatstop') ) then
          stop 'stop after eloop'
        end if


         end if


      end if
    end do !iatom
    call timing_stop('vpot->tmat')

    if ( config_testflag('stopcalctmat') ) then
      stop '[energyloop] stop after calctmat because of test flag'
    end if

    !##############################
    ! calculate the Greens functions
    !##############################
    call timing_start('gref->gmat')
    write(1337,*) 'call gdyson'
    call gdyson(3434560,ie,ispin,nspin,natom,lmaxatom,tmat,use_fullgmat,gmat,& 
                gmatonsite,ielast,mpi_iebounds(:,my_rank),ITSCF,saveGmat)
    write(1337,*) 'end call gdyson'
    if(config_testflag('gtest')) write(30000+my_rank,'(832E)') gmat%gmat

    call timing_stop('gref->gmat')


       if (config%calcJijmat==1) then
         call timing_start('gmat->Jij')
         call calcJijmatrix(gmat,tmat,natom,wez(ie),Jijmatrix,Aimatrix)
         call timing_stop('gmat->Jij')
       end if


        if ( config%ncoll==1 .and. config_testflag('calctmatfirstIter') ) then
          if (config%ncoll==1) then
            do iatom=1,natom
              call rotatematrix(tmat(iatom,ispin)%tmat,density(iatom)%theta, density(iatom)%phi,(lmaxatom(iatom)+1)**2,'glob->loc')
            end do !iatom
          end if
        end if

        if (config%ncoll==1) then
          do iatom=1,natom
            call rotatematrix(gmatonsite(iatom,ispin)%gmat,density(iatom)%theta, density(iatom)%phi,(lmaxatom(iatom)+1)**2,'glob->loc')
          end do !iatom
        end if


    if ( config_testflag('write_gmatonsite') ) then
      do iatom=1,natom
        write(73467345,*) '# iatom' ,iatom,'ispin',ispin,'dimension GF' ,ubound(gmatonsite(iatom,ispin)%gmat)
        write(73467345,'(3I4,5000g24.16)') iatom,ispin,ie,gmatonsite(iatom,ispin)%gmat
      end do !iatom
    end if
    !##############################333
    ! calculate the density out of the onsite greens function
    !##############################333
    call timing_start('gonsite->density')
    do iatom=1,natom

       if ( config_testflag('Gref=0') ) then
          if (ITSCF==1.and. my_rank==0) print *, 'GMAT=0'
          gmatonsite(iatom,ispin)%gmat=(0.0D0,0.0D0)
       end if


       if ( .not. config_testflag('tmatnew') ) then

          if (use_fullgmat==0) then

             write(1337,*) 'call RHOVAL'
             call RHOVAL(ie,ielast,ez(ie) ,WEZ(IE), gmatonsite(iatom,ispin)%gmat,ISPIN, NSPIN, &
                  IATOM,CELL(iatom),VPOT(:,:,ispin,iatom),SHAPEFUN(iatom), GAUNTCOEFF(lmaxatom(iatom)), ZATOM(iatom),&
                  density(iatom), LMAXATOM(iatom), (LMAXATOM(iatom)+1)**2,config ,lmaxd,energyparts,efermi)
             write(1337,*) 'end RHOVAL'

          else

             if (ispin==2) stop' [eloop] use_fullmat==1 but nspin==2'
             write(1337,*) 'call RHOVAL'
             call RHOVALFULL(ie,ielast,ez(ie) ,WEZ(IE), gmatonsite(iatom,ispin)%gmat, NSPIN, &
                  IATOM,CELL(iatom),VPOT(:,:,1,iatom),VPOT(:,:,2,iatom),SHAPEFUN(iatom), GAUNTCOEFF(lmaxatom(iatom)), ZATOM(iatom),&
                  density(iatom), LMAXATOM(iatom), (LMAXATOM(iatom)+1)**2,config ,lmaxd,energyparts,efermi)
             write(1337,*) 'end RHOVAL'

          end if

          if ( config_testflag('tmatdebug') ) then
             do lm1=1,(2*LMAXATOM(iatom)+1)**2
                write(9004,'(50000F)') density(iatom)%rho2ns(:,lm1,ispin)
             end do
          end if

       else ! ( .not. config_testflag('tmatnew') )


          if ( wavefunction(iatom,ispin)%deallocate==1 ) then 
             if ( config_testflag('rlltodisc') ) then
                call wavefunctodisc_read(wavefunction(iatom,ispin),cellnew(iatom),iatom,ispin)
             else
                write(1337,*) 'Single site wavefunctions not allocated for atom',iatom
                write(1337,*) 'Recalculating wave functions ...'
                call calctmat_bauernew (cell(IATOM),tmat(iatom,ispin),lmaxatom(iatom),ez(IE),ZATOM(iatom), &
                     cellnew(iatom),wavefunction(iatom,ispin),ispin,nspin,config%kspinorbit, &
                     use_fullgmat,density(IATOM)%theta,density(IATOM)%phi,config%ncoll,config%nsra,config,idotime, &
                     ie,ldau(iatom))        ! lda+u

             end if
          end if

          write(*,*) 'entering rhovalnew iatom',iatom
          call rhoval_new(ez(ie),ie,wez(ie),cellnew(iatom),wavefunction(iatom,ispin), &                  
               cell(iatom),gmatonsite(iatom,ispin)%gmat,iatom,ispin,nspin,SHAPEFUN(iatom), &
               GAUNTCOEFF(lmaxatom(iatom)), ZATOM(iatom), DENSITY(iatom), &
               LMAXATOM(iatom),(LMAXATOM(iatom)+1)**2,config,lmaxd,energyparts,config%kspinorbit,use_fullgmat,nspinden,efermi, &
               ldau(iatom))        ! lda+u
          write(*,*) 'exited rhovalnew iatom',iatom


          if ( wavefunction(iatom,ispin)%deallocate==1 ) then 
             deallocate(wavefunction(iatom,ispin)%rll,     wavefunction(iatom,ispin)%sll)
             if (config%kspinorbit==1) then
                deallocate(wavefunction(iatom,ispin)%rllleft, wavefunction(iatom,ispin)%sllleft)
             end if
          end if



          if ( config_testflag('tmatdebug') ) then
             write(9002,'(50000F)') cell(iatom)%rmesh(:)
             do lm1=1,(2*LMAXATOM(iatom)+1)**2
                write(9004+my_rank,'(50000E)') density(iatom)%rho2ns_complex(:,lm1,1)
                write(9004+my_rank,'(50000E)') density(iatom)%rho2ns_complex(:,lm1,2)
             end do
          end if




       end if

       if ( config_testflag('rhotest') ) then
          write(9902,'(50000F)') cell(iatom)%rmesh(:)
          do lm1=1,(2*LMAXATOM(iatom)+1)**2
             write(9904+my_rank,'(50000E)') density(iatom)%rho2ns(:,lm1,ispin)
          end do
       end if

       if (config_testflag('write_rho2nscompnew')) then
          do lm1=1,(2*LMAXATOM(iatom)+1)**2
             write(9604+my_rank,'(50000E)') density(iatom)%rho2ns_complexnew(:,lm1,1)
             write(9704+my_rank,'(50000E)') density(iatom)%rho2ns_complexnew(:,lm1,2)
          end do
       end if


       if ( config_testflag('checknan') ) then
          call checknan(DENSITY(iatom)%rho2ns_complex,ierror)
          if (ierror==1) then 
             write(*,*) '[energyloop] error scf = ',itscf,'ie ',ie,'ispin',ispin
             stop
          end if
       end if


    end do !iatom

    ! ###################################
    ! # Start the spin loop
    ! ###################################
    call timing_stop('gonsite->density')
 end do !ispin

! stop 'after spinloop'
!####################################################################
!#################################################################### 
!#                      Start the energy loop
!####################################################################
!####################################################################

if ( config_testflag('eloopstop') ) then
  stop 'stop after eloop'
end if

end do !ie


#ifdef MPI

call mpi_barrier(MPI_COMM_WORLD,ierror)

call log_write('>>>>>>>>>>>>>>>>>>> MPI comm >>>>>>>>>>>>>>>>>>>')
call timing_start('MPI comm : charge density')
do iatom=1,natom
       if (.not. config_testflag('tmatnew') ) then

        allocate(rho2ns_temp(cell(iatom)%nrmax,(2*lmaxatom(iatom)+1)**2,2))

!         write(*,*) ubound(density(iatom)%rho2ns)
!         write(*,*) ubound(rho2ns_temp)
!         write(*,*) cell(iatom)%nrmax*(2*lmaxatom(iatom)+1)**2*2

  
        CALL MPI_REDUCE(density(iatom)%rho2ns, rho2ns_temp, cell(iatom)%nrmax*(2*lmaxatom(iatom)+1)**2*2, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
        density(iatom)%rho2ns=rho2ns_temp
  
        deallocate(rho2ns_temp)!(cell(iatom)%nrmax,(2*lmaxatom(iatom)+1)**2,2))

       else if ( config_testflag('tmatnew') ) then

        allocate(rho2ns_complex_temp(cell(iatom)%nrmax,(2*lmaxatom(iatom)+1)**2,nspinden))
!        write(*,*) '1','iatom= ',iatom
        CALL MPI_REDUCE(density(iatom)%rho2ns_complex, rho2ns_complex_temp, cell(iatom)%nrmax*(2*lmaxatom(iatom)+1)**2*nspinden, &
                        MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
        density(iatom)%rho2ns_complex=rho2ns_complex_temp
!        write(*,*) '2','iatom= ',iatom
  
  !        if (my_rank==0) write(*,*) density(iatom)%espv 
        deallocate(rho2ns_complex_temp)!(cell(iatom)%nrmax,(2*lmaxatom(iatom)+1)**2,2))
!        write(*,*) '3','iatom= ',iatom

!-----------------------------------------------------------------------------------     ! lda+u
!  mpi-reduce green-function parts                                                       ! lda+u
       lmsmax = 2*(lmaxatom(iatom)+1)**2  ! (2x2 in spin space)                          ! lda+u
       allocate(gfint_temp(lmsmax,lmsmax))                                               ! lda+u
       CALL MPI_REDUCE(density(iatom)%gfint, gfint_temp, lmsmax**2, &                    ! lda+u
                        MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)           ! lda+u
       density(iatom)%gfint=gfint_temp                                                   ! lda+u
       deallocate(gfint_temp)                                                            ! lda+u
!         write(*,*) '4','iatom= ',iatom

       allocate(gflle_temp(lmsmax,lmsmax,ielast))                                        ! lda+u
       CALL MPI_REDUCE(density(iatom)%gflle, gflle_temp, ielast*lmsmax**2, &             ! lda+u
                        MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)           ! lda+u
       density(iatom)%gflle=gflle_temp                                                   ! lda+u
       deallocate(gflle_temp)                                                            ! lda+u
!-----------------------------------------------------------------------------------     ! lda+u
!         write(*,*) '5','iatom= ',iatom

       end if



       CALL MPI_REDUCE(density(iatom)%rho2ns_integrated, rho2ns_integrated_temp, 4, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%rho2ns_integrated=rho2ns_integrated_temp


       CALL MPI_REDUCE(density(iatom)%rho2ns_integrated_scattering, rho2ns_integrated_temp, 4, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%rho2ns_integrated_scattering=rho2ns_integrated_temp

      CALL MPI_REDUCE(density(iatom)%orbitalmom, orbitalmom_temp, 3, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%orbitalmom=orbitalmom_temp

      CALL MPI_REDUCE(density(iatom)%orbitalmom_ns, orbitalmom_temp, 3, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%orbitalmom_ns=orbitalmom_temp

      CALL MPI_REDUCE(density(iatom)%orbitalmom_lm, orbitalmom_temp2, 30, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%orbitalmom_lm=orbitalmom_temp2

      CALL MPI_REDUCE(density(iatom)%orbitalmom_sp, orbitalmom_temp3, 6, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%orbitalmom_sp=orbitalmom_temp3





! 
! 
! 
! 
       do ispin=1,nspin
         ncharge_temp=0.0D0
         CALL MPI_REDUCE(density(iatom)%ncharge, ncharge_temp,nspin*(lmaxatom(iatom)+2), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
         density(iatom)%ncharge=ncharge_temp(0:lmaxatom(iatom)+1,:)
       end do !ispin=1,nspin
! 
! write(*,*) 'include den and denlm'
! !         write(my_rank+1000,'(5000F)') dimag(density(iatom)%den)
       allocate(den_temp(0:lmaxatom(iatom)+2,nspin,ielast))
       CALL MPI_REDUCE(density(iatom)%den, den_temp, (1+lmaxatom(iatom)+2)*nspin*ielast, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%den=den_temp
       deallocate(den_temp)
! 
       allocate(denlm_temp((lmaxatom(iatom)+1)**2,nspin,ielast))
       CALL MPI_REDUCE(density(iatom)%denlm, denlm_temp, (lmaxatom(iatom)+1)**2*nspin*ielast, &
                       MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
       density(iatom)%denlm=denlm_temp
       deallocate(denlm_temp)
! !         write(my_rank+2000,'(5000F)') dimag(density(iatom)%den)

end do

       if (config%calcJijmat==1) then
         allocate(Jijmatrix_temp(3,3,natom,natom))
         CALL MPI_REDUCE(Jijmatrix, Jijmatrix_temp, 9*natom**2, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
         Jijmatrix=Jijmatrix_temp
         deallocate(Jijmatrix_temp)
         allocate(Aimatrix_temp(3,natom))
         CALL MPI_REDUCE(Aimatrix, Aimatrix_temp, 3*natom, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierror)
         Aimatrix=Aimatrix_temp
         deallocate(Aimatrix_temp)


       end if
  
!         write(*,*) '7','iatom= ',iatom

allocate(espv_temp(0:lmaxd+1,NSPIN,NATOM))
CALL MPI_REDUCE(energyparts%espv, espv_temp, (lmaxd+2)*NSPIN*NATOM, &
                MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD, ierror)

energyparts%espv=espv_temp
deallocate(espv_temp)
call timing_stop('MPI comm : charge density')
call log_write('<<<<<<<<<<<<<<<<<<< MPI comm <<<<<<<<<<<<<<<<<<<')
#endif

if (my_rank==0) then 


if (config%calcJijmat==1) then

call log_write('>>>>>>>>>>>>>>>>>>> calcJij >>>>>>>>>>>>>>>>>>>')
call calccouplingconstants_writeoutJij(natom,Jijmatrix,Aimatrix,density,ITSCF)
call log_write('<<<<<<<<<<<<<<<<<<< calcJij <<<<<<<<<<<<<<<<<<<')





end if

!      do lm1=1,(2*lmaxatom(1)+1)**2
!        Do ispin=1,ubound(density(1)%RHO2NS_COMPLEX,3)
!          write(7990+ispin,'(5000E)') density(1)%RHO2NS_COMPLEX(:,lm1,ispin)
!        END DO
!      end do


if (config%ncoll==1) then
  do iatom=1,natom
    call rotatespin(density(iatom),cell(iatom),lmaxatom(iatom),config,iatom)
  end do !iatom

! density%theta
  if (config%imixspin>0 .and. ITSCF>1) then
    call mixbroydenspin (natom,density,20,ITSCF)
  else if (config%imixspin<0 .and. ITSCF>1) then
    call mixbroydenspinangle (natom,density,20,ITSCF)
  end if

  do iatom=1,natom
    call rotatevector(density(iatom)%rho2ns_complex,density(iatom)%rho2ns, &
                      cell(iatom)%nrmax,(2*lmaxatom(iatom)+1)**2, &
                      density(iatom)%theta,density(iatom)%phi,density(iatom)%thetaold,density(iatom)%phiold)

    write(23452326,'(5000F)') density(iatom)%theta*180/pi,density(iatom)%phi*180/pi
    write(23452327,'(5000F)') density(iatom)%magmoment
  end do !iatom



else

  if ( config_testflag('tmatnew') ) then
    do iatom=1,natom
      do ispin=1,nspin 
        do lm1=1,(2*lmaxatom(iatom)+1)**2
          density(iatom)%rho2ns(:,lm1,ispin)= dimag ( density(iatom)%rho2ns_complex(:,lm1,ispin) )
        end do
      end do 
    end do
  end if
end if






! write(*,*) 'test'
!c
!c---> gather the parts in the last energy loop
!c     be carefull!!! this is for the temperature!!!

! do ispin=1,nspin
!   do iatom=1,natom
!     do ilval=0,lmaxatom(iatom)+1
!             energyparts%espv(ilval,ispin,iatom) = energyparts%espv(ilval,ispin,iatom) - &
!                                dreal(ez(ielast))*dimag(density(iatom)%ncharge(ilval,ispin))
! 
!     end do
!   end do
! end do 





! c           write(6,'(13x,2e12.6)') e
! c           write(6,'(13x,i3,1x,e12.6)') l,espv(l,i1,ispin)



!   if(config_runflag('ldos')) then
!     do iatom=1,natom
!       do ispin=1,nspin
!         OPEN(UNIT=30, &
!             FILE="out_ldos.atom="//char(48+IATOM/10)//char(48+mod(IATOM,10))//"_spin"//char(48+ISPIN)//".dat")
!         do ie=1,ielast
!           WRITE(30,'(300G24.16)') DREAL(EZ(IE)),(-DIMAG(DENSITY(IATOM)%DEN(ilval,ISPIN,IE))/PI,ilval=0,LMAXATOM(IATOM)+1)
!         end do !ie
!         close(30)
!       end do !ispin
!     end do !iatom
!   end if !lmdos
! 
! 
! 
!   if(config_runflag('lmdos')) then
!     do iatom=1,natom
!       do ispin=1,nspin
!         OPEN(UNIT=30, &
!             FILE="out_lmdos.atom="//char(48+IATOM/10)//char(48+mod(IATOM,10))//"_spin"//char(48+ISPIN)//".dat")
!         do ie=1,ielast
!           WRITE(30,'(300G24.16)') DREAL(EZ(IE)),(-DIMAG(DENSITY(IATOM)%DENLM(LM1,ISPIN,IE))/PI,LM1=1,(LMAXATOM(IATOM)+1)**2)
!         end do !ie
!         close(30)
!       end do !ispin
!     end do !iatom
!   end if !lmdos




  if(config_runflag('ldos').or.config_runflag('lmdos')) then

! calculate l-summed dos and store in den(lmax+2,ispin,ie)
    do iatom=1,natom
      do ispin=1,nspin
        do ie=1,ielast
          DENSITY(IATOM)%DEN(lmaxatom(iatom)+2,ISPIN,IE) = (0.d0,0.d0)
                  do ilval=0,lmaxatom(iatom)+1
            DENSITY(IATOM)%DEN(lmaxatom(iatom)+2,ISPIN,IE) =                      &
            DENSITY(IATOM)%DEN(lmaxatom(iatom)+2,ISPIN,IE) + DENSITY(IATOM)%DEN(ilval,ISPIN,IE)
                  enddo ! ilval
        end do !ie
      end do !ispin
    end do !iatom

! write out dos
    do iatom=1,natom
      do ispin=1,nspin

      call complexdos3(lmaxatom(iatom),ielast,iatom,ispin,nspin,density(iatom)%den,density(iatom)%denlm,ez)

! l-dos
        OPEN(UNIT=30, &
            FILE="out_ldos.atom="//char(48+IATOM/10)//char(48+mod(IATOM,10))//"_spin"//char(48+ISPIN)//".dat")
        do ie=1,ielast
          WRITE(30,'(300G24.16)') DREAL(EZ(IE)),-DIMAG(DENSITY(IATOM)%DEN(LMAXATOM(IATOM)+2,ISPIN,IE))/PI,    & ! e, l-summed dos,
                                (-DIMAG(DENSITY(IATOM)%DEN(ilval,ISPIN,IE))/PI,ilval=0,LMAXATOM(IATOM)+1)       ! l-resolved dos
        end do !ie
        close(30)
! lm-dos
        OPEN(UNIT=30, &
            FILE="out_lmdos.atom="//char(48+IATOM/10)//char(48+mod(IATOM,10))//"_spin"//char(48+ISPIN)//".dat")
        do ie=1,ielast
          WRITE(30,'(300G24.16)') DREAL(EZ(IE)),(-DIMAG(DENSITY(IATOM)%DENLM(LM1,ISPIN,IE))/PI,LM1=1,(LMAXATOM(IATOM)+1)**2)
        end do !ie
        close(30)
! complex lmdos
        OPEN(UNIT=30, &
            FILE="out_clmdos.atom="//char(48+IATOM/10)//char(48+mod(IATOM,10))//"_spin"//char(48+ISPIN)//".dat")
                WRITE(30,*) ielast, 'IELAST'
                WRITE(30,*) lmaxatom(iatom), 'LMAX'
        do ie=1,ielast
          WRITE(30,'(300G24.16)') EZ(IE), -DENSITY(IATOM)%DEN(LMAXATOM(IATOM)+2,ISPIN,IE)/PI,     &
                                  (-(DENSITY(IATOM)%DENLM(LM1,ISPIN,IE))/PI,LM1=1,(LMAXATOM(IATOM)+1)**2)
        end do !ie
                write(30,*) '    '
        close(30)

      end do !ispin
    end do !iatom
  end if ! ldos or lmdos

  if(config%calcorbitalmoment==1) then
  do iatom=1,natom
    write(*,*       ) 'Orbital moments are: '
    write(*, '(100F9.5)') density(iatom)%orbitalmom
    write(88943362,'(10F)') density(iatom)%orbitalmom
    write(88943363,'(10F)') density(iatom)%orbitalmom_ns
  
    write(88943364,*) '#atom',iatom
    write(88943365,*) '#atom',iatom
    do ialpha=1,3
      write(88943364,'(I,100F9.5)') ialpha,density(iatom)%orbitalmom_lm(0:lmaxatom(iatom),ialpha)
      write(88943365,'(I,100F9.5)') ialpha,density(iatom)%orbitalmom_sp(1:2,ialpha)
    end do
  end do !iatom
  end if



!-------------------------------------------------------------------------
!-- charge density for up and down to spin density
!-------------------------------------------------------------------------
! before     : rho2ns(1/2) = spin up/spin down
! afterwards : rho2ns(1) = charge density
!              rho2ns(2) = spin density
  if (nspin==2) then
    do iatom=1,natom
      density(iatom)%rho2ns(:,:,2) = density(iatom)%rho2ns(:,:,2) & 
                                          - density(iatom)%rho2ns(:,:,1)
      density(iatom)%rho2ns(:,:,1) = density(iatom)%rho2ns(:,:,2) & 
                                          + 2.0D0 * density(iatom)%rho2ns(:,:,1)
    end do !iatom
  end if


  if ( config_testflag('write_density') ) then
    if (nspin==1) then
      do iatom=1,natom
        write(2342345,*) '# atom ',iatom
        write(2342345,'(50000g24.16)') density(iatom)%rho2ns(:,:,1)
      end do !iatom
    elseif (nspin==2) then
      do iatom=1,natom
        write(2342345,*) '# atom ',iatom,' charge density'
        write(2342345,'(50000g24.16)') density(iatom)%rho2ns(:,:,1)
        write(2342345,*) '# atom ',iatom,' spin density'
        write(2342345,'(50000g24.16)') density(iatom)%rho2ns(:,:,2)
      end do !iatom
    else
      stop '[energyloop] nspin error'
    end if
  end if ! config_testflag('write_gmatonsite')

end if !my_rank=0

call timing_stop('energyloop')

end subroutine energyloop

end module mod_energyloop
