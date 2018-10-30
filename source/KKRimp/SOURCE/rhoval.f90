MODULE MOD_RHOVAL
  CONTAINS
!-------------------------------------------------------------------------
!> Summary: Main driver for valence charge density, old solver
!> Category: physical-observables, KKRimp
!>
!> @todo Code looks very messy. Clean up required. @endtodo
!> @warning Contains ToDos for LDA+U case. @endwarning
!-------------------------------------------------------------------------

      SUBROUTINE RHOVAL(IE,IELAST, ERYD ,WEZ, GMATLL,ISPIN, NSPIN, &
                        IATOM,CELL,VPOT,SHAPEFUN,GAUNTCOEFF, ZATOM, DENSITY, &
                         LMAXATOM,LMMAXATOM,config,lmaxd,energyparts,efermi)
 !  DEN, &
 !                         RHO2NS, ESPV, &
! C
! C **********************************************************************
! C * For KREL = 1 (relativistic mode)                                   *
! C *                                                                    *
! C *  NPOTD = 2 * NATYPD                                                *
! C *  LMMAXD = 2 * (LMAXD+1)^2                                          *
! C *  NSPIND = 1                                                        *
! C *                                                                    *
! C *  LDA+U implementation     Mar. 2002-Dec.2004                       *
! C *                           ph.mavropoulos, h. ebert, v. popescu     *
! C * Notes:                                                             *
! C *  average WLDAU for spherical wavefunctions:                        *
! C *  The spherical part of the d or f wavefunction is found by adding  *
! C *  the average interaction potential WLDAUAV to the spherical        *
! C *  potential. Then the non-spherical parts are found by using only   *
! C *  the deviation of WLDAU from the average. This speeds up the       *
! C *  convergence of the Born series. See also subroutines              *
! C *  regsol, pnstmat and pnsqns                                        *
! C *                                                                    *
! C **********************************************************************
      use nrtype
      use mod_physic_params, only: cvlight
!       use configparams, only: icst,ins,nsra
      use type_cell
      use type_gauntcoeff
      use type_shapefun
      use type_density
      use type_config
      use type_energyparts
!       use mod_gauntharmonics, only: gauntcoeff
      use mod_wfmesh
      use mod_cradwf
      use mod_pnsqns
      use mod_rholm
      use mod_rhons
      use mod_timing
      use mod_config, only: config_runflag,config_testflag
      implicit none
!interface variables
      integer,intent(in)                    ::   lmaxd
      integer,intent(in)                    ::   ie
      integer                               ::   ielast
      complex(kind=dpc),intent(in)          ::   eryd
      complex(kind=dpc),intent(in)          ::   wez
      complex(kind=dpc),intent(in)          ::   gmatll(lmmaxatom,lmmaxatom)
      integer,intent(in)                    ::   ispin
      integer,intent(in)                    ::   nspin
      integer,intent(in)                    ::   iatom
      type(cell_type),intent(in)            ::   cell
      real(kind=dp),intent(in)              ::   vpot(1:cell%nrmaxd,(2*lmaxd+1)**2)
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

!local variables
  ! LDA+U Stuff
      real(kind=dp),allocatable             ::   vpottemp(:,:)
      integer                               ::   idoldau, lopt,lmlo ,lmhi
      real(kind=dp)                         ::   wldauav
      real(kind=dp),allocatable             ::   cutoff(:) ,wldau(:,:,:)

      complex(kind=dpc),parameter           ::   czero=(0.0D0,0.0D0)
      integer                               ::   lm1,lm2,lval
      complex(kind=dpc)                     ::   df,ek
      logical,parameter                     ::   lirrsol=.true.

      complex(kind=dpc)                     ::   alpha(0:lmaxatom),ar(lmmaxatom,lmmaxatom), &
                                                 cr(lmmaxatom,lmmaxatom), &
                                                 dr(lmmaxatom,lmmaxatom), &
                                                 ekl(0:lmaxatom),fz(cell%nrmax,0:lmaxatom), &
                                                 pns(lmmaxatom,lmmaxatom,cell%nrmin_ns:cell%nrmax,2), & 
                                                 pz(cell%nrmax,0:lmaxatom), &
                                                 qns(lmmaxatom,lmmaxatom,cell%nrmin_ns:cell%nrmax,2), &
                                                 qz(cell%nrmax,0:lmaxatom), &
                                                 sz(cell%nrmax,0:lmaxatom),tmat(0:lmaxatom)
      real(kind=dp)                         ::   rs(cell%nrmax,0:lmaxatom), s(0:lmaxatom)

      external daxpy,dscal,drvrho 

      intrinsic atan,dble,dimag,sqrt

!---------------------------------------------------------------
!---   LDA+U to implement in the future 
!---------------------------------------------------------------
!       IF ( IDOLDAU.EQ.1 ) THEN
!          WLDAUAV = 0.D0
!          LMLO = LOPT*LOPT + 1
!          LMHI = (LOPT+1)*(LOPT+1)
!          MMAX = LMHI - LMLO + 1
!          DO M1 = 1,MMAX                                        
!             WLDAUAV = WLDAUAV + WLDAU(M1,M1,ISPIN)         
!          ENDDO                                                 
!          WLDAUAV = WLDAUAV/DBLE(MMAX)                        
!
! -> Note: Application if WLDAU makes the potential discontinuous.
!    A cutoff can be used if necessary to make the potential continuous
!    for example (array bounds should be adjusted):
!
!ccc            CUTOFF(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) *
!ccc     &                   ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
!ccc            CUTOFF(IR) = 1D0/CUTOFF(IR)
!---------------------------------------------------------------
!          DO M1 = 1,IRMD
!             CUTOFF(M1) = 1.D0
!          END DO
!       END IF
!---------------------------------------------------------------

DF = WEZ/DBLE(NSPIN)

IDOLDAU=0
if (idoldau==0) then 
! allocate(wldau(1,1,2))
  allocate(cutoff(cell%nrmax), wldau(2*lmaxatom+1,2*lmaxatom+1,nspin))
else
 stop 'ldau not implemented'
end if
!   call timing_start()

allocate(vpottemp(cell%nrmin_ns:cell%nrmax,(2*lmaxatom+1)**2))
!-----------------------------------------------------------------------
! calculate the spherical wavefunctions
!-----------------------------------------------------------------------
call wfmesh(eryd,ek,cvlight,config%nsra,zatom,cell%rmesh,s,rs, &
            cell%nrmax,lmaxatom)
call cradwf(eryd,ek,config%nsra,alpha,cell%npan,cell%nrcut,cvlight,rs,s, &
            pz,fz,qz,sz,tmat,vpot(:,1),cell%drmeshdi,cell%rmesh,zatom,lirrsol, &
            idoldau,lopt,wldauav,cutoff,lmaxatom,lmaxatom+1,cell%nrmax)


! if (maxval(vpot(:,1))>10D10) stop 'error'
! write(*,*) maxval(vpot(:,1))
! stop

!-----------------------------------------------------------------------
! transform to non-spherical wavefunctions
!-----------------------------------------------------------------------


IF (config%INS.GT.0) vpottemp=vpot(cell%nrmin_ns:cell%nrmax,1:(2*lmaxatom+1)**2)

! write(88888,*) AR,CR,DR,cell%DRMESHDI,EK,config%ICST,PZ,QZ,FZ,SZ, &
!                           PNS,QNS,config%NSRA,vpottemp,cell%NPAN,cell%NRCUT, &
!                           gauntcoeff%CLEB,gauntcoeff%ICLEB,gauntcoeff%IEND,gauntcoeff%LOFLM,LMAXATOM,&
! !                           IDOLDAU,LOPT,LMLO,LMHI,&
!                           WLDAU(1,1,ISPIN),WLDAUAV,CUTOFF,&
!                           cell%nrmax, cell%nrmin_ns,gauntcoeff%ncleb, &
!                           lmaxatom,2*lmaxatom+1,lmmaxatom,(2*lmaxatom+1)**2
! stop


IF (config%INS.GT.0) CALL PNSQNS(AR,CR,DR,cell%DRMESHDI,EK,config%ICST,PZ,QZ,FZ,SZ, &
                          PNS,QNS,config%NSRA,vpottemp,cell%NPAN,cell%NRCUT, &
                          gauntcoeff%CLEB,gauntcoeff%ICLEB,gauntcoeff%IEND,gauntcoeff%LOFLM,LMAXATOM,&
                          IDOLDAU,LOPT,LMLO,LMHI,&
                          WLDAU(:,:,ISPIN),WLDAUAV,CUTOFF,&
                          cell%nrmax, cell%nrmin_ns,gauntcoeff%ncleb, &
                          lmaxatom,2*lmaxatom+1,lmmaxatom,(2*lmaxatom+1)**2)



! write(88889,*) AR,CR,DR,cell%DRMESHDI,EK,config%ICST,PZ,QZ,FZ,SZ, &
!                           PNS,QNS,config%NSRA,vpottemp,cell%NPAN,cell%NRCUT, &
!                           gauntcoeff%CLEB,gauntcoeff%ICLEB,gauntcoeff%IEND,gauntcoeff%LOFLM,LMAXATOM,&
! !                           IDOLDAU,LOPT,LMLO,LMHI,&
!                           WLDAU(1,1,ISPIN),WLDAUAV,CUTOFF,&
!                           cell%nrmax, cell%nrmin_ns,gauntcoeff%ncleb, &
!                           lmaxatom,2*lmaxatom+1,lmmaxatom,(2*lmaxatom+1)**2

! stop

!     write(*,*) 'time wf', timing_stop()

do lval = 0,lmaxatom
    ekl(lval) = ek*dble(2*lval+1)
end do
!       write(*,*) 'test'
!  density%DEN = CZERO

!-----------------------------------------------------------------------
! calculate the charge density for spherical/non-spherical input potential
!-----------------------------------------------------------------------

!    if (ie==1) then
!      if(config_runflag('lmdos')) then
!       OPEN(UNIT=30, &
!            FILE="out_lmdos.atom="//char(48+IATOM/10)//char(48+mod(IATOM,10))//"_spin"//char(48+ISPIN)//".dat")
!       WRITE (30,*) ' '   ! lm-dos
!       WRITE (30,'(a8,I3,a4,I5)') '# ISPIN=',ISPIN,' IATOM=',IATOM   ! lm-dos
! 
!       WRITE(30,9000) DREAL(ENERG),(-DIMAG(DENLM(LM))/PI,LM=1,LMMAXD)
!  9000 FORMAT(30E12.4)
! 
!     end if
!   end if !ie==1

            if ( config%ins.eq.0 ) then
!              spherical
!                write(10000+ispin,*) 'ispin',ispin
!                write(10000+ispin,*) density%den(:,ispin,ie),df,gmatll,config%nsra, &
!                     density%rho2ns(:,:,ispin),cell%drmeshdi,cell%npan,cell%nrcut,pz,fz,qz,sz, &
!                     gauntcoeff%cleb,gauntcoeff%icleb,gauntcoeff%iend,gauntcoeff%jend,ekl, &
!                     cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2

               call rholm(density%den(:,ispin,ie),df,gmatll,config%nsra, &
                    density%rho2ns(:,:,ispin),cell%drmeshdi,cell%npan,cell%nrcut,pz,fz,qz,sz, &
                    gauntcoeff%cleb,gauntcoeff%icleb,gauntcoeff%iend,gauntcoeff%jend,ekl, &
                    cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2)

!                write(20000+ispin,*) 'ispin',ispin
!                write(20000+ispin,*) density%den(:,ispin,ie),df,gmatll,config%nsra, &
!                     density%rho2ns(:,:,ispin),cell%drmeshdi,cell%npan,cell%nrcut,pz,fz,qz,sz, &
!                     gauntcoeff%cleb,gauntcoeff%icleb,gauntcoeff%iend,gauntcoeff%jend,ekl, &
!                     cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2

            else
!              non-spherical
!       write(*,*) 'test2'
!                write(66666,*) density%den(:,ie),df,cell%drmeshdi,gmatll,ek, &
!                     density%rho2ns(:,:,ispin),cell%npan,cell%nrcut,shapefun%thetas,shapefun%lm2index,shapefun%lmused, &
!                     config%nsra,qns,pns,ar,cr,pz,fz,qz,sz,gauntcoeff%cleb(:,1),gauntcoeff%icleb, &
!                     gauntcoeff%jend,gauntcoeff%iend,ekl, &
!                     shapefun%nrshaped,shapefun%nlmshaped, cell%nrmin_ns, &
!                     cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2



               call rhons(density%den(:,ispin,ie),density%denlm(:,ispin,ie),df,cell%drmeshdi,gmatll,ek, &
                    density%rho2ns(:,:,ispin),cell%npan,cell%nrcut,shapefun%thetas,shapefun%lm2index,shapefun%lmused, &
                    config%nsra,qns,pns,ar,cr,pz,fz,qz,sz,gauntcoeff%cleb(:,1),gauntcoeff%icleb, &
                    gauntcoeff%jend,gauntcoeff%iend,ekl, &
                    shapefun%nrshaped,shapefun%nlmshaped, cell%nrmin_ns, &
                    cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2)
!       write(*,*) 'test'

!                write(66667,*) density%den(:,ie),df,cell%drmeshdi,gmatll,ek, &
!                     density%rho2ns(:,:,ispin),cell%npan,cell%nrcut,shapefun%thetas,shapefun%lm2index,shapefun%lmused, &
!                     config%nsra,qns,pns,ar,cr,pz,fz,qz,sz,gauntcoeff%cleb(:,1),gauntcoeff%icleb, &
!                     gauntcoeff%jend,gauntcoeff%iend,ekl, &
!                     shapefun%nrshaped,shapefun%nlmshaped, cell%nrmin_ns, &
!                     cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2

!                  write(*,*) density%rho2ns(:,:,ispin)
            end if
! stop
!     write(*,*) 'time rho', timing_stop()

! if (ie==ielast) then
!   if(config_runflag('lmdos')) then
!     close(30)
!   end if
! end if



if ( config_testflag('tmatdebug') ) then

write(5001,'(5000F)') cell%rmesh(cell%nrmin_ns:cell%nrmax)
write(5002,'(5000F)') cell%rmesh(:)

do lm1=1,(lmaxatom+1)**2
  do lm2=1,(lmaxatom+1)**2
    write(4002,'(5000E)') pns(lm2,lm1,:,1)
  end do
end do

do lm1=1,(lmaxatom+1)**2
  do lm2=1,(lmaxatom+1)**2
    write(4003,'(5000F)') qns(lm2,lm1,:,1)
  end do
end do

do lm1=0,lmaxatom
    write(4004,'(5500E)') pz(:,lm1)
end do

do lm1=0,lmaxatom
    write(4005,'(5500E)') qz(:,lm1)
end do

do lm1=0,lmaxatom
    write(4006,'(5500E)') fz(:,lm1)
end do

do lm1=0,lmaxatom
    write(4007,'(5500E)') sz(:,lm1)
end do


do lm1=1,(lmaxatom+1)**2
  do lm2=1,(lmaxatom+1)**2
    write(4008,'(5000E)') pns(lm2,lm1,:,2)
  end do
end do

do lm1=1,(lmaxatom+1)**2
  do lm2=1,(lmaxatom+1)**2
    write(4009,'(5000F)') qns(lm2,lm1,:,2)
  end do
end do


end if


! stop 'rhoval end'


            do lval = 0,lmaxatom+1
          density%ncharge(lval,ispin) = density%ncharge(lval,ispin) + DIMAG(density%den(lval,ispin,ie)*df)
            end do


            do lval = 0,lmaxatom+1
               energyparts%espv(lval,ispin,iatom) = energyparts%espv(lval,ispin,iatom) + dimag( (eryd-efermi)*density%den(lval,ispin,ie)*df)
            end do




!            end if

!----------------------------------------------------------------------
!--- LDA+U Calculation to implement in the future
!----------------------------------------------------------------------
! 
!          IF ( ( IDOLDAU.EQ.1 ).AND.( LOPT.GE.0 ) ) 
!      &        CALL DENSITYMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL(1,1,IE),
!      &                        IPAN,IRCUT,DRDI,EK,
!      &                        IRMIN,LOPT,MMAX,LMLO,LMHI,PHILDAU,DENMATC
!      &        ,den,ie) 
!         END DO

end subroutine rhoval

end module mod_rhoval
