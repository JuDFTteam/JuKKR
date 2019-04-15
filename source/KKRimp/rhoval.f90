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
                         LMAXATOM,LMMAXATOM,config,lmaxd_global,energyparts,efermi)
! C
! C **********************************************************************
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
      use nrtype, only: dp
      use mod_physic_params, only: cvlight
      use global_variables, only: irid, nfund, lmpotd, irmind, ipand, ncleb, lmaxd, irmd, lmmaxd
      use type_cell
      use type_gauntcoeff
      use type_shapefun
      use type_density
      use type_config
      use type_energyparts
      use mod_wfmesh
      use mod_cradwf
      use mod_pnsqns
      use mod_rholm
      use mod_rhons, only: rhons
      use mod_timing
      use mod_config, only: config_runflag,config_testflag
      implicit none
!interface variables
      integer,intent(in)                    ::   lmaxd_global
      integer,intent(in)                    ::   ie
      integer                               ::   ielast
      complex(kind=dp),intent(in)          ::   eryd
      complex(kind=dp),intent(in)          ::   wez
      complex(kind=dp),intent(in)          ::   gmatll(lmmaxatom,lmmaxatom)
      integer,intent(in)                    ::   ispin
      integer,intent(in)                    ::   nspin
      integer,intent(in)                    ::   iatom
      type(cell_type),intent(in)            ::   cell
      real(kind=dp),intent(in)              ::   vpot(1:cell%nrmaxd,(2*lmaxd_global+1)**2)
      type(shapefun_type),intent(in)        ::   shapefun
      type(gauntcoeff_type),intent(in)      ::   gauntcoeff
      real(kind=dp),intent(in)              ::   zatom
      type(density_type)                    ::   density
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

      complex(kind=dp),parameter           ::   czero=(0.0D0,0.0D0)
      integer                               ::   lm1,lm2,lval
      complex(kind=dp)                     ::   df,ek
      logical,parameter                     ::   lirrsol=.true.

      complex(kind=dp)                     ::   alpha(0:lmaxatom),ar(lmmaxatom,lmmaxatom), &
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
!          WLDAUAV = WLDAUAV/real(MMAX, kind=dp)
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

DF = WEZ/real(NSPIN, kind=dp)

IDOLDAU=0
if (idoldau==0) then 
  allocate(cutoff(cell%nrmax), wldau(2*lmaxatom+1,2*lmaxatom+1,nspin))
else
 stop 'ldau not implemented'
end if

allocate(vpottemp(cell%nrmin_ns:cell%nrmax,(2*lmaxatom+1)**2))
!-----------------------------------------------------------------------
! calculate the spherical wavefunctions
!-----------------------------------------------------------------------
call wfmesh(eryd,ek,cvlight,config%nsra,zatom,cell%rmesh,s,rs, &
  cell%nrmax,lmaxatom)
call cradwf(eryd,ek,config%nsra,alpha,cell%npan,cell%nrcut,cvlight,rs,s, &
  pz,fz,qz,sz,tmat,vpot(:,1),cell%drmeshdi,cell%rmesh,zatom,lirrsol, &
  idoldau,lopt,wldauav,cutoff,lmaxatom,lmaxatom+1,cell%nrmax)

!-----------------------------------------------------------------------
! transform to non-spherical wavefunctions
!-----------------------------------------------------------------------


IF (config%INS.GT.0) then
  vpottemp=vpot(cell%nrmin_ns:cell%nrmax,1:(2*lmaxatom+1)**2)

  CALL PNSQNS(AR,CR,DR,cell%DRMESHDI,EK,config%ICST,PZ,QZ,FZ,SZ, &
    PNS,QNS,config%NSRA,vpottemp,cell%NPAN,cell%NRCUT, &
    gauntcoeff%CLEB,gauntcoeff%ICLEB,gauntcoeff%IEND,gauntcoeff%LOFLM,LMAXATOM,&
    IDOLDAU,LOPT,LMLO,LMHI,&
    WLDAU(:,:,ISPIN),WLDAUAV,CUTOFF,&
    cell%nrmax, cell%nrmin_ns,gauntcoeff%ncleb, &
    lmaxatom,2*lmaxatom+1,lmmaxatom,(2*lmaxatom+1)**2)
end if ! config%ins>0

do lval = 0,lmaxatom
    ekl(lval) = ek*real(2*lval+1, kind=dp)
end do

!-----------------------------------------------------------------------
! calculate the charge density for spherical/non-spherical input potential
!-----------------------------------------------------------------------

if ( config%ins.eq.0 ) then
  ! spherical
  call rholm(density%den(:,ispin,ie),df,gmatll,config%nsra, &
    density%rho2ns(:,:,ispin),cell%drmeshdi,cell%npan,cell%nrcut,pz,fz,qz,sz, &
    gauntcoeff%cleb,gauntcoeff%icleb,gauntcoeff%iend,gauntcoeff%jend,ekl, &
    cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2)

else
  ! non-spherical

  !call rhons(density%den(:,ispin,ie),density%denlm(:,ispin,ie),df,cell%drmeshdi,gmatll,ek, &
  !     density%rho2ns(:,:,ispin),cell%npan,cell%nrcut,shapefun%thetas,shapefun%lm2index,shapefun%lmused, &
  !     config%nsra,qns,pns,ar,cr,pz,fz,qz,sz,gauntcoeff%cleb(:,1),gauntcoeff%icleb, &
  !     gauntcoeff%jend,gauntcoeff%iend,ekl, &
  !     shapefun%nrshaped,shapefun%nlmshaped, cell%nrmin_ns, &
  !     cell%nrmax,gauntcoeff%ncleb,lmaxatom,lmmaxatom,(2*lmaxatom+1)**2)

  !den = density%den(:,ispin,ie)
  !denlm = density%denlm(:,ispin,ie)
  !drdi = cell%drmeshdi
  !rho2ns = density%rho2ns(:,:,ispin)
  !ircut = cell%nrcut
  !thetas = shapefun%thetas
  !ifunm = shapefun%lm2index
  !lmsp = shapefun%lmused
  !nsra = config%nsra
  !cleb = gauntcoeff%cleb(:,1)
  !icleb = gauntcoeff%icleb
  !jend = gauntcoeff%jend
  !iend = gauntcoeff%iend
  !ipan = cell%npan

  ! set values of global_variables, used inside of rhons
  irid = shapefun%nrshaped
  nfund = shapefun%nlmshaped
  lmpotd = (2*lmaxatom+1)**2
  irmind = cell%nrmin_ns
  ipand = cell%npan
  ncleb = gauntcoeff%ncleb
  irmd = cell%nrmax
  lmmaxd = lmmaxatom
  lmaxd = lmaxatom

  call rhons(density%den(:,ispin,ie), df, cell%drmeshdi, gmatll, ek, density%rho2ns(:,:,ispin), cell%npan, cell%nrcut, cell%nrmin_ns, shapefun%thetas, shapefun%lm2index, shapefun%lmused, & 
    config%nsra, qns, pns, ar, cr, pz, fz, qz, sz, gauntcoeff%cleb(:,1), gauntcoeff%icleb, gauntcoeff%jend, gauntcoeff%iend, ekl, density%denlm(:,ispin,ie))

end if !( config%ins.eq.0 ) 



if ( config_testflag('tmatdebug') ) then

  write(5001,'(5000E25.14)') cell%rmesh(cell%nrmin_ns:cell%nrmax)
  write(5002,'(5000E25.14)') cell%rmesh(:)

  do lm1=1,(lmaxatom+1)**2
    do lm2=1,(lmaxatom+1)**2
      write(4002,'(5000E25.14)') pns(lm2,lm1,:,1)
    end do
  end do

  do lm1=1,(lmaxatom+1)**2
    do lm2=1,(lmaxatom+1)**2
      write(4003,'(5000E25.14)') qns(lm2,lm1,:,1)
    end do
  end do

  do lm1=0,lmaxatom
    write(4004,'(5500E25.14)') pz(:,lm1)
  end do

  do lm1=0,lmaxatom
    write(4005,'(5500E25.14)') qz(:,lm1)
  end do

  do lm1=0,lmaxatom
    write(4006,'(5500E25.14)') fz(:,lm1)
  end do

  do lm1=0,lmaxatom
    write(4007,'(5500E25.14)') sz(:,lm1)
  end do

  do lm1=1,(lmaxatom+1)**2
    do lm2=1,(lmaxatom+1)**2
      write(4008,'(5000E25.14)') pns(lm2,lm1,:,2)
    end do
  end do

  do lm1=1,(lmaxatom+1)**2
    do lm2=1,(lmaxatom+1)**2
      write(4009,'(5000E25.14)') qns(lm2,lm1,:,2)
    end do
  end do

end if


do lval = 0,lmaxatom+1
  density%ncharge(lval,ispin) = density%ncharge(lval,ispin) + aIMAG(density%den(lval,ispin,ie)*df)
end do


do lval = 0,lmaxatom+1
  energyparts%espv(lval,ispin,iatom) = energyparts%espv(lval,ispin,iatom) + aimag( (eryd-efermi)*density%den(lval,ispin,ie)*df)
end do

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
