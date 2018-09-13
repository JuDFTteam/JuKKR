!-------------------------------------------------------------------------------
! MODULE: MOD_MAIN2
!> @brief Wrapper module for the calculation of the DFT quantities for the JM-KKR package
!> @details The code uses the information obtained in the main0 module, this is
!> mostly done via the get_params_2() call, that obtains parameters of the type
!> t_params and passes them to local variables
!> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
module MOD_MAIN2

   use mod_Profiling
   use Constants
   use global_variables
   Use mod_datatypes, Only: dp

   use mod_brydbm
   use mod_convol
   use mod_ecoub
   use mod_epathtb
   use mod_epotinb
   use mod_espcb
   use mod_etotb1
   use mod_force
   use mod_forceh
   use mod_forcxc
   use mod_mtzero
   use mod_mdirnewang
   use mod_mixstr
   use mod_rhosymm
   use mod_relpotcvt
   use mod_rhototb
   use mod_vmadelblk
   use mod_vintras
   use mod_vinterface
   use mod_scfiterang
   use mod_rites
   use mod_writekkrflex
   use mod_vxcdrv
  use mod_rinit

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: main2
   !> @brief Main wrapper routine dealing with the calculation of the DFT quantities
   !> @details Calculates the potential from density, exc-potential, calculate total energy, ...
   !> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
   !> and many others ...
   !----------------------------------------------------------------------------
   subroutine main2()

      use mod_types, only: t_inc
      use mod_wunfiles
#ifdef CPP_TIMING
      use mod_timing
#endif
      use mod_version_info

      use mod_main0

      implicit none

      real (kind=dp), parameter :: eps=1.0D-12

      integer :: IOBROY
      parameter ( IOBROY = 20 )
      integer :: NMVECMAX
      parameter (NMVECMAX = 4)
      ! ..
      ! .. Local scalars
      integer :: NK  ! > ITERMDIR variables
      integer :: IRC1
      integer :: IPOT
      integer :: NMVEC  ! > ITERMDIR variables
      integer :: ICONT
      integer :: ISPIN
      integer :: IRMIN1
      integer :: LMAXD1
      integer :: LSMEAR
      integer :: LRECABMAD
      integer :: i_stat,i_all
      integer :: I,J,IE,I1,I2,IH,IT,IO,LM,IR,IREC
      real (kind=dp) :: DF
      real (kind=dp) :: RV
      real (kind=dp) :: MIX
      real (kind=dp) :: SUM
      real (kind=dp) :: FPI
      real (kind=dp) :: RFPI
      real (kind=dp) :: EFOLD
      real (kind=dp) :: EFNEW
      real (kind=dp) :: DENEF
      real (kind=dp) :: FSOLD
      real (kind=dp) :: VSHIFT  ! fxf
      real (kind=dp) :: RMSAVM
      real (kind=dp) :: RMSAVQ
      real (kind=dp) :: RMSAV0
      real (kind=dp) :: CHRGNT
      real (kind=dp) :: CHRGOLD
      real (kind=dp) :: EXCDIFF   ! > Scale magn. part of xc-potential
      real (kind=dp) :: E2SHIFT
      real (kind=dp) :: ERRAVANG  ! > ITERMDIR variables
      real (kind=dp) :: CHRGSEMICORE
      ! .. Local Arrays
      integer, dimension(NATYPD)                    :: LCOREMAX
      integer, dimension(NATYPD,NAEZD)               :: ITOQ
      integer, dimension(20,NATYPD)                 :: NKCORE
      integer, dimension(20,NPOTD)                 :: KAPCORE
      real (kind=dp), dimension(NATYPD)           :: EU  ! > LDA+U
      real (kind=dp), dimension(NATYPD)           :: EDC ! > LDA+U
      real (kind=dp), dimension(LMPOTD)           :: C00
      real (kind=dp), dimension(LMPOTD)           :: BVMAD
      real (kind=dp), dimension(NATYPD)           :: DENEFAT
      real (kind=dp), dimension(2)               :: VMT_INIT
      real (kind=dp), dimension(IRMD,NPOTD)       :: RHOC     ! > core charge density
      real (kind=dp), dimension(LMPOTD,LMPOTD)     :: AVMAD
      real (kind=dp), dimension(0:LPOTD,NATYPD)    :: EXCNM    ! > Scale magn. part of xc-potential
      real (kind=dp), dimension(LMPOTD,NAEZD)      :: VINTERS
      real (kind=dp), dimension(IRMD,NPOTD)       :: VSPSMDUM
      logical, dimension(NATYPD,LMPOT)              :: LPOTSYMM
      !-------------------------------------------------------------------------
      ! ITERMDIR variables
      !-------------------------------------------------------------------------
      real (kind=dp), dimension(NATYPD,NMVECMAX) :: MVGAM
      real (kind=dp), dimension(NATYPD,NMVECMAX) :: MVPHI
      real (kind=dp), dimension(NATYPD,NMVECMAX) :: MVTET
      complex (kind=dp), dimension(NATYPD,3,NMVECMAX) :: MVEVI
      complex (kind=dp), dimension(NATYPD,3,NMVECMAX) :: MVEVIEF
      !-------------------------------------------------------------------------
      !   ECOU(0:LPOT,NATYP)      ! Coulomb energy
      !   EPOTIN(NATYP),          ! energy of input potential (EPOTINB
      !   ESPC(0:3,NPOTD),        ! energy single particle core
      !   ESPV(0:LMAXD1,NPOTD)    ! energy single particle valence
      !                           ! both changed for the relativistic case
      !   EXC(0:LPOT,NATYP),      ! E_xc
      !-------------------------------------------------------------------------
      real (kind=dp), dimension(NATYPD)           :: EPOTIN   ! > energy of input potential (EPOTINB
      real (kind=dp), dimension(0:3,NPOTD)       :: ESPC     ! > energy single particle core
      real (kind=dp), dimension(0:LPOTD,NATYPD)    :: EXC      ! > exchange correlation energy
      real (kind=dp), dimension(0:LPOTD,NATYPD)    :: ECOU     ! > Coulomb energy
      real (kind=dp), dimension(0:LMAXD+1,NPOTD)  :: ESPV     ! > energy single particle valence both changed for the relativistic case
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD)  :: RHOORB
      real (kind=dp), dimension(KREL*20+(1-KREL),NPOTD)   :: ECOREREL
      !-------------------------------------------------------------------------
      !  CMINST(LMPOT,NATYP)            ! charge moment of interstitial
      !  CMOM(LMPOT,NATYP)              ! LM moment of total charge
      !  CHRGATOM(NATYP,
      !           2*KREL+(1-KREL)*NSPIND) ! total charge per atom
      !-------------------------------------------------------------------------
      real (kind=dp), dimension(LMPOTD,NATYPD)                    :: CMOM        ! > LM moment of total charge
      real (kind=dp), dimension(LMPOTD,NATYPD)                    :: CMINST      ! > charge moment of interstitial
      real (kind=dp), dimension(NATYPD,2*KREL+(1-KREL)*NSPIND)   :: CHRGATOM    ! > total charge per atom
      !-------------------------------------------------------------------------
      ! FORCES
      !-------------------------------------------------------------------------
      real (kind=dp), dimension(-1:1,NATYPD) :: FLM  ! > Forces
      real (kind=dp), dimension(-1:1,NATYPD) :: FLMC ! > Forces
      !-------------------------------------------------------------------------
      ! For SIMULASA
      !-------------------------------------------------------------------------
      integer :: IPOS,ILM_MAPP,IAS

      ! .. Allocatable arrays
      real (kind=dp), dimension(:,:,:), allocatable :: VONS ! > output potential (nonspherical VONS)

      !-------------------------------------------------------------------------
      !  R2NEF (IRMD,LMPOT,NATYP,2)  ! rho at FERMI energy
      !  RHO2NS(IRMD,LMPOT,NATYP,2)  ! radial density
      !   nspin=1            : (*,*,*,1) radial charge density
      !   nspin=2 or krel=1  : (*,*,*,1) rho(2) + rho(1) -> charge
      !                               (*,*,*,2) rho(2) - rho(1) -> mag. moment
      !  RHOC(IRMD,NPOTD)              ! core charge density
      !-------------------------------------------------------------------------
      real (kind=dp), dimension(:,:,:,:), allocatable :: R2NEF  ! > rho at FERMI energy
      real (kind=dp), dimension(:,:,:,:), allocatable :: RHO2NS ! > radial density
      !-------------------------------------------------------------------------
      ! Scale magn. part of xc-potential:
      real (kind=dp), dimension(:,:,:), allocatable    :: VXCM
      real (kind=dp), dimension(:,:,:), allocatable    :: VXCNM
      real (kind=dp), dimension(:,:,:,:), allocatable  :: RHO2NSNM

      ! .. External Subroutines
      logical :: OPT
      logical :: TEST


      LMAXD1=LMAX+1

      ! Allocations
      !allocate(THETAS(IRID,NFUND,NCELLD),stat=i_stat)
      !call memocc(i_stat,product(shape(THETAS))*kind(THETAS),'THETAS','main2')
      !allocate(THESME(IRID,NFUND,NCELLD),stat=i_stat)
      !call memocc(i_stat,product(shape(THESME))*kind(THESME),'THESME','main2')
      allocate(VONS(IRMD,LMPOT,NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(VONS))*kind(VONS),'VONS','main2')
      !allocate(VINS(IRMIND:IRMD,LMPOT,NSPOTD),stat=i_stat)
      !call memocc(i_stat,product(shape(VINS))*kind(VINS),'VINS','main2')
      allocate(VXCM(IRMD,LMPOT,NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(VXCM))*kind(VXCM),'VXCM','main2')
      allocate(VXCNM(IRMD,LMPOT,NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(VXCNM))*kind(VXCNM),'VXCNM','main2')
      allocate(R2NEF(IRMD,LMPOT,NATYP,2),stat=i_stat)
      call memocc(i_stat,product(shape(R2NEF))*kind(R2NEF),'R2NEF','main2')
      allocate(RHO2NS(IRMD,LMPOT,NATYP,2),stat=i_stat)
      call memocc(i_stat,product(shape(RHO2NS))*kind(RHO2NS),'RHO2NS','main2')
      allocate(RHO2NSNM(IRMD,LMPOT,NATYP,2),stat=i_stat)
      call memocc(i_stat,product(shape(RHO2NSNM))*kind(RHO2NSNM),'RHO2NSNM','main2')

      ! Consistency check
      if ( (KREL.lt.0) .or. (KREL.gt.1) ) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
      if ( (KREL.eq.1) .and. (NSPIND.eq.2) ) stop ' set NSPIN = 1 for KREL = 1 in the inputcard'
      !-------------------------------------------------------------------------
      ! This routine previously used to read from unformatted files created by
      ! the main0 module, now  instead of unformatted files take parameters from
      ! types defined in wunfiles.F90
      !-------------------------------------------------------------------------
      call get_params_2(t_params,KREL,NATYP,IPAND,NPOTD,NATOMIMPD,LMXSPD,NFUND,  &
         LMPOT,NCELLD,IRMD,NEMBD1,NEMBD,IRMIND,NSRA,INS,NSPIN,IPAN,IRCUT,LCORE,    &
         NCORE,LMAX,NTCELL,LPOT,NLBASIS,NRBASIS,NRIGHT,NLEFT,NATOMIMP,ATOMIMP,   &
         IMIX,QBOUND,FCM,ITDBRY,IRNS,KPRE,KSHAPE,KTE,KVMAD,KXC, ICC,ISHIFT,      &
         IXIPOL,KFORCE,IFUNM,LMSP,IMT,IRC,IRMIN,IRWS,LLMSP,ITITLE,NFU,HOSTIMP,   &
         ILM_MAP,IMAXSH,IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,IESEMICORE,     &
         KAOEZ,IQAT,NOQ,LLY,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,ZREL,JWSREL,IRSHIFT,   &
         MIXING,LAMBDA_XC,A,B,THETAS, DRDI,RMESH,ZAT,RMT,RMTNEW,RWS,EMIN,EMAX,TK,    &
         ALAT,EFOLD,CHRGOLD,CMOMHOST,CONC,GSH,EBOTSEMI,EMUSEMI,TKSEMI,VINS,VISP, &
         RMREL,DRDIREL,VBC,FSOLD,R2DRDIREL,ECORE,EZ,WEZ,TXC,LINTERFACE,LRHOSYM,  &
         NGSHD,NAEZ,IRID,NSPOTD,IEMXD)

      !-------------------------------------------------------------------------
      ! Reading the density parameters stored in t_params
      !-------------------------------------------------------------------------
      call read_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,ESPV,ECORE,  &
         IDOLDAU,LOPT,EU,EDC,CHRGSEMICORE,RHOORB,ECOREREL,NKCORE,KAPCORE,KREL,&
         NATYP,NPOTD,IRMD,LMPOT,LMAXD1)
      !-------------------------------------------------------------------------
      ! End read in variables
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Setting up constants
      !-------------------------------------------------------------------------
      FPI = 4.0D0*PI
      RFPI = SQRT(FPI)
      RMSAV0 = 1.0D10
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! Setting dummy argument LSMEAR to allow compatibility with IMPURITY
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      LSMEAR=0
      !
      ICONT = 1
      IPF = 1337
      NSPIN = 2*KREL + (1-KREL)*NSPIN
      IDOSEMICORE = 0
      if ( OPT('SEMICORE') ) IDOSEMICORE = 1
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ITERATION BEGIN  ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ITSCF = ITSCF + 1         ! initialised to 0 in main0
      t_inc%i_iteration = ITSCF
      !
      write(1337,'(/,79("*"))')
      write(1337,'(19X,A,I3,A,I3,A)') '****** ITERATION : ',ITSCF,' OUT OF ',SCFSTEPS,' ******'
      write(1337,'(79("*"),/)')
      !
      if (IMIX.ge.3) then
         open (IOBROY,FILE='broy_io.unformatted',FORM='unformatted',STATUS='unknown')
         open (IOBROY+2,FILE='broy_io2.unformatted',FORM='unformatted',STATUS='unknown')
      endif
      !-------------------------------------------------------------------------
      ! the next four lines may not always work
      !-------------------------------------------------------------------------
      NSHELL(0) = NATYP
      do I1 = 1,NATYP
         NSHELL(I1) = 1
      end do
      !-------------------------------------------------------------------------
      ! Determine total charge density expanded in spherical harmonics
      !-------------------------------------------------------------------------
      if(TEST('flow    ')) write(1337,*) '>>> RHOTOTB'
      call RHOTOTB(IPF,NATYP,NAEZ,NSPIN,RHO2NS,RHOC,RHOORB,ZAT,DRDI,IRWS,IRCUT,  &
         NFU,LLMSP,THETAS,NTCELL,KSHAPE,IPAN,CHRGNT,ITSCF,NSHELL,NOQ,CONC,  &
         KAOEZ,CHRGATOM,IRMD,NEMB,LMPOT)

      if(TEST('flow    ')) write(1337,*) '<<< RHOTOTB'

      if ( TEST('RHOVALTW') ) then !Bauer
         do I1 = 1,NATYP
            open(unit=324234,file='out_rhotot')
            write(324234,*) '#IATOM',I1
            write(324234,'(50000F14.7)') RHO2NS(:,:,I1,1)
            if (NSPIN==2) write(324234,'(50000F14.7)') RHO2NS(:,:,I1,2)
         end do
      end if

      !-------------------------------------------------------------------------
      ! Determine new Fermi level due to valence charge up to old Fermi level
      ! EMAX and density of states DENEF
      !-------------------------------------------------------------------------
      if (ITSCF.gt.1.and.CHRGNT*CHRGOLD.lt.0.D0.and.ABS(CHRGNT).gt.5.D-2) then
         E2SHIFT = CHRGNT/(CHRGNT-CHRGOLD)*(EMAX-EFOLD)
      else
         E2SHIFT = CHRGNT/DENEF
      end if
      !
      E2SHIFT = MIN(ABS(E2SHIFT),0.05_dp)*SIGN(1.0_dp,E2SHIFT)
      EFOLD = EMAX
      CHRGOLD = CHRGNT
      if (TEST('no-neutr').OR.OPT('no-neutr')) then
         write(1337,*) 'test-opt no-neutr: Setting FERMI level shift to zero'
         E2SHIFT = 0.d0
      endif
      if (TEST('slow-neu').OR.OPT('slow-neu')) then
         write(1337,*) 'test-opt slow-neu: FERMI level shift * STRMIX'
         E2SHIFT = E2SHIFT * MIXING
      endif
      if (ISHIFT.EQ.0) EMAX = EMAX - E2SHIFT
      !-------------------------------------------------------------------------
      FSEMICORE=0d0
      if ( IDOSEMICORE.eq.1 ) then
         !----------------------------------------------------------------------
         ! Semicore treatment, recalculate the normalisation factor
         !----------------------------------------------------------------------
         if ( CHRGSEMICORE.lt.1D-10 ) CHRGSEMICORE = 1D-10
         !  Number of semicore bands
         I1 = NINT(CHRGSEMICORE)
         FSEMICORE = real(I1, kind=dp)/CHRGSEMICORE * FSOLD
         write(1337,'(6X,"< SEMICORE > : ",/,21X,"charge found in semicore :",F10.6,/,21X,"new normalisation factor :",F20.16,/)')&
            CHRGSEMICORE,FSEMICORE
      end if
      !-------------------------------------------------------------------------
      !write (6,FMT=9020) EFOLD,E2SHIFT
      write (1337,FMT=9020) EFOLD,E2SHIFT
      !-------------------------------------------------------------------------
      ! Divided by NAEZ because the weight of each atom has been already
      !     taken into account in 1c
      !-------------------------------------------------------------------------
      write (1337,FMT=9030) EMAX,DENEF/real(NAEZ, kind=dp)
      write (6,FMT=9030) EMAX,DENEF/real(NAEZ, kind=dp)
      write(1337,'(79("+"),/)')
      !-------------------------------------------------------------------------
      DF = 2.0D0/PI*E2SHIFT/real(NSPIN, kind=dp)
      !-------------------------------------------------------------------------
      ! ISPIN LOOP
      !-------------------------------------------------------------------------
      do ISPIN = 1,NSPIN
         !----------------------------------------------------------------------
         if (KTE.EQ.1) then
            do I1 = 1,NATYP
               IPOT = (I1-1)*NSPIN + ISPIN
               ESPV(0,IPOT) = ESPV(0,IPOT) -EFOLD*CHRGNT/real(NSPIN*NAEZ, kind=dp)
            end do
         end if                 ! (kte.eq.1)
         !----------------------------------------------------------------------
         ! Get correct density
         !----------------------------------------------------------------------
         if (.not.(OPT('DECIMATE'))) then
            do I1 = 1,NATYP
               do LM = 1,LMPOT
                  call DAXPY(IRC(I1),DF,R2NEF(1,LM,I1,ISPIN),1,&
                     RHO2NS(1,LM,I1,ISPIN),1)
               end do
            end do
         end if
         !----------------------------------------------------------------------
      end do
      !-------------------------------------------------------------------------
      ! End of ISPIN loop
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! ITERMDIR
      !-------------------------------------------------------------------------
      if ((KREL.eq.1).and.(OPT('ITERMDIR'))) then
         MVEVI   = t_params%MVEVI
         MVEVIEF = t_params%MVEVIEF
         !
         call RINIT(NAEZ,QMGAM)
         do I1 = 1,NAEZ
            ITOQ(1,I1) = KAOEZ(1,I1)
         end do
         NK = 2 * LMAX + 1
         NMVEC = 2
         !
         FACT(0) = 1.0D0
         do I = 1,100
            FACT(I) = FACT(I-1)*real(I, kind=dp)
         end do
         !----------------------------------------------------------------------
         if (.not.(OPT('DECIMATE'))) then
            do I1 = 1,NATYP
               do LM = 1, NMVEC
                  do IT = 1,3
                     MVEVI(I1,IT,LM) = MVEVI(I1,IT,LM)+ E2SHIFT*MVEVIEF(I1,IT,LM)
                  end do
               end do
            end do
         end if
         !----------------------------------------------------------------------
         do I1 = 1,NATYP
            call MDIRNEWANG(I1,NMVEC,MVEVI,MVPHI,MVTET,MVGAM,NATYP,LMAX,NMVECMAX)
         end do
         !
         open(67,FILE='itermdir.unformatted',FORM='unformatted')
         call SCFITERANG(ITSCF,ITOQ,FACT,MVPHI,MVTET,MVGAM,QMPHI,QMTET,QMGAM,&
            NAEZ,NK,ERRAVANG,NAEZ,NATYP,NMVECMAX,LMMAXD)
         t_params%MVEVI   = MVEVI
         t_params%MVEVIEF = MVEVIEF
      end if
      !-------------------------------------------------------------------------
      ! End of ITERMDIR
      !-------------------------------------------------------------------------

      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! POTENTIAL PART
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      if (LRHOSYM) then
         call RHOSYMM(LMPOT,NSPIN,1,NATYP,RHO2NS,IXIPOL,IRWS,IRCUT,IPAN,KSHAPE,&
            NATYP,IRMD)
      endif

      CMINST(:,:) = 0.d0
      CMOM(:,:) =0.d0
      VONS(:,:,:) = 0.d0
      call VINTRAS(CMOM,CMINST,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,       &
         RMESH,DRDI,IRWS,IRCUT,IPAN,KSHAPE,NTCELL,ILM_MAP,IFUNM,IMAXSH,GSH,  &
         THETAS,LMSP,LMPOT,NATYP)

      if ( TEST('vintrasp') ) then !Bauer
         open(unit=786785,file='test_vintraspot')
         do i1=1,nspin*natyp
            write(786785,*) '# atom/spin index ',i1
            write(786785,'(50000E25.16)') vons(:,:,i1)
         end do !iatom
         close(786785)
      end if !

      !-------------------------------------------------------------------------
      !fivos     IF ( .NOT.TEST('NoMadel ').AND. ( SCFSTEPS.GT.1 )
      !fivos     &     .OR. (ICC .GT. 0 ) )THEN
      !-------------------------------------------------------------------------
      if ( LINTERFACE ) then
         call VINTERFACE(CMOM,CMINST,LPOT,NSPIN,NAEZ,NATYP,VONS,ZAT,RMESH,IRWS,IRCUT,&
            IPAN,KSHAPE,NOQ,KAOEZ,IQAT,CONC,CHRGATOM(1,1),ICC,HOSTIMP,NLBASIS,   &
            NLEFT,NRBASIS,NRIGHT,CMOMHOST,CHRGNT,VINTERS,NAEZ,LMPOT)
         !----------------------------------------------------------------------
      else
         !----------------------------------------------------------------------
         call VMADELBLK(CMOM,CMINST,LPOT,NSPIN,NAEZ,VONS,ZAT,RMESH,IRWS,IRCUT,IPAN,  &
            KSHAPE,NOQ,KAOEZ,CONC,CHRGATOM(1,1),ICC,HOSTIMP,VINTERS,NEMB,    &
            LMPOT,NATYP)
      end if

      if (OPT('KKRFLEX ')) then
         call WRITEKKRFLEX(NATOMIMP,NSPIN,IELAST,(LPOT+1)**2,ALAT,NATYP,  &
            KSHAPE,VBC,ATOMIMP,HOSTIMP,NOQ,ZAT,KAOEZ,CONC,CMOM,CMINST,VINTERS,   &
            NEMB,NAEZ)
      end if

      !-------------------------------------------------------------------------
      !fivos      END IF
      !-------------------------------------------------------------------------
      if ( TEST('Vspher  ') ) VONS(1:IRMD,2:LMPOT,1:NPOTD) = 0.D0
      if ( TEST('vpotout ') ) then !bauer
         open(unit=54633163,file='test_vpotout_inter')
         do i1=1,natyp*nspin
            write(54633163,*) '# atom ',i1
            write(54633163,'(50000E25.16)') vons(:,:,i1)
         end do !iatom
         close(54633163)
      end if ! config_testflag('write_gmatonsite')

      !-------------------------------------------------------------------------
      ! Write the CMOMS to a file
      !-------------------------------------------------------------------------
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! In case of DECIMATION output, we store ONLY information connected
      ! with the effective CPA-medium (for the case of NO-CPA NAEZ=NATYP)
      ! hence the CMOMS are calculated site-dependent. In the same format
      ! are read in by <MAIN0> -- < CMOMSREAD >     v.popescu 01/02/2002
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      if (OPT('deci-out').and.(ITSCF.eq.1)) then
         open(37, file='decifile', form='formatted', position='append')
         write(37,1080) NAEZ,LMPOT
         do IH=1,NAEZ
            write(37,*) IH
            do LM=1,LMPOT
               C00(LM) = 0.0D0
               !----------------------------------------------------------------
               ! Store the charge on SITE IH
               !----------------------------------------------------------------
               do IO=1,NOQ(IH)
                  IT = KAOEZ(IO,IH)
                  C00(LM) = C00(LM) + CMOM(LM,IT) * CONC(IT)
                  if (INS.NE.0) C00(LM) = C00(LM)+CMINST(LM,IT) * CONC(IT)
                  if (LM.eq.1) C00(1) = C00(1) - ZAT(IT)/RFPI*CONC(IT)
               end do
               !----------------------------------------------------------------
            end do
            write(37,1090) (C00(LM),LM=1,LMPOT)
         end do
         close(37)
      end if
      !-------------------------------------------------------------------------
      ! FORCES
      !-------------------------------------------------------------------------
      if ( (KFORCE.eq.1).and.(KREL.ne.1) ) then
         if (INS.eq.0) then
            call FORCEH(CMOM,FLM,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,RMESH,DRDI,IRWS,ZAT)
            call FORCE(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,RMESH,DRDI,IRWS)
         else
            call FORCEH(CMOM,FLM,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,RMESH,DRDI,IMT,ZAT)
            call FORCE(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,RMESH,DRDI,IMT)
         end if
      end if
      !-------------------------------------------------------------------------
      ! Force Calculation stops here look after VXCDRV
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! ENERGIES
      !-------------------------------------------------------------------------
      if (KTE.eq.1) then
         ! Single-particle core energy
         call ESPCB(ESPC,NSPIN,NATYP,ECORE,LCORE,LCOREMAX,NCORE)
         ! "Energy of the input potential"
         ! Int V(r) rho(r) d^3r
         call EPOTINB(EPOTIN,NSPIN,NATYP,RHO2NS,VISP,RMESH,DRDI,INS,IRMIN,IRWS,   &
            LPOT,VINS,IRCUT,IPAN,ZAT)
         ! Coulomb hartree energy
         call ECOUB(CMOM,ECOU,LPOT,NSPIN,NATYP,RHO2NS,VONS,ZAT,RMESH,DRDI,IRWS,   &
            KVMAD,KSHAPE,IRCUT,IPAN,IMAXSH,IFUNM,ILM_MAP,NTCELL,GSH,THETAS,LMSP,LPOT)

      end if
      !-------------------------------------------------------------------------
      ! End of calculation of the energy
      !-------------------------------------------------------------------------
      VXCM(:,:,:) = 0.D0
      call VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,1,NATYP,RHO2NS,VXCM,RMESH,DRDI, &
         A,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,ILM_MAP,IMAXSH,IFUNM,THETAS,LMSP)

      if ( TEST('Vspher  ') ) VONS(1:IRMD,2:LMPOT,1:NPOTD) = 0.D0

      ! Recalculate XC-potential with zero spin density for magn. moment scaling
      VXCNM(:,:,:) = 0.D0                 ! Initialize
      EXCNM(:,:) = 0.D0
      if (abs(LAMBDA_XC-1.D0)>eps.and.NSPIN.eq.2) then
         RHO2NSNM(:,:,:,1) = RHO2NS(:,:,:,1) ! Copy charge density
         RHO2NSNM(:,:,:,2) = 0.D0            ! Set spin density to zero
         call VXCDRV(EXCNM,KTE,KXC,LPOT,NSPIN,1,NATYP,RHO2NSNM,VXCNM,&
            RMESH,DRDI,A,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,              &
            ILM_MAP,IMAXSH,IFUNM,THETAS,LMSP)
         ! Compute the EXC-difference
         EXCDIFF = 0.D0
         do I1 = 1,NATYP
            do LM = 0,LPOT
               EXCDIFF = EXCDIFF + EXC(LM,I1) - EXCNM(LM,I1)
            end do
         end do
         write(1337,*) 'LAMBDA_XC=',LAMBDA_XC,'EXCDIF=',EXCDIFF
      endif

      VONS(:,:,:) = VONS(:,:,:) +   & ! Add xc-potential with magn. part weighted by lambda_xc
         LAMBDA_XC*VXCM(:,:,:) + (1.D0-LAMBDA_XC)*VXCNM(:,:,:)
      EXC(:,:) = LAMBDA_XC*EXC(:,:)     + (1.D0-LAMBDA_XC)*EXCNM(:,:)


      if ( TEST('vpotout ') ) then !bauer
         open(unit=57633263,file='test_vpotout_xc')
         do i1=1,natyp*nspin
            write(57633263,*) '# atom ',i1
            write(57633263,'(50000E25.16)') vons(:,:,i1)
         end do !iatom
         close(57633263)
      end if ! config_testflag('write_gmatonsite')

      !-------------------------------------------------------------------------
      ! FORCES
      !-------------------------------------------------------------------------
      ! Force calculation continues here
      !-------------------------------------------------------------------------
      if ( (KFORCE.eq.1).and.(KREL.ne.1) ) then
         if (KSHAPE.eq.0) THEN
            call FORCXC(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,RMESH,ALAT,DRDI,IRWS,0)
         else
            call FORCXC(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,RMESH,ALAT,DRDI,IMT,0)
         end if
      end if
      !-------------------------------------------------------------------------
      ! Force calculation ends
      !-------------------------------------------------------------------------
      if (ISHIFT.eq.2) then                                                      ! fxf
         !        OPEN (67,FILE='vmtzero_init',FORM='formatted')                 ! fxf
         !        READ (67,*) VMT_INIT(1)                                        ! fxf
         !        CLOSE(67)                                                      ! fxf
         VMT_INIT(1) = 0.d0
         VMT_INIT(2) = VMT_INIT(1) ! read initial muffin-tin zero                ! fxf
         ! vmt_init is needed as a common reference for potential mixing if more
         ! iterations are to be remembered, e.g., in Anderson mixing.
         !----------------------------------------------------------------------
         ! Shift new potential to initial muffin-tin zero                        ! fxf
         call MTZERO(LMPOT,NATYP,CONC,NSPIN,VONS,VMT_INIT,ZAT,RMESH,DRDI,   &        ! fxf
            IMT,IRCUT,IPAN,NTCELL,LMSP,IFUNM,THETAS,IRWS,E2SHIFT,ISHIFT,&        ! fxf
            NSHELL,LINTERFACE)                                                   ! fxf
         open (67,FILE='vmtzero',FORM='formatted')                               ! fxf
         write (67,*) VMT_INIT(1)                                                ! fxf
         close(67)                                                               ! fxf
         ! Shift old potential to initial muffin-tin zero for correct mixing     ! fxf
         VSHIFT = - VBC(1)                                                       ! fxf
         call POTENSHIFT(VISP,VINS,NATYP,NSPIN,IRCUT,IRC,IRMIN,NTCELL,  &        ! fxf
            IMAXSH,ILM_MAP,IFUNM,LMSP,LMPOT,GSH,THETAS,THESME,RFPI,RMESH,KSHAPE,&        ! fxf
            VSHIFT,IRMD,NPOTD,IRMIND,LMXSPD)                                      ! fxf
      else if (ISHIFT.EQ.1) then
         ! Shift new potential to old MT-zero for correct mixing
         ! (convolution with shapes is done later)
         do ISPIN = 1,NSPIN
            do IH = 1,NATYP
               IPOT = NSPIN * (IH-1) + ISPIN
               IRC1 = IRC(IH)
               VSHIFT = RFPI * VBC(ISPIN)
               VONS(1:IRC1,1,IPOT) = VONS(1:IRC1,1,IPOT) + VSHIFT
            enddo
         enddo

      else                                                                       ! fxf
         ! Before fxf, only the following call was present.
         call MTZERO(LMPOT,NATYP,CONC,NSPIN,VONS,VBC,ZAT,RMESH,DRDI,  &
            IMT,IRCUT,IPAN,NTCELL,LMSP,IFUNM,THETAS,IRWS,E2SHIFT, &
            ISHIFT,NSHELL,LINTERFACE)
      END IF                                                                     ! fxf
      !-------------------------------------------------------------------------
      WRITE(1337,'(79("="),/)')
      !-------------------------------------------------------------------------
      ! Convolute potential with shape function for next iteration
      !-------------------------------------------------------------------------

      if ( TEST('vpotout ') ) then !bauer
         open(unit=12633269,file='test_vpotout_shift')
         do i1=1,natyp*nspin
            write(12633269,*) '# atom ',i1
            write(12633269,'(50000E25.16)') vons(:,:,i1)
         end do !iatom
         close(12633269)
      end if ! config_testflag('write_gmatonsite')

      if (KSHAPE.ne.0) then
         do ISPIN = 1,NSPIN
            do I1 = 1,NATYP
               IPOT = NSPIN* (I1-1) + ISPIN

               if ( TEST('vpotout ') ) then !bauer
                  open(unit=12642269,file='test_convol')
                  write(12642269,*) '# atom ',i1

                  write(12642269,*) IRCUT(1,I1),IRC(I1),IMAXSH(LMPOT),ILM_MAP, &
                     IFUNM,LMPOT,GSH,THETAS,ZAT(I1),RFPI,RMESH(:,I1),          &
                     VONS(:,:,IPOT),LMSP
                  close(12642269)
               end if           ! config_testflag('write_gmatonsite')

               call CONVOL(IRCUT(1,I1),IRC(I1),NTCELL(I1),IMAXSH(LMPOT),   &
                  ILM_MAP,IFUNM,LMPOT,GSH,THETAS,THESME,ZAT(I1),RFPI,RMESH(:,I1),  &
                  VONS(:,:,IPOT),VSPSMDUM(1,1),LMSP)
            end do
         end do
      end if

      if ( TEST('vpotout ') ) then !bauer
         open(unit=57633269,file='test_vpotout_conv')
         do i1=1,natyp*nspin
            write(57633269,*) '# atom ',i1
            write(57633269,'(50000E25.16)') vons(:,:,i1)
         end do !iatom
         close(57633269)
      end if ! config_testflag('write_gmatonsite')

      !-------------------------------------------------------------------------
      ! Symmetrisation of the potentials
      !-------------------------------------------------------------------------
      ! Keep only symmetric part of the potential
      if (TEST('potcubic')) then
         write(1337,*) 'Keeping only symmetric part of potential:'
         write(1337,*) 'Components L = 1, 11, 21, 25, 43, 47.'
         do IPOT = 1,NPOTD
            do LM = 1, LMPOT
               if (LM.ne.1.and.LM.ne.11.and.LM.ne.21  &
                  .and.LM.ne.25.and.LM.ne.43.and.LM.ne.47) then
                  do I = 1,IRMD
                     VONS(I,LM,IPOT) = 0.d0
                  enddo
               endif
            enddo
         enddo
      endif

      if(TEST('potsymm ')) then
         ! declarations needed:
         !     real (kind=dp) AVMAD(LMPOT,LMPOT),BVMAD(LMPOT)
         !     INTEGER LRECABMAD,I2,IREC
         !     LOGICAL LPOTSYMM(NATYP,LMPOT)
         LRECABMAD = WLENGTH*2*LMPOT*LMPOT + WLENGTH*2*LMPOT
         open (69,ACCESS='direct',RECL=LRECABMAD,FILE='abvmad.unformatted',   &
            FORM='unformatted')
         do I1 = 1,NATYP
            do LM = 1,LMPOT
               LPOTSYMM(I1,LM)=.FALSE.
            end do
            do I2 = 1,NATYP
               IREC = I2 + NAEZ*(I1-1)
               read (69,REC=IREC) AVMAD,BVMAD
               do LM = 1,LMPOT
                  if(ABS(BVMAD(LM)).GT.1D-10) LPOTSYMM(I1,LM)=.TRUE.
               end do
            end do
            do LM = 1,LMPOT
               if(LPOTSYMM(I1,LM)) then
                  write(1337,*) 'atom ',I1,'lm = ',LM,' contribution used'
               else
                  do ISPIN=1,NSPIN
                     IPOT = NSPIN* (I1-1) + ISPIN
                     do IR = 1,IRMD
                        VONS(IR,LM,IPOT) = 0.0D0
                     end do
                  end do
               endif
            end do
         end do
         close (69)
      end if

      if ( TEST('vpotout ') ) then !bauer
         open(unit=54633563,file='test_vpotout')
         do i1=1,natyp*nspin
            write(54633563,*) '# atom ',i1
            write(54633563,'(50000E25.16)') vons(:,:,i1)
         end do !iatom
         close(54633563)
      end if ! config_testflag('write_gmatonsite')

      ! for simulasa:
      if(opt('simulasa')) then
         do IAS = 1, NPOTD
            do ILM_MAPP = 1, LMPOT
               do IPOS = 1, IRMD
                  if (ILM_MAPP .NE. 1) then
                     VONS(IPOS, ILM_MAPP, IAS) = 0.D0
                  endif
               enddo
            enddo
         enddo
      endif

      !-------------------------------------------------------------------------
      ! Final construction of the potentials (straight mixing)
      !-------------------------------------------------------------------------
      MIX = MIXING
      if (TEST('alt mix ')) MIX = MIXING/real(1+MOD(ITSCF,2), kind=dp)
      if (TEST('spec mix')) then
         MIX = MIXING/(1.0D0 + 1.0D+3 * ABS(CHRGNT)/real(NAEZ*NSPIN, kind=dp))
      endif
      write(1337,*) 'MIXSTR',MIX
      call MIXSTR(RMSAVQ,RMSAVM,INS,LPOT,LMPOT,0,NSHELL, &
         1,NATYP,CONC,NSPIN,ITSCF,RFPI,FPI,IPF,MIX,FCM,  &
         IRC,IRMIN,RMESH,DRDI,VONS,VISP,VINS,VSPSMDUM,VSPSMDUM,LSMEAR)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! End of  POTENTIAL PART
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !
      !-------------------------------------------------------------------------
      if (ITSCF.NE.1) RMSAV0 = 1.0d2*MAX(RMSAVQ,RMSAVM)
      !
      write(1337,FMT=9160) MIX
      write(1337,'(79("="),/)')
      !-------------------------------------------------------------------------
      if (MAX(RMSAVQ,RMSAVM).LT.QBOUND) then
         t_inc%i_iteration = t_inc%N_iteration
      else
         !----------------------------------------------------------------------
         ! Potential mixing procedures: Broyden or Andersen updating schemes
         !----------------------------------------------------------------------
         if (IMIX.GE.3) then
            call BRYDBM(VISP,VONS,VINS,VSPSMDUM,VSPSMDUM,INS,  &
               LMPOT,RMESH,DRDI,MIX,CONC,IRC,IRMIN,NSPIN,1,NATYP,  &
               ITDBRY,IMIX,IOBROY,IPF,LSMEAR)
         endif
         !----------------------------------------------------------------------
         ! Reset to start new iteration
         !----------------------------------------------------------------------
         do I = 1,NSPIN*NATYP
            IT = I
            if (NSPIN.EQ.2) IT = (I+1)/2
            !
            IRC1 = IRC(IT)
            call DCOPY(IRC1,VONS(1,1,I),1,VISP(1,I),1)
            !
            if ( ( INS.NE.0 ).AND.( LPOT.GT.0 ) ) then
               IRMIN1 = IRMIN(IT)
               do LM = 2,LMPOT
                  do J = IRMIN1,IRC1
                     VINS(J,LM,I) = VONS(J,LM,I)
                  end do
                  SUM = 0.0D0
                  do IR = IRMIN1,IRC1
                     RV = VINS(IR,LM,I)*RMESH(IR,IT)
                     SUM = SUM + RV*RV*DRDI(IR,IT)
                  end do
                  if ( SQRT(SUM).LT.QBOUND ) then
                     do J = IRMIN1,IRC1
                        VINS(J,LM,I) = 0.0D0
                     end do
                  end if
               end do
            end if
         end do
         !----------------------------------------------------------------------
      end if
      !-------------------------------------------------------------------------
      rewind 11
      !
      EFNEW = EMAX
      if (OPT('rigid-ef').OR.OPT('DECIMATE')) EFNEW = EFOLD
      !
      if (ISHIFT.EQ.2) then ! Shift mixed potential to new muffin-tin zero       ! fxf
         VBC(1) = VBC(1) + E2SHIFT                                               ! fxf
         VBC(2) = VBC(1)                                                         ! fxf
         VSHIFT = VBC(1)                                                         ! fxf
         call POTENSHIFT(VISP,VINS,NATYP,NSPIN,IRCUT,IRC,IRMIN,NTCELL,  &        ! fxf
            IMAXSH,ILM_MAP,IFUNM,LMSP,LMPOT,GSH,THETAS,THESME,RFPI,RMESH,KSHAPE,&        ! fxf
            VSHIFT,IRMD,NPOTD,IRMIND,LMXSPD)                                      ! fxf
         write(1337,*) 'New VMT ZERO:',VBC(1)                                    ! fxf
      end if                                                                     ! fxf
      !
      call RITES(11,1,NATYP,NSPIN,ZAT,ALAT,RMT,RMTNEW,RWS,ITITLE, &
         RMESH,DRDI,VISP,IRWS,A,B,TXC,KXC,INS,IRNS,LPOT,VINS,         &
         QBOUND,IRC,KSHAPE,EFNEW,VBC,ECORE,LCORE,NCORE,           &
         ECOREREL,NKCORE,KAPCORE,LMPOT)
      close (11)
      !-------------------------------------------------------------------------
      ! ENERGIES calculation
      !-------------------------------------------------------------------------
      if ((KTE.EQ.1 .AND. ICC.EQ.0) .or. OPT('KKRFLEX ')) then
         call ETOTB1(ECOU,EPOTIN,ESPC,ESPV,EXC,KPRE,LMAX,LPOT,   &
            LCOREMAX,NSPIN,NATYP,NSHELL(1),CONC,IDOLDAU,LOPT,EU,EDC)
      endif
      !-------------------------------------------------------------------------
      ! End of ENERGIES calculation
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! CONVERGENCY TESTS
      !-------------------------------------------------------------------------
      if ( OPT('SEARCHEF') .and. (ABS(E2SHIFT).LT.1D-8)) then
         t_inc%i_iteration = t_inc%N_iteration
         ICONT = 0
         goto 260
      end if
      !-------------------------------------------------------------------------
      if (MAX(RMSAVQ,RMSAVM).LT.QBOUND) then
         write(6,'(17X,A)') '++++++ SCF ITERATION CONVERGED ++++++'
         write(6,'(79("*"))')
         ICONT = 0
         go to 260
         !----------------------------------------------------------------------
      else
         !---------------------------------------------------------------------
         if (MAX(RMSAVQ,RMSAVM).gt.RMSAV0) then
            write(6,*) 'ITERATION DIVERGED ---'
            ICONT = 0
            go to 260
         end if
         !----------------------------------------------------------------------
         if (ITSCF.ge.SCFSTEPS) then
            t_inc%i_iteration = t_inc%N_iteration
            write(6,'(12X,A)') '++++++ NUMBER OF SCF STEPS EXHAUSTED ++++++'
            write(6,'(79("*"))')
            ICONT = 0
            goto 260
         end if
      end if
      !-------------------------------------------------------------------------
      260  continue                  ! jump mark
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! ! ! ! ! ! ! ! ! ! ! ! ! !    ITERATION END    ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !
      !-------------------------------------------------------------------------
      ! Update energy contour
      !-------------------------------------------------------------------------
      if ( ICONT.eq.1 ) then
         call EPATHTB(EZ,DEZ,EMAX,IELAST,IESEMICORE,IDOSEMICORE,EMIN,EMAX,TK, &
            NPOL,NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,                   &
            NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IEMXD)
         do IE = 1,IELAST
            WEZ(IE) = -2.D0/PI*DEZ(IE)
            if ( IE.LE.IESEMICORE ) WEZ(IE) = WEZ(IE)*FSEMICORE
         end do
         write(1337,'(79("="))')
      end if
      !-------------------------------------------------------------------------
      ! Convert VISP potential to the relativistic form VTREL,BTREL.
      !-------------------------------------------------------------------------
      if ( KREL.eq.1 ) then
         call RELPOTCVT(2,VISP,ZAT,RMESH,DRDI,IRCUT,VTREL,BTREL,ZREL, &
            RMREL,JWSREL,DRDIREL,R2DRDIREL,IRSHIFT,IPAND,IRMD,NPOTD,NATYP)
      endif
      !
      !-------------------------------------------------------------------------
      ! Write out information for the next iteration
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! New_energy_mesh
      !-------------------------------------------------------------------------
      call save_emesh(IELAST,EZ,WEZ,EMIN,EMAX,IESEMICORE,FSEMICORE,NPOL,TK,   &
         NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,                           &
         NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IEMXD,t_params)
      !-------------------------------------------------------------------------
      ! Output_potential
      !-------------------------------------------------------------------------
      call save_scfinfo(t_params,VINS,VISP,ECORE,VBC,RMREL,DRDIREL,  &
         R2DRDIREL,ZREL,JWSREL,IRSHIFT,VTREL,BTREL,                  &
         ITSCF,SCFSTEPS,EFOLD,CHRGOLD,CMOMHOST,KREL,                 &
         IRMIND,IRMD,LMPOT,NSPOTD,NATYP,NPOTD,                     &
         NEMBD1)
      !-------------------------------------------------------------------------
      9020 format ('                old',' E Fermi ',F14.10,' Delta E_F = ',E16.8)
      9030 format ('                new',' E FERMI ',F14.10,'  DOS(E_F) = ',F12.6)
      9100 format (19X,' TIME IN ITERATION : ',f9.2,/,79('*'))
      9160 format(20X,'mixing factor used : ',1P,D12.2)
      1080 format('CMOMC',2I6)
      1090 format(4D22.14)

      ! Deallocate arrays
      i_all=-product(shape(vxcm))*kind(vxcm)
      deallocate(vxcm,stat=i_stat)
      call memocc(i_stat,i_all,'vxcm','main2')
      i_all=-product(shape(vxcnm))*kind(vxcnm)
      deallocate(vxcnm,stat=i_stat)
      call memocc(i_stat,i_all,'vxcnm','main2')
      i_all=-product(shape(thetas))*kind(thetas)
      deallocate(R2NEF,stat=i_stat)
      call memocc(i_stat,i_all,'R2NEF','main2')
      i_all=-product(shape(RHO2NS))*kind(RHO2NS)
      deallocate(RHO2NS,stat=i_stat)
      call memocc(i_stat,i_all,'RHO2NS','main2')
      i_all=-product(shape(RHO2NSNM))*kind(RHO2NSNM)
      deallocate(RHO2NSNM,stat=i_stat)
      call memocc(i_stat,i_all,'RHO2NSNM','main2')
      deallocate(VONS,stat=i_stat)
      i_all = product(shape(VONS))*kind(VONS)
      call memocc(i_stat,i_all,'VONS','main2')

   end subroutine main2


   !----------------------------------------------------------------------------
   ! SUBROUTINE: POTENSHIFT
   !> @brief Adds a constant (=VSHIFT) to the potentials of atoms
   !----------------------------------------------------------------------------
   subroutine POTENSHIFT(VISP,VINS,NATYP,NSPIN, IRCUT,IRC,IRMIN,NTCELL,IMAXSH,   &
      ILM_MAP,IFUNM,LMSP,LMPOT,GSH,THETAS,THESME,RFPI,RMESH,KSHAPE,VSHIFT,IRMD,NPOTD, &
      IRMIND,LMXSPD)

      implicit none
      !
      ! .. Input variables
      integer, intent(in) :: IRMD       ! > Maximum number of radial points
      integer, intent(in) :: LMPOT     ! > (LPOT+1)**2
      integer, intent(in) :: NATYP     ! > Number of kinds of atoms in unit cell
      integer, intent(in) :: NSPIN     ! > Counter for spin directions
      integer, intent(in) :: NPOTD     ! > (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      integer, intent(in) :: LMXSPD    ! > (2*LPOT+1)**2
      integer, intent(in) :: IRMIND    ! > IRMD-IRNSD
      integer, intent(in) :: KSHAPE    ! > Exact treatment of WS cell
      real (kind=dp), intent(in) :: RFPI
      real (kind=dp), intent(in) :: VSHIFT
      integer, dimension(NATYP), intent(in)           :: IRC      ! > R point for potential cutting
      integer, dimension(NATYP), intent(in)           :: IRMIN    ! > Max R for spherical treatment
      integer, dimension(NATYP), intent(in)           :: NTCELL   ! > Index for WS cell
      integer, dimension(0:LMPOT), intent(in)         :: IMAXSH
      integer, dimension(NGSHD,3), intent(in)         :: ILM_MAP
      integer, dimension(NATYP,LMXSPD), intent(in)    :: LMSP
      integer, dimension(NATYP,LMXSPD), intent(in)    :: IFUNM
      integer, dimension(0:IPAND,NATYP), intent(in)   :: IRCUT    ! > R points of panel borders
      real (kind=dp), dimension(NGSHD), intent(in)              :: GSH
      real (kind=dp), dimension(IRMD,NATYP), intent(in)          :: RMESH    ! > Radial mesh ( in units a Bohr)
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(in)  :: THESME
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(in)  :: THETAS   ! > shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      ! .. Input/Output:
      real (kind=dp), dimension(IRMD,NPOTD), intent(inout) :: VISP  ! > Spherical part of the potential
      real (kind=dp), dimension(IRMIND:IRMD,LMPOT,NSPOTD), intent(inout) :: VINS   ! > Non-spherical part of the potential
      ! .. Local variables
      integer :: ISPIN,IH,IPOT,IR,LM,IMT1,IRC1,IRMIN1
      real (kind=dp), dimension(IRMD) :: PSHIFTR
      real (kind=dp), dimension(IRMD,LMPOT) :: PSHIFTLMR

      do IH = 1,NATYP
         IMT1 = IRCUT(1,IH)
         IRC1 = IRC(IH)
         IRMIN1 = IRMIN(IH)
         do ISPIN = 1,NSPIN
            write (1337,*) 'SHIFTING OF THE POTENTIALS OF ATOM',IH,' BY', VSHIFT, 'RY.'
            IPOT = NSPIN * (IH-1) + ISPIN
            !
            call RINIT(IRMD*LMPOT,PSHIFTLMR)
            call RINIT(IRMD,PSHIFTR)
            do IR = 1,IRC1
               PSHIFTLMR(IR,1) = VSHIFT
            enddo
            !
            if (KSHAPE.EQ.0) then ! ASA
               do IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               end do
            else                ! Full-potential
               call CONVOL(IMT1,IRC1,NTCELL(IH),IMAXSH(LMPOT),ILM_MAP,IFUNM,   &
                  LMPOT,GSH,THETAS,THESME,0.0_dp,RFPI,RMESH(1,IH),PSHIFTLMR, &
                  PSHIFTR,LMSP)
               !
               do IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               enddo
               !
               do LM = 2,LMPOT
                  do IR = IRMIN1,IRC1
                     VINS(IR,LM,IPOT)=VINS(IR,LM,IPOT)+PSHIFTLMR(IR,LM)*RFPI
                  enddo
               enddo
            end if              ! (kshape.eq.0)
         end do
      end do

   end subroutine potenshift

end module MOD_MAIN2
