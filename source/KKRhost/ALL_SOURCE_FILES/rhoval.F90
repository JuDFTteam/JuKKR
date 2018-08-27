module mod_rhoval

contains

!-------------------------------------------------------------------------------
! SUBROUTINE: RHOVAL
!
!> @note - Average WLDAU for spherical wavefunctions:
!> The spherical part of the d or f wavefunction is found by adding
!> the average interaction potential WLDAUAV to the spherical
!> potential. Then the non-spherical parts are found by using only
!> the deviation of WLDAU from the average. This speeds up the
!> convergence of the Born series. See also subroutines
!> regsol, pnstmat and pnsqns
!>
!> @note -LDA+U implementation Mar. 2002-Dec.2004 Ph. Mavropoulos, H. Ebert, V. Popescu
!> @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine RHOVAL(IHOST,LDORHOEF,ICST,INS,IELAST,NSRA,ISPIN,NSPIN,NSPINPOT,I1,EZ,&
   WEZ,DRDI,R,VINS,VISP,ZAT,IPAN,IRCUT,IRMIN,THETAS,IFUNM,LMSP,RHO2NS,R2NEF,     &
   RHOORB,DEN,DENLM,MUORB,ESPV,CLEB,LOFLM,ICLEB,IEND,JEND,SOLVER,SOCTL,CTL,VTREL,&
   BTREL,RMREL,DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,ITERMVDIR,QMTET,QMPHI,      &
   MVEVIL,MVEVILEF,NMVECMAX,IDOLDAU,LOPT,PHILDAU,WLDAU,DENMATC,NATYP,NQDOS,LMAX)
   !
#ifdef CPP_MPI
   use mpi
#endif
   use mod_mympi, only: myrank, master
#ifdef CPP_MPI
   use mod_mympi, only: distribute_linear_on_tasks
#endif
   use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid, init_tgmat
   use Constants
   use mod_Profiling
   use mod_version_info
   use global_variables
   use mod_DataTypes
   use mod_pnsqns
   use mod_cradwf
   use mod_rholm
   use mod_rhons
   use mod_wfmesh
   use mod_drvrho

   implicit none

   ! .. Input variables
   integer, intent(in) :: I1
   integer, intent(in) :: INS
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
   integer, intent(in) :: IPAN      !< Number of panels in non-MT-region
   integer, intent(in) :: ICST      !< Number of Born approximation
   integer, intent(in) :: NSRA
   integer, intent(in) :: ZREL      !< atomic number (cast integer)
   integer, intent(in) :: LOPT      !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
   integer, intent(in) :: ISPIN
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: IHOST
   integer, intent(in) :: IRMIN     !< Max R for spherical treatment
   integer, intent(in) :: IELAST
   integer, intent(in) :: JWSREL    !< index of the WS radius
   integer, intent(in) :: IRSHIFT   !< shift of the REL radial mesh with respect no NREL
   integer, intent(in) :: IDOLDAU   !< flag to perform LDA+U
   integer, intent(in) :: NSPINPOT
   integer, intent(in) :: NMVECMAX
   integer, intent(inout) :: NQDOS
   real (kind=dp), intent(in) :: ZAT !< Nuclear charge
   logical, intent(in) :: LDORHOEF
   character(len=10), intent(in) :: SOLVER
   real (kind=dp), dimension(IRMD), intent(in)          :: R
   real (kind=dp), dimension(KREL*LMAX+1), intent(in)  :: CTL
   real (kind=dp), dimension(IRMD), intent(in)          :: DRDI  !< Derivative dr/di
   real (kind=dp), dimension(IRMD), intent(in)          :: VISP  !< Spherical part of the potential
   real (kind=dp), dimension(IRMD*KREL+(1-KREL)), intent(in) :: VTREL       !< potential (spherical part)
   real (kind=dp), dimension(IRMD*KREL+(1-KREL)), intent(in) :: BTREL       !< magnetic field
   real (kind=dp), dimension(IRMD*KREL+(1-KREL)), intent(in) :: RMREL       !< radial mesh
   real (kind=dp), dimension(KREL*LMAX+1), intent(in)       :: SOCTL
   real (kind=dp), dimension(IRMD*KREL+(1-KREL)), intent(in) :: DRDIREL     !< derivative of radial mesh
   real (kind=dp), dimension(IRMD*KREL+(1-KREL)), intent(in) :: R2DRDIREL   !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
   real (kind=dp), dimension(IRMIND:IRMD,LMPOTD), intent(in) :: VINS        !< Non-spherical part of the potential
   real (kind=dp), dimension(NCLEB,2), intent(in)           :: CLEB        !< GAUNT coefficients (GAUNT)
   real (kind=dp), dimension(IRID,NFUND), intent(in)        :: THETAS      !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   complex (kind=dp), dimension(IEMXD), intent(in) :: EZ
   complex (kind=dp), dimension(IEMXD), intent(in) :: WEZ
   complex (kind=dp), dimension(IRMD), intent(in)  :: PHILDAU
   complex (kind=dp), dimension(0:LMAX+1,IEMXD*(1+KREL),NQDOS), intent(in) :: DEN
   complex (kind=dp), dimension(LMMAXD,IEMXD*(1+KREL),NQDOS), intent(in)   :: DENLM
   ! .. In/Out variables
   real (kind=dp), dimension(MMAXD,MMAXD,NSPIND), intent(inout) :: WLDAU !< potential matrix
   !---------------------------------------------------------------------------
   !     IHOST = 1   < -- this routine is called by the HOST tbkkr-program
   !     IHOST <> 1  < --                 called by the IMPURITY program
   !---------------------------------------------------------------------------
   ! .. Output variables
   real (kind=dp), dimension(IRMD*KREL+(1-KREL)), intent(out) :: RHOORB
   real (kind=dp), dimension(0:LMAX+1+1,3), intent(out)      :: MUORB    !< orbital magnetic moment
   real (kind=dp), dimension(0:LMAX+1,2), intent(out)        :: ESPV     !< changed for REL case
   real (kind=dp), dimension(IRMD,LMPOTD,2), intent(out)       :: R2NEF    !< rho at FERMI energy
   real (kind=dp), dimension(IRMD,LMPOTD,2), intent(out)       :: RHO2NS   !< radial density
   complex (kind=dp), dimension(MMAXD,MMAXD), intent(out) :: DENMATC
   !----------------------------------------------------------------------------
   !      ITERMDIR variables
   !----------------------------------------------------------------------------
   logical, intent(in) :: ITERMVDIR
   real (kind=dp), intent(in) :: QMTET  !< \f$ \theta\f$ angle of the agnetization with respect to the z-axis
   real (kind=dp), intent(in) :: QMPHI  !< \f$ \phi\f$ angle of the agnetization with respect to the z-axis
   complex (kind=dp), dimension(0:LMAX,3,NMVECMAX), intent(out) :: MVEVIL ! OUTPUT
   complex (kind=dp), dimension(0:LMAX,3,NMVECMAX), intent(out) :: MVEVILEF ! OUTPUT
   !----------------------------------------------------------------------------
   !      ITERMDIR variables
   !----------------------------------------------------------------------------
   integer, dimension(LMXSPD), intent(in)    :: LMSP
   integer, dimension(LMXSPD), intent(in)    :: IFUNM
   integer, dimension(0:IPAND), intent(in)   :: IRCUT    !< R points of panel borders
   integer, dimension(LM2D), intent(in)      :: LOFLM    !< l of lm=(l,m) (GAUNT)
   integer, dimension(NCLEB,4), intent(in)   :: ICLEB    !< Pointer array
   integer, dimension(LMPOTD,0:LMAX,0:LMAX), intent(in) :: JEND !< Pointer array for icleb()
   ! .. Local Scalars
   ! .. Parameters
   integer :: LMAXD1
   integer :: i_stat, i_all
   real (kind=dp) :: WLDAUAV
   complex (kind=dp) :: DF,ERYD,EK
#ifndef CPP_MPI
   complex (kind=dp) :: DENTOT ! qdos
#endif
   integer :: IDIM,IE,IR,L,LM1,LM2,LMHI,LMLO,IREC,ISPINPOT,LASTEZ,M1,MMAX
   integer :: IQ  ! NQDOS ! qdos number of qdos points
   integer :: IX  ! qdos
   integer :: LRECGFLLE ! lmlm-dos
   integer, dimension(4) :: LMSHIFT1 ! lmlm-dos
   integer, dimension(4) :: LMSHIFT2 ! lmlm-dos

   ! .. Local Arrays
   real (kind=dp), dimension(0:LMAX)    :: S
   real (kind=dp), dimension(IRMD)       :: CUTOFF
   real (kind=dp), dimension(IRMD,0:LMAX) :: RS
   complex (kind=dp), dimension(0:LMAX)   :: EKL
   complex (kind=dp), dimension(0:LMAX)   :: TMAT
   complex (kind=dp), dimension(0:LMAX)   :: ALPHA
   complex (kind=dp), dimension(0:LMAX+1) :: DENDUM
   complex (kind=dp), dimension(LMMAXD)   :: DUM_DENLM
   complex (kind=dp), dimension(IRMD,0:LMAX)     :: QZ
   complex (kind=dp), dimension(IRMD,0:LMAX)     :: SZ
   complex (kind=dp), dimension(IRMD,0:LMAX)     :: PZ
   complex (kind=dp), dimension(IRMD,0:LMAX)     :: FZ
   complex (kind=dp), dimension(LMMAXD,LMMAXD)  :: AR
   complex (kind=dp), dimension(LMMAXD,LMMAXD)  :: CR
   complex (kind=dp), dimension(LMMAXD,LMMAXD)  :: DR
   complex (kind=dp), dimension(LMMAXD,LMMAXD)  :: GMAT0
   complex (kind=dp), dimension(LMMAXD,LMMAXD,IEMXD) :: GMATLL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,IRMIND:IRMD,2) :: PNS
   complex (kind=dp), dimension(LMMAXD,LMMAXD,IRMIND:IRMD,2) :: QNS
   !     .. first 2 indices in dmuorb are the spin-resolved contributions,
   !     .. the 3rd one should be the sum of them
   complex (kind=dp), dimension(0:KREL*LMAX+(1-KREL),3) :: DMUORB
   ! .. Local allocatable arrays
   complex (kind=dp), dimension(:,:), allocatable :: QVEC   !< qdos, q-vectors for qdos
   complex (kind=dp), dimension(:,:), allocatable :: GLDAU
   complex (kind=dp), dimension(:,:), allocatable :: DUM_GFLLE ! lmlm-dos
   complex (kind=dp), dimension(:,:,:,:), allocatable :: GFLLE ! qdos

   ! This routine needs irregular wavefunctions
   logical :: LIRRSOL
   parameter ( LIRRSOL = .TRUE. )

#ifdef CPP_MPI
   integer :: ie_start
#endif
   integer :: ie_end, ie_num
   ! .. External Functions
   logical :: OPT,TEST                          ! qdos
   external :: OPT,TEST                         ! qdos

   LMAXD1= LMAX+1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LDAU
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ( IDOLDAU.EQ.1 ) then
      WLDAUAV = 0.D0
      LMLO = LOPT*LOPT + 1
      LMHI = (LOPT+1)*(LOPT+1)
      MMAX = LMHI - LMLO + 1
      do M1 = 1,MMAX
         WLDAUAV = WLDAUAV + WLDAU(M1,M1,ISPIN)
      enddo
      WLDAUAV = WLDAUAV/DBLE(MMAX)
      !-------------------------------------------------------------------------
      ! Note: Application if WLDAU makes the potential discontinuous.
      ! A cutoff can be used if necessary to make the potential continuous
      ! for example (array bounds should be adjusted):
      !
      ! CUTOFF(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) * &
      !  ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
      ! CUTOFF(IR) = 1D0/CUTOFF(IR)
      !-------------------------------------------------------------------------
      do M1 = 1,IRMD
         CUTOFF(M1) = 1.D0
      end do
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LDAU
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Initialise variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ( KREL.EQ.0 ) then
      do LM1 = 1,LMPOTD
         do IR = 1,IRMD
            RHO2NS(IR,LM1,ISPIN) = 0.0D0
         end do
      end do
      !
      do L = 0,LMAXD1
         ESPV(L,ISPIN) = 0.0D0
      end do
      !-------------------------------------------------------------------------
   else
      !-------------------------------------------------------------------------
      do ISPINPOT = 1,2
         do LM1 = 1,LMPOTD
            do IR = 1,IRMD
               RHO2NS(IR,LM1,ISPINPOT) = 0.0D0
               R2NEF(IR,LM1,ISPINPOT)  = 0.0D0
            end do
         end do
         !
         do L = 0,LMAXD1
            ESPV(L,ISPINPOT) = 0.0D0
         end do
         !
      end do
      !
      do IR = 1,IRMD
         RHOORB(IR) = 0.0D0
      end do
      !
      do IR = 1,3
         do L = 0,LMAX
            DMUORB(L,IR) = CZERO
         end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ITERMDIR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ITERMVDIR) then
         do LM1 = 1,3
            do LM2 = 1, NMVECMAX
               do L = 0, LMAX
                  MVEVIL(L,LM1,LM2) = CZERO
                  MVEVILEF(L,LM1,LM2) = CZERO
               end do
            end do
         end do
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ITERMDIR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end if ! KREL = 0/1
   !----------------------------------------------------------------------------
   LASTEZ = IELAST
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! End initialise variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef CPP_MPI
   if (OPT('qdos    ')) then                                                     ! qdos
      if(NATYP.ge.100) then                                                      ! qdos
         open(31,FILE="qdos."//char(48+I1/100)//char(48+mod(I1/10,10))//   &     ! qdos
            char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")                    ! qdos
      else                                                                       ! qdos
         open(31,FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//  &     ! qdos
            char(48+ISPIN)//".dat")                                              ! qdos
      end if                                                                     ! qdos
      call version_print_header(31)                                              ! qdos
      write (31,'(7(A,3X))') '#   Re(E)','Im(E)',&                               ! qdos
         'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'                               ! qdos
   endif                                                                         ! qdos
   !
   open(30,FILE="lmdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//       &     ! lmdos
      char(48+ISPIN)//".dat")                                                    ! lmdos
   call version_print_header(31)                                                 ! lmdos
   write (30,*) ' '                                                              ! lmdos
   write (30,8600) '# ISPIN=',ISPIN,' I1=',I1                                    ! lmdos
   8600 format (a8,I3,a4,I5)                                                     ! lmdos

   ! write out complex qdos for interpolation to the real axis                   ! complex qdos
   if(test('compqdos')) then                                                     ! complex qdos
      if(NATYP.ge.100) then                                                      ! complex qdos
         open(3031,FILE="cqdos."//char(48+I1/100)//char(48+mod(I1/10,10))//   &  ! complex qdos
            char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")                    ! complex qdos
      else                                                                       ! complex qdos
         open(3031,FILE="cqdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//  &  ! complex qdos
            char(48+ISPIN)//".dat")                                              ! complex qdos
      end if                                                                     ! complex qdos
      call version_print_header(3031)                                            ! complex qdos
      write (3031,'(A)') '# lmax, natyp, nspin, nqdos, ielast:'                  ! complex qdos
      write (3031,'(5I9)') lmax, natyp, nspin, nqdos, ielast                    ! complex qdos
      write (3031,'(7(A,3X))') '#   Re(E)','Im(E)',&                             ! complex qdos
         'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'                               ! complex qdos
   endif                                                                         ! qdos
#endif
   !
   ! set LMSHIFT value which is need to construct dentmp
   LMSHIFT1(1)=0
   LMSHIFT1(2)=LMMAXD
   LMSHIFT1(3)=0
   LMSHIFT1(4)=LMMAXD
   LMSHIFT2(1)=0
   LMSHIFT2(2)=LMMAXD
   LMSHIFT2(3)=LMMAXD
   LMSHIFT2(4)=0
   NQDOS = 1                                                                     ! qdos
   if (OPT('qdos    ')) then                                                     ! qdos
      ! Read BZ path for qdos calculation:                                       ! qdos
      open(67,FILE='qvec.dat')                                                   ! qdos
      read(67,*) NQDOS                                                           ! qdos
      allocate(QVEC(3,NQDOS),stat=i_stat)                                        ! qdos
      call memocc(i_stat,product(shape(QVEC))*kind(QVEC),'QVEC','RHOVAL')        ! qdos
      do IQ = 1,NQDOS                                                            ! qdos
         read(67,*) (QVEC(IX,IQ),IX=1,3)                                         ! qdos
      enddo                                                                      ! qdos
      close(67)                                                                  ! qdos
   end if
   !
   allocate(GFLLE(LMMAXD,LMMAXD,IELAST,NQDOS),stat=i_stat)
   call memocc(i_stat,product(shape(GFLLE))*kind(GFLLE),'GFLLE','RHOVAL')
   allocate(DUM_GFLLE(LMMAXD,LMMAXD),stat=i_stat)
   call memocc(i_stat,product(shape(DUM_GFLLE))*kind(DUM_GFLLE),'DUM_GFLLE','RHOVAL')
   if (IDOLDAU.EQ.1) then
      allocate(GLDAU(LMMAXD,LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(GLDAU))*kind(GLDAU),'GLDAU','RHOVAL')
      GLDAU=CZERO
   endif
   !
   if(myrank==master) write(1337,*) 'atom',I1
#ifdef CPP_MPI
   ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
   ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

   do ie_num=1,ie_end
      IE = ie_start+ie_num
#else
! ! omp: start parallel region here
! !$omp parallel do default(none)
! !$omp& private(eryd,ie,ir,irec,lm1,lm2)
! !$omp& private(jlk_index,tmatll,ith)
! !$omp& shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmax)
! !$omp& shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew,i1)
! !$omp& shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop)
! !$omp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)
! !$omp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,cdentemp)
! !$omp& shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft)
! !$omp& private(iq,df,ek,tmattemp,gmatll,gmat0,iorb,dentemp)
! !$omp& private(rho2ns_temp,dentot,dentmp,rho2,temp1)
! !$omp& shared(ldorhoef,nqdos,lmshift1,lmshift2,wez,lmsp,imt1,ifunm)
! !$omp& shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop)
! !$omp& reduction(+:rho2int,espv) reduction(-:muorb)
! !$omp& reduction(-:denorbmom,denorbmomsp,denorbmomlm,denorbmomns)
! !$omp& shared(t_tgmat)
! !$omp& private(alphasph,alphall)
   do IE=1,IELAST
         ie_num = ie
         ie_end = ielast
#endif
      if(t_inc%i_write>0) write(1337,*) 'energy',IE,EZ(IE)
      !
      ERYD = EZ(IE)
      DF = WEZ(IE)/DBLE(NSPINPOT)
      !-------------------------------------------------------------------------
      ! non/scalar-relativistic OR relativistic
      !-------------------------------------------------------------------------
      if ( KREL.EQ.0 ) then
         call WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,R,S,RS,IRCUT(IPAN),IRMD,LMAX)
         call CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S,PZ,FZ,QZ,SZ,&
            TMAT,VISP,DRDI,R,ZAT,LIRRSOL,IDOLDAU,LOPT,WLDAUAV,CUTOFF)
         !-----------------------------------------------------------------------
         ! Non-spherical
         !-----------------------------------------------------------------------
         if (INS.GT.0) then
            call PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ,  &
               PNS,QNS,NSRA,VINS,IPAN,IRMIN,IRCUT,          & ! Added IRMIN 1.7.2014
               CLEB,ICLEB,IEND,LOFLM,LMAX,                  &
               IDOLDAU,LOPT,LMLO,LMHI,                      &
               WLDAU(1,1,ISPIN),WLDAUAV,CUTOFF, LMAX)
         endif
         !
         do L = 0,LMAX
            EKL(L) = EK*DBLE(2*L+1)
         end do
         !
         do IQ = 1,NQDOS                                                         ! qdos
            !-------------------------------------------------------------------
            ! Read in Green function
            !-------------------------------------------------------------------
            IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST * (ISPIN-1) + &          ! qdos (without qdos, IQ=NQDOS=1)
               NQDOS * IELAST * NSPIN * (I1-1)                                   ! qdos
            if (t_tgmat%gmat_to_file) then
               read(69,REC=IREC) GMAT0
            else
               IREC = IQ + NQDOS * (ie_num-1) + NQDOS *  &
                  t_mpi_c_grid%ntot2 * (ISPIN-1) +       &
                  NQDOS * t_mpi_c_grid%ntot2 * NSPIN * (I1-1)
               GMAT0(:,:) = t_tgmat%gmat(:,:,irec)
            end if
            if (TEST('GMAT=0  ')) then
               write(*,*) 'TEST GMAT=0, setting GMAT to zero'
               GMAT0 =CZERO
            endif
            !-------------------------------------------------------------------
            ! Spherical/non-spherical input potential
            !-------------------------------------------------------------------
            if ( INS.EQ.0 ) then
               call RHOLM(DEN(0,IE,IQ),DF,GMAT0,NSRA,             &
                  RHO2NS(1,1,ISPIN),DRDI,IPAN,IRCUT,PZ,FZ,QZ,SZ,  &
                  CLEB(1,1),ICLEB,IEND,JEND,EKL)
            else
               call RHONS(DEN(0,IE,IQ),DF,DRDI,GMAT0,EK,                &
                  RHO2NS(1,1,ISPIN),IPAN,IRCUT,IRMIN,THETAS,IFUNM,LMSP, & ! Added IRMIN 1.7.2014
                  NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,       &
                  JEND,IEND,EKL,DENLM(1,IE,IQ),GFLLE(:,:,IE,IQ))
            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! LDA+U
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (IDOLDAU.EQ.1) then
               do LM1=1,LMMAXD
                  do LM2=1,LMMAXD
                     GLDAU(LM1,LM2)=GLDAU(LM1,LM2)+ DF*GFLLE(LM1,LM2,IE,IQ)
                  enddo
               enddo
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! LDA+U
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef CPP_MPI
            ! Write out qdos:
            if (OPT('qdos    ')) then                                            ! qdos
               DENTOT = CMPLX(0.D0,0.D0, kind=dp)                                        ! qdos
               do L = 0,LMAXD1                                                   ! qdos
                  DENTOT = DENTOT + DEN(L,IE,IQ)                                 ! qdos
               enddo                                                             ! qdos
               write(30,9000) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),       &     ! lmdos
                  -aimag(DENTOT)/PI,(-aimag(DENLM(L,IE,IQ))/PI,L=1,LMMAXD)       ! lmdos
               write(31,9000) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),       &     ! qdos
                  -aimag(DENTOT)/PI,(-aimag(DENLM(L,IE,IQ))/PI,L=1,LMMAXD)       ! qdos
               9000 format(5F10.6,40E16.8)                                       ! qdos
               ! writeout complex qdos for interpolation                         ! complex qdos
               if(test('compqdos')) then                                         ! complex qdos
                  write(3031,9001) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),  &     ! complex qdos
                     DENTOT,(DENLM(L,IE,IQ),L=1,LMMAXD)                          ! complex qdos
               end if                                                            ! complex qdos
               9001 format(6F10.6,80E16.8)                                       ! qdos
            endif                                                                ! qdos
#endif
         end do ! IQ                                                             ! qdos
         !----------------------------------------------------------------------
         do L = 0,LMAXD1
            ESPV(L,ISPIN) = ESPV(L,ISPIN)+aimag(ERYD*DEN(L,IE,1)*DF)
         end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! get the charge at the Fermi energy (IELAST)
         ! call RHOLM/RHONS with the energy weight CONE --> not overwrite DF
         !                  with the dummy DENDUM       --> not overwrite DEN
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if ( (IE.EQ.IELAST) .AND. LDORHOEF ) then
            if (INS.EQ.0) then
               call RHOLM(DENDUM,CONE,GMAT0,NSRA,              &
                  R2NEF(1,1,ISPIN),DRDI,IPAN,IRCUT,PZ,FZ,QZ,SZ,&
                  CLEB(1,1),ICLEB,IEND,JEND,EKL)
            else
               call RHONS(DENDUM,CONE,DRDI,GMAT0,EK,                    &
                  R2NEF(1,1,ISPIN),IPAN,IRCUT,IRMIN,THETAS,IFUNM,LMSP,  & ! Added IRMIN 1.7.2014
                  NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,       &
                  JEND,IEND,EKL,DUM_DENLM,DUM_GFLLE)
            end if
         end if
         !----------------------------------------------------------------------
      else ! ( KREL.EQ.0 )
         !----------------------------------------------------------------------
         IQ = 1 ! reset IQ to zero, problem with qdos!!!

         !              !!!! PROBLEM WITH ARRAY DIMENSIONS FOR VTREL ETC. !!!!
         ! #ifdef CPP_MPI
         !              call MPI_FINALIZE(L)
         ! #endif
         !              stop '[rhoval] ERROR array dimensions need to be checked!'
         call DRVRHO_QDOS(LDORHOEF,RHO2NS,R2NEF,DEN,DMUORB,RHOORB,&
            IE,ERYD,DF,LASTEZ,                                    &
            GMATLL,VTREL,BTREL,RMREL,DRDIREL,                     &
            R2DRDIREL,ZREL,JWSREL,IRSHIFT,SOLVER,SOCTL,CTL,       &
            QMTET,QMPHI,ITERMVDIR,MVEVIL,MVEVILEF,LMMAXD,         &
            LMAX,IRMD,LMPOTD,IEMXD,NMVECMAX,                        &
            I1,NQDOS)                                                            ! qdos
         !
         do L = 0,LMAXD1
            ESPV(L,1) = ESPV(L,1) + aimag(ERYD*DEN(L,IE,IQ)*DF)
            ESPV(L,2) = ESPV(L,2) + aimag(ERYD*DEN(L,IE+IEMXD,IQ)*DF)
         end do
         !
         do IR = 1,3
            do L = 0,LMAX
               MUORB(L,IR) = MUORB(L,IR) + aimag(DMUORB(L,IR)*DF)
            end do
         end do
      end if
      !-------------------------------------------------------------------------
      ! Non/scalar-relativistic OR relativistic
      !-------------------------------------------------------------------------
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDAU
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDA+U calculation
      !         IF ( ( IDOLDAU.EQ.1 ).AND.( LOPT.GE.0 ) )
      !     &        CALL DENSITYMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL(1,1,IE),
      !     &                        IPAN,IRCUT,DRDI,EK,
      !     &                        IRMIN,LOPT,MMAX,LMLO,LMHI,PHILDAU,DENMATC
      !     &        ,den,ie) ! test fivos 19.9.08
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDAU
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end do ! IE = 1,IELAST

   ! LDA+U
   if (IDOLDAU.EQ.1) then
      do M1=1,MMAX
         LM1=LMLO-1+M1
         DENMATC(1:MMAX,M1)=(1.0/(2.0*CI))*&
            (GLDAU(LMLO:LMHI,LM1)-CONJG(GLDAU(LM1,LMLO:LMHI)))
      enddo
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write out gflle
   if (OPT('lmlm-dos')) then                                                     ! lmlm-dos
      if (ISPIN.EQ.1) then                                                       ! lmlm-dos
         ! 4 words = 16 bytes / complex number                                   ! lmlm-dos
         LRECGFLLE = WLENGTH*4*LMMAXD*LMMAXD*IELAST*NQDOS                        ! lmlm-dos
         open(96,ACCESS='direct',RECL=LRECGFLLE,FILE="gflle",FORM='unformatted') ! lmlm-dos
      endif                                                                      ! lmlm-dos
      IREC = I1 + NATYP * (ISPIN-1)                                             ! lmlm-dos
      write(96,REC=IREC) GFLLE(:,:,:,:)                                          ! lmlm-dos
   endif
   !
   i_all=-product(shape(GFLLE))*kind(GFLLE)
   deallocate(GFLLE,stat=i_stat)
   call memocc(i_stat,i_all,'GFLLE','RHOVAL')
   i_all=-product(shape(DUM_GFLLE))*kind(DUM_GFLLE)
   deallocate(DUM_GFLLE,stat=i_stat)
   call memocc(i_stat,i_all,'DUM_GFLLE','RHOVAL')
   !
   if (IDOLDAU.EQ.1) then
      i_all=-product(shape(GLDAU))*kind(GLDAU)
      deallocate(GLDAU,stat=i_stat)
      call memocc(i_stat,i_all,'GLDAU','RHOVAL')
   endif
   !
   if ( IHOST.NE.1) return
   !
   ! Transformation of ISPIN=1,2 from (spin-down,spin-up) to (charge-density,spin-density)
   if (ISPIN.EQ.2) then
      IDIM = IRMD*LMPOTD
      call DSCAL(IDIM,2.D0,RHO2NS(1,1,1),1)
      call DAXPY(IDIM,-0.5D0,RHO2NS(1,1,1),1,RHO2NS(1,1,2),1)
      call DAXPY(IDIM,1.0D0,RHO2NS(1,1,2),1,RHO2NS(1,1,1),1)
      !-------------------------------------------------------------------------
      ! Do the same at the Fermi energy
      !-------------------------------------------------------------------------
      call DSCAL(IDIM,2.D0,R2NEF(1,1,1),1)
      call DAXPY(IDIM,-0.5D0,R2NEF(1,1,1),1,R2NEF(1,1,2),1)
      call DAXPY(IDIM,1.0D0,R2NEF(1,1,2),1,R2NEF(1,1,1),1)
   end if

end subroutine RHOVAL

end module mod_rhoval
