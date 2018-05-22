!-------------------------------------------------------------------------------
! MODULE: memoryhandling
!
! DESCRIPTION:
!> @brief Subroutine to handle allocation/deallocation of arrays
!> @details Module to handle the allocation of arrays to be later distributed
!> it aims to bring modularity to the memory management
!
!> @author
!> Jonathan Chico
!> @date 14.11.2017
!> @todo The number of arrays in the misc section should be reduced, and they should
!> be located in the appropriate routines
!-------------------------------------------------------------------------------
module memoryhandling

   use Profiling

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_cell
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the unit cell.
   !
   !> @author
   !> Jonathan Chico
   !> @date 14.11.2017
   !----------------------------------------------------------------------------
   subroutine allocate_cell(flag,NAEZ,NEMB,NATYP,CLS,IMT,IRWS,IRNS,NTCELL,REFPOT,&
      KFG,KAOEZ,RMT,ZAT,RWS,MTFAC,RMTREF,RMTREFAT,RMTNEW,RBASIS,LMXC)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: NEMB !< number of 'embedding' positions
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      integer, dimension(:), allocatable, intent(inout) :: CLS !< Cluster around atomic sites
      integer, dimension(:), allocatable, intent(inout) :: IMT !< R point at MT radius
      integer, dimension(:), allocatable, intent(inout) :: IRWS !< R point at WS radius
      integer, dimension(:), allocatable, intent(inout) :: IRNS !< Position of atoms in the unit cell in units of bravais vectors
      integer, dimension(:), allocatable, intent(inout) :: LMXC
      integer, dimension(:), allocatable, intent(inout) :: NTCELL !< Index for WS cell
      integer, dimension(:), allocatable, intent(inout) :: REFPOT !< Ref. pot. card  at position
      integer, dimension(:,:), allocatable, intent(inout) :: KFG
      integer, dimension(:,:), allocatable, intent(inout) :: KAOEZ !< atom types located at a given site
      double precision, dimension(:), allocatable, intent(inout) :: RMT !< Muffin-tin radius of true system
      double precision, dimension(:), allocatable, intent(inout) :: ZAT !< Nuclear charge
      double precision, dimension(:), allocatable, intent(inout) :: RWS !< Wigner Seitz radius
      double precision, dimension(:), allocatable, intent(inout) :: MTFAC  !< Scaling factor for radius MT
      double precision, dimension(:), allocatable, intent(inout) :: RMTREF !< Muffin-tin radius of reference system
      double precision, dimension(:), allocatable, intent(inout) :: RMTREFAT
      double precision, dimension(:), allocatable, intent(inout) :: RMTNEW !< Adapted muffin-tin radius
      double precision, dimension(:,:), allocatable, intent(inout) :: RBASIS !< Position of atoms in the unit cell in units of bravais vectors

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(REFPOT(NAEZ+NEMB),stat=i_stat)
         call memocc(i_stat,product(shape(REFPOT))*kind(REFPOT),'REFPOT','allocate_cell')
         REFPOT = 0
         allocate(RMTREFAT(NAEZ+NEMB),stat=i_stat)
         call memocc(i_stat,product(shape(RMTREFAT))*kind(RMTREFAT),'RMTREFAT','allocate_cell')
         RMTREFAT = -1.D0 ! Signals the need for later calculation
         allocate(KAOEZ(NATYP,NAEZ+NEMB),stat=i_stat)
         call memocc(i_stat,product(shape(KAOEZ))*kind(KAOEZ),'KAOEZ','allocate_cell')
         KAOEZ = 0
         allocate(CLS(NAEZ+NEMB),stat=i_stat)
         call memocc(i_stat,product(shape(CLS))*kind(CLS),'CLS','allocate_cell')
         CLS = 1
         allocate(RBASIS(3,NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(RBASIS))*kind(RBASIS),'RBASIS','allocate_cell')
         RBASIS = 0.D0
         allocate(MTFAC(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(MTFAC))*kind(MTFAC),'MTFAC','allocate_cell')
         MTFAC = 0.D0
         allocate(ZAT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(ZAT))*kind(ZAT),'ZAT','allocate_cell')
         ZAT = -1.D0      ! Negative value signals read-in from pot-file
         allocate(NTCELL(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(NTCELL))*kind(NTCELL),'NTCELL','allocate_cell')
         NTCELL = 0
         allocate(RMTREF(NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(RMTREF))*kind(RMTREF),'RMTREF','allocate_cell')
         RMTREF = -1.D0
         allocate(IRNS(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IRNS))*kind(IRNS),'IRNS','allocate_cell')
         IRNS = -1        ! Negative value signals to use FPRADIUS
         allocate(RWS(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(RWS))*kind(RWS),'RWS','allocate_cell')
         RWS = 0.D0
         allocate(IRWS(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IRWS))*kind(IRWS),'IRWS','allocate_cell')
         IRWS = 0
         allocate(RMT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(RMT))*kind(RMT),'RMT','allocate_cell')
         RMT = 0.D0
         allocate(RMTNEW(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(RMTNEW))*kind(RMTNEW),'RMTNEW','allocate_cell')
         RMTNEW = 0.D0
         allocate(IMT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IMT))*kind(IMT),'IMT','allocate_cell')
         IMT = 0
         allocate(KFG(4,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(KFG))*kind(KFG),'KFG','allocate_cell')
         KFG = 0
         allocate(LMXC(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(LMXC))*kind(LMXC),'LMXC','allocate_cell')
         LMXC = 0
      else
         if (allocated(ZAT)) then
            i_all=-product(shape(ZAT))*kind(ZAT)
            deallocate(ZAT,stat=i_stat)
            call memocc(i_stat,i_all,'ZAT','allocate_cell')
         endif
         if (allocated(NTCELL)) then
            i_all=-product(shape(NTCELL))*kind(NTCELL)
            deallocate(NTCELL,stat=i_stat)
            call memocc(i_stat,i_all,'NTCELL','allocate_cell')
         endif
         if (allocated(NTCELL)) then
            i_all=-product(shape(NTCELL))*kind(NTCELL)
            deallocate(NTCELL,stat=i_stat)
            call memocc(i_stat,i_all,'NTCELL','allocate_cell')
         endif
         if (allocated(RMTREF)) then
            i_all=-product(shape(RMTREF))*kind(RMTREF)
            deallocate(RMTREF,stat=i_stat)
            call memocc(i_stat,i_all,'RMTREF','allocate_cell')
         endif
         if (allocated(IRNS)) then
            i_all=-product(shape(IRNS))*kind(IRNS)
            deallocate(IRNS,stat=i_stat)
            call memocc(i_stat,i_all,'IRNS','allocate_cell')
         endif
         if (allocated(IRWS)) then
            i_all=-product(shape(IRWS))*kind(IRWS)
            deallocate(IRWS,stat=i_stat)
            call memocc(i_stat,i_all,'IRWS','allocate_cell')
         endif
         if (allocated(RWS)) then
            i_all=-product(shape(RWS))*kind(RWS)
            deallocate(RWS,stat=i_stat)
            call memocc(i_stat,i_all,'RWS','allocate_cell')
         endif
         if (allocated(RMT)) then
            i_all=-product(shape(RMT))*kind(RMT)
            deallocate(RMT,stat=i_stat)
            call memocc(i_stat,i_all,'RMT','allocate_cell')
         endif
         if (allocated(RMTNEW)) then
            i_all=-product(shape(RMTNEW))*kind(RMTNEW)
            deallocate(RMTNEW,stat=i_stat)
            call memocc(i_stat,i_all,'RMTNEW','allocate_cell')
         endif
         if (allocated(IMT)) then
            i_all=-product(shape(IMT))*kind(IMT)
            deallocate(IMT,stat=i_stat)
            call memocc(i_stat,i_all,'IMT','allocate_cell')
         endif
         if (allocated(KFG)) then
            i_all=-product(shape(KFG))*kind(KFG)
            deallocate(KFG,stat=i_stat)
            call memocc(i_stat,i_all,'KFG','allocate_cell')
         endif
         if(allocated(KAOEZ)) then
            i_all=-product(shape(KAOEZ))*kind(KAOEZ)
            deallocate(KAOEZ,stat=i_stat)
            call memocc(i_stat,i_all,'KAOEZ','allocate_cell')
         endif
         if (allocated(REFPOT)) then
            i_all=-product(shape(REFPOT))*kind(REFPOT)
            deallocate(REFPOT,stat=i_stat)
            call memocc(i_stat,i_all,'REFPOT','allocate_cell')
         endif
         if (allocated(RMTREFAT)) then
            i_all=-product(shape(RMTREFAT))*kind(RMTREFAT)
            deallocate(RMTREFAT,stat=i_stat)
            call memocc(i_stat,i_all,'RMTREFAT','allocate_cell')
         endif
         if (allocated(LMXC)) then
            i_all=-product(shape(LMXC))*kind(LMXC)
            deallocate(LMXC,stat=i_stat)
            call memocc(i_stat,i_all,'LMXC','allocate_cell')
         endif
      endif

   end subroutine allocate_cell

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_semi_inf_host
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the left and right host for the calculation of slabs
   !
   !> @author
   !> Jonathan Chico
   !> @date 20.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_semi_inf_host(flag,NEMB,TLEFT,TRIGHT)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NEMB !< number of 'embedding' positions
      double precision, dimension(:,:), allocatable, intent(inout) :: TLEFT  !< Vectors of the basis for the left host
      double precision, dimension(:,:), allocatable, intent(inout) :: TRIGHT !< vectors of the basis for the right host

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         allocate(TLEFT(3,NEMB+1),stat=i_stat)
         call memocc(i_stat,product(shape(TLEFT))*kind(TLEFT),'TLEFT','allocate_cell')
         TLEFT = 0.D0
         allocate(TRIGHT(3,NEMB+1),stat=i_stat)
         call memocc(i_stat,product(shape(TRIGHT))*kind(TRIGHT),'TRIGHT','allocate_cell')
         TRIGHT = 0.D0
      else
         if (allocated(TLEFT)) then
            i_all=-product(shape(TLEFT))*kind(TLEFT)
            deallocate(TLEFT,stat=i_stat)
            call memocc(i_stat,i_all,'TLEFT','allocate_cell')
         endif
         if (allocated(TRIGHT)) then
            i_all=-product(shape(TRIGHT))*kind(TRIGHT)
            deallocate(TRIGHT,stat=i_stat)
            call memocc(i_stat,i_all,'TRIGHT','allocate_cell')
         endif
      endif


   end subroutine allocate_semi_inf_host

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_potential
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the potential
   !
   !> @author
   !> Jonathan Chico
   !> @date 14.11.2017
   !----------------------------------------------------------------------------
   subroutine allocate_potential(flag,NAEZ,NEMB,IRM,NATYP,NPOTD,IPAND,NFUND,LMXSPD, &
      LMPOT,IRMIND,NSPOTD,NFU,IRC,NCORE,IRMIN,LMSP,LMSP1,IRCUT,LCORE,LLMSP,ITITLE,  &
      FPRADIUS,VISP,ECORE,VINS)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: NEMB !< number of 'embedding' positions
      integer, intent(in) :: IRM
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      integer, intent(in) :: NPOTD !< 2*NATYP
      integer, intent(in) :: IPAND
      integer, intent(in) :: NFUND
      integer, intent(in) :: LMXSPD
      integer, intent(in) :: LMPOT
      integer, intent(in) :: IRMIND
      integer, intent(in) :: NSPOTD
      integer, dimension(:), allocatable, intent(inout) :: NFU
      integer, dimension(:), allocatable, intent(inout) :: IRC !< R point for potential cutting
      integer, dimension(:), allocatable, intent(inout) :: NCORE !< Number of core states
      integer, dimension(:), allocatable, intent(inout) :: IRMIN !< Max R for spherical treatment
      integer, dimension(:,:), allocatable, intent(inout) :: LMSP !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      integer, dimension(:,:), allocatable, intent(inout) :: LMSP1
      integer, dimension(:,:), allocatable, intent(inout) :: IRCUT !< R points of panel borders
      integer, dimension(:,:), allocatable, intent(inout) :: LCORE !< Angular momentum of core states
      integer, dimension(:,:), allocatable, intent(inout) :: LLMSP !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
      integer, dimension(:,:), allocatable, intent(inout) :: ITITLE
      double precision, dimension(:), allocatable, intent(inout) :: FPRADIUS !< R point at which full-potential treatment starts
      double precision, dimension(:,:), allocatable, intent(inout) :: VISP   !< Spherical part of the potential
      double precision, dimension(:,:), allocatable, intent(inout) :: ECORE  !< Core energies
      double precision, dimension(:,:,:), allocatable, intent(inout) :: VINS !< Non-spherical part of the potential

      ! .. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         allocate(FPRADIUS(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(FPRADIUS))*kind(FPRADIUS),'FPRADIUS','allocate_potential')
         FPRADIUS = -1.D0 ! Negative value signals to use IRNS from pot-file (sub. startb1)
         allocate(IRC(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IRC))*kind(IRC),'IRC','allocate_potential')
         IRC = 0
         allocate(IRCUT(0:IPAND,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IRCUT))*kind(IRCUT),'IRCUT','allocate_potential')
         IRCUT = 0
         allocate(IRMIN(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IRMIN))*kind(IRMIN),'IRMIN','allocate_potential')
         IRMIN = 0
         allocate(LCORE(20,NPOTD),stat=i_stat)
         call memocc(i_stat,product(shape(LCORE))*kind(LCORE),'LCORE','allocate_potential')
         LCORE = 0
         allocate(ECORE(20,NPOTD),stat=i_stat)
         call memocc(i_stat,product(shape(ECORE))*kind(ECORE),'ECORE','allocate_potential')
         ECORE = 0.D0
         allocate(NCORE(NPOTD),stat=i_stat)
         call memocc(i_stat,product(shape(NCORE))*kind(NCORE),'NCORE','allocate_potential')
         NCORE = 0
         allocate(LMSP1(LMXSPD,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(LMSP1))*kind(LMSP1),'LMSP1','allocate_potential')
         LMSP1 = 0
         allocate(LLMSP(NATYP,NFUND),stat=i_stat)
         call memocc(i_stat,product(shape(LLMSP))*kind(LLMSP),'LLMSP','allocate_potential')
         LLMSP = 0
         allocate(LMSP(NATYP,LMXSPD),stat=i_stat)
         call memocc(i_stat,product(shape(LMSP))*kind(LMSP),'LMSP','allocate_potential')
         LMSP = 0
         allocate(ITITLE(20,NPOTD),stat=i_stat)
         call memocc(i_stat,product(shape(ITITLE))*kind(ITITLE),'ITITLE','allocate_potential')
         ITITLE = 0
         allocate(NFU(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(NFU))*kind(NFU),'NFU','allocate_potential')
         NFU = 0
         allocate(VINS(IRMIND:IRM,LMPOT,NSPOTD), stat=i_stat)
         call memocc(i_stat,product(shape(VINS))*kind(VINS),'VINS','allocate_misc')
         VINS = 0.D0
         allocate(VISP(IRM,NPOTD),stat=i_stat)
         call memocc(i_stat,product(shape(VISP))*kind(VISP),'VISP','allocate_misc')
         VISP = 0.D0

      else
         if (allocated(FPRADIUS)) then
            i_all=-product(shape(FPRADIUS))*kind(FPRADIUS)
            deallocate(FPRADIUS,stat=i_stat)
            call memocc(i_stat,i_all,'FPRADIUS','allocate_potential')
         endif
         if (allocated(IRC)) then
            i_all=-product(shape(IRC))*kind(IRC)
            deallocate(IRC,stat=i_stat)
            call memocc(i_stat,i_all,'IRC','allocate_potential')
         endif
         if (allocated(IRCUT)) then
            i_all=-product(shape(IRCUT))*kind(IRCUT)
            deallocate(IRCUT,stat=i_stat)
            call memocc(i_stat,i_all,'IRCUT','allocate_potential')
         endif
         if (allocated(IRMIN)) then
            i_all=-product(shape(IRMIN))*kind(IRMIN)
            deallocate(IRMIN,stat=i_stat)
            call memocc(i_stat,i_all,'IRMIN','allocate_potential')
         endif
         if (allocated(LCORE)) then
            i_all=-product(shape(LCORE))*kind(LCORE)
            deallocate(LCORE,stat=i_stat)
            call memocc(i_stat,i_all,'LCORE','allocate_potential')
         endif
         if (allocated(LMSP1)) then
            i_all=-product(shape(LMSP1))*kind(LMSP1)
            deallocate(LMSP1,stat=i_stat)
            call memocc(i_stat,i_all,'LMSP1','allocate_potential')
         endif
         if (allocated(LLMSP)) then
            i_all=-product(shape(LLMSP))*kind(LLMSP)
            deallocate(LLMSP,stat=i_stat)
            call memocc(i_stat,i_all,'LLMSP','allocate_potential')
         endif
         if (allocated(LMSP)) then
            i_all=-product(shape(LMSP))*kind(LMSP)
            deallocate(LMSP,stat=i_stat)
            call memocc(i_stat,i_all,'LMSP','allocate_potential')
         endif
         if (allocated(ITITLE)) then
            i_all=-product(shape(ITITLE))*kind(ITITLE)
            deallocate(ITITLE,stat=i_stat)
            call memocc(i_stat,i_all,'ITITLE','allocate_potential')
         endif
         if (allocated(NFU)) then
            i_all=-product(shape(NFU))*kind(NFU)
            deallocate(NFU,stat=i_stat)
            call memocc(i_stat,i_all,'NFU','allocate_potential')
         endif
         if (allocated(VINS)) then
            i_all=-product(shape(VINS))*kind(VINS)
            deallocate(VINS,stat=i_stat)
            call memocc(i_stat,i_all,'VINS','allocate_misc')
         endif
         if (allocated(VISP)) then
            i_all=-product(shape(VISP))*kind(VISP)
            deallocate(VISP,stat=i_stat)
            call memocc(i_stat,i_all,'VISP','allocate_misc')
         endif

      endif

   end subroutine allocate_potential

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_cpa
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the CPA treatment
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_cpa(flag,NAEZ,NEMB,NATYP,NOQ,ICPA,IQAT,HOSTIMP,CONC)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: NEMB !< number of 'embedding' positions
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell

      integer, dimension(:), allocatable, intent(inout) :: NOQ !< Number of diff. atom types located
      integer, dimension(:), allocatable, intent(inout) :: ICPA !< ICPA = 0/1 site-dependent CPA flag
      integer, dimension(:), allocatable, intent(inout) :: IQAT !< the site on which an atom is located on a given site
      integer, dimension(:), allocatable, intent(inout) :: HOSTIMP
      double precision, dimension(:), allocatable, intent(inout) :: CONC !< concentration of a given atom

      ! .. Local variables
      integer :: i_stat, i_all


      if (flag>0) then

         allocate(NOQ(NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(NOQ))*kind(NOQ),'NOQ','allocate_cpa')
         NOQ = 1
         allocate(ICPA(NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(ICPA))*kind(ICPA),'ICPA','allocate_cpa')
         ICPA = 0
         allocate(IQAT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IQAT))*kind(IQAT),'IQAT','allocate_cpa')
         IQAT = 0
         allocate(CONC(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(CONC))*kind(CONC),'CONC','allocate_cpa')
         CONC = 1.D0
         allocate(HOSTIMP(0:NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(HOSTIMP))*kind(HOSTIMP),'HOSTIMP','allocate_cpa')
         HOSTIMP = 0

      else

         if (allocated(NOQ)) then
            i_all=-product(shape(NOQ))*kind(NOQ)
            deallocate(NOQ,stat=i_stat)
            call memocc(i_stat,i_all,'NOQ','allocate_cpa')
         endif
         if (allocated(ICPA)) then
            i_all=-product(shape(ICPA))*kind(ICPA)
            deallocate(ICPA,stat=i_stat)
            call memocc(i_stat,i_all,'ICPA','allocate_cpa')
         endif
         if (allocated(IQAT)) then
            i_all=-product(shape(IQAT))*kind(IQAT)
            deallocate(IQAT,stat=i_stat)
            call memocc(i_stat,i_all,'IQAT','allocate_cpa')
         endif
         if(allocated(CONC)) then
            i_all=-product(shape(CONC))*kind(CONC)
            deallocate(CONC,stat=i_stat)
            call memocc(i_stat,i_all,'CONC','allocate_cpa')
         endif
         if(allocated(HOSTIMP)) then
            i_all=-product(shape(HOSTIMP))*kind(HOSTIMP)
            deallocate(HOSTIMP,stat=i_stat)
            call memocc(i_stat,i_all,'HOSTIMP','allocate_cpa')
         endif

      endif

   end subroutine allocate_cpa

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_ldau
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the LDA+U approach
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_ldau(flag,NATYP,LOPT,UEFF,JEFF,EREFLDAU)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell

      integer, dimension(:), allocatable, intent(inout) :: LOPT !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      double precision, dimension(:), allocatable, intent(inout) :: UEFF !< input U parameter for each atom
      double precision, dimension(:), allocatable, intent(inout) :: JEFF !< input J parameter for each atom
      double precision, dimension(:), allocatable, intent(inout) :: EREFLDAU !< the energies of the projector's wave functions (REAL)

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(LOPT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(LOPT))*kind(LOPT),'LOPT','allocate_ldau')
         LOPT = -1        !  not perform lda+u (default)
         allocate(UEFF(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(UEFF))*kind(UEFF),'UEFF','allocate_ldau')
         UEFF = 0.D0
         allocate(JEFF(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(JEFF))*kind(JEFF),'JEFF','allocate_ldau')
         JEFF = 0.D0
         allocate(EREFLDAU(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(EREFLDAU))*kind(EREFLDAU),'EREFLDAU','allocate_ldau')
         EREFLDAU = 0.5D0
      else
         if (allocated(LOPT)) then
            i_all=-product(shape(LOPT))*kind(LOPT)
            deallocate(LOPT,stat=i_stat)
            call memocc(i_stat,i_all,'LOPT','allocate_ldau')
         endif
         if (allocated(UEFF)) then
            i_all=-product(shape(UEFF))*kind(UEFF)
            deallocate(UEFF,stat=i_stat)
            call memocc(i_stat,i_all,'UEFF','allocate_ldau')
         endif
         if (allocated(JEFF)) then
            i_all=-product(shape(JEFF))*kind(JEFF)
            deallocate(JEFF,stat=i_stat)
            call memocc(i_stat,i_all,'JEFF','allocate_ldau')
         endif
         if (allocated(EREFLDAU)) then
            i_all=-product(shape(EREFLDAU))*kind(EREFLDAU)
            deallocate(EREFLDAU,stat=i_stat)
            call memocc(i_stat,i_all,'EREFLDAU','allocate_ldau')
         endif

      endif

   end subroutine allocate_ldau

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_ldau_potential
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the potentials for the LDA+U approach
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_ldau_potential(flag,IRM,NATYP,MMAXD,NSPIND,ITLDAU,WLDAU,&
      ULDAU,PHILDAU)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: IRM
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      integer, intent(in) :: MMAXD
      integer, intent(in) :: NSPIND !< Counter for spin directions (KREL+(1-KREL)*(KSP+1))
      integer, dimension(:), allocatable, intent(inout) :: ITLDAU !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
      double precision, dimension(:,:,:,:), allocatable, intent(inout) :: WLDAU !< potential matrix
      double precision, dimension(:,:,:,:,:), allocatable, intent(inout) :: ULDAU !< calculated Coulomb matrix elements (EREFLDAU)
      double complex, dimension(:,:), allocatable, intent(inout) :: PHILDAU

      ! .. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         allocate(ITLDAU(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(ITLDAU))*kind(ITLDAU),'ITLDAU','allocate_ldau_potential')
         ITLDAU = 0
         allocate(ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(ULDAU))*kind(ULDAU),'ULDAU','allocate_ldau_potential')
         ULDAU = 0.D0
         allocate(WLDAU(MMAXD,MMAXD,NSPIND,NATYP) ,stat=i_stat)
         call memocc(i_stat,product(shape(WLDAU))*kind(WLDAU),'WLDAU','allocate_ldau_potential')
         WLDAU =  0.D0
         allocate(PHILDAU(IRM,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(PHILDAU))*kind(PHILDAU),'PHILDAU','allocate_ldau_potential')
         PHILDAU=(0.D0,0.D0)

      else
         if (allocated(ITLDAU)) then
            i_all=-product(shape(ITLDAU))*kind(ITLDAU)
            deallocate(ITLDAU,stat=i_stat)
            call memocc(i_stat,i_all,'ITLDAU','allocate_ldau_potential')
         endif
         if (allocated(ULDAU)) then
            i_all=-product(shape(ULDAU))*kind(ULDAU)
            deallocate(ULDAU,stat=i_stat)
            call memocc(i_stat,i_all,'ULDAU','allocate_ldau_potential')
         endif
         if (allocated(WLDAU)) then
            i_all=-product(shape(WLDAU))*kind(WLDAU)
            deallocate(WLDAU,stat=i_stat)
            call memocc(i_stat,i_all,'WLDAU','allocate_ldau_potential')
         endif
         if (allocated(PHILDAU)) then
            i_all=-product(shape(PHILDAU))*kind(PHILDAU)
            deallocate(PHILDAU,stat=i_stat)
            call memocc(i_stat,i_all,'PHILDAU','allocate_ldau_potential')
         endif

      endif

   end subroutine allocate_ldau_potential

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_magnetization
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the magnetisation
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_magnetization(flag,NAEZ,NATYP,LMMAXD,INIPOL,IXIPOL,QMTET,&
      QMPHI,DROTQ)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      integer, intent(in) :: LMMAXD
      integer, dimension(:), allocatable, intent(inout) :: INIPOL !< Initial spin polarisation
      integer, dimension(:), allocatable, intent(inout) :: IXIPOL !< Constraint of spin pol.
      double precision, dimension(:), allocatable, intent(inout) :: QMTET
      double precision, dimension(:), allocatable, intent(inout) :: QMPHI
      double complex, dimension(:,:,:), allocatable, intent(inout) :: DROTQ !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity

      !.. Local variables
      integer :: i_stat, i_all

      if(flag>0) then
         allocate(QMTET(NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(QMTET))*kind(QMTET),'QMTET','allocate_magnetization')
         QMTET = 0.D0
         allocate(QMPHI(NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(QMPHI))*kind(QMPHI),'QMPHI','allocate_magnetization')
         QMPHI = 0.D0
         allocate(INIPOL(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(INIPOL))*kind(INIPOL),'INIPOL','allocate_magnetization')
         INIPOL = 0
         allocate(IXIPOL(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IXIPOL))*kind(IXIPOL),'IXIPOL','allocate_magnetization')
         IXIPOL = 0
         allocate(DROTQ(LMMAXD,LMMAXD,NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(DROTQ))*kind(DROTQ),'DROTQ','allocate_magnetization')
         DROTQ = (0.D0,0.D0)
      else
         if (allocated(QMTET)) then
            i_all=-product(shape(QMTET))*kind(QMTET)
            deallocate(QMTET,stat=i_stat)
            call memocc(i_stat,i_all,'QMTET','allocate_magnetization')
         endif
         if (allocated(QMPHI)) then
            i_all=-product(shape(QMPHI))*kind(QMPHI)
            deallocate(QMPHI,stat=i_stat)
            call memocc(i_stat,i_all,'QMPHI','allocate_magnetization')
         endif
         if (allocated(INIPOL)) then
            i_all=-product(shape(INIPOL))*kind(INIPOL)
            deallocate(INIPOL,stat=i_stat)
            call memocc(i_stat,i_all,'INIPOL','allocate_magnetization')
         endif
         if (allocated(IXIPOL)) then
            i_all=-product(shape(IXIPOL))*kind(IXIPOL)
            deallocate(IXIPOL,stat=i_stat)
            call memocc(i_stat,i_all,'IXIPOL','allocate_magnetization')
         endif
         if (allocated(DROTQ)) then
            i_all=-product(shape(DROTQ))*kind(DROTQ)
            deallocate(DROTQ,stat=i_stat)
            call memocc(i_stat,i_all,'DROTQ','allocate_magnetization')
         endif

      endif

   end subroutine allocate_magnetization

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_SOC
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the spin-orbit coupling (SOC)
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_SOC(flag,KREL,NATYP,LMAX,SOCSCALE,CSCL,SOCSCL)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: KREL
      integer, intent(in) :: LMAX !< Maximum l component in wave function expansion
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      double precision, dimension(:), allocatable, intent(inout) :: SOCSCALE !< Spin-orbit scaling
      double precision, dimension(:,:), allocatable, intent(inout) :: CSCL !< Speed of light scaling
      double precision, dimension(:,:), allocatable, intent(inout) :: SOCSCL

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(SOCSCL(KREL*LMAX+1,KREL*NATYP+(1-KREL)),stat=i_stat)
         call memocc(i_stat,product(shape(SOCSCL))*kind(SOCSCL),'SOCSCL','allocate_SOC')
         SOCSCL = 1.D0
         allocate(CSCL(KREL*LMAX+1,KREL*NATYP+(1-KREL)),stat=i_stat)
         call memocc(i_stat,product(shape(CSCL))*kind(CSCL),'CSCL','allocate_SOC')
         CSCL = 0.D0
         allocate(SOCSCALE(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(SOCSCALE))*kind(SOCSCALE),'SOCSCALE','allocate_SOC')
         SOCSCALE = 1.D0  ! Spin-orbit scaling
      else
         if (allocated(SOCSCL)) then
            i_all=-product(shape(SOCSCL))*kind(SOCSCL)
            deallocate(SOCSCL,stat=i_stat)
            call memocc(i_stat,i_all,'SOCSCL','allocate_SOC')
         endif
         if (allocated(SOCSCALE)) then
            i_all=-product(shape(SOCSCALE))*kind(SOCSCALE)
            deallocate(SOCSCALE,stat=i_stat)
            call memocc(i_stat,i_all,'SOCSCALE','allocate_SOC')
         endif
         if (allocated(CSCL)) then
            i_all=-product(shape(CSCL))*kind(CSCL)
            deallocate(CSCL,stat=i_stat)
            call memocc(i_stat,i_all,'CSCL','allocate_SOC')
         endif
      endif

   end subroutine allocate_SOC

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_energies
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe energies
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_energies(flag,IEMXD,EZ,DEZ,WEZ)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: IEMXD

      double complex, dimension(:), allocatable, intent(inout) :: EZ
      double complex, dimension(:), allocatable, intent(inout) :: DEZ
      double complex, dimension(:), allocatable, intent(inout) :: WEZ

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(EZ(IEMXD),stat=i_stat)
         call memocc(i_stat,product(shape(EZ))*kind(EZ),'EZ','allocate_energies')
         EZ = (0.D0,0.D0)
         allocate(DEZ(IEMXD),stat=i_stat)
         call memocc(i_stat,product(shape(DEZ))*kind(DEZ),'DEZ','allocate_energies')
         DEZ = (0.D0,0.D0)
         allocate(WEZ(IEMXD),stat=i_stat)
         call memocc(i_stat,product(shape(WEZ))*kind(WEZ),'WEZ','allocate_energies')
         WEZ = (0.D0,0.D0)
      else
         if (allocated(EZ)) then
            i_all=-product(shape(EZ))*kind(EZ)
            deallocate(EZ,stat=i_stat)
            call memocc(i_stat,i_all,'EZ','allocate_energies')
         endif
         if (allocated(DEZ)) then
            i_all=-product(shape(DEZ))*kind(DEZ)
            deallocate(DEZ,stat=i_stat)
            call memocc(i_stat,i_all,'DEZ','allocate_energies')
         endif
         if (allocated(WEZ)) then
            i_all=-product(shape(WEZ))*kind(WEZ)
            deallocate(WEZ,stat=i_stat)
            call memocc(i_stat,i_all,'WEZ','allocate_energies')
         endif

      endif

   end subroutine allocate_energies

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_relativistic
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe relativistic corrections
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_relativistic(flag,KREL,IRM,NAEZ,NATYP,ZREL,JWSREL,IRSHIFT,&
      VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,QMGAM,QMGAMTAB,QMPHITAB,QMTETTAB)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: KREL
      integer, intent(in) :: IRM
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      integer, dimension(:), allocatable, intent(inout) :: ZREL !< atomic number (cast integer)
      integer, dimension(:), allocatable, intent(inout) :: JWSREL !< index of the WS radius
      integer, dimension(:), allocatable, intent(inout) :: IRSHIFT !< shift of the REL radial mesh with respect no NREL
      double precision, dimension(:), allocatable, intent(inout) :: QMGAM
      double precision, dimension(:,:), allocatable, intent(inout) :: VTREL !< potential (spherical part)
      double precision, dimension(:,:), allocatable, intent(inout) :: BTREL !< magnetic field
      double precision, dimension(:,:), allocatable, intent(inout) :: RMREL !< radial mesh
      double precision, dimension(:,:), allocatable, intent(inout) :: DRDIREL !< derivative of radial mesh
      double precision, dimension(:,:), allocatable, intent(inout) :: R2DRDIREL !< r**2 * drdi
      double precision, dimension(:,:), allocatable, intent(inout) :: QMGAMTAB
      double precision, dimension(:,:), allocatable, intent(inout) :: QMPHITAB
      double precision, dimension(:,:), allocatable, intent(inout) :: QMTETTAB

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(VTREL(IRM*KREL+(1-KREL),NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(VTREL))*kind(VTREL),'VTREL','allocate_relativistic')
         VTREL = 0.D0
         allocate(BTREL(IRM*KREL+(1-KREL),NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(BTREL))*kind(BTREL),'BTREL','allocate_relativistic')
         BTREL= 0.D0
         allocate(DRDIREL(IRM*KREL+(1-KREL),NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(DRDIREL))*kind(DRDIREL),'DRDIREL','allocate_relativistic')
         DRDIREL = 0.D0
         allocate(R2DRDIREL(IRM*KREL+(1-KREL),NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(R2DRDIREL))*kind(R2DRDIREL),'R2DRDIREL','allocate_relativistic')
         R2DRDIREL =0.D0
         allocate(RMREL(IRM*KREL+(1-KREL),NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(RMREL))*kind(RMREL),'RMREL','allocate_relativistic')
         RMREL = 0.D0
         allocate(IRSHIFT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IRSHIFT))*kind(IRSHIFT),'IRSHIFT','allocate_relativistic')
         IRSHIFT = 0
         allocate(JWSREL(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(JWSREL))*kind(JWSREL),'JWSREL','allocate_relativistic')
         JWSREL = 0
         allocate(ZREL(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(ZREL))*kind(ZREL),'ZREL','allocate_relativistic')
         ZREL = 0
         allocate(QMGAM(NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(QMGAM))*kind(QMGAM),'QMGAM','allocate_relativistic')
         QMGAM = 0.D0
         allocate(QMGAMTAB(NAEZ,3),stat=i_stat)
         call memocc(i_stat,product(shape(QMGAMTAB))*kind(QMGAMTAB),'QMGAMTAB','allocate_relativistic')
         QMGAMTAB = 0.D0
         allocate(QMPHITAB(NAEZ,3),stat=i_stat)
         call memocc(i_stat,product(shape(QMPHITAB))*kind(QMPHITAB),'QMPHITAB','allocate_relativistic')
         QMPHITAB = 0.D0
         allocate(QMTETTAB(NAEZ,3),stat=i_stat)
         call memocc(i_stat,product(shape(QMTETTAB))*kind(QMTETTAB),'QMTETTAB','allocate_relativistic')
         QMTETTAB = 0.D0

      else
         if (allocated(VTREL)) then
            i_all=-product(shape(VTREL))*kind(VTREL)
            deallocate(VTREL,stat=i_stat)
            call memocc(i_stat,i_all,'VTREL','allocate_relativistic')
         endif
         if (allocated(BTREL)) then
            i_all=-product(shape(BTREL))*kind(BTREL)
            deallocate(BTREL,stat=i_stat)
            call memocc(i_stat,i_all,'BTREL','allocate_relativistic')
         endif
         if (allocated(DRDIREL)) then
            i_all=-product(shape(DRDIREL))*kind(DRDIREL)
            deallocate(DRDIREL,stat=i_stat)
            call memocc(i_stat,i_all,'DRDIREL','allocate_relativistic')
         endif
         if (allocated(R2DRDIREL)) then
            i_all=-product(shape(R2DRDIREL))*kind(R2DRDIREL)
            deallocate(R2DRDIREL,stat=i_stat)
            call memocc(i_stat,i_all,'R2DRDIREL','allocate_relativistic')
         endif
         if (allocated(RMREL)) then
            i_all=-product(shape(RMREL))*kind(RMREL)
            deallocate(RMREL,stat=i_stat)
            call memocc(i_stat,i_all,'RMREL','allocate_relativistic')
         endif
         if (allocated(IRSHIFT)) then
            i_all=-product(shape(IRSHIFT))*kind(IRSHIFT)
            deallocate(IRSHIFT,stat=i_stat)
            call memocc(i_stat,i_all,'IRSHIFT','allocate_relativistic')
         endif
         if (allocated(JWSREL)) then
            i_all=-product(shape(JWSREL))*kind(JWSREL)
            deallocate(JWSREL,stat=i_stat)
            call memocc(i_stat,i_all,'JWSREL','allocate_relativistic')
         endif
         if (allocated(ZREL)) then
            i_all=-product(shape(ZREL))*kind(ZREL)
            deallocate(ZREL,stat=i_stat)
            call memocc(i_stat,i_all,'ZREL','allocate_relativistic')
         endif
         if (allocated(QMGAM)) then
            i_all=-product(shape(QMGAM))*kind(QMGAM)
            deallocate(QMGAM,stat=i_stat)
            call memocc(i_stat,i_all,'QMGAM','allocate_relativistic')
         endif
         if (allocated(QMGAMTAB)) then
            i_all=-product(shape(QMGAMTAB))*kind(QMGAMTAB)
            deallocate(QMGAMTAB,stat=i_stat)
            call memocc(i_stat,i_all,'QMGAMTAB','allocate_relativistic')
         endif
         if (allocated(QMPHITAB)) then
            i_all=-product(shape(QMPHITAB))*kind(QMPHITAB)
            deallocate(QMPHITAB,stat=i_stat)
            call memocc(i_stat,i_all,'QMPHITAB','allocate_relativistic')
         endif
         if (allocated(QMTETTAB)) then
            i_all=-product(shape(QMTETTAB))*kind(QMTETTAB)
            deallocate(QMTETTAB,stat=i_stat)
            call memocc(i_stat,i_all,'QMTETTAB','allocate_relativistic')
         endif

      endif

   end subroutine allocate_relativistic

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_rel_transformations
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe relativistic transformations
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_rel_transformations(flag,LMMAXD,NRREL,IRREL,RC,CREL,RREL,SRREL)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: LMMAXD
      integer, dimension(:,:), allocatable, intent(inout) :: NRREL
      integer, dimension(:,:,:), allocatable, intent(inout) :: IRREL
      double complex, dimension(:,:), allocatable, intent(inout) :: RC !< NREL REAL spher. harm. >  CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
      double complex, dimension(:,:), allocatable, intent(inout) :: CREL !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
      double complex, dimension(:,:), allocatable, intent(inout) :: RREL !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
      double complex, dimension(:,:,:), allocatable,intent(inout) :: SRREL

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         allocate(RREL(LMMAXD,LMMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(RREL))*kind(RREL),'RREL','allocate_rel_transformations')
         RREL = (0.D0,0.D0)
         allocate(SRREL(2,2,LMMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(SRREL))*kind(SRREL),'SRREL','allocate_rel_transformations')
         SRREL = (0.D0,0.D0)
         allocate(IRREL(2,2,LMMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(IRREL))*kind(IRREL),'IRREL','allocate_rel_transformations')
         IRREL = 0
         allocate(NRREL(2,LMMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(NRREL))*kind(NRREL),'NRREL','allocate_rel_transformations')
         NRREL = 0
         allocate(CREL(LMMAXD,LMMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(CREL))*kind(CREL),'CREL','allocate_rel_transformations')
         CREL = (0.D0,0.D0)
         allocate(RC(LMMAXD,LMMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(RC))*kind(RC),'RC','allocate_rel_transformations')
         RC = (0.D0,0.D0)
      else
         if (allocated(RREL)) then
            i_all=-product(shape(RREL))*kind(RREL)
            deallocate(RREL,stat=i_stat)
            call memocc(i_stat,i_all,'RREL','allocate_rel_transformations')
         endif
         if (allocated(SRREL)) then
            i_all=-product(shape(SRREL))*kind(SRREL)
            deallocate(SRREL,stat=i_stat)
            call memocc(i_stat,i_all,'SRREL','allocate_rel_transformations')
         endif
         if (allocated(IRREL)) then
            i_all=-product(shape(IRREL))*kind(IRREL)
            deallocate(IRREL,stat=i_stat)
            call memocc(i_stat,i_all,'IRREL','allocate_rel_transformations')
         endif
         if (allocated(NRREL)) then
            i_all=-product(shape(NRREL))*kind(NRREL)
            deallocate(NRREL,stat=i_stat)
            call memocc(i_stat,i_all,'NRREL','allocate_rel_transformations')
         endif
         if (allocated(CREL)) then
            i_all=-product(shape(CREL))*kind(CREL)
            deallocate(CREL,stat=i_stat)
            call memocc(i_stat,i_all,'CREL','allocate_rel_transformations')
         endif
         if (allocated(RC)) then
            i_all=-product(shape(RC))*kind(RC)
            deallocate(RC,stat=i_stat)
            call memocc(i_stat,i_all,'RC','allocate_rel_transformations')
         endif

      endif

   end subroutine allocate_rel_transformations

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_clusters
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe clusters
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_clusters(flag,NAEZ,LMAX,NCLEB,NCLSD,NEMBD1,NSHELD,NACLSD,&
      LMPOT,NATOMIMPD,NSH1,NSH2,NACLS,NSHELL,ATOMIMP,ATOM,EZOA,ICLEB,JEND,RATOM,&
      RCLSIMP,CMOMHOST,RCLS)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: LMAX !< Maximum l component in wave function expansion
      integer, intent(in) :: NCLEB
      integer, intent(in) :: NCLSD
      integer, intent(in) :: NEMBD1
      integer, intent(in) :: NSHELD
      integer, intent(in) :: NACLSD
      integer, intent(in) :: LMPOT
      integer, intent(in) :: NATOMIMPD
      integer, dimension(:), allocatable, intent(inout) :: NSH1 !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(:), allocatable, intent(inout) :: NSH2 !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(:), allocatable, intent(inout) :: NACLS !< Number of atoms in cluster
      integer, dimension(:), allocatable, intent(inout) :: NSHELL !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
      integer, dimension(:), allocatable, intent(inout) :: ATOMIMP
      integer, dimension(:,:), allocatable, intent(inout) :: ATOM !< Atom at site in cluster
      integer, dimension(:,:), allocatable, intent(inout) :: EZOA !< EZ of atom at site in cluster
      integer, dimension(:,:), allocatable, intent(inout) :: ICLEB !< Pointer array
      integer, dimension(:,:,:), allocatable, intent(inout) :: JEND !< Pointer array for icleb()
      double precision, dimension(:,:), allocatable, intent(inout) :: RATOM
      double precision, dimension(:,:), allocatable, intent(inout) :: RCLSIMP
      double precision, dimension(:,:), allocatable, intent(inout) :: CMOMHOST !< Charge moments of each atom of the (left/right) host
      double precision, dimension(:,:,:), allocatable, intent(inout) :: RCLS !< Real space position of atom in cluster

      integer :: i_stat, i_all

      if (flag>0) then

         allocate(ATOM(3,NSHELD),stat=i_stat)
         call memocc(i_stat,product(shape(ATOM))*kind(ATOM),'ATOM','allocate_clusters')
         ATOM = 0
         allocate(RATOM(3,NSHELD),stat=i_stat)
         call memocc(i_stat,product(shape(RATOM))*kind(RATOM),'RATOM','allocate_clusters')
         RATOM = 0.D0
         allocate(RCLS(3,NACLSD,NCLSD),stat=i_stat)
         call memocc(i_stat,product(shape(RCLS))*kind(RCLS),'RCLS','allocate_clusters')
         RCLS = 0.D0
         allocate(RCLSIMP(3,NATOMIMPD),stat=i_stat)
         call memocc(i_stat,product(shape(RCLSIMP))*kind(RCLSIMP),'RCLSIMP','allocate_clusters')
         RCLSIMP = 0.D0
         allocate(NACLS(NCLSD),stat=i_stat)
         call memocc(i_stat,product(shape(NACLS))*kind(NACLS),'NACLS','allocate_clusters')
         NACLS = 0
         allocate(EZOA(NACLSD,NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(EZOA))*kind(EZOA),'EZOA','allocate_clusters')
         EZOA = 0
         allocate(ATOMIMP(NATOMIMPD),stat=i_stat)
         call memocc(i_stat,product(shape(ATOMIMP))*kind(ATOMIMP),'ATOMIMP','allocate_clusters')
         ATOMIMP = 0
         allocate(ICLEB(NCLEB,4),stat=i_stat)
         call memocc(i_stat,product(shape(ICLEB))*kind(ICLEB),'ICLEB','allocate_clusters')
         ICLEB = 0
         allocate(NSH1(NSHELD),stat=i_stat)
         call memocc(i_stat,product(shape(NSH1))*kind(NSH1),'NSH1','allocate_clusters')
         NSH1 = 0
         allocate(NSH2(NSHELD),stat=i_stat)
         call memocc(i_stat,product(shape(NSH2))*kind(NSH2),'NSH2','allocate_clusters')
         NSH2 = 0
         allocate(NSHELL(0:NSHELD),stat=i_stat)
         call memocc(i_stat,product(shape(NSHELL))*kind(NSHELL),'NSHELL','allocate_clusters')
         NSHELL = 0
         allocate(CMOMHOST(LMPOT,NEMBD1),stat=i_stat)
         call memocc(i_stat,product(shape(CMOMHOST))*kind(CMOMHOST),'CMOMHOST','allocate_clusters')
         CMOMHOST = 0.D0
         allocate(JEND(LMPOT,0:LMAX,0:LMAX),stat=i_stat)
         call memocc(i_stat,product(shape(JEND))*kind(JEND),'JEND','allocate_clusters')
         JEND = 0

      else
         if (allocated(ATOM)) then
            i_all=-product(shape(ATOM))*kind(ATOM)
            deallocate(ATOM,stat=i_stat)
            call memocc(i_stat,i_all,'ATOM','allocate_clusters')
         endif
         if (allocated(RATOM)) then
            i_all=-product(shape(RATOM))*kind(RATOM)
            deallocate(RATOM,stat=i_stat)
            call memocc(i_stat,i_all,'RATOM','allocate_clusters')
         endif
         if (allocated(RCLS)) then
            i_all=-product(shape(RCLS))*kind(RCLS)
            deallocate(RCLS,stat=i_stat)
            call memocc(i_stat,i_all,'RCLS','allocate_clusters')
         endif
         if (allocated(RCLSIMP)) then
            i_all=-product(shape(RCLSIMP))*kind(RCLSIMP)
            deallocate(RCLSIMP,stat=i_stat)
            call memocc(i_stat,i_all,'RCLSIMP','allocate_clusters')
         endif
         if (allocated(NACLS)) then
            i_all=-product(shape(NACLS))*kind(NACLS)
            deallocate(NACLS,stat=i_stat)
            call memocc(i_stat,i_all,'NACLS','allocate_clusters')
         endif
         if (allocated(EZOA)) then
            i_all=-product(shape(EZOA))*kind(EZOA)
            deallocate(EZOA,stat=i_stat)
            call memocc(i_stat,i_all,'EZOA','allocate_clusters')
         endif
         if (allocated(ATOMIMP)) then
            i_all=-product(shape(ATOMIMP))*kind(ATOMIMP)
            deallocate(ATOMIMP,stat=i_stat)
            call memocc(i_stat,i_all,'ATOMIMP','allocate_clusters')
         endif
         if (allocated(ICLEB)) then
            i_all=-product(shape(ICLEB))*kind(ICLEB)
            deallocate(ICLEB,stat=i_stat)
            call memocc(i_stat,i_all,'ICLEB','allocate_clusters')
         endif
         if (allocated(NSH1)) then
            i_all=-product(shape(NSH1))*kind(NSH1)
            deallocate(NSH1,stat=i_stat)
            call memocc(i_stat,i_all,'NSH1','allocate_clusters')
         endif
         if (allocated(NSH2)) then
            i_all=-product(shape(NSH2))*kind(NSH2)
            deallocate(NSH2,stat=i_stat)
            call memocc(i_stat,i_all,'NSH2','allocate_clusters')
         endif
         if (allocated(NSHELL)) then
            i_all=-product(shape(NSHELL))*kind(NSHELL)
            deallocate(NSHELL,stat=i_stat)
            call memocc(i_stat,i_all,'NSHELL','allocate_clusters')
         endif
         if (allocated(CMOMHOST)) then
            i_all=-product(shape(CMOMHOST))*kind(CMOMHOST)
            deallocate(CMOMHOST,stat=i_stat)
            call memocc(i_stat,i_all,'CMOMHOST','allocate_clusters')
         endif
         if (allocated(JEND)) then
            i_all=-product(shape(JEND))*kind(JEND)
            deallocate(JEND,stat=i_stat)
            call memocc(i_stat,i_all,'JEND','allocate_clusters')
         endif

      endif

   end subroutine allocate_clusters

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_expansion
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the functions for the expansion of the Green function
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_expansion(flag,LM2D,IRID,NFUND,NTOTD,NCLEB,LASSLD,NCELLD,&
      NCHEBD,LOFLM,WG,CLEB,YRG,THETAS,THETASNEW)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: LM2D
      integer, intent(in) :: IRID
      integer, intent(in) :: NFUND
      integer, intent(in) :: NTOTD
      integer, intent(in) :: NCLEB
      integer, intent(in) :: LASSLD
      integer, intent(in) :: NCELLD
      integer, intent(in) :: NCHEBD
      integer, dimension(:), allocatable, intent(inout) :: LOFLM !< l of lm=(l,m) (GAUNT)
      double precision, dimension(:), allocatable, intent(inout) :: WG !< Integr. weights for Legendre polynomials
      double precision, dimension(:,:), allocatable, intent(inout) :: CLEB !< GAUNT coefficients (GAUNT)
      double precision, dimension(:,:,:), allocatable, intent(inout) :: YRG !< Spherical harmonics (GAUNT2)
      double precision, dimension(:,:,:), allocatable, intent(inout) :: THETAS !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      double precision, dimension(:,:,:), allocatable, intent(inout) :: THETASNEW

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         allocate(WG(LASSLD),stat=i_stat)
         call memocc(i_stat,product(shape(WG))*kind(WG),'WG','allocate_expansion')
         WG = 0.D0
         allocate(YRG(LASSLD,0:LASSLD,0:LASSLD),stat=i_stat)
         call memocc(i_stat,product(shape(YRG))*kind(YRG),'YRG','allocate_expansion')
         YRG = 0.D0
         allocate(THETAS(IRID,NFUND,NCELLD),stat=i_stat)
         call memocc(i_stat,product(shape(THETAS))*kind(THETAS),'THETAS','allocate_expansion')
         THETAS = 0.D0
         allocate(THETASNEW(NTOTD*(NCHEBD+1),NFUND,NCELLD),stat=i_stat)
         call memocc(i_stat,product(shape(THETASNEW))*kind(THETASNEW),'THETASNEW','allocate_expansion')
         THETASNEW = 0.D0
         allocate(CLEB(NCLEB,2),stat=i_stat)
         call memocc(i_stat,product(shape(CLEB))*kind(CLEB),'CLEB','allocate_expansion')
         CLEB = 0.D0
         allocate(LOFLM(LM2D),stat=i_stat)
         call memocc(i_stat,product(shape(LOFLM))*kind(LOFLM),'LOFLM','allocate_expansion')
         LOFLM = 0

      else
         if (allocated(WG)) then
            i_all=-product(shape(WG))*kind(WG)
            deallocate(WG,stat=i_stat)
            call memocc(i_stat,i_all,'WG','allocate_expansion')
         endif
         if (allocated(YRG)) then
            i_all=-product(shape(YRG))*kind(YRG)
            deallocate(YRG,stat=i_stat)
            call memocc(i_stat,i_all,'YRG','allocate_expansion')
         endif
         if (allocated(THETAS)) then
            i_all=-product(shape(THETAS))*kind(THETAS)
            deallocate(THETAS,stat=i_stat)
            call memocc(i_stat,i_all,'THETAS','allocate_expansion')
         endif
         if (allocated(THETASNEW)) then
            i_all=-product(shape(THETASNEW))*kind(THETASNEW)
            deallocate(THETASNEW,stat=i_stat)
            call memocc(i_stat,i_all,'THETASNEW','allocate_expansion')
         endif
         if (allocated(CLEB)) then
            i_all=-product(shape(CLEB))*kind(CLEB)
            deallocate(CLEB,stat=i_stat)
            call memocc(i_stat,i_all,'CLEB','allocate_expansion')
         endif
         if (allocated(LOFLM)) then
            i_all=-product(shape(LOFLM))*kind(LOFLM)
            deallocate(LOFLM,stat=i_stat)
            call memocc(i_stat,i_all,'LOFLM','allocate_expansion')
         endif

      endif

   end subroutine allocate_expansion

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_mesh
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the integration mesh
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_mesh(flag,IRM,NATYP,A,B,R,DRDI)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: IRM
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      double precision, dimension(:), allocatable, intent(inout) :: A !< Constants for exponential R mesh
      double precision, dimension(:), allocatable, intent(inout) :: B
      double precision, dimension(:,:), allocatable, intent(inout) :: R !< Radial mesh ( in units a Bohr)
      double precision, dimension(:,:), allocatable, intent(inout) :: DRDI !< Derivative dr/di

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         allocate(DRDI(IRM,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(DRDI))*kind(DRDI),'DRDI','allocate_mesh')
         DRDI = 0.D0
         allocate(R(IRM,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(R))*kind(R),'R','allocate_mesh')
         R = 0.D0
         allocate(A(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(A))*kind(A),'A','allocate_mesh')
         A = 0.D0
         allocate(B(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(B))*kind(B),'B','allocate_mesh')
         B = 0.D0

      else
         if (allocated(DRDI)) then
            i_all=-product(shape(DRDI))*kind(DRDI)
            deallocate(DRDI,stat=i_stat)
            call memocc(i_stat,i_all,'DRDI','allocate_mesh')
         endif
         if (allocated(R)) then
            i_all=-product(shape(R))*kind(R)
            deallocate(R,stat=i_stat)
            call memocc(i_stat,i_all,'R','allocate_mesh')
         endif
         if (allocated(A)) then
            i_all=-product(shape(A))*kind(A)
            deallocate(A,stat=i_stat)
            call memocc(i_stat,i_all,'A','allocate_mesh')
         endif
         if (allocated(B)) then
            i_all=-product(shape(B))*kind(B)
            deallocate(B,stat=i_stat)
            call memocc(i_stat,i_all,'B','allocate_mesh')
         endif

      endif

   end subroutine allocate_mesh

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_pannels
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays that
   !> describe the pannels
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_pannels(flag,NATYP,NTOTD,IPAN,NPAN_TOT,NPAN_EQ_AT,NPAN_LOG_AT,&
      IPAN_INTERVALL,RPAN_INTERVALL)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      integer, intent(in) :: NTOTD
      integer, dimension(:), allocatable, intent(inout) :: IPAN !< Number of panels in non-MT-region
      integer, dimension(:), allocatable, intent(inout) :: NPAN_TOT
      integer, dimension(:), allocatable, intent(inout) :: NPAN_EQ_AT
      integer, dimension(:), allocatable, intent(inout) :: NPAN_LOG_AT
      integer, dimension(:,:), allocatable, intent(inout) :: IPAN_INTERVALL
      double precision, dimension(:,:), allocatable, intent(inout) :: RPAN_INTERVALL

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(IPAN(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IPAN))*kind(IPAN),'IPAN','allocate_pannels')
         IPAN = 0
         allocate(NPAN_TOT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(NPAN_TOT))*kind(NPAN_TOT),'NPAN_TOT','allocate_pannels')
         NPAN_TOT = 0
         allocate(NPAN_EQ_AT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(NPAN_EQ_AT))*kind(NPAN_EQ_AT),'NPAN_EQ_AT','allocate_pannels')
         NPAN_EQ_AT = 0
         allocate(NPAN_LOG_AT(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(NPAN_LOG_AT))*kind(NPAN_LOG_AT),'NPAN_LOG_AT','allocate_pannels')
         NPAN_LOG_AT = 0
         allocate(RPAN_INTERVALL(0:NTOTD,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(RPAN_INTERVALL))*kind(RPAN_INTERVALL),'RPAN_INTERVALL','allocate_pannels')
         RPAN_INTERVALL = 0.D0
         allocate(IPAN_INTERVALL(0:NTOTD,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IPAN_INTERVALL))*kind(IPAN_INTERVALL),'IPAN_INTERVALL','allocate_pannels')
         IPAN_INTERVALL = 0

      else
         if (allocated(IPAN)) then
            i_all=-product(shape(IPAN))*kind(IPAN)
            deallocate(IPAN,stat=i_stat)
            call memocc(i_stat,i_all,'IPAN','allocate_pannels')
         endif
         if (allocated(NPAN_TOT)) then
            i_all=-product(shape(NPAN_TOT))*kind(NPAN_TOT)
            deallocate(NPAN_TOT,stat=i_stat)
            call memocc(i_stat,i_all,'NPAN_TOT','allocate_pannels')
         endif
         if (allocated(NPAN_EQ_AT)) then
            i_all=-product(shape(NPAN_EQ_AT))*kind(NPAN_EQ_AT)
            deallocate(NPAN_EQ_AT,stat=i_stat)
            call memocc(i_stat,i_all,'NPAN_EQ_AT','allocate_pannels')
         endif
         if (allocated(NPAN_LOG_AT)) then
            i_all=-product(shape(NPAN_LOG_AT))*kind(NPAN_LOG_AT)
            deallocate(NPAN_LOG_AT,stat=i_stat)
            call memocc(i_stat,i_all,'NPAN_LOG_AT','allocate_pannels')
         endif
         if (allocated(RPAN_INTERVALL)) then
            i_all=-product(shape(RPAN_INTERVALL))*kind(RPAN_INTERVALL)
            deallocate(RPAN_INTERVALL,stat=i_stat)
            call memocc(i_stat,i_all,'RPAN_INTERVALL','allocate_pannels')
         endif
         if (allocated(IPAN_INTERVALL)) then
            i_all=-product(shape(IPAN_INTERVALL))*kind(IPAN_INTERVALL)
            deallocate(IPAN_INTERVALL,stat=i_stat)
            call memocc(i_stat,i_all,'IPAN_INTERVALL','allocate_pannels')
         endif

      endif

   end subroutine allocate_pannels

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_misc
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of misc arrays
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_misc(flag,NR,IRM,IRID,LMAX,NAEZ,NATYP,NFUND,NREFD,IEMXD,&
      NTOTD,NSHELD,LMMAXD,NEMBD1,NCHEBD,NCELLD,LMXSPD,NSPINDD,NSYMAXD,NPRINCD,IFUNM,&
      IFUNM1,ICHECK,VREF,S,RR,DROR,RNEW,RS,RROT,THESME,DSYMLL,DSYMLL1,LEFTTINVLL,&
      RIGHTTINVLL)

      implicit none

      integer, intent(in) :: NR
      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: IRM
      integer, intent(in) :: IRID
      integer, intent(in) :: LMAX !< Maximum l component in wave function expansion
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: NATYP !< number of kinds of atoms in unit cell
      integer, intent(in) :: NFUND
      integer, intent(in) :: NREFD
      integer, intent(in) :: IEMXD
      integer, intent(in) :: NTOTD
      integer, intent(in) :: LMMAXD
      integer, intent(in) :: NSHELD
      integer, intent(in) :: NEMBD1
      integer, intent(in) :: NCHEBD
      integer, intent(in) :: NCELLD
      integer, intent(in) :: LMXSPD
      integer, intent(in) :: NSPINDD
      integer, intent(in) :: NSYMAXD
      integer, intent(in) :: NPRINCD
      integer, dimension(:,:), allocatable, intent(inout) :: IFUNM
      integer, dimension(:,:), allocatable, intent(inout) :: IFUNM1
      integer, dimension(:,:), allocatable, intent(inout) :: ICHECK
      double precision, dimension(:), allocatable, intent(inout) :: VREF
      double precision, dimension(:,:), allocatable, intent(inout) :: S
      double precision, dimension(:,:), allocatable, intent(inout) :: RR !< Set of real space vectors (in a.u.)
      double precision, dimension(:,:), allocatable, intent(inout) :: DROR
      double precision, dimension(:,:), allocatable, intent(inout) :: RNEW
      double precision, dimension(:,:,:), allocatable, intent(inout) :: RS
      double precision, dimension(:,:,:), allocatable, intent(inout) :: RROT
      double precision, dimension(:,:,:), allocatable, intent(inout) :: THESME
      double complex, dimension(:,:,:), allocatable, intent(inout) :: DSYMLL
      double complex, dimension(:,:,:), allocatable, intent(inout) :: DSYMLL1
      double complex, dimension(:,:,:,:,:), allocatable, intent(inout) :: LEFTTINVLL
      double complex, dimension(:,:,:,:,:), allocatable, intent(inout) :: RIGHTTINVLL

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(RR(3,0:NR),stat=i_stat)
         call memocc(i_stat,product(shape(RR))*kind(RR),'RR','allocate_misc')
         RR = 0.D0
         allocate(VREF(NREFD),stat=i_stat)
         call memocc(i_stat,product(shape(VREF))*kind(VREF),'VREF','allocate_misc')
         VREF = 0.D0
         allocate(DROR(IRM,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(DROR))*kind(DROR),'DROR','allocate_misc')
         DROR = 0.D0
         allocate(S(0:LMAX,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(S))*kind(S),'S','allocate_misc')
         S = 0.D0
         allocate(RROT(48,3,NSHELD),stat=i_stat)
         call memocc(i_stat,product(shape(RROT))*kind(RROT),'RROT','allocate_misc')
         RROT = 0.D0
         allocate(IFUNM(NATYP,LMXSPD),stat=i_stat)
         call memocc(i_stat,product(shape(IFUNM))*kind(IFUNM),'IFUNM','allocate_misc')
         IFUNM = 0
         allocate(IFUNM1(LMXSPD,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(IFUNM1))*kind(IFUNM1),'IFUNM1','allocate_misc')
         IFUNM1 = 0
         allocate(RS(IRM,0:LMAX,NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(RS))*kind(RS),'RS','allocate_misc')
         RS = 0.D0
         allocate(THESME(IRID,NFUND,NCELLD),stat=i_stat)
         call memocc(i_stat,product(shape(THESME))*kind(THESME),'THESME','allocate_misc')
         THESME = 0.D0
         allocate(RNEW(NTOTD*(NCHEBD+1),NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(RNEW))*kind(RNEW),'RNEW','allocate_misc')
         RNEW = 0.D0
         allocate(DSYMLL(LMMAXD,LMMAXD,NSYMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(DSYMLL))*kind(DSYMLL),'DSYMLL','allocate_misc')
         DSYMLL = (0.D0,0.D0)
         allocate(DSYMLL1(LMMAXD,LMMAXD,NSYMAXD),stat=i_stat)
         call memocc(i_stat,product(shape(DSYMLL1))*kind(DSYMLL1),'DSYMLL1','allocate_misc')
         DSYMLL1 = (0.D0,0.D0)
         allocate(ICHECK(NAEZ/NPRINCD,NAEZ/NPRINCD),stat=i_stat)
         call memocc(i_stat,product(shape(ICHECK))*kind(ICHECK),'ICHECK','allocate_misc')
         ICHECK = 0
         allocate(LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), stat=i_stat)
         call memocc(i_stat,product(shape(LEFTTINVLL))*kind(LEFTTINVLL),'LEFTTINVLL','allocate_misc')
         LEFTTINVLL = (0.D0,0.D0)
         allocate(RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), stat=i_stat)
         call memocc(i_stat,product(shape(RIGHTTINVLL))*kind(RIGHTTINVLL),'RIGHTTINVLL','allocate_misc')
         RIGHTTINVLL = (0.D0,0.D0)

      else
         if (allocated(RR)) then
            i_all=-product(shape(RR))*kind(RR)
            deallocate(RR,stat=i_stat)
            call memocc(i_stat,i_all,'RR','allocate_misc')
         endif
         if (allocated(VREF)) then
            i_all=-product(shape(VREF))*kind(VREF)
            deallocate(VREF,stat=i_stat)
            call memocc(i_stat,i_all,'VREF','allocate_misc')
         endif
         if (allocated(DROR)) then
            i_all=-product(shape(DROR))*kind(DROR)
            deallocate(DROR,stat=i_stat)
            call memocc(i_stat,i_all,'DROR','allocate_misc')
         endif
         if (allocated(S)) then
            i_all=-product(shape(S))*kind(S)
            deallocate(S,stat=i_stat)
            call memocc(i_stat,i_all,'S','allocate_misc')
         endif
         if (allocated(RROT)) then
            i_all=-product(shape(RROT))*kind(RROT)
            deallocate(RROT,stat=i_stat)
            call memocc(i_stat,i_all,'RROT','allocate_misc')
         endif
         if (allocated(IFUNM)) then
            i_all=-product(shape(IFUNM))*kind(IFUNM)
            deallocate(IFUNM,stat=i_stat)
            call memocc(i_stat,i_all,'IFUNM','allocate_misc')
         endif
         if (allocated(IFUNM1)) then
            i_all=-product(shape(IFUNM1))*kind(IFUNM1)
            deallocate(IFUNM1,stat=i_stat)
            call memocc(i_stat,i_all,'IFUNM1','allocate_misc')
         endif
         if (allocated(RS)) then
            i_all=-product(shape(RS))*kind(RS)
            deallocate(RS,stat=i_stat)
            call memocc(i_stat,i_all,'RS','allocate_misc')
         endif
         if (allocated(THESME)) then
            i_all=-product(shape(THESME))*kind(THESME)
            deallocate(THESME,stat=i_stat)
            call memocc(i_stat,i_all,'THESME','allocate_misc')
         endif
         if (allocated(RNEW)) then
            i_all=-product(shape(RNEW))*kind(RNEW)
            deallocate(RNEW,stat=i_stat)
            call memocc(i_stat,i_all,'RNEW','allocate_misc')
         endif
         if (allocated(DSYMLL)) then
            i_all=-product(shape(DSYMLL))*kind(DSYMLL)
            deallocate(DSYMLL,stat=i_stat)
            call memocc(i_stat,i_all,'DSYMLL','allocate_misc')
         endif
         if (allocated(DSYMLL1)) then
            i_all=-product(shape(DSYMLL1))*kind(DSYMLL1)
            deallocate(DSYMLL1,stat=i_stat)
            call memocc(i_stat,i_all,'DSYMLL1','allocate_misc')
         endif
         if (allocated(ICHECK)) then
            i_all=-product(shape(ICHECK))*kind(ICHECK)
            deallocate(ICHECK,stat=i_stat)
            call memocc(i_stat,i_all,'ICHECK','allocate_misc')
         endif
         if (allocated(LEFTTINVLL)) then
            i_all=-product(shape(LEFTTINVLL))*kind(LEFTTINVLL)
            deallocate(LEFTTINVLL,stat=i_stat)
            call memocc(i_stat,i_all,'LEFTTINVLL','allocate_misc')
         endif
         if (allocated(RIGHTTINVLL)) then
            i_all=-product(shape(RIGHTTINVLL))*kind(RIGHTTINVLL)
            deallocate(RIGHTTINVLL,stat=i_stat)
            call memocc(i_stat,i_all,'RIGHTTINVLL','allocate_misc')
         endif

      endif

   end subroutine allocate_misc

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_green
   !
   ! DESCRIPTION:
   !> @brief subroutine handling the allocation/deallocation of arrays handling
   !> the Green functions
   !
   !> @author
   !> Jonathan Chico
   !> @date 19.12.2017
   !----------------------------------------------------------------------------
   subroutine allocate_green(flag,NAEZ,IEMXD,NGSHD,NSHELD,LMPOT,NOFGIJD,ISH,JSH,&
      KMESH,IMAXSH,IQCALC,IOFGIJ,JOFGIJ,IJTABSH,IJTABSYM,IJTABCALC,IJTABCALC_I,ILM_MAP,GSH)

      implicit none

      integer, intent(in) :: flag ! Allocate/deallocate (1/-1) arrays
      integer, intent(in) :: NAEZ !< number of atoms in unit cell
      integer, intent(in) :: IEMXD
      integer, intent(in) :: NGSHD
      integer, intent(in) :: NSHELD
      integer, intent(in) :: LMPOT
      integer, intent(in) :: NOFGIJD
      integer, dimension(:,:), allocatable, intent(inout) :: ISH
      integer, dimension(:,:), allocatable, intent(inout) :: JSH
      integer, dimension(:), allocatable, intent(inout) :: KMESH
      integer, dimension(:), allocatable, intent(inout) :: IMAXSH
      integer, dimension(:), allocatable, intent(inout) :: IQCALC
      integer, dimension(:), allocatable, intent(inout) :: IOFGIJ !< Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
      integer, dimension(:), allocatable, intent(inout) :: JOFGIJ !< Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
      integer, dimension(:), allocatable, intent(inout) :: IJTABSH !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
      integer, dimension(:), allocatable, intent(inout) :: IJTABSYM !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
      integer, dimension(:), allocatable, intent(inout) :: IJTABCALC !< Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
      integer, dimension(:), allocatable, intent(inout) :: IJTABCALC_I
      integer, dimension(:,:), allocatable, intent(inout) :: ILM_MAP
      double precision, dimension(:), allocatable, intent(inout) :: GSH

      integer :: i_stat, i_all

      if (flag>0) then

         allocate(GSH(NGSHD),stat=i_stat)
         call memocc(i_stat,product(shape(GSH))*kind(GSH),'GSH','allocate_green')
         GSH = 0.D0
         allocate(KMESH(IEMXD),stat=i_stat)
         call memocc(i_stat,product(shape(KMESH))*kind(KMESH),'KMESH','allocate_green')
         KMESH = 0
         allocate(ILM_MAP(NGSHD,3),stat=i_stat)
         call memocc(i_stat,product(shape(ILM_MAP))*kind(ILM_MAP),'ILM_MAP','allocate_green')
         ILM_MAP = 0
         allocate(IQCALC(NAEZ),stat=i_stat)
         call memocc(i_stat,product(shape(IQCALC))*kind(IQCALC),'IQCALC','allocate_green')
         IQCALC = 0
         allocate(JOFGIJ(NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(JOFGIJ))*kind(JOFGIJ),'JOFGIJ','allocate_green')
         JOFGIJ = 0
         allocate(IOFGIJ(NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(IOFGIJ))*kind(IOFGIJ),'IOFGIJ','allocate_green')
         IOFGIJ = 0
         allocate(IMAXSH(0:LMPOT),stat=i_stat)
         call memocc(i_stat,product(shape(IMAXSH))*kind(IMAXSH),'IMAXSH','allocate_green')
         IMAXSH = 0
         allocate(IJTABSH(NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(IJTABSH))*kind(IJTABSH),'IJTABSH','allocate_green')
         IJTABSH = 0
         allocate(IJTABSYM(NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(IJTABSYM))*kind(IJTABSYM),'IJTABSYM','allocate_green')
         IJTABSYM = 0
         allocate(IJTABCALC(NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(IJTABCALC))*kind(IJTABCALC),'IJTABCALC','allocate_green')
         IJTABCALC = 0
         allocate(ISH(NSHELD,NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(ISH))*kind(ISH),'ISH','allocate_green')
         ISH =0
         allocate(JSH(NSHELD,NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(JSH))*kind(JSH),'JSH','allocate_green')
         JSH = 0
         allocate(IJTABCALC_I(NOFGIJD),stat=i_stat)
         call memocc(i_stat,product(shape(IJTABCALC_I))*kind(IJTABCALC_I),'IJTABCALC_I','allocate_green')
         IJTABCALC_I = 0

      else

         if (allocated(GSH)) then
            i_all=-product(shape(GSH))*kind(GSH)
            deallocate(GSH,stat=i_stat)
            call memocc(i_stat,i_all,'GSH','allocate_misc')
         endif
         if (allocated(KMESH)) then
            i_all=-product(shape(KMESH))*kind(KMESH)
            deallocate(KMESH,stat=i_stat)
            call memocc(i_stat,i_all,'KMESH','allocate_misc')
         endif
         if (allocated(ILM_MAP)) then
            i_all=-product(shape(ILM_MAP))*kind(ILM_MAP)
            deallocate(ILM_MAP,stat=i_stat)
            call memocc(i_stat,i_all,'ILM_MAP','allocate_misc')
         endif
         if (allocated(IQCALC)) then
            i_all=-product(shape(IQCALC))*kind(IQCALC)
            deallocate(IQCALC,stat=i_stat)
            call memocc(i_stat,i_all,'IQCALC','allocate_misc')
         endif
         if (allocated(JOFGIJ)) then
            i_all=-product(shape(JOFGIJ))*kind(JOFGIJ)
            deallocate(JOFGIJ,stat=i_stat)
            call memocc(i_stat,i_all,'JOFGIJ','allocate_misc')
         endif
         if (allocated(IOFGIJ)) then
            i_all=-product(shape(IOFGIJ))*kind(IOFGIJ)
            deallocate(IOFGIJ,stat=i_stat)
            call memocc(i_stat,i_all,'IOFGIJ','allocate_misc')
         endif
         if (allocated(IMAXSH)) then
            i_all=-product(shape(IMAXSH))*kind(IMAXSH)
            deallocate(IMAXSH,stat=i_stat)
            call memocc(i_stat,i_all,'IMAXSH','allocate_misc')
         endif
         if (allocated(IJTABSH)) then
            i_all=-product(shape(IJTABSH))*kind(IJTABSH)
            deallocate(IJTABSH,stat=i_stat)
            call memocc(i_stat,i_all,'IJTABSH','allocate_misc')
         endif
         if (allocated(IJTABSYM)) then
            i_all=-product(shape(IJTABSYM))*kind(IJTABSYM)
            deallocate(IJTABSYM,stat=i_stat)
            call memocc(i_stat,i_all,'IJTABSYM','allocate_misc')
         endif
         if (allocated(IJTABCALC)) then
            i_all=-product(shape(IJTABCALC))*kind(IJTABCALC)
            deallocate(IJTABCALC,stat=i_stat)
            call memocc(i_stat,i_all,'IJTABCALC','allocate_misc')
         endif
         if (allocated(ISH)) then
            i_all=-product(shape(ISH))*kind(ISH)
            deallocate(ISH,stat=i_stat)
            call memocc(i_stat,i_all,'ISH','allocate_misc')
         endif
         if (allocated(JSH)) then
            i_all=-product(shape(JSH))*kind(JSH)
            deallocate(JSH,stat=i_stat)
            call memocc(i_stat,i_all,'JSH','allocate_misc')
         endif
         if (allocated(IJTABCALC_I)) then
            i_all=-product(shape(IJTABCALC_I))*kind(IJTABCALC_I)
            deallocate(IJTABCALC_I,stat=i_stat)
            call memocc(i_stat,i_all,'IJTABCALC_I','allocate_misc')
         endif
      endif

   end subroutine allocate_green

end module memoryhandling
