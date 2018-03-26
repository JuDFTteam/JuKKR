!-------------------------------------------------------------------------------
! SUBROUTINE: VXCGGA
!> @brief Add the exchange-correlation-potential given by GGA to the potential
!> and if total energies should be calculated (KTE=1) the
!> exchange-correlation-energies are calculated.
!> @details Use as input the charge density times \f$r^2\f$ (rho2ns(...,1)) and
!> in the spin-polarized case (NSPIN=2) the spin density times \f$r^2\f$
!> (rho2ns(...,2)) .
!> The density times \f$4\pi\f$ is generated at an angular mesh.
!> The exchange-correlation potential and the exchange-correlation
!> energy are calculated at those mesh points with a subroutine.
!> In the paramagnetic case the "spin-density" is set equal zero.
!> After that the exchange-correlation potential and in the case of
!> total energies (KTE=1) the exchange-correlation energy are
!> expanded into spherical harmonics.
!> The ex.-cor. potential is added to the given potential.
!> The expansion into spherical harmonics uses the orthogonality
!> of these harmonics.
!> - Therefore a gauss-legendre integration for \f$\theta\f$ and a
!> gauss-tschebyscheff integration for \f$\phi\f$ is used.
!>
!> All needed values for the angular mesh and angular integration
!> are generate in the subroutine sphere.
!> The ex.-cor. potential is extrapolated to the origin only
!> for the lm=1 value .
!> @author B. Drittler, R. Zeller
!> @note
!> - B. Drittler Oct. 1989: Modified for shape functions
!> - R. Zeller Nov. 1993: simplified and modified for Paragon X/PS
!> - R. Zeller 23/6/1996: cor error
!> - Jonathan Chico: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine VXCGGA(EXC,KTE,KXC,LMAX,NSPIN,IATYP,RHO2NS,V,R,DRDI,A,&
      IRWS,IRCUT,IPAN,KSHAPE,GSH,ILM,IMAXSH,&
      IFUNM,THETAS,WTYR,IJEND,LMSP,THET,YLM,DYLMT1,&
      DYLMT2,DYLMF1,DYLMF2,DYLMTF,NGSHD,LMPOT,NFUND,LMXSPD)

   use Constants

   implicit none

   ! .. Input variables
   integer, intent(in) :: KTE    !< Calculation of the total energy On/Off (1/0)
   integer, intent(in) :: KXC    !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
   integer, intent(in) :: LMAX   !< Maximum l component in wave function expansion
   integer, intent(in) :: IRWS   !< IATYP Entry in the IRWS array with the R point at WS radius
   integer, intent(in) :: IPAN   !< IATYP Entry in the IPAN array with the number of panels in non-MT-region
   integer, intent(in) :: IATYP
   integer, intent(in) :: IJEND
   integer, intent(in) :: NGSHD  !< Shape functions parameters in non-spherical part
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: NFUND
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: IPAND  !< Number of panels in non-spherical part
   integer, intent(in) :: LMMAX  !< (LMAX+1)^2
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   integer, intent(in) :: LMXSPD !< (2*LPOT+1)**2
   double precision, intent(in):: A !< IATYP entry for the array A with the constants for exponential R mesh
   ! .. Array Arguments
   integer, dimension(LMXSPD), intent(in)    :: LMSP  !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(0:IPAND), intent(in)   :: IRCUT !< R points of panel borders
   integer, dimension(LMXSPD), intent(in)    :: IFUNM
   integer, dimension(0:LMPOTD), intent(in)  :: IMAXSH
   integer, dimension(NGSHD,3), intent(in)   :: ILM
   double precision, dimension(IRM), intent(in)    :: R  !< IATYP entry of the radial mesh ( in units a Bohr)
   double precision, dimension(NGSHD), intent(in)  :: GSH
   double precision, dimension(IRM), intent(in)    :: DRDI !< IATYP entry of the derivative dr/di
   double precision, dimension(IJEND), intent(in)  :: THET
   double precision, dimension(IJEND,LMPOT), intent(in)  :: YLM
   double precision, dimension(IJEND,LMPOT), intent(in)  :: WTYR
   double precision, dimension(IRID,NFUND), intent(in)   :: THETAS  !< IATYP entry of the shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   double precision, dimension(IJEND,LMPOT), intent(in)  :: DYLMF1
   double precision, dimension(IJEND,LMPOT), intent(in)  :: DYLMF2
   double precision, dimension(IJEND,LMPOT), intent(in)  :: DYLMT1
   double precision, dimension(IJEND,LMPOT), intent(in)  :: DYLMT2
   double precision, dimension(IJEND,LMPOT), intent(in)  :: DYLMTF
   double precision, dimension(IRM,LMPOT,2), intent(in)  :: RHO2NS   !< radial density
   ! .. Input/Output variables
   double precision, dimension(0:LPOT,NATYP), intent(inout) :: EXC   !< exchange correlation energy
   double precision, dimension(IRM,LMPOT,2), intent(inout) :: V
   ! .. Local Scalars
   double precision :: VXC1,VXC2,VXC3,ZERO,ZERO1
   double precision :: CHGDEN,DX,ELMXC,FPI,R1,R2,RPOINT,SPIDEN,VLMXC
   integer :: LM2,LMMAX,M,MESH,NSPIN2
   integer :: IFUN,IPAN1,IPOT,IR,IRC0,IRC1,IRH,IRS1,ISPIN,J,L,L1MAX,LM
   ! .. Local Arrays
   double precision, dimension(IJEND) :: EXCIJ
   double precision, dimension(IRM,0:LPOT)   :: ER
   double precision, dimension(IJEND,2)      :: VXC
   double precision, dimension(IRM,LMPOT)    :: DRRL
   double precision, dimension(2:3,2)        :: VXCR
   double precision, dimension(IRM,LMPOT)    :: ESTOR
   double precision, dimension(LMPOT,2)      :: RHOLM
   double precision, dimension(IRM,LMPOT)    :: DRRUL
   double precision, dimension(IRM,LMPOT)    :: DDRRL
   double precision, dimension(IRM,LMPOT)    :: DDRRUL
   double precision, dimension(IRM,2,LMPOT)  :: RHOL
   ! .. External Functions
   double precision :: DDOT
   external DDOT
   ! .. External Subroutines ..
   external :: GRADRL,MKXCPE,SIMP3,SIMPK,MKXCPE2
   ! .. Intrinsic Functions ..
   intrinsic :: ABS,ATAN,MOD
   ! .. Data statements ..
   data ZERO,ZERO1/0.d0,1.d-12/
   !     ..
   write (1337,FMT=*) ' GGA CALCULATION '
   FPI = 4.0D0*PI
   !----------------------------------------------------------------------------
   ! Loop over given representive atoms
   !----------------------------------------------------------------------------
   if (KSHAPE.NE.0) then
      IPAN1 = IPAN
      IRC1 = IRCUT(IPAN)
      IRS1 = IRCUT(1)
      IRC0 = 2
      if (KREL.EQ.1) stop ' REL + FULL POTENTIAL N/A '
   else
      IRC1 = IRWS
      IRS1 = IRC1
      IPAN1 = 1
      IRC0 = 2
      if (KREL.EQ.1) IRC0 = 2 + MOD(IRCUT(1),2)
   end if

   do ISPIN = 1,NSPIN
      VXCR(2,ISPIN) = 0.0D0
      VXCR(3,ISPIN) = 0.0D0
   enddo
   !----------------------------------------------------------------------------
   ! Initialize for ex.-cor. energy
   !----------------------------------------------------------------------------
   if (KTE.EQ.1) then
      do L = 0,LMAX
         EXC(L,IATYP) = 0.0D0
         do IR = 1,IRC1
            ER(IR,L) = 0.0D0
         enddo ! IR
      enddo ! L
      !
      do LM = 1,LMMAX
         do IR = 1,IRC1
            ESTOR(IR,LM) = 0.0D0
         enddo ! IR
      enddo ! LM
   end if
   !
   L1MAX = LMAX + 1
   MESH = IRWS
   DX = A
   !
   if (NSPIN.EQ.2) then
      do LM = 1,LMMAX
         do IR = 2,MESH
            R1 = R(IR)
            R2 = R1*R1
            CHGDEN = RHO2NS(IR,LM,1)/R2
            SPIDEN = RHO2NS(IR,LM,2)/R2
            if (ABS(CHGDEN).LE.ZERO1) CHGDEN = ZERO
            if (ABS(SPIDEN).LE.ZERO1) SPIDEN = ZERO
            RHOL(IR,2,LM) = (CHGDEN+SPIDEN)/2.d0
            RHOL(IR,1,LM) = (CHGDEN-SPIDEN)/2.d0
         enddo ! IR
         ! extrapolate
         RHOL(1,1,LM) = RHOL(2,1,LM)
         RHOL(1,2,LM) = RHOL(2,2,LM)
      enddo ! LM
   else
      !
      do LM = 1,LMMAX
         do IR = 2,MESH
            R1 = R(IR)
            R2 = R1*R1
            !
            CHGDEN = RHO2NS(IR,LM,1)/R2
            if (ABS(CHGDEN).LE.ZERO1) CHGDEN = ZERO
            RHOL(IR,1,LM) = CHGDEN/2.d0
            RHOL(IR,2,LM) = CHGDEN/2.d0
         enddo ! IR
         ! extrapolate
         RHOL(1,1,LM) = RHOL(2,1,LM)
         RHOL(1,2,LM) = RHOL(2,2,LM)
      enddo ! LM
   end if

   call GRADRL(NSPIN,MESH,L1MAX,DX,RHOL,R,DRDI,IPAN1,IPAND,IRCUT, &
      DRRL,DDRRL,DRRUL,DDRRUL,IRMD,LMPOTD)

   !----------------------------------------------------------------------------
   ! Loop over radial mesh
   !----------------------------------------------------------------------------

   do IR = IRC0,IRC1
      RPOINT = R(IR)
      !-------------------------------------------------------------------------
      ! Calculate the ex.-cor. potential
      !-------------------------------------------------------------------------
      NSPIN2 = 2

      do ISPIN = 1,NSPIN2
         do  LM = 1,LMMAX
            RHOLM(LM,ISPIN) = RHOL(IR,ISPIN,LM)
         enddo
      enddo
      !    only for spin-polarized
      !
      ! PW91 functional
      if(KXC.EQ.3)then
         call MKXCPE(NSPIN2,IR,IJEND,L1MAX,RPOINT,RHOLM,VXC,EXCIJ,   &
            THET,YLM,DYLMT1,DYLMT2,DYLMF1,DYLMF2,DYLMTF,DRRL,        &
            DDRRL,DRRUL,DDRRUL,IRMD,LMPOTD)
      ! PBE functional
      elseif(KXC.EQ.4)then
         CALL MKXCPE2(IR,IJEND,RPOINT,RHOLM,VXC,EXCIJ,YLM,DYLMT1, &
            DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,DRRUL,DDRRUL,         &
            IRMD,LMPOTD,LMMAX,.false.)
      ! PBEsol functional
      elseif(KXC.EQ.5)then
         call MKXCPE2(IR,IJEND,RPOINT,RHOLM,VXC,EXCIJ,YLM,DYLMT1, &
            DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,DRRUL,DDRRUL,         &
            IRMD,LMPOTD,LMMAX,.true.)
      else
         write(1337,*) ' KXC ???'
         stop
      endif
      !-------------------------------------------------------------------------
      ! Expand the ex.-cor. potential into spherical harmonics ,
      !   using the orthogonality
      !-------------------------------------------------------------------------
      do ISPIN = 1,NSPIN
         !----------------------------------------------------------------------
         ! Determine the corresponding potential number
         !----------------------------------------------------------------------
         IPOT = ISPIN
         do LM = 1,LMMAX
            VLMXC = DDOT(IJEND,VXC(1,ISPIN),1,WTYR(1,LM),1)
            V(IR,LM,IPOT) = V(IR,LM,IPOT) + VLMXC
            !-------------------------------------------------------------------
            ! Store the ex.-c. potential of ir=2 and =3 for the extrapolation
            !-------------------------------------------------------------------
            if (LM.EQ.1 .AND. (IR.EQ.2.OR.IR.EQ.3)) VXCR(IR,ISPIN) = VLMXC
         enddo ! LM
      enddo ! ISPIN
      !-------------------------------------------------------------------------
      ! File er in case of total energies
      !-------------------------------------------------------------------------
      if (KTE.EQ.1) then
         !----------------------------------------------------------------------
         ! Expand ex.-cor. energy into spherical harmonics
         !   using the orthogonality
         !----------------------------------------------------------------------
         do L = 0,LMAX
            do M = -L,L
               LM = L*L + L + M + 1
               ELMXC = DDOT(IJEND,EXCIJ,1,WTYR(1,LM),1)
               !----------------------------------------------------------------
               ! Multiply the lm-component of the ex.-cor. energy with the same
               ! lm-component of the charge density times r**2 and sum over lm
               ! this corresponds to a integration over the angular .
               !----------------------------------------------------------------
               if ((KSHAPE.NE.0) .AND. (IR.GT.IRS1)) then
                  ESTOR(IR,LM) = ELMXC
               else
                  ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*ELMXC
               end if
            enddo ! M
         enddo ! L
      end if
   enddo !IR
   !----------------------------------------------------------------------------
   ! Integrate er in case of total energies to get exc
   !----------------------------------------------------------------------------
   if (KTE.EQ.1) then
      if (KSHAPE.EQ.0) then
         do L = 0,LMAX
            call SIMP3(ER(1,L),EXC(L,IATYP),1,IRS1,DRDI)
         enddo
      else
         do L = 0,LMAX
            do M = -L,L
               LM = L*L + L + M + 1
               !----------------------------------------------------------------
               ! Convolute with shape function
               !----------------------------------------------------------------
               do J = IMAXSH(LM-1) + 1,IMAXSH(LM)
                  LM2 = ILM(J,2)
                  if (LMSP(ILM(J,3)).GT.0) then
                     IFUN = IFUNM(ILM(J,3))
                     do IR = IRS1 + 1,IRC1
                        IRH = IR - IRS1
                        ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*GSH(J)*   &
                                    THETAS(IRH,IFUN)*ESTOR(IR,LM2)
                     enddo ! IR
                  end if
               enddo ! J
            enddo ! M
            call SIMPK(ER(1,L),EXC(L,IATYP),IPAN1,IRCUT,DRDI)
         enddo ! L
      end if
   end if
   !----------------------------------------------------------------------------
   ! Extrapolate ex.-cor potential to the origin only for lm=1
   !----------------------------------------------------------------------------
   do ISPIN = 1,NSPIN
      IPOT = ISPIN
      !
      VXC2 = VXCR(2,ISPIN)
      VXC3 = VXCR(3,ISPIN)
      VXC1 = VXC2 - R(2)* (VXC3-VXC2)/ (R(3)-R(2))
      !
      V(1,1,IPOT) = V(1,1,IPOT) + VXC1
   enddo
   !
end subroutine VXCGGA
