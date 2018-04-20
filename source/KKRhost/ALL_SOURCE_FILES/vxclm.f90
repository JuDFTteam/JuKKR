!-------------------------------------------------------------------------------
! SUBROUTINE: VXCLM
!> @brief Add the exchange-correlation-potential to the given potential
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
subroutine VXCLM(EXC,KTE,KXC,LMAX,NSPIN,IATYP,RHO2NS,V,R,DRDI,IRWS,IRCUT,IPAN,   &
   KSHAPE,GSH,ILM,IMAXSH,IFUNM,THETAS,YR,WTYR,IJEND,LMSP,LMPOT,LMXSPD,LMMAX,IRM, &
   LPOT,NATYP)

   use Constants
   use global_variables

   implicit none

   ! .. Scalar Arguments
   integer, intent(in) :: IRM    !< Maximum number of radial points
   integer, intent(in) :: KTE    !< Calculation of the total energy On/Off (1/0)
   integer, intent(in) :: KXC    !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
   integer, intent(in) :: LMAX   !< Maximum l component in wave function expansion
   integer, intent(in) :: IRWS   !< IATYP Entry in the IRWS array with the R point at WS radius
   integer, intent(in) :: IPAN   !< IATYP Entry in the IPAN array with the number of panels in non-MT-region
   integer, intent(in) :: LPOT   !< Maximum l component in potential expansion
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: IATYP
   integer, intent(in) :: IJEND
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: LMMAX  !< (LMAX+1)^2
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   integer, intent(in) :: LMXSPD !< (2*LPOT+1)**2
   ! .. Array Arguments
   integer, dimension(LMXSPD), intent(in)    :: LMSP  !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(0:IPAND), intent(in)   :: IRCUT !< R points of panel borders
   integer, dimension(LMXSPD), intent(in)    :: IFUNM
   integer, dimension(0:LMPOT), intent(in)   :: IMAXSH
   integer, dimension(NGSHD,3), intent(in)   :: ILM
   double precision, dimension(IRM), intent(in)          :: R        !< Radial mesh ( in units a Bohr)
   double precision, dimension(IJEND,LMPOT), intent(in)  :: YR
   double precision, dimension(NGSHD), intent(in)        :: GSH
   double precision, dimension(IRM), intent(in)          :: DRDI     !< Derivative dr/di
   double precision, dimension(IJEND,LMPOT), intent(in)  :: WTYR
   double precision, dimension(IRID,NFUND), intent(in)   :: THETAS   !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   double precision, dimension(IRM,LMPOT,2), intent(in)  :: RHO2NS   !< radial density
   ! .. Input/Output variables
   double precision, dimension(0:LPOT,NATYP), intent(inout) :: EXC !< exchange correlation energy
   double precision, dimension(IRM,LMPOT,2), intent(inout)  :: V
   ! .. Local Scalars
   integer :: IFUN,IJ,IPOT,IR,IRC1,IRH,IRS1,IS,ISPIN,J,L,LM,LM2,M
   double precision :: ELMXC,FPI,FPIPR2,VLMXC,VXC1,VXC2,VXC3,factor
   ! .. Local Arrays
   double precision, dimension(IJEND)        :: EXCIJ
   double precision, dimension(IRM,0:LPOT)   :: ER
   double precision, dimension(IJEND,2)      :: VXC
   double precision, dimension(2:3,2)        :: VXCR
   double precision, dimension(IJEND,2)      :: FPRHO
   double precision, dimension(IRM,LMPOT)    :: ESTOR
   ! .. External Functions
   double precision :: DDOT
   external :: DDOT
   ! .. External Subroutines ..
   external :: DAXPY,SIMP3,SIMPK,VOSKO,VXCSPO
   ! ..
   write(1337,*) 'Including cutoff of vxc for small density'
   FPI = 4.0D0*PI

   !----------------------------------------------------------------------------
   ! Loop over given representive atoms
   !----------------------------------------------------------------------------
   if (KSHAPE.NE.0) then
      IRC1 = IRCUT(IPAN)
      IRS1 = IRCUT(1)
   else
      IRC1 = IRWS
      IRS1 = IRC1
   end if

   do ISPIN = 1,NSPIN
      VXCR(2,ISPIN) = 0.0D0
      VXCR(3,ISPIN) = 0.0D0
   enddo ! ISPIN
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
   !----------------------------------------------------------------------------
   ! Loop over radial mesh
   !----------------------------------------------------------------------------
   do IR = 2,IRC1
      !-------------------------------------------------------------------------
      ! Generate the densities on an angular mesh
      !-------------------------------------------------------------------------
      do IS = 1,2
         do IJ = 1,IJEND
            FPRHO(IJ,IS) = 0.D0
         enddo ! IJ
      enddo ! IS

      FPIPR2 = FPI/R(IR)**2
      do ISPIN = 1,NSPIN
         do LM = 1,LMMAX
            call DAXPY(IJEND,RHO2NS(IR,LM,ISPIN)*FPIPR2,YR(1,LM),1,FPRHO(1,ISPIN),1)
         enddo ! LM
      enddo
      !-------------------------------------------------------------------------
      ! Calculate the ex.-cor. potential
      !-------------------------------------------------------------------------
      if (KXC.LE.1) then
         call VXCSPO(EXCIJ,FPRHO,VXC,KXC,IJEND,IJEND)
      else
         call VOSKO(EXCIJ,FPRHO,VXC,IJEND,IJEND)
      end if

      do IJ=1,IJEND
         factor = (1.d0-dexp(-dabs(fprho(IJ,1))*1000.d0))
         do ISPIN=1,NSPIN
            VXC(IJ,ISPIN) =VXC(IJ,ISPIN) * factor  !cutoff
         enddo ! ISPIN
      enddo !IJ
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
            ! store the ex.-c. potential of ir=2 and =3 for the extrapolation
            !-------------------------------------------------------------------
            if (LM.EQ.1 .AND. (IR.EQ.2.OR.IR.EQ.3)) VXCR(IR,ISPIN) = VLMXC
         enddo ! LM
      enddo ! ISPIN
      !-------------------------------------------------------------------------
      ! File er in case of total energies
      !-------------------------------------------------------------------------
      if (KTE.EQ.1) then
         !----------------------------------------------------------------------
         ! expand ex.-cor. energy into spherical harmonics
         !   using the orthogonality
         !----------------------------------------------------------------------
         do L = 0,LMAX
            do M = -L,L
               LM = L*L + L + M + 1
               ELMXC = DDOT(IJEND,EXCIJ,1,WTYR(1,LM),1)
               !----------------------------------------------------------------
               ! multiply the lm-component of the ex.-cor. energy with the same
               ! lm-component of the charge density times r**2 and sum over lm
               ! this corresponds to a integration over the angular .
               !----------------------------------------------------------------
               if ((KSHAPE.NE.0) .AND. (IR.GT.IRS1)) then
                  ESTOR(IR,LM) = ELMXC
               else
                  ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*ELMXC
               end if
            enddo !M
         enddo !L
      end if
   enddo !IR

   !----------------------------------------------------------------------------
   ! Integrate er in case of total energies to get exc
   !----------------------------------------------------------------------------
   if (KTE.EQ.1) then
      if (KSHAPE.EQ.0) then
         do L = 0,LMAX
            call SIMP3(ER(1,L),EXC(L,IATYP),1,IRS1,DRDI)
         enddo !L
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
            call SIMPK(ER(1,L),EXC(L,IATYP),IPAN,IRCUT,DRDI)
         enddo ! L
      end if
   end if
   !----------------------------------------------------------------------------
   ! Extrapolate ex.-cor potential to the origin only for lm=1
   !----------------------------------------------------------------------------
   do ISPIN = 1,NSPIN
      IPOT = ISPIN
      VXC2 = VXCR(2,ISPIN)
      VXC3 = VXCR(3,ISPIN)
      VXC1 = VXC2 - R(2)* (VXC3-VXC2)/ (R(3)-R(2))
      V(1,1,IPOT) = V(1,1,IPOT) + VXC1
   enddo ! ISPIN
   !
end subroutine VXCLM
