!-------------------------------------------------------------------------------
! SUBROUTINE: VINTRAS
!> @brief Calculate the electron-intracell-potentials and the charge-moments
!> of given charge densities. ( For each spin-direction the potential is the
!> same in the polarized case.)
!> @details Initialize the potential \f$ V\f$ with the electron-intracell-potentials
!> the intracell-potential is expanded into spherical harmonics .
!> the lm-term of the intracell-potential of the representive atom i is given by
!>
!> \f$V\left(r,lm,i\right)=\frac{8\pi}{2l+1} \left(\int_{0}^{r}dr'\frac{r'^{l}}{r^{l+1}} rho2ns(r',lm,i,1) + \int_{r}^{r_{cut}}dr' \frac{r^l}{r'^{l+1}}rho2ns(r',lm,i,1) ) \right)\f$
!>
!> the lm contribution of the charge moment of the representive atom i is given by
!>
!> \f$ cmom\left(lm,i\right)=\int_{0}^{r_{cut}} dr' r'^{l}rho2ns(r',lm,i,1)\f$
!>
!>  (see notes by b.drittler and u.klemradt) \f$ r_{cut}\f$ is muffin tin or
!> Wigner-Seitz sphere radius, depending on kshape turned on or off
!
!> @note Attention : \f$ rho2ns(...,1)\f$ is the real charge density times \f$r^2\f$
!> developed into spherical harmonics . (see deck rholm)
!>
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!> @author B. Drittler
!> @date May 1987
!-----------------------------------------------------------------------

subroutine VINTRAS(CMOM,CMINST,LMAX,NSPIN,NSTART,NEND,RHO2NS,V,R, &
      DRDI,IRWS,IRCUT,IPAN,KSHAPE,NTCELL,ILM,IFUNM,IMAXSH,GSH,THETAS,LMSP)

   use Constants

   implicit none

   ! .. Input Variables
   integer, intent(in) :: IRM   !< Maximum number of radial points
   integer, intent(in) :: LMAX   !< Maximum l component in wave function expansion
   integer, intent(in) :: NEND
   integer, intent(in) :: NFUND
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: IPAND  !< Number of panels in non-spherical part
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NGSHD  !< Shape functions parameters in non-spherical part
   integer, intent(in) :: NSTART
   integer, intent(in) :: NCELLD !< Number of cells (shapes) in non-spherical part
   integer, intent(in) :: LMXSPD !< (2*LPOT+1)**2
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   integer, dimension(NATYP), intent(in)     :: IRWS !< R point at WS radius
   integer, dimension(NATYP), intent(in)     :: IPAN   !< Number of panels in non-MT-region
   integer, dimension(NATYP), intent(in)     :: NTCELL !< Index for WS cell
   integer, dimension(0:LMPOT), intent(in)   :: IMAXSH
   integer, dimension(NGSHD,3), intent(in)         :: ILM
   integer, dimension(NATYP,LMXSPD), intent(in)    :: LMSP !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(NATYP,LMXSPD), intent(in)    :: IFUNM
   integer, dimension(0:IPAND,NATYP), intent(in)   :: IRCUT   !< R points of panel borders

   double precision, dimension(NGSHD), intent(in) :: GSH
   double precision, dimension(IRM,NATYP), intent(in) :: R !< Radial mesh ( in units a Bohr)
   double precision, dimension(IRM,NATYP), intent(in) :: DRDI !< Derivative dr/di
   double precision, dimension(IRID,NFUND,NCELLD), intent(in) :: THETAS  !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   double precision, dimension(IRM,LMPOT,NATYP,2), intent(in) :: RHO2NS   !< radial density

   ! .. Output variables
   double precision, dimension(LMPOT,NATYP), intent(out) :: CMOM   !< LM moment of total charge
   double precision, dimension(LMPOT,NATYP), intent(out) :: CMINST
   double precision, dimension(IRM,LMPOT,NPOTD), intent(out) :: V

   ! .. Local Variables
   double precision :: FAC,RL
   integer :: I,IATYP,ICELL,IEND,IFUN,IPOT,IRC1,IRS1,ISTART,J,L,LM,LM2,LM3,M
   ! .. Local Arrays
   integer, dimension(0:IPAND) :: IRCUTM
   double precision, dimension(IRM) :: V1
   double precision, dimension(IRM) :: V2
   double precision, dimension(IRM) :: VINT1
   double precision, dimension(IRM) :: VINT2
   ! .. External Subroutines ..
   external SINWK,SOUTK
   ! .. Intrinsic Functions ..
   intrinsic ATAN,REAL
   !     ..
   do IATYP = NSTART,NEND
      if (KSHAPE.NE.0) then
         IRS1 = IRCUT(1,IATYP)
         IRC1 = IRCUT(IPAN(IATYP),IATYP)
         ICELL = NTCELL(IATYP)
         do I = 0,IPAN(IATYP)
            IRCUTM(I) = IRCUT(I,IATYP)
         enddo ! I
      else
         IRS1 = IRWS(IATYP)
         IRC1 = IRS1
         IRCUTM(0) = IRCUT(0,IATYP)
         IRCUTM(1) = IRC1
      end if
      !-------------------------------------------------------------------------
      ! Determine the right potential numbers
      !-------------------------------------------------------------------------
      IPOT = NSPIN*IATYP

      do L = 0,LMAX
         FAC = 8.0D0*PI/real(2*L+1)
         do M = -L,L
            LM = L*L + L + M + 1
            !-------------------------------------------------------------------
            ! Set up of the integrands v1 and v2
            !-------------------------------------------------------------------
            V1(1) = 0.0D0
            V2(1) = 0.0D0
            do I = 2,IRS1
               RL = R(I,IATYP)**L
               V1(I) = RHO2NS(I,LM,IATYP,1)*RL*DRDI(I,IATYP)
               V2(I) = RHO2NS(I,LM,IATYP,1)/R(I,IATYP)/RL*DRDI(I,IATYP)
            enddo ! I
            !-------------------------------------------------------------------
            ! Convolute charge density of interstial with shape function if kshape.gt.0
            !-------------------------------------------------------------------
            if (KSHAPE.NE.0) then
               do I = IRS1 + 1,IRC1
                  V1(I) = 0.0D0
               enddo ! I
               ISTART = IMAXSH(LM-1) + 1
               IEND = IMAXSH(LM)
               do J = ISTART,IEND
                  LM2 = ILM(J,2)
                  LM3 = ILM(J,3)
                  if (LMSP(ICELL,LM3).GT.0) then
                     IFUN = IFUNM(ICELL,LM3)
                     do I = IRS1 + 1,IRC1
                        V1(I) = V1(I) + GSH(J)*RHO2NS(I,LM2,IATYP,1)*THETAS(I-IRS1,IFUN,ICELL)
                     enddo ! I
                  end if
               enddo ! J

               do I = IRS1 + 1,IRC1
                  RL = R(I,IATYP)**L
                  V2(I) = V1(I)/R(I,IATYP)/RL*DRDI(I,IATYP)
                  V1(I) = V1(I)*RL*DRDI(I,IATYP)
               enddo ! I
            end if
            !-------------------------------------------------------------------
            ! Now integrate v1 and v2
            !-------------------------------------------------------------------
            call SOUTK(V1,VINT1,IPAN(IATYP),IRCUTM)
            call SINWK(V2,VINT2,IPAN(IATYP),IRCUTM)
            !-------------------------------------------------------------------
            ! Gather all parts
            !-------------------------------------------------------------------
            if (LM.EQ.1) then
               V(1,LM,IPOT) = FAC*VINT2(1)
            else
               V(1,LM,IPOT) = 0.0D0
            end if

            do I = 2,IRC1
               RL = R(I,IATYP)**L
               V(I,LM,IPOT) = FAC* (VINT1(I)/R(I,IATYP)/RL+VINT2(I)*RL)
            enddo ! I
            !-------------------------------------------------------------------
            ! Store charge moment - in case of kshape.gt.0 this is the moment
            !      of the charge in the muffin tin sphere
            !-------------------------------------------------------------------
            CMOM(LM,IATYP) = VINT1(IRS1)
            !-------------------------------------------------------------------
            ! Store charge moment of interstial in case of kshape.gt.0
            !-------------------------------------------------------------------
            if (KSHAPE.NE.0) CMINST(LM,IATYP) = VINT1(IRC1) -VINT1(IRS1)
            !
            if (NSPIN.EQ.2) then
               do I = 1,IRC1
                  V(I,LM,IPOT-1) = V(I,LM,IPOT)
               enddo ! I
            end if
         enddo ! M
      enddo ! L
   enddo ! IATYP
   return

end subroutine VINTRAS
