!-------------------------------------------------------------------------------
! SUBROUTINE: VXCDRV
!> @brief Wrapper for the calculation of the exchange-correlation energy, for
!> the different treatments of the xc-potential.
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,NSTART,NEND,RHO2NS,VONS, &
      R,DRDI,A,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,ILM, &
      IMAXSH,IFUNM,THETAS,LMSP,NPOTD,LMMAX,LMPOT,NGSHD,NFUND,LMXSPD)

   implicit none

   ! .. Input variables
   integer, intent(in) :: KTE    !< Calculation of the total energy On/Off (1/0)
   integer, intent(in) :: KXC    !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
   integer, intent(in) :: NEND
   integer, intent(in) :: LPOT   !< Maximum l component in potential expansion
   integer, intent(in) :: NPOTD  !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
   integer, intent(in) :: LMMAX  !< (LMAX+1)^2
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: NGSHD  !< Shape functions parameters in non-spherical part
   integer, intent(in) :: NFUND
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: NSTART
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   integer, intent(in) :: LMXSPD !< (2*LPOT+1)**2
   ! .. Array Arguments ..
   double precision, dimension(NATYP), intent(in)  :: A   !< Constants for exponential R mesh
   double precision, dimension(NGSHD), intent(in)  :: GSH
   double precision, dimension(IRM,NATYP), intent(in)       :: R !< Radial mesh ( in units a Bohr)
   double precision, dimension(0:LPOT,NATYP), intent(inout) :: EXC   !< exchange correlation energy
   double precision, dimension(IRM,NATYP), intent(in)       :: DRDI  !< Derivative dr/di
   double precision, dimension(IRID,NFUND,NCELLD), intent(in) :: THETAS !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   double precision, dimension(IRM,LMPOT,NATYP,2), intent(in) :: RHO2NS !< radial density
   double precision, dimension(IRMD,LMPOT,NPOTD), intent(out) :: VONS !< output potential (nonspherical VONS)
   integer, dimension(NATYP), intent(in)           :: IRWS   !< R point at WS radius
   integer, dimension(NATYP), intent(in)           :: IPAN   !< Number of panels in non-MT-region
   integer, dimension(NATYP), intent(in)           :: NTCELL !< index for WS cell
   integer, dimension(0:LMPOT), intent(in)         :: IMAXSH
   integer, dimension(NGSHD,3), intent(in)         :: ILM
   integer, dimension(NATYP,LMXSPD), intent(in)    :: LMSP   !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(NATYP,LMXSPD), intent(in)    :: IFUNM
   integer, dimension(0:IPAND,NATYP), intent(in)   :: IRCUT  !< r points of panel borders
   ! ..
   ! .. External Subroutines
   external DCOPY,SPHERE_GGA,SPHERE_NOGGA,VXCGGA,VXCLM
   ! .. Local Scalars
   integer :: IATYP,ICELL,IPOT,LMX1
   integer :: IJD
   ! .. Parameters
   parameter (IJD = 434)
   ! .. Local Arrays ..
   integer, dimension(LMXSPD) :: LMSPIAT
   integer, dimension(LMXSPD) :: IFUNMIAT
   double precision, dimension(IJD)          :: THET
   double precision, dimension(IJD,LMPOT)    :: YR
   double precision, dimension(IJD,3)        :: RIJ
   double precision, dimension(IJD,LMPOT)    :: YLM
   double precision, dimension(IJD,LMPOT)    :: WTYR
   double precision, dimension(IJD,LMPOT)    :: DYLMF1
   double precision, dimension(IJD,LMPOT)    :: DYLMF2
   double precision, dimension(IJD,LMPOT)    :: DYLMT1
   double precision, dimension(IJD,LMPOT)    :: DYLMT2
   double precision, dimension(IJD,LMPOT)    :: DYLMTF
   double precision, dimension(IRM,LMPOT,2)  :: RHO2IAT
   !     ..
   if (KXC.lt.3) then
      call SPHERE_NOGGA(LPOT,YR,WTYR,RIJ,IJD)
   else
      call SPHERE_GGA(LPOT,YR,WTYR,RIJ,IJD,LMPOTD,THET,YLM,DYLMT1, &
         DYLMT2,DYLMF1,DYLMF2,DYLMTF)
   end if
   do IATYP = NSTART,NEND
      ICELL = NTCELL(IATYP)
      IPOT = NSPIN* (IATYP-1) + 1
      do LMX1 = 1,LMXSPD
         IFUNMIAT(LMX1) = IFUNM(ICELL,LMX1)
         LMSPIAT(LMX1) = LMSP(ICELL,LMX1)
      end do
      call DCOPY(IRMD*LMPOTD,RHO2NS(1,1,IATYP,1),1,RHO2IAT(1,1,1),1)
      if (NSPIN.eq.2 .or. KREL.eq.1) then
         call DCOPY(IRMD*LMPOTD,RHO2NS(1,1,IATYP,2),1,RHO2IAT(1,1,2),1)
      end if
      if (KXC.LT.3) then
         call VXCLM(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2IAT, &
            VONS(1,1,IPOT),R(1,IATYP),DRDI(1,IATYP), &
            IRWS(IATYP),IRCUT(0,IATYP),IPAN(IATYP), &
            KSHAPE,GSH,ILM,IMAXSH,IFUNMIAT,THETAS(1,1,ICELL), &
            YR,WTYR,IJD,LMSPIAT,NGSHD,LMPOT,NFUND,LMXSPD)
      else
         !----------------------------------------------------------------------
         ! GGA EX-COR POTENTIAL
         !----------------------------------------------------------------------
         call VXCGGA(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2IAT, &
            VONS(1,1,IPOT),R(1,IATYP),DRDI(1,IATYP),A(IATYP),&
            IRWS(IATYP),IRCUT(0,IATYP),IPAN(IATYP),&
            KSHAPE,GSH,ILM,IMAXSH,IFUNMIAT,THETAS(1,1,ICELL),&
            WTYR,IJD,LMSPIAT,THET,YLM,DYLMT1,DYLMT2,&
            DYLMF1,DYLMF2,DYLMTF,NGSHD,LMPOT,NFUND,LMXSPD)
      end if
   end do
end subroutine VXCDRV
