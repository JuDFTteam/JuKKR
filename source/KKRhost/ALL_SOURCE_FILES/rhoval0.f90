!-------------------------------------------------------------------------------
! SUBROUTINE: RHOVAL0
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine RHOVAL0(EZ,DRDI,RMESH,IPAN,IRCUT,IRWS,THETAS,DOS0,DOS1,IRM,LMAX)
   !
   use Constants
   use global_variables

   implicit none

   ! .. Input variables
   integer, intent(in) :: IRM       !< Maximum number of radial points
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: IPAN  !< Number of panels in non-MT-region
   integer, intent(in) :: IRWS  !< R point at WS radius
   double complex, intent(in) :: EZ
   integer, dimension(0:IPAND), intent(in) :: IRCUT  !< R points of panel borders
   double precision, dimension(IRM), intent(in) :: DRDI !< Derivative dr/di
   double precision, dimension(IRM), intent(in) :: RMESH
   double precision, dimension(IRID,NFUND), intent(in) :: THETAS  !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   ! .. Output variables
   double complex, intent(out) :: DOS0
   double complex, intent(out) :: DOS1
   ! .. Local Scalars
   integer :: IR,L,L1,IMT1
   integer :: LMAXD1
   double precision :: C0LL
   double complex :: EK,CIEK,DENL
   ! .. Local Arrays ..
   double complex, dimension(0:LMAX+1) :: BESSJW
   double complex, dimension(0:LMAX+1) :: BESSYW
   double complex, dimension(0:LMAX+1) :: HANKWS
   double complex, dimension(IRM,0:LMAX)   :: PZ
   double complex, dimension(IRM,0:LMAX)   :: QZ
   double complex, dimension(IRM,0:LMAX+1)  :: CDEN0
   double complex, dimension(IRM,0:LMAX+1)  :: CDEN1

   ! .. Intrinsic Functions
   intrinsic :: ATAN,DBLE,SQRT
   !
   LMAXD1= LMAX+1
   EK = SQRT(EZ)
   C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
   CIEK=CI*EK
   !
   !----------------------------------------------------------------------------
   do IR = 2,IRWS
      call BESHAN(HANKWS,BESSJW,BESSYW,RMESH(IR)*EK,LMAXD1)
      do L = 0,LMAX
         PZ(IR,L) = BESSJW(L)*RMESH(IR)
         QZ(IR,L) = (BESSYW(L) - CI*BESSJW(L))*RMESH(IR)
      enddo
   enddo
   IMT1=IRCUT(1)
   do L1 = 0,LMAXD1
      CDEN0(1,L1) = CZERO
      CDEN1(1,L1) = CZERO
   end do
   do IR = 2,IRWS
      CDEN0(IR,0) = EK*PZ(IR,0)*QZ(IR,0)
      CDEN1(IR,0) = EK*PZ(IR,0)**2*(0.D0,-1.D0)
      CDEN1(IR,LMAXD1) = CIEK*RMESH(IR)**2
   end do
   do L1 = 1,LMAX
      do IR = 2,IRWS
         CDEN0(IR,L1) = EK*PZ(IR,L1)*QZ(IR,L1)*(L1+L1+1)
         CDEN1(IR,L1) = EK*PZ(IR,L1)**2*(0.D0,-1.D0)*(L1+L1+1)
      end do
   end do
   !
   do L1 = 0,LMAXD1 !LMAXD1
      if (IPAN.GT.1) then
         do IR = IMT1 + 1,IRWS
            CDEN0(IR,L1) = CDEN0(IR,L1)*THETAS(IR-IMT1,1)*C0LL
            CDEN1(IR,L1) = CDEN1(IR,L1)*THETAS(IR-IMT1,1)*C0LL
         end do
      end if
      call CSIMPK(CDEN0(1,L1),DENL,IPAN,IRCUT,DRDI)
      DOS0 = DOS0 + DENL
      call CSIMPK(CDEN1(1,L1),DENL,IPAN,IRCUT,DRDI)
      DOS1 = DOS1 + DENL
   end do
   !
end subroutine RHOVAL0
