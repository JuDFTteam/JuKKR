!-------------------------------------------------------------------------------
! SUBROUTINE: VLLMAT
!> @brief
!-------------------------------------------------------------------------------
subroutine VLLMAT(IRMIN,NRMAXD,IRC,LMMAX,LMMAXSO,VNSPLL0,VINS,LMPOT, &
      CLEB,ICLEB,IEND,NSPIN,Z,RNEW,USE_SRATRICK,NCLEB)

   implicit none

   !.. Input variables
   integer, intent(in) :: IRC    !< r point for potential cutting
   integer, intent(in) :: IEND
   integer, intent(in) :: NCLEB  !< Number of Clebsch-Gordon coefficients
   integer, intent(in) :: IRMIN  !< max r for spherical treatment
   integer, intent(in) :: LMMAX  !< (LMAX+1)^2
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: NRMAXD !< NTOTD*(NCHEBD+1)
   integer, intent(in) :: LMMAXSO
   integer, intent(in) :: USE_SRATRICK
   double precision, intent(in) :: Z

   integer, dimension(NCLEB,4), intent(in) :: ICLEB
   double precision, dimension(*), intent(in) :: CLEB !< GAUNT coefficients (GAUNT)
   double precision, dimension(IRMIN:IRC,LMPOT,NSPIN), intent(in) :: VINS  !< Non-spherical part of the potential
   double precision, dimension(IRMIN:NRMAXD), intent(in) ::  RNEW
   double complex, dimension(LMMAXSO,LMMAXSO,IRMIN:IRC), intent(out) :: VNSPLL0

   !.. Local variables
   integer :: ISP
   integer :: I,IR,J,LM1,LM2,LM3
   double precision, dimension(LMMAX,LMMAX,IRMIN:IRC,NSPIN) ::  VNSPLL

   do ISP=1,NSPIN
      do LM1 = 1,LMMAX
         do LM2 = 1,LM1
            do IR = IRMIN,IRC
               VNSPLL(LM1,LM2,IR,ISP) = 0.0D0
            enddo ! IR
         enddo ! LM2
      enddo ! LM11

      do J = 1,IEND
         LM1 = ICLEB(J,1)
         LM2 = ICLEB(J,2)
         LM3 = ICLEB(J,3)
         do I = IRMIN,IRC
            VNSPLL(LM1,LM2,I,ISP) = VNSPLL(LM1,LM2,I,ISP) +CLEB(J)*VINS(I,LM3,ISP)
         enddo ! I
      enddo ! J
      !-------------------------------------------------------------------------
      ! Use symmetry of the gaunt coef.
      !-------------------------------------------------------------------------
      do LM1 = 1,LMMAX
         do LM2 = 1,LM1 - 1
            do I = IRMIN,IRC
               VNSPLL(LM2,LM1,I,ISP) = VNSPLL(LM1,LM2,I,ISP)
            enddo ! I
         enddo ! LM2
      enddo ! LM1

      if (USE_SRATRICK.EQ.0) then
         do LM1=1,LMMAX
            do I=IRMIN,IRC
               VNSPLL(LM1,LM1,I,ISP)=VNSPLL(LM1,LM1,I,ISP)+VINS(I,1,ISP)-2d0*Z/RNEW(I)
            enddo
         enddo
      endif

   end do !NSPIN

   ! Set vnspll as twice as large
   VNSPLL0(1:LMMAX,1:LMMAX,IRMIN:IRC)=&
      cmplx(VNSPLL(1:LMMAX,1:LMMAX,IRMIN:IRC,1),0d0)

   if(NSPIN==2)then! hack to make routine work for Bxc-field
      VNSPLL0(LMMAX+1:LMMAXSO,LMMAX+1:LMMAXSO,IRMIN:IRC)=&
         cmplx(VNSPLL(1:LMMAX,1:LMMAX,IRMIN:IRC,NSPIN),0d0)
   end if

end subroutine VLLMAT
