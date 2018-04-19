!-------------------------------------------------------------------------------
! SUBROUTINE: VLLNS
!> @brief Transformation of the wavefunctions for non spherical potentials.
!> @details To determine the non - spherical wavefunctions the potential
!> has to be lm1 and lm2 dependent . the potential is stored
!> only as lm dependent , therefore a transformation in the
!> following way has to be done :
!> \f$ vnsll(r,lm1,lm2)   =  \sum_{lm3} \left\{  c(lm1,lm2,lm3) *vins(r,lm3)  \right\}\f$
!> where c(lm1,lm2,lm3) are the gaunt coeffients. (see notes by B. Drittler)
!> @author B. Drittler
!> @date July 1988
!> @note attention : The gaunt coeffients are stored in an index array only for lm1.gt.lm2
!> (see subroutine gaunt)
!> - R. Zeller Sep. 2000: modified
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,IRM,NCLEB,LMPOT,IRMIND,LMMAXD)

   implicit none

   ! .. Input variables
   integer, intent(in) :: IRM       !< Maximum number of radial points
   integer, intent(in) :: IEND
   integer, intent(in) :: NCLEB     !< Number of Clebsch-Gordon coefficients
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: IRMIND    !< IRM-IRNSD
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   ! .. Array Arguments
   integer, dimension(NCLEB,4), intent(in) :: ICLEB  !< Pointer array
   double precision, dimension(NCLEB,2), intent(in) :: CLEB   !< GAUNT coefficients (GAUNT)
   double precision, dimension(IRMIND:IRM,LMPOT), intent(in) :: VINS !< Non-spherical part of the potential
   ! .. Output variables
   double precision, dimension(LMMAXD,LMMAXD,IRMIND:IRM), intent(out) :: VNSPLL
   ! .. Local Scalars
   integer :: IR,J,LM1,LM2,LM3
   ! ..
   do LM1 = 1,LMMAXD
      do LM2 = 1,LM1
         do IR = IRMIND,IRM
            VNSPLL(LM1,LM2,IR) = 0.0D0
         enddo
      enddo
   enddo
   !
   do J = 1,IEND
      LM1 = ICLEB(J,1)
      LM2 = ICLEB(J,2)
      LM3 = ICLEB(J,3)
      do  IR = IRMIND,IRM
         VNSPLL(LM1,LM2,IR) = VNSPLL(LM1,LM2,IR) +CLEB(J,1)*VINS(IR,LM3)
      enddo
   enddo
   !----------------------------------------------------------------------------
   ! Use symmetry of the gaunt coef.
   !----------------------------------------------------------------------------
   do LM1 = 1,LMMAXD
      do LM2 = 1,LM1 - 1
         do IR = IRMIND,IRM
            VNSPLL(LM2,LM1,IR) = VNSPLL(LM1,LM2,IR)
         enddo
      enddo
   enddo

   do LM1 = 1,LMMAXD
      do IR = IRMIND,IRM
         VNSPLL(LM1,LM1,IR) = VNSPLL(LM1,LM1,IR) + VINS(IR,1)
      end do
   end do

end subroutine VLLNS
