!-------------------------------------------------------------------------------
! SUBROUTINE: SINWK
!> @brief This subroutine does an outwards integration of a function with kinks
!> @detail To integrate the function \f$ fint\left(r\right)=\int_{r}^{r_c}f\left(r'\right)dr' \f$
!> at each kink the integration is restarted the starting value for this integration is determined by
!> a 4 point lagrangian integration, coefficients given by M. Abramowitz and I.A. Stegun,
!> handbook of mathematical functions, nbs applied mathematics series 55 (1968).
!> the weights \f$drdi\f$ have to be multiplied before calling this subroutine.
!> @author B. Drittler
!> @date Oct. 1989
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
SUBROUTINE SINWK(F,FINT,IPAN,IRCUT)

   use global_variables
   ! .. Scalar Arguments
   integer, intent(in) :: IPAN   !< Number of panels in non-MT-region
   ! .. Array Arguments
   INTEGER, dimension(0:IPAND), intent(in) :: IRCUT   !< R points of panel borders
   double precision, dimension(*), intent(in) :: F
   ! .. Output variables
   double precision, dimension(*), intent(out) :: FINT
   ! .. Local Scalars
   integer :: I,IEN,IP,IST
   double precision :: A1,A2

   A1 = 1.0D0/3.0D0
   A2 = 4.0D0/3.0D0
   !----------------------------------------------------------------------------
   ! Loop over kinks
   !----------------------------------------------------------------------------
   do IP = IPAN,1,-1
      IST = IRCUT(IP)
      IEN = IRCUT(IP-1) + 1

      if (IP.EQ.IPAN) then
         FINT(IST) = 0.0D0
         !----------------------------------------------------------------------
         ! Integrate fint(ist-1) with a 4 point lagrangian
         !----------------------------------------------------------------------
         FINT(IST-1) = (F(IST-3)-5.0D0*F(IST-2)+19.0D0*F(IST-1)+9.0D0*F(IST))/24.0D0

      else
         FINT(IST) = FINT(IST+1)
         !----------------------------------------------------------------------
         ! Integrate fint(ist-1) with a 4 point lagrangian
         !----------------------------------------------------------------------
         FINT(IST-1) = FINT(IST+1) + (F(IST-3)-5.0D0*F(IST-2)+19.0D0*F(IST-1)+9.0D0*F(IST))/24.0D0
      end if
      !-------------------------------------------------------------------------
      ! Calculate fint with an extended 3-point-simpson
      !-------------------------------------------------------------------------
      do I = IST - 2,IEN,-1
         FINT(I) = ((FINT(I+2)+F(I+2)*A1)+F(I+1)*A2) + F(I)*A1
      enddo ! I
   enddo ! IP

end subroutine SINWK
