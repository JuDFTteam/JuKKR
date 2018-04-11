!-------------------------------------------------------------------------------
! SUBROUTINE: SINWK
!> @brief This subroutine does an integration up to \f$ r_{cut}\f$ of an real function
!> \f$f\f$ with an extended 3-point-simpson
!> @details The integration of the function
!> \f$ fint =\int_{0}^{r_{cut}} f\left(r'\right)dr'\f$ has been modified
!> for functions with kinks - at each kink the integration is restarted.
!> @note Attention : Input \f$f\f$ is destroyed !
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine SIMPK(F,FINT,IPAN,IRCUT,DRDI)

   use global_variables

   integer, intent(in) :: IPAN !< Number of panels in non-MT-region
   integer, dimension(0:IPAND), intent(in) :: IRCUT !< R points of panel borders
   double precision, dimension(*), intent(in) :: DRDI !< Derivative dr/di
   ! .. Output variables
   double precision, intent(out) :: FINT
   ! .. In/Out variables
   double precision, dimension(*), intent(inout) :: F
   ! .. Local Scalars
   integer :: I,IEN,IP,IST,N
   double precision :: A1,A2

   ! .. External Functions
   double precision :: SSUM
   external :: SSUM
   !     ..
   A1 = 4.0D0/3.0D0
   A2 = 2.0D0/3.0D0
   FINT = 0.0D0
   !
   do IP = 1,IPAN
      !-------------------------------------------------------------------------
      ! Loop over kinks
      !-------------------------------------------------------------------------
      IST = IRCUT(IP-1) + 1
      IEN = IRCUT(IP)
      !
      do I = IST,IEN
         F(I) = F(I)*DRDI(I)
      enddo ! I
      if (MOD(IEN-IST,2).EQ.0) then
         FINT = FINT + (F(IST)-F(IEN))/3.0D0
         IST = IST + 1
         N = (IEN-IST+1)/2
      else
         !----------------------------------------------------------------------
         ! Four point Lagrange integration for the first step
         !----------------------------------------------------------------------
         FINT = FINT + (9.0D0*F(IST)+19.0D0*F(IST+1)-5.0D0*F(IST+2)+ &
            F(IST+3))/24.0D0 + (F(IST+1)-F(IEN))/3.0D0
         IST = IST + 2
         N = (IEN-IST+1)/2
      end if
      !-------------------------------------------------------------------------
      ! Calculate with an extended 3-point-simpson
      !-------------------------------------------------------------------------
      FINT = FINT + A1*SSUM(N,F(IST),2) + A2*SSUM(N,F(IST+1),2)
   enddo ! IP
end subroutine SIMPK
