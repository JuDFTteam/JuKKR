!-------------------------------------------------------------------------------
!> @brief Generate an angular mesh and spherical harmonics at those
!> mesh points. For an angular integration the weights are generated .
!
!> @author R. Zeller
!> @date Feb. 1996
!> @note - Jonathan Chico: Rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine SPHERE_NOGGA(LMAX,YR,WTYR,RIJ,IJD)

   use Constants
   ! ..
   ! .. Scalar Arguments
   integer, intent(in) :: IJD
   integer, intent(in) :: LMAX !< Maximum l component in wave function expansion
   ! .. Output variables
   double precision, dimension(IJD,*), intent(out) :: YR
   double precision, dimension(IJD,3), intent(out) :: RIJ
   double precision, dimension(IJD,*), intent(out) :: WTYR
   ! .. Local variables
   integer :: IJ,LM1
   double precision :: WGHT
   double precision :: R,R1,R2,R3
   double precision, dimension(1000) :: Y
   ! .. External Subroutines
   external :: YMY
   !
   write (1337,*) ' SPHERE : read LEBEDEV mesh'
   if (IJD.GT.1000) stop ' SPHERE '
   !
   do IJ = 1,IJD
      call LEBEDEV(IJ,R1,R2,R3,WGHT)
      RIJ(IJ,1) = R1
      RIJ(IJ,2) = R2
      RIJ(IJ,3) = R3
      call YMY(R1,R2,R3,R,Y,LMAX)
      do LM1 = 1, (LMAX+1)**2
         YR(IJ,LM1) = Y(LM1)
      enddo ! LM1
      !-------------------------------------------------------------------------
      ! Multiply the spherical harmonics with the weights
      !-------------------------------------------------------------------------
      do LM1 = 1, (LMAX+1)**2
         WTYR(IJ,LM1) = YR(IJ,LM1)*WGHT*PI*4.D0
      enddo ! LM1
   enddo ! IJ

end subroutine SPHERE_NOGGA
