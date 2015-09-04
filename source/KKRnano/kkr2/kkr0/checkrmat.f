!*==checkrmat.f    processed by SPAG 6.05Rc at 12:32 on  3 Oct 2021
      logical function checkrmat(rmat,co1,si1,co2,si2,co3,si3,i,j)
!   ********************************************************************
!   *                                                                  *
!   *  check whether the values of the cosinus and sinus found for the *
!   *  euler angles tet1, tet2, tet3 are consistent with the           *
!   *  rotation matrix   rmat                                          *
!   *                                                                  *
!   ********************************************************************
      implicit none
!
!*** start of declarations rewritten by spag
!
! dummy arguments
!
      double precision, intent(in) :: co1,co2,co3,si1,si2,si3
      integer, intent(in) :: i,j
      double precision, intent(in) :: rmat(3,3)
!
! local variables
!
      double precision :: a, b
      logical equal
!
!*** end of declarations rewritten by spag
!
      equal(a,b) = (abs(a-b) < 1.d-7)
!
      checkrmat = .false.
!
      if (i == 1) then
         if (j == 1) then
            checkrmat = equal(rmat(1,1),co3*co2*co1-si3*si1)
         else if (j == 2) then
            checkrmat = equal(rmat(1,2),co3*co2*si1+si3*co1)
         else if (j == 3) then
            checkrmat = equal(rmat(1,3),-co3*si2)
         end if
      else if (i == 2) then
         if (j == 1) then
            checkrmat = equal(rmat(2,1),-si3*co2*co1-co3*si1)
         else if (j == 2) then
            checkrmat = equal(rmat(2,2),-si3*co2*si1+co3*co1)
         else if (j == 3) then
            checkrmat = equal(rmat(2,3),si3*si2)
         end if
      else if (j == 1) then
         checkrmat = equal(rmat(3,1),si2*co1)
      else if (j == 2) then
         checkrmat = equal(rmat(3,2),si2*si1)
      else if (j == 3) then
         checkrmat = equal(rmat(3,3),co2)
      end if
      end
