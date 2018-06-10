SUBROUTINE dinv33(matrix,iopt,invers,det)
!- Inverts 3X3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   matrix:input matrix
!i   iopt  :if 0, usual inverse
!i             1, transpose of inverse
!i             2, 2*pi*inverse
!i             3, 2*pi*transpose of inverse
!o Outputs:
!o   invers:as modified according to iopt
!o   det   :determinant, (or det/2*pi if iopt=2,3)
!r Remarks:
!r  To generate reciprocal lattice vectors, call dinv33(plat,3,plat)
! ----------------------------------------------------------------------
      implicit none 
! Passed parameters:                                                    
      integer iopt 
      double precision matrix(3,3),invers(3,3),det 
! Local parameters:                                                     
      integer i,j 
      double precision ddot1,twopi,pi
      parameter (pi=3.141592653589793d0)
      parameter(twopi=2.d0*pi)

! external calls:                                                       
      external cross,ddot1,dscal1,dswap1 
                                                                        
CALL cross(matrix(1,2),matrix(1,3),invers(1,1))
CALL cross(matrix(1,3),matrix(1,1),invers(1,2))
CALL cross(matrix(1,1),matrix(1,2),invers(1,3))
det = ddot1(3,matrix,1,invers,1)
IF (iopt >= 2) det = det/twopi
IF (MOD(iopt,2) == 0) THEN
  DO i = 1, 3
    DO j = i+1, 3
      CALL dswap1(1,invers(i,j),1,invers(j,i),1)
    END DO
  END DO
END IF
CALL dscal1(9,1.d0/det,invers,1)
END SUBROUTINE dinv33
