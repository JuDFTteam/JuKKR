    Subroutine dinv33(matrix, iopt, invers, det)
      Use mod_datatypes, Only: dp
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
      Implicit None
! Passed parameters:                                                    
      Integer :: iopt
      Real (Kind=dp) :: matrix(3, 3), invers(3, 3), det
! Local parameters:                                                     
      Integer :: i, j
      Real (Kind=dp) :: ddot1, twopi, pi
      Parameter (pi=3.141592653589793E0_dp)
      Parameter (twopi=2.E0_dp*pi)

! external calls:                                                       
      External :: cross, ddot1, dscal1, dswap1

      Call cross(matrix(1,2), matrix(1,3), invers(1,1))
      Call cross(matrix(1,3), matrix(1,1), invers(1,2))
      Call cross(matrix(1,1), matrix(1,2), invers(1,3))
      det = ddot1(3, matrix, 1, invers, 1)
      If (iopt>=2) det = det/twopi
      If (mod(iopt,2)==0) Then
        Do i = 1, 3
          Do j = i + 1, 3
            Call dswap1(1, invers(i,j), 1, invers(j,i), 1)
          End Do
        End Do
      End If
      Call dscal1(9, 1.E0_dp/det, invers, 1)
    End Subroutine
