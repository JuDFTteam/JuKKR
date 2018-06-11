    Subroutine nrmliz(n, r, rn)
      Use mod_datatypes, Only: dp
!-  normalizes a vector
! ----------------------------------------------------------------------
!i Inputs
!i   n     :number of vectors
!i   r     :vector
!o Outputs:
!o   rn    :normalized vector
! ----------------------------------------------------------------------
      Implicit None
      real(kind=dp), parameter :: eps=1E-14_dp
! Passed parameters:                                                    
      Integer :: n
      Real (Kind=dp) :: r(3, *), rn(3, *)
! Local parameters                                                      
      Integer :: i
      Real (Kind=dp) :: d, d2
! External calls                                                        
      External :: dcopy, dscal

      Call dcopy(3*n, r, 1, rn, 1)
      Do i = 1, n
        d2 = r(1, i)*r(1, i) + r(2, i)*r(2, i) + r(3, i)*r(3, i)
        d = sqrt(d2)
        If (abs(d)<eps) Call dscal(3, 1.E0_dp/d, rn(1,i), 1)
      End Do

    End Subroutine
