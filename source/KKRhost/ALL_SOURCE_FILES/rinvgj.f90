SUBROUTINE rinvgj(ainv,a,arraydim,n)
!   ********************************************************************
!   *                                                                  *
!   *                      AINV = A**(-1)                              *
!   *                                                                  *
!   *  invert A using the GAUSS-JORDAN - algorithm                     *
!   *  the 1- matrix is not set up and use is made of its structure    *
!   *                                                                  *
!   *                    REAL*8 VERSION                                *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER ARRAYDIM,N
REAL*8 A(ARRAYDIM,ARRAYDIM),AINV(ARRAYDIM,ARRAYDIM)

! Local variables
INTEGER ICOL,L,LL
REAL*8 T,T1

ainv(1,1) = 0D0
!                                                        scan columns
DO icol = 1,n
  
!                                               make A(ICOL,ICOL) = 1
  t1 = 1.0D0/a(icol,icol)
  DO l = (icol+1),n
    a(icol,l) = a(icol,l)*t1
  END DO
  
  DO l = 1,(icol-1)
    ainv(icol,l) = ainv(icol,l)*t1
  END DO
  ainv(icol,icol) = t1
  
!                                    make A(LL,ICOL) = 0 for LL<>ICOL
  DO ll = 1,n
    IF ( ll /= icol ) THEN
      t = a(ll,icol)
      DO l = (icol+1),n
        a(ll,l) = a(ll,l) - a(icol,l)*t
      END DO
      
      DO l = 1,(icol-1)
        ainv(ll,l) = ainv(ll,l) - ainv(icol,l)*t
      END DO
      ainv(ll,icol) = -t1*t
    END IF
  END DO
END DO

END SUBROUTINE rinvgj
