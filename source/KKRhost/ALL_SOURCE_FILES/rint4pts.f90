SUBROUTINE rint4pts(y,jtop,z)
!   ********************************************************************
!   *                                                                  *
!   *      perform the integral  Z(i)   =  INT   Y(i') di'             *
!   *                                    R=0..R(i)                     *
!   *                                                                  *
!   *      via a 4-point integration formula                           *
!   *                                                                  *
!   *      JTOP:     Y is tabulated form 1 .. JTOP                     *
!   *      Y(i):     function to be integrated                         *
!   *                                                                  *
!   *                       REAL    - VERSION                          *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER JTOP
REAL*8 Y(JTOP),Z(JTOP)
! Local variables
INTEGER I,IG,J,K,M,N1,N2
REAL*8 Q(5,5),Q5(5,5),S,SVN
DATA Q5/0.D0,251.D0,232.D0,243.D0,224.D0,0.D0,646.D0,992.D0, &
     918.D0,1024.D0,0.D0, - 264.D0,192.D0,648.D0,384.D0,0.D0, &
     106.D0,32.D0,378.D0,1024.D0,0.D0, - 19.D0, - 8.D0, - 27.D0, &
     224.D0/

DO i = 1,5
  DO j = 1,5
    q(i,j) = q5(i,j)/720.0D0
  END DO
END DO

z(1) = 0.0D0
svn = z(1)

DO ig = 1,jtop - 4,4
  n1 = ig
  n2 = ig + 4
  DO m = n1 + 1,n2
    i = m - n1 + 1
    s = svn
    DO k = n1,n2
      j = k - n1 + 1
      s = s + q(i,j)*y(k)
    END DO
    z(m) = s
  END DO
  svn = z(n2)
END DO

IF ( n2 /= jtop ) THEN
  n1 = jtop - 4
  n2 = jtop
  svn = z(n1)
  DO m = n1 + 1,n2
    i = m - n1 + 1
    s = svn
    DO k = n1,n2
      j = k - n1 + 1
      s = s + q(i,j)*y(k)
    END DO
    z(m) = s
  END DO
END IF

END SUBROUTINE rint4pts
