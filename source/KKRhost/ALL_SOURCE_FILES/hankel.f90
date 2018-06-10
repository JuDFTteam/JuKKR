SUBROUTINE hankel(h,l,arg)
!  this subroutine uses the explicit formulas for the hankel
!  functions. for higher l-values these formulas may lead to
!  loss of significant figures. This subroutine should be used
!  only for core states.
implicit none
!.. Scalar Arguments ..
      DOUBLE COMPLEX ARG
      INTEGER L
!..
!.. Array Arguments ..
      DOUBLE COMPLEX H(*)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX A1,A2,A3,A4
!..
!.. Intrinsic Functions ..
      INTRINSIC EXP
!..
!.. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
!     ..
h(1) = -EXP(arg*ci)/arg
IF (l /= 1) THEN
  a1 = (1.d0,0.d0) - arg*ci
  h(2) = h(1)*a1/arg
  IF (l /= 2) THEN
    a1 = 3.d0*a1
    a2 = arg*arg
    h(3) = h(1)* (a1-a2)/a2
    IF (l /= 3) THEN
      a1 = 5.d0*a1
      a3 = a2*arg*ci
      a4 = a2*arg
      a2 = 6.d0*a2
      h(4) = h(1)* (a1-a2+a3)/a4
      IF (l /= 4) THEN
        a1 = 7.d0*a1
        a2 = 7.5D0*a2
        a3 = 10.d0*a3
        a4 = a4*arg
        h(5) = h(1)* (a1-a2+a3+a4)/a4
        IF (l /= 5) THEN
          h(6) = (9.0D0,0.0D0)*h(5)/arg - h(4)
          IF (l /= 6) THEN
            WRITE (6,FMT=9000) l
            STOP 'HANKEL'
            
          END IF
          
        END IF
        
      END IF
      
    END IF
    
  END IF
  
END IF

RETURN


9000 FORMAT (2X,' hankel :  l=',i2,' is too large')
END SUBROUTINE hankel
