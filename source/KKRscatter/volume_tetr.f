      SUBROUTINE VOLUME_TETR(a,b,c,volume)

      implicit none

      integer :: i 
      double precision,intent(in)  :: a(3), b(3), c(3)
      double precision             :: d(3)            
      double precision,intent(out) :: volume

      d(1)=a(2)*b(3)-a(3)*b(2)
      d(2)=a(3)*b(1)-a(1)*b(3)
      d(3)=a(1)*b(2)-a(2)*b(1)

      volume=1/6.d0*abs(c(1)*d(1)+c(2)*d(2)+c(3)*d(3))

      END SUBROUTINE VOLUME_TETR
