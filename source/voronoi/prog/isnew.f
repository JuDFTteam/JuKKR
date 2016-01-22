      LOGICAL FUNCTION ISNEW(x,y,z,rp,np,isame)
      implicit none
      real*8 x,y,z,rp(3,*),diff
      integer np,i,isame
      the_same=.false.
      do i=1,np
         diff = sqrt((rp(1,i)-x)**2+(rp(2,i)-y)**2+(rp(3,i)-z)**2)
         if (diff.lt.1.d-4) THEN 
         THE_SAME=.true.
         isame = i
         end if
      end do
      end
