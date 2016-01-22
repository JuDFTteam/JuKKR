       SUBROUTINE COMMONPLANE(IP1,IP2,LOCALPLANES,NPOINTS,POI2PLANE,
     &                        POI2POI)
       implicit none
       integer npoimax,nneimax,nplanemax
       parameter (npoimax=300,nneimax=100,nplanemax=100)
c
c This looks at the line connecting 2 points and returns the index 
c of common planes
c
      integer ip1,ip2,localplanes(0:nplanemax),npoints,
     &      poi2plane(0:nplanemax,npoimax),poi2poi(0:nneimax,npoimax)
c local
      integer i,i1,i2,ic,nl,index(nplanemax)
c     ---------
      nl = 0
      DO I1 = 1,POI2PLANE(0,IP1)

         DO I2 = 1,POI2PLANE(0,ip2)

            ic = poi2plane(i1,ip1) - poi2plane(i2,ip2)
c       write(6,*) 'comm',i1,i2,poi2plane(i1,ip1),poi2plane(i2,ip2)
            if (ic.eq.0) then
c     common plane found
               nl = nl + 1
               index(nl) = poi2plane(i1,ip1)
            end if
         end do
      end do

c      if (nl.ne.2.and.nl.ne.0) write(6,*) 
c     &     'COMMONPLANE found more than 2 common planes',nl

      LOCALPLANES(0) = nl

      do i=1,nl
         localplanes(i) = index(i) 
      end do
      end 
