      SUBROUTINE FINDEDGES(npolypoi,poi2poi,polyedge1,polyedge2,
     &     nedge)
      implicit none
      integer npoimax,nneimax,nplanemax
      parameter (npoimax=300,nneimax=100,nplanemax=100)
      integer polyedge1(1000),polyedge2(1000)
      integer npolypoi,poi2poi(0:nneimax,npoimax)
      integer i,nn,np,j,ii,it,it1,it2,pair1,pair2,nedge
c
      nedge  = 0
      do i=1,npoimax 
         np = poi2poi(0,i)
         do ii=1,np
            pair1 = i
            pair2 = poi2poi(ii,i)
            nn = 0
            do j=1,nedge
              it1 = abs(polyedge1(j)-pair1) + abs(polyedge2(j)-pair2)
              it2 = abs(polyedge1(j)-pair2) + abs(polyedge2(j)-pair1)  
              it = it1*it2
              if (it.ne.0) then
c new edge found
              nn = nn + 1
              polyedge1(nedge + nn) = pair1
              polyedge2(nedge + nn) = pair2                    
              end if
            end do
            nedge = nedge + 1
         end do
      end do
      end      
