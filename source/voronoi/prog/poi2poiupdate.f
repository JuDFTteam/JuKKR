      SUBROUTINE POI2POIUPDATE(poi2poi,poi2plane,npolyplan,npolypoi,
     &     keepit,polyplaneA3,polyplaneB3,polyplaneC3,polyplaneD3,
     &     nedge,poiofplane)
      implicit none  
      integer npoimax,nneimax,nplanemax
      parameter (npoimax=300,nneimax=100,nplanemax=100)    
c ***** 
c    This sub updates the neighbours array poi2poi, it checks all atoms
c    some are newly defined, and looks if a pair of atoms belong
c    to the same 2 planes
c
      integer npolyplan,npolypoi
      integer poi2poi(0:nneimax,npoimax),poi2plane(0:nplanemax,npoimax)
      integer poiofplane(nplanemax,npoimax),nedge(*)  
      logical keepit(npoimax)
      real*8 polyplaneA3(nplanemax),polyplaneB3(nplanemax),
     &     polyplaneC3(nplanemax),polyplaneD3(nplanemax)
      real*8 at,bt,ct,dt
c local
      integer ip1,ip2,nnei,localplanes(0:nplanemax),points(npoimax)
      integer iplane,np,index,np1,i1,i
c ---------------
      do ip1=1,npolypoi
            nnei = 0
            do ip2=1,npolypoi
               if (ip1.ne.ip2) then
                  CALL COMMONPLANE(IP1,IP2,LOCALPLANES,NPOLYPOI,
     &                 POI2PLANE,POI2POI)
                  if (localplanes(0).eq.2) then
c     this pair are neighbours, add it to the list of poi2poi
                     nnei = nnei + 1
                     poi2poi(nnei,ip1) = ip2 
                  end if
               end if
         end do
         poi2poi(0,ip1) = nnei
      end do
c     
c     Poi2poi updated now look if all planes have 3 or more atoms in them
c     
      do ip1=1,npoimax
         points(ip1) = 0
      end do
      do ip1=1,npolypoi
         np =  poi2plane(0,ip1)
         do iplane =1,np 
            index = poi2plane(iplane,ip1)
            if (index.gt.0) then
            points(index) = points(index) + 1
            
            poiofplane(index,points(index)) = ip1
            end if
         end do
      end do
c     
c     
      np = 0
      do iplane =1,npolyplan
         i1 = 0
c         write(6,*) 'Plane ',iplane,' has ',points(iplane),' points'
          
         NEDGE(IPLANE) = points(iplane)
c         write(6,*) (poiofplane(iplane,i),i=1,nedge(iplane))
         

      end do      

      end
      






