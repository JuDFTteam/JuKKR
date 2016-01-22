      SUBROUTINE CUTPLANE(npolypoi,npolyplan,polypoints,
     &     poi2poi,poi2plane,
     &     polyplaneA3,polyplaneB3,polyplaneC3,polyplaneD3,
     &     A3,B3,C3,D3,NEDGE,poiofplane)
      implicit none    
      integer npoimax,nneimax,nplanemax
      parameter (npoimax=300,nneimax=100,nplanemax=100)
       
c * This sub takes a polyhedron and the equation of a plane 
c   A3*x+B3*y+C3*z=D3       and produces  new polyhedron cuted 
c   by this plane.
c                                                     v. 4.9.01
      integer npolypoi,npolyplan
      integer poi2poi(0:nneimax,npoimax),poi2plane(0:nplanemax,npoimax)
      integer nedge(*),poiofplane(nplanemax,npoimax)
      real*8 polyplaneA3(nplanemax),polyplaneB3(nplanemax),
     &     polyplaneC3(nplanemax),polyplaneD3(nplanemax)
      real*8 polypoints(3,npoimax)
      real*8 a3,b3,c3,d3
c --------- Local 
      real*8 xcut,ycut,zcut,x1,y1,z1,x2,y2,z2,a1
      integer localplanes(0:nplanemax)
      logical keepit(npoimax),lerrase
      integer ip,nnew,i,isimilar(nneimax),ip0,ip2,num,in,i1,
     &   isimilar1(nneimax),ip1
      real*8 newpoint(3,nneimax)
      integer newplanes(0:nneimax,npoimax),nnew1 
c external 
      logical halfspace
c ===================================================     
c       write(6,*) '**********************************'
c      do ip=1,npolypoi
c        write(6,fmt='(i5,3f8.3)') ip,(polypoints(i,ip),i=1,3)
c      end do  


      DO IP=1,npolypoi
         xcut = polypoints(1,ip)
         ycut = polypoints(2,ip)
         zcut = polypoints(3,ip)
         keepit(ip) =  (HALFSPACE(A3,B3,C3,D3,XCUT,YCUT,ZCUT))
c          write(6,*) 'ip, halfplane ',ip,keepit(ip) 
      end do
      nnew = 0      

      DO i=1,nneimax
         isimilar(i) = 0
         isimilar1(i) = 0
      end do
      
      do ip=1,npolypoi
         if (.not.keepit(ip)) then ! one point is left outside, have to redefine 
                                   ! polyhedron
           ! write(6,*) '%%%%%%%%%%%%%%%%%%%',ip 
            
c     
c     Introduce new points
c     
            x1 = polypoints(1,ip)
            y1 = polypoints(2,ip)
            z1 = polypoints(3,ip)
c           write(6,*) 'Old point: ',ip,'  has ',poi2poi(0,ip),
c     &           ' neigh ',(poi2poi(i,ip),i=1,poi2poi(0,ip))
c     loop in all neighbours of ip
            do ip0 = 1,poi2poi(0,ip)
               ip2 = poi2poi(ip0,ip)
c     
c     Now find crossing of line connecting point ip with all neighbors ip2
c     and plane a3,b3,c3,d3
c     
               if (keepit(ip2)) THEN
                  x2 = polypoints(1,ip2)
                  y2 = polypoints(2,ip2)
                  z2 = polypoints(3,ip2)
                  
                  CALL CROSSPOIPLANE(x1,y1,z1,x2,y2,z2,a3,b3,c3,d3,
     &                 xcut,ycut,zcut,a1)
c                  write(6,888) 
c     &                 xcut,ycut,zcut,a1,ip,ip2
 888          format('new ',3F8.3,' a= ',f5.3,' in line :',2I3)
c     new point introduced, keep it until all points are found
                  nnew = nnew + 1
                  newpoint(1,nnew) = xcut
                  newpoint(2,nnew) = ycut
                  newpoint(3,nnew) = zcut
                  
                  if (abs(a1).lt.1.d-8)  ISIMILAR(nnew) = ip
                  if (abs(a1-1.d0).lt.1.d-8) ISIMILAR1(nnew) = ip2
c     find planes that go through ip,ip2 and keep
                  CALL COMMONPLANE(IP,IP2,LOCALPLANES,NPOLYPOI,
     &                 POI2PLANE,POI2POI)
c                 
c                  write(6,*) 'Common planes ',
c     &                 (localplanes(i),i=1,LOCALPLANES(0)) 
                  newplanes(0,nnew) =  LOCALPLANES(0)
                  do i=1, LOCALPLANES(0)
                     newplanes(i,nnew) = localplanes(i)
                  end do     
c     
               END IF           
            end do   
c     ---------- all neigbours of excluded point are done now decide if you
c     will keep the plane, and if you will keep the excluded point...
            

         end if                 ! if (.NOT.keepit(ip)) THEN
      end do                    ! ip=1,npolypoi
c
c New points are now stored see if we have a cut in the polyhedron
c      
      if (nnew.gt.2) then  
c     add plane in the list of planes and update points
         npolyplan = npolyplan + 1
         if (npolyplan.gt.nplanemax) STOP 'Increase nplanemax '
         polyplaneA3(npolyplan) = A3
         polyplaneB3(npolyplan) = B3
         polyplaneC3(npolyplan) = C3 
         polyplaneD3(npolyplan) = D3

         do in=1,nnew

            ip1 = isimilar(in)
            ip2 = isimilar1(in) ! aa
c if all ok they can never be both zero              
            ip = 0
            if (ip1.ne.0) ip = ip1
            if (ip2.ne.0) ip = ip2
            if (ip1.ne.0.and.ip2.ne.0) STOP ' check ip1,ip2'
            if (ip.ne.0) then
c     point already exists add plane in list of planes and find 
c     neighbors in this plane
c               write(6,*) ' Now updating existing point',ip,in
               poi2plane(0,ip) = poi2plane(0,ip) + 1
              ! num = poi2plane(0,ip) + 1     ! changed now!
                num = poi2plane(0,ip)
               poi2plane(num,ip) = npolyplan
               keepit(ip) = .true. ! change the flag ...             
            else   
c     point is new, just append it to the list 
c     
               

                  npolypoi = npolypoi + 1
                  if (npolypoi.gt.npoimax) STOP 'Increase npoimax'
c      write(6,*) ' Now adding new point '
                  do i=1,3
                     polypoints(i,npolypoi) = newpoint(i,in)
                  end do
                  
                  poi2plane(0,npolypoi)  = newplanes(0,in) + 1
                  poi2plane(1,npolypoi) = npolyplan
                  do i1=2,poi2plane(0,npolypoi)
                     poi2plane(i1,npolypoi)=newplanes(i1-1,in)
                  end do
c     
                  keepit(npolypoi) = .true.

               
            end if
c     
c     Now the plane list is updated, nead to update 
c     the poi2poi array (neighbors list) 
c     
         end do
c         do i=1,npolypoi
c         
c         write(6,*) i,' pl ',(poi2plane(i1,i),i1=1,poi2plane(0,i))
c         end do 
c         write(6,*) ' Points and Planes updated 1:',npolypoi,npolyplan
c now erase points outside polyhedron
         nnew1 = npolypoi
 
         call erasepoints(npolypoi,keepit,poi2plane,poi2poi,polypoints)


c         write(6,*) ' Points updated 2:',npolypoi,npolyplan

         CALL poi2poiupdate(poi2poi,poi2plane,
     &        npolyplan,npolypoi,keepit,
     &        polyplaneA3,polyplaneB3,polyplaneC3,polyplaneD3,NEDGE,
     &        poiofplane)
            
         
      end if                    ! nnew.gt.2
      

      end
      
      











