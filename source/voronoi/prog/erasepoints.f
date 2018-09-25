      SUBROUTINE ERASEPOINTS(npolypoi,keepit,poi2plane,poi2poi,
     &     polypoints)
      implicit none
c#@# KKRtags: VORONOI    
      integer npoimax,nneimax,nplanemax
      parameter (npoimax=300,nneimax=100,nplanemax=100)
      integer npolypoi
      integer poi2poi(0:nneimax,npoimax),poi2plane(0:nplanemax,npoimax)
      real*8 polypoints(3,npoimax)
      logical keepit(npoimax)
      integer nnew1,n,i,ip,i0
      integer t_poi2plane(0:nplanemax),t_poi2poi(0:nneimax)
      real*8 tr(3)
c 
      nnew1 = npolypoi 
c
      n = 0 
      i0 = 0 
 
      do ip=1,npolypoi
        ! write(6,*) 'ip = ',ip
         if (.not.keepit(ip)) then
            ! write(6,*) 'Erasing point ',ip
            i0 = i0 + 1
            ! !!!!!!!! erase(i0) = ip
            nnew1 = nnew1 - 1
            
         else
          !   write(6,*) 'change ',n+1
            n = n + 1
c     copy to tmp          
            t_poi2poi(0) = poi2poi(0,ip)
           ! write(6,*) 'tpoi',t_poi2poi(0)
            do i=1,t_poi2poi(0)
               t_poi2poi(i) = poi2poi(i,ip)
              ! write(6,*) 'tpoi',poi2poi(i,ip)
            end do
            t_poi2plane(0) = poi2plane(0,ip)
            do i=1,t_poi2plane(0)
               t_poi2plane(i) = poi2plane(i,ip)
            end do
            do i=1,3
               tr(i) = polypoints(i,ip)
            end do
c     now map back
            ! write(6,*) 'lalala',n,poi2poi(0,n),poi2plane(0,n)

            poi2poi(0,n) = t_poi2poi(0)
            do i=1,t_poi2poi(0)
               ! write(6,*) 'lalala',t_poi2poi(i)
               poi2poi(i,n) = t_poi2poi(i) 
            end do 
            poi2plane(0,n) = t_poi2plane(0) 
            do i=1,t_poi2plane(0)
              ! write(6,*) 'lalala',t_poi2plane(i)
               poi2plane(i,n) = t_poi2plane(i) 

            end do         
            do i=1,3
               polypoints(i,n) = tr(i)
            end do
            ! write (6,*) 'poloioihh',ip
         end if
         
      end do
c  now we know which have to be errased...
      npolypoi = n
         
      ! write(6,*) 'Finished ',npolypoi              

      end



