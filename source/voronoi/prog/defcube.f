      SUBROUTINE defcube(alfa,npoints,nplanes,polypoints,poi2poi,
     &     poi2plane,
     &     polyplaneA3,polyplaneB3,polyplaneC3,polyplaneD3)
      implicit none
      integer npoimax,nneimax,nplanemax
      parameter (npoimax=300,nneimax=100,nplanemax=100)
c
      real*8 polypoints(3,npoimax)
      real*8 polyplaneA3(nplanemax),polyplaneB3(nplanemax),
     &     polyplaneC3(nplanemax),polyplaneD3(nplanemax)
      integer poi2poi(0:nneimax,npoimax),
     &     poi2plane(0:nplanemax,npoimax)
      integer npoints,nplanes 
      real*8 alfa
c local
      integer ipoint1,i,iplane
      real*8 ao2
c
c This just defines all arrays for a cube of edge alfa
c
c Initialize added on 26.11.2001
      do ipoint1=0,nneimax
         do i=1,npoimax
            poi2poi(ipoint1,i)  = 0
         end do
      end do
      do ipoint1=0,nplanemax
         do i=1,npoimax
            poi2plane(ipoint1,i)  = 0
         end do
      end do
c added on 26.11.2001

      do ipoint1=1,npoimax
         do i=1,3
            polypoints(i,ipoint1) = 0.d0       
         end do
      end do
      npoints = 8
      ao2 = alfa/2.d0
      polypoints(1,1) =  ao2
      polypoints(2,1) =  ao2
      polypoints(3,1) =  ao2
c  
      polypoints(1,2) =  ao2
      polypoints(2,2) = -ao2
      polypoints(3,2) =  ao2
c
      polypoints(1,3) = -ao2
      polypoints(2,3) = -ao2
      polypoints(3,3) =  ao2
c  
      polypoints(1,4) = -ao2
      polypoints(2,4) =  ao2
      polypoints(3,4) =  ao2
c
      polypoints(1,5) =  ao2
      polypoints(2,5) =  ao2
      polypoints(3,5) = -ao2
c  
      polypoints(1,6) =  ao2
      polypoints(2,6) = -ao2
      polypoints(3,6) = -ao2
c
      polypoints(1,7) = -ao2
      polypoints(2,7) = -ao2
      polypoints(3,7) = -ao2
c  
      polypoints(1,8) = -ao2
      polypoints(2,8) =  ao2
      polypoints(3,8) = -ao2 
c ************************************************      
      poi2poi(0,1) = 3 ! number of neighbors of point 1 
      poi2poi(1,1) = 2 ! point indeces
      poi2poi(2,1) = 4
      poi2poi(3,1) = 5
      poi2plane(0,1) = 3  ! number of planes passing from point 1
      poi2plane(1,1) = 1  ! plane indeces
      poi2plane(2,1) = 3
      poi2plane(3,1) = 5 
c
      poi2poi(0,2) = 3 ! number of neighbors of point 
      poi2poi(1,2) = 1  ! point indeces
      poi2poi(2,2) = 3
      poi2poi(3,2) = 6
      poi2plane(0,2) = 3  ! number of planes passing from point 
      poi2plane(1,2) = 1   ! plane indeces
      poi2plane(2,2) = 4
      poi2plane(3,2) = 5         
c
      poi2poi(0,3) = 3 ! number of neighbors of point 
      poi2poi(1,3) = 2  ! point indeces
      poi2poi(2,3) = 4
      poi2poi(3,3) = 7
      poi2plane(0,3) = 3  ! number of planes passing from point 
      poi2plane(1,3) = 2   ! plane indeces
      poi2plane(2,3) = 4
      poi2plane(3,3) = 5 
c
      poi2poi(0,4) = 3 ! number of neighbors of point 
      poi2poi(1,4) = 1  ! point indeces
      poi2poi(2,4) = 3
      poi2poi(3,4) = 8
      poi2plane(0,4) = 3  ! number of planes passing from point 
      poi2plane(1,4) = 2    ! plane indeces
      poi2plane(2,4) = 3 
      poi2plane(3,4) = 5
c
      poi2poi(0,5) = 3 ! number of neighbors of point 
      poi2poi(1,5) = 1  ! point indeces
      poi2poi(2,5) = 6
      poi2poi(3,5) = 8
      poi2plane(0,5) = 3  ! number of planes passing from point 
      poi2plane(1,5) = 1   ! plane indeces
      poi2plane(2,5) = 3
      poi2plane(3,5) = 6
c
      poi2poi(0,6) = 3 ! number of neighbors of point 
      poi2poi(1,6) = 2  ! point indeces
      poi2poi(2,6) = 5
      poi2poi(3,6) = 7
      poi2plane(0,6) = 3  ! number of planes passing from point 
      poi2plane(1,6) = 1   ! plane indeces
      poi2plane(2,6) = 4
      poi2plane(3,6) = 6
c
      poi2poi(0,7) = 3 ! number of neighbors of point 
      poi2poi(1,7) = 3  ! point indeces
      poi2poi(2,7) = 6
      poi2poi(3,7) = 8
      poi2plane(0,7) = 3  ! number of planes passing from point 
      poi2plane(1,7) = 2   ! plane indeces
      poi2plane(2,7) = 4
      poi2plane(3,7) = 6
c
      poi2poi(0,8) = 3 ! number of neighbors of point 
      poi2poi(1,8) = 4  ! point indeces
      poi2poi(2,8) = 5
      poi2poi(3,8) = 7
      poi2plane(0,8) = 3  ! number of planes passing from point 
      poi2plane(1,8) = 2   ! plane indeces
      poi2plane(2,8) = 3
      poi2plane(3,8) = 6
c ---------------------------
      nplanes = 6
      do iplane=1,nplanes
         polyplaneA3(iplane) = 0.d0
         polyplaneB3(iplane) = 0.d0
         polyplaneC3(iplane) = 0.d0
         polyplaneD3(iplane) = 0.d0
      end do
      polyplaneA3(1) = 1.d0
      polyplaneD3(1) = alfa/2.d0

      polyplaneA3(2) = 1.d0
      polyplaneD3(2) = -alfa/2.d0

      polyplaneB3(3) = 1.d0
      polyplaneD3(3) = alfa/2.d0

      polyplaneB3(4) = 1.d0
      polyplaneD3(4) = -alfa/2.d0
      
      polyplaneC3(5) = 1.d0
      polyplaneD3(5) = alfa/2.d0

      polyplaneC3(6) = 1.d0
      polyplaneD3(6) = -alfa/2.d0
c
c
      return 
      end       











