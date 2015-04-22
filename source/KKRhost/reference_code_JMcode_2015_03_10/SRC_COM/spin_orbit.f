C ************************************************************************
      SUBROUTINE SPIN_ORBIT_ONE_L(lmax,L_S)

      implicit none
C ************************************************************************
c      in this subroutine the matrix L*S is calculated for the basis of
c      real spherical harmonics 

c      schematically it has the form
c      (  -L_z    L_+  )
c      (  L_-     L_z  )
c

       integer, intent(in)         ::    lmax
       double complex, intent(out) ::    L_S((2*lmax+1)*2,(2*lmax+1)*2)

c  local variables 
       integer                     ::    i1,i2,i1l 
       double complex              ::    icompl
       double complex,allocatable  ::    l_min(:,:)
       double complex,allocatable  ::    l_up(:,:)
       double precision            ::    lfac



       icompl=(0d0,1d0) 


       allocate(l_min(-lmax:lmax,-lmax:lmax))
       allocate(l_up(-lmax:lmax,-lmax:lmax))

c  initialize the matrix 

       do i1=1,(2*lmax+1)*2
         do i2=1,(2*lmax+1)*2
           L_S(i2,i1)=0d0
         end do
       end do

       do i1=-lmax,lmax
         do i2=-lmax,lmax
           L_min(i2,i1)=0d0
           L_up(i2,i1)=0d0
         end do
       end do

c  fill the second and the forth quadrant with L_z
c (-L_z,respectively) 


       do i1=1,2*lmax+1
         i1l=i1-lmax-1       ! the value of m (varies from -l to +l)  
         i2=2*lmax+1-(i1-1)  

c         L_S(i2,i1)=icompl*i1l 
         L_S(i2,i1)=-icompl*i1l 

       end do 

       do i1=2*lmax+2,(2*lmax+1)*2
         i1l=i1-lmax-1-(2*lmax+1)       ! the value of m (varies from -l to +l)  
         i2=(2*lmax+1)*2-(i1-(2*lmax+2))  

c         L_S(i2,i1)=-icompl*i1l 
         L_S(i2,i1)=icompl*i1l 

       end do 


c  implement now L_- in the third quadrant

       IF (lmax>0) then

         lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
         l_min(0,-1)=-icompl*lfac
c         l_min(0,-1)=icompl*lfac
         l_min(0,1)=lfac
         l_min(-1,0)=icompl*lfac
         l_min(1,0)=-lfac
       
         IF (lmax > 1) then

            do i1=2,lmax

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-i1*(i1-1d0))
              l_min(-i1,-i1+1)=-lfac
              l_min(-i1,i1-1)=icompl*lfac
              l_min(i1,-i1+1)=-icompl*lfac
              l_min(i1,i1-1)=-lfac

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-(i1-1)*(i1))
              l_min(-i1+1,-i1)=lfac
              l_min(-i1+1,i1)=icompl*lfac
              l_min(i1-1,-i1)=-icompl*lfac
              l_min(i1-1,i1)=lfac

            end do

         END IF
       END IF

 
        do i1=-lmax,lmax
         do i2=-lmax,lmax
           L_S(i2+3*lmax+2,i1+lmax+1)=l_min(i1,i2)
         end do
       end do


c  implement now L_+ in the   quadrant

       IF (lmax>0) then

         lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
         l_up(0,-1)=-icompl*lfac
         l_up(0,1)=-lfac
         l_up(-1,0)=icompl*lfac
         l_up(1,0)=lfac
       
         IF (lmax > 1) then

            do i1=2,lmax

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-i1*(i1-1d0))
              l_up(-i1,-i1+1)=lfac
              l_up(-i1,i1-1)=icompl*lfac
              l_up(i1,-i1+1)=-icompl*lfac
              l_up(i1,i1-1)=lfac

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-(i1-1)*(i1))
              l_up(-i1+1,-i1)=-lfac
              l_up(-i1+1,i1)=icompl*lfac
              l_up(i1-1,-i1)=-icompl*lfac
              l_up(i1-1,i1)=-lfac

            end do

         END IF
       END IF

c
        do i1=-lmax,lmax
         do i2=-lmax,lmax
           L_S(i2+lmax+1,i1+3*lmax+2)=l_up(i1,i2)
         end do
       end do



       deallocate(l_min)
       deallocate(l_up)


       END
