MODULE mod_mixbroydenspin

contains

!====================================================================================================================
       subroutine mixbroydenspin (natom,density,max_iter,iter) 
use type_density
implicit none
integer            ::  natom
type(density_type) ::  density(natom)
!local

double precision,allocatable :: vector(:,:) !(mvlen,2)
double precision             :: alpha
integer                      :: mvlen,vlen
double precision             :: rms
integer                      :: iter,n_init
integer                      :: mbroylen
integer                      :: max_iter
double precision             :: totmagmoment,totxymagmoment
double precision             :: magmoment(3)
integer                      :: iatom,ipos
double precision             :: totmagmoment_temp,totmagmoment_nomix(natom)

allocate( vector(3*natom,2) )
mvlen=natom*3
vlen=mvlen
alpha=0.01
mbroylen=max_iter
n_init=1 !number of linear mixing steps


if (iter==1) then
  do iatom=1,natom
    density(iatom)%magmomentold(1)=density(iatom)%magmoment(1)
    density(iatom)%magmomentold(2)=density(iatom)%magmoment(2)
    density(iatom)%magmomentold(3)=density(iatom)%magmoment(3)
  end do
end if


  do iatom=1,natom
    totmagmoment_nomix(iatom) = sqrt ( density(iatom)%magmoment(1)**2+ &
                                       density(iatom)%magmoment(2)**2+ &
                                       density(iatom)%magmoment(3)**2 )
  end do


! print *, vector
do iatom=1,natom
  ipos=(iatom-1)*3+1
  vector(ipos,  2)=density(iatom)%magmoment(1)
  vector(ipos+1,2)=density(iatom)%magmoment(2)
  vector(ipos+2,2)=density(iatom)%magmoment(3)
  vector(ipos,  1)=density(iatom)%magmomentold(1)
  vector(ipos+1,1)=density(iatom)%magmomentold(2)
  vector(ipos+2,1)=density(iatom)%magmomentold(3)
end do !iatom
print *,vector

rms=0.0d0
do ipos=1,3*natom
    rms=rms+(vector(ipos,2)-vector(ipos,1) )**2
end do
rms=sqrt(rms)


call broyden (vector, vlen, alpha, rms, iter,  &
                         n_init,mbroylen,mvlen) 
print *,vector

write(*,*) 'Broyden spin mixing'
write(*,*) 'rms for spin is ',rms
do iatom=1,natom
  ipos=(iatom-1)*3+1

!     density(iatom)%magmomentold(1)    = density(iatom)%magmoment(1)
!     density(iatom)%magmomentold(2)    = density(iatom)%magmoment(2)
!     density(iatom)%magmomentold(3)    = density(iatom)%magmoment(3)



  if (density(iatom)%magmomentfixed/=1) then
    density(iatom)%magmoment(1)    = vector(ipos,  2)
    density(iatom)%magmoment(2)    = vector(ipos+1,2)
    density(iatom)%magmoment(3)    = vector(ipos+2,2)

    
    totmagmoment_temp = sqrt ( density(iatom)%magmoment(1)**2+ &
                               density(iatom)%magmoment(2)**2+ &
                               density(iatom)%magmoment(3)**2 )
! totmagmoment_nomix
    density(iatom)%magmoment(1)    = vector(ipos,  2)/totmagmoment_temp*totmagmoment_nomix(iatom)
    density(iatom)%magmoment(2)    = vector(ipos+1,2)/totmagmoment_temp*totmagmoment_nomix(iatom)
    density(iatom)%magmoment(3)    = vector(ipos+2,2)/totmagmoment_temp*totmagmoment_nomix(iatom)
    print *,totmagmoment_temp,totmagmoment_nomix(iatom)
!     density(iatom)%magmomentold(1) = vector(ipos,  1)
!     density(iatom)%magmomentold(2) = vector(ipos+1,1)
!     density(iatom)%magmomentold(3) = vector(ipos+2,1)
  end if
end do !iatom

do iatom=1,natom
  magmoment=density(iatom)%magmoment
  totmagmoment=SQRT(magmoment(1)**2+magmoment(2)**2+magmoment(3)**2)
  totxymagmoment=SQRT(magmoment(1)**2+magmoment(2)**2)

  density(iatom)%theta= acos(magmoment(3)/totmagmoment)
  density(iatom)%phi  = datan2(magmoment(2),magmoment(1))

  write(1337,*) 'Mixing of angles with Broyden mixing'
  write(*,*)    'Mixing of angles with Broyden mixing'
  write(*,*)   'new theta1 [deg]',density(iatom)%theta*180/pi
  write(1337,*)'new theta1 [deg]',density(iatom)%theta*180/pi
  write(*,*)   'new phi1 [deg]',density(iatom)%phi*180/pi
  write(1337,*)'new phi1 [deg]',density(iatom)%phi*180/pi
!   write(23452324,'(5000F)') density(iatom)%theta*180/pi,density(iatom)%phi*180/pi

end do

end subroutine mixbroydenspin


      subroutine mixbroydenspinangle (natom,density,max_iter,iter) 
use type_density
implicit none
integer            ::  natom
type(density_type) ::  density(natom)
!local

double precision,allocatable :: vector(:,:) !(mvlen,2)
double precision             :: alpha
integer                      :: mvlen,vlen
double precision             :: rms
integer                      :: iter,n_init
integer                      :: mbroylen
integer                      :: max_iter
integer                      :: iatom,ipos

allocate( vector(2*natom,2) )
mvlen=natom*2
vlen=mvlen
alpha=0.01
mbroylen=max_iter
n_init=1 !number of linear mixing steps

! print *, vector
do iatom=1,natom
  ipos=(iatom-1)*2+1
  vector(ipos,  2)=density(iatom)%theta
  vector(ipos+1,2)=density(iatom)%phi
  vector(ipos,  1)=density(iatom)%thetaold
  vector(ipos+1,1)=density(iatom)%phiold
end do !iatom
print *,vector

rms=0.0d0
do ipos=1,2*natom
    rms=rms+(vector(ipos,2)-vector(ipos,1) )**2
end do
rms=sqrt(rms)


call broyden (vector, vlen, alpha, rms, iter,  &
                         n_init,mbroylen,mvlen) 
print *,vector

write(*,*) 'Broyden spin mixing'
write(*,*) 'rms for spin is ',rms
do iatom=1,natom
  ipos=(iatom-1)*2+1

  if (density(iatom)%magmomentfixed/=1) then
    density(iatom)%theta    = vector(ipos,  2)
    density(iatom)%phi    = vector(ipos+1,2)
  end if
end do !iatom

do iatom=1,natom
!   magmoment=density(iatom)%magmoment
!   totmagmoment=SQRT(magmoment(1)**2+magmoment(2)**2+magmoment(3)**2)
!   totxymagmoment=SQRT(magmoment(1)**2+magmoment(2)**2)

!   density(iatom)%theta= acos(magmoment(3)/totmagmoment)
!   density(iatom)%phi  = datan2(magmoment(2),magmoment(1))

  write(1337,*) 'Mixing of angles with Broyden mixing'
  write(*,*)    'Mixing of angles with Broyden mixing'
  write(*,*)   'new theta1 [deg]',density(iatom)%theta*180/pi
  write(1337,*)'new theta1 [deg]',density(iatom)%theta*180/pi
  write(*,*)   'new phi1 [deg]',density(iatom)%phi*180/pi
  write(1337,*)'new phi1 [deg]',density(iatom)%phi*180/pi
!   write(23452324,'(5000F)') density(iatom)%theta*180/pi,density(iatom)%phi*180/pi

end do

end subroutine mixbroydenspinangle



       subroutine broyden (vector, vlen, alpha, rms, iter,  &
                         n_init,mbroylen,mvlen) 

! ,ntasks)

!

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     History: Original code written by D.D. Johnson (see PRB 38, 12807)

!                  note:  there are a few typos in that paper but 

!                  the code is working!

!              Rewritten by W. A. Shelton for LSMS code 6/21/94

!                  this version is easy to read (no goto!!!! more comments ...)

!                  and is setup for MPP machines (not tested)

!              Rewritten by T. C. Schulthess, ORNL, March 97

!                  this version should work for any code (see comments below)

!

!     Bug fixes:   TCS, 8/5/97 see comments below 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!

!     further comments on how to use this subroutine:

!     (if nothing useful stands here I had no time yet to update these

!     comments, please consult usage in lkkr code version 0.6.3 or later,

!     or call Thomas Schulthess (423) 5768067)

!

!     vector(r,i) -> i=1: old vector (input), scratch (ouput)

!                 -> i=2: new vector (input), mixed vector (output)

!     vlen        -> length of vector

!     alpha       -> linear mixing factor

!     rms         -> RMS difference between old and new vector

!     iter        -> iteration number (if 1, linear mixing, broyden reset)

!     broylen     -> number of iterations that are used for mixing (<=mbroylen)

!     u, vt, f, df, vold, and w are working arrays that need to be saved

!                   between call to this subroutine

!     a, b, d, cm, and ipiv are working arrays that need not be saved

!     mbroylen    -> maximum number of iterations that can be saved

!     mvlen       -> maximum length of vectors

!

!     See declaration for exact dimentsions and types

!

!     There are two options for matrix inversions, a Gaussian

!     elimination routine called invert1 and calls to lapack routines

!     with pivoting (see comments "using invert1" and "using lapack").

!     Obviously only one can be used, comment out the other one.

!

!     When using this subroutine in a parallel code in which only parts

!     of the vectors are known on every node, make sure that the calls

!     to gldsum (global sum) are correct (LKKR and LSMS codes have

!     different calls).

!

!     In a serial code, either comment out the calls to glbsum or

!     provide a dummy subroutine

!

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       implicit none 

!


!       

       integer mbroylen,mvlen ! used for dimensioning 

!     

       double precision vector(mvlen,2)

       integer vlen

       double precision alpha

       double precision rms

       integer iter

       integer n_init,broylen

       double precision, allocatable, save:: u(:,:)

       double precision, allocatable, save:: vt(:,:)

       double precision, allocatable, save:: f(:)

       double precision, allocatable, save:: df(:)

       double precision, allocatable, save:: vold(:)

       double precision a(mbroylen,mbroylen)

       double precision b(mbroylen,mbroylen)

       double precision d(mbroylen,mbroylen)

       double precision cm(mbroylen)

       double precision, allocatable, save:: w(:)

       integer ipiv(mbroylen) 

!

       integer i,j,k,info

       double precision fac1,fac2,fnorm,dfnorm,w0,aij,gmi,cmj,wtmp 

!

       integer lastit,lastm1,nn

       double precision zero,one

       parameter (zero=0.d0,one=1.d0)

       double precision amix

       save lastit,amix 

      

       broylen=mbroylen

!

!      if (broylen > mbroylen) then

!         write(6&

!      &        ,'('' broyden: broylen='',i5,'' exeeds mbroylen='',i5)'&

!      &        ) broylen,mbroylen

!         stop

!      endif 

!



       if (.not. allocated(u)) then

         allocate(u(mvlen,mbroylen))

         allocate(vt(mvlen,mbroylen))

         allocate(f(mvlen))

         allocate(df(mvlen))

         allocate(vold(mvlen))

         allocate(w(mbroylen))

       end if



!

       if( iter <=  n_init)then 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     first n_init iteration: preform linear mixing, load f and vold, set

!                      different pointers and variables

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         lastit = iter-1          ! initialize pointers

         lastm1 = lastit-1 



         amix = alpha           ! for safety reasons

         do k = 1,vlen

            f(k) = vector(k,2) - vector(k,1)

            vold(k) = vector(k,1)

         enddo 



         do k = 1,vlen          ! this is the linear mixing

            vector(k,2) = vector(k,1) + amix * f(k)

         enddo 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       else 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     iter > n_init: this is where the non-linear mixing is done

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         lastit = lastit+1      ! update pointers

         lastm1 = lastit-1 



         if( iter > broylen) then ! set current lenght of broyden cycle

            nn = broylen

         else

            nn = lastit         !lastm1

         endif 



         w0=.01d0               ! set weighting factor for the zeroth iteration 



!---- find: f[i] := vector(2)[i] - vector(1)[i]





         do k = 1,vlen

            df(k) = vector(k,2) - vector(k,1) - f(k)

         enddo

         do k = 1,vlen

            f(k) = vector(k,2) - vector(k,1)

         enddo 



!---- find: fnorm  := |f|



         dfnorm = zero

         fnorm = zero

         do k = 1,vlen

            dfnorm = dfnorm + df(k)*df(k)

            fnorm  = fnorm  + f(k)*f(k)

         enddo 



         dfnorm = sqrt( dfnorm )

         fnorm  = sqrt( fnorm ) 



!---- set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|



         fac2 = one/dfnorm

         fac1 = amix*fac2

         do k = 1,vlen

            vector(k,2) = fac1*df(k) + fac2*(vector(k,1) - vold(k))

            vold(k) = vector(k,1)

            vector(k,1) = fac2*df(k)

         enddo 



!---- store vector(1) and vector(2) in the stacks u and vt restpectively



         

                 call broy_sav(u,vt,vector,iter-1,broylen,vlen,mvlen) 

         



!---- calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:



         do j=1,nn - 1          ! off diagonal elements of a(i,j)

           do i = j+1,nn

              aij = zero

              do k = 1,vlen

                 aij = aij + vt(k,j)*vt(k,i)

              enddo 

              a(i,j) = aij

              a(j,i) = aij

           enddo

        enddo 



        do i = 1,nn             ! diagonal elements a(i,i) and cm(i)

           aij = zero

           cmj = zero

           do k=1,vlen

              cmj = cmj + vt(k,i)*f(k)

              aij = aij + vt(k,i)*vt(k,i)

           enddo 

           a(i,i) = aij

           cm(i) = cmj

        enddo 



!---- shift down weights in stack



        if(iter-1 > broylen)then

           do i=1,broylen-1

              w(i)=w(i+1)

           enddo

        endif

        wtmp = zero

        if( rms > 1.0d-09 ) wtmp=2.0*sqrt(0.010d+00/rms)

        if( wtmp < one )    wtmp=1.00d+00

        if(iter > broylen)then

           w(broylen)=wtmp

        else

           w(lastit)=wtmp       !w(lastm1)=wtmp

        endif 



!---- now calculate the b-matrix:

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> uses lapack

        do i=1,nn

           do j=1,nn

              b(j,i)= a(j,i)*w(j)*w(i)

           enddo

          b(i,i)= w0**2 + a(i,i)*w(i)*w(i)

        enddo

        call dgetrf(nn,nn,b,mbroylen,ipiv,info)

        call dgetri(nn,b,mbroylen, ipiv, d, nn, info )

      !  write(6,*) ' optimum lwork', d(1,1)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< uses lapack

!---- mix vectors: 



        do k=1,vlen

           vector(k,2)= vold(k) + amix * f(k)

        enddo

        do i=1,nn

           gmi = zero

           do j=1,nn

              gmi = gmi + cm(j)*b(j,i)*w(j)

           enddo

           do k=1,vlen

              vector(k,2) = vector(k,2) - gmi*u(k,i)*w(i)

           enddo

        enddo 



       endif 

      

      




      

      end subroutine



!******************************************************************************


!     ==================================================================

      subroutine broy_sav(fins,fots,vector,itscf,istore,ivsiz,mivsiz) 

!     ==================================================================

      implicit none

      integer itscf

      integer ivsiz

      integer mivsiz

      integer istore 

!     ==================================================================

      double precision  vector(mivsiz,2)

      double precision  fins(mivsiz,istore)

      double precision  fots(mivsiz,istore) 
      
!     ==================================================================

      integer i, j

!      write(6,'('' IN BROY_SAV: istore,itscf '',2i5)') istore,itscf

!     ==================================================================

      if( itscf <= istore ) then 

!     ==================================================================

!     Load the first istore iterations in increasing iteration count

!     ==================================================================

        do i = 1,ivsiz

          fins(i,itscf) = vector(i,2)

        enddo 

!     ==================================================================

        do i = 1,ivsiz

          fots(i,itscf) = vector(i,1)

        enddo 

!     ==================================================================

      else 

!     ==================================================================

!     Re-load so that the ordering is in increasing iteration count

!     ==================================================================

        do j = 1,istore - 1 

!          write(6,'('' IN BROY_SAV: j,j+1 '',2i5)') j,j+1

          do i = 1,ivsiz

            fins(i,j) = fins(i,j+1)

          enddo 

!     ==================================================================

          do i = 1,ivsiz

            fots(i,j) = fots(i,j+1)

          enddo 

!     ==================================================================

        enddo 

!     ==================================================================

!     Load current charge densities in the last storage location

!     ==================================================================

        do i = 1,ivsiz

          fins(i,istore) = vector(i,2)

        enddo 

!     ==================================================================

        do i = 1,ivsiz

          fots(i,istore) = vector(i,1)

        enddo 

!     ==================================================================

      endif 

!     ==================================================================

      

      end subroutine






END  MODULE mod_mixbroydenspin