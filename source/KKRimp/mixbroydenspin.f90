!------------------------------------------------------------------------------------
!> Summary:
!> Author:
!>
!------------------------------------------------------------------------------------
MODULE mod_mixbroydenspin

contains

   !====================================================================================================================
   !-------------------------------------------------------------------------------
   !> Summary: Broyden mixing: spin
   !> Author:
   !> Category: KKRimp, potential
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !>
   !-------------------------------------------------------------------------------
   subroutine mixbroydenspin (natom,density,max_iter,iter) 
      use type_density
      use mod_types, only: t_inc
      use mod_broyden, only: broyden
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

      if (t_inc%i_write>0) write(1337,*) 'Mixing of angles with Broyden mixing'
      write(*,*)    'Mixing of angles with Broyden mixing'
      write(*,*)   'new theta1 [deg]',density(iatom)%theta*180/pi
      if (t_inc%i_write>0) write(1337,*)'new theta1 [deg]',density(iatom)%theta*180/pi
      write(*,*)   'new phi1 [deg]',density(iatom)%phi*180/pi
      if (t_inc%i_write>0) write(1337,*)'new phi1 [deg]',density(iatom)%phi*180/pi
      !   write(23452324,'(5000F)') density(iatom)%theta*180/pi,density(iatom)%phi*180/pi

      end do

   end subroutine mixbroydenspin


   !-------------------------------------------------------------------------------
   !> Summary: Broyden mixing: spinangle
   !> Author:
   !> Category: KKRimp, potential
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !>
   !-------------------------------------------------------------------------------
   subroutine mixbroydenspinangle (natom,density,max_iter,iter) 
      use type_density
      use mod_types, only: t_inc
      use mod_broyden, only: broyden
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

      if (t_inc%i_write>0) write(1337,*) 'Mixing of angles with Broyden mixing'
      write(*,*)    'Mixing of angles with Broyden mixing'
      write(*,*)   'new theta1 [deg]',density(iatom)%theta*180/pi
      if (t_inc%i_write>0) write(1337,*)'new theta1 [deg]',density(iatom)%theta*180/pi
      write(*,*)   'new phi1 [deg]',density(iatom)%phi*180/pi
      if (t_inc%i_write>0) write(1337,*)'new phi1 [deg]',density(iatom)%phi*180/pi
      !   write(23452324,'(5000F)') density(iatom)%theta*180/pi,density(iatom)%phi*180/pi

      end do

   end subroutine mixbroydenspinangle

END  MODULE mod_mixbroydenspin
