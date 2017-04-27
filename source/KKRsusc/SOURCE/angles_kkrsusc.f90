      subroutine angles_kkrsusc(natom,theta,phi)
!     This routine passes angles from kkrflex to kkrsusc

      use global

      implicit none 
!     --------------------------------------------------------------
      real(kind=r8b), intent(in)  :: theta(natom), phi(natom)
      real(kind=r8b)              :: theta1, phi1
      integer(kind=i4b)           :: natom, ia, ia2
!     --------------------------------------------------------------


      write(*,*) nasusc      
      ! Loop over kkrsusc atoms      
      do ia = 1, nasusc 
        ia2 = iasusc(ia) 
        theta1 = theta(ia2)
        phi1   = phi(ia2)
        magdir(1,ia) = cos(phi1)*sin(theta1)
        magdir(2,ia) = sin(phi1)*sin(theta1)
        magdir(3,ia) = cos(theta1) 
      end do ! iatoms

      end subroutine angles_kkrsusc
