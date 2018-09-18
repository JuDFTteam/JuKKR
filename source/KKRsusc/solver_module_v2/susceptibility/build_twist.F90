  subroutine build_twist(twist)
! assembles the xc ALDA kernel in the density basis
! input:  kxc multipoles up to lmmax2 (as GS density)
! output: kxc multipoles up to lmmax0 (as susceptibility)
  use global

  implicit none

  complex(kind=c8b), intent(out) :: twist(4,4,nasusc2,nasusc2)
! -----------------------------------------------------------------
  real(kind=r8b), parameter :: uz(3) = (/0.d0, 0.d0, 1.d0/)
  integer(kind=i4b) :: ia2, ja2, ia, ja, mu, nu, is1, is2, js1, js2
  real(kind=r8b)    :: rotmati(3,3), rotmatj(3,3), twistij(3,3), fac

! ***********************
  atomj: do ja2=1,nasusc2
! ***********************
    ja = iasusc2(ja2)
    if (lrot) then
!     rotation to local frame
      call rotvec(magdir(:,ja),uz,rotmatj)
    else
!     unit matrix
      call rotvec(uz,uz,rotmatj)
    end if
!   ***********************
    atomi: do ia2=1,nasusc2
!   ***********************
      ia = iasusc2(ia2)
      if (lrot) then
!       rotation to local frame
        call rotvec(magdir(:,ia),uz,rotmati)
      else
!       unit matrix
        call rotvec(uz,uz,rotmati)
      end if
!     ------------------------------------------------------------------
!              assemble the twist matrix in the cartesian basis
!     ------------------------------------------------------------------
      twistij(:,:) = 0.d0
      do nu=1,3
        do mu=1,3
          fac = sum(rotmati(3,:)*rotmatj(3,:))*sum(rotmati(mu,:)*rotmatj(nu,:))
          fac = fac - sum(rotmati(3,:)*rotmatj(nu,:))*sum(rotmati(mu,:)*rotmatj(3,:))
          twistij(mu,nu) = fac
        end do
      end do
      if(loutsusc) then 
        write(*,'("ia2,ja2=",2i4)') ia2, ja2
        do mu=1,4
          write(*,'(12f8.4)') twistij(mu,:)
        end do
      end if
!     ------------------------------------------------------------------
!                   convert from cartesian to spin basis
!     ------------------------------------------------------------------
      do js2=1,nsmax
      do js1=1,nsmax
        nu = is2i(js1,js2)
        do is2=1,nsmax
        do is1=1,nsmax
          mu = is2i(is1,is2)
          twist(mu,nu,ia2,ja2) = 0.5d0*sum(pc2s(is1,is2,1:3)*matmul(twistij(:,:),ds2c(1:3,js1,js2)))
        end do
        end do
      end do
      end do
      if(loutsusc) then
        write(*,'("ia2,ja2=",2i4)') ia2, ja2
        do mu=1,4
          write(*,'(12f8.4)') twist(mu,:,ia2,ja2)
        end do
      end if 
!   ************
    end do atomi
!   ************
! ************
  end do atomj
! ************
! All done
  end subroutine build_twist
