  subroutine charge_density(ia,gfsum)
! charge, spin and orbital densities
! xc kernel
  use global

  implicit none

! Which atom
  integer(kind=i4b), intent(in) :: ia
! Matrix elements of energy integrated GF
  complex(kind=c8b), intent(in) :: gfsum(nlmsb,nlmsb)
! ----------------------------------------------------------------------
  integer(kind=i4b) :: i2(2), i3(3)
  integer(kind=i4b) :: j, jlms, js, jlm, jl, jm, jb
  integer(kind=i4b) :: i, ilms, is, ilm, il, im, ib
  integer(kind=i4b) :: nr, k, klm, ir
  real(kind=r8b),    allocatable :: phi1(:), phi2(:), r(:)
  complex(kind=c8b), allocatable :: rho(:,:,:,:), work(:)
  real(kind=r8b)    :: fac
  complex(kind=c8b) :: znorm
  complex(kind=c8b), external :: radint


  allocate(phi1(nrmax),phi2(nrmax),r(nrmax),work(nrmax),rho(nrmax,0:6,lmmax,lmmax))
! ----------------------------------------------------------------
! density matrix (COMPLEX)
  nr = nrpts(ia)
  r = rmesh(1:nr,ia)
  rho = 0.d0
  do j=1,nlmsba(ia)
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    jlms = lms2i(jlm,js)
    phi2(1:nr) = phiref(1:nr,jb,i2lm(2,jlm),js,ia)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      ilms = lms2i(ilm,is)
      phi1(1:nr) = phiref(1:nr,ib,i2lm(2,ilm),is,ia)
      do k=1,4
!        rho(1:nr,mod(k,4),ilm,jlm) = rho(1:nr,mod(k,4),ilm,jlm) + gfsum(i,j)*pauli(js,is,k)*phi1(1:nr)*phi2(1:nr)
        rho(1:nr,mod(k,4),ilm,jlm) = rho(1:nr,mod(k,4),ilm,jlm) + gfsum(i,j)*ds2c(k,js,is)*phi1(1:nr)*phi2(1:nr)   ! dot product: gfsum and ds2c as vectors
      end do
    end do
  end do
! orbital densities - only charge part (spin trace kills the rest)
! ----------------------------------------------------------------
  do jlm=1,lmmax
    do ilm=1,lmmax
      do k=4,6
        do klm=1,lmmax
          rho(1:nr,k,ilm,jlm) = rho(1:nr,k,ilm,jlm) + lorb(ilm,klm,k-3)*rho(1:nr,0,klm,jlm)
        end do
      end do
    end do
  end do
! ----------------------------------------------------------------
! fill up new rho2ns; project the spin density on the spin axis
!  write(iodb,'("charge_density: ia=",i4," magdir=",3f10.6)') ia, magdir(:,ia)
  new_rho2ns(:,:,:,ia) = 0.d0
  rho_lm(:,:,:,ia) = 0.d0
  do klm=1,lmmax2
    do jlm=1,lmmax
      do ilm=1,lmmax
        fac = rgaunt(ilm,jlm,klm)
        if (abs(fac) > ylmtol) then
!         charge density
          new_rho2ns(1:nr,klm,1,ia) = new_rho2ns(1:nr,klm,1,ia) + fac*rho(1:nr,0,ilm,jlm)
          rho_lm(1:nr,klm,0,ia) = rho_lm(1:nr,klm,0,ia) + fac*rho(1:nr,0,ilm,jlm)/r(1:nr)**2
!         spin density
          do i=1,3
!           project on output spin direction
!            new_rho2ns(1:nr,klm,2,ia) = new_rho2ns(1:nr,klm,2,ia) + fac*newdir(i,ia)*rho(1:nr,i,ilm,jlm)
!           project on input spin direction
            new_rho2ns(1:nr,klm,2,ia) = new_rho2ns(1:nr,klm,2,ia) + fac*magdir(i,ia)*rho(1:nr,i,ilm,jlm)
!           not projected form
            rho_lm(1:nr,klm,i,ia) = rho_lm(1:nr,klm,i,ia) + fac*rho(1:nr,i,ilm,jlm)/r(1:nr)**2
          end do
        end if
      end do
    end do
  end do
! ----------------------------------------------------------------

! kill spin density
  if (lnobxc .and. inobxc(ia) == 1) new_rho2ns(:,:,2,ia) = 0.d0
! ----------------------------------------------------------------
! compute xc potential and kernel
  if (lkxc .and. ikxc(ia) > 0) call get_xc2(ia)
! ----------------------------------------------------------------
! save spherical part of valence densities
!  open(100+ia)
!  open(200+ia)
!    nr = nrpts(ia)
!    r = rmesh(1:nr,ia)
!    nrv(:,:,ia) = 0.d0
!    do ilm=1,lmmax
!      nrv(1:nr,:,ia) = nrv(1:nr,:,ia) + real(rho(1:nr,:,ilm,ilm))
!    end do
!    do ir=1,nr
!      write(100+ia,'(8es16.8)') r(ir), nrv(ir,:,ia)!/r(ir)**2,
!      write(200+ia,'(8es16.8)') r(ir), old_rho2ns(ir,3,:,ia), new_rho2ns(ir,3,:,ia)
!    end do
!  close(100+ia)
!  close(200+ia)
! ----------------------------------------------------------------
  deallocate(phi1,phi2,r,work,rho)
! All done!
  end subroutine charge_density
