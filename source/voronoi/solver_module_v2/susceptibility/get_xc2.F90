  subroutine get_xc2(ia,lmmax,kxclm)
! xc potential and energy for atom ia
! Uses the density from the impurity program
! The density matrix and the core densities are multiplied by 4pi
! This was factored out in the vxc_vwn subroutine

  use global

  implicit none

! which atom, angular momentum cutoff
  integer(kind=i4b), intent(in)    :: ia, lmmax
! full xc kernel
  real(kind=r8b),    intent(inout) :: kxclm(nrmax,lmmax,4,4,nasusc)
! ---------------------------------------------------
  real(kind=r8b), parameter :: fourpi = 16.d0*atan(1.d0)
  real(kind=r8b), parameter :: dx = 1.d-6  ! optimized by hand - NOT COOL
  real(kind=r8b), parameter :: magtol = 1.d-8
  real(kind=r8b), parameter :: tol = 1.d-8
  integer(kind=i4b) :: nr, ir, i2(2), i3(3), k, n, ill
  integer(kind=i4b) :: i, ilm, il, im, ib, is
  integer(kind=i4b) :: j, jlm, jl, jm, jb, js
  real(kind=r8b)    :: phi1, phi2, den, maglen, magdir(3), vxc, bxc, exc
  real(kind=r8b)    :: rho(4,nll), vxclm(nrmax,lmmax), bxclm(nrmax,3,lmmax)
  real(kind=r8b)    :: exc2(1), vxc2(2), fpirho(2), qcore, mcore
  real(kind=r8b)    :: dvdn, dvdm, dbdn, dbdm, trans(3,3), longt(3,3)
  real(kind=r8b)    :: kxc(nrmax)
  complex(kind=c8b) :: work(nrmax)



! loop over radial points
  nr = nrpts(ia); vxclm = 0.d0; bxclm = 0.d0
  kxc = 0.d0; kxclm(:,:,:,:,ia) = 0.d0
! Test
  work = nrc(:,ia)
  qcore = real(radint(nr,work(1:nr),drmesh(1:nr,ia)),npanat(ia),ircutat(:,ia))
  work = mrc(:,ia)
  mcore = real(radint(nr,work(1:nr),drmesh(1:nr,ia)),npanat(ia),ircutat(:,ia))
  write(iodb,'("ia=",i4," qcore, mcore=",2f12.8)') ia, qcore, mcore
! Test
  open(file='xctest.dat',unit=iofile)
  do ir=1,nr
!   ++++++++++++++++++++++++++++++++++++++++++
!   loop over components of the density matrix
!   ++++++++++++++++++++++++++++++++++++++++++
    rho = 0.d0
    do ilm=1,lmmax2
      if (abs(gs_qlm(ilm,ia)) > tol .or. abs(gs_mlm(ilm,ia)) > tol) then
        rho(3,:) = rho(3,:) + old_rho2ns(ir,ilm,2,ia)*ylm(:,ilm)
        rho(4,:) = rho(4,:) + old_rho2ns(ir,ilm,1,ia)*ylm(:,ilm)
      end if
    end do
!    write(iodb,'("get_xc: rho=",4es16.8)') rho(:,1)*rmesh(ir,ia)**2
!   ++++++++++++++++++++++++++++++++++++++++++
!   Add the core charge and magnetization densities
    rho(3,:) = rho(3,:) + mrc(ir,ia)
    rho(4,:) = rho(4,:) + nrc(ir,ia)
!   Divide by r^2; the 4pi was taken care of in vxc_vwn
    rho = rho/rmesh(ir,ia)**2
!   Expand the xc potentials in spherical harmonics
    do ill=1,nll
      den = rho(4,ill)
      maglen = sqrt(dot_product(rho(1:3,ill),rho(1:3,ill)))
      if (maglen > magtol) then
        magdir = rho(1:3,ill)/maglen
!        magdir = (/0.d0,0.d0,1.d0/)
      else
        magdir = (/0.d0,0.d0,1.d0/)
        maglen = 0.d0
      end if
      fpirho(1) = den
      fpirho(2) = maglen
      call vxc_vwn(den,maglen,vxc,bxc,exc)
!      call vosko(exc2,fpirho,vxc2,1,1)
!      write(iodb,'("get_xc: diff vosko=",6es16.8)') vxc-bxc, vxc+bxc, exc, vxc2, exc2
      vxclm(ir,1:lmmax) = vxclm(ir,1:lmmax) + vxc*wll(ill)*ylm(ill,1:lmmax)
      do k=1,3
        bxclm(ir,k,1:lmmax) = bxclm(ir,k,1:lmmax) + magdir(k)*bxc*wll(ill)*ylm(ill,1:lmmax)
      end do
!     *************************************************
!     ****************    xc kernel    ****************
!     *************************************************
!     transverse and longitudinal parts
      trans = 0.d0; longt = 0.d0
      do n=1,3
        do k=1,3
          longt(k,n) = magdir(k)*magdir(n)
          trans(k,n) = -longt(k,n)
        end do
        trans(n,n) = trans(n,n) + 1.d0
      end do
!     transverse part of the xc kernel
      if (maglen > magtol .and. ikxc /= 2) then
        kxc(ir) = kxc(ir) + wll(ill)*bxc/maglen
        do n=1,3
          do k=1,3
            kxclm(ir,1:lmmax,k,n,ia) = kxclm(ir,1:lmmax,k,n,ia) + trans(k,n)*(bxc/maglen)*wll(ill)*ylm(ill,1:lmmax)
          end do
        end do
      end if
!     derivatives
      dvdn = vxc; dvdm = vxc; dbdn = bxc; dbdm = bxc
      call vxc_vwn(den+dx,maglen,vxc,bxc,exc)
      dvdn = (vxc - dvdn)/dx; dbdn = (bxc - dbdn)/dx
      call vxc_vwn(den,maglen+dx,vxc,bxc,exc)
      dvdm = (vxc - dvdm)/dx; dbdm = (bxc - dbdm)/dx
!     symmetry: comment this out for debugging
      dvdm = 0.5d0*(dvdm + dbdn)
      dbdn = dvdm
!     longitudinal part of the xc kernel
      if (ikxc /= 1) then
        do n=1,3
          do k=1,3
            kxclm(ir,1:lmmax,k,n,ia) = kxclm(ir,1:lmmax,k,n,ia) + longt(k,n)*dbdm*wll(ill)*ylm(ill,1:lmmax)
          end do
          kxclm(ir,1:lmmax,n,4,ia) = kxclm(ir,1:lmmax,n,4,ia) + magdir(n)*dbdn*wll(ill)*ylm(ill,1:lmmax)
          kxclm(ir,1:lmmax,4,n,ia) = kxclm(ir,1:lmmax,4,n,ia) + magdir(n)*dvdm*wll(ill)*ylm(ill,1:lmmax)
        end do
        kxclm(ir,1:lmmax,4,4,ia) = kxclm(ir,1:lmmax,4,4,ia) + dvdn*wll(ill)*ylm(ill,1:lmmax)
      end if
!     *************************************************
    end do
    kxc(ir) = kxc(ir)/rmesh(ir,ia)**2  ! why?
    kxclm(ir,:,:,:,ia) = kxclm(ir,:,:,:,ia)/rmesh(ir,ia)**2  ! why?
    write(iofile,'(100es16.8)') rmesh(ir,ia), vr(ir,ia), br(ir,ia), vxclm(ir,1), bxclm(ir,:,1), kxc(ir)
  end do
  close(iofile)
! Test: multipoles
  do ilm=1,lmmax
    do n=1,4
      do k=1,4
        work(1:nr) = kxclm(1:nr,ilm,k,n,ia)
        rho(k,1) = radint(nr,work(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        if (abs(rho(k,1)) < 1.d-6) then
          rho(k,1) = 0.d0
          kxclm(:,ilm,k,n,ia) = 0.d0
        end if
      end do
      write(iodb,'("get_xc: multipoles of kxc for ia,ilm=",3i4,4es16.8)') ia, i2lm(:,ilm), rho(:,1)
    end do
  end do
  if (ikxc == 0) kxclm(:,:,:,:,ia) = 0.d0
! All done!
  end subroutine get_xc2

