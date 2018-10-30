  subroutine get_xc_basis(ia,lmmax,gfsum,kxc,kxcbasis)
! xc potential and energy for atom ia
! The density matrix and the core densities are multiplied by 4pi
! This was factored out in the vxc_vwn subroutine

  implicit none

! which atom, angular momentum cutoff
  integer(kind=i4b), intent(in)  :: ia, lmmax
! density matrix
  complex(kind=c8b), intent(in)  :: gfsum(nlmsb,nlmsb)
! transverse kernel, spherical average
  real(kind=r8b),    intent(out) :: kxc(nrmax)
! full xc kernel, in the density basis
  complex(kind=c8b), intent(out) :: kxcbasis(ndenmax,ndenmax)
! ---------------------------------------------------
  real(kind=r8b), parameter :: fourpi = 16.d0*atan(1.d0)
  real(kind=r8b), parameter :: dx = 1.d-4
  real(kind=r8b), parameter :: magtol = 1.d-8
  integer(kind=i4b) :: ia2, nr, ir, i2(2), i3(3), i4(4), k, n, ill, ngf0, ngf1, nden0, nden1
  integer(kind=i4b) :: i, ilm, il, im, ib, is, iq
  integer(kind=i4b) :: j, jlm, jl, jm, jb, js, jq
  real(kind=r8b)    :: phi1, phi2, den, maglen, magdir(3), vxc, bxc, exc, dr(nrmax)
  real(kind=r8b)    :: rho(4,nll), drho(4,nll)
  real(kind=r8b)    :: vxclm(nrmax,lmmax), bxclm(nrmax,3,lmmax)
  real(kind=r8b)    :: dvxclm(nrmax,lmmax0,ndenmax,4,4)
  real(kind=r8b)    :: exc2(1), vxc2(2), fpirho(2), qcore, mcore
  complex(kind=c8b) :: work(nrmax), norm, rhoden(ndenmax), multipoles(nsmax2,lmmax0)


  if (lmmax < lmmax0) stop 'get_xc_basis: lmmax < lmmax0'

! Check susc
  ia2 = 0
  do i=1,nasusc2
    if (ia == iasusc2(i)) then
      ia2 = i
      exit
    end if
  end do
! Atom not used for susc
  if (ia2 == 0) return

  nr = nrpts(ia)
  dr(1:nr) = drmesh(1:nr,ia)
  vxclm = 0.d0; bxclm = 0.d0; kxc = 0.d0
  dvxclm = 0.d0; kxcbasis = 0.d0
! ----------------------------------------------------------------------
! check the core
  work = nrc(:,ia)
  qcore = real(radint(nr,work,dr),npanat(ia),ircutat(:,ia))
  work = mrc(:,ia)
  mcore = real(radint(nr,work,dr),npanat(ia),ircutat(:,ia))
  write(iodb,'("ia=",i4," qcore, mcore=",2f12.8)') ia, qcore, mcore
! ----------------------------------------------------------------------
! Convert density matrix to the density basis
  rhoden = 0.d0; multipoles = 0.d0
  nden0 = sum(nalmsbden(1:ia-1))
  nden1 = sum(nalmsbden(1:ia))
  do jq=1,nden1-nden0
    ngf0 = sum(nalmsbgf(1:ia-1))
    ngf1 = sum(nalmsbgf(1:ia))
    do iq=1,ngf1-ngf0
      if (abs(dengaunt(iq,jq,ia)) > ylmtol) then
        i3 = i2almsbgf(:,iq+ngf0)
        i = i3(1); j = i3(2)!; ia = i3(3)
        rhoden(jq) = rhoden(jq) + dengaunt(iq,jq,ia)*gfsum(i,j)
      end if
    end do
!   Test: multipoles
    i4 = i2almsbden(:,jq+nden0)
    ib = i4(1); ilm = i4(2); is = i4(3)
    multipoles(is,ilm) = multipoles(is,ilm) + rhoden(jq)*suscnorm(jq+nden0)
  end do
  do ilm=1,lmmax0
    write(iodb,'("get_xc_basis: multipoles=",2i4,8f12.8)') i2lm(:,ilm), multipoles(:,ilm)
  end do
! ----------------------------------------------------------------------
  open(file='xctest.dat',unit=iofile)
! loop over radial points
! **********
  do ir=1,nr
! **********
!   loop over components of the density matrix
    rho = 0.d0
    do iq=1,nden1-nden0
      i4 = i2almsbden(:,iq+nden0)
      ib = i4(1); ilm = i4(2); is = i4(3)
!     from complex to real here
      do k=1,4
        rho(k,:) = rho(k,:) + ds2c(k,i2is(1,is),i2is(2,is))*rhoden(iq)*suscbasis(ir,ib,ilm,is,ia2)*ylm(:,ilm)
      end do
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
      else
        maglen = 0.d0
        magdir = (/0.d0,0.d0,1.d0/)
        rho(1:3,ill) = 0.d0
      end if
      fpirho(1) = den
      fpirho(2) = maglen
      call vxc_vwn(den,maglen,vxc,bxc,exc)
!      call vosko(exc2,fpirho,vxc2,1,1)
!      write(iodb,'("get_xc: diff vosko=",6es16.8)') vxc-bxc, vxc+bxc, exc, vxc2, exc2
      vxclm(ir,:) = vxclm(ir,:) + vxc*wll(ill)*ylm(ill,1:lmmax)
      do k=1,3
        bxclm(ir,k,:) = bxclm(ir,k,:) + magdir(k)*bxc*wll(ill)*ylm(ill,1:lmmax)
      end do
      kxc(ir) = kxc(ir) + wll(ill)*bxc/(maglen + magtol)
    end do
    kxc(ir) = kxc(ir)/rmesh(ir,ia)**2
    write(iofile,'(100es16.8)') rmesh(ir,ia), vr(ir,ia), br(ir,ia), vxclm(ir,1), bxclm(ir,:,1), kxc(ir)
!   Now for the full kernel, in the density basis
!   Do one cartesian component at a time
    do n=1,2
    do jq=1,nden1-nden0
      i4 = i2almsbden(:,jq+nden0)
      jb = i4(1); jlm = i4(2); js = i4(3)
!     perturbation proportional to a basis function
      if (abs(ds2c(n,i2is(1,js),i2is(2,js))) > ylmtol) then
        drho = rho
        drho(n,:) = rho(n,:) + dx*suscbasis(ir,jb,jlm,js,ia)*ylm(:,jlm)/rmesh(ir,ia)**2
!       Expand the xc potentials in spherical harmonics
        do ill=1,nll
          den = drho(4,ill)
          maglen = sqrt(dot_product(drho(1:3,ill),drho(1:3,ill)))
          if (maglen > magtol) then
            magdir = drho(1:3,ill)/maglen
          else
            maglen = 0.d0
            magdir = (/0.d0,0.d0,1.d0/)
            drho(1:3,ill) = 0.d0
          end if
          fpirho(1) = den
          fpirho(2) = maglen
          call vxc_vwn(den,maglen,vxc,bxc,exc)
!          call vosko(exc2,fpirho,vxc2,1,1)
!          write(iodb,'("get_xc: diff vosko=",6es16.8)') vxc-bxc, vxc+bxc, exc, vxc2, exc2
!          dvxclm(ir,:,jq,4,n) = dvxclm(ir,:,jq,4,n) + vxc*wll(ill)*ylm(ill,1:lmmax0)
          do k=1,2
            dvxclm(ir,:,jq,k,n) = dvxclm(ir,:,jq,k,n) + magdir(k)*bxc*wll(ill)*ylm(ill,1:lmmax0)
          end do
        end do
!        dvxclm(ir,:,jq,4,n) = (dvxclm(ir,:,jq,4,n) - vxclm(ir,1:lmmax0))/dx
        do k=1,2
          dvxclm(ir,:,jq,k,n) = (dvxclm(ir,:,jq,k,n) - bxclm(ir,k,1:lmmax0))/dx
        end do
      end if
    end do
    end do
! ******
  end do
! ******
! ----------------------------------------------------------------------
! Back to the density basis
!  dr(1:nr) = dr(1:nr)/rmesh(1:nr,ia)**2
  do n=1,2
  do jq=1,nden1-nden0
    i4 = i2almsbden(:,jq+nden0)
    jb = i4(1); jlm = i4(2); js = i4(3)
    do iq=1,nden1-nden0
      i4 = i2almsbden(:,iq+nden0)
      ib = i4(1); ilm = i4(2); is = i4(3)
      work = 0.d0
      do k=1,2
        work(1:nr) = work(1:nr) + pc2s(i2is(1,is),i2is(2,is),k)*dvxclm(1:nr,ilm,jq,k,n)*ds2c(n,i2is(1,js),i2is(2,js))
      end do
      work(1:nr) = work(1:nr)*suscbasis(1:nr,ib,ilm,is,ia)!*suscbasis(1:nr,jb,jlm,js,ia)
      norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
      if (abs(norm) < 1.d-3) then
        norm = 0.d0
      else
        write(iodb,'("get_xc_basis: kxcbasis=",5i4,2es16.8)') ia, is, js, iq, jq, norm
      end if
      kxcbasis(iq,jq) = kxcbasis(iq,jq) + norm
    end do
  end do
  end do
  close(iofile)
! All done!
  end subroutine get_xc_basis

