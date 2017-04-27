  subroutine build_vscfb2(bxctwist)
! assembles the SCF potential in the projection basis
! converts from cartesian to spin representation
! 
  use global

  implicit none

! twist Bxc
  complex(kind=c8b), intent(out) :: bxctwist(nlmsb,nlmsb,3,nasusc)
! -----------------------------------------------------------------
  real(kind=r8b), parameter :: uz(3) = (/0.d0,0.d0,1.d0/)
  integer(kind=i4b) :: ia, is, js, ib, jb, nbi, nbj, i3(3), i2(2), i4(4), ia2
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j, nr, k, ir, s1, s2, iq
  real(kind=r8b)    :: dr(nrmax), maxnorm, uvec(3)
  complex(kind=c8b) :: work(nrmax), norm
  real(kind=r8b)    :: potcart(nrmax,4)
  complex(kind=c8b) :: potspin(nrmax,nsmax,nsmax), vdiff(nrmax,4)
  integer(kind=i4b) :: imax, jmax
  complex(kind=c8b), external :: radint
  real(kind=r8b)    :: rotmat(3,3,nasusc), magdirnew(3)

  vlmsbgf = 0.d0; bxctwist = 0.d0
! use the ASA
! use the fact that the basis was constructed per l channel
! ******************
  do ia=1,nasusc
! ******************
!   rotation matrices to local frame
!    call rotvec(uz,magdir(:,ia),rotmat(:,:,ia))
    call rotvec(uz,uz,rotmat(:,:,ia))
    magdirnew(:) = matmul(transpose(rotmat(:,:,ia)),magdir(:,ia))
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)
    maxnorm = 0.d0
!   ----------------------------------------------------------------------
!   potential in cartesian labels
!   this is the real space orientation of the potential
!    uvec = (/ 0.d0, 1.d0, 0.d0/)
    potcart = 0.d0
!   magnetic
    do k=1,3
      potcart(1:nr,k) = br(1:nr,ia)*magdir(k,ia)
    end do
!   charge
!    potcart(1:nr,4) = vr(1:nr,ia) - 2.d0*zat(ia)*rmesh(1:nr,ia)
!   potential from cartesian to spin labels
    potspin = 0.d0
    do js=1,2
      do is=1,2
        do k=1,4
          potspin(1:nr,is,js) = potspin(1:nr,is,js) + pc2s(is,js,k)*potcart(1:nr,k)
        end do
      end do
    end do
!   ----------------------------------------------------------------------
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
!     ----------------------------------------------------------------
!     selection rules
        if (il == jl .and. im == jm) then
!       revise the basis
          work(1:nr) = phiref(1:nr,ib,il,is,ia)*potspin(1:nr,is,js)*phiref(1:nr,jb,jl,js,ia)
          norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
          vlmsbgf(i,j,ia) = norm
!         components of \vec{B} \times \vec{\sigma}
          work(1:nr) = phiref(1:nr,ib,il,1,ia)*br(1:nr,ia)*phiref(1:nr,jb,jl,1,ia)
          norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
          bxctwist(i,j,1,ia) = bxctwist(i,j,1,ia) + norm*(magdirnew(2)*pauli(is,js,3) - magdirnew(3)*pauli(is,js,2))
          bxctwist(i,j,2,ia) = bxctwist(i,j,2,ia) + norm*(magdirnew(3)*pauli(is,js,1) - magdirnew(1)*pauli(is,js,3))
          bxctwist(i,j,3,ia) = bxctwist(i,j,3,ia) + norm*(magdirnew(1)*pauli(is,js,2) - magdirnew(2)*pauli(is,js,1))
!          write(iodb,'(4i4,2es16.8)') il, im, is, js, norm
          if (abs(norm) > maxnorm) then
            maxnorm = abs(norm)
            imax = i; jmax = j
          end if
        end if
!     ----------------------------------------------------------------
      end do
    end do
    if (lhdio) write(iodb,'(" SCF potential, ia=",i4)') ia
    if (lhdio) write(iodb,'("vscfb norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!    write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
!    do k=1,4
!    do jlm=1,lmmax
!    do ilm=1,lmmax
!      do js=1,nsmax
!      do is=1,nsmax
!        do jb=1,iwsusc(i2lm(2,jlm),js,ia)
!        do ib=1,iwsusc(i2lm(2,ilm),is,ia)
!          j = lmsb2i(jb,jlm,js,ia)
!          i = lmsb2i(ib,ilm,is,ia)
!          if (abs(vscfb(i,j,ia)) > 1.d-8) then
!            write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vscfb(i,j,ia)
!          end if
!        end do
!        end do
!      end do
!      end do
!    end do
!    end do
!    end do
! ******************
  end do
! ******************
! ----------------------------------------------------------------------
  if (.not.lsusc) return
! ----------------------------------------------------------------------
! Now construct the potential in the susceptibility basis
  vlmsbden = 0.d0
! ******************
  do ia2=1,nasusc2
! ******************
    ia = iasusc2(ia2)
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)!/rmesh(1:nr,ia)**2
!   --------------------------------------------------------------------
!   potential in cartesian labels
!   this is the real space orientation of the potential
!    uvec = (/ 0.d0, 1.d0, 0.d0/)
    potcart = 0.d0
!   magnetic
    do k=1,3
      potcart(1:nr,k) = br(1:nr,ia)*magdir(k,ia)
    end do
!   charge
!    potcart(1:nr,4) = vr(1:nr,ia) - 2.d0*zat(ia)*rmesh(1:nr,ia)
!   potential from cartesian to spin labels
    potspin = 0.d0
    do js=1,2
      do is=1,2
        do k=1,4
          potspin(1:nr,is,js) = potspin(1:nr,is,js) + pc2s(is,js,k)*potcart(1:nr,k)
        end do
      end do
    end do
    do is=1,4
      s1 = i2is(1,is); s2 = i2is(2,is)
      vdiff(1:nr,is) = potspin(1:nr,s1,s2)
    end do
!   --------------------------------------------------------------------
    do iq=1+sum(nalmsbden(1:ia2-1)),sum(nalmsbden(1:ia2))
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
      s1 = i2is(1,is); s2 = i2is(2,is)
!     spherical input potential
      if (ilm == 1) then
        work(1:nr) = suscbasis(1:nr,ib,ilm,is,ia2)*potspin(1:nr,s1,s2)
        norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
        vlmsbden(iq) = norm
        vdiff(1:nr,is) = vdiff(1:nr,is) - vlmsbden(iq)*suscbasis(1:nr,ib,ilm,is,ia2)
      end if
    end do
    do is=1,4
      work(1:nr) = abs(vdiff(1:nr,is))
      norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
      if (lhdio) write(iodb,'("build_vscfb: ia,is=",2i4,"  vdiff=",es8.1)') ia, is, real(norm)
    end do
!   test output
!    work = 0.d0
!    do iq=1+sum(nalmsbden(1:ia-1)),sum(nalmsbden(1:ia))
!      i4 = i2almsbden(:,iq)
!      ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
!      if (ilm == 1) then
!        work(1:nr) = work(1:nr) + vlmsbden(iq)*suscbasis(1:nr,ib,ilm,is,ia2)
!      end if
!    end do
!    do ir=1,nr
!      write(777,'(4es16.8)') rmesh(ir,ia), br(ir,ia), work(ir)
!    end do
! ******************
  end do
! ******************
! All done
  end subroutine build_vscfb2
