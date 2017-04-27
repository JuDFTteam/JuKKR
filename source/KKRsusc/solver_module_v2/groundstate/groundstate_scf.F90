  subroutine groundstate_scf()
! Ground state properties without fitting GFs
! The matrix elements of the operators are filled in subroutine observables
  use global

  implicit none

! temporary parameters
!  integer(kind=i4b), parameter :: ne = 1001, nesusc = 300
! i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
!  real(kind=r8b),    parameter :: omegamax = 4.d-1
! auxiliary
  complex(kind=c8b) :: de, work(nrmax), ze, za, zb, tmp, block(4,4)
  complex(kind=c8b) :: gf(nlmsb,nlmsb), gfsum(nlmsb,nlmsb)
  integer(kind=i4b) :: nr, imap(nlmsb,nlmsb)
  real(kind=r8b)    :: r(nrmax), zre, zim, norm, start, finish
! moments of charge density, spin and orbital moments
  real(kind=r8b)    :: phi1(nrmax), phi2(nrmax)
  complex(kind=c8b) :: rho(nrmax,0:6,lmmax,lmmax)
  complex(kind=c8b) :: rhoint(0:6,lmmax,lmmax)
  complex(kind=c8b) :: zdos(lmmax,lmmax,nsmax,nsmax,nasusc,nedos)
  complex(kind=c8b) :: zdosl(0:nlmax,nsmax,nasusc,nedos)
  real(kind=r8b)    :: rhotot(0:6,0:nlmax)
  real(kind=r8b)    :: rhoylm(0:6,lmmax2)
  real(kind=r8b)    :: new_len, new_dir(3), fac
  character*14      :: filename
! looping
  integer(kind=i4b) :: ia, ja, ie, k, l, ir, klm
! decoding
  integer(kind=i4b) :: i2(2), i3(3), dsusc
  integer(kind=i4b) :: j, jlms, js, jlm, jl, jm, jb, jalmsb, jq
  integer(kind=i4b) :: i, ilms, is, ilm, il, im, ib, ialmsb, iq
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call cpu_time(start)
  nrv = 0.d0; zdosl = 0.d0; zdos = 0.d0
  do ia=1,nasusc
    nr = nrpts(ia)
    gfsum = 0.d0
! use the energy mesh from SCF calculation here
! construct the energy integrated site-diagonal GF
    do ie=1,nescf
      de = descf(ie)
!     ****  site diagonal part for local quantities  ****
      call projected_gf(ie+nesusc-nescf,ia,ia,gf,.true.,.true.)
!       (G^dagger - G)/(i 2pi)
      gf = (conjg(transpose(de*gf)) - de*gf)/i2pi
      gfsum = gfsum + gf
!     ***************************************************
    end do
! ================================================================
! Compute DOS
    if (idos == 1) then
!   ===================
    za = cmplx(e0dos,eimdos); zb = cmplx(e1dos,eimdos)
    do ie=1,nesusc
      ze = esusc(ie)
      call projected_gf(ie,ia,ia,gf,.true.,.true.)
!     DOS
      do j=1,nlmsba(ia)
        i3 = i2lmsb(:,j,ia)
        jb = i3(1); jlm = i3(2); js = i3(3)
        i2 = i2lm(:,jlm)
        jm = i2(1); jl = i2(2)
        phi1(1:nr) = phiref(1:nr,jb,jl,js,ia)
        do i=1,nlmsba(ia)
          i3 = i2lmsb(:,i,ia)
          ib = i3(1); ilm = i3(2); is = i3(3)
          i2 = i2lm(:,ilm)
          im = i2(1); il = i2(2)
          zdos(ilm,jlm,is,js,ia,ie) = zdos(ilm,jlm,is,js,ia,ie) + gf(i,j)*overlap(i,j,ia)
          if (ilm == jlm .and. is == js) zdosl(il,is,ia,ie) = zdosl(il,is,ia,ie) + gf(i,j)*overlap(i,j,ia)
        end do
      end do
    end do
! DOS to file
    write(filename,'("ldos",i4.4,".dat")') ia
    open(file=filename,unit=iofile,status='replace')
    do ie=1,nesusc
      ze = esusc(ie)
      write(iofile,'(100es16.8)') ze, sum(zdosl(:,1,ia,ie)), sum(zdosl(:,2,ia,ie)), zdosl(:,:,ia,ie)
    end do
    close(iofile)
    write(filename,'("lmdos",i4.4,".dat")') ia
    open(file=filename,unit=iofile,status='replace')
    do ie=1,nesusc
      ze = esusc(ie)
      write(iofile,'(100es16.8)') ze, ((zdos(ilm,ilm,is,is,ia,ie),ilm=1,lmmax),is=1,nsmax)
    end do
    close(iofile)
!   ======
    end if
! ================================================================
! filter small elements of the GF away
    imap = 0
    do j=1,nlmsba(ia)
      do i=1,nlmsba(ia)
        zre = real(gfsum(i,j)); zim = aimag(gfsum(i,j))
        if (abs(zre) < gstol) zre = 0.d0
        if (abs(zim) < gstol) zim = 0.d0
        gfsum(i,j) = cmplx(zre,zim)
        if (abs(gfsum(i,j)) > gstol) imap(i,j) = 1
      end do
    end do
!    write(*,*) "before density matrix"
! ----------------------------------------------------------------
! xc-potential and energy
!    call overlaps_susc(ia,gfsum)
!    call get_xc(ia,lmmax4,gfsum,kxclm)
!    call get_rhoden(ia,lmmax2,gfsum,rhoden)
!    call get_xc2(ia,lmmax4,kxclm)
!    call get_xc_basis(ia,lmmax2,gfsum,kxc,kxcbasis(:,:,ia))
! ----------------------------------------------------------------
! density matrix and its radial integral (COMPLEX)
    rho = 0.d0; rhoint = 0.d0; rhomat(:,:,ia) = 0.d0
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
          rho(1:nr,mod(k,4),ilm,jlm) = rho(1:nr,mod(k,4),ilm,jlm)   &
                                       + gfsum(i,j)*pauli(js,is,k)*phi1(1:nr)*phi2(1:nr)
!          rhoint(mod(k,4),ilm,jlm)   = rhoint(mod(k,4),ilm,jlm)     &
!                                       + overlap(i,j,ia)*gfsum(i,j)*pauli(js,is,k)
          rhoint(mod(k,4),ilm,jlm)   = rhoint(mod(k,4),ilm,jlm)     &
                                       + overlap(i,j,ia)*gfsum(i,j)*ds2c(k,js,is)
        end do
        rhomat(ilms,jlms,ia) = rhomat(ilms,jlms,ia) + gfsum(i,j)*overlap(i,j,ia)
      end do
    end do
! extra: orbital densities and moments - only charge part (spin trace)
! ----------------------------------------------------------------
    do jlm=1,lmmax
    do ilm=1,lmmax
      do k=4,6
      do klm=1,lmmax
        rho(1:nr,k,ilm,jlm) = rho(1:nr,k,ilm,jlm) + lorb(ilm,klm,k-3)*rho(1:nr,0,klm,jlm)
        rhoint(k,ilm,jlm)   = rhoint(k,ilm,jlm)   + lorb(ilm,klm,k-3)*rhoint(0,klm,jlm)
      end do
      end do
    end do
    end do
!   ----------------------------------------------------------------
!    do i=1,nlmsb
!      write(*,'(1000i2)') imap(i,1:nlmsb)
!    end do
!    do i=1,nlmsb
!      write(*,'("gfsum ",1000es10.1)') gfsum(i,1:nlmsb)
!    end do
! this is the lm decomposition of the spherical part
    rhotot = 0.d0
    do ilm=1,lmmax
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      do k=0,6
        rhotot(k,il) = rhotot(k,il) + rhoint(k,ilm,ilm)
      end do
    end do
! this is the spherical harmonics resummation
    rhoylm = 0.d0
    do klm=1,lmmax2
      do jlm=1,lmmax
      do ilm=1,lmmax
        norm = rgaunt(ilm,jlm,klm)
        if (abs(norm) > ylmtol) then
          rhoylm(:,klm) = rhoylm(:,klm) + rhoint(:,ilm,jlm)*norm
        end if
      end do
      end do
    end do
! and the winner is...
    write(*,'("GS quantities for ia =",i4,": sum, spdf")') ia
    write(*,'("  qe =",100f16.8)') sum(rhotot(0,:)), rhotot(0,:)
    write(*,'("  msx=",100f16.8)') sum(rhotot(1,:)), rhotot(1,:)
    write(*,'("  msy=",100f16.8)') sum(rhotot(2,:)), rhotot(2,:)
    write(*,'("  msz=",100f16.8)') sum(rhotot(3,:)), rhotot(3,:)
    write(*,'("  mox=",100f16.8)') sum(rhotot(4,:)), rhotot(4,:)
    write(*,'("  moy=",100f16.8)') sum(rhotot(5,:)), rhotot(5,:)
    write(*,'("  moz=",100f16.8)') sum(rhotot(6,:)), rhotot(6,:)
! output density matrix
    write(iodb,'(" density matrix for ia=",i4)') ia
    do ilms=1,nlms
      write(iodb,'(1000f8.4)') rhomat(ilms,:,ia)
    end do
! magnitude and direction of spin moment
    new_dir = sum(rhotot(1:3,:),dim=2)
    new_len = sqrt(dot_product(new_dir,new_dir))
    new_dir = new_dir/(new_len + 1.d-12)
    write(*,'("  new_len=",f16.8,",  new_dir=",3f16.8)') new_len, new_dir
    if (lrot) then
      if (ispinrot == 0) then
        magdir(:,ia) = new_dir
      else
        new_dir = (1.d0 - dirmix)*magdir(:,ia) + dirmix*new_dir
        magdir(:,ia) = new_dir/sqrt(dot_product(new_dir,new_dir))
      end if
    end if
! LDA+U energy: U/2 Tr rho (1 - rho)
    norm = 0.d0
    do ilms=1,nlms
      i2 = i2lms(:,ilms)
      norm = norm + 0.5d0*(ueff(i2(2),ia) - jeff(i2(2),ia))*rhomat(ilms,ilms,ia)*(1.d0 - rhomat(ilms,ilms,ia))
    end do
    write(*,'("U=",4f8.4,"  J=",4f8.4)') ueff(:,ia), jeff(:,ia)
    write(*,'(/," LDA+U energy: ia=",i4," E=",f16.8," Ry = ",f16.8," eV")') ia, real(norm), real(norm)*13.61d0
! output decomposition
    write(iodb,'("decomposition into spherical harmonics")')
    write(iodb,'("  qeylm=")')
    do ilm=1,lmmax2
      norm = rhoylm(0,ilm)
      if (abs(norm) > atol) then
        write(iodb,'(2i4,f16.8," | ",f16.8)') i2lm(:,ilm), norm, gs_qlm(ilm,ia)
      end if
    end do
    write(iodb,'("  msylm=")')
    do ilm=1,lmmax2
      norm = sqrt(dot_product(rhoylm(1:3,ilm),rhoylm(1:3,ilm)))
      if (norm > atol) then
        write(iodb,'(2i4,4f16.8," | ",f16.8)') i2lm(:,ilm), norm, rhoylm(1:3,ilm), gs_mlm(ilm,ia)
      end if
    end do
    write(iodb,'("  moylm=")')
    do ilm=1,lmmax2
      norm = sqrt(dot_product(rhoylm(4:6,ilm),rhoylm(4:6,ilm)))
      if (norm > atol) then
        write(iodb,'(2i4,8f16.8)') i2lm(:,ilm), norm, rhoylm(4:6,ilm)
      end if
    end do
! fill up new rho2ns
    new_rho2ns(:,:,:,ia) = 0.d0
    if (.not.lrot) new_dir = (/0.d0,0.d0,1.d0/)
    do klm=1,lmmax2
      do jlm=1,lmmax
      do ilm=1,lmmax
        fac = rgaunt(ilm,jlm,klm)
        if (abs(fac) > ylmtol) then
!       charge density
          new_rho2ns(1:nr,klm,1,ia) = new_rho2ns(1:nr,klm,1,ia) + fac*rho(1:nr,0,ilm,jlm)
!       magnetisation density
          do i=1,3
            new_rho2ns(1:nr,klm,2,ia) = new_rho2ns(1:nr,klm,2,ia) + fac*new_dir(i)*rho(1:nr,i,ilm,jlm)
          end do
        end if
      end do
      end do
    end do
! save spherical part of valence densities
    open(100+ia)
      nr = nrpts(ia)
      r = rmesh(1:nr,ia)
      do ilm=1,lmmax
        nrv(1:nr,:,ia) = nrv(1:nr,:,ia) + real(rho(1:nr,:,ilm,ilm))
      end do
      do ir=1,nr
        write(100+ia,'(8es16.8)') r(ir), nrv(ir,:,ia)!/r(ir)**2,
        write(200+ia,'(8es16.8)') r(ir), old_rho2ns(ir,7,:,ia), new_rho2ns(ir,7,:,ia)
      end do
    close(100+ia)
  end do
  call cpu_time(finish)
  write(*,'(/," GS quantities  time=",f10.3," s",/)') finish - start
! All done
  end subroutine groundstate_scf

