  subroutine integrated_charge(ia,gfsum,egfsum,tgfsum,charge,mspin,morb)
! Cell-integrated valence charge, spin and orbital moments
! Density matrix for LDA+U
  use global

  implicit none

! Which atom
  integer(kind=i4b), intent(in)  :: ia
! Matrix elements of energy integrated GF
  complex(kind=c8b), intent(in)  :: gfsum(nlmsb,nlmsb), egfsum(nlmsb,nlmsb), tgfsum(nlmsb,nlmsb)
! Valence charge, spin and orbital moments
  real(kind=r8b),    intent(out) :: charge, mspin(3), morb(3)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: small = 1.d-6
  complex(kind=c8b), parameter :: cone = (1.d0,0.d0), czero = (0.d0,0.d0)
  integer(kind=i4b) :: j, jlms, js, jlm, jl, jm, jb
  integer(kind=i4b) :: i, ilms, is, ilm, il, im, ib
  integer(kind=i4b) :: nr, k, klm, i2(2), i3(3)
  real(kind=r8b)    :: fac, norm, newlen
  character*14      :: filename
  complex(kind=c8b), allocatable :: rhoint(:,:,:), rhomat2(:,:)
  real(kind=r8b),    allocatable :: rhotot(:,:), rhoylm(:,:)
  integer(kind=i4b) :: info, lwork, ieval(nlms)
  complex(kind=c8b), allocatable :: work(:), vsocxy(:,:)
  real(kind=r8b),    allocatable :: rwork(:), evals(:)
  complex(kind=c8b) :: zetorque(3,0:nlmax,nsmax), znorm
  complex(kind=c8b), external :: radint


  allocate(rhoint(0:6,lmmax,lmmax),rhotot(0:6,0:nlmax),rhoylm(0:6,lmmax2),rhomat2(nlms,nlms),work(nrmax))
!  allocate(vsocx(nlmsba(ia),nlmsba(ia)),vsocy(nlmsba(ia),nlmsba(ia)))
! ----------------------------------------------------------------
! density matrix and its radial integral (COMPLEX)
  rhoint = 0.d0; rhomat(:,:,ia) = 0.d0
  ebandv(:,:,ia) = 0.d0; etorque(:,:,ia) = 0.d0; eldau(:,:,ia) = 0.d0
  zetorque = 0.d0
  nr = nrpts(ia)
  do j=1,nlmsba(ia)
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    jlms = lms2i(jlm,js)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      ilms = lms2i(ilm,is)
      do k=1,4
!        rhoint(mod(k,4),ilm,jlm) = rhoint(mod(k,4),ilm,jlm) + overlap(i,j,ia)*gfsum(i,j)*pauli(js,is,k)
        rhoint(mod(k,4),ilm,jlm) = rhoint(mod(k,4),ilm,jlm) + overlap(i,j,ia)*gfsum(i,j)*ds2c(k,js,is)  !  dot product: gfsum and ds2c as vectors
      end do
      if (ilm == jlm) then
        il = i2lm(2,ilm)
!       band energy, spin dn projection
        ebandv(il,1,ia)     = ebandv(il,1,ia)     + 0.5d0*overlap(i,j,ia)*egfsum(i,j)*(pauli(js,is,4) - sum(pauli(js,is,1:3)*magdir(:,ia)))
!       band energy, spin up projection
        ebandv(il,nsmax,ia) = ebandv(il,nsmax,ia) + 0.5d0*overlap(i,j,ia)*egfsum(i,j)*(pauli(js,is,4) + sum(pauli(js,is,1:3)*magdir(:,ia)))
!       torque between mspin and bxc
        work(1:nr) = br(1:nr,ia)*phiref(:,jb,i2lm(2,ilm),js,ia)*phiref(:,ib,i2lm(2,ilm),is,ia)
        znorm = radint(nr,work(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        do k=1,3
          etorque(k,il,ia) = etorque(k,il,ia) - znorm*pauli(js,is,k)*gfsum(i,j)
!          etorque(k,il,is,ia) = etorque(k,il,is,ia) + overlap(i,j,ia)*tgfsum(i,j,k)
!          zetorque(k,il,is) = zetorque(k,il,is) + tgfsum(i,j,k)
        end do
      end if
      rhomat(ilms,jlms,ia) = rhomat(ilms,jlms,ia) + gfsum(i,j)*overlap(i,j,ia)
    end do
  end do
! magnetic torque: transverse part
  do il=0,nlmax
    etorque(:,il,ia) = etorque(:,il,ia) - magdir(:,ia)*dot_product(magdir(:,ia),etorque(:,il,ia))
  end do
!  write(iodb,'(" band energy for ia=",i4,100f12.8)') ia, ebandv(:,:,ia)
!  write(iodb,'(" torque      for ia=",i4,100f12.8)') ia, zetorque(:,:,:)
!  etorque(:,:,:,ia) = real(zetorque)
! ----------------------------------------------------------------
! test: diagonalize rhomat
  if (lrhodiag) then
    lwork = 2*nlms
    allocate(work(lwork),rwork(2*lwork),evals(nlms))
    rhomat2 = rhomat(:,:,ia)
    call zheev('V','U',nlms,rhomat2,nlms,evals,work,lwork,rwork,info)
    if (info /= 0) stop 'integrated_charge: failure in zheev'
!   sort the eigenvalues according to ilms
    do ilms=1,nlms
      norm = 0.d0
!     search for largest element for ilms
      do jlms=1,nlms
        if (abs(rhomat2(ilms,jlms)) > norm) then
!         if this eigenvector was already picked skip it
          if (any(ieval(1:ilms-1) == jlms)) cycle
          norm = abs(rhomat2(ilms,jlms))
          ieval(ilms) = jlms
        end if
      end do
!     normalize the eigenvalues according to ilms
!      znorm = rhomat2(ilms,ieval(ilms))
!      znorm = conjg(znorm)/sqrt(abs(znorm))
      znorm = 1.d0/sqrt(dot_product(rhomat2(1:nlms,ieval(ilms)),rhomat2(1:nlms,ieval(ilms))))
      rhomat2(1:nlms,ieval(ilms)) = rhomat2(1:nlms,ieval(ilms))*znorm
    end do
    if (lhdio) write(iodb,'("Eigendecomposition of rhomat for ia=",i4)') ia
    do ilms=1,nlms
      if (lhdio) write(iodb,'(f8.4," | ",100f6.3)') evals(ieval(ilms)), rhomat2(1:nlms,ieval(ilms))
    end do
    deallocate(work,rwork,evals)
  end if
! ----------------------------------------------------------------
! LDA+U energy
  if (lldau .and. ildau(ia) == 1) then
    call zgemm('N','N',nlms,nlms,nlms,cone,rhomat(:,:,ia),nlms,rhomat(:,:,ia),nlms,czero,rhomat2,nlms)
    do ilms=1,nlms
      i2 = i2lms(:,ilms)
      ilm = i2(1); is = i2(2)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      eldau(il,is,ia) = eldau(il,is,ia) + 0.5d0*(ueff(il,ia) - jeff(il,ia))*(rhomat(ilms,ilms,ia) - rhomat2(ilms,ilms))
    end do
!    write(iodb,'(" U,J=",10f8.4)') ueff(:,ia), jeff(:,ia)
    if (lhdio) write(iodb,'(" LDA+U energy for ia=",i4,100f12.8)') ia, eldau(:,:,ia)
  end if
! ----------------------------------------------------------------
! output density matrix
!  write(iodb,'(" rhomat for ia=",i4)') ia
!  do ilms=1,nlms
!    write(iodb,'(1000f8.4)') real(rhomat(ilms,:,ia))
!  end do
!  write(iodb,'(" rhomat2 for ia=",i4)') ia
!  do ilms=1,nlms
!    write(iodb,'(1000f8.4)') real(rhomat2(ilms,:))
!  end do
! ----------------------------------------------------------------
! orbital densities and moments - only charge part (spin trace)
  do jlm=1,lmmax
    do ilm=1,lmmax
      do k=4,6
        do klm=1,lmmax
          rhoint(k,ilm,jlm) = rhoint(k,ilm,jlm) + lorb(ilm,klm,k-3)*rhoint(0,klm,jlm)
        end do
      end do
    end do
  end do
! ----------------------------------------------------------------
! this is the lm decomposition of the spherical part
  rhotot = 0.d0
  do ilm=1,lmmax
    i2 = i2lm(:,ilm)
    im = i2(1); il = i2(2)
    rhotot(:,il) = rhotot(:,il) + rhoint(:,ilm,ilm)
  end do
! ----------------------------------------------------------------
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
! ----------------------------------------------------------------
! total charge, spin and orbital moments
  charge = sum(rhotot(0,:))
  mspin  = sum(rhotot(1:3,:),dim=2)
  morb   = sum(rhotot(4:6,:),dim=2)
! ----------------------------------------------------------------
! magnitude and direction of spin moment
  newdir(:,ia) = mspin
  newlen = sqrt(dot_product(newdir(:,ia),newdir(:,ia)))
  if (newlen < small) then
    newdir(:,ia) = (/0.d0,0.d0,1.d0/)
  else
    newdir(:,ia) = newdir(:,ia)/newlen
  end if
!  newdir(:,ia) = newdir(:,ia)/(newlen + 1.d-12)
!  write(iodb,'(7f12.8)') mspin, newlen, newdir(:,ia)
! ----------------------------------------------------------------
! Write info to file
!  write(filename,'("charge",i4.4,".dat")') ia
!  open(file=filename,unit=iofile,status='replace')
  write(filename,'("chargespdf.dat")')
  if (ia == 1) then
    open(file=filename,unit=iofile,status='replace')
  else
    open(file=filename,unit=iofile,status='old',access='append')
  end if
  write(iofile,'("GS quantities for ia =",i4,": sum, spdf")') ia
  write(iofile,'("  qe =",100f16.8)') sum(rhotot(0,:)), rhotot(0,:)
  write(iofile,'("  msx=",100f16.8)') sum(rhotot(1,:)), rhotot(1,:)
  write(iofile,'("  msy=",100f16.8)') sum(rhotot(2,:)), rhotot(2,:)
  write(iofile,'("  msz=",100f16.8)') sum(rhotot(3,:)), rhotot(3,:)
  write(iofile,'("  mox=",100f16.8)') sum(rhotot(4,:)), rhotot(4,:)
  write(iofile,'("  moy=",100f16.8)') sum(rhotot(5,:)), rhotot(5,:)
  write(iofile,'("  moz=",100f16.8)') sum(rhotot(6,:)), rhotot(6,:)
! ----------------------------------------------------------------
  write(iofile,'("decomposition into spherical harmonics")')
  write(iofile,'("  qeylm=")')
  do ilm=1,lmmax2
    norm = rhoylm(0,ilm)
    if (abs(norm) > atol) then
      write(iofile,'(2i4,f16.8," | ",f16.8)') i2lm(:,ilm), norm, gs_qlm(ilm,ia)
    end if
    gs_qlm(ilm,ia) = norm
  end do
! ----------------------------------------------------------------
  write(iofile,'("  msylm=")')
  do ilm=1,lmmax2
    norm = sqrt(dot_product(rhoylm(1:3,ilm),rhoylm(1:3,ilm)))
    if (norm > atol) then
      write(iofile,'(2i4,4f16.8," | ",f16.8)') i2lm(:,ilm), norm, rhoylm(1:3,ilm), gs_mlm(ilm,ia)
    end if
    gs_mlm(ilm,ia) = norm
  end do
! ----------------------------------------------------------------
  write(iofile,'("  moylm=")')
  do ilm=1,lmmax2
    norm = sqrt(dot_product(rhoylm(4:6,ilm),rhoylm(4:6,ilm)))
    if (norm > atol) then
      write(iofile,'(2i4,8f16.8)') i2lm(:,ilm), rhoylm(4:6,ilm), norm
    end if
  end do
  close(iofile)
! ----------------------------------------------------------------
  deallocate(rhoint,rhotot,rhoylm,rhomat2)
! All done!
  end subroutine integrated_charge
