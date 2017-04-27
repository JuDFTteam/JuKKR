  subroutine anisotropic_jij(lldautoo)
! Magnetic exchange interactions
! Eq. 9 in PRB 79, 045209 (2009)
! One factor of 1/2 is taken into the definition of the magnetic Hamiltonian:
! H = -1/2 \sum_{ij} e_i J_{ij} e_j

  use global

  implicit none

! Whether to add LDA+U correction to spin splitting
  logical,  intent(in) :: lldautoo
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
  complex(kind=c8b), parameter :: i2pi = cmplx(0.d0,8.d0*atan(1.d0))
  real(kind=r8b),    parameter :: uz(3) = (/0.d0,0.d0,1.d0/)
! Which atoms for Jij's
  integer(kind=i4b) :: najij
  integer(kind=i4b) :: iajij(nasusc)
! Exchange interactions
  complex(kind=c8b), allocatable :: jijdos(:,:,:,:,:)
  real(kind=r8b),    allocatable :: jij(:,:,:,:), jiout(:), chi0(:,:,:,:), linres(:,:), rhs(:)
  real(kind=r8b),    allocatable :: bxcef(:,:), msef(:,:), dosef(:), msbxc(:,:), rotmat(:,:,:), dev(:,:)
! Storage for angular momentum matrices
  complex(kind=c8b), allocatable :: dti(:,:,:,:)
  complex(kind=c8b), allocatable :: gij(:,:), gji(:,:), work1(:,:), work2(:,:)
! ----------------------------------------------------------------------
  real(kind=r8b), save :: jiinp(1000) = 0.d0
  integer(kind=i4b) :: ie, i2(2), ia2, ia, ja2, ja, jalms, ialms, jlms, ilms, i, j, istart, iend, k
  complex(kind=c8b) :: de, trace, jijtemp(3,3), jitemp, ze0, dy, fac
  real(kind=r8b)    :: re, im, avgdev, scaling, jiavg, bin(3), bout(3), udir(3)
  character*17      :: filename
  integer(kind=i4b) :: info
  integer(kind=i4b), allocatable :: ipiv(:)

! Use ijij to find which atoms to compute Jij for
  najij = 0
  do ia=1,nasusc
    if (ijij(ia) /= 0) then
      najij = najij + 1
      iajij(najij) = ia
    end if
  end do
  write(*,'("anisotropic_jij: RAM=",f16.3," MB",/)') 16.d0*(najij+4)*nlmsb*nlmsb/1024.d0**2
  allocate(dti(nlmsb,nlmsb,6,najij),jij(3,3,najij,najij),jiout(najij),jijdos(3,3,najij,najij,nesusc))
  allocate(chi0(3,3,najij,najij),linres(2*najij,2*najij),rhs(2*najij),ipiv(2*najij))
  allocate(bxcef(3,najij),msef(3,najij),dosef(najij),msbxc(3,najij),rotmat(3,3,najij),dev(3,najij))
  allocate(gij(nlmsb,nlmsb),gji(nlmsb,nlmsb),work1(nlmsb,nlmsb),work2(nlmsb,nlmsb))
!  write(*,'("anisotropic_jij: RAM=",f16.3," MB")') 16.d0*(najij+4)*nlms*nlms/1024.d0**2
!  allocate(dti(nlms,nlms,3,najij),jij(3,3,najij,najij))
!  allocate(gij(nlms,nlms),gji(nlms,nlms),work1(nlms,nlms),work2(nlms,nlms))
  jij = czero; bxcef = 0.d0; msef = 0.d0; dosef = 0.d0; msbxc = 0.d0; chi0 = 0.d0
! --------------------------------------------------------------------
!                       t-matrix corrections
! --------------------------------------------------------------------
  write(*,'(" Magnetic exchange interactions: J onsite")')
  dti = czero
  do ia2=1,najij
    ia = iajij(ia2)
!   dti = Ri Bi Ri \vec\sigma
    call build_dti(1,ia,lldautoo,dti(:,:,:,ia2),jiout(ia2))
    write(*,'(" ia=",i4,"  Jiiout=",es16.8)') iajij(ia2), jiout(ia2)
!    do i=1,3
!      trace = czero
!      do ilms=1,nlms
!        trace = trace + dti(ilms,ilms,i,ia2)
!      end do
!      write(iodb,'(" tr dti=",2i4,2es16.8)') ie, ia, trace
!    end do
!   rotations to local frame
    if (lrot .and. ispinrot == 1) then
      call rotvec(uz,magdir(:,ia),rotmat(:,:,ia2))
    else
      call rotvec(uz,uz,rotmat(:,:,ia2))
    end if
!    do i=1,3
!      write(*,'("rotmat, ia2=",i4,3es16.8)') ia2, rotmat(i,:,ia2)
!    end do
  end do
! **************
! Energy loop
  do ie=1,nesusc
! **************
    if (ldos) then
      de = desusc(ie)
    else
      if (ie > nescf) then
        de = czero
      else
        de = descf(ie)
      end if
    end if
!   --------------------------------------------------------------------
!                         t-matrix corrections
!   --------------------------------------------------------------------
    do ia2=1,najij
      ia = iajij(ia2)
!     dti = Ri Bi Ri \vec\sigma
!      call build_dti_v2(ie,ia,lldautoo,dti(:,:,:,ia2))
!      do i=1,3
!        trace = czero
!        do ilms=1,nlms
!          trace = trace + dti(ilms,ilms,i,ia2)
!        end do
!        write(iodb,'(" tr dti=",2i4,2es16.8)') ie, ia, trace
!      end do
    end do
!   --------------------------------------------------------------------
!                           Jij for all pairs    
!   --------------------------------------------------------------------
    do ja2=1,najij
      ja = iajij(ja2)
      do ia2=1,najij
        ia = iajij(ia2)
!       collect the relevant blocks of the structural GF
!        gij = czero; gji = czero
!        do jlms=1,nlms
!          jalms = alms2i(jlms,ja)
!          do ilms=1,nlms
!            ialms = alms2i(ilms,ia)
!            gij(ilms,jlms) = gstruct(ialms,jalms,ie)
!            gji(jlms,ilms) = gstruct(jalms,ialms,ie)
!          end do
!        end do
        call projected_gf(ie,ia,ja,gij,.true.,.true.)
        call projected_gf(ie,ja,ia,gji,.true.,.true.)
!       ---------------------------
!        trace = czero
!        do ilms=1,nlms
!          trace = trace + gij(ilms,ilms)
!        end do
!        write(iodb,'(" tr gij=",3i4,2es16.8)') ie, ia, ja, trace
!        trace = czero
!        do ilms=1,nlms
!          trace = trace + gji(ilms,ilms)
!        end do
!        write(iodb,'(" tr gji=",3i4,2es16.8)') ie, ja, ia, trace
!       ---------------------------
!       cartesian components of Jij tensor
        jijtemp = czero
        do j=1,3
          do i=1,3
!           ------------------------------------------------------------
!!           dti * gij --> work1
!            call zgemm('N','N',nlms,nlms,nlms,cone,dti(:,:,i,ia2),nlms,gij,nlms,czero,work1,nlms)
!!           (dti * gij) * dtj --> work2
!            call zgemm('N','N',nlms,nlms,nlms,cone,work1,nlms,dti(:,:,j,ja2),nlms,czero,work2,nlms)
!!           de * (dti * gij * dtj) * gji --> work1
!            call zgemm('N','N',nlms,nlms,nlms,descf(ie),work2,nlms,gji,nlms,czero,work1,nlms)
!            trace = czero
!            do ilms=1,nlms
!              trace = trace + work1(ilms,ilms)
!            end do
!            jij(i,j,ia2,ja2) = jij(i,j,ia2,ja2) + 4.d0*trace/i2pi
!!           dti * gij --> work1
!            call zgemm('C','C',nlms,nlms,nlms,cone,dti(:,:,i,ia2),nlms,gji,nlms,czero,work1,nlms)
!!           (dti * gij) * dtj --> work2
!            call zgemm('N','C',nlms,nlms,nlms,cone,work1,nlms,dti(:,:,j,ja2),nlms,czero,work2,nlms)
!!           de * (dti * gij * dtj) * gji --> work1
!            call zgemm('N','C',nlms,nlms,nlms,conjg(descf(ie)),work2,nlms,gij,nlms,czero,work1,nlms)
!            trace = czero
!            do ilms=1,nlms
!              trace = trace + work1(ilms,ilms)
!            end do
!            jij(i,j,ia2,ja2) = jij(i,j,ia2,ja2) - 4.d0*trace/i2pi
!           ------------------------------------------------------------
!           Generalized Lichtenstein formula
!           dVi * Gij(E+i0) --> work1
            call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ia),cone,dti(:,:,i,ia2),nlmsb,gij,nlmsb,czero,work1,nlmsb)
!           (dVi * Gij(E+i0)) * dVj --> work2
            call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,work1,nlmsb,dti(:,:,j,ja2),nlmsb,czero,work2,nlmsb)
!           (dVi * Gij(E+i0) * dVj) * Gji(E+i0) --> work1
            call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ja),cone,work2,nlmsb,gji,nlmsb,czero,work1,nlmsb)
!           onsite correction -- REVISE !!!
            if (ljionsite .and. ia2 == ja2 .and. i == j) then
!             - dVi * Gii(E+i0) --> work1
              do k=1,3
                call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ia),-cone*magdir(k,ia),dti(:,:,k,ia2),nlmsb,gij,nlmsb,cone,work1,nlmsb)
              end do
            end if
            trace = czero
            do ilms=1,nlmsba(ia)
              trace = trace + work1(ilms,ilms)
            end do
            jij(i,j,ia2,ja2) = jij(i,j,ia2,ja2) + real(2.d0*de*trace/i2pi)
!           spectral density for Jij
            re = real(trace); im = aimag(trace)
            if (abs(re) < 1.d-10) re = 0.d0
            if (abs(im) < 1.d-10) im = 0.d0
            jijtemp(i,j) = cmplx(re,im)
!           ------------------------------------------------------------
!           Spherical average of KS susceptibility
!           pauli_i * Gij(E+i0) --> work1
            call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ia),cone,dti(:,:,i+3,ia2),nlmsb,gij,nlmsb,czero,work1,nlmsb)
!           (pauli_i * Gij(E+i0)) * pauli_j --> work2
            call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,work1,nlmsb,dti(:,:,j+3,ja2),nlmsb,czero,work2,nlmsb)
!           (pauli_i * Gij(E+i0) * pauli_j) * Gji(E+i0) --> work1
            call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ja),cone,work2,nlmsb,gji,nlmsb,czero,work1,nlmsb)
            trace = czero
            do ilms=1,nlmsba(ia)
              trace = trace + work1(ilms,ilms)
            end do
            chi0(i,j,ia2,ja2) = chi0(i,j,ia2,ja2) - real(2.d0*de*trace/i2pi)
!           ------------------------------------------------------------
!           dVi * Gij(E-i0)  = dVi * Gji(E+i0)^\dagger--> work1
!            call zgemm('N','C',nlmsba(ia),nlmsba(ja),nlmsba(ia),cone,dti(:,:,i,ia2),nlmsb,gji,nlmsb,czero,work1,nlmsb)
!           (dVi * Gij(E-i0)) * dVj = (dVi * Gji(E+i0)^\dagger) * dVj --> work2
!            call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,work1,nlmsb,dti(:,:,j,ja2),nlmsb,czero,work2,nlmsb)
!           (dVi * Gij(E-i0) * dVj) * Gji(E-i0) = (dVi * Gji(E+i0)^\dagger * dVj) * Gij(E+i0)^\dagger --> work1
!            call zgemm('N','C',nlmsba(ia),nlmsba(ia),nlmsba(ja),cone,work2,nlmsb,gij,nlmsb,czero,work1,nlmsb)
!            trace = czero
!            do ilms=1,nlmsba(ia)
!              trace = trace + work1(ilms,ilms)
!            end do
!            jij(i,j,ia2,ja2) = jij(i,j,ia2,ja2) - real(conjg(de)*trace/i2pi)
!           ------------------------------------------------------------
          end do  ! i loop
!         ------------------------------------------------------------
!         DOS at EF from atoms in Jij <=> static susceptibility sum rule
!         Gji(E+i0) * Gij(E+i0) --> work1
          call zgemm('N','N',nlmsba(ja),nlmsba(ja),nlmsba(ia),cone,gji,nlmsb,gij,nlmsb,czero,work2,nlmsb)
          trace = czero
          do ilms=1,nlmsba(ja)
            trace = trace + work2(ilms,ilms)
          end do
          dosef(ja2) = dosef(ja2) + real(2.d0*de*trace/i2pi)/3.d0  ! it's inside the j loop
!         averaged Bxc at EF from atoms in Jij
!         dVj * ( Gji(E+i0) * Gij(E+i0)) --> work1
          call zgemm('N','N',nlmsba(ja),nlmsba(ja),nlmsba(ja),cone,dti(:,:,j,ja2),nlmsb,work2,nlmsb,czero,work1,nlmsb)
          trace = czero
          do ilms=1,nlmsba(ja)
            trace = trace + work1(ilms,ilms)
          end do
          bxcef(j,ja2) = bxcef(j,ja2) + real(2.d0*de*trace/i2pi)
!         averaged mspin at EF from atoms in Jij
!         pauli_j * ( Gji(E+i0) * Gij(E+i0)) --> work1
          call zgemm('N','N',nlmsba(ja),nlmsba(ja),nlmsba(ja),cone,dti(:,:,j+3,ja2),nlmsb,work2,nlmsb,czero,work1,nlmsb)
          trace = czero
          do ilms=1,nlmsba(ja)
            trace = trace + work1(ilms,ilms)
          end do
          msef(j,ja2) = msef(j,ja2) - real(2.d0*de*trace/i2pi)
!         averaged ms * Bxc
          if (ia2 == ja2) then
!           dVj * Gii(E+i0) --> work1
            call zgemm('N','N',nlmsba(ja),nlmsba(ja),nlmsba(ja),cone,dti(:,:,j,ja2),nlmsb,gij,nlmsb,czero,work1,nlmsb)
            trace = czero
            do ilms=1,nlmsba(ja)
              trace = trace + work1(ilms,ilms)
            end do
            msbxc(j,ja2) = msbxc(j,ja2) + real(2.d0*de*trace/i2pi)
          end if
!         ------------------------------------------------------------
        end do  ! j loop
!        jij(:,:,ia2,ja2) = jij(:,:,ia2,ja2) + (jijtemp - conjg(transpose(jijtemp)))/i2pi
!       ---------------------------
        jijdos(:,:,ia2,ja2,ie) = jijtemp
!       ---------------------------------------------------------------- 
      end do  ! ia2 loop
!     Fermi surface contributions
!      if (ie == nesusc) then
!       averaged Bxc at EF
!        do i=1,3
!         dVi * Gij(E+i0) --> work1
!          call zgemm('N','N',nlmsba(ja),nlmsba(ja),nlmsba(ja),cone,dti(:,:,i,ja2),nlmsb,gij,nlmsb,czero,work1,nlmsb)
!          trace = czero
!          do ilms=1,nlmsba(ja)
!            trace = trace + work1(ilms,ilms)
!          end do
!          bxcef(i,ja2) = -real(2.d0*trace/i2pi)
!        end do
!       DOS at EF
!        trace = czero
!        do ilms=1,nlmsba(ja)
!          trace = trace + gij(ilms,ilms)
!        end do
!        dosef(ja2) = -real(2.d0*trace/i2pi)
!      end if
!     ---------------------------------------------------------------- 
    end do  ! ja2 loop
!   --------------------------------------------------------------------
! ******
  end do  ! ie loop
! ******
! Rotate onsite contributions to local frame
  do ja2=1,najij
    bxcef(:,ja2) = matmul(transpose(rotmat(:,:,ja2)),bxcef(:,ja2))
    msef(:,ja2)  = matmul(transpose(rotmat(:,:,ja2)),msef(:,ja2))
    msbxc(:,ja2) = matmul(transpose(rotmat(:,:,ja2)),msbxc(:,ja2))
  end do
! Transform the tensors too
  write(*,'(/," Magnetic exchange interactions: J tensor")')
  do ja2=1,najij
    jiavg = 0.d0
    do ia2=1,najij
!     rotate to local frame
!     R_i^T J_ij R_j
      jij(:,:,ia2,ja2) = matmul(transpose(rotmat(:,:,ia2)),matmul(jij(:,:,ia2,ja2),rotmat(:,:,ja2)))
      chi0(:,:,ia2,ja2) = matmul(transpose(rotmat(:,:,ia2)),matmul(chi0(:,:,ia2,ja2),rotmat(:,:,ja2)))
!     ------------------------------------------------------------------
!     correction from EF due to constant Ne
      if (ljijef) then
        do j=1,3
          do i=1,3
            jij(i,j,ia2,ja2)  = jij(i,j,ia2,ja2)  - bxcef(i,ia2)*bxcef(j,ja2)/sum(dosef)
            chi0(i,j,ia2,ja2) = chi0(i,j,ia2,ja2) + msef(i,ia2)*msef(j,ja2)/sum(dosef)
          end do
        end do
      end if
!     ------------------------------------------------------------------
      jij(:,:,ia2,ja2) = 2.d0*jij(:,:,ia2,ja2)
!     chop small numbers
      where (abs(jij(:,:,ia2,ja2)) < 1.d-10) jij(:,:,ia2,ja2) = 0.d0
!      jiavg = jiavg + 0.5d0*(jij(1,1,ia2,ja2) + jij(2,2,ia2,ja2) + jij(1,2,ia2,ja2) - jij(2,1,ia2,ja2))
      jiavg = jiavg + 0.5d0*(jij(1,1,ia2,ja2) + jij(2,2,ia2,ja2))
      write(*,'(" ia,ja=",2i4)') iajij(ia2), iajij(ja2)
      do i=1,3
        write(*,'(100es16.8)') jij(i,1:3,ia2,ja2)
      end do
    end do
    write(*,'(" ja2=",i4,"  dosef=",es16.8," bxcef=",3es16.8," msef=",3es16.8)') ja2, dosef(ja2), bxcef(:,ja2), msef(:,ja2)
    write(*,'(" ja= ",i4,"  Ji=   ",es16.8," msbxc=",3es16.8)') iajij(ja2), jiavg, msbxc(:,ja2)
    jiinp(ja2) = jiavg
    jiout(ja2) = jiavg
  end do
! **********************************************************************
! Constraining fields
  if (lrot .and. ispinrot == 1 .and. lbconst) then
    linres(:,:) = 0.d0
    do ja2=1,najij
    do j=1,2
      ja = iajij(ja2)
      do ia2=1,najij
      do i=1,2
!       keep transverse part only
        linres(i+2*(ia2-1),j+2*(ja2-1)) = chi0(i,j,ia2,ja2)
      end do
      end do
!     M_j R_j^T e_j^out
      rhs(j+2*(ja2-1)) = gs_mlm(1,ja)*dot_product(rotmat(:,j,ja2),newdir(:,ja))
    end do
    end do
    call dgesv(2*najij,1,linres,2*najij,ipiv,rhs,2*najij,info)
    if (info /= 0) stop 'anisotropic_jij: failure in dgesv'
!   rotate back to global frame: dev = R.rhs
    do ia2=1,najij
      dev(:,ia2) = rotmat(:,1,ia2)*rhs(2*ia2-1) + rotmat(:,2,ia2)*rhs(2*ia2)
    end do
!   keep deviations small
!    avgdev = sqrt(sum(dev(:,:)*dev(:,:)))/najij
!    if (avgdev < 1.d-1) then
!      scaling = 1.d0
!    else
!      scaling = 1.d-1/avgdev
!    end if
!    dev(:,:) = scaling*dev(:,:)
    where (abs(dev) < 1.d-12) dev = 0.d0
!   new spin direction
    write(*,'("Constraining fields:")')
    do ia2=1,najij
      ia = iajij(ia2)
      bin(:) = bconlen(ia)*bcondir(:,ia)
      bout(:) = dirmix2*dev(:,ia)
!     make sure the field is perpendicular to the target spin direction
      udir(:) = magdir(:,ia)
!      udir(:) = newdir(:,ia)
!      udir(:) = magdir(:,ia) + newdir(:,ia)
!      udir(:) = udir(:)/sqrt(dot_product(udir,udir))
      bout(:) = bout(:) - udir(:)*dot_product(udir(:),bout(:))
      bout(:) = bout(:) + bin(:)
      bconlen(ia) = sqrt(dot_product(bout,bout))
      if (bconlen(ia) > 1.d-8) then
        bcondir(:,ia) = bout(:)/bconlen(ia)
      else
        bconlen(ia) = 0.d0
        bcondir(:,ia) = (/0.d0,0.d0,1.d0/)
      end if
      write(*,'("  ia2=",i4,"  bconlen, bcondir=",4f16.8,"  b.e=",f16.8)') ia2, bconlen(ia), bcondir(:,ia), dot_product(bcondir(:,ia),newdir(:,ia))
    end do
  end if
! **********************************************************************
! Linear response for spin directions
  if (lrot .and. ispinrot == 1 .and. .not.lbconst) then
    linres(:,:) = 0.d0
    do ja2=1,najij
    do j=1,2
      ja = iajij(ja2)
      do ia2=1,najij
      do i=1,2
!       keep transverse part only
        linres(i+2*(ia2-1),j+2*(ja2-1)) = -jij(i,j,ia2,ja2)
      end do
      end do
!     first value
!      if (abs(jiinp(ja2)) < 1.d-8) jiinp(ja2) = jiout(ja2)
      linres(j+2*(ja2-1),j+2*(ja2-1)) = linres(j+2*(ja2-1),j+2*(ja2-1)) + 1.01d0*jiinp(ja2)  ! avoid a singular matrix
!     R_j^T e_j^out
      rhs(j+2*(ja2-1)) = dirmix2*jiout(ja2)*dot_product(rotmat(:,j,ja2),newdir(:,ja))
!     save current value for next iteration
!      jiinp(ja2) = jiout(ja2)
    end do
    end do
    call dgesv(2*najij,1,linres,2*najij,ipiv,rhs,2*najij,info)
    if (info /= 0) stop 'anisotropic_jij: failure in dgesv'
!   rotate back to global frame: dev = R.rhs
    do ia2=1,najij
      dev(:,ia2) = rotmat(:,1,ia2)*rhs(2*ia2-1) + rotmat(:,2,ia2)*rhs(2*ia2)
    end do
!   keep deviations small
    avgdev = sqrt(sum(dev(:,:)*dev(:,:)))/najij
    if (avgdev < 1.d-1) then
      scaling = 1.d0
    else
      scaling = 1.d-1/avgdev
    end if
    dev(:,:) = scaling*dev(:,:)
    where (abs(dev) < 1.d-6) dev = 0.d0
!   new spin direction
    write(*,'("Deviations in orientations:")')
    do ia2=1,najij
      ia = iajij(ia2)
      magdir(:,ia) = magdir(:,ia) + dev(:,ia2)
      magdir(:,ia) = magdir(:,ia)/sqrt(dot_product(magdir(:,ia),magdir(:,ia)))
      write(*,'("  ia2=",i4,"  dev=",3f16.8,"  magdir=",3f16.8)') ia2, dev(:,ia2), magdir(:,ia)
    end do
  end if
! **********************************************************************
! --------------
! Print Jij DOS
  if (ldos) then
! --------------
    do ja2=1,najij
      do ia2=1,najij
        write(filename,'("ia",i4.4,"ja",i4.4,".jije")') iajij(ia2), iajij(ja2)
        open(file=filename,unit=iofile,status='replace')
        write(iofile,'("# Energy then (Jij(i,j),j=1,3),i=1,3)")')
        do ie=1,nesusc
!         analytical continuation using rational function
          if (ldosacon) then
            ze0 = real(esusc(ie))
            if (ie < nedosacon/2) then
              istart = 1
              iend   = nedosacon
            else if (ie + nedosacon/2 + mod(nedosacon,2) > nesusc) then
              istart = nesusc - nedosacon + 1
              iend   = nesusc
            else
              istart = ie - nedosacon/2 + 1
              iend   = ie + nedosacon/2 + mod(nedosacon,2)
            end if
            do j=1,3
              do i=1,3
                call zratint(nedosacon,esusc(istart:iend),jijdos(i,j,ia2,ja2,istart:iend),ze0,jijtemp(i,j),dy)
              end do
            end do
            write(iofile,'(1000es16.8)') ze0, ((jijtemp(i,j),j=1,3),i=1,3)
!         computed values for complex energy
          else
            write(iofile,'(1000es16.8)') esusc(ie), ((jijdos(i,j,ia2,ja2,ie),j=1,3),i=1,3)
          end if
        end do
        close(iofile)
      end do
    end do
! ------
  end if
! ------
  deallocate(jij,jiout,jijdos,work1,work2,dti,gij,gji,chi0,linres,ipiv,rhs,bxcef,msef,dosef,msbxc,rotmat)
! All done!
  end subroutine anisotropic_jij
