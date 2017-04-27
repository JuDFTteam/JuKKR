  subroutine density_matrix(onsite,struct)
  use global

  implicit none

! Onsite GF
  logical,           intent(in) :: onsite
! Structural GF
  logical,           intent(in) :: struct
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: pi = 4.d0*atan(1.d0)
  complex(kind=c8b), parameter :: cone = (1.d0,0.d0), czero = (0.d0,0.d0)
  complex(kind=c8b), allocatable :: ze(:), gf(:,:), zdmat(:,:,:)
  complex(kind=c8b) :: za, zb, ze0, dy, zms(1:4), zmo(1:3), zmso(1:3,1:3)
  character*14      :: filename
  integer(kind=i4b) :: ie, ne, istart, iend, natoms, ndim
  integer(kind=i4b) :: i2(2), i3(3)
  integer(kind=i4b), allocatable :: i2dmalms(:,:), ipiv(:)
  integer(kind=i4b) :: j, jlmsb, ja, jlms, js, jlm, jl, jm, jb
  integer(kind=i4b) :: i, ilmsb, ia, ilms, is, ilm, il, im, ib
  integer(kind=i4b) :: q, qlmsb, qa, qlms, qs, qlm, ql, qm, qb
  integer(kind=i4b) :: p, plmsb, pa, plms, ps, plm, pl, pm, pb
  real(kind=r8b),    allocatable :: pzef(:,:), rotmat(:,:,:,:)
  real(kind=r8b)    :: norm, u0(3), u1(3)
  integer(kind=i4b) :: info
  complex(kind=c8b), allocatable :: work(:,:)


! Energy mesh for DOS
  if (idos == 1) then
    ne = nesusc
    allocate(ze(ne))
    ze = esusc
  else
    ne = nedos
    allocate(ze(ne))
    za = dose0; zb = dose1
    do ie=1,nedos
      ze(ie) = za + (zb-za)*(ie-1)/(ne-1.d0)
    end do
  end if
! Atoms for output: determine dimension
  natoms = 0; ndim = 0
  do ia=1,nasusc
    if (iadmat(ia) /= 0) then
      natoms = natoms + 1
      do ilmsb=1,nlmsba(ia)
        i3 = i2lmsb(:,ilmsb,ia)
        ib = i3(1); il = i2lm(2,i3(2)); is = i3(3)
        if (ib == 1 .and. ildmat(il,is,ia) /= 0) ndim = ndim + 1
      end do
    end if
  end do
  write(*,'(/,"density_matrix: natoms=",i4," matrix dimension=",i8)') natoms, ndim
! Storage
  allocate(zdmat(ndim,ndim,ne),work(ndim,ndim),ipiv(ndim),gf(nlmsb,nlmsb))
  allocate(i2dmalms(2,ndim),pzef(nbmax,ndim),rotmat(lmmax,lmmax,natoms,natoms))
! Closest point to EF
  ie = minloc(abs(ze(:)-dosezero),dim=1)
  write(*,'("density_matrix: EF near ie=",i4,"  E=",2f12.8)') ie, ze(ie)
! Pointers and projection wfn coefficients
  i = 0; pzef(:,:) = czero
  do ia=1,nasusc
    if (iadmat(ia) /= 0) then
      do ilmsb=1,nlmsba(ia)
        i3 = i2lmsb(:,ilmsb,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        ilms = lms2i(ilm,is); il = i2lm(2,ilm)
        if (ildmat(il,is,ia) /= 0) then
          if (ib == 1) then
            i = i + 1
            i2dmalms(:,i) = (/ ilms, ia /)
          end if
!         right regular wfn at EF
          pzef(ib,i) = abs(pzr(ilmsb,ilms,ia,ie))
          if (il == i2lm(1,ilm)) write(*,'("density_matrix: pzef=",3i4,2f12.8)') i2dmalms(2,i), i2dmalms(1,i), ib, pzef(ib,i)
        end if
      end do
    end if
  end do
! normalize projection coefficients
  do i=1,ndim
    norm = dot_product(pzef(1:nbmax,i),pzef(1:nbmax,i))
    pzef(1:nbmax,i) = pzef(1:nbmax,i)/sqrt(norm)
  end do
! orbital rotation matrices for C3v symmetry
  u0(:) = (/ 1.d0, 0.d0, 0.d0 /)
  do ia=1,natoms
    u1(:) = (/ cos(2.d0*ia*pi/3.d0), sin(2.d0*ia*pi/3.d0), 0.d0 /)
    call orb_rotation(u0,u1,rotmat(:,:,ia,ia))
  end do
! for rotating the bonds to the z direction
  u1(:) = (/ 0.d0, 0.d0, 1.d0 /)
! 1 2
  u0(:) = (/ 1.d0, 0.d0, 0.d0 /)
  call orb_rotation(u0,u1,rotmat(:,:,1,2))
  u0(:) = -u0(:)
  call orb_rotation(u0,u1,rotmat(:,:,2,1))
! 1 3
  u0(:) = (/ cos(pi/3.d0), sin(pi/3.d0), 0.d0 /)
  call orb_rotation(u0,u1,rotmat(:,:,1,3))
  u0(:) = -u0(:)
  call orb_rotation(u0,u1,rotmat(:,:,3,1))
! 2 3
  u0(:) = (/ cos(2.d0*pi/3.d0), sin(2.d0*pi/3.d0), 0.d0 /)
  call orb_rotation(u0,u1,rotmat(:,:,2,3))
  u0(:) = -u0(:)
  call orb_rotation(u0,u1,rotmat(:,:,3,2))
! ----------------------------------------------------------------------
! **********
  do ie=1,ne
! **********
!   fill in matrix to invert
    work(:,:) = czero
    do ja=1,nasusc
      do ia=1,nasusc
!       +++++++++++++++++++++++++++++++++++++++++++++++
        if (iadmat(ia) /= 0 .and. iadmat(ja) /= 0) then
!       +++++++++++++++++++++++++++++++++++++++++++++++
!         GF at projection energies
          if (.not.lfit) then
            call projected_gf(ie,ia,ja,gf,onsite,struct)
!         GF from interpolation
          else if (ifit == 1) then
            call biratint_gf(ze(ie),ia,ja,gf,onsite,struct)
!         GF from fit
          else if (ifit == 2) then
            call ratval_gf(ia,ja,ze(ie),gf)
          end if
!         fill in GF for selected atoms and l-channels
          do jlmsb=1,nlmsba(ja)
            i3 = i2lmsb(:,jlmsb,ja)
            jb = i3(1); jlm = i3(2); js = i3(3)
!           ------------------------------------------------------------
!           do we want this element?
            do j=1,ndim
              if (i2dmalms(1,j) == lms2i(jlm,js) .and. i2dmalms(2,j) == ja) exit
            end do
            if (j > ndim) cycle
!           ------------------------------------------------------------
            do ilmsb=1,nlmsba(ia)
              i3 = i2lmsb(:,ilmsb,ia)
              ib = i3(1); ilm = i3(2); is = i3(3)
!             ----------------------------------------------------------
!             do we want this element?
              do i=1,ndim
                if (i2dmalms(1,i) == lms2i(ilm,is) .and. i2dmalms(2,i) == ia) exit
              end do
              if (i > ndim) cycle
!             ----------------------------------------------------------
!             projection on wfn at fermi energy
              work(i,j) = work(i,j) + pzef(ib,i)*gf(ilmsb,jlmsb)*pzef(jb,j)
            end do
          end do
!       ++++++
        end if
!       ++++++
      end do
    end do
!   apply symmetry transformations
!    zdmat(:,:,ie) = czero
!    do j=1,ndim
!      jlms = i2dmalms(1,j); ja = i2dmalms(2,j)
!      jlm = i2lms(1,jlms); js = i2lms(2,jlms)
!      do i=1,ndim
!        ilms = i2dmalms(1,i); ia = i2dmalms(2,i)
!        ilm = i2lms(1,ilms); is = i2lms(2,ilms)
!       ----------------------------------------------------------------
!        do q=1,ndim
!          qlms = i2dmalms(1,q); qa = i2dmalms(2,q)
!          qlm = i2lms(1,qlms); qs = i2lms(2,qlms)
!          if (qa == ja .and. qs == js) then
!            do p=1,ndim
!              plms = i2dmalms(1,p); pa = i2dmalms(2,p)
!              plm = i2lms(1,plms); ps = i2lms(2,plms)
!              if(pa == ia .and. ps == is) then
!                zdmat(i,j,ie) = zdmat(i,j,ie) + rotmat(ilm,plm,ia,ia)*work(p,q)*rotmat(jlm,qlm,ja,ja)
!              end if
!            end do
!          end if
!        end do
!       ----------------------------------------------------------------
!      end do
!    end do
!    work(:,:) = zdmat(:,:,ie)
    zdmat(:,:,ie) = work(:,:)
!   now invert G(E)^-1 = E - H(E)
    zdmat(:,:,ie) = czero
    do i=1,ndim
      zdmat(i,i,ie) = cone
    end do
    call zgesv(ndim,ndim,work,ndim,ipiv,zdmat(:,:,ie),ndim,info)
    if (info /= 0) stop 'density_matrix: failure in zgesv'
!   H(E) = E - G(E)^-1
    zdmat(:,:,ie) = -zdmat(:,:,ie)
    do i=1,ndim
      zdmat(i,i,ie) = zdmat(i,i,ie) + ze(ie)
    end do
!   rotate bonds to z axis
!    work(:,:) = czero
!    do j=1,ndim
!      jlms = i2dmalms(1,j); ja = i2dmalms(2,j)
!      jlm = i2lms(1,jlms); js = i2lms(2,jlms)
!      do i=1,ndim
!        ilms = i2dmalms(1,i); ia = i2dmalms(2,i)
!        ilm = i2lms(1,ilms); is = i2lms(2,ilms)
!       ----------------------------------------------------------------
!        do q=1,ndim
!          qlms = i2dmalms(1,q); qa = i2dmalms(2,q)
!          qlm = i2lms(1,qlms); qs = i2lms(2,qlms)
!          if (qa == ja .and. qs == js) then
!            do p=1,ndim
!              plms = i2dmalms(1,p); pa = i2dmalms(2,p)
!              plm = i2lms(1,plms); ps = i2lms(2,plms)
!              if(pa == ia .and. ps == is) then
!                work(i,j) = work(i,j) + rotmat(ilm,plm,ia,ia)*zdmat(p,q,ie)*rotmat(jlm,qlm,ja,ja)
!              end if
!            end do
!          end if
!        end do
!       ----------------------------------------------------------------
!      end do
!    end do
!    zdmat(:,:,ie) = work(:,:)
! ******
  end do
! ******
! ----------------------------------------------------------------------
! output to file
  open(file='density_matrix.dat',unit=iofile,status='replace')
  write(iofile,'("# Density matrix: first indices, then data for each energy")')
  write(iofile,'("# ", i8,"  ",100i4)') ndim, i2dmalms(1,:), i2dmalms(2,:)
! **********
  do ie=1,ne
! **********
!   analytical continuation using rational function
    if (ldosacon) then
      ze0 = real(ze(ie))
      if (ie < nedosacon/2) then
        istart = 1
        iend   = nedosacon
      else if (ie + nedosacon/2 + mod(nedosacon,2) > ne) then
        istart = ne - nedosacon + 1
        iend   = ne
      else
        istart = ie - nedosacon/2 + 1
        iend   = ie + nedosacon/2 + mod(nedosacon,2)
      end if
      do j=1,ndim
        do i=1,ndim
          call zratint(nedosacon,ze(istart:iend),zdmat(i,j,istart:iend),ze0,work(i,j),dy)
        end do
      end do
      ze0 = dosefac*(ze0 - dosezero)
      do i=1,ndim
        work(i,i) = work(i,i) - dosezero
      end do
      write(iofile,'(1000es12.4)') ze0, dosefac*work(:,:)
!   computed values for complex energy
    else
      ze0 = dosefac*(ze(ie) - dosezero)
      do i=1,ndim
        zdmat(i,i,ie) = zdmat(i,i,ie) - dosezero
      end do
      write(iofile,'(1000es12.4)') ze0, dosefac*zdmat(:,:,ie)
    end if
! ******
  end do
! ******
  close(iofile)
  deallocate(ze,zdmat,work,ipiv,gf,i2dmalms,pzef,rotmat)
! All done!
  end subroutine density_matrix
