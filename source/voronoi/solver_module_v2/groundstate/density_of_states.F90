  subroutine density_of_states(ia,onsite,struct)
  use global

  implicit none

! Atom of interest
  integer(kind=i4b), intent(in) :: ia
! Onsite GF
  logical,           intent(in) :: onsite
! Structural GF
  logical,           intent(in) :: struct
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: invpi = 0.25d0/atan(1.d0)
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi), cone = (1.d0,0.d0), czero = (0.d0,0.d0)
  complex(kind=c8b), allocatable :: zdos(:,:,:,:,:), zdosl(:,:,:), zdoslm(:,:,:)
  complex(kind=c8b), allocatable :: ze(:), gf(:,:), zdosal(:,:), zdosalm(:,:)
  real(kind=r8b),    allocatable :: center(:,:), norm(:,:), centerl(:,:)
  complex(kind=c8b) :: za, zb, ze0, dy, zms(1:4), zmo(1:3), zmso(1:3,1:3)
  real(kind=r8b)    :: norm2
  character*14      :: filename
  integer(kind=i4b) :: ie, ne, istart, iend
  integer(kind=i4b) :: i2(2), i3(3)
  integer(kind=i4b) :: j, jlms, js, jlm, jl, jm, jb
  integer(kind=i4b) :: i, ilms, is, ilm, il, im, ib
  integer(kind=i4b) :: info, lwork, ieval(nlms)
  complex(kind=c8b), allocatable :: work(:), rhomat2(:,:)
  real(kind=r8b),    allocatable :: rwork(:), evals(:)
  complex(kind=c8b) :: znorm

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
! ----------------------------------------------------------------------
! Storage
  allocate(zdos(lmmax,lmmax,nsmax,nsmax,ne),zdosl(0:nlmax,nsmax,ne),zdoslm(lmmax,nsmax,ne),gf(nlmsb,nlmsb))
  allocate(center(lmmax,nsmax),norm(lmmax,nsmax),centerl(0:nlmax,nsmax))
  allocate(zdosal(0:nlmax,nsmax),zdosalm(lmmax,nsmax))
  zdosl = czero; zdoslm = czero; zdos = czero
  if (lrhodiag) then
    lwork = 2*nlms
    allocate(work(lwork),rwork(2*lwork),evals(nlms),rhomat2(nlms,nlms))
  end if
! ----------------------------------------------------------------------
! **********
  do ie=1,ne
! **********
!   GF at projection energies
    if (.not.lfit) then
      call projected_gf(ie,ia,ia,gf,onsite,struct)
!   GF from interpolation
    else if (ifit == 1) then
!      call baryint_gf2(nesusc,esusc,ia,ia,ze(ie),numd,gf,onsite,struct)
!      call ratint_gf2(nesusc,esusc,ia,ia,ze(ie),numd,gf,onsite,struct)
      call biratint_gf(ze(ie),ia,ia,gf,onsite,struct)
!   GF from fit
    else if (ifit == 2) then
      call ratval_gf(ia,ia,ze(ie),gf)
    end if
!   -----------------------------------------------
!   lm-DOS calculation
    if (lrhodiag) then
      rhomat2 = czero
      do j=1,nlmsba(ia)
        i3 = i2lmsb(:,j,ia)
        jb = i3(1); jlm = i3(2); js = i3(3)
        jlms = lms2i(jlm,js)
        do i=1,nlmsba(ia)
          i3 = i2lmsb(:,i,ia)
          ib = i3(1); ilm = i3(2); is = i3(3)
          ilms = lms2i(ilm,is)
          rhomat2(ilms,jlms) = rhomat2(ilms,jlms) + (conjg(gf(j,i)) - gf(i,j))/i2pi
        end do
      end do
      call zheev('V','U',nlms,rhomat2,nlms,evals,work,lwork,rwork,info)
      if (info /= 0) stop 'integrated_charge: failure in zheev'
!     sort the eigenvalues according to ilms
      do ilms=1,nlms
        norm2 = 0.d0
!       search for largest element for ilms
        do jlms=1,nlms
          if (abs(rhomat2(ilms,jlms)) > norm2) then
!           if this eigenvector was already picked skip it
            if (any(ieval(1:ilms-1) == jlms)) cycle
            norm2 = abs(rhomat2(ilms,jlms))
            ieval(ilms) = jlms
          end if
        end do
!       normalize the eigenvalues according to ilms
!        znorm = rhomat2(ilms,ieval(ilms))
!        znorm = conjg(znorm)/sqrt(abs(znorm))
!        rhomat2(1:nlms,ieval(ilms)) = rhomat2(1:nlms,ieval(ilms))*znorm
        i2 = i2lms(:,ilms)
        ilm = i2(1); is = i2(2)
        zdos(ilm,ilm,is,is,ie) = evals(ieval(ilms))
      end do
    else
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
          zdos(ilm,jlm,is,js,ie) = zdos(ilm,jlm,is,js,ie) + gf(i,j)*overlap(i,j,ia)
!          if (ilm == jlm .and. is == js) zdosl(il,is,ie) = zdosl(il,is,ie) + gf(i,j)*overlap(i,j,ia)
        end do
      end do
    end if
!   l-DOS and lm-DOS
    do ilm=1,lmmax
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      do js=1,nsmax
        do is=1,nsmax
!         spin-down projection
          zdoslm(ilm,1,ie)     = zdoslm(ilm,1,ie)     + 0.5d0*zdos(ilm,ilm,is,js,ie)*(pauli(js,is,4) - sum(pauli(js,is,1:3)*magdir(:,ia)))
!         spin up projection
          zdoslm(ilm,nsmax,ie) = zdoslm(ilm,nsmax,ie) + 0.5d0*zdos(ilm,ilm,is,js,ie)*(pauli(js,is,4) + sum(pauli(js,is,1:3)*magdir(:,ia)))
        end do
      end do
      zdosl(il,1,ie)     = zdosl(il,1,ie)     + zdoslm(ilm,1,ie)
      zdosl(il,nsmax,ie) = zdosl(il,nsmax,ie) + zdoslm(ilm,nsmax,ie)
    end do
! ******
  end do
! ******
  if (.not.lrhodiag) zdos = -invpi*zdos; zdosl = -invpi*zdosl; zdoslm = -invpi*zdoslm
! ----------------------------------------------------------------------
! center of mass (only works for straight line DOS)
  center = 0.d0; norm = 0.d0; centerl = 0.d0
  do ie=1,ne
    do is=1,nsmax
      do ilm=1,lmmax
        center(ilm,is) = center(ilm,is) + real(ze(ie))*aimag(zdos(ilm,ilm,is,is,ie))
        norm(ilm,is)   = norm(ilm,is)   + aimag(zdos(ilm,ilm,is,is,ie))
      end do
    end do
  end do
  center = center/(norm + 1.d-12)
  do is=1,nsmax
    do ilm=1,lmmax
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      centerl(il,is) = centerl(il,is) + center(ilm,is)/(2.d0*il + 1.d0)
    end do
  end do
! ----------------------------------------------------------------------
! l-DOS to file
  write(filename,'("ldos",i4.4,".dat")') ia
  open(file=filename,unit=iofile,status='replace')
  write(iofile,'("# Centers of mass")')
  write(iofile,'("#",100f8.4)') ((dosefac*(centerl(il,is)-dosezero),il=0,nlmax),is=1,nsmax)
  write(iofile,'("# ze  zdosl(0:nlmax,1:nsmax)")')
  do ie=1,ne
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
      do is=1,nsmax
        do il=0,nlmax
          call zratint(nedosacon,ze(istart:iend),zdosl(il,is,istart:iend),ze0,zdosal(il,is),dy)
        end do
      end do
      ze0 = dosefac*(ze0 - dosezero)
      write(iofile,'(100es16.8)') ze0, sum(zdosal(:,1))/dosefac, sum(zdosal(:,nsmax))/dosefac, zdosal(:,:)/dosefac
!   computed values for complex energy
    else
      ze0 = dosefac*(ze(ie) - dosezero)
      write(iofile,'(100es16.8)') ze0, sum(zdosl(:,1,ie))/dosefac, sum(zdosl(:,nsmax,ie))/dosefac, zdosl(:,:,ie)/dosefac
    end if
  end do
  close(iofile)
! ----------------------------------------------------------------------
! lm-DOS to file
  write(filename,'("lmdos",i4.4,".dat")') ia
  open(file=filename,unit=iofile,status='replace')
  write(iofile,'("# Centers of mass")')
  write(iofile,'("#",100f8.4)') ((dosefac*(center(ilm,is)-dosezero),ilm=1,lmmax),is=1,nsmax)
  write(iofile,'("# ze  zdos(1:lmmax,1:nsmax)")')
  do ie=1,ne
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
      do is=1,nsmax
        do ilm=1,lmmax
!          call zratint(nedosacon,ze(istart:iend),zdos(ilm,ilm,is,is,istart:iend),ze0,zdosalm(ilm,is),dy)
          call zratint(nedosacon,ze(istart:iend),zdoslm(ilm,is,istart:iend),ze0,zdosalm(ilm,is),dy)
        end do
      end do
      ze0 = dosefac*(ze0 - dosezero)
      write(iofile,'(100es16.8)') ze0, zdosalm(:,:)/dosefac
!   computed values for complex energy
    else
      ze0 = dosefac*(ze(ie) - dosezero)
      write(iofile,'(100es16.8)') ze0, ((zdoslm(ilm,is,ie)/dosefac,ilm=1,lmmax),is=1,nsmax)
    end if
  end do
  close(iofile)
! ----------------------------------------------------------------------
! observables to file
  write(filename,'("lobs",i4.4,".dat")') ia
  open(file=filename,unit=iofile,status='replace')
  write(iofile,'("# Q, ms_x, ms_y, ms_z, mo_x, mo_y, mo_z, mso_ab (a,b=x,y,z)")')
  do ie=1,ne
!   orbital moment
    zmo(:) = (0.d0,0.d0)
    do is=1,nsmax
      do jlm=1,lmmax
      do ilm=1,lmmax
        zmo(1:3) = zmo(1:3) + zdos(ilm,jlm,is,is,ie)*lorb(jlm,ilm,1:3)
      end do
      end do
    end do
!   spin moment
    zms(:) = (0.d0,0.d0)
    do js=1,nsmax
    do is=1,nsmax
      do jlm=1,lmmax
        zms(1:4) = zms(1:4) + zdos(jlm,jlm,is,js,ie)*pauli(js,is,1:4)
      end do
    end do
    end do
!   spin-orbital moment
    zmso(:,:) = (0.d0,0.d0)
    do i=1,3
      do js=1,nsmax
      do is=1,nsmax
        do jlm=1,lmmax
        do ilm=1,lmmax
          zmso(1:3,i) = zmso(1:3,i) + zdos(ilm,jlm,is,js,ie)*lorb(jlm,ilm,1:3)*pauli(js,is,i)
        end do
        end do
      end do
      end do
    end do
    ze0 = dosefac*(ze(ie) - dosezero)
    write(iofile,'(100es16.8)') ze0, zms(4)/dosefac, zms(1:3)/dosefac, zmo(1:3)/dosefac, zmso(1:3,1:3)/dosefac
  end do
  close(iofile)
! ----------------------------------------------------------------------
! density matrix for d-block
!  if (ldosdmat .and. iadmat(ia) == 1) then
!    write(filename,'("dmate",i4.4,".dat")') ia
!    open(file=filename,unit=iofile,status='replace')
!    write(iofile,'("# spin down: (-2,-2) (-2,+1) (+1,-2) (+1,+1) (0,0) (-1,-1) (-1,+2) (+2,-1) (+2,+2)")')
!    do ie=1,ne
!      ze0 = dosefac*(ze(ie) - dosezero)
!      write(iofile,'(100es16.8)') ze0, zdos(5,5,1,1,ie)/dosefac, zdos(5,8,1,1,ie)/dosefac, zdos(8,5,1,1,ie)/dosefac, zdos(8,8,1,1,ie)/dosefac, &
!          zdos(7,7,1,1,ie)/dosefac, zdos(6,6,1,1,ie)/dosefac, zdos(6,9,1,1,ie)/dosefac, zdos(9,6,1,1,ie)/dosefac, zdos(9,9,1,1,ie)/dosefac
!    end do
!    close(iofile)
!  end if
! ----------------------------------------------------------------------
  deallocate(ze,zdos,zdosl,zdoslm,zdosal,zdosalm,gf,center,norm,centerl)
  if (lrhodiag) deallocate(work,rwork,evals,rhomat2)
! All done!
  end subroutine density_of_states
