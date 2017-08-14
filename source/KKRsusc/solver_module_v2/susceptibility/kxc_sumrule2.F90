  subroutine kxc_sumrule2()
! assembles the xc ALDA kernel in the density basis from the sum rule
! only transverse kernel at present
! diagonal in the density basis
  use global

  implicit none

! this is the magnetization calculated in ms_from_bxc3
!  complex(kind=c8b) :: mssusc(nbmax,nbmax,lmmax,lmmax,nasusc2)
! -----------------------------------------------------------------
  real(kind=r8b),    parameter   :: tol = 1.d-7
  complex(kind=c8b), allocatable :: bxcden(:), mxcden(:), mtotden(:), kxcden(:), mxclm(:,:), mtotlm(:,:), twist(:,:,:,:)
  integer(kind=i4b) :: i, i3(3), i4(4)
  integer(kind=i4b) :: ia, ia2, iq, is, ilm, ib, q1, ib1, ilm1, is1
  integer(kind=i4b) :: ja, ja2, jq, js, jlm, jb, q2, ib2, ilm2, is2
  integer(kind=i4b) :: ipiv(ndensum), info, nden0, iqs(nsmax2), iq2qs(ndensum), nqmax
  real(kind=r8b)    :: re, im
  logical           :: file_exists


! ----------------------------------------------------------------------
! if the kernel file exists just read in values
  inquire(file='kernel.dat',exist=file_exists)
  if (file_exists) then
    write(*,'("kernel read from file")')
    open(file='kernel.dat',unit=iofile,status='old')
    kernel = 0.d0
    read(iofile,*) nqmax
    do i=1,nqmax
      read(iofile,*) iq, jq, re, im
!     +-
      kernel(iq,jq) = cmplx(re,im)
!     -+
      kernel(jq,iq) = cmplx(re,im)
    end do
    close(iofile)
    return
  end if
! ----------------------------------------------------------------------
! construct pointers to spin-independent basis index
  iqs(:) = 0
  do iq=1,ndensum
    is = i2almsbden(3,iq)
    iqs(is) = iqs(is) + 1
    iq2qs(iq) = iqs(is)
  end do
  write(*,'("kxc_sumrule2: iqs(1:nsmax2)=",4i8)') iqs(:)
  if (any(iqs(:) - iqs(1) /= 0)) stop 'kxcsumrule2: inconsistency in iqs!'
  nqmax = iqs(1)
  allocate(bxcden(nqmax),mxcden(nqmax),mtotden(nqmax),kxcden(nqmax))
  allocate(mxclm(lmmax0,nasusc2),mtotlm(lmmax0,nasusc2),twist(nsmax2,nsmax2,nasusc2,nasusc2))
! ----------------------------------------------------------------------
! transform the magnetization from the GF to the density basis
!  do ilm1=1,lmmax
!    do ib1=1,nbmax
!      write(*,'("mssusc for ilm1,ib1=",2i4,2es16.8)') ilm1, ib1, mssusc(ib1,ib1,ilm1,ilm1,1)
!    end do
!  end do
!  do iq=1,ndensum
!    write(*,'(10000es10.2)') kssusc0(iq,:)
!  end do
  mxcden = 0.d0; mxclm = 0.d0; mtotden = 0.d0; mtotlm = 0.d0
  do iq=1,ndensum
    i4 = i2almsbden(:,iq)
    ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
    ia = iasusc2(ia2)
    nden0 = sum(nalmsbden(1:ia2-1))
!   construct the magnetization (only once)
    if (is == 1) then
      i = iq2qs(iq)
      do q2=1,nlmsba(ia)
        i3 = i2lmsb(:,q2,ia)
        ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
        do q1=1,nlmsba(ia)
          i3 = i2lmsb(:,q1,ia)
          ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
          jq = q1 + (q2-1)*nlmsba(ia)
          mxcden(i)  = mxcden(i)  + dengaunt(jq,iq-nden0,ia2)*mxcsusc(ib1,ib2,ilm1,ilm2,ia2)
          mtotden(i) = mtotden(i) + dengaunt(jq,iq-nden0,ia2)*mtotsusc(ib1,ib2,ilm1,ilm2,ia2)
        end do
      end do
      call zchop(mxcden(i),tol)
!      write(*,'(4i4,2es16.8)') ia2, i2lm(:,ilm), ib, mxcden(i)
      mxclm(ilm,ia2) = mxclm(ilm,ia2) + suscnorm(iq)*mxcden(i)
      call zchop(mtotden(i),tol)
!      write(*,'(4i4,2es16.8)') ia2, i2lm(:,ilm), ib, mxcden(i)
      mtotlm(ilm,ia2) = mtotlm(ilm,ia2) + suscnorm(iq)*mtotden(i)
    end if
  end do
  write(*,'("mssusc in density basis: tot and xc")')
  do ia2=1,nasusc2
    do ilm=1,lmmax0
      write(*,'(3i4,4es16.8)') ia2, i2lm(:,ilm), mtotlm(ilm,ia2), mxclm(ilm,ia2)
    end do
  end do
! ----------------------------------------------------------------------
! twist the static susceptibility
  call build_twist(twist)
  kernel = 0.d0
  do jq=1,ndensum
    js = i2almsbden(3,jq); ja2 = i2almsbden(4,jq)
    q2 = iq2qs(jq)
    do iq=1,ndensum
      is = i2almsbden(3,iq); ia2 = i2almsbden(4,iq)
      q1 = iq2qs(iq)
!      write(*,'("iq,jq,ia2,ja2,is,js,q1,q2=",8i4)') iq, jq, ia2, ja2, is, js, q1, q2
      kernel(q1,q2) = kernel(q1,q2) + 0.5d0*twist(is,js,ia2,ja2)*kssusc0(iq,jq)
    end do
  end do
  do q1=1,nqmax
    write(*,'(10000es10.2)') kernel(q1,1:nqmax)
  end do
! ----------------------------------------------------------------------
! Local frame: m_z = (twisted \chi_0)_zz B_z
  bxcden(:) = mxcden(:)
  call zgesv(nqmax,1,kernel,ndensum,ipiv,bxcden,nqmax,info)
  if (info /= 0) stop 'kxc_sumrule: failure in zgesv'
! divide bxc by the magnetization
  mxclm = 0.d0
  do iq=1,ndensum
    i4 = i2almsbden(:,iq)
    ilm = i4(2); is = i4(3); ia2 = i4(4)
    ia = iasusc2(ia2)
    if (is == 1) then
      q1 = iq2qs(iq)
      call zchop(bxcden(q1),tol)
      kxcden(q1) = 0.d0
      if (abs(mtotden(q1)) > tol .and. abs(bxcden(q1)) > tol) kxcden(q1) = bxcden(q1)/mtotden(q1)
!      write(*,'(4i4,8es16.8)') ia2, is, ilm, ib, mtotden(q1), mxcden(q1), bxcden(q1), kxcden(q1)
      mxclm(ilm,ia2) = mxclm(ilm,ia2) + suscnorm(iq)*bxcden(q1)
    end if
  end do
  write(*,'("bxc in density basis:")')
  do ia2=1,nasusc2
    do ilm=1,lmmax0
      write(*,'(3i4,2es16.8)') ia2, i2lm(:,ilm), mxclm(ilm,ia2)
    end do
  end do
! ----------------------------------------------------------------------
! copy results to the kernel and output
  write(*,'("kernel written to file")')
  open(file='kernel.dat',unit=iofile,status='new')
  kernel = 0.d0
  write(iofile,'(i8)') nqmax
  do jq=1,ndensum
    js = i2almsbden(3,jq); q2 = iq2qs(jq)
    do iq=1,ndensum
      is = i2almsbden(3,iq); q1 = iq2qs(iq)
!     kernel assumed diagonal in spatial basis
      if (q1 == q2) then
!       +-
        if (is == 1 .and. js == 2) then
          kernel(iq,jq) = kxcden(q1)
          write(iofile,'(2i8,2es16.8)') iq, jq, kxcden(q1)
        end if
!       -+
        if (is == 2 .and. js == 1) kernel(iq,jq) = kxcden(q1)
      end if
    end do
  end do
  close(iofile)
! ----------------------------------------------------------------------
  deallocate(bxcden,mxcden,mtotden,kxcden)
  deallocate(mxclm,mtotlm,twist)
! All done
  end subroutine kxc_sumrule2
