  subroutine kxc_sumrule()
! assembles the xc ALDA kernel in the density basis from the sum rule
! only transverse kernel at present
! diagonal in the density basis
  use global

  implicit none

! this is the magnetization calculated in ms_from_bxc3
!  complex(kind=c8b) :: mssusc(nbmax,nbmax,lmmax,lmmax,nasusc2)
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-7
  integer(kind=i4b) :: i3(3), i4(4), ia, ia2, ja, ja2, iq, jq, is, js, ilm, jlm, ib, jb, q1, ib1, ilm1, is1, q2, ib2, ilm2, is2, nden0
  complex(kind=c8b) :: bxcden(ndensum), mxcden(ndensum), mtotden(ndensum), kxcden(ndensum), mxclm(lmmax0,nasusc2), mtotlm(lmmax0,nasusc2)
  real(kind=r8b)    :: re, im
  integer(kind=i4b) :: iqmax, jqmax, ipiv(ndensum), info, iqp(ndensum), iqm(ndensum), iq2ilmb(3,ndensum)
  logical           :: exists

! if the kernel file exists just read in values
  inquire(file='kernel.dat',exist=exists)
  if (exists) then
    write(*,'("kernel read from file")')
    open(file='kernel.dat',unit=iofile,status='old')
    kernel = 0.d0
    read(iofile,*) iqmax
    do q1=1,iqmax
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
  iqmax = 0; jqmax = 0
  do iq=1,ndensum
    i4 = i2almsbden(:,iq)
    ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
    ia = iasusc2(ia2)
    nden0 = sum(nalmsbden(1:ia2-1))
!   m+ only:  m+ = \chi0^+- B-
    if (is == 1) then
      iqmax = iqmax + 1
      iqp(iqmax) = iq
      do q2=1,nlmsba(ia)
        i3 = i2lmsb(:,q2,ia)
        ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
        do q1=1,nlmsba(ia)
          i3 = i2lmsb(:,q1,ia)
          ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
          jq = q1 + (q2-1)*nlmsba(ia)
!          if (is1 == 2 .and. is2 == 1) then
!            write(*,'(5i4,4es16.8)') iq, ia, ia2, q1, q2, dengaunt(jq,iq,ia2), mssusc(ib1,ib2,ilm1,ilm2,ia2)
          mxcden(iqmax)  = mxcden(iqmax)  + dengaunt(jq,iq-nden0,ia2)*mxcsusc(ib1,ib2,ilm1,ilm2,ia2)
          mtotden(iqmax) = mtotden(iqmax) + dengaunt(jq,iq-nden0,ia2)*mtotsusc(ib1,ib2,ilm1,ilm2,ia2)
!          end if
        end do
      end do
      re = real(mxcden(iqmax)); im = aimag(mxcden(iqmax))
      if (abs(re) < tol) re = 0.d0
      if (abs(im) < tol) im = 0.d0
      mxcden(iqmax) = cmplx(re,im)
!      write(*,'(4i4,2es16.8)') ia2, i2lm(:,ilm), ib, mxcden(iqmax)
      mxclm(ilm,ia2) = mxclm(ilm,ia2) + suscnorm(iq)*mxcden(iqmax)
      re = real(mtotden(iqmax)); im = aimag(mtotden(iqmax))
      if (abs(re) < tol) re = 0.d0
      if (abs(im) < tol) im = 0.d0
      mtotden(iqmax) = cmplx(re,im)
!      write(*,'(4i4,2es16.8)') ia2, i2lm(:,ilm), ib, mxcden(iqmax)
      mtotlm(ilm,ia2) = mtotlm(ilm,ia2) + suscnorm(iq)*mtotden(iqmax)
    else if (is == 2) then
      jqmax = jqmax + 1
      iqm(jqmax) = iq
    end if
  end do
  write(*,'("mssusc in density basis: tot and xc")')
  do ia2=1,nasusc2
    do ilm=1,lmmax0
      write(*,'(3i4,4es16.8)') ia2, i2lm(:,ilm), mtotlm(ilm,ia2), mxclm(ilm,ia2)
    end do
  end do
  write(*,'("kxc_sumrule: iqmax,jqmax=",2i8)') iqmax, jqmax
!  do q1=1,iqmax
!    write(*,'(8i4)') i2almsbden(1,iqp(q1)), i2almsbden(1,iqm(q1)), i2almsbden(2,iqp(q1)), i2almsbden(2,iqm(q1)),  &
!                     i2almsbden(3,iqp(q1)), i2almsbden(3,iqm(q1)), i2almsbden(4,iqp(q1)), i2almsbden(4,iqm(q1))
!  end do
! ----------------------------------------------------------------------
! symmetrize the static susceptibility
  kernel = 0.d0
  do q2=1,jqmax
    do q1=1,iqmax
      kernel(q1,q2) = 0.5d0*(kssusc0(iqp(q1),iqm(q2)) + kssusc0(iqm(q1),iqp(q2)))
      kernel(q1,q2) = kernel(q1,q2) + 0.5d0*(kssusc0(iqp(q1),iqp(q2)) + kssusc0(iqm(q1),iqm(q2)))
    end do
  end do
  do iq=1,iqmax
    write(*,'(10000es10.2)') kernel(iq,1:jqmax)
  end do
! NEW: in local frame m_z = (twisted \chi_0)_zz B_z
  bxcden(1:iqmax) = mxcden(1:iqmax)
  call zgesv(iqmax,1,kernel,ndensum,ipiv,bxcden,ndensum,info)
  if (info /= 0) stop 'kxc_sumrule: failure in zgesv'
! divide bxc by the magnetization
  mxclm = 0.d0
  do q1=1,iqmax
    i4 = i2almsbden(:,iqp(q1))
    ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
    ia = iasusc2(ia2)
    re = real(bxcden(q1)); im = aimag(bxcden(q1))
    if (abs(re) < tol) re = 0.d0
    if (abs(im) < tol) im = 0.d0
    bxcden(q1) = cmplx(re,im)
    kxcden(q1) = 0.d0
    if (abs(mtotden(q1)) > tol .and. abs(bxcden(q1)) > tol) kxcden(q1) = bxcden(q1)/mtotden(q1)
    write(444,'(4i4,8es16.8)') ia2, is, ilm, ib, mtotden(q1), mxcden(q1), bxcden(q1), kxcden(q1)
    mxclm(ilm,ia2) = mxclm(ilm,ia2) + suscnorm(iqp(q1))*bxcden(q1)
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
  write(iofile,'(i8)') iqmax
  do q1=1,iqmax
    kernel(iqp(q1),iqm(q1)) = kxcden(q1)
    kernel(iqm(q1),iqp(q1)) = kxcden(q1)
    write(iofile,'(2i8,2es16.8)') iqp(q1), iqm(q1), kxcden(q1)
  end do
  close(iofile)
! All done
  end subroutine kxc_sumrule
