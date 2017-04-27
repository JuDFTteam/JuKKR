  subroutine orbmoment()
! Construction of the orbital and spin angular momentum matrices
! Also matrix elements in the lmsb basis
  use global

  implicit none

  real(kind=r8b),    parameter :: invr2 = 1.d0/sqrt(2.d0)
  complex(kind=c8b), parameter :: iu = (0.d0,1.d0)
  integer(kind=i4b) :: l1, m1, lm1, l2, m2, lm2, i, j, i2(2)
  integer(kind=i4b) :: ib, ilm, is, jb, jlm, js, k, ia, i3(3)
  complex(kind=c8b) :: u(lmmax4,lmmax4), fac
  complex(kind=c8b) :: lp(lmmax4,lmmax4), lm(lmmax4,lmmax4)


  if (allocated(ldots)) deallocate(ldots)

! storage for matrix elements in lmsb basis
  allocate(ldots(nlmsb,nlmsb,nasusc))
! initialize pauli matrices
! x
  pauli(1,1,1) = ( 0.0d0, 0.0d0); pauli(1,2,1) = ( 1.0d0, 0.0d0)
  pauli(2,1,1) = ( 1.0d0, 0.0d0); pauli(2,2,1) = ( 0.0d0, 0.0d0)
! y
  pauli(1,1,2) = ( 0.0d0, 0.0d0); pauli(1,2,2) = ( 0.0d0, 1.0d0)
  pauli(2,1,2) = ( 0.0d0,-1.0d0); pauli(2,2,2) = ( 0.0d0, 0.0d0)
! z
  pauli(1,1,3) = (-1.0d0, 0.0d0); pauli(1,2,3) = ( 0.0d0, 0.0d0)
  pauli(2,1,3) = ( 0.0d0, 0.0d0); pauli(2,2,3) = ( 1.0d0, 0.0d0)
! n
  pauli(1,1,4) = ( 1.0d0, 0.0d0); pauli(1,2,4) = ( 0.0d0, 0.0d0)
  pauli(2,1,4) = ( 0.0d0, 0.0d0); pauli(2,2,4) = ( 1.0d0, 0.0d0)
! transformations from cartesian to spin labels and back
! potential cartesian to spin
  pc2s(2,1,1) = ( 1.0d0, 0.0d0); pc2s(2,1,2) = ( 0.0d0,-1.0d0); pc2s(2,1,3) = ( 0.0d0, 0.0d0); pc2s(2,1,4) = ( 0.0d0, 0.0d0) 
  pc2s(1,2,1) = ( 1.0d0, 0.0d0); pc2s(1,2,2) = ( 0.0d0, 1.0d0); pc2s(1,2,3) = ( 0.0d0, 0.0d0); pc2s(1,2,4) = ( 0.0d0, 0.0d0) 
  pc2s(2,2,1) = ( 0.0d0, 0.0d0); pc2s(2,2,2) = ( 0.0d0, 0.0d0); pc2s(2,2,3) = ( 1.0d0, 0.0d0); pc2s(2,2,4) = ( 1.0d0, 0.0d0) 
  pc2s(1,1,1) = ( 0.0d0, 0.0d0); pc2s(1,1,2) = ( 0.0d0, 0.0d0); pc2s(1,1,3) = (-1.0d0, 0.0d0); pc2s(1,1,4) = ( 1.0d0, 0.0d0) 
! potential spin to cartesian
  ps2c(1,2,1) = ( 0.5d0, 0.0d0); ps2c(2,2,1) = ( 0.0d0, 0.5d0); ps2c(3,2,1) = ( 0.0d0, 0.0d0); ps2c(4,2,1) = ( 0.0d0, 0.0d0) 
  ps2c(1,1,2) = ( 0.5d0, 0.0d0); ps2c(2,1,2) = ( 0.0d0,-0.5d0); ps2c(3,1,2) = ( 0.0d0, 0.0d0); ps2c(4,1,2) = ( 0.0d0, 0.0d0) 
  ps2c(1,2,2) = ( 0.0d0, 0.0d0); ps2c(2,2,2) = ( 0.0d0, 0.0d0); ps2c(3,2,2) = ( 0.5d0, 0.0d0); ps2c(4,2,2) = ( 0.5d0, 0.0d0) 
  ps2c(1,1,1) = ( 0.0d0, 0.0d0); ps2c(2,1,1) = ( 0.0d0, 0.0d0); ps2c(3,1,1) = (-0.5d0, 0.0d0); ps2c(4,1,1) = ( 0.5d0, 0.0d0) 
! density cartesian to spin
  dc2s(2,1,1) = ( 0.5d0, 0.0d0); dc2s(2,1,2) = ( 0.0d0, 0.5d0); dc2s(2,1,3) = ( 0.0d0, 0.0d0); dc2s(2,1,4) = ( 0.0d0, 0.0d0) 
  dc2s(1,2,1) = ( 0.5d0, 0.0d0); dc2s(1,2,2) = ( 0.0d0,-0.5d0); dc2s(1,2,3) = ( 0.0d0, 0.0d0); dc2s(1,2,4) = ( 0.0d0, 0.0d0) 
  dc2s(2,2,1) = ( 0.0d0, 0.0d0); dc2s(2,2,2) = ( 0.0d0, 0.0d0); dc2s(2,2,3) = ( 0.5d0, 0.0d0); dc2s(2,2,4) = ( 0.5d0, 0.0d0) 
  dc2s(1,1,1) = ( 0.0d0, 0.0d0); dc2s(1,1,2) = ( 0.0d0, 0.0d0); dc2s(1,1,3) = (-0.5d0, 0.0d0); dc2s(1,1,4) = ( 0.5d0, 0.0d0) 
! density spin to cartesian
  ds2c(1,2,1) = ( 1.0d0, 0.0d0); ds2c(2,2,1) = ( 0.0d0,-1.0d0); ds2c(3,2,1) = ( 0.0d0, 0.0d0); ds2c(4,2,1) = ( 0.0d0, 0.0d0) 
  ds2c(1,1,2) = ( 1.0d0, 0.0d0); ds2c(2,1,2) = ( 0.0d0, 1.0d0); ds2c(3,1,2) = ( 0.0d0, 0.0d0); ds2c(4,1,2) = ( 0.0d0, 0.0d0) 
  ds2c(1,2,2) = ( 0.0d0, 0.0d0); ds2c(2,2,2) = ( 0.0d0, 0.0d0); ds2c(3,2,2) = ( 1.0d0, 0.0d0); ds2c(4,2,2) = ( 1.0d0, 0.0d0) 
  ds2c(1,1,1) = ( 0.0d0, 0.0d0); ds2c(2,1,1) = ( 0.0d0, 0.0d0); ds2c(3,1,1) = (-1.0d0, 0.0d0); ds2c(4,1,1) = ( 1.0d0, 0.0d0) 

! orbital angular momentum
! 1 -> columns
  do lm1=1,lmmax4
    i2 = i2lm(:,lm1)
    m1 = i2(1); l1 = i2(2)
! 2 -> rows
    do lm2=1,lmmax4
      i2 = i2lm(:,lm2)
      m2 = i2(1); l2 = i2(2)
!   -----------------------------------------------------------
!   L+ in complex Ylm basis
!   row index larger than column index
      if (l2 == l1 .and. m2 == m1+1) then
        fac = l1*(l1+1) - m1*m2
        lp(lm2,lm1) = sqrt(fac)
!      write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, lp(lm2,lm1)
      else
        lp(lm2,lm1) = 0.d0
      end if
!   -----------------------------------------------------------
!   L- in complex Ylm basis
!   row index smaller than column index
      if (l2 == l1 .and. m2 == m1-1) then
        fac = l1*(l1+1) - m1*m2
        lm(lm2,lm1) = sqrt(fac)
!      write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, lm(lm2,lm1)
      else
        lm(lm2,lm1) = 0.d0
      end if
!   -----------------------------------------------------------
!   Lz in complex Ylm basis
!   column index equal to row index
      if (l2 == l1 .and. m2 == m1) then
        lorb(lm2,lm1,3) = m1
      else
        lorb(lm2,lm1,3) = 0.d0
      end if
!   -----------------------------------------------------------
!   Transformation matrix from complex to real Ylm
!   -----------------------------------------------------------
      if (l2 == l1 .and. m2 == 0 .and. m1 == 0) then
        u(lm2,lm1) = 1.d0
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 ==  m1 .and. m1 < 0) then
        u(lm2,lm1) = iu*invr2
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 == -m1 .and. m1 < 0) then
        u(lm2,lm1) = invr2
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 ==  m1 .and. m1 > 0) then
        u(lm2,lm1) = invr2*(-1)**m1
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 == -m1 .and. m1 > 0) then
        u(lm2,lm1) = -iu*invr2*(-1)**m1
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else
        u(lm2,lm1) = 0.d0
      end if
!   -----------------------------------------------------------
    end do
  end do
! -----------------------------------------------------------------
! Lx
  lorb(:,:,1) = 0.5d0*(lp + lm)
! Ly
  lorb(:,:,2) = 0.5d0*iu*(lm - lp)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (lhdio) then
    write(iodb,*) "orbmom", nlmax, lmmax4
    write(iodb,*) "Angular momentum matrices in complex Ylm basis"
    write(iodb,*) "L+"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    write(iodb,*) "L-"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lm(lm1,1:lmmax4)
    end do
    write(iodb,*) "Lx"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lorb(lm1,1:lmmax4,1)
    end do
    write(iodb,*) "Lx = i(LzLy-LyLz)"
    lp = iu*(matmul(lorb(:,:,3),lorb(:,:,2)) - matmul(lorb(:,:,2),lorb(:,:,3)))
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    write(iodb,*) "Ly"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lorb(lm1,1:lmmax4,2)
    end do
    write(iodb,*) "Ly = i(LxLz - LzLx)"
    lp = iu*(matmul(lorb(:,:,1),lorb(:,:,3)) - matmul(lorb(:,:,3),lorb(:,:,1)))
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    write(iodb,*) "Lz"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lorb(lm1,1:lmmax4,3)
    end do
    write(iodb,*) "Lz = i(LyLx - LxLy)"
    lp = iu*(matmul(lorb(:,:,2),lorb(:,:,1)) - matmul(lorb(:,:,1),lorb(:,:,2)))
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    lp = 0.d0
    do k=1,3
      lp = lp + matmul(lorb(:,:,k),lorb(:,:,k))
    end do
    write(iodb,*) "L^2"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
  end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Change to real Ylm
  do k=1,3
    lorb(:,:,k) = matmul(conjg(u),matmul(lorb(:,:,k),transpose(u)))
  end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (lhdio) then
    write(iodb,*) "Angular momentum matrices in real Ylm basis"
    write(iodb,*) "U"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') u(lm1,1:lmmax4)
    end do
    write(iodb,*) "U^dagger"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') conjg(u(1:lmmax4,lm1))
    end do
    write(iodb,*) "U^dagger U"
    lp = matmul(conjg(transpose(u)),u)
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    write(iodb,*) "Lx"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lorb(lm1,1:lmmax4,1)
    end do
    write(iodb,*) "Lx = i(LzLy-LyLz)"
    lp = iu*(matmul(lorb(:,:,3),lorb(:,:,2)) - matmul(lorb(:,:,2),lorb(:,:,3)))
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    write(iodb,*) "Ly"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lorb(lm1,1:lmmax4,2)
    end do
    write(iodb,*) "Ly = i(LxLz - LzLx)"
    lp = iu*(matmul(lorb(:,:,1),lorb(:,:,3)) - matmul(lorb(:,:,3),lorb(:,:,1)))
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    write(iodb,*) "Lz"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lorb(lm1,1:lmmax4,3)
    end do
    write(iodb,*) "Lz = i(LyLx - LxLy)"
    lp = iu*(matmul(lorb(:,:,2),lorb(:,:,1)) - matmul(lorb(:,:,1),lorb(:,:,2)))
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
    lp = 0.d0
    do k=1,3
      lp = lp + matmul(lorb(:,:,k),lorb(:,:,k))
    end do
    write(iodb,*) "L^2"
    do lm1=1,lmmax4
      write(iodb,'(100f4.1)') lp(lm1,1:lmmax4)
    end do
  end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Now put together the matrix elements of L.S
  ldots = 0.d0
  do ia=1,nasusc
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
!        if (ib == jb) then
        do k=1,3
          ldots(i,j,ia) = ldots(i,j,ia) + lorb(ilm,jlm,k)*pauli(is,js,k)
        end do
!        end if
      end do
    end do
  end do
! All done
  end subroutine orbmoment
