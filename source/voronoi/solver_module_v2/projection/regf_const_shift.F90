  subroutine regf_const_shift()
! shifts the real part of the GF by a constant

  implicit none

  integer(kind=i4b) :: ne, i, ie, ia, lmax
  integer(kind=i4b) :: i3(3), i2(2), ib, ilm, is, im, il
  real(kind=r8b), allocatable :: dummy(:,:,:)
  complex(kind=c8b) :: ze, gf(nlmsb,nlmsb), zdos(0:nlmax,nsmax), avgdev(0:nlmax,nsmax)
  real(kind=r8b)    :: zere, zeim

! read file with complex dos from impurity program
  open(file='imp.zdos',unit=iofile,status='old')
  read(iofile,*) ne, lmax
  allocate(dummy(2,-1:lmax,nsmax))
  do ia=1,nasusc
    avgdev = 0.d0
    call ratfit_gf(ia,ia,.true.,.true.)
    do ie=1,ne
      zdos = 0.d0
!     total zdos put into il = -1 component
      read(iofile,*) i, zere, zeim, dummy
      do is=1,2
        do il=0,nlmax
          zdos(il,is) = cmplx(dummy(1,il,is),dummy(2,il,is))
        end do
      end do
!      write(*,'(i4,100f8.3)') i, zere, zeim, zdos
      ze = cmplx(zere,zeim)
      call ratval_gf(ia,ia,ze,gf)
!     **** zdos from projected GF ****
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
        zdos(il,is) = zdos(il,is) - gf(i,i)*overlap(i,i,ia)
      end do
!     ********************************
!     running averaged difference
      avgdev = avgdev + zdos(0:nlmax,:)
    end do
    avgdev = avgdev/ne
!   done with the differences; print them
    write(*,'("ia, il, avgdev(is=1:2)")')
    do il=0,nlmax
      write(*,'(2i4,4f16.8)') ia, il, avgdev(il,:)
    end do
!   now implement correction
    do ie=1,nesusc
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
!     the location of the shift is a bit arbitrary
        if (ib == 1) then
!        if (ib <= iwsusc(il,is,ia)) then
!          gfpq(i,i,ia,ie) = gfpq(i,i,ia,ie) + fudge*real(avgdev(il,is))/((2*il+1)*iwsusc(il,is,ia))
          gfpq(i,i,ia,ie) = gfpq(i,i,ia,ie) + fudge*real(avgdev(il,is))/((2*il+1))
        end if
      end do
    end do
  end do
  close(iofile)
  deallocate(dummy)
! All done!
  end subroutine regf_const_shift

