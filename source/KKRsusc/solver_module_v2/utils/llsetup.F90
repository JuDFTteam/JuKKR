  subroutine llsetup(lmaxll)
! Setup the LL mesh for integration on a sphere

  use global, only: i4b, r8b, c8b, nll, ull, wll

  implicit none

! maximum degree of Ylm that will integrate to zero
  integer(kind=i4b), intent(in)  :: lmaxll
!----------------------------------------------------------------------
  real(kind=r8b), allocatable :: xll(:), yll(:), zll(:)
  integer(kind=i4b) :: ill

  select case (lmaxll)
    case (0)
      nll = 2
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      xll(1) = 0.d0; yll(1) = 0.d0; zll(1) =  1.d0; wll(1) = 0.5d0
      xll(2) = 0.d0; yll(2) = 0.d0; zll(2) = -1.d0; wll(2) = 0.5d0
    case (1)
      nll = 6
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0006(xll,yll,zll,wll,nll)
    case (2)
      nll = 14
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0014(xll,yll,zll,wll,nll)
    case (3)
      nll = 26
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0026(xll,yll,zll,wll,nll)
    case (4)
      nll = 38
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0038(xll,yll,zll,wll,nll)
    case (5)
      nll = 50
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0050(xll,yll,zll,wll,nll)
    case (6)
      nll = 74
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0074(xll,yll,zll,wll,nll)
    case (7)
      nll = 86
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0086(xll,yll,zll,wll,nll)
    case (8)
      nll = 110
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0110(xll,yll,zll,wll,nll)
    case (9)
      nll = 146
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0146(xll,yll,zll,wll,nll)
    case (10)
      nll = 170
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0170(xll,yll,zll,wll,nll)
    case (11)
      nll = 194
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0194(xll,yll,zll,wll,nll)
    case (12)
      nll = 230
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0230(xll,yll,zll,wll,nll)
    case (13)
      nll = 266
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0266(xll,yll,zll,wll,nll)
    case (14)
      nll = 302
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0302(xll,yll,zll,wll,nll)
    case (15)
      nll = 350
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0350(xll,yll,zll,wll,nll)
    case (16)
      nll = 434
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0434(xll,yll,zll,wll,nll)
    case (17)
      nll = 590
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0590(xll,yll,zll,wll,nll)
    case (18)
      nll = 770
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0770(xll,yll,zll,wll,nll)
    case (19)
      nll = 974
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD0974(xll,yll,zll,wll,nll)
    case (20)
      nll = 1202
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD1202(xll,yll,zll,wll,nll)
    case (21)
      nll = 1454
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD1454(xll,yll,zll,wll,nll)
    case (22)
      nll = 1730
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD1730(xll,yll,zll,wll,nll)
    case (23)
      nll = 2030
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD2030(xll,yll,zll,wll,nll)
    case (24)
      nll = 2354
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD2354(xll,yll,zll,wll,nll)
    case (25)
      nll = 2702
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD2702(xll,yll,zll,wll,nll)
    case (26)
      nll = 3074
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD3074(xll,yll,zll,wll,nll)
    case (27)
      nll = 3470
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD3470(xll,yll,zll,wll,nll)
    case (28)
      nll = 3890
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD3890(xll,yll,zll,wll,nll)
    case (29)
      nll = 4334
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD4334(xll,yll,zll,wll,nll)
    case (30)
      nll = 4802
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD4802(xll,yll,zll,wll,nll)
    case (31)
      nll = 5294
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD5294(xll,yll,zll,wll,nll)
    case (32)
      nll = 5810
      allocate(xll(nll),yll(nll),zll(nll),ull(3,nll),wll(nll))
      call LD5810(xll,yll,zll,wll,nll)
    case default
      stop 'llsetup: unknown lmaxll'
  end select
  do ill=1,nll
    ull(:,ill) = (/ xll(ill), yll(ill), zll(ill) /)
  end do
! All done
  end subroutine llsetup 
