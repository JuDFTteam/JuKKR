  subroutine overlaps_gf()
! Compute and store overlaps between basis functions
! Also constructs the pointers for the product basis
  use global

  implicit none

  integer(kind=i4b) :: ia, i3(3), nr, ia2
  integer(kind=i4b) :: i, j, is, js, ilm, jlm, il, jl, im, jm, ib, jb, q1, q2
  real(kind=r8b)    :: ram
  complex(kind=c8b) :: work(nrmax)
  complex(kind=c8b), external :: radint

! Overlaps between products of basis functions
  overlap = 0.d0
  do ia=1,nasusc
    nr = nrpts(ia)
!    write(iodb,*) "Overlaps for ia=", ia
!    write(iodb,'("  ib  im  il  is  jb  jm  jl  js  overlap")')
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        work = phiref(:,jb,i2lm(2,jlm),js,ia)*phiref(:,ib,i2lm(2,ilm),is,ia)
        overlap(i,j,ia) = radint(nr,work(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
!        write(iodb,'(8i4,es16.6)') ib, i2lm(:,ilm), is, jb, i2lm(:,jlm), js, overlap(i,j,ia)
      end do
    end do
  end do
! ----------------------------------------------------------------------
  if (allocated(nalmsbgf)) deallocate(nalmsbgf)
  allocate(nalmsbgf(nasusc))
  nalmsbgf = nlmsba*nlmsba
  ngfsum = sum(nalmsbgf); ngfmax = maxval(nlmsba*nlmsba)
  write(*,'("overlaps_gf:    ngfsum =",i8,"; ngfmax =",i8)') ngfsum, ngfmax
! combined indices (decoding only)
  if (allocated(i2almsbgf)) deallocate(i2almsbgf)
  if (allocated(almsbgf2i)) deallocate(almsbgf2i) !added by Sascha
  allocate(i2almsbgf(3,ngfsum))
  allocate(almsbgf2i(ngfmax,ngfmax,nasusc))
  i = 0
  do ia=1,nasusc
    do q2=1,nlmsba(ia)
      do q1=1,nlmsba(ia)
        i = i + 1
        almsbgf2i(q1,q2,ia) = i
        i2almsbgf(:,i) = (/q1,q2,ia/)
      end do
    end do
  end do
  write(*,'("overlaps_gf:    ngfsum  should be ",i8)') i
! Added by Sascha:
! Gradient basis
  ngradmax = 0
  ngradsum = 0
  if(allocated(nalmsbgrad)) deallocate(nalmsbgrad)
  allocate(nalmsbgrad(nasusc2))
  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    ngradsum = ngradsum + nalmsbgf(ia)
    nalmsbgrad(ia2) = nalmsbgf(ia)
    if(nalmsbgf(ia) > ngradmax) ngradmax = nalmsbgf(ia)
  end do
  if (allocated(i2almsbgrad)) deallocate(i2almsbgrad)
  if (allocated(almsbgrad2i)) deallocate(almsbgrad2i)
  allocate(i2almsbgrad(3,ngradsum))
  allocate(almsbgrad2i(ngradmax,ngradmax,nasusc2))
  i = 0
  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    do q2=1,nlmsba(ia)
      do q1=1,nlmsba(ia)
        i = i + 1
        almsbgrad2i(q1,q2,ia2) = i
        i2almsbgrad(:,i) = (/q1,q2,ia2/)
      end do
    end do
  end do
! All done!
  end subroutine overlaps_gf
