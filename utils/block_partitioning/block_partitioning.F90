  program block_partitioning
! Given a banded matrix sparsity pattern,
! an optimal partition into block tridiagonal structure is found

  implicit none

  integer, allocatable :: links(:,:), bsizes(:)
  integer :: i, nc, nb

! Read sparsity pattern
  open(file='couplings.dat',unit=10,status='old')
  read(10,*)  ! skip line
  read(10,*) nc
  allocate(links(nc,nc),bsizes(nc))
  do i=1,nc
    read(10,'(1000i1)') links(1:nc,i)
  end do
  close(10)
  call block_partitioner(nc,links,nb,bsizes)
  deallocate(links)


  contains


  subroutine block_partitioner(nc,links,nblocks,blocksizes)
! Find an optimal block tridiagonal partition
! for sparsity pattern given in links

! Number of columns of matrix
  integer, intent(in)  :: nc
! Location of non-zero entries
  integer, intent(in)  :: links(nc,nc)
! Optimal number of blocks
  integer, intent(out) :: nblocks
! Sizes of each block, in (1:nblocks)
  integer, intent(out) :: blocksizes(nc)
! ----------------------------------------------------------------------
  integer :: trialsizes(nc)
  integer :: i, j, iup, idn, maxblock, minsize, sum2, minsum2, imax
  logical :: lcheck
  integer :: bandmin, bandmax, bandavg, banddev

! Initial values
  nblocks = 1; blocksizes(1) = nc; blocksizes(2:nc) = 0
! Find bandwidth of each column
  bandmax = 0; bandmin = nc; bandavg = 0; banddev = 0
  do j=1,nc
    iup = 0
    do i=j,1,-1
      if (links(i,j) == 1) iup = i
    end do
    idn = 0
    do i=j,nc
      if (links(i,j) == 1) idn = i
    end do
    if (idn-iup < bandmin) bandmin = idn-iup
    if (idn-iup > bandmax) bandmax = idn-iup
    bandavg = bandavg + idn - iup
    banddev = banddev + (idn - iup)**2
!    write(*,'("bandwidth: j, up, dn, diff=",4i8)') j, iup, idn, idn-iup
  end do
  bandavg = bandavg/nc
  banddev = ceiling(sqrt(real(banddev/nc - bandavg**2)))
  write(*,'("bandwidth: min, max, avg, dev=",4i8)') bandmin, bandmax, bandavg, banddev
  call brute_partitioner(nc,links,bandmin,bandmax,bandavg,banddev,nblocks,blocksizes)
! All done!
  end subroutine block_partitioner


  subroutine brute_partitioner(nc,links,bandmin,bandmax,bandavg,banddev,nblocks,blocksizes)
! Brute force search for optimal block partition

  implicit none

! Size of full matrix
  integer, intent(in)  :: nc
! Sparsity pattern
  integer, intent(in)  :: links(nc,nc)
! Bandwidth statistics
  integer, intent(in)  :: bandmin, bandmax, bandavg, banddev
! Optimal number of blocks
  integer, intent(out) :: nblocks
! Sizes of each block, in (1:nblocks)
  integer, intent(out) :: blocksizes(nc)
! ----------------------------------------------------------------------
  integer :: trialsizes(nc)
  integer :: i, k, maxblock, minsize, sum2, minsum2, imax
  logical :: lcheck
  integer*8 :: j, jdiv, jrem, imax2i

! Minimum block size is half of the minimum bandwith
!  minsize = bandmin/2
! Allowance for uneven bandwidth
  minsize = (bandmin+banddev)/2
! Maximum number of blocks
  maxblock = nc/minsize
  write(*,'("minsize, maxblock=",2i8)') minsize, maxblock
! Find combinations of block sizes that add to nc
  minsum2 = nc**2
  do i=2,maxblock
!   Sum of squares of block dimensions is bound by nc**2
    imax = ceiling(sqrt(1.d0*(nc**2 - (i-1)*minsize**2))) - minsize
!   Sum of block dimensions is bound by nc
    imax = min(imax,max(nc-i*minsize,1))
!   Empirical reduction
    imax = min(imax,max(bandmax/2-minsize,1))
!   Make sure the big integers don't explode
    imax2i = imax**(i-1)
    write(*,'("nblocks,imax=",2i8,i16)') i, imax, imax2i
!   Go over all combinations, encoded in j
!   Digits are decoded into jdiv
    do j=0,imax2i
      jrem = j
      do k=i-1,1,-1
        jdiv = jrem/(imax**(k-1))
        jrem = mod(jrem,imax**(k-1))
        trialsizes(k) = minsize + jdiv
      end do
      trialsizes(i) = nc - sum(trialsizes(1:i-1))
      if (trialsizes(i) < minsize) cycle
!     Do the blocks leave out any matrix elements?
      call check_blocks(nc,links,i,trialsizes(1:i),lcheck)
      if (lcheck) then
        sum2 = sum(trialsizes(1:i)**2)
        if (sum2 < minsum2) then
          minsum2 = sum2
          nblocks = i
          blocksizes(1:i) = trialsizes(1:i)
        end if
      end if
    end do
    write(*,'("minsum2, nblocks, blocksizes=",1000i8)') minsum2, nblocks, blocksizes(1:nblocks)
  end do
! All done!
  end subroutine brute_partitioner


  pure subroutine check_blocks(nc,links,nblocks,blocksizes,lcheck)
! Check if the given blocks are a tridiagonal partition

  implicit none

! Number of columns
  integer, intent(in)  :: nc
! Sparsity pattern
  integer, intent(in)  :: links(nc,nc)
! Number of blocks
  integer, intent(in)  :: nblocks
! Block sizes
  integer, intent(in)  :: blocksizes(1:nblocks)
! Is this a correct partition?
  logical, intent(out) :: lcheck
! ----------------------------------------------------------------------
  integer :: ib, i0, i1, j0, j1
  integer :: zeroed(nc,nc)

  lcheck = .true.
  zeroed(:,:) = links(:,:)
  i0 = 1; i1 = blocksizes(1)
  do ib=1,nblocks-1
    i0 = sum(blocksizes(1:ib-1)) + 1
    i1 = sum(blocksizes(1:ib))
    j0 = i1 + 1
    j1 = i1 + blocksizes(ib+1)
    zeroed(i0:i1,i0:i1) = 0
    zeroed(i0:i1,j0:j1) = 0
    zeroed(j0:j1,i0:i1) = 0
  end do
  zeroed(j0:j1,j0:j1) = 0
  if (any(zeroed(:,:) == 1)) lcheck = .false.
! All done!
  end subroutine check_blocks


! All done!
  end program block_partitioning
