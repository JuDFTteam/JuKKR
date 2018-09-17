module godfrin
  use :: mod_datatypes, only: dp
  ! Implementation of Godfrin's algorithm:
  ! EM Godfrin, J Phys Condens Matter 3 (1991) 7843-7848
  ! Extended to handle periodic boundary conditions
  ! Option to use standard sparse solver added

  implicit none

  type :: type_godfrin             ! GODFRIN Flaviano
    ! For the input parameters of the godfrin inversion scheme
    ! na: number of atoms, nb: number of blocks
    integer :: na, nb
    ! ldiag: diagonal part of inverse only, lper: periodic system, lpardiso:
    ! use PARDISO solver instead
    logical :: ldiag, lper, lpardiso
    ! bdims: block dimensions (how many atoms in each block)
    integer, allocatable :: bdims(:)
  end type type_godfrin

  type (type_godfrin), save :: t_godfrin ! GODFRIN Flaviano

  ! Only the driver is made public
  private
  public :: sparse_inverse, t_godfrin
#ifdef CPP_MPI
  public :: bcast_params_godfrin
#endif

  ! -----------------------------------------------------------------------
  ! Variables set by sparse_inverse
  ! Number of blocks of matrix A, maximum size of each block
  integer :: nb, mb
  ! Dimension of each block of matrix A
  integer, allocatable :: bdims(:)
  ! Storage for the non-zero elements of matrix A
  complex (kind=dp), allocatable :: aii(:, :, :), aij(:, :, :), aji(:, :, :)
  ! Storage for coefficients of recursion
  complex (kind=dp), allocatable :: x(:, :, :), y(:, :, :), tmp(:, :), c(:, :), d(:, :)
  ! Extra storage for periodic case
  complex (kind=dp), allocatable :: borders(:, :, :, :), delta(:, :)
  ! ----------------------------------------------------------------------
  ! Variables needed to call LAPACK routines
  integer :: info, lwork
  integer, allocatable :: ipiv(:)
  complex (kind=dp), allocatable :: work(:)
  ! ----------------------------------------------------------------------
  ! Some parameters
  complex (kind=dp), parameter :: czero = (0.e0_dp, 0.e0_dp), cone = (1.e0_dp, 0.e0_dp)
  complex (kind=dp), parameter :: cminus = (-1.e0_dp, 0.e0_dp)
  ! ----------------------------------------------------------------------
  ! Switch to print or not the running time
  logical, parameter :: ltiming = .false.


contains


  subroutine sparse_inverse(a, na, nblocks, blockdims, diagonal, periodic, use_pardiso)
    ! Main driver for sparse matrix inversion
    ! Typical dimensions:
    ! --> na = nlayers*nlms, nlayers number of layers, nlms is
    ! nspin*(lmax+1)^2
    ! If screened structure constants truncated after nscreen layers:
    ! --> na = (nlayers/nscreen)*(nscreen*nlms) = nblocks*blockdim
    ! assuming nlayers is divisible by nscreen
    ! In general na >= sum(blockdims(1:nblocks))

    ! ************************************************************************************
    ! WARNING: if the matrix a is not diagonally dominant this routine might
    ! be inaccurate
    ! ************************************************************************************

    implicit none


    ! Size of square matrix a in memory
    integer, intent (in) :: na
    ! Block tridiagonal matrix a to be inverted
    complex (kind=dp), intent (inout) :: a(na, na)
    ! Number of blocks of square matrix a
    integer, intent (in) :: nblocks
    ! Dimensions of each block
    integer, intent (in) :: blockdims(nblocks)
    ! Only diagonal blocks of inverse?
    logical, intent (in) :: diagonal
    ! Periodic boundary conditions?
    logical, intent (in) :: periodic
    ! Use MKL sparse solver instead?
    logical, intent (in) :: use_pardiso
    ! ----------------------------------------------------------------------
    ! Reallocate memory?
    logical :: reallocate
    ! Variables needed to run the PARDISO solver
    integer *8 :: pt(64)           ! this is 4 bytes for 32-bit OS and 8 bytes
    ! for 64-bit OS
    integer :: iparm(64), maxfct, mnum, phase, msglvl, mtype, n, ntot
    integer, allocatable :: ia(:), ja(:), perm(:)
    complex (kind=dp), allocatable :: asparse(:), b(:, :)
    ! Timing
    real (kind=dp) :: start, finish


    ! Initialize data structures
    nb = nblocks
    reallocate = .false.
    ! bdims not allocated => might be first call, allocate everything
    if (.not. allocated(bdims)) then
      allocate (bdims(nb))
      reallocate = .true.
      ! size of bdims has changed => new matrix size, reallocate everything
    else if (size(bdims)/=nblocks) then
      deallocate (bdims)
      allocate (bdims(nb))
      reallocate = .true.
    end if
    bdims = blockdims
    n = sum(bdims(1:nb))
    ! Inconsistency check
    if (na<sum(bdims(1:nb))) then
      write (*, '("sparse_inverse: na=",i8," sum(bdims)=",i8 )') na, sum(bdims(1:nb))
      stop
    end if
    ! ======================================================================
    ! GODFRIN ALGORITHM
#ifdef __INTEL_COMPILER
    if (.not. use_pardiso) then
#else
      ! deactivate pardiso if no intel compiler is used (see rinput13.f90 for
      ! further information)
      if (.true.) then
#endif
        ! ======================================================================
        ! Largest block size
        mb = maxval(bdims(1:nb))
        ! Storage for the non-zero blocks of A
        if (reallocate) allocate (aii(mb,mb,nb), aij(mb,mb,nb), aji(mb,mb,nb))
        ! Storage for recursion coefficients
        if (reallocate) allocate (x(mb,mb,nb), y(mb,mb,nb), tmp(mb,mb), c(mb,mb), d(mb,mb))
        ! Extra storage for periodic case
        if (reallocate .and. periodic) allocate (delta(mb,mb), borders(mb,mb,4,nb))
        ! ----------------------------------------------------------------------
        ! For the inverses: determine optimum size of work array
        if (reallocate) then
          lwork = -1
          allocate (ipiv(mb), work(1))
          call zgetri(mb, tmp, mb, ipiv, work, lwork, info)
          lwork = int(real(work(1)))
          deallocate (work)
          allocate (work(lwork))
        end if
        ! Save blocks of matrix in memory
        call cpu_time(start)
        call tridiag_save(a, na, periodic)
        call cpu_time(finish)
        if (ltiming) write (*, '("tridiag_save           time=",f10.3," s")') finish - start
        ! This computes the diagonal elements
        call cpu_time(start)
        call get_diagonal(a, na, periodic)
        call cpu_time(finish)
        if (ltiming) write (*, '("get_diagonal           time=",f10.3," s")') finish - start
        ! This fills in the rest
        if (.not. diagonal) then
          call cpu_time(start)
          call get_offdiagonal(a, na, periodic)
          call cpu_time(finish)
          if (ltiming) write (*, '("get_offdiagonal        time=",f10.3," s")') finish - start
        end if
        ! Corrections for periodic blocks
        if (periodic) then
          call cpu_time(start)
          call diagonal_correction(a, na)
          call cpu_time(finish)
          if (ltiming) write (*, '("diagonal_correction    time=",f10.3," s")') finish - start
          if (.not. diagonal) then
            call cpu_time(start)
            call offdiagonal_correction(a, na)
            call cpu_time(finish)
            if (ltiming) write (*, '("offdiagonal_correction time=",f10.3," s")') finish - start
          end if
        end if
        ! Free memory
        ! deallocate(aii,aij,aji,x,y,tmp,c,d)
        ! if (periodic) deallocate(delta,borders)
        ! ======================================================================
      else
        ! PARDISO SOLVER
        ! ======================================================================
        mtype = 3                  ! complex structurally symmetric
        maxfct = 1                 ! number of sparse factors in memory
        mnum = 1                   ! matrix to factorize
        phase = 13                 ! analysis, numerical factorization, solve,
        ! iterative refinement
        msglvl = 0                 ! write messages if 1
        ! Total number of non-zero elements (expected)
        ntot = sum(bdims(1:nb)**2) + sum(bdims(1:nb-1)**2) + sum(bdims(2:nb)**2)
        if (periodic) ntot = ntot + bdims(1)**2 + bdims(nb)**2
        call cpu_time(start)
        allocate (asparse(ntot), b(n,n), tmp(n,n), ia(n+1), ja(ntot), perm(n))
        call sparse_fill(n, a, na, asparse, b, nb, bdims, ntot, ia, ja, periodic)
        call pardisoinit(pt, mtype, iparm)
        call pardiso(pt, maxfct, mnum, mtype, phase, n, asparse, ia, ja, perm, n, iparm, msglvl, b, tmp, info)
        if (info/=0) stop 'fail in PARDISO'
        a(1:n, 1:n) = tmp
        deallocate (asparse, b, tmp, ia, ja, perm)
        call cpu_time(finish)
        if (ltiming) write (*, '("PARDISO sparse solver  time=",f10.3," s")') finish - start
        ! ======================================================================
      end if
      ! ======================================================================
      ! All done!
    end subroutine sparse_inverse
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine sparse_fill(k, a, na, asparse, b, nb, bdims, ntot, ia, ja, periodic)
      ! Copies the matrix A to storage in compressed sparse row format
      ! Fills in matrix B with unit matrix

      implicit none

      ! Number of rows/columns of A, number of non-zero blocks
      integer, intent (in) :: k
      ! The size of the matrix A in storage
      integer, intent (in) :: na
      ! The target matrix
      complex (kind=dp), intent (inout) :: a(na, na)
      ! Total number of non-zero elements
      integer, intent (in) :: ntot
      ! The non-zero blocks stored in compressed sparse row format
      complex (kind=dp), intent (out) :: asparse(ntot)
      ! Unit matrix
      complex (kind=dp), intent (out) :: b(k, k)
      ! Number of blocks per row/column
      integer, intent (in) :: nb
      ! Size of each block
      integer, intent (in) :: bdims(nb)
      ! Non-zero rows
      integer, intent (out) :: ia(k+1)
      ! Non-zero columns in each row
      integer, intent (out) :: ja(ntot)
      ! Periodic blocks
      logical, intent (in) :: periodic
      ! ----------------------------------------------------------------------
      integer :: i, ib, i0, j, jb, j0, m, n, itot

      itot = 0
      ! Rows
      ! write(*,*) "Debug: ib, i, jb, j, ia, ja, itot"
      do ib = 1, nb
        i0 = sum(bdims(1:ib-1))
        n = bdims(ib)
        do i = 1, n
          ia(i0+i) = itot + 1
          ! ------------------------------------------------------------------
          ! Columns
          ! A(N,1)
          if (periodic .and. ib==nb) then
            jb = 1
            j0 = sum(bdims(1:jb-1))
            m = bdims(jb)
            do j = 1, m
              ! compressed sparse row storage
              itot = itot + 1
              asparse(itot) = a(i0+i, j0+j)
              ja(itot) = j0 + j
              ! write(*,'(10i8)') ib, i, jb, j, ia(i0+i), ja(itot), itot
            end do
          end if
          ! A(i,j)
          do jb = max(1, ib-1), min(nb, ib+1)
            j0 = sum(bdims(1:jb-1))
            m = bdims(jb)
            do j = 1, m
              ! compressed sparse row storage
              itot = itot + 1
              asparse(itot) = a(i0+i, j0+j)
              ja(itot) = j0 + j
              ! write(*,'(10i8)') ib, i, jb, j, ia(i0+i), ja(itot), itot
            end do
          end do
          ! A(1,N)
          if (periodic .and. ib==1) then
            jb = nb
            j0 = sum(bdims(1:jb-1))
            m = bdims(jb)
            do j = 1, m
              ! compressed sparse row storage
              itot = itot + 1
              asparse(itot) = a(i0+i, j0+j)
              ja(itot) = j0 + j
              ! write(*,'(10i8)') ib, i, jb, j, ia(i0+i), ja(itot), itot
            end do
          end if
          ! ------------------------------------------------------------------
        end do
      end do
      ia(k+1) = itot + 1
      ! Goodbye A
      a = czero
      ! rhs for inversion
      b = czero
      do i = 1, k
        b(i, i) = cone
      end do
      ! write(*,*) itot
      ! write(*,'(1000i6)') ia(1:k+1)
      ! write(*,'(1000i6)') ja(1:ntot)
      ! write(*,'(1000es10.3)') asparse(1:ntot)
      ! write(*,'(1000es10.3)') b(1:k,1:k)
      ! All done!
    end subroutine sparse_fill
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine tridiag_save(a, na, periodic)
      ! Copy tridiagonals of matrix A

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      complex (kind=dp), intent (inout) :: a(na, na)
      ! Corner blocks?
      logical, intent (in) :: periodic
      ! -------------------------------
      integer :: i, i0, j0, m, n


      ! ------------------------------------------
      do i = 1, nb - 1
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        j0 = i0 + n
        m = bdims(i+1)
        ! diagonal blocks A(i,i)
        aii(1:n, 1:n, i) = a(i0+1:i0+n, i0+1:i0+n)
        ! superdiagonal blocks A(i,i+1)
        aij(1:n, 1:m, i) = a(i0+1:i0+n, j0+1:j0+n)
        ! subdiagonal blocks A(i+1,i)
        aji(1:m, 1:n, i) = a(j0+1:j0+m, i0+1:i0+n)
      end do
      ! Last diagonal block
      i0 = sum(bdims(1:nb-1))
      n = bdims(nb)
      aii(1:n, 1:n, nb) = a(i0+1:i0+n, i0+1:i0+n)
      ! ------------------------------------------
      if (periodic) then
        ! Blocks connecting first and last layers
        i0 = 0
        n = bdims(1)
        j0 = sum(bdims(1:nb-1))
        m = bdims(nb)
        ! superdiagonal block A(1,N)
        aij(1:n, 1:m, nb) = a(i0+1:i0+n, j0+1:j0+m)
        ! subdiagonal block A(N,1)
        aji(1:m, 1:n, nb) = a(j0+1:j0+m, i0+1:i0+n)
        ! First and last blocks of the diagonal are modified
        ! A(1,1) = A(1,1) - 1
        do i = 1, n
          aii(i, i, 1) = aii(i, i, 1) + cminus
        end do
        ! A(N,N) = A(N,N) - A(N,1) A(1,N)
        call zgemm('N', 'N', m, m, n, cminus, aji(:,:,nb), mb, aij(:,:,nb), mb, cone, aii(:,:,nb), mb)
      end if
      ! Goodbye input matrix
      a = czero
      ! All done
    end subroutine tridiag_save
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine get_diagonal(a, na, periodic)
      ! Fills in the diagonal of the matrix inverse
      ! Also computes the first/last row/column for periodic blocks

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      complex (kind=dp), intent (inout) :: a(na, na)
      ! Corner blocks?
      logical, intent (in) :: periodic
      ! -------------------------------
      integer :: i, i0, n, m


      ! ----------------------------------------------------------------------
      ! Recursion for Y matrices
      ! Y(1) = 0
      n = bdims(1)
      y(1:n, 1:n, 1) = czero
      do i = 2, nb
        m = bdims(i-1)
        n = bdims(i)
        ! Y(i) = (A(i-1,i-1) - Y(i-1))^-1
        y(1:m, 1:m, i) = aii(1:m, 1:m, i-1) - y(1:m, 1:m, i-1)
        ! ***** Inverse *****
        call zgetrf(m, m, y(:,:,i), mb, ipiv, info)
        if (info/=0) stop 'get_diagonal: LU fail in Y'
        call zgetri(m, y(:,:,i), mb, ipiv, work, lwork, info)
        if (info/=0) stop 'get_diagonal: inverse fail in Y'
        ! ******************
        ! tmp = Y(i) * A(i-1,i)  -- superdiagonal
        call zgemm('N', 'N', m, n, m, cone, y(:,:,i), mb, aij(:,:,i-1), mb, czero, tmp, mb)
        ! Y(i) = A(i,i-1) * tmp  -- subdiagonal
        call zgemm('N', 'N', n, n, m, cone, aji(:,:,i-1), mb, tmp, mb, czero, y(:,:,i), mb)
        ! call out_mat(y(1:n,1:n,i),n,n,i,0,'Y(i)  ')
      end do
      ! ----------------------------------------------------------------------
      ! Recursion for X matrices
      ! X(n) = 0
      n = bdims(nb)
      x(1:n, 1:n, nb) = czero
      do i = nb - 1, 1, -1
        m = bdims(i+1)
        n = bdims(i)
        ! X(i) = (A(i+1,i+1) - X(i+1))^-1
        x(1:m, 1:m, i) = aii(1:m, 1:m, i+1) - x(1:m, 1:m, i+1)
        ! ***** Inverse *****
        call zgetrf(m, m, x(:,:,i), mb, ipiv, info)
        if (info/=0) stop 'get_diagonal: LU fail in X'
        call zgetri(m, x(:,:,i), mb, ipiv, work, lwork, info)
        if (info/=0) stop 'get_diagonal: inverse fail in X'
        ! ******************
        ! tmp = X(i) * A(i+1,i)  -- subdiagonal
        call zgemm('N', 'N', m, n, m, cone, x(:,:,i), mb, aji(:,:,i), mb, czero, tmp, mb)
        ! X(i) = A(i,i+1) * tmp  -- superdiagonal
        call zgemm('N', 'N', n, n, m, cone, aij(:,:,i), mb, tmp, mb, czero, x(:,:,i), mb)
        ! call out_mat(x(1:n,1:n,i),n,n,i,0,'X(i)  ')
      end do
      ! ----------------------------------------------------------------------
      ! Diagonal elements of inverse
      do i = 1, nb
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        ! B(i,i) = (A(i,i) - X(i) - Y(i))^-1
        tmp(1:n, 1:n) = aii(1:n, 1:n, i) - x(1:n, 1:n, i) - y(1:n, 1:n, i)
        ! ***** Inverse *****
        call zgetrf(n, n, tmp, mb, ipiv, info)
        if (info/=0) stop 'get_diagonal: LU fail in B(i,i)'
        call zgetri(n, tmp, mb, ipiv, work, lwork, info)
        if (info/=0) stop 'get_diagonal: inverse fail in B(i,i)'
        ! *******************
        ! call out_mat(tmp(1:n,1:n),n,n,i,i,'B(i,i)')
        a(i0+1:i0+n, i0+1:i0+n) = tmp(1:n, 1:n)
      end do
      ! ----------------------------------------------------------------------
      ! All done!
    end subroutine get_diagonal
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine get_offdiagonal(a, na, periodic)
      ! Fills in the off-diagonal elements of the matrix inverse
      ! It includes the correction for periodic blocks

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      complex (kind=dp), intent (inout) :: a(na, na)
      ! Corner blocks?
      logical, intent (in) :: periodic
      ! -------------------------------
      logical, parameter :: columns = .true.
      integer :: i, j, istop, jstop


      ! Off-diagonal elements of inverse without corner blocks
      if (columns) then
        jstop = 0
        ! First and last columns were already calculated
        if (periodic) jstop = 1
        ! Build one column at a time
        do j = 1 + jstop, nb - jstop
          call build_col(a, na, j, jstop)
        end do
      else
        istop = 0
        ! First and last rows were already calculated
        if (periodic) istop = 1
        ! Build one row at a time
        do i = 1 + istop, nb - istop
          call build_row(a, na, i, istop)
        end do
      end if
      ! All done!
    end subroutine get_offdiagonal
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine diagonal_correction(a, na)
      ! Corrects the diagonal blocks for the periodic case

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      complex (kind=dp), intent (inout) :: a(na, na)
      ! -------------------------------
      integer :: i, i0, n, k


      ! ----------------------------------------------------------------------
      ! Save borders of inverse matrix with block prefactors
      call save_borders(a, na)
      ! -----------------------------------------------------------------------
      ! Update the diagonals of the inverse matrix with the correction
      ! Using rescaled border rows and columns (see save_borders)
      k = bdims(1)
      do i = 1, nb
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        ! tmp = B(i,i)
        tmp(1:n, 1:n) = a(i0+1:i0+n, i0+1:i0+n)
        ! tmp = tmp - B(i,1) * B(1,i)
        call zgemm('N', 'N', n, n, k, cminus, borders(:,:,1,i), mb, borders(:,:,3,i), mb, cone, tmp, mb)
        ! tmp = tmp - B(i,1) * B(N,i)
        call zgemm('N', 'N', n, n, k, cminus, borders(:,:,1,i), mb, borders(:,:,4,i), mb, cone, tmp, mb)
        ! tmp = tmp - B(i,N) * B(1,i)
        call zgemm('N', 'N', n, n, k, cminus, borders(:,:,2,i), mb, borders(:,:,3,i), mb, cone, tmp, mb)
        ! tmp = tmp - B(i,N) * B(N,i)
        call zgemm('N', 'N', n, n, k, cminus, borders(:,:,2,i), mb, borders(:,:,4,i), mb, cone, tmp, mb)
        ! B(i,i) = tmp
        a(i0+1:i0+n, i0+1:i0+n) = tmp(1:n, 1:n)
      end do
      ! All done!
    end subroutine diagonal_correction
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine offdiagonal_correction(a, na)
      ! Corrects the off-diagonal blocks for the periodic case

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      complex (kind=dp), intent (inout) :: a(na, na)
      ! -------------------------------
      integer :: i, i0, n, j, j0, m, k


      ! Update the off-diagonals of the inverse matrix with the correction
      ! Using rescaled border rows and columns (see save_borders)
      k = bdims(1)
      do j = 1, nb
        j0 = sum(bdims(1:j-1))
        m = bdims(j)
        do i = 1, nb
          if (i==j) cycle          ! diagonals were taken care of already
          i0 = sum(bdims(1:i-1))
          n = bdims(i)
          ! tmp = B(i,j)
          tmp(1:n, 1:m) = a(i0+1:i0+n, j0+1:j0+m)
          ! tmp = tmp - B(i,1) * B(1,j)
          call zgemm('N', 'N', n, m, k, cminus, borders(:,:,1,i), mb, borders(:,:,3,j), mb, cone, tmp, mb)
          ! tmp = tmp - B(i,1) * B(N,j)
          call zgemm('N', 'N', n, m, k, cminus, borders(:,:,1,i), mb, borders(:,:,4,j), mb, cone, tmp, mb)
          ! tmp = tmp - B(i,N) * B(1,j)
          call zgemm('N', 'N', n, m, k, cminus, borders(:,:,2,i), mb, borders(:,:,3,j), mb, cone, tmp, mb)
          ! tmp = tmp - B(i,N) * B(N,j)
          call zgemm('N', 'N', n, m, k, cminus, borders(:,:,2,i), mb, borders(:,:,4,j), mb, cone, tmp, mb)
          ! B(i,j) = tmp
          a(i0+1:i0+n, j0+1:j0+m) = tmp(1:n, 1:m)
        end do
      end do
      ! All done!
    end subroutine offdiagonal_correction
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine save_borders(a, na)
      ! Saves the borders of the inverse matrix without corner blocks
      ! These are multiplied by appropriate blocks for later use
      ! The column blocks have dimensions bdims(i) * bdims(1)
      ! The row blocks have dimensions bdims(1) * bdims(j)

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      ! Should have the diagonals of the inverse already stored
      complex (kind=dp), intent (inout) :: a(na, na)
      ! -------------------------------
      integer :: i, j, i0, j0, m, n, k


      ! ----------------------------------------------------------------------
      ! Generate the borders of the matrix inverse
      call build_col(a, na, 1, 0)
      call build_col(a, na, nb, 0)
      call build_row(a, na, 1, 0)
      call build_row(a, na, nb, 0)
      ! ----------------------------------------------------------------------
      ! Central part of the Sherman-Morrison-Woodbury formula
      i0 = 0
      n = bdims(1)
      j0 = sum(bdims(1:nb-1))
      m = bdims(nb)
      ! Delta = A(1,N) * B(N,N) * A(N,1)
      delta(1:m, 1:m) = a(j0+1:j0+m, j0+1:j0+m)
      call zgemm('N', 'N', m, n, m, cone, delta, mb, aji(:,:,nb), mb, czero, tmp, mb)
      call zgemm('N', 'N', n, n, m, cone, aij(:,:,nb), mb, tmp, mb, czero, delta, mb)
      ! Delta = Delta + A(1,N) * B(N,1)
      tmp = a(j0+1:j0+m, i0+1:i0+n)
      call zgemm('N', 'N', n, n, m, cone, aij(:,:,nb), mb, tmp, mb, cone, delta, mb)
      ! Delta = Delta + B(1,N) * A(N,1)
      tmp = a(i0+1:i0+n, j0+1:j0+m)
      call zgemm('N', 'N', n, n, m, cone, tmp, mb, aji(:,:,nb), mb, cone, delta, mb)
      ! Delta = Delta + B(1,1) + 1
      delta(1:n, 1:n) = delta(1:n, 1:n) + a(i0+1:i0+n, i0+1:i0+n)
      do i = 1, n
        delta(i, i) = delta(i, i) + cone
      end do
      ! ***** Inverse *****
      call zgetrf(n, n, delta, mb, ipiv, info)
      if (info/=0) stop 'save_borders: LU fail in Delta'
      call zgetri(n, delta, mb, ipiv, work, lwork, info)
      if (info/=0) stop 'save_borders: inverse fail in Delta'
      ! *******************
      ! call out_mat(delta(1:n,1:n),n,n,1,1,'Delta ')
      ! ----------------------------------------------------------------------
      ! First column: B(i,1)
      j0 = 0
      m = bdims(1)
      do i = 1, nb
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        borders(1:n, 1:m, 1, i) = a(i0+1:i0+n, j0+1:j0+m) ! B(i,1)
        ! call out_mat(borders(1:n,1:m,1,i),n,m,i,1,'B(i,1)')
      end do
      ! ----------------------------------------------------------------------
      ! Last column: B(i,N) * A(N,1)
      j0 = sum(bdims(1:nb-1))
      m = bdims(nb)
      k = bdims(1)
      do i = 1, nb
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        tmp(1:n, 1:m) = a(i0+1:i0+n, j0+1:j0+m) ! B(i,N)
        call zgemm('N', 'N', n, k, m, cone, tmp, mb, aji(:,:,nb), mb, czero, borders(:,:,2,i), mb)
        ! call out_mat(borders(1:n,1:k,2,i),n,k,i,nb,'B(i,N)')
      end do
      ! ----------------------------------------------------------------------
      ! First row: Delta * B(1,j)
      i0 = 0
      n = bdims(1)
      do j = 1, nb
        j0 = sum(bdims(1:j-1))
        m = bdims(j)
        tmp(1:n, 1:m) = a(i0+1:i0+n, j0+1:j0+m) ! B(1,j)
        call zgemm('N', 'N', n, m, n, cone, delta, mb, tmp, mb, czero, borders(:,:,3,j), mb)
        ! call out_mat(borders(1:n,1:m,3,j),n,m,1,j,'B(1,j)')
      end do
      ! ----------------------------------------------------------------------
      ! Last row: Delta * A(1,N) * B(N,j)
      i0 = sum(bdims(1:nb-1))
      n = bdims(nb)
      k = bdims(1)
      do j = 1, nb
        j0 = sum(bdims(1:j-1))
        m = bdims(j)
        borders(1:n, 1:m, 4, j) = a(i0+1:i0+n, j0+1:j0+m) ! B(N,j)
        call zgemm('N', 'N', k, m, n, cone, aij(:,:,nb), mb, borders(:,:,4,j), mb, czero, tmp, mb)
        call zgemm('N', 'N', k, m, k, cone, delta, mb, tmp, mb, czero, borders(:,:,4,j), mb)
        ! call out_mat(borders(1:k,1:m,4,j),k,m,j,nb,'B(N,j)')
      end do
      ! ----------------------------------------------------------------------
      ! All done!
    end subroutine save_borders
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine build_col(a, na, j, jstop)
      ! Builds the j-th column of the inverse starting from the diagonal

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      complex (kind=dp), intent (inout) :: a(na, na)
      ! Desired column, where to stop
      integer, intent (in) :: j, jstop
      ! -------------------------------
      integer :: i, i0, n, j0, m, k


      ! ----------------------------------------------------------------------
      ! Subdiagonal blocks (i > j)
      ! B(i,j) = - C(i) * A(i,i-1) * B(i-1,j)
      ! ----------------------------------------------------------------------
      j0 = sum(bdims(1:j-1))
      m = bdims(j)
      k = m
      ! B(j,j)  -->  saved in array D
      d(1:m, 1:m) = a(j0+1:j0+m, j0+1:j0+m)
      do i = j + 1, nb - jstop
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        ! tmp = - A(i,i-1) * B(i-1,j)  --> subdiagonal
        call zgemm('N', 'N', n, m, k, cminus, aji(:,:,i-1), mb, d, mb, czero, tmp, mb)
        ! C(i) = (A(i,i) - X(i))^-1
        c(1:n, 1:n) = aii(1:n, 1:n, i) - x(1:n, 1:n, i)
        ! ***** Inverse *****
        call zgetrf(n, n, c, mb, ipiv, info)
        if (info/=0) stop 'build_col: LU fail in C(i)'
        call zgetri(n, c, mb, ipiv, work, lwork, info)
        if (info/=0) stop 'build_col: inverse fail in C(i)'
        ! *******************
        ! B(i,j) = C(i) * tmp  -->  D used as temporary storage
        call zgemm('N', 'N', n, m, n, cone, c, mb, tmp, mb, czero, d, mb)
        ! call out_mat(d(1:n,1:m),n,m,i,j,'B(i,j)')
        a(i0+1:i0+n, j0+1:j0+m) = d(1:n, 1:m)
        k = n                      ! save old row dimension
      end do
      ! ----------------------------------------------------------------------
      ! Superdiagonal blocks (i < j)
      ! B(i,j) = - D(i) * A(i,i+1) * B(i+1,j)
      ! ----------------------------------------------------------------------
      k = m
      ! B(j,j)  -->  saved in array C
      c(1:m, 1:m) = a(j0+1:j0+m, j0+1:j0+m)
      do i = j - 1, 1 + jstop, -1
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        ! tmp = - A(i,i+1) * B(i+1,j)  -->  superdiagonal
        call zgemm('N', 'N', n, m, k, cminus, aij(:,:,i), mb, c, mb, czero, tmp, mb)
        ! D(i) = (A(i,i) - Y(i))^-1
        d(1:n, 1:n) = aii(1:n, 1:n, i) - y(1:n, 1:n, i)
        ! ***** Inverse *****
        call zgetrf(n, n, d, mb, ipiv, info)
        if (info/=0) stop 'build_col: LU fail in D(i)'
        call zgetri(n, d, mb, ipiv, work, lwork, info)
        if (info/=0) stop 'build_col: inverse fail in D(i)'
        ! *******************
        ! B(i,j) = D(i) * tmp  -->  C used as temporary storage
        call zgemm('N', 'N', n, m, n, cone, d, mb, tmp, mb, czero, c, mb)
        ! call out_mat(c(1:n,1:m),n,m,i,j,'B(i,j)')
        a(i0+1:i0+n, j0+1:j0+m) = c(1:n, 1:m)
        k = n                      ! save old row dimension
      end do
      ! All done!
    end subroutine build_col
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine build_row(a, na, i, istop)
      ! Builds the i-th row of the inverse starting from the diagonal

      implicit none

      ! Size of square matrix a in memory
      integer, intent (in) :: na
      ! Block tridiagonal matrix a to be inverted
      complex (kind=dp), intent (inout) :: a(na, na)
      ! Desired column, where to stop
      integer, intent (in) :: i, istop
      ! -------------------------------
      integer :: i0, n, j, j0, m, k


      ! ----------------------------------------------------------------------
      ! Subdiagonal blocks (i > j)
      ! B(i,j) = - B(i,j+1) * A(j+1,j) * D(j)
      ! ----------------------------------------------------------------------
      i0 = sum(bdims(1:i-1))
      n = bdims(i)
      k = n
      ! B(i,i)  -->  saved in array C
      c(1:n, 1:n) = a(i0+1:i0+n, i0+1:i0+n)
      do j = i - 1, 1 + istop, -1
        j0 = sum(bdims(1:j-1))
        m = bdims(j)
        ! tmp = - B(i,j+1) A(j+1,j)  -->  subdiagonal
        call zgemm('N', 'N', n, m, k, cminus, c, mb, aji(:,:,j), mb, czero, tmp, mb)
        ! D(i) = (A(i,i) - Y(i))^-1
        d(1:m, 1:m) = aii(1:m, 1:m, j) - y(1:m, 1:m, j)
        ! ***** Inverse *****
        call zgetrf(m, m, d, mb, ipiv, info)
        if (info/=0) stop 'build_row: LU fail in D(i)'
        call zgetri(m, d, mb, ipiv, work, lwork, info)
        if (info/=0) stop 'build_row: inverse fail in D(i)'
        ! *******************
        ! B(i,j) = tmp * D(j)
        call zgemm('N', 'N', n, m, m, cone, tmp, mb, d, mb, czero, c, mb)
        ! call out_mat(c(1:n,1:m),n,m,i,j,'B(i,j)')
        a(i0+1:i0+n, j0+1:j0+m) = c(1:n, 1:m)
        k = m                      ! save old column dimension
      end do
      ! ----------------------------------------------------------------------
      ! Superdiagonal blocks (i < j)
      ! B(i,j) = - B(i,j-1) * A(j-1,j) * C(j)
      ! ----------------------------------------------------------------------
      k = n
      ! B(i,i)  -->  saved in array D
      d(1:n, 1:n) = a(i0+1:i0+n, i0+1:i0+n)
      do j = i + 1, nb - istop
        j0 = sum(bdims(1:j-1))
        m = bdims(j)
        ! tmp = - B(i,j-1) * A(j-1,j)  -->  superdiagonal
        call zgemm('N', 'N', n, m, k, cminus, d, mb, aij(:,:,j-1), mb, czero, tmp, mb)
        ! C(j) = (A(j,j) - X(j))^-1
        c(1:m, 1:m) = aii(1:m, 1:m, j) - x(1:m, 1:m, j)
        ! ***** Inverse *****
        call zgetrf(m, m, c, mb, ipiv, info)
        if (info/=0) stop 'build_row: LU fail in C(i)'
        call zgetri(m, c, mb, ipiv, work, lwork, info)
        if (info/=0) stop 'build_row: inverse fail in C(i)'
        ! *******************
        ! B(i,j) = tmp * C(j) -->  D used as temporary storage
        call zgemm('N', 'N', n, m, m, cone, tmp, mb, c, mb, czero, d, mb)
        ! call out_mat(d(1:n,1:m),n,m,i,j,'B(i,j)')
        a(i0+1:i0+n, j0+1:j0+m) = d(1:n, 1:m)
        k = m                      ! save old column dimension
      end do
      ! All done!
    end subroutine build_row
    ! -----------------------------------------------------------------------


    ! -----------------------------------------------------------------------
    subroutine out_mat(a, n, m, ib, jb, label)

      implicit none

      ! Number of rows and columns
      integer, intent (in) :: n, m
      ! Matrix to output
      complex (kind=dp), intent (in) :: a(n, m)
      ! Info: which block, text label
      integer, intent (in) :: ib, jb
      character (len=6), intent (in) :: label
      ! ------------------------------------
      integer :: i

      write (*, *) label, ib, jb
      do i = 1, n
        write (*, '(1000es10.3)') a(i, 1:m)
      end do
      ! All done!
    end subroutine out_mat
    ! -----------------------------------------------------------------------

#ifdef CPP_MPI
    subroutine bcast_params_godfrin(t_godfrin) ! GODFRIN Flaviano
      ! broadcast parameters for godfrin matrix inversion scheme
      use :: mpi
      implicit none

      type (type_godfrin), intent (inout) :: t_godfrin
      integer :: ierr

      call mpi_bcast(t_godfrin%na, 1, mpi_integer, 0, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting maxmem_number'
      call mpi_bcast(t_godfrin%nb, 1, mpi_integer, 0, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting maxmem_units'
      call mpi_bcast(t_godfrin%ldiag, 1, mpi_logical, 0, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_rll'
      call mpi_bcast(t_godfrin%lper, 1, mpi_logical, 0, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_sll'
      call mpi_bcast(t_godfrin%lpardiso, 1, mpi_logical, 0, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_rllleft'

      if (.not. allocated(t_godfrin%bdims)) allocate (t_godfrin%bdims(t_godfrin%nb))
      call mpi_bcast(t_godfrin%bdims, t_godfrin%nb, mpi_integer, 0, mpi_comm_world, ierr)
      if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_sllleft'

    end subroutine bcast_params_godfrin
#endif

  end module godfrin
