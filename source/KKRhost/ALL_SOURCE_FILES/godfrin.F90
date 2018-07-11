Module godfrin
  Use mod_datatypes, Only: dp
  ! Implementation of Godfrin's algorithm:
  ! EM Godfrin, J Phys Condens Matter 3 (1991) 7843-7848
  ! Extended to handle periodic boundary conditions
  ! Option to use standard sparse solver added

  Implicit None

  Type :: type_godfrin   ! GODFRIN Flaviano
    ! For the input parameters of the godfrin inversion scheme
    ! na: number of atoms, nb: number of blocks
    Integer :: na, nb
    ! ldiag: diagonal part of inverse only, lper: periodic system, lpardiso:
    ! use PARDISO solver instead
    Logical :: ldiag, lper, lpardiso
    ! bdims: block dimensions (how many atoms in each block)
    Integer, Allocatable :: bdims(:)
  End Type

  Type (type_godfrin), Save :: t_godfrin   ! GODFRIN Flaviano

  ! Only the driver is made public
  Private
  Public :: sparse_inverse, t_godfrin
#ifdef CPP_MPI
  Public :: bcast_params_godfrin
#endif

  ! -----------------------------------------------------------------------
  ! Variables set by sparse_inverse
  ! Number of blocks of matrix A, maximum size of each block
  Integer :: nb, mb
  ! Dimension of each block of matrix A
  Integer, Allocatable :: bdims(:)
  ! Storage for the non-zero elements of matrix A
  Complex (Kind=dp), Allocatable :: aii(:, :, :), aij(:, :, :), &
    aji(:, :, :)
  ! Storage for coefficients of recursion
  Complex (Kind=dp), Allocatable :: x(:, :, :), y(:, :, :), tmp(:, :), &
    c(:, :), d(:, :)
  ! Extra storage for periodic case
  Complex (Kind=dp), Allocatable :: borders(:, :, :, :), delta(:, :)
  ! ----------------------------------------------------------------------
  ! Variables needed to call LAPACK routines
  Integer :: info, lwork
  Integer, Allocatable :: ipiv(:)
  Complex (Kind=dp), Allocatable :: work(:)
  ! ----------------------------------------------------------------------
  ! Some parameters
  Complex (Kind=dp), Parameter :: czero = (0.E0_dp, 0.E0_dp), &
    cone = (1.E0_dp, 0.E0_dp)
  Complex (Kind=dp), Parameter :: cminus = (-1.E0_dp, 0.E0_dp)
  ! ----------------------------------------------------------------------
  ! Switch to print or not the running time
  Logical, Parameter :: ltiming = .False.


Contains


  Subroutine sparse_inverse(a, na, nblocks, blockdims, diagonal, periodic, &
    use_pardiso)
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

    Implicit None


    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! Number of blocks of square matrix a
    Integer, Intent (In) :: nblocks
    ! Dimensions of each block
    Integer, Intent (In) :: blockdims(nblocks)
    ! Only diagonal blocks of inverse?
    Logical, Intent (In) :: diagonal
    ! Periodic boundary conditions?
    Logical, Intent (In) :: periodic
    ! Use MKL sparse solver instead?
    Logical, Intent (In) :: use_pardiso
    ! ----------------------------------------------------------------------
    ! Reallocate memory?
    Logical :: reallocate
    ! Variables needed to run the PARDISO solver
    Integer *8 :: pt(64)   ! this is 4 bytes for 32-bit OS and 8 bytes
      ! for 64-bit OS
    Integer :: iparm(64), maxfct, mnum, phase, msglvl, mtype, n, ntot
    Integer, Allocatable :: ia(:), ja(:), perm(:)
    Complex (Kind=dp), Allocatable :: asparse(:), b(:, :)
    ! Timing
    Real (Kind=dp) :: start, finish


    ! Initialize data structures
    nb = nblocks
    reallocate = .False.
    ! bdims not allocated => might be first call, allocate everything
    If (.Not. allocated(bdims)) Then
      Allocate (bdims(nb))
      reallocate = .True.
    ! size of bdims has changed => new matrix size, reallocate everything
    Else If (size(bdims)/=nblocks) Then
      Deallocate (bdims)
      Allocate (bdims(nb))
      reallocate = .True.
    End If
    bdims = blockdims
    n = sum(bdims(1:nb))
    ! Inconsistency check
    If (na<sum(bdims(1:nb))) Then
      Write (*, '("sparse_inverse: na=",i8," sum(bdims)=",i8 )') na, &
        sum(bdims(1:nb))
      Stop
    End If
    ! ======================================================================
    ! GODFRIN ALGORITHM
#ifdef __INTEL_COMPILER
    If (.Not. use_pardiso) Then
#else
    ! deactivate pardiso if no intel compiler is used (see rinput13.f90 for
    ! further information)
    If (.true.) then
#endif
      ! ======================================================================
      ! Largest block size
      mb = maxval(bdims(1:nb))
      ! Storage for the non-zero blocks of A
      If (reallocate) Allocate (aii(mb,mb,nb), aij(mb,mb,nb), &
        aji(mb,mb,nb))
      ! Storage for recursion coefficients
      If (reallocate) Allocate (x(mb,mb,nb), y(mb,mb,nb), tmp(mb,mb), &
        c(mb,mb), d(mb,mb))
      ! Extra storage for periodic case
      If (reallocate .And. periodic) Allocate (delta(mb,mb), &
        borders(mb,mb,4,nb))
      ! ----------------------------------------------------------------------
      ! For the inverses: determine optimum size of work array
      If (reallocate) Then
        lwork = -1
        Allocate (ipiv(mb), work(1))
        Call zgetri(mb, tmp, mb, ipiv, work, lwork, info)
        lwork = int(real(work(1)))
        Deallocate (work)
        Allocate (work(lwork))
      End If
      ! Save blocks of matrix in memory
      Call cpu_time(start)
      Call tridiag_save(a, na, periodic)
      Call cpu_time(finish)
      If (ltiming) Write (*, '("tridiag_save           time=",f10.3," s")' &
        ) finish - start
      ! This computes the diagonal elements
      Call cpu_time(start)
      Call get_diagonal(a, na, periodic)
      Call cpu_time(finish)
      If (ltiming) Write (*, '("get_diagonal           time=",f10.3," s")' &
        ) finish - start
      ! This fills in the rest
      If (.Not. diagonal) Then
        Call cpu_time(start)
        Call get_offdiagonal(a, na, periodic)
        Call cpu_time(finish)
        If (ltiming) Write (*, &
          '("get_offdiagonal        time=",f10.3," s")') finish - start
      End If
      ! Corrections for periodic blocks
      If (periodic) Then
        Call cpu_time(start)
        Call diagonal_correction(a, na)
        Call cpu_time(finish)
        If (ltiming) Write (*, &
          '("diagonal_correction    time=",f10.3," s")') finish - start
        If (.Not. diagonal) Then
          Call cpu_time(start)
          Call offdiagonal_correction(a, na)
          Call cpu_time(finish)
          If (ltiming) Write (*, &
            '("offdiagonal_correction time=",f10.3," s")') finish - start
        End If
      End If
      ! Free memory
      ! deallocate(aii,aij,aji,x,y,tmp,c,d)
      ! if (periodic) deallocate(delta,borders)
      ! ======================================================================
    Else
      ! PARDISO SOLVER
      ! ======================================================================
      mtype = 3   ! complex structurally symmetric
      maxfct = 1   ! number of sparse factors in memory
      mnum = 1   ! matrix to factorize
      phase = 13   ! analysis, numerical factorization, solve,
      ! iterative refinement
      msglvl = 0   ! write messages if 1
      ! Total number of non-zero elements (expected)
      ntot = sum(bdims(1:nb)**2) + sum(bdims(1:nb-1)**2) + &
        sum(bdims(2:nb)**2)
      If (periodic) ntot = ntot + bdims(1)**2 + bdims(nb)**2
      Call cpu_time(start)
      Allocate (asparse(ntot), b(n,n), tmp(n,n), ia(n+1), ja(ntot), &
        perm(n))
      Call sparse_fill(n, a, na, asparse, b, nb, bdims, ntot, ia, ja, &
        periodic)
      Call pardisoinit(pt, mtype, iparm)
      Call pardiso(pt, maxfct, mnum, mtype, phase, n, asparse, ia, ja, &
        perm, n, iparm, msglvl, b, tmp, info)
      If (info/=0) Stop 'fail in PARDISO'
      a(1:n, 1:n) = tmp
      Deallocate (asparse, b, tmp, ia, ja, perm)
      Call cpu_time(finish)
      If (ltiming) Write (*, '("PARDISO sparse solver  time=",f10.3," s")' &
        ) finish - start
      ! ======================================================================
    End If
    ! ======================================================================
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine sparse_fill(k, a, na, asparse, b, nb, bdims, ntot, ia, ja, &
    periodic)
    ! Copies the matrix A to storage in compressed sparse row format
    ! Fills in matrix B with unit matrix

    Implicit None

    ! Number of rows/columns of A, number of non-zero blocks
    Integer, Intent (In) :: k
    ! The size of the matrix A in storage
    Integer, Intent (In) :: na
    ! The target matrix
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! Total number of non-zero elements
    Integer, Intent (In) :: ntot
    ! The non-zero blocks stored in compressed sparse row format
    Complex (Kind=dp), Intent (Out) :: asparse(ntot)
    ! Unit matrix
    Complex (Kind=dp), Intent (Out) :: b(k, k)
    ! Number of blocks per row/column
    Integer, Intent (In) :: nb
    ! Size of each block
    Integer, Intent (In) :: bdims(nb)
    ! Non-zero rows
    Integer, Intent (Out) :: ia(k+1)
    ! Non-zero columns in each row
    Integer, Intent (Out) :: ja(ntot)
    ! Periodic blocks
    Logical, Intent (In) :: periodic
    ! ----------------------------------------------------------------------
    Integer :: i, ib, i0, j, jb, j0, m, n, itot

    itot = 0
    ! Rows
    ! write(*,*) "Debug: ib, i, jb, j, ia, ja, itot"
    Do ib = 1, nb
      i0 = sum(bdims(1:ib-1))
      n = bdims(ib)
      Do i = 1, n
        ia(i0+i) = itot + 1
        ! ------------------------------------------------------------------
        ! Columns
        ! A(N,1)
        If (periodic .And. ib==nb) Then
          jb = 1
          j0 = sum(bdims(1:jb-1))
          m = bdims(jb)
          Do j = 1, m
            ! compressed sparse row storage
            itot = itot + 1
            asparse(itot) = a(i0+i, j0+j)
            ja(itot) = j0 + j
            ! write(*,'(10i8)') ib, i, jb, j, ia(i0+i), ja(itot), itot
          End Do
        End If
        ! A(i,j)
        Do jb = max(1, ib-1), min(nb, ib+1)
          j0 = sum(bdims(1:jb-1))
          m = bdims(jb)
          Do j = 1, m
            ! compressed sparse row storage
            itot = itot + 1
            asparse(itot) = a(i0+i, j0+j)
            ja(itot) = j0 + j
            ! write(*,'(10i8)') ib, i, jb, j, ia(i0+i), ja(itot), itot
          End Do
        End Do
        ! A(1,N)
        If (periodic .And. ib==1) Then
          jb = nb
          j0 = sum(bdims(1:jb-1))
          m = bdims(jb)
          Do j = 1, m
            ! compressed sparse row storage
            itot = itot + 1
            asparse(itot) = a(i0+i, j0+j)
            ja(itot) = j0 + j
          ! write(*,'(10i8)') ib, i, jb, j, ia(i0+i), ja(itot), itot
          End Do
        End If
        ! ------------------------------------------------------------------
      End Do
    End Do
    ia(k+1) = itot + 1
    ! Goodbye A
    a = czero
    ! rhs for inversion
    b = czero
    Do i = 1, k
      b(i, i) = cone
    End Do
    ! write(*,*) itot
    ! write(*,'(1000i6)') ia(1:k+1)
    ! write(*,'(1000i6)') ja(1:ntot)
    ! write(*,'(1000es10.3)') asparse(1:ntot)
    ! write(*,'(1000es10.3)') b(1:k,1:k)
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine tridiag_save(a, na, periodic)
    ! Copy tridiagonals of matrix A

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! Corner blocks?
    Logical, Intent (In) :: periodic
    ! -------------------------------
    Integer :: i, i0, j0, m, n


    ! ------------------------------------------
    Do i = 1, nb - 1
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
    End Do
    ! Last diagonal block
    i0 = sum(bdims(1:nb-1))
    n = bdims(nb)
    aii(1:n, 1:n, nb) = a(i0+1:i0+n, i0+1:i0+n)
    ! ------------------------------------------
    If (periodic) Then
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
      Do i = 1, n
        aii(i, i, 1) = aii(i, i, 1) + cminus
      End Do
      ! A(N,N) = A(N,N) - A(N,1) A(1,N)
      Call zgemm('N', 'N', m, m, n, cminus, aji(:,:,nb), mb, aij(:,:,nb), &
        mb, cone, aii(:,:,nb), mb)
    End If
    ! Goodbye input matrix
    a = czero
    ! All done
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine get_diagonal(a, na, periodic)
    ! Fills in the diagonal of the matrix inverse
    ! Also computes the first/last row/column for periodic blocks

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! Corner blocks?
    Logical, Intent (In) :: periodic
    ! -------------------------------
    Integer :: i, i0, n, m


    ! ----------------------------------------------------------------------
    ! Recursion for Y matrices
    ! Y(1) = 0
    n = bdims(1)
    y(1:n, 1:n, 1) = czero
    Do i = 2, nb
      m = bdims(i-1)
      n = bdims(i)
      ! Y(i) = (A(i-1,i-1) - Y(i-1))^-1
      y(1:m, 1:m, i) = aii(1:m, 1:m, i-1) - y(1:m, 1:m, i-1)
      ! ***** Inverse *****
      Call zgetrf(m, m, y(:,:,i), mb, ipiv, info)
      If (info/=0) Stop 'get_diagonal: LU fail in Y'
      Call zgetri(m, y(:,:,i), mb, ipiv, work, lwork, info)
      If (info/=0) Stop 'get_diagonal: inverse fail in Y'
      ! ******************
      ! tmp = Y(i) * A(i-1,i)  -- superdiagonal
      Call zgemm('N', 'N', m, n, m, cone, y(:,:,i), mb, aij(:,:,i-1), mb, &
        czero, tmp, mb)
      ! Y(i) = A(i,i-1) * tmp  -- subdiagonal
      Call zgemm('N', 'N', n, n, m, cone, aji(:,:,i-1), mb, tmp, mb, &
        czero, y(:,:,i), mb)
      ! call out_mat(y(1:n,1:n,i),n,n,i,0,'Y(i)  ')
    End Do
    ! ----------------------------------------------------------------------
    ! Recursion for X matrices
    ! X(n) = 0
    n = bdims(nb)
    x(1:n, 1:n, nb) = czero
    Do i = nb - 1, 1, -1
      m = bdims(i+1)
      n = bdims(i)
      ! X(i) = (A(i+1,i+1) - X(i+1))^-1
      x(1:m, 1:m, i) = aii(1:m, 1:m, i+1) - x(1:m, 1:m, i+1)
      ! ***** Inverse *****
      Call zgetrf(m, m, x(:,:,i), mb, ipiv, info)
      If (info/=0) Stop 'get_diagonal: LU fail in X'
      Call zgetri(m, x(:,:,i), mb, ipiv, work, lwork, info)
      If (info/=0) Stop 'get_diagonal: inverse fail in X'
      ! ******************
      ! tmp = X(i) * A(i+1,i)  -- subdiagonal
      Call zgemm('N', 'N', m, n, m, cone, x(:,:,i), mb, aji(:,:,i), mb, &
        czero, tmp, mb)
      ! X(i) = A(i,i+1) * tmp  -- superdiagonal
      Call zgemm('N', 'N', n, n, m, cone, aij(:,:,i), mb, tmp, mb, czero, &
        x(:,:,i), mb)
    ! call out_mat(x(1:n,1:n,i),n,n,i,0,'X(i)  ')
    End Do
    ! ----------------------------------------------------------------------
    ! Diagonal elements of inverse
    Do i = 1, nb
      i0 = sum(bdims(1:i-1))
      n = bdims(i)
      ! B(i,i) = (A(i,i) - X(i) - Y(i))^-1
      tmp(1:n, 1:n) = aii(1:n, 1:n, i) - x(1:n, 1:n, i) - y(1:n, 1:n, i)
      ! ***** Inverse *****
      Call zgetrf(n, n, tmp, mb, ipiv, info)
      If (info/=0) Stop 'get_diagonal: LU fail in B(i,i)'
      Call zgetri(n, tmp, mb, ipiv, work, lwork, info)
      If (info/=0) Stop 'get_diagonal: inverse fail in B(i,i)'
      ! *******************
      ! call out_mat(tmp(1:n,1:n),n,n,i,i,'B(i,i)')
      a(i0+1:i0+n, i0+1:i0+n) = tmp(1:n, 1:n)
    End Do
    ! ----------------------------------------------------------------------
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine get_offdiagonal(a, na, periodic)
    ! Fills in the off-diagonal elements of the matrix inverse
    ! It includes the correction for periodic blocks

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! Corner blocks?
    Logical, Intent (In) :: periodic
    ! -------------------------------
    Logical, Parameter :: columns = .True.
    Integer :: i, j, istop, jstop


    ! Off-diagonal elements of inverse without corner blocks
    If (columns) Then
      jstop = 0
      ! First and last columns were already calculated
      If (periodic) jstop = 1
      ! Build one column at a time
      Do j = 1 + jstop, nb - jstop
        Call build_col(a, na, j, jstop)
      End Do
    Else
      istop = 0
      ! First and last rows were already calculated
      If (periodic) istop = 1
      ! Build one row at a time
      Do i = 1 + istop, nb - istop
        Call build_row(a, na, i, istop)
      End Do
    End If
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine diagonal_correction(a, na)
    ! Corrects the diagonal blocks for the periodic case

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! -------------------------------
    Integer :: i, i0, n, k


    ! ----------------------------------------------------------------------
    ! Save borders of inverse matrix with block prefactors
    Call save_borders(a, na)
    ! -----------------------------------------------------------------------
    ! Update the diagonals of the inverse matrix with the correction
    ! Using rescaled border rows and columns (see save_borders)
    k = bdims(1)
    Do i = 1, nb
      i0 = sum(bdims(1:i-1))
      n = bdims(i)
      ! tmp = B(i,i)
      tmp(1:n, 1:n) = a(i0+1:i0+n, i0+1:i0+n)
      ! tmp = tmp - B(i,1) * B(1,i)
      Call zgemm('N', 'N', n, n, k, cminus, borders(:,:,1,i), mb, &
        borders(:,:,3,i), mb, cone, tmp, mb)
      ! tmp = tmp - B(i,1) * B(N,i)
      Call zgemm('N', 'N', n, n, k, cminus, borders(:,:,1,i), mb, &
        borders(:,:,4,i), mb, cone, tmp, mb)
      ! tmp = tmp - B(i,N) * B(1,i)
      Call zgemm('N', 'N', n, n, k, cminus, borders(:,:,2,i), mb, &
        borders(:,:,3,i), mb, cone, tmp, mb)
      ! tmp = tmp - B(i,N) * B(N,i)
      Call zgemm('N', 'N', n, n, k, cminus, borders(:,:,2,i), mb, &
        borders(:,:,4,i), mb, cone, tmp, mb)
      ! B(i,i) = tmp
      a(i0+1:i0+n, i0+1:i0+n) = tmp(1:n, 1:n)
    End Do
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine offdiagonal_correction(a, na)
    ! Corrects the off-diagonal blocks for the periodic case

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! -------------------------------
    Integer :: i, i0, n, j, j0, m, k


    ! Update the off-diagonals of the inverse matrix with the correction
    ! Using rescaled border rows and columns (see save_borders)
    k = bdims(1)
    Do j = 1, nb
      j0 = sum(bdims(1:j-1))
      m = bdims(j)
      Do i = 1, nb
        If (i==j) Cycle   ! diagonals were taken care of already
        i0 = sum(bdims(1:i-1))
        n = bdims(i)
        ! tmp = B(i,j)
        tmp(1:n, 1:m) = a(i0+1:i0+n, j0+1:j0+m)
        ! tmp = tmp - B(i,1) * B(1,j)
        Call zgemm('N', 'N', n, m, k, cminus, borders(:,:,1,i), mb, &
          borders(:,:,3,j), mb, cone, tmp, mb)
        ! tmp = tmp - B(i,1) * B(N,j)
        Call zgemm('N', 'N', n, m, k, cminus, borders(:,:,1,i), mb, &
          borders(:,:,4,j), mb, cone, tmp, mb)
        ! tmp = tmp - B(i,N) * B(1,j)
        Call zgemm('N', 'N', n, m, k, cminus, borders(:,:,2,i), mb, &
          borders(:,:,3,j), mb, cone, tmp, mb)
        ! tmp = tmp - B(i,N) * B(N,j)
        Call zgemm('N', 'N', n, m, k, cminus, borders(:,:,2,i), mb, &
          borders(:,:,4,j), mb, cone, tmp, mb)
        ! B(i,j) = tmp
        a(i0+1:i0+n, j0+1:j0+m) = tmp(1:n, 1:m)
      End Do
    End Do
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine save_borders(a, na)
    ! Saves the borders of the inverse matrix without corner blocks
    ! These are multiplied by appropriate blocks for later use
    ! The column blocks have dimensions bdims(i) * bdims(1)
    ! The row blocks have dimensions bdims(1) * bdims(j)

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    ! Should have the diagonals of the inverse already stored
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! -------------------------------
    Integer :: i, j, i0, j0, m, n, k


    ! ----------------------------------------------------------------------
    ! Generate the borders of the matrix inverse
    Call build_col(a, na, 1, 0)
    Call build_col(a, na, nb, 0)
    Call build_row(a, na, 1, 0)
    Call build_row(a, na, nb, 0)
    ! ----------------------------------------------------------------------
    ! Central part of the Sherman-Morrison-Woodbury formula
    i0 = 0
    n = bdims(1)
    j0 = sum(bdims(1:nb-1))
    m = bdims(nb)
    ! Delta = A(1,N) * B(N,N) * A(N,1)
    delta(1:m, 1:m) = a(j0+1:j0+m, j0+1:j0+m)
    Call zgemm('N', 'N', m, n, m, cone, delta, mb, aji(:,:,nb), mb, czero, &
      tmp, mb)
    Call zgemm('N', 'N', n, n, m, cone, aij(:,:,nb), mb, tmp, mb, czero, &
      delta, mb)
    ! Delta = Delta + A(1,N) * B(N,1)
    tmp = a(j0+1:j0+m, i0+1:i0+n)
    Call zgemm('N', 'N', n, n, m, cone, aij(:,:,nb), mb, tmp, mb, cone, &
      delta, mb)
    ! Delta = Delta + B(1,N) * A(N,1)
    tmp = a(i0+1:i0+n, j0+1:j0+m)
    Call zgemm('N', 'N', n, n, m, cone, tmp, mb, aji(:,:,nb), mb, cone, &
      delta, mb)
    ! Delta = Delta + B(1,1) + 1
    delta(1:n, 1:n) = delta(1:n, 1:n) + a(i0+1:i0+n, i0+1:i0+n)
    Do i = 1, n
      delta(i, i) = delta(i, i) + cone
    End Do
    ! ***** Inverse *****
    Call zgetrf(n, n, delta, mb, ipiv, info)
    If (info/=0) Stop 'save_borders: LU fail in Delta'
    Call zgetri(n, delta, mb, ipiv, work, lwork, info)
    If (info/=0) Stop 'save_borders: inverse fail in Delta'
    ! *******************
    ! call out_mat(delta(1:n,1:n),n,n,1,1,'Delta ')
    ! ----------------------------------------------------------------------
    ! First column: B(i,1)
    j0 = 0
    m = bdims(1)
    Do i = 1, nb
      i0 = sum(bdims(1:i-1))
      n = bdims(i)
      borders(1:n, 1:m, 1, i) = a(i0+1:i0+n, j0+1:j0+m)   ! B(i,1)
      ! call out_mat(borders(1:n,1:m,1,i),n,m,i,1,'B(i,1)')
    End Do
    ! ----------------------------------------------------------------------
    ! Last column: B(i,N) * A(N,1)
    j0 = sum(bdims(1:nb-1))
    m = bdims(nb)
    k = bdims(1)
    Do i = 1, nb
      i0 = sum(bdims(1:i-1))
      n = bdims(i)
      tmp(1:n, 1:m) = a(i0+1:i0+n, j0+1:j0+m)   ! B(i,N)
      Call zgemm('N', 'N', n, k, m, cone, tmp, mb, aji(:,:,nb), mb, czero, &
        borders(:,:,2,i), mb)
      ! call out_mat(borders(1:n,1:k,2,i),n,k,i,nb,'B(i,N)')
    End Do
    ! ----------------------------------------------------------------------
    ! First row: Delta * B(1,j)
    i0 = 0
    n = bdims(1)
    Do j = 1, nb
      j0 = sum(bdims(1:j-1))
      m = bdims(j)
      tmp(1:n, 1:m) = a(i0+1:i0+n, j0+1:j0+m)   ! B(1,j)
      Call zgemm('N', 'N', n, m, n, cone, delta, mb, tmp, mb, czero, &
        borders(:,:,3,j), mb)
      ! call out_mat(borders(1:n,1:m,3,j),n,m,1,j,'B(1,j)')
    End Do
    ! ----------------------------------------------------------------------
    ! Last row: Delta * A(1,N) * B(N,j)
    i0 = sum(bdims(1:nb-1))
    n = bdims(nb)
    k = bdims(1)
    Do j = 1, nb
      j0 = sum(bdims(1:j-1))
      m = bdims(j)
      borders(1:n, 1:m, 4, j) = a(i0+1:i0+n, j0+1:j0+m)   ! B(N,j)
      Call zgemm('N', 'N', k, m, n, cone, aij(:,:,nb), mb, &
        borders(:,:,4,j), mb, czero, tmp, mb)
      Call zgemm('N', 'N', k, m, k, cone, delta, mb, tmp, mb, czero, &
        borders(:,:,4,j), mb)
      ! call out_mat(borders(1:k,1:m,4,j),k,m,j,nb,'B(N,j)')
    End Do
    ! ----------------------------------------------------------------------
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine build_col(a, na, j, jstop)
    ! Builds the j-th column of the inverse starting from the diagonal

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! Desired column, where to stop
    Integer, Intent (In) :: j, jstop
    ! -------------------------------
    Integer :: i, i0, n, j0, m, k


    ! ----------------------------------------------------------------------
    ! Subdiagonal blocks (i > j)
    ! B(i,j) = - C(i) * A(i,i-1) * B(i-1,j)
    ! ----------------------------------------------------------------------
    j0 = sum(bdims(1:j-1))
    m = bdims(j)
    k = m
    ! B(j,j)  -->  saved in array D
    d(1:m, 1:m) = a(j0+1:j0+m, j0+1:j0+m)
    Do i = j + 1, nb - jstop
      i0 = sum(bdims(1:i-1))
      n = bdims(i)
      ! tmp = - A(i,i-1) * B(i-1,j)  --> subdiagonal
      Call zgemm('N', 'N', n, m, k, cminus, aji(:,:,i-1), mb, d, mb, &
        czero, tmp, mb)
      ! C(i) = (A(i,i) - X(i))^-1
      c(1:n, 1:n) = aii(1:n, 1:n, i) - x(1:n, 1:n, i)
      ! ***** Inverse *****
      Call zgetrf(n, n, c, mb, ipiv, info)
      If (info/=0) Stop 'build_col: LU fail in C(i)'
      Call zgetri(n, c, mb, ipiv, work, lwork, info)
      If (info/=0) Stop 'build_col: inverse fail in C(i)'
      ! *******************
      ! B(i,j) = C(i) * tmp  -->  D used as temporary storage
      Call zgemm('N', 'N', n, m, n, cone, c, mb, tmp, mb, czero, d, mb)
      ! call out_mat(d(1:n,1:m),n,m,i,j,'B(i,j)')
      a(i0+1:i0+n, j0+1:j0+m) = d(1:n, 1:m)
      k = n   ! save old row dimension
    End Do
    ! ----------------------------------------------------------------------
    ! Superdiagonal blocks (i < j)
    ! B(i,j) = - D(i) * A(i,i+1) * B(i+1,j)
    ! ----------------------------------------------------------------------
    k = m
    ! B(j,j)  -->  saved in array C
    c(1:m, 1:m) = a(j0+1:j0+m, j0+1:j0+m)
    Do i = j - 1, 1 + jstop, -1
      i0 = sum(bdims(1:i-1))
      n = bdims(i)
      ! tmp = - A(i,i+1) * B(i+1,j)  -->  superdiagonal
      Call zgemm('N', 'N', n, m, k, cminus, aij(:,:,i), mb, c, mb, czero, &
        tmp, mb)
      ! D(i) = (A(i,i) - Y(i))^-1
      d(1:n, 1:n) = aii(1:n, 1:n, i) - y(1:n, 1:n, i)
      ! ***** Inverse *****
      Call zgetrf(n, n, d, mb, ipiv, info)
      If (info/=0) Stop 'build_col: LU fail in D(i)'
      Call zgetri(n, d, mb, ipiv, work, lwork, info)
      If (info/=0) Stop 'build_col: inverse fail in D(i)'
      ! *******************
      ! B(i,j) = D(i) * tmp  -->  C used as temporary storage
      Call zgemm('N', 'N', n, m, n, cone, d, mb, tmp, mb, czero, c, mb)
      ! call out_mat(c(1:n,1:m),n,m,i,j,'B(i,j)')
      a(i0+1:i0+n, j0+1:j0+m) = c(1:n, 1:m)
      k = n   ! save old row dimension
    End Do
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine build_row(a, na, i, istop)
    ! Builds the i-th row of the inverse starting from the diagonal

    Implicit None

    ! Size of square matrix a in memory
    Integer, Intent (In) :: na
    ! Block tridiagonal matrix a to be inverted
    Complex (Kind=dp), Intent (Inout) :: a(na, na)
    ! Desired column, where to stop
    Integer, Intent (In) :: i, istop
    ! -------------------------------
    Integer :: i0, n, j, j0, m, k


    ! ----------------------------------------------------------------------
    ! Subdiagonal blocks (i > j)
    ! B(i,j) = - B(i,j+1) * A(j+1,j) * D(j)
    ! ----------------------------------------------------------------------
    i0 = sum(bdims(1:i-1))
    n = bdims(i)
    k = n
    ! B(i,i)  -->  saved in array C
    c(1:n, 1:n) = a(i0+1:i0+n, i0+1:i0+n)
    Do j = i - 1, 1 + istop, -1
      j0 = sum(bdims(1:j-1))
      m = bdims(j)
      ! tmp = - B(i,j+1) A(j+1,j)  -->  subdiagonal
      Call zgemm('N', 'N', n, m, k, cminus, c, mb, aji(:,:,j), mb, czero, &
        tmp, mb)
      ! D(i) = (A(i,i) - Y(i))^-1
      d(1:m, 1:m) = aii(1:m, 1:m, j) - y(1:m, 1:m, j)
      ! ***** Inverse *****
      Call zgetrf(m, m, d, mb, ipiv, info)
      If (info/=0) Stop 'build_row: LU fail in D(i)'
      Call zgetri(m, d, mb, ipiv, work, lwork, info)
      If (info/=0) Stop 'build_row: inverse fail in D(i)'
      ! *******************
      ! B(i,j) = tmp * D(j)
      Call zgemm('N', 'N', n, m, m, cone, tmp, mb, d, mb, czero, c, mb)
      ! call out_mat(c(1:n,1:m),n,m,i,j,'B(i,j)')
      a(i0+1:i0+n, j0+1:j0+m) = c(1:n, 1:m)
      k = m   ! save old column dimension
    End Do
    ! ----------------------------------------------------------------------
    ! Superdiagonal blocks (i < j)
    ! B(i,j) = - B(i,j-1) * A(j-1,j) * C(j)
    ! ----------------------------------------------------------------------
    k = n
    ! B(i,i)  -->  saved in array D
    d(1:n, 1:n) = a(i0+1:i0+n, i0+1:i0+n)
    Do j = i + 1, nb - istop
      j0 = sum(bdims(1:j-1))
      m = bdims(j)
      ! tmp = - B(i,j-1) * A(j-1,j)  -->  superdiagonal
      Call zgemm('N', 'N', n, m, k, cminus, d, mb, aij(:,:,j-1), mb, &
        czero, tmp, mb)
      ! C(j) = (A(j,j) - X(j))^-1
      c(1:m, 1:m) = aii(1:m, 1:m, j) - x(1:m, 1:m, j)
      ! ***** Inverse *****
      Call zgetrf(m, m, c, mb, ipiv, info)
      If (info/=0) Stop 'build_row: LU fail in C(i)'
      Call zgetri(m, c, mb, ipiv, work, lwork, info)
      If (info/=0) Stop 'build_row: inverse fail in C(i)'
      ! *******************
      ! B(i,j) = tmp * C(j) -->  D used as temporary storage
      Call zgemm('N', 'N', n, m, m, cone, tmp, mb, c, mb, czero, d, mb)
      ! call out_mat(d(1:n,1:m),n,m,i,j,'B(i,j)')
      a(i0+1:i0+n, j0+1:j0+m) = d(1:n, 1:m)
      k = m   ! save old column dimension
    End Do
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------


  ! -----------------------------------------------------------------------
  Subroutine out_mat(a, n, m, ib, jb, label)

    Implicit None

    ! Number of rows and columns
    Integer, Intent (In) :: n, m
    ! Matrix to output
    Complex (Kind=dp), Intent (In) :: a(n, m)
    ! Info: which block, text label
    Integer, Intent (In) :: ib, jb
    Character (Len=6), Intent (In) :: label
    ! ------------------------------------
    Integer :: i

    Write (*, *) label, ib, jb
    Do i = 1, n
      Write (*, '(1000es10.3)') a(i, 1:m)
    End Do
    ! All done!
  End Subroutine
  ! -----------------------------------------------------------------------

#ifdef CPP_MPI
  Subroutine bcast_params_godfrin(t_godfrin)   ! GODFRIN Flaviano
    ! broadcast parameters for godfrin matrix inversion scheme
    use :: mpi
    Implicit None

    Type (type_godfrin), Intent (Inout) :: t_godfrin
    Integer :: ierr

    Call mpi_bcast(t_godfrin%na, 1, mpi_integer, 0, mpi_comm_world, ierr)
    If (ierr/=mpi_success) Stop &
      '[bcast_params_savewf] Error broadcasting maxmem_number'
    Call mpi_bcast(t_godfrin%nb, 1, mpi_integer, 0, mpi_comm_world, ierr)
    If (ierr/=mpi_success) Stop &
      '[bcast_params_savewf] Error broadcasting maxmem_units'
    Call mpi_bcast(t_godfrin%ldiag, 1, mpi_logical, 0, mpi_comm_world, ierr)
    If (ierr/=mpi_success) Stop &
      '[bcast_params_savewf] Error broadcasting save_rll'
    Call mpi_bcast(t_godfrin%lper, 1, mpi_logical, 0, mpi_comm_world, ierr)
    If (ierr/=mpi_success) Stop &
      '[bcast_params_savewf] Error broadcasting save_sll'
    Call mpi_bcast(t_godfrin%lpardiso, 1, mpi_logical, 0, mpi_comm_world, ierr)
    If (ierr/=mpi_success) Stop &
      '[bcast_params_savewf] Error broadcasting save_rllleft'

    If (.Not. allocated(t_godfrin%bdims)) Allocate (t_godfrin%bdims( &
      t_godfrin%nb))
    Call mpi_bcast(t_godfrin%bdims, t_godfrin%nb, mpi_integer, 0, &
      mpi_comm_world, ierr)
    If (ierr/=mpi_success) Stop &
      '[bcast_params_savewf] Error broadcasting save_sllleft'

  End Subroutine
#endif

End Module
