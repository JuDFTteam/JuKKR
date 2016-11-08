!> Build coefficient matrix for solution of Dyson equation.
!>
!> @author Elias Rabel

module fillKKRMatrix_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: getKKRMatrixStructure, buildKKRCoeffMatrix, buildRightHandSide, solveFull
  public :: dump
  public :: getKKRSolutionStructure
  public :: convertBSRToFullMatrix, convertToFullMatrix, convertFullMatrixToBSR

  double complex, parameter :: ZERO=(0.d0, 0.d0), CONE=(1.d0, 0.d0)
  
  interface dump
    module procedure dumpDenseMatrix, dumpSparseMatrixData
  endinterface

  contains

  !------------------------------------------------------------------------------
  !> Setup of the sparsity pattern of the KKR-Solution.
  subroutine getKKRSolutionStructure(lmax_a, bsr_X)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, create

    integer(kind=1), intent(in) :: lmax_a(:,:) !< lmax of each interaction dim(naez_trc,num_local_atoms)
    type(SparseMatrixDescription), intent(inout) :: bsr_X

    integer :: iRow, jCol, Xind

    call create(bsr_X, nRows=size(lmax_a, 1), nnzb=count(lmax_a >= 0), nCols=size(lmax_a, 2))

    assert( size(bsr_X%RowStart) == 1 + bsr_X%nRows )
    assert( size(bsr_X%ColIndex) == bsr_X%nnzb )
    
    Xind = 0
    do iRow = 1, bsr_X%nRows
      bsr_X%RowStart(iRow) = Xind + 1 ! start indices into ColIndex
      do jCol = 1, bsr_X%nCols

        if (lmax_a(iRow,jCol) >= 0) then ! sub-optimal indexing here
          Xind = Xind + 1
          bsr_X%ColIndex(Xind) = jCol
        endif ! block is non-zero

      enddo ! jCol
    enddo ! iRow
    bsr_X%RowStart(bsr_X%nRows + 1) = Xind + 1 ! final, important since the ranges are always [RowStart(i) ... RowStart(i+1)-1]
    assert( Xind == bsr_X%nnzb ) ! check

  endsubroutine ! getKKRSolutionStructure
  
  
  !------------------------------------------------------------------------------
  !> Setup of the sparsity pattern of the KKR-Matrix.
  subroutine getKKRMatrixStructure(numn0, indn0, bsr_A)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, create

    integer, intent(in) :: numn0(:) !< dim(nRows) number of non-zero elements in each row
    integer(kind=2), intent(in) :: indn0(:,:) !< dim(maxval(numn0),nRows) column indices of non-zero elements per row
    type(SparseMatrixDescription), intent(out) :: bsr_A

    integer :: Aind, iRow, icol, jCol

    call create(bsr_A, nRows=size(numn0), nnzb=sum(numn0)) ! nCols==nRows by default

    assert( size(indn0, 2) >= bsr_A%nRows )
    assert( size(bsr_A%RowStart) == 1 + bsr_A%nRows )
    assert( size(bsr_A%ColIndex) == bsr_A%nnzb )

    Aind = 0
    do iRow = 1, bsr_A%nRows
      bsr_A%RowStart(iRow) = Aind + 1 ! start indices into ColIndex
      do icol = 1, numn0(iRow)
        assert( Aind < bsr_A%nnzb )

        assert( icol <= size(indn0, 1) )
        jCol = indn0(icol,iRow)
        assert( 1 <= jCol .and. jCol <= bsr_A%nCols )

        Aind = Aind + 1
        bsr_A%ColIndex(Aind) = jCol

      enddo ! icol
    enddo ! iRow
    bsr_A%RowStart(bsr_A%nRows + 1) = Aind + 1 ! final, important since the ranges are always [RowStart(i) ... RowStart(i+1)-1]
    assert( Aind == bsr_A%nnzb ) ! check

  endsubroutine ! getKKRMatrixStructure

  !---------------------------------------------------------------------
  !> Given G_ref we build (T*G_ref - 1) -> coefficent matrix.
  !>
  !> @param smat  on input:  sparse block matrix containing Gref
  !>              on output: coefficient matrix (1 - TG_ref)
  !> @param sparse  sparse matrix description
  !> RowStart    for each row give index of first non-zero block in ColIndex
  !> ColIndex    column index array of non-zero blocks
  subroutine buildKKRCoeffMatrix(smat, tmatLL, bsr_A)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription

    double complex, intent(inout) :: smat(:,:,:) !< dim(lmsa,lmsd,bsr_A%nnzb)
    double complex, intent(in) :: tmatLL(:,:,:) !< dim(lmsd,lmsd,bsr_A%nRows)
    type(SparseMatrixDescription), intent(in) :: bsr_A

    double complex :: temp(size(smat, 1))
    integer :: lmsa, lmsd, lm2, lm1
    integer :: iRow, jCol, Aind

    lmsa = size(smat, 1) ! may be memory-aligned
    lmsd = size(tmatLL, 2)
    assert( size(tmatLL, 1) == lmsd )
    assert( lmsa >= lmsd )
    assert( bsr_A%nRows == size(tmatLL, 3) )

    do iRow = 1, bsr_A%nRows
      do Aind = bsr_A%RowStart(iRow), bsr_A%RowStart(iRow + 1) - 1
        jCol = bsr_A%ColIndex(Aind) ! ColIndex gives the block-column indices of non-zero blocks

#ifndef NDEBUG
        if (1 > iRow .or. iRow > bsr_A%nRows) then
          write (*,*) "buildKKRCoeffMatrix: invalid iRow", iRow
          stop
        endif

        if (1 > jCol .or. jCol > bsr_A%nCols) then
          write (*,*) "buildKKRCoeffMatrix: invalid jCol", jCol
          stop
        endif
#endif
        ! multiply the T matrix onto 

        do lm2 = 1, lmsd

          temp(:) = ZERO
          do lm1 = 1, lmsd
            temp(:) = temp(:) + tmatLL(:,lm1,iRow) * smat(lm1,lm2,Aind) ! T*G
          enddo ! lm1

          if (iRow == jCol) temp(lm2) = temp(lm2) - CONE ! subtract 1.0 from the diagonal

          smat(:,lm2,Aind) = temp(:)
        enddo ! lm2

      enddo ! Aind ! block columns
    enddo ! block rows

  endsubroutine ! buildKKRCoeffMatrix

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine buildRightHandSide(mat_B, bsr_B, atom_indices, tmatLL)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, exists
  
    double complex, intent(out) :: mat_B(:,:,:) !> dim(lmsa,lmsd,nRHSs)
    type(SparseMatrixDescription), intent(in) :: bsr_B
    integer(kind=2), intent(in) :: atom_indices(:) !> truncation zone indices of the local atoms
    double complex, intent(in), optional :: tmatLL(:,:,:) !< dim(lmsd,lmsd,nRows)

    integer :: iRHS, atom_index, lm2, lmsd

    mat_B = ZERO
    lmsd  = size(tmatLL, 2)
    assert( size(tmatLL, 1) == lmsd ) ! assume square shape
    assert( size(mat_B, 1)  >= lmsd )
!   assert( size(mat_B, 2)  == lmsd ) ! <future>
    assert( modulo(size(mat_B, 2), lmsd) == 0 ) ! the number of columns is a multiple of the block dimension

    do iRHS = 1, size(atom_indices)
      atom_index = atom_indices(iRHS)
#ifdef useBSR
      atom_index = exists(bsr_B, row=atom_index, col=iRHS)
      if (atom_index < 1) stop "fatal! bsr_B cannot find diagonal element"
      
      if (present(tmatLL)) then
        mat_B(:lmsd,:,atom_index) = tmatLL(:,:,atom_indices(iRHS))
      else
        do lm2 = 1, lmsd
          mat_B(lm2,lm2,atom_index) = CONE ! set the block to a unity matrix
        enddo ! lm2
      endif
#else
      if (present(tmatLL)) then
        mat_B(:lsmd,lmsd*(iRHS - 1) + 1:lmsd*iRHS,atom_index) = tmatLL(:,:,atom_index)
      else
        do lm2 = 1, lmsd
          mat_B(lm2,lm2 + lmsd*(iRHS - 1),atom_index) = CONE ! set the block to a unity matrix
        enddo ! lm2
      endif
#endif
    enddo ! iRHS

  endsubroutine ! buildRightHandSide


  
  subroutine convertBSRToFullMatrix(smat, bsr, full)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
  
    double complex, intent(in) :: smat(:,:,:) !< dim(lmsa,lmsd,nnzb)
    type(SparseMatrixDescription), intent(in) :: bsr
    double complex, intent(out) :: full(:,:)
    
    call convertToFullMatrix(smat, bsr%RowStart, bsr%ColIndex, full)
    
  endsubroutine
  
  !----------------------------------------------------------------------------
  !> Given the sparse matrix data 'smat' and the sparsity information,
  !> create the dense matrix representation of the matrix.
  subroutine convertToFullMatrix(smat, RowStart, ColIndex, full)
  
    double complex, intent(in) :: smat(:,:,:) !< dim(lmsa,lmsd,nnzb)
    integer, intent(in) :: RowStart(:) !> dim(nRows + 1)
    integer, intent(in) :: ColIndex(:) !> dim(nnzb)
    double complex, intent(out) :: full(:,:)

    integer :: iRow, jCol, Aind, lmsd

    full = ZERO

    lmsd = size(smat, 2)
    assert( size(smat, 1) == lmsd )
    assert( size(ColIndex) == size(smat, 3) )

    do iRow = 1, size(RowStart) - 1
      do Aind = RowStart(iRow), RowStart(iRow+1) - 1
        jCol = ColIndex(Aind)

        full(lmsd*(iRow - 1) + 1:lmsd*iRow,lmsd*(jCol - 1) + 1:lmsd*jCol) = smat(1:lmsd,:,Aind)

      enddo ! Aind
    enddo ! iRow

  endsubroutine ! convertToFullMatrix
  
  
  !----------------------------------------------------------------------------
  !> Solution of a system of linear equations with multiple right hand sides,
  !> using standard dense matrix LAPACK routines.
  integer function solveFull(full_A, full_X) result(info)
    double complex, intent(inout) :: full_A(:,:)
    double complex, intent(inout) :: full_X(:,:) ! on entry this contains full_B, on exit the solution

    integer, allocatable :: ipvt(:)
    integer :: ndim, nRHSs
    external :: zgetrf, zgetrs ! LAPACK

    ndim  = size(full_A, 1)
    nRHSs = size(full_X, 2)
    assert( size(full_A, 2) == ndim ) ! must be square
    assert( size(full_X, 1) == ndim ) ! must match the dims of A
    
    allocate(ipvt(ndim))

    ! factorization
    call zgetrf(ndim, ndim, full_A, ndim, ipvt, info) ! LU-factorize
!   if (info /= 0) ! ToDo: warn

    call zgetrs('n', ndim, nRHSs, full_A, ndim, ipvt, full_X, ndim, info) ! solve the system of linear equations
!   if (info /= 0) ! ToDo: warn

    deallocate(ipvt, stat=info)
  endfunction ! solveFull
  
  
  subroutine convertFullMatrixToBSR(smat, bsr, full)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
  
    double complex, intent(out) :: smat(:,:,:) !< dim(lmsa,lmsd,nnzb)
    type(SparseMatrixDescription), intent(in) :: bsr
    double complex, intent(in) :: full(:,:)
    
    call convertFullMatrixTo(smat, bsr%RowStart, bsr%ColIndex, full)
    
  endsubroutine
  
  !----------------------------------------------------------------------------
  !> Given the sparse matrix data 'smat' and the sparsity information,
  !> create the dense matrix representation of the matrix.
  subroutine convertFullMatrixTo(smat, RowStart, ColIndex, full)
  
    double complex, intent(out) :: smat(:,:,:) !< dim(lmsa,lmsd,nnzb)
    integer, intent(in) :: RowStart(:) !> dim(nRows + 1)
    integer, intent(in) :: ColIndex(:) !> dim(nnzb)
    double complex, intent(in) :: full(:,:)

    integer :: iRow, jCol, Aind, lmsd

    smat = ZERO

    lmsd = size(smat, 2)
    assert( size(smat, 1) == lmsd )
    assert( size(ColIndex) == size(smat, 3) )

    do iRow = 1, size(RowStart) - 1
      do Aind = RowStart(iRow), RowStart(iRow+1) - 1
        jCol = ColIndex(Aind)

        smat(1:lmsd,:,Aind) = full(lmsd*(iRow - 1) + 1:lmsd*iRow,lmsd*(jCol - 1) + 1:lmsd*jCol)

      enddo ! Aind
    enddo ! iRow

  endsubroutine ! convertToFullMatrix
  
  
  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to a formatted or an unformatted file
  subroutine dumpSparseMatrixData(smat, filename, formatted)
    double complex, intent(in) :: smat(:,:,:)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: formatted

    if (formatted) then
      call dumpSparseMatrixDataFormatted(smat, filename)
    else
      call dumpSparseMatrixDataBinary(smat, filename)
    endif
  endsubroutine ! dump
  
  !----------------------------------------------------------------------------
  !> Write dense matrix to unformatted file
  subroutine dumpDenseMatrix(mat, filename, formatted)
    double complex, intent(in) :: mat(:,:)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: formatted

    if (formatted) then
      call dumpDenseMatrixFormatted(mat, filename)
    else
      call dumpDenseMatrixBinary(mat, filename)
    endif
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to unformatted file
  subroutine dumpSparseMatrixDataBinary(smat, filename)
    double complex, intent(in) :: smat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) smat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read sparse matrix data from unformatted file
  subroutine loadSparseMatrixData(smat, filename)
    double complex, intent(out) :: smat(:,:,:)
    character(len=*), intent(in) :: filename
    
    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) smat
    close(fu)
  endsubroutine ! load

  !----------------------------------------------------------------------------
  !> Write dense matrix to unformatted file
  subroutine dumpDenseMatrixBinary(mat, filename)
    double complex, intent(in) :: mat(:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) mat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read dense matrix data from unformatted file
  !> - useful for testing.
  subroutine loadDenseMatrix(mat, filename)
    double complex, intent(out) :: mat(:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) mat
    close(fu)
  endsubroutine ! load

  !----------------------------------------------------------------------------
  !> Write GLLh to unformatted file
  !> - useful for testing.
  subroutine dumpGLLhBinary(mat, filename)
    double complex, intent(in) :: mat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) mat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read GLLh from unformatted file
  !> - useful for testing.
  subroutine loadGLLh(mat, filename)
    double complex, intent(out) :: mat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) mat
    close(fu)
  endsubroutine ! load

  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to formatted file
  !> - useful for testing.
  subroutine dumpSparseMatrixDataFormatted(smat, filename)
    double complex, intent(in) :: smat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97
    integer :: fi, si, bi

    open(fu, file=filename, form='formatted', action='write')
    do bi = 1, size(smat, 3) ! block index
      do si = 1, size(smat, 2) !  slow index
        do fi = 1, size(smat, 1) !  fast index
          write(fu, *) real(smat(fi,si,bi)), aimag(smat(fi,si,bi))
        enddo ! fi
      enddo ! si
    enddo ! bi
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Write dense matrix to formatted file
  !> First line gives matrix dimension: rows cols
  !> - useful for testing.
  subroutine dumpDenseMatrixFormatted(mat, filename)
    double complex, intent(in) :: mat(:,:)
    character(len=*), intent(in) :: filename
    
    integer, parameter :: fu = 97
    integer :: ii, jj

    open(fu, file=filename, form='formatted', action='write')
    write(fu, *) size(mat, 1), size(mat, 2)
    do jj = 1, size(mat, 2)
      do ii = 1, size(mat, 1)
         write(fu, *) real(mat(ii,jj)), aimag(mat(ii,jj))
      enddo ! ii
    enddo ! jj
    close(fu)
  endsubroutine ! dump

endmodule ! fillKKRMatrix_mod
