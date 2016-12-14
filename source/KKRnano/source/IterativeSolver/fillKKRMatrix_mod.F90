!> Build coefficient matrix for solution of Dyson equation.
!>
!> @author Elias Rabel

module fillKKRMatrix_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: getKKRMatrixStructure, getKKRSolutionStructure, getRightHandSideStructure
  public :: buildKKRCoeffMatrix, buildRightHandSide
  public :: getGreenDiag
  public :: dump
  public :: convertBSRToFullMatrix, convertToFullMatrix, convertFullMatrixToBSR

  double complex, parameter :: ZERO=(0.d0, 0.d0), CONE=(1.d0, 0.d0)
  
  interface dump
    module procedure dumpDenseMatrix, dumpSparseMatrixData
  endinterface

!   interface convert
!     module procedure convertBSRToFullMatrix, convertToFullMatrix, convertFullMatrixToBSR, convertFullMatrixTo
!   endinterface

  contains

  !------------------------------------------------------------------------------
  !> Setup of the sparsity pattern of the KKR-Solution.
  subroutine getKKRSolutionStructure(bsr_X, lmax_a)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, create
    type(SparseMatrixDescription), intent(out) :: bsr_X
    integer(kind=1), intent(in) :: lmax_a(:,:) !< lmax of each interaction dim(naez_trc,num_local_atoms), -1: truncated

    integer, parameter :: GROUPING = 0 ! -1:never, 0:auto, 1:always
    integer :: iRow, jCol, Xind, nnzb, group
 
    nnzb = count(lmax_a >= 0) ! number of non-zero blocks in X

    if (nnzb < size(lmax_a)) then ! truncation can be applied
      group = max(0, GROUPING)
      if (GROUPING > 0) warn(6, "X_mat will be stored rectangular!") ! always group
    else
      group = 1 + GROUPING
      if (GROUPING < 0) warn(6, "X_mat will NOT be made rectangular, although it could!") ! never group
    endif

    if (group > 0) then

      call create(bsr_X, nRows=size(lmax_a, 1), nnzb=size(lmax_a, 1), nCols=1)

      ! generate BSR descriptor of a flat structure, fuse RHS atoms into rectangular blocks
      bsr_X%RowStart(:) = [(iRow, iRow = 1, bsr_X%nRows + 1)] ! start indices into ColIndex
      bsr_X%ColIndex(:) = 1

      return
    endif

    call create(bsr_X, nRows=size(lmax_a, 1), nnzb=nnzb, nCols=size(lmax_a, 2))

    assert( size(bsr_X%RowStart) == 1 + bsr_X%nRows )
    assert( size(bsr_X%ColIndex) == bsr_X%nnzb )

    Xind = 0
    do iRow = 1, bsr_X%nRows
      bsr_X%RowStart(iRow) = Xind + 1 ! start indices into ColIndex
      do jCol = 1, bsr_X%nCols

        if (lmax_a(iRow,jCol) >= 0) then ! sub-optimal indexing into lmax_a here
          Xind = Xind + 1 ! create a new entry
          bsr_X%ColIndex(Xind) = jCol
        endif ! block is non-zero

      enddo ! jCol
! #ifndef NDEBUG
!       write(0, '(i6,a,999i1)') iRow,"  ",lmax_a(iRow,:)+1 ! beware: the iRow index is in internal representation
! #endif
    enddo ! iRow
    bsr_X%RowStart(bsr_X%nRows + 1) = Xind + 1 ! final, important since the ranges are always [RowStart(i) ... RowStart(i+1)-1]
    assert( Xind == bsr_X%nnzb ) ! check

  endsubroutine ! getKKRSolutionStructure
  
  
  !------------------------------------------------------------------------------
  !> Setup of the sparsity pattern of the KKR-Matrix.
  subroutine getKKRMatrixStructure(bsr_A, numn0, indn0)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, create
    type(SparseMatrixDescription), intent(out) :: bsr_A
    integer, intent(in) :: numn0(:) !< dim(nRows) number of non-zero elements in each row
    integer(kind=2), intent(in) :: indn0(:,:) !< dim(maxval(numn0),nRows) column indices of non-zero elements per row

    integer :: Aind, iRow, iacls, jCol

    call create(bsr_A, nRows=size(numn0), nnzb=sum(numn0)) ! nCols==nRows by default

    assert( size(indn0, 2) >= bsr_A%nRows )
    assert( size(bsr_A%RowStart) == 1 + bsr_A%nRows )
    assert( size(bsr_A%ColIndex) == bsr_A%nnzb )

    Aind = 0
    do iRow = 1, bsr_A%nRows
      bsr_A%RowStart(iRow) = Aind + 1 ! start indices into ColIndex
      do iacls = 1, numn0(iRow)
        assert( Aind < bsr_A%nnzb )

        assert( iacls <= size(indn0, 1) )
        jCol = indn0(iacls,iRow)
        assert( 1 <= jCol .and. jCol <= bsr_A%nCols )

        Aind = Aind + 1 ! create a new entry
        bsr_A%ColIndex(Aind) = jCol

      enddo ! iacls
    enddo ! iRow
    bsr_A%RowStart(bsr_A%nRows + 1) = Aind + 1 ! final, important since the ranges are always [RowStart(i) ... RowStart(i+1)-1]
    assert( Aind == bsr_A%nnzb ) ! check

  endsubroutine ! getKKRMatrixStructure

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine getRightHandSideStructure(bsr_B, row_indices, nRHS_group)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, create
    type(SparseMatrixDescription), intent(out) :: bsr_B
    integer(kind=2), intent(in) :: row_indices(:) !> truncation zone indices of the local atoms
    integer, intent(in) :: nRHS_group ! 1:no grouping, 2+:grouping of RHS columns

    integer :: iRow, jCol, row_index, ist, Bind, nRows
    integer, allocatable :: ColIndices(:) ! temporary structure

    nRows = maxval(row_indices)
    call create(bsr_B, nRows, nnzb=size(row_indices), nCols=size(row_indices))

    allocate(ColIndices(bsr_B%nRows))
    ColIndices = -1
    
    do jCol = 1, bsr_B%nCols ! loop over right hand sides
      row_index = row_indices(jCol)
      ColIndices(row_index) = jCol
    enddo ! jCol == iRHS
    
    Bind = 0
    do iRow = 1, bsr_B%nRows
      bsr_B%RowStart(iRow) = Bind + 1 ! start indices into ColIndex
      jCol = ColIndices(iRow)
      if (jCol > 0) then
  
        Bind = Bind + 1 ! create a new entry
        bsr_B%ColIndex(Bind) = (jCol - 1)/nRHS_group + 1 ! when grouping is active, these need to be reduced by the group size

      endif ! in this row one non-zero entry exists
    enddo ! iRow
    bsr_B%RowStart(bsr_B%nRows + 1) = Bind + 1 ! final, important since the ranges are always [RowStart(i) ... RowStart(i+1)-1]
    assert( Bind == bsr_B%nnzb ) ! check
    
    deallocate(ColIndices, stat=ist) ! ignore status   
  endsubroutine ! getRightHandSideStructure

  
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
        ! multiply the T matrix onto smat

#ifdef  TRANSPOSE_TO_ROW_MAJOR
#define smat(a,b,c) SMAT(b,a,c)
#endif

        do lm2 = 1, lmsd

          temp(:) = ZERO
          do lm1 = 1, lmsd
            temp(:) = temp(:) + tmatLL(:,lm1,iRow) * smat(lm1,lm2,Aind) ! T*G
          enddo ! lm1

          if (iRow == jCol) temp(lm2) = temp(lm2) - CONE ! subtract 1.0 from the diagonal

          smat(:,lm2,Aind) = temp(:)
        enddo ! lm2

#ifdef  TRANSPOSE_TO_ROW_MAJOR
#undef  smat
#endif

      enddo ! Aind ! block columns
    enddo ! block rows

  endsubroutine ! buildKKRCoeffMatrix

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine buildRightHandSide(mat_B, bsr_B, row_indices, tmatLL)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, exists
  
    double complex, intent(out) :: mat_B(:,:,:) !> dim(lmsa,lmsd,nRHSs)
    type(SparseMatrixDescription), intent(in) :: bsr_B
    integer(kind=2), intent(in) :: row_indices(:) !> truncation zone indices of the local atoms
    double complex, intent(in), optional :: tmatLL(:,:,:) !< dim(lmsd,lmsd,nRows)

    integer :: iRHS, row_index, lms, lmsd, Bind, nRHS_group, iRHS_group, ing

    mat_B = ZERO
    lmsd  = size(tmatLL, 2)
    assert( size(tmatLL, 1) == lmsd ) ! assume square shape
    assert( size(mat_B, 1)  >= lmsd )
!   assert( size(mat_B, 2)  == lmsd ) ! <future>
#ifdef  TRANSPOSE_TO_ROW_MAJOR
#define mat_B(a,b,c) MAT_B(b,a,c)
    assert( modulo(size(mat_B, 1), lmsd) == 0 ) ! the number of columns is a multiple of the block dimension
    nRHS_group = size(mat_B, 1) / lmsd ! integer division should be exact
#else
    assert( modulo(size(mat_B, 2), lmsd) == 0 ) ! the number of columns is a multiple of the block dimension
    nRHS_group = size(mat_B, 2) / lmsd ! integer division should be exact
#endif
    assert( nRHS_group >= 1 )

    do iRHS = 1, size(row_indices)
      row_index = row_indices(iRHS)
      iRHS_group = (iRHS - 1)/nRHS_group  + 1
      ing = modulo((iRHS - 1),nRHS_group) + 1
      Bind = exists(bsr_B, row=row_index, col=iRHS_group) ! find index in BSR structure
      if (Bind < 1) stop "fatal! bsr_B cannot find diagonal element"

      if (present(tmatLL)) then
        mat_B(:lmsd,lmsd*(ing - 1) + 1:lmsd*ing,Bind) = tmatLL(:,:,row_index)
      else
        do lms = 1, lmsd
          mat_B(lms,lms + lmsd*(ing - 1),Bind) = CONE ! set the block to a unity matrix
        enddo ! lms
      endif

#ifdef  TRANSPOSE_TO_ROW_MAJOR
#undef  mat_B
#endif

    enddo ! iRHS
    
  endsubroutine ! buildRightHandSide


  
  !------------------------------------------------------------------------------
  !> Copy the diagonal elements G_{LL'}^{nn'} of the Green's-function,
  !> dependent on (k,E) into matrix G_diag
  subroutine getGreenDiag(G_diag, mat_X, bsr_X, atom_index, iRHS)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, exists
    double complex, intent(out) :: G_diag(:,:) ! dim(lmsd,lmsd)
    double complex, intent(in) :: mat_X(:,:,:) !> dim(lmsa,lmsd*nRHS_group,nnzb)
    type(SparseMatrixDescription), intent(in) :: bsr_X
    integer(kind=2), intent(in) :: atom_index
    integer, intent(in) :: iRHS

    integer :: lmsd, Xind, row, nRHS_group, iRHS_group, ing
    !                                      nn
    !         Copy the diagonal elements G_LL' of the Green's-function,
    !         dependent on (k,E) into matrix G_diag
    !         (n = n' = atom_index)
    
    row = atom_index ! convert from integer(kind=2) to default integer

    lmsd = size(G_diag, 2)
    assert( lmsd == size(G_diag, 1) )
!   assert( lmsd == size(mat_X, 1) ) ! fails of alignment padding is non-zero
    
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    assert( modulo(size(mat_X, 1), lmsd) == 0 ) ! the number of columns is a multiple of the block dimension
    nRHS_group = size(mat_X, 1) / lmsd ! integer division should be exact
#else
    assert( modulo(size(mat_X, 2), lmsd) == 0 ) ! the number of columns is a multiple of the block dimension
    nRHS_group = size(mat_X, 2) / lmsd ! integer division should be exact
#endif
    assert( nRHS_group >= 1 )
    
      iRHS_group = (iRHS - 1)/nRHS_group  + 1
      ing = modulo((iRHS - 1),nRHS_group) + 1
      Xind = exists(bsr_X, row, col=iRHS_group) ! find index in BSR structure
      if (Xind < 1) die_here("diagonal element not contained in X, row="-row-", col="-iRHS)

      G_diag = ZERO
#ifdef  TRANSPOSE_TO_ROW_MAJOR
      G_diag(:,:) = transpose(mat_X(:lmsd,(ing - 1)*lmsd + 1:lmsd*ing,Xind))
#else
      G_diag(:,:) = mat_X(:lmsd,(ing - 1)*lmsd + 1:lmsd*ing,Xind)
#endif
  endsubroutine ! getGreenDiag

  
  subroutine convertBSRToFullMatrix(full, bsr, smat)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
  
    double complex, intent(out) :: full(:,:)
    type(SparseMatrixDescription), intent(in) :: bsr
    double complex, intent(in) :: smat(:,:,:) !< dim(lmsa,lmsd,nnzb)
    
    call convertToFullMatrix(full, bsr%RowStart, bsr%ColIndex, smat)
    
  endsubroutine
  
  !----------------------------------------------------------------------------
  !> Given the sparse matrix data 'smat' and the sparsity information,
  !> create the dense matrix representation of the matrix.
  subroutine convertToFullMatrix(full, RowStart, ColIndex, smat)
  
    double complex, intent(out) :: full(:,:) !< dim(lmsa*nRows,lmsd*nCols)
    integer, intent(in) :: RowStart(:) !> dim(nRows + 1)
    integer, intent(in) :: ColIndex(:) !> dim(nnzb)
    double complex, intent(in) :: smat(:,:,:) !< dim(lmsa,lmsd,nnzb)

    integer :: iRow, jCol, Aind, lmsd, lmsa

    full = ZERO

    lmsa = size(smat, 1)
    lmsd = size(smat, 2)
!   assert( lmsd == size(smat, 1) ) ! does not hold if we use this method also to convert X_mat with rectangular (non-square) blocks that originate from grouping
    assert( size(ColIndex) == size(smat, 3) )

    do iRow = 1, size(RowStart) - 1
      do Aind = RowStart(iRow), RowStart(iRow+1) - 1
        jCol = ColIndex(Aind)

        full(lmsa*(iRow - 1) + 1:lmsa*iRow,lmsd*(jCol - 1) + 1:lmsd*jCol) = smat(1:lmsa,:,Aind)

      enddo ! Aind
    enddo ! iRow

  endsubroutine ! convertToFullMatrix
  
  
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

    integer :: iRow, jCol, Aind, lmsd, lmsa

    smat = ZERO

    lmsa = size(smat, 1)
    lmsd = size(smat, 2)
    assert( size(ColIndex) == size(smat, 3) )

    do iRow = 1, size(RowStart) - 1
      do Aind = RowStart(iRow), RowStart(iRow+1) - 1
        jCol = ColIndex(Aind)

        smat(1:lmsa,:,Aind) = full(lmsa*(iRow - 1) + 1:lmsa*iRow,lmsd*(jCol - 1) + 1:lmsd*jCol)

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
