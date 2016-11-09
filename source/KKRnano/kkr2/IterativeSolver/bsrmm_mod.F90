!> Multiplication of two BSR-matrices (block sparse row matrices)
#define DEBUG

module bsrmm_mod
  implicit none
  private
  public :: bsr_times_bsr

  
#ifdef DEBUG
  logical, private :: dbg = .true.
#define cDBG if(dbg)
#else
#define cDBG !
#endif

  contains

  subroutine bsr_times_bsr(Y, ia, ja, A, ix, jx, X, nFlops)
    double complex, intent(out) :: Y(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    integer, intent(in) :: ia(:) !> dim(A%nRows + 1) !  start indices
    integer, intent(in) :: ja(:) !> dim(A%nnzb)      ! column indices
    double complex, intent(in)  :: A(:,:,:) ! dim(lmsa,lmsd,A%nnzb)
    integer, intent(in) :: ix(:) !> dim(X%nRows + 1) !  start indices
    integer, intent(in) :: jx(:) !> dim(X%nnzb)      ! column indices
    double complex, intent(in)  :: X(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    integer(kind=8), intent(inout) :: nFlops
    ! for the result Y we assume the same BSR structure as X
#define iy ix
#define jy jx

    external :: ZGEMM ! from BLAS

    ! local variables
    integer :: XnRows, leadDim_Y, leadDim_X, leadDim_A, lmsd, nRHS
    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)
    ! private variables
    integer :: iRow, Yind, jCol, Aind, kCol, Xind
    double complex :: beta

    leadDim_A = size(A, 1)
    leadDim_X = size(X, 1)
    leadDim_Y = size(Y, 1)
    
    lmsd = size(A, 2)
    nRHS = size(X, 2)

    if (leadDim_Y /= leadDim_X) stop __LINE__
    if (leadDim_A /= leadDim_X) stop __LINE__
    if (size(X, 2) /= size(Y, 2)) stop __LINE__
    if (size(X, 3) /= size(Y, 3)) stop __LINE__

    XnRows = size(ix) - 1
#define YnRows XnRows
    if (XnRows /= size(ia) - 1) stop __LINE__

    do iRow = 1, YnRows
      ! reuse elements of Y, i.e. keep the accumulator in the cache
      do Yind = iy(iRow), iy(iRow + 1) - 1 
        jCol = jy(Yind) ! update   matrix block element Y_full[iRow,jCol]

        beta = zero
        do Aind = ia(iRow), ia(iRow + 1) - 1
          kCol = ja(Aind) ! for each matrix block element A_full[iRow,kCol]
#define kRow kCol          
          Xind = BSR_entry_exists(ix, jx, row=kRow, col=jCol) ! find out if X_full[kRow,jCol] exists
          if (Xind > 0) then ! yes

#ifdef DEBUG
#define show(X) #X,"=",X,", "
!! cDBG write(*, "(99(2a,i0,a))") show(iRow),show(jCol),show(kCol),show(Yind),show(Aind),show(Xind)
            if (Xind > size(X, 3)) then
              write(0,*) __FILE__,__LINE__," ERROR: ",show(iRow),show(jCol),show(kCol),show(Yind),show(Aind),show(Xind),&
                show(size(Y,3)),show(size(A,3)),show(size(X,3)),show(jx(ix(kRow):ix(kRow+1)-1))
            endif
            if (Xind > size(X, 3)) stop __LINE__
            if (Aind > size(A, 3)) stop __LINE__
            if (Yind > size(Y, 3)) stop __LINE__
#endif
          
            ! now: Y[:,:,Yind] += A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
            !                    M     N     K          A                       B                             C
            call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), leadDim_A, X(:,1,Xind), leadDim_X, beta, Y(:,1,Yind), leadDim_Y)

            nFlops = nFlops + (8_8*lmsd)*(nRHS*lmsd)

            beta = one
          endif ! X_full[kRow,jCol] exists
#undef kRow
        enddo ! Aind
        
      enddo ! Yind
    enddo ! iRow
cDBG iRow = show_BSR_structure(6, ia, ja, ix, jx)
cDBG dbg = .false. !! switch off after 1st iteration
  endsubroutine ! bsr_times_bsr

  integer function BSR_entry_exists(RowStart, Colindex, row, col) result(Ind)
    integer, intent(in) :: RowStart(:), ColIndex(:), row, col

    ! ToDo: the ColIndex list should be sorted ascendingly, so bisection search will be faster
    do Ind = RowStart(row), RowStart(row + 1) - 1
      if (ColIndex(Ind) == col) return ! Ind
    enddo ! Ind
    Ind = -1 ! not found

  endfunction ! exists

#ifdef DEBUG
  integer function show_BSR_structure(unit, ia, ja, ix, jx) result(ios)
    integer, intent(in) :: unit ! output unit to write to
    integer, intent(in) :: ia(:), ja(:), ix(:), jx(:)
    
    integer :: nRows, iRow, Aind, Xind, jCol
    character(len=192) :: line
    nRows = size(ia) - 1
    do iRow = 1, nRows
      line = "" ! clear line
      do Aind = ia(iRow), ia(iRow + 1) - 1 ; jCol = ja(Aind)                     ; if(jCol > len(line)) cycle
        line(jCol:jCol) = char(48 + mod(Aind, 10))
      enddo ! Aind                                                 ! add space
      do Xind = ix(iRow), ix(iRow + 1) - 1 ; jCol = jx(Xind)       + (nRows + 4) ; if(jCol > len(line)) cycle
        line(jCol:jCol) = char(48 + mod(Xind, 10))
      enddo ! Xind
      write(unit, fmt='(i6,9a)', iostat=ios) iRow,"    ",trim(line)
    enddo ! iRow
  endfunction ! show
#endif
  
endmodule ! vbrmv_mat_mod
