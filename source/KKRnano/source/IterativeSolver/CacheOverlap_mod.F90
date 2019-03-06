!> For the multiplication of two BSR-matrices (block sparse row matrices)
!> we want to permute the order of rows such that a high data-reuse is given

module CacheOverlap_mod
  implicit none
  private
  public :: maximize_overlap

  contains

  subroutine maximize_overlap(perm, ia, ja, weights) ! scales as nRows^2
    integer, intent(out) :: perm(:) ! dim(nRows) new permutation
    integer, intent(in) :: ia(:) ! dim(nRows+1)
    integer, intent(in) :: ja(:) ! dim(nnzb)
    integer, intent(in) :: weights(:) ! dim(nCols), nCols==nRows
     
    integer :: nRows, nCols, iRow, jRow, start_pair(2), inxt(1), ist, Aind, ijovl, ovlsum, sums(0:1)
    integer, allocatable :: ovl(:,:), rowi(:), iperm(:)
    
    nRows = size(ia) - 1
    if (size(ja) /= ia(nRows + 1) - 1) stop __LINE__
    nCols = nRows
    if (size(weights) /= nCols) stop __LINE__

    if (size(perm) /= nRows) stop __LINE__
    
    allocate(ovl(nRows,nRows)) ! shared
    allocate(rowi(nCols)) ! private
    
    ovlsum = 0
    do iRow = 1, nRows
      rowi(:) = 0
      do Aind = ia(iRow), ia(iRow + 1) - 1
        rowi(ja(Aind)) = weights(ja(Aind))
      enddo ! Aind

      ovl(iRow,iRow) = 0 ! diagonal elements are meaningless
      do jRow = 1, iRow - 1 ! loop over off-diagonal elements

        ijovl = 0
        do Aind = ia(jRow), ia(jRow + 1) - 1
          ijovl = ijovl + rowi(ja(Aind))
        enddo ! Aind
        ovl(jRow,iRow) = ijovl ! weighted overlap
        ovl(iRow,jRow) = ijovl ! symmetric

        if (jRow == iRow - 1) ovlsum = ovlsum + ijovl
      enddo ! jRow
    enddo ! iRow
    sums(0) = ovlsum
    
    deallocate(rowi, stat=ist) ! ignore status

#ifdef DEBUG
    do iRow = 1, nRows
      write(*, '(a,999I4)') "overlap: ",ovl(:,iRow) ! DEBUG
    enddo ! iRow
#endif

    allocate(iperm(nRows))
    iperm = 0 ! init
    
    start_pair = maxloc(ovl)
    jRow = start_pair(2) ! start where the overlap is largest
    
    ovlsum = 0
    do iRow = 1, nRows
      perm(iRow) = jRow
      iperm(jRow) = iRow ! inverse permutation
      inxt = maxloc(ovl(:,jRow), mask=(iperm == 0))
      ovlsum = ovlsum + ovl(inxt(1),jRow)
      ovl(:,jRow) = 0
      ovl(jRow,:) = 0
      jRow = inxt(1)
    enddo ! iRow
    sums(1) = ovlsum

    if (any(ovl /= 0)) stop 'DEBUG: maximize_overlap: not all elements were eleminated!'
    
    write(*, '(9(a,f0.3))') "improve overlap from ",sums(0)/dble(nRows)," to ",sums(1)/dble(nRows)
#ifdef DEBUG
    write(*, '(a,999(" ",i0))') "overlap permutation     ", perm  ! DEBUG
    write(*, '(a,999(" ",i0))') "overlap inv. permutation", iperm ! DEBUG
#endif
    deallocate(ovl, iperm, stat=ist) ! ignore status
  endsubroutine ! max

endmodule ! CacheOverlap_mod
