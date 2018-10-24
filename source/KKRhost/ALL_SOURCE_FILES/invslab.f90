module mod_invslab

  private
  public :: invslab

contains

  !-------------------------------------------------------------------------------
  !> Summary: Matrix inversion for slab geometry
  !> Author: 
  !> Category: KKRhost, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> ---> ALGORITM FOR SLAB GEOMETRY
  !>
  !> ---> factorization D ^-1 = (prod L) * M * (prod U)
  !>
  !> see notes R. Zeller
  !-------------------------------------------------------------------------------
  subroutine invslab(gdi, gup, gdow, gin, icheck)

    use :: global_variables, only: ndim_slabinv, nlayerd, alm
    use :: mod_datatypes, only: dp
    use :: mod_bofm, only: bofm
    use :: mod_btom, only: btom
    use :: mod_cinit, only: cinit
    use :: mod_constants, only: czero,cone
    implicit none

    integer :: ipvt(ndim_slabinv)
    complex (kind=dp) :: gup(ndim_slabinv, ndim_slabinv, nlayerd)
    complex (kind=dp) :: gdow(ndim_slabinv, ndim_slabinv, nlayerd)
    complex (kind=dp) :: gdi(ndim_slabinv, ndim_slabinv, nlayerd)
    complex (kind=dp) :: dmat(ndim_slabinv, ndim_slabinv, nlayerd)
    complex (kind=dp) :: dinver(ndim_slabinv, ndim_slabinv, nlayerd)
    complex (kind=dp) :: f(ndim_slabinv, ndim_slabinv)
    complex (kind=dp) :: e(ndim_slabinv, ndim_slabinv)
    complex (kind=dp) :: g(ndim_slabinv, ndim_slabinv)
    complex (kind=dp) :: cunit(ndim_slabinv, ndim_slabinv)
    complex (kind=dp) :: gin(alm, alm)
    complex (kind=dp) :: gdiold(ndim_slabinv, ndim_slabinv, nlayerd)

    integer :: n, lm, info, i, j, irow
    integer :: icheck(nlayerd, nlayerd)

    ! ---> to be changed: set up the triangular matrix.
    ! first version:

    call cinit(ndim_slabinv*ndim_slabinv, e)
    call cinit(ndim_slabinv*ndim_slabinv, f)
    call cinit(ndim_slabinv*ndim_slabinv, g)
    call cinit(ndim_slabinv*ndim_slabinv, cunit)
    call cinit(ndim_slabinv*ndim_slabinv*nlayerd, dmat)

    do n = 1, ndim_slabinv
      cunit(n, n) = cone
    end do

    do n = 1, nlayerd
      call zcopy(ndim_slabinv*ndim_slabinv, gdi(1,1,n), 1, gdiold(1,1,n), 1)
    end do

    ! ---> calculate D_1 = (M_11)**(-1)

    call zcopy(ndim_slabinv*ndim_slabinv, gdi(1,1,1), 1, e(1,1), 1)
    call zcopy(ndim_slabinv*ndim_slabinv, cunit, 1, dmat(1,1,1), 1)

    call zgetrf(ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, info)
    call zgetrs('N', ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, dmat(1,1,1), ndim_slabinv, info)

    ! ---> claculate D_N (2 <= N <= NLAYERD)
    do n = 2, nlayerd
      ! ---> F = D(N-1) * M1(N-1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, dmat(1,1,n-1), ndim_slabinv, gup(1,1,n-1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
      
      ! ---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)
      call zcopy(ndim_slabinv*ndim_slabinv, gdi(1,1,n), 1, e(1,1), 1)
      call zcopy(ndim_slabinv*ndim_slabinv, cunit(1,1), 1, dmat(1,1,n), 1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, gdow(1,1,n-1), ndim_slabinv, f(1,1), ndim_slabinv, cone, e(1,1), ndim_slabinv)
      call zcopy(ndim_slabinv*ndim_slabinv, e(1,1), 1, dinver(1,1,n), 1)
      call zgetrf(ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, info)
      call zgetrs('N', ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, dmat(1,1,n), ndim_slabinv, info)
    end do

    ! At this point the matrix DMAT(ndim_slabinv,ndim_slabinv,nlayerd) contains the
    ! matrices [of dimension (ndim_slabinv,ndim_slabinv)]  D^n, n=1,..,nlayerd

    ! ---> calculate Z_n for 1 =< n <= n-1
    do n = nlayerd, 1, (-1)
      if (n==nlayerd) then
        call zcopy(ndim_slabinv*ndim_slabinv, dmat(1,1,nlayerd), 1, e(1,1), 1)
        call btom(nlayerd, nlayerd, e, ndim_slabinv, gin, alm, .false.)
        call zcopy(ndim_slabinv*ndim_slabinv, dmat(1,1,nlayerd), 1, gdi(1,1,nlayerd), 1)
      else
        call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, gdow(1,1,n), ndim_slabinv, dmat(1,1,n), ndim_slabinv, czero, f(1,1), ndim_slabinv)
        call bofm(n+1, n+1, e, ndim_slabinv, gin, alm)
        call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, e(1,1), ndim_slabinv, f(1,1), ndim_slabinv, czero, g(1,1), ndim_slabinv)
        call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, gup(1,1,n), ndim_slabinv, g(1,1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
        do lm = 1, ndim_slabinv
          f(lm, lm) = cone + f(lm, lm)
        end do
        call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, dmat(1,1,n), ndim_slabinv, f(1,1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
        call btom(n, n, e, ndim_slabinv, gin, alm, .false.)
        call zcopy(ndim_slabinv*ndim_slabinv, e(1,1), 1, gdi(1,1,n), 1)
      end if

      ! here start the two loops on the row index,
      ! in order to span all the matrix
      ! and to calculate just the blocks that are needed
      ! for the construction of the cluster of green's function

      if (icheck(n,n)==0) go to 110

      if (n==1) go to 100

      ! this is the loop for element G_ij with i<j

      do irow = (n-1), 1, (-1)
        if (icheck(irow,n)==1) then
          call bofm(irow+1, n, e, ndim_slabinv, gin, alm)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, gup(1,1,irow), ndim_slabinv, e(1,1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, dmat(1,1,irow), ndim_slabinv, f(1,1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
          call btom(irow, n, e, ndim_slabinv, gin, alm, .false.)
        end if
      end do

100   continue

      if (n==nlayerd) go to 110

      ! this is the loop for element G_ij with i>j

      do irow = n + 1, nlayerd, 1
        if (icheck(irow,n)==1) then
          call zcopy(ndim_slabinv*ndim_slabinv, cunit(1,1), 1, e(1,1), 1)
          call bofm(irow, irow, f, ndim_slabinv, gin, alm)
          call zgetrf(ndim_slabinv, ndim_slabinv, f(1,1), ndim_slabinv, ipvt, info)
          call zgetrs('N', ndim_slabinv, ndim_slabinv, f(1,1), ndim_slabinv, ipvt, e(1,1), ndim_slabinv, info)
          do i = 1, ndim_slabinv
            do j = 1, ndim_slabinv
              f(i, j) = gdiold(i, j, irow) - (dinver(i,j,irow)-e(i,j))
            end do
          end do
          call zcopy(ndim_slabinv*ndim_slabinv, cunit(1,1), 1, e(1,1), 1)
          call zgetrf(ndim_slabinv, ndim_slabinv, f(1,1), ndim_slabinv, ipvt, info)
          call zgetrs('N', ndim_slabinv, ndim_slabinv, f(1,1), ndim_slabinv, ipvt, e(1,1), ndim_slabinv, info)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, e(1,1), ndim_slabinv, gdow(1,1,irow-1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
          call bofm(irow-1, n, e, ndim_slabinv, gin, alm)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, f(1,1), ndim_slabinv, e(1,1), ndim_slabinv, czero, g(1,1), ndim_slabinv)
 
          ! corrected 15.3.2000
          call btom(irow, n, g, ndim_slabinv, gin, alm, .false.)
 
        end if
      end do

110   continue

    end do

  end subroutine invslab

end module mod_invslab
