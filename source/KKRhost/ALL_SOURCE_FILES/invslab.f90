module mod_invslab

contains

  ! ************************************************************************
  subroutine invslab(gdi, gup, gdow, gin, icheck)
    ! ************************************************************************
    ! ************************************************************************

    ! ---> ALGORITM FOR SLAB GEOMETRY

    ! ------------------------------------------------------------------------

    ! ---> factorization D ^-1 = (prod L) * M * (prod U)

    ! see notes R. Zeller

    ! ------------------------------------------------------------------------
    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_bofm
    use :: mod_btom
    use :: mod_cinit
    implicit none

    complex (kind=dp), parameter :: czero = (0.e0_dp, 0.e0_dp)
    complex (kind=dp), parameter :: cone = (1.e0_dp, 0.e0_dp)
    integer :: ipvt(ndim_slabinv)
    complex (kind=dp) :: gup(ndim_slabinv, ndim_slabinv, nlayerd), gdow(ndim_slabinv, ndim_slabinv, nlayerd), gdi(ndim_slabinv, ndim_slabinv, nlayerd), &
      dmat(ndim_slabinv, ndim_slabinv, nlayerd), dinver(ndim_slabinv, ndim_slabinv, nlayerd), f(ndim_slabinv, ndim_slabinv), e(ndim_slabinv, ndim_slabinv), &
      g(ndim_slabinv, ndim_slabinv), cunit(ndim_slabinv, ndim_slabinv), gin(alm, alm), gdiold(ndim_slabinv, ndim_slabinv, nlayerd)

    integer :: n, lm, info, i, j, irow
    ! ---> to be changed: set up the triangular matrix.
    integer :: icheck(nlayerd, nlayerd)
    ! first version:

    call cinit(ndim_slabinv*ndim_slabinv, e)
    call cinit(ndim_slabinv*ndim_slabinv, f)
    call cinit(ndim_slabinv*ndim_slabinv, g)
    call cinit(ndim_slabinv*ndim_slabinv, cunit)
    call cinit(ndim_slabinv*ndim_slabinv*nlayerd, dmat)
    ! ---> calculate D_1 = (M_11)**(-1)
    do n = 1, ndim_slabinv
      cunit(n, n) = cone
    end do

    do n = 1, nlayerd
      call zcopy(ndim_slabinv*ndim_slabinv, gdi(1,1,n), 1, gdiold(1,1,n), 1)
    end do

    ! ---> claculate D_N (2 <= N <= NLAYERD)
    call zcopy(ndim_slabinv*ndim_slabinv, gdi(1,1,1), 1, e(1,1), 1)
    call zcopy(ndim_slabinv*ndim_slabinv, cunit, 1, dmat(1,1,1), 1)

    call zgetrf(ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, info)
    call zgetrs('N', ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, dmat(1,1,1), ndim_slabinv, info)

    ! ---> F = D(N-1) * M1(N-1)
    do n = 2, nlayerd
      ! ---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, dmat(1,1,n-1), ndim_slabinv, gup(1,1,n-1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
      call zcopy(ndim_slabinv*ndim_slabinv, gdi(1,1,n), 1, e(1,1), 1)
      call zcopy(ndim_slabinv*ndim_slabinv, cunit(1,1), 1, dmat(1,1,n), 1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, gdow(1,1,n-1), ndim_slabinv, f(1,1), ndim_slabinv, cone, e(1,1), ndim_slabinv)
      ! At this point the matrix DMAT(ndim_slabinv,ndim_slabinv,nlayerd) contains the
      call zcopy(ndim_slabinv*ndim_slabinv, e(1,1), 1, dinver(1,1,n), 1)
      ! matrices [of dimension (ndim_slabinv,ndim_slabinv)]  D^n, n=1,..,nlayerd
      call zgetrf(ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, info)
      call zgetrs('N', ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, dmat(1,1,n), ndim_slabinv, info)
    end do

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
        ! here start the two loops on the row index,
        call btom(n, n, e, ndim_slabinv, gin, alm, .false.)
        ! in order to span all the matrix
        call zcopy(ndim_slabinv*ndim_slabinv, e(1,1), 1, gdi(1,1,n), 1)
        ! and to calculate just the blocks that are needed
      end if

      ! for the construction of the cluster of green's function
      ! this is the loop for element G_ij with i<j
      if (icheck(n,n)==0) go to 110
      if (n==1) go to 100
      do irow = (n-1), 1, (-1)
        if (icheck(irow,n)==1) then
          call bofm(irow+1, n, e, ndim_slabinv, gin, alm)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, gup(1,1,irow), ndim_slabinv, e(1,1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, dmat(1,1,irow), ndim_slabinv, f(1,1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
          call btom(irow, n, e, ndim_slabinv, gin, alm, .false.)
        end if
        ! this is the loop for element G_ij with i>j
      end do

100   continue

      if (n==nlayerd) go to 110

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
          ! corrected 15.3.2000
          call bofm(irow-1, n, e, ndim_slabinv, gin, alm)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, f(1,1), ndim_slabinv, e(1,1), ndim_slabinv, czero, g(1,1), ndim_slabinv)
          call btom(irow, n, g, ndim_slabinv, gin, alm, .false.)
        end if
      end do

110   continue

    end do

  end subroutine invslab

end module mod_invslab
