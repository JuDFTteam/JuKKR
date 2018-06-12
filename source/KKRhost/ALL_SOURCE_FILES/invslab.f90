! ************************************************************************
subroutine invslab(gdi, gup, gdow, gin, icheck)
  use :: mod_datatypes, only: dp
  ! ************************************************************************
  ! ************************************************************************

  ! ---> ALGORITM FOR SLAB GEOMETRY

  ! ------------------------------------------------------------------------

  ! ---> factorization D ^-1 = (prod L) * M * (prod U)

  ! see notes R. Zeller

  ! ------------------------------------------------------------------------
  implicit none

  ! .. parameters ..
  include 'inc.p'
  ! *  NPOTD = 2 * NATYPD                                               *
  ! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
  ! *  NSPIND = 1                                                       *
  ! *                                                                   *
  ! *********************************************************************

  ! .. local scalars ..
  ! .. data statements ..

  ! .. external subroutines ..
  integer, parameter :: lmmaxd = (krel+korbit+1)*(lmaxd+1)**2
  integer, parameter :: ndim = nprincd*lmmaxd
  integer, parameter :: almd = naezd*lmmaxd
  complex (kind=dp), parameter :: ci = (0.e0_dp, 1.e0_dp)
  complex (kind=dp), parameter :: czero = (0.e0_dp, 0.e0_dp)
  complex (kind=dp), parameter :: cone = (1.e0_dp, 0.e0_dp)
  integer :: ipvt(ndim)
  complex (kind=dp) :: gup(ndim, ndim, nlayerd), gdow(ndim, ndim, nlayerd), &
    gdi(ndim, ndim, nlayerd), dmat(ndim, ndim, nlayerd), &
    dinver(ndim, ndim, nlayerd), f(ndim, ndim), e(ndim, ndim), g(ndim, ndim), &
    cunit(ndim, ndim), gin(almd, almd), gdiold(ndim, ndim, nlayerd)

  integer :: n, lm, info, i, j, irow
  ! ---> to be changed: set up the triangular matrix.
  integer :: icheck(nlayerd, nlayerd)
  ! first version:

  external :: cinit, zcopy, zgemm, zgetrf, zgetrs, btom





  call cinit(ndim*ndim, e)
  call cinit(ndim*ndim, f)
  call cinit(ndim*ndim, g)
  call cinit(ndim*ndim, cunit)
  call cinit(ndim*ndim*nlayerd, dmat)
  ! ---> calculate D_1 = (M_11)**(-1)
  do n = 1, ndim
    cunit(n, n) = cone
  end do

  do n = 1, nlayerd
    call zcopy(ndim*ndim, gdi(1,1,n), 1, gdiold(1,1,n), 1)
  end do



  ! ---> claculate D_N (2 <= N <= NLAYERD)

  call zcopy(ndim*ndim, gdi(1,1,1), 1, e(1,1), 1)
  call zcopy(ndim*ndim, cunit, 1, dmat(1,1,1), 1)

  call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
  call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, dmat(1,1,1), ndim, info)
  ! ---> F = D(N-1) * M1(N-1)


  do n = 2, nlayerd


    ! ---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)

    call zgemm('N', 'N', ndim, ndim, ndim, cone, dmat(1,1,n-1), ndim, &
      gup(1,1,n-1), ndim, czero, f(1,1), ndim)




    call zcopy(ndim*ndim, gdi(1,1,n), 1, e(1,1), 1)
    call zcopy(ndim*ndim, cunit(1,1), 1, dmat(1,1,n), 1)

    call zgemm('N', 'N', ndim, ndim, ndim, -cone, gdow(1,1,n-1), ndim, f(1,1), &
      ndim, cone, e(1,1), ndim)
    ! At this point the matrix DMAT(ndim,ndim,nlayerd) contains the
    call zcopy(ndim*ndim, e(1,1), 1, dinver(1,1,n), 1)
    ! matrices [of dimension (ndim,ndim)]  D^n, n=1,..,nlayerd
    call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
    call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, dmat(1,1,n), ndim, info)

  end do
  ! ---> calculate Z_n for 1 =< n <= n-1





  do n = nlayerd, 1, (-1)

    if (n==nlayerd) then

      call zcopy(ndim*ndim, dmat(1,1,nlayerd), 1, e(1,1), 1)

      call btom(nlayerd, nlayerd, e, ndim, gin, almd, .false.)

      call zcopy(ndim*ndim, dmat(1,1,nlayerd), 1, gdi(1,1,nlayerd), 1)
    else

      call zgemm('N', 'N', ndim, ndim, ndim, cone, gdow(1,1,n), ndim, &
        dmat(1,1,n), ndim, czero, f(1,1), ndim)

      call bofm(n+1, n+1, e, ndim, gin, almd)

      call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, f(1,1), ndim, &
        czero, g(1,1), ndim)
      call zgemm('N', 'N', ndim, ndim, ndim, cone, gup(1,1,n), ndim, g(1,1), &
        ndim, czero, f(1,1), ndim)

      do lm = 1, ndim
        f(lm, lm) = cone + f(lm, lm)
      end do

      call zgemm('N', 'N', ndim, ndim, ndim, cone, dmat(1,1,n), ndim, f(1,1), &
        ndim, czero, e(1,1), ndim)
      ! here start the two loops on the row index,
      call btom(n, n, e, ndim, gin, almd, .false.)
      ! in order to span all the matrix
      call zcopy(ndim*ndim, e(1,1), 1, gdi(1,1,n), 1)
      ! and to calculate just the blocks that are needed
    end if
    ! for the construction of the cluster of green's function



    ! this is the loop for element G_ij with i<j

    if (icheck(n,n)==0) go to 110

    if (n==1) go to 100



    do irow = (n-1), 1, (-1)

      if (icheck(irow,n)==1) then

        call bofm(irow+1, n, e, ndim, gin, almd)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, gup(1,1,irow), ndim, &
          e(1,1), ndim, czero, f(1,1), ndim)

        call zgemm('N', 'N', ndim, ndim, ndim, -cone, dmat(1,1,irow), ndim, &
          f(1,1), ndim, czero, e(1,1), ndim)

        call btom(irow, n, e, ndim, gin, almd, .false.)

      end if

      ! this is the loop for element G_ij with i>j
    end do

100 continue

    if (n==nlayerd) go to 110



    do irow = n + 1, nlayerd, 1

      if (icheck(irow,n)==1) then

        call zcopy(ndim*ndim, cunit(1,1), 1, e(1,1), 1)

        call bofm(irow, irow, f, ndim, gin, almd)

        call zgetrf(ndim, ndim, f(1,1), ndim, ipvt, info)
        call zgetrs('N', ndim, ndim, f(1,1), ndim, ipvt, e(1,1), ndim, info)

        do i = 1, ndim
          do j = 1, ndim
            f(i, j) = gdiold(i, j, irow) - (dinver(i,j,irow)-e(i,j))
          end do
        end do

        call zcopy(ndim*ndim, cunit(1,1), 1, e(1,1), 1)

        call zgetrf(ndim, ndim, f(1,1), ndim, ipvt, info)
        call zgetrs('N', ndim, ndim, f(1,1), ndim, ipvt, e(1,1), ndim, info)

        call zgemm('N', 'N', ndim, ndim, ndim, -cone, e(1,1), ndim, &
          gdow(1,1,irow-1), ndim, czero, f(1,1), ndim)
        ! corrected 15.3.2000
        call bofm(irow-1, n, e, ndim, gin, almd)


        call zgemm('N', 'N', ndim, ndim, ndim, cone, f(1,1), ndim, e(1,1), &
          ndim, czero, g(1,1), ndim)



        call btom(irow, n, g, ndim, gin, almd, .false.)


      end if


    end do

110 continue

  end do


  ! ************************************************************************
  ! ************************************************************************
  ! ************************************************************************
  return

end subroutine invslab
