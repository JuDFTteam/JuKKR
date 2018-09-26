module mod_surfgf

contains

  ! -------------------------------------------------------------------------------
  ! SUBROUTINE: SURFGF
  !> @brief Solve surface green's function: \f$
  ! f(x)=ml\left(m0-x\right)^{\left(-1\right)*mr} \f$
  !> @details method: decimation technique
  !> input:  ml,m0,mr - complex rectangular matrices
  !> output: x        - result, matrix of same type as before
  !> @note NEW VERSION (speeded up) by V.Bellini (march,1999)
  !>
  !> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to
  ! Fortran90
  ! -------------------------------------------------------------------------------
  subroutine surfgf(ndim, ml, m0, mr, x, itermax, errmax, ichck)

    use :: mod_constants
    use :: mod_profiling
    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_cinit

    implicit none

    ! -------------------------------------------------------------------------------
    ! For KREL = 1 (relativistic mode)

    ! NPOTD = 2 * NATYPD
    ! LMMAXD = 2 * (LMAXD+1)^2
    ! NSPIND = 1
    ! LMGF0D = (LMAXD+1)^2 dimension of the reference system Green
    ! function, set up in the spin-independent non-relativstic
    ! (l,m_l)-representation

    ! -------------------------------------------------------------------------------
    ! .. Input variables
    integer, intent (in) :: ndim
    integer, intent (in) :: ichck
    integer, intent (in) :: itermax
    real (kind=dp), intent (in) :: errmax
    ! .. Input arrays
    complex (kind=dp), dimension (ndim, ndim), intent (in) :: m0
    complex (kind=dp), dimension (ndim, ndim), intent (in) :: ml
    complex (kind=dp), dimension (ndim, ndim), intent (in) :: mr
    ! .. Output variables
    complex (kind=dp), dimension (ndim, ndim), intent (out) :: x
    ! .. Local Scalars
    integer :: i, info, iter, j, n
    real (kind=dp) :: err, sum, xim, xre
    ! .. Local Arrays
    integer, dimension (ndim) :: ipvt
    complex (kind=dp), dimension (ndim, ndim) :: aa, alfa, bb, beta, cc, cunit, eps, tempin
    complex (kind=dp), dimension (ndim, ndim) :: tempout, y1, y2
    ! ..

    call cinit(ndim*ndim, cunit)
    do n = 1, ndim
      cunit(n, n) = cone
    end do

    call zcopy(ndim*ndim, m0, 1, eps, 1)
    call zcopy(ndim*ndim, ml, 1, alfa, 1)
    call zcopy(ndim*ndim, mr, 1, beta, 1)
    call zcopy(ndim*ndim, m0, 1, x, 1)

    iter = 1
100 continue

    call zcopy(ndim*ndim, eps, 1, y1, 1)
    call zcopy(ndim*ndim, y1, 1, tempin, 1)
    call zgetrf(ndim, ndim, tempin, ndim, ipvt, info)

    ! aa = eps^-1 * alfa
    call zcopy(ndim*ndim, alfa, 1, tempout, 1)
    call zgetrs('N', ndim, ndim, tempin, ndim, ipvt, tempout, ndim, info)
    call zcopy(ndim*ndim, tempout, 1, aa, 1)

    ! bb = eps^-1 * beta

    call zcopy(ndim*ndim, beta, 1, tempout, 1)
    call zgetrs('N', ndim, ndim, tempin, ndim, ipvt, tempout, ndim, info)
    call zcopy(ndim*ndim, tempout, 1, bb, 1)

    ! alfa_new = alfa * aa

    call zgemm('N', 'N', ndim, ndim, ndim, cone, alfa, ndim, aa, ndim, czero, y1, ndim)

    ! beta_new = beta * bb

    call zgemm('N', 'N', ndim, ndim, ndim, cone, beta, ndim, bb, ndim, czero, y2, ndim)

    ! cc = - alfa * bb

    call zgemm('N', 'N', ndim, ndim, ndim, -cone, alfa, ndim, bb, ndim, czero, cc, ndim)

    ! x_new = x + cc

    call zaxpy(ndim*ndim, cone, cc, 1, x, 1)

    ! cc = eps + cc

    call zaxpy(ndim*ndim, cone, cc, 1, eps, 1)

    ! eps_new = cc - beta * aa

    call zgemm('N', 'N', ndim, ndim, ndim, -cone, beta, ndim, aa, ndim, cone, eps, ndim)

    call zcopy(ndim*ndim, y1, 1, alfa, 1)
    call zcopy(ndim*ndim, y2, 1, beta, 1)

    sum = 0.e0_dp
    do i = 1, ndim
      do j = 1, ndim
        xre = real(alfa(i,j))
        xim = aimag(alfa(i,j))
        sum = sum + xre*xre + xim*xim
      end do
    end do

    err = sqrt(sum)
    if (err<errmax .or. iter>itermax) go to 110
    iter = iter + 1
    go to 100

110 continue

    call zcopy(ndim*ndim, x, 1, tempin, 1)
    call zcopy(ndim*ndim, cunit, 1, tempout, 1)
    call zgetrf(ndim, ndim, tempin, ndim, ipvt, info)
    call zgetrs('N', ndim, ndim, tempin, ndim, ipvt, tempout, ndim, info)
    call zcopy(ndim*ndim, tempout, 1, x, 1)

    call zgemm('N', 'N', ndim, ndim, ndim, cone, x, ndim, mr, ndim, czero, tempin, ndim)
    call zgemm('N', 'N', ndim, ndim, ndim, cone, ml, ndim, tempin, ndim, czero, x, ndim)

    if (iter>itermax) then
      write (6, fmt='('' itermax too small.  iter='',i3)') iter
      write (6, '('' Surfgf:  iter='',i4,''  error='',d14.7)') iter, err
    end if
    if (ichck==0) return
    ! write(6,'('' Surfgf:  iter='',i4,''  error='',d12.7)') iter,err
    ! write(6,'(/'' X matrix'')')
    ! write(6,*)

    return
  end subroutine surfgf

end module mod_surfgf
