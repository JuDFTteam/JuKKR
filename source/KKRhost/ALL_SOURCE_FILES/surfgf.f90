!-------------------------------------------------------------------------------
! SUBROUTINE: SURFGF
!> @brief Solve surface green's function: \f$ f(x)=ml\left(m0-x\right)^{\left(-1\right)*mr} \f$
!> @details method: decimation technique
!> input:  ml,m0,mr - complex rectangular matrices
!> output: x        - result, matrix of same type as before
!> @note NEW VERSION (speeded up) by V.Bellini (march,1999)
!>
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine surfgf(ndim, ml, m0, mr, x, itermax, errmax, ichck, lmmaxd)

      Use constants
      Use profiling
      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

!-------------------------------------------------------------------------------
! For KREL = 1 (relativistic mode)
!
!  NPOTD = 2 * NATYPD
!  LMMAXD = 2 * (LMAXD+1)^2
!  NSPIND = 1
!  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green
!          function, set up in the spin-independent non-relativstic
!          (l,m_l)-representation
!
!-------------------------------------------------------------------------------
! .. Input variables
      Integer, Intent (In) :: ndim
      Integer, Intent (In) :: ichck
      Integer, Intent (In) :: lmmaxd !< (KREL+KORBIT+1)(LMAX+1)^2
      Integer, Intent (In) :: itermax
      Real (Kind=dp), Intent (In) :: errmax
! .. Input arrays
      Complex (Kind=dp), Dimension (ndim, ndim), Intent (In) :: m0
      Complex (Kind=dp), Dimension (ndim, ndim), Intent (In) :: ml
      Complex (Kind=dp), Dimension (ndim, ndim), Intent (In) :: mr
! .. Output variables
      Complex (Kind=dp), Dimension (ndim, ndim), Intent (Out) :: x
! .. Local Scalars
      Integer :: i, info, iter, j, n
      Real (Kind=dp) :: err, sum, xim, xre
! .. Local Arrays
      Integer, Dimension (ndim) :: ipvt
      Complex (Kind=dp), Dimension (ndim, ndim) :: aa, alfa, bb, beta, cc, &
        cunit, eps, tempin
      Complex (Kind=dp), Dimension (ndim, ndim) :: tempout, y1, y2
! .. External Subroutines ..
      External :: cinit, zaxpy, zcopy, zgemm, zgetrf, zgetrs
! .. Intrinsic Functions ..
      Intrinsic :: real, aimag
!     ..

      Call cinit(ndim*ndim, cunit)
      Do n = 1, ndim
        cunit(n, n) = cone
      End Do

      Call zcopy(ndim*ndim, m0, 1, eps, 1)
      Call zcopy(ndim*ndim, ml, 1, alfa, 1)
      Call zcopy(ndim*ndim, mr, 1, beta, 1)
      Call zcopy(ndim*ndim, m0, 1, x, 1)

      iter = 1
100   Continue

      Call zcopy(ndim*ndim, eps, 1, y1, 1)
      Call zcopy(ndim*ndim, y1, 1, tempin, 1)
      Call zgetrf(ndim, ndim, tempin, ndim, ipvt, info)

!     aa = eps^-1 * alfa
      Call zcopy(ndim*ndim, alfa, 1, tempout, 1)
      Call zgetrs('N', ndim, ndim, tempin, ndim, ipvt, tempout, ndim, info)
      Call zcopy(ndim*ndim, tempout, 1, aa, 1)

!     bb = eps^-1 * beta

      Call zcopy(ndim*ndim, beta, 1, tempout, 1)
      Call zgetrs('N', ndim, ndim, tempin, ndim, ipvt, tempout, ndim, info)
      Call zcopy(ndim*ndim, tempout, 1, bb, 1)

!     alfa_new = alfa * aa

      Call zgemm('N', 'N', ndim, ndim, ndim, cone, alfa, ndim, aa, ndim, &
        czero, y1, ndim)

!     beta_new = beta * bb

      Call zgemm('N', 'N', ndim, ndim, ndim, cone, beta, ndim, bb, ndim, &
        czero, y2, ndim)

!     cc = - alfa * bb

      Call zgemm('N', 'N', ndim, ndim, ndim, -cone, alfa, ndim, bb, ndim, &
        czero, cc, ndim)

!     x_new = x + cc

      Call zaxpy(ndim*ndim, cone, cc, 1, x, 1)

!     cc = eps + cc

      Call zaxpy(ndim*ndim, cone, cc, 1, eps, 1)

!     eps_new = cc - beta * aa

      Call zgemm('N', 'N', ndim, ndim, ndim, -cone, beta, ndim, aa, ndim, &
        cone, eps, ndim)

      Call zcopy(ndim*ndim, y1, 1, alfa, 1)
      Call zcopy(ndim*ndim, y2, 1, beta, 1)

      sum = 0.E0_dp
      Do i = 1, ndim
        Do j = 1, ndim
          xre = real(alfa(i,j))
          xim = aimag(alfa(i,j))
          sum = sum + xre*xre + xim*xim
        End Do
      End Do

      err = sqrt(sum)
      If (err<errmax .Or. iter>itermax) Go To 110
      iter = iter + 1
      Go To 100

110   Continue

      Call zcopy(ndim*ndim, x, 1, tempin, 1)
      Call zcopy(ndim*ndim, cunit, 1, tempout, 1)
      Call zgetrf(ndim, ndim, tempin, ndim, ipvt, info)
      Call zgetrs('N', ndim, ndim, tempin, ndim, ipvt, tempout, ndim, info)
      Call zcopy(ndim*ndim, tempout, 1, x, 1)

      Call zgemm('N', 'N', ndim, ndim, ndim, cone, x, ndim, mr, ndim, czero, &
        tempin, ndim)
      Call zgemm('N', 'N', ndim, ndim, ndim, cone, ml, ndim, tempin, ndim, &
        czero, x, ndim)

      If (iter>itermax) Then
        Write (6, Fmt='('' itermax too small.  iter='',i3)') iter
        Write (6, '('' Surfgf:  iter='',i4,''  error='',d14.7)') iter, err
      End If
      If (ichck==0) Return
!      write(6,'('' Surfgf:  iter='',i4,''  error='',d12.7)') iter,err
!      write(6,'(/'' X matrix'')')
!      write(6,*)
!
      Return
    End Subroutine
