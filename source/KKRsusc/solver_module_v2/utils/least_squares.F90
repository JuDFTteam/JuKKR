  subroutine least_squares(m,n,a,b,chi2)
! Solution of Ax = b by SVD: A = U sigma V^dagger, U and V unitary
! A is destroyed, x returned in b
! m >= n
  use global, only: i4b, r8b, c8b

  implicit none

! m data points, n unknowns
  integer(kind=i4b), intent(in)    :: m, n
! design matrix, rhs then x
  complex(kind=c8b), intent(inout) :: a(m,n), b(m)
! approximate chi-square
  real(kind=r8b),    intent(out)   :: chi2
! ------------------------------------------
! tolerance for small singular values
  real(kind=r8b),    parameter :: reltol = 1.d-8, abstol = 1.d-6
  complex(kind=c8b), parameter :: one = (1.d0,0.d0), zero = (0.d0,0.d0), minus = (-1.d0,0.d0)
  complex(kind=c8b) :: v(n,n), udummy, work(2*m+n), bsave(m)
  real(kind=r8b)    :: sigma(n), maxsigma, minsigma, rwork(5*n)
  integer(kind=i4b) :: info

  if (m < n) stop 'least_squares: m < n'
! computing the SVD of matrix A
!  write(iodb,'("least_squares: SVD")')
  call zgesvd('O','S',m,n,a,m,sigma,udummy,m,v,n,work,2*m+n,rwork,info)
  if (info /= 0) stop 'least_squares: fail in SVD'
!  write(iodb,'("singular values:",100es16.8)') sigma
! backsubstitution
!  write(iodb,'("least_squares: backsubstitution")')
! save rhs and singular values
  bsave = b; rwork(1:n) = sigma
! set small singular values to zero
  maxsigma = maxval(sigma)
  minsigma = minval(sigma)
  where (rwork(1:n) < reltol*maxsigma) rwork = 0.d0
!  write(iodb,'("least_squares: inv cond number=",es10.1)') minsigma/maxsigma
  where (sigma > reltol*maxsigma)
    sigma = 1.d0/sigma
  elsewhere
    sigma = 0.d0
  end where
!  write(iodb,'("singular values:",100es16.8)') sigma
! multiply b by U^dagger, put result in work
  call zgemv('C',m,n,one,a,m,b,1,zero,work,1)
! now rescale with singular values
  work(1:n) = sigma*work(1:n)
! multiply by V, put result in b
  call zgemv('C',n,n,one,v,n,work,1,zero,b,1)
  b(n+1:m) = zero
! Now chi-square residual
! multiply x by V^dagger, put result in work
  call zgemv('N',n,n,one,v,n,b,1,zero,work,1)
! now rescale with singular values
  work(1:n) = rwork(1:n)*work(1:n)
! multiply by U, subtract from b
  call zgemv('N',m,n,one,a,m,work,1,minus,bsave,1)
! final result
  chi2 = dot_product(bsave,bsave)
!  write(iodb,'("least_squares: chi2=",es10.1)') chi2
! All done!
  end subroutine least_squares

