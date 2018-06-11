!module mod_cheb

!contains 

subroutine getcmatrix(ncheb, cmatrix)
! calculates the C matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
  use mod_DataTypes
  implicit none
  integer, intent (in) :: ncheb
  real (kind=dp), intent (out) :: cmatrix(0:ncheb, 0:ncheb)
  real (kind=dp) :: pi
!local
  integer :: icheb1, icheb2

  pi = 4d0*datan(1d0)
  do icheb1 = 0, ncheb
    do icheb2 = 0, ncheb
! maybe incorrect
      cmatrix(icheb2, icheb1) = dcos(icheb1*pi*((ncheb- &
        icheb2)+0.5d0)/(ncheb+1))
    end do
  end do
end subroutine


subroutine getcinvmatrix(ncheb, cinvmatrix)
! calculates the C**-1 matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
  use mod_DataTypes
  implicit none
  integer, intent (in) :: ncheb
  real (kind=dp), intent (out) :: cinvmatrix(0:ncheb, 0:ncheb)
!local
  real (kind=dp) :: pi
  integer :: icheb1, icheb2
  real (kind=dp) :: fac

  pi = 4d0*datan(1d0)
  fac = 1.0d0/(ncheb+1)
  do icheb1 = 0, ncheb
    do icheb2 = 0, ncheb
      cinvmatrix(icheb1, icheb2) = fac*dcos(icheb1*pi*((ncheb- &
        icheb2)+0.5d0)/(ncheb+1))
    end do
    fac = 2.0d0/(ncheb+1)
  end do
end subroutine


subroutine getccmatrix(ncheb, rmesh, nrmesh, cmatrix)
! calculates the C matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
  use mod_DataTypes
  implicit none
  integer, intent (in) :: ncheb, nrmesh
  real (kind=dp), intent (in) :: rmesh(nrmesh)
  real (kind=dp), intent (out) :: cmatrix(1:nrmesh, 0:ncheb)
  integer :: icheb, ir

  do ir = 1, nrmesh
    do icheb = 0, ncheb
      cmatrix(ir, icheb) = cos(real(icheb, kind=dp)*dacos(rmesh(ir)))
    end do
  end do
end subroutine


subroutine getlambda(ncheb, lambda)
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion 
  use mod_DataTypes
  implicit none
  integer, intent (in) :: ncheb
  real (kind=dp), intent (out) :: lambda(0:ncheb, 0:ncheb)
!local
  integer :: icheb, icheb2

  do icheb2 = 1, ncheb, 2
    lambda(0, icheb2) = icheb2
  end do
  do icheb = 1, ncheb
    do icheb2 = icheb + 1, ncheb, 2
      lambda(icheb, icheb2) = icheb2*2
    end do
  end do
end subroutine


subroutine getclambdacinv(ncheb, clambdacinv)
  use mod_DataTypes
  implicit none
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion
  integer :: ncheb
  real (kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)
!local
  real (kind=dp) :: lambda(0:ncheb, 0:ncheb)
  real (kind=dp) :: cmatrix(0:ncheb, 0:ncheb)
  real (kind=dp) :: cinvmatrix(0:ncheb, 0:ncheb)
  real (kind=dp) :: temp1(0:ncheb, 0:ncheb)
  integer :: n

  lambda = (0.0d0, 0.0d0)
  cmatrix = (0.0d0, 0.0d0)
  cinvmatrix = (0.0d0, 0.0d0)
  lambda = (0.0d0, 0.0d0)
  temp1 = (0.0d0, 0.0d0)

  call getlambda(ncheb, lambda)
  call getcinvmatrix(ncheb, cinvmatrix)
  call getcmatrix(ncheb, cmatrix)
  n = ncheb + 1
  call dgemm('N', 'N', n, n, n, 1d0, lambda, n, cinvmatrix, n, 0d0, temp1, n)
  call dgemm('N', 'N', n, n, n, 1d0, cmatrix, n, temp1, n, 0d0, clambdacinv, n)
end subroutine


subroutine getclambda2cinv(ncheb, clambda2cinv)
  use mod_DataTypes
  implicit none
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion
  integer :: ncheb
  real (kind=dp) :: clambda2cinv(0:ncheb, 0:ncheb)
!local
  real (kind=dp) :: lambda(0:ncheb, 0:ncheb)
  real (kind=dp) :: cmatrix(0:ncheb, 0:ncheb)
  real (kind=dp) :: cinvmatrix(0:ncheb, 0:ncheb)
  real (kind=dp) :: temp1(0:ncheb, 0:ncheb)
  real (kind=dp) :: temp2(0:ncheb, 0:ncheb)

  lambda = (0.0d0, 0.0d0)
  cmatrix = (0.0d0, 0.0d0)
  cinvmatrix = (0.0d0, 0.0d0)
  lambda = (0.0d0, 0.0d0)
  temp1 = (0.0d0, 0.0d0)

  call getlambda(ncheb, lambda)
  call getcinvmatrix(ncheb, cinvmatrix)
  call getcmatrix(ncheb, cmatrix)

  call matmat_dmdm(lambda, lambda, ncheb, temp1)
  call matmat_dmdm(temp1, cinvmatrix, ncheb, temp2)
  call matmat_dmdm(cmatrix, temp2, ncheb, clambda2cinv)
end subroutine


subroutine diffcheb(fn, ncheb, dfndr)
  use mod_DataTypes
  implicit none
  integer :: ncheb
  real (kind=dp) :: fn(0:ncheb)
  real (kind=dp) :: dfndr(0:ncheb)
  real (kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)

!needs to be checked!!!!!!1
  call getclambdacinv(ncheb, clambdacinv(0:ncheb,0:ncheb))
  call matvec_dmdm(ncheb, clambdacinv(0:ncheb,0:ncheb), fn(0:ncheb), dfndr(0:ncheb))
end subroutine


! helper functions

subroutine matvec_dmdm(ncheb, mat1, vec1, outvec)
  use mod_DataTypes
  implicit none
  integer, intent (in) :: ncheb
  real (kind=dp), intent (in) :: mat1(0:ncheb, 0:ncheb), vec1(0:ncheb)
  real (kind=dp), intent (out) :: outvec(0:ncheb)
  integer :: n, m

  m = size(mat1, 1)
  n = size(mat1, 2)
  if (size(vec1,1)/=n) stop &
    'matvec_dmdm: dimensions of first input array differ.'
  call dgemv('N', m, n, 1.0d0, mat1, m, vec1, 1, 0.0d0, outvec, 1)
end subroutine matvec_dmdm

subroutine matmat_dmdm(mat1, mat2, ncheb, outmat)
  use mod_DataTypes
  implicit none
  integer, intent(in) :: ncheb
  real (kind=dp), intent (in) :: mat1(0:ncheb, 0:ncheb), mat2(0:ncheb, 0:ncheb)
  real (kind=dp), intent (out) :: outmat(0:ncheb, 0:ncheb)
  
  integer :: n

  n = ncheb + 1
  call dgemm('N', 'N', n, n, n, 1d0, mat1, n, mat2, n, 0d0, outmat, n)
end subroutine matmat_dmdm


!end module mod_cheb
