!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: Functions and subroutines for Chebychev radial mesh
!> Author: 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> Module providing functionality needed to used Chebychev radial mesh, which 
!> is used in the rllsll-solver (see PhD thesis of David Bauer for details)
!-------------------------------------------------------------------------------
module mod_cheb

  private
  public :: getcmatrix, getcinvmatrix, getccmatrix, getlambda, getclambdacinv, getclambda2cinv, diffcheb

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate C matrix of Chebychev mesh
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the C matrix according to:
  !> Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
  !-------------------------------------------------------------------------------
  subroutine getcmatrix(ncheb, cmatrix)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: pi
    implicit none
    integer, intent (in) :: ncheb
    real (kind=dp), intent (out) :: cmatrix(0:ncheb, 0:ncheb)
    ! local
    integer :: icheb1, icheb2

    do icheb1 = 0, ncheb
      do icheb2 = 0, ncheb
        ! maybe incorrect
        cmatrix(icheb2, icheb1) = cos(icheb1*pi*((ncheb-icheb2)+0.5e0_dp)/(ncheb+1))
      end do
    end do
  end subroutine getcmatrix


  !-------------------------------------------------------------------------------
  !> Summary: Calculate inverse C matrix of Chebychev mesh
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the C**-1 matrix according to:
  !> Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
  !-------------------------------------------------------------------------------
  subroutine getcinvmatrix(ncheb, cinvmatrix)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: pi
    implicit none
    integer, intent (in) :: ncheb
    real (kind=dp), intent (out) :: cinvmatrix(0:ncheb, 0:ncheb)
    ! local
    integer :: icheb1, icheb2
    real (kind=dp) :: fac

    fac = 1.0e0_dp/(ncheb+1)
    do icheb1 = 0, ncheb
      do icheb2 = 0, ncheb
        cinvmatrix(icheb1, icheb2) = fac*cos(icheb1*pi*((ncheb-icheb2)+0.5e0_dp)/(ncheb+1))
      end do
      fac = 2.0e0_dp/(ncheb+1)
    end do
  end subroutine getcinvmatrix


  !-------------------------------------------------------------------------------
  !> Summary: Get CC-matrix of Chebychev mesh
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the C matrix according to:
  !> Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
  !>
  !> @note
  !> Similar to getcmatrix but has form cos(x)*acos(r) instead of cos(x)
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine getccmatrix(ncheb, rmesh, nrmesh, cmatrix)
    use :: mod_datatypes, only: dp
    implicit none
    integer, intent (in) :: ncheb, nrmesh
    real (kind=dp), intent (in) :: rmesh(nrmesh)
    real (kind=dp), intent (out) :: cmatrix(1:nrmesh, 0:ncheb)
    integer :: icheb, ir

    do ir = 1, nrmesh
      do icheb = 0, ncheb
        cmatrix(ir, icheb) = cos(real(icheb,kind=dp)*acos(rmesh(ir)))
      end do
    end do
  end subroutine getccmatrix


  !-------------------------------------------------------------------------------
  !> Summary: Get Lambda-matrix
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Set up the Lambda matrix which differentiates the coefficients of an
  ! Chebyshev expansion
  !-------------------------------------------------------------------------------
  subroutine getlambda(ncheb, lambda)
    use :: mod_datatypes, only: dp
    implicit none
    integer, intent (in) :: ncheb
    real (kind=dp), intent (out) :: lambda(0:ncheb, 0:ncheb)
    ! local
    integer :: icheb, icheb2

    do icheb2 = 1, ncheb, 2
      lambda(0, icheb2) = icheb2
    end do
    do icheb = 1, ncheb
      do icheb2 = icheb + 1, ncheb, 2
        lambda(icheb, icheb2) = icheb2*2
      end do
    end do
  end subroutine getlambda


  !-------------------------------------------------------------------------------
  !> Summary: Computes C.Lambda.C^-1 matrix for chebycheb differentiation
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Set up the product of C-matrix, Lambda-matrix and C^-1-matrix
  !-------------------------------------------------------------------------------
  subroutine getclambdacinv(ncheb, clambdacinv)
    use :: mod_datatypes, only: dp
    implicit none
    integer :: ncheb
    real (kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)
    ! local
    real (kind=dp) :: lambda(0:ncheb, 0:ncheb)
    real (kind=dp) :: cmatrix(0:ncheb, 0:ncheb)
    real (kind=dp) :: cinvmatrix(0:ncheb, 0:ncheb)
    real (kind=dp) :: temp1(0:ncheb, 0:ncheb)
    integer :: n

    lambda = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    cmatrix = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    cinvmatrix = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    lambda = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    temp1 = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)

    call getlambda(ncheb, lambda)
    call getcinvmatrix(ncheb, cinvmatrix)
    call getcmatrix(ncheb, cmatrix)
    n = ncheb + 1
    call dgemm('N', 'N', n, n, n, 1e0_dp, lambda, n, cinvmatrix, n, 0e0_dp, temp1, n)
    call dgemm('N', 'N', n, n, n, 1e0_dp, cmatrix, n, temp1, n, 0e0_dp, clambdacinv, n)
  end subroutine getclambdacinv


  !-------------------------------------------------------------------------------
  !> Summary: Computes C.Lambda^2.C^-1 matrix for chebycheb differentiation
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Set up the product of C-matrix, Lambda-matrix^2 and C^-1-matrix
  !> @note
  !> Similar to getclambdacinv but for squared Lambda-matrix
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine getclambda2cinv(ncheb, clambda2cinv)
    use :: mod_datatypes, only: dp
    implicit none
    integer :: ncheb
    real (kind=dp) :: clambda2cinv(0:ncheb, 0:ncheb)
    ! local
    real (kind=dp) :: lambda(0:ncheb, 0:ncheb)
    real (kind=dp) :: cmatrix(0:ncheb, 0:ncheb)
    real (kind=dp) :: cinvmatrix(0:ncheb, 0:ncheb)
    real (kind=dp) :: temp1(0:ncheb, 0:ncheb)
    real (kind=dp) :: temp2(0:ncheb, 0:ncheb)

    lambda = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    cmatrix = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    cinvmatrix = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    lambda = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)
    temp1 = cmplx(0.0e0_dp, 0.0e0_dp, kind=dp)

    call getlambda(ncheb, lambda)
    call getcinvmatrix(ncheb, cinvmatrix)
    call getcmatrix(ncheb, cmatrix)

    call matmat_dmdm(lambda, lambda, ncheb, temp1)
    call matmat_dmdm(temp1, cinvmatrix, ncheb, temp2)
    call matmat_dmdm(cmatrix, temp2, ncheb, clambda2cinv)
  end subroutine getclambda2cinv


  !-------------------------------------------------------------------------------
  !> Summary: Chebychev differentiation
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Computed differenction in Chebychev mesh by computing Lambda matrix 
  !> and multiply this to input matrix
  !-------------------------------------------------------------------------------
  subroutine diffcheb(fn, ncheb, dfndr)
    use :: mod_datatypes, only: dp
    implicit none
    integer :: ncheb
    real (kind=dp) :: fn(0:ncheb)
    real (kind=dp) :: dfndr(0:ncheb)
    real (kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)

    !> needs to be checked
    call getclambdacinv(ncheb, clambdacinv(0:ncheb,0:ncheb))
    call matvec_dmdm(ncheb, clambdacinv(0:ncheb,0:ncheb), fn(0:ncheb), dfndr(0:ncheb))
  end subroutine diffcheb


  ! helper functions:

  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for double-precision matrix-vector multiplication 
  !> Author: 
  !> Category: KKRhost, sanity-check, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Helper function that wraps dgemv after sanity check of input
  !-------------------------------------------------------------------------------
  subroutine matvec_dmdm(ncheb, mat1, vec1, outvec)
    use :: mod_datatypes, only: dp
    implicit none
    integer, intent (in) :: ncheb
    real (kind=dp), intent (in) :: mat1(0:ncheb, 0:ncheb), vec1(0:ncheb)
    real (kind=dp), intent (out) :: outvec(0:ncheb)
    integer :: n, m

    m = size(mat1, 1)
    n = size(mat1, 2)
    if (size(vec1,1)/=n) stop 'matvec_dmdm: dimensions of first input array differ.'
    call dgemv('N', m, n, 1.0e0_dp, mat1, m, vec1, 1, 0.0e0_dp, outvec, 1)
  end subroutine matvec_dmdm

  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for double-precision matrix-matrix multiplication
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Helper function wrapping dgemm
  !-------------------------------------------------------------------------------
  subroutine matmat_dmdm(mat1, mat2, ncheb, outmat)
    use :: mod_datatypes, only: dp
    implicit none
    integer, intent (in) :: ncheb
    real (kind=dp), intent (in) :: mat1(0:ncheb, 0:ncheb), mat2(0:ncheb, 0:ncheb)
    real (kind=dp), intent (out) :: outmat(0:ncheb, 0:ncheb)

    integer :: n

    n = ncheb + 1
    call dgemm('N', 'N', n, n, n, 1e0_dp, mat1, n, mat2, n, 0e0_dp, outmat, n)
  end subroutine matmat_dmdm

end module mod_cheb
