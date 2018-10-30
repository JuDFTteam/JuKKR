!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dmpy
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Matrix-matrix multiplication
  !> Author: 
  !> Category: KKRhost, numerical-tools
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !> @note only used by symlat which seems to be unused @endnote
  !>
  !> Matrix multiplication:  c = a * b
  !> 
  !> i Inputs:
  !> i   a     :left matrix to multiply
  !> i   nca   :spacing between elements in adjacent columns in matrix a
  !> i   nra   :spacing between elements in adjacent rows in matrix a
  !> i   b     :right matrix to multiply
  !> i   ncb   :spacing between elements in adjacent columns in matrix b
  !> i   nrb   :spacing between elements in adjacent rows in matrix b
  !> i   ncc   :spacing between elements in adjacent columns in matrix c
  !> i   nrc   :spacing between elements in adjacent rows in matrix c
  !> i   n     :the number of rows to calculate
  !> i   m     :the number of columns to calculate
  !> i   l     :length of vector for matrix multiply
  !> o Outputs:
  !> o   c     :product matrix
  !> r Remarks:
  !> r   This is a general-purpose matrix multiplication routine,
  !> r   multiplying a subblock of matrix a by a subblock of matrix b.
  !> r   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
  !> r   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
  !> r   Arrays are Locally one-dimensional so as to optimize inner loop,
  !> r   which is executed n*m*l times.  No attempt is made to optimize
  !> r   the outer loops, executed n*m times.
  !> r     Examples: product of (n,l) subblock of a into (l,m) subblock of b
  !> r   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
  !> r     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
  !> r     To generate the tranpose of that product, use:
  !> r   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
  !-------------------------------------------------------------------------------
  subroutine dmpy(a, nca, nra, b, ncb, nrb, c, ncc, nrc, n, m, l)

    implicit none
    ! Passed parameters
    integer :: nca, nra, ncb, nrb, ncc, nrc, n, m, l
    real (kind=dp) :: a(0:*), b(0:*), c(0:*)
    ! Local parameters
    real (kind=dp) :: sum
    integer :: i, j, k, nakpi, nbjpk

    ! #ifdefC CRAY
    ! CALL MXMA(A,NRA,NCA,B,NRB,NCB,C,NRC,NCC,N,L,M)
    ! #else
    do i = n - 1, 0, -1
      do j = m - 1, 0, -1
        sum = 0.e0_dp
        nakpi = nra*i
        nbjpk = ncb*j
        do k = l - 1, 0, -1
          sum = sum + a(nakpi)*b(nbjpk)
          nakpi = nakpi + nca
          nbjpk = nbjpk + nrb
        end do
        c(i*nrc+j*ncc) = sum
      end do
    end do
    ! #endif
  end subroutine dmpy

end module mod_dmpy
