SUBROUTINE dmpy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
!- Matrix multiplication:  c = a * b
! ----------------------------------------------------------------
!i Inputs:
!i   a     :left matrix to multiply
!i   nca   :spacing between elements in adjacent columns in matrix a
!i   nra   :spacing between elements in adjacent rows in matrix a
!i   b     :right matrix to multiply
!i   ncb   :spacing between elements in adjacent columns in matrix b
!i   nrb   :spacing between elements in adjacent rows in matrix b
!i   ncc   :spacing between elements in adjacent columns in matrix c
!i   nrc   :spacing between elements in adjacent rows in matrix c
!i   n     :the number of rows to calculate
!i   m     :the number of columns to calculate
!i   l     :length of vector for matrix multiply
!o Outputs:
!o   c     :product matrix
!r Remarks:
!r   This is a general-purpose matrix multiplication routine,
!r   multiplying a subblock of matrix a by a subblock of matrix b.
!r   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
!r   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
!r   Arrays are Locally one-dimensional so as to optimize inner loop,
!r   which is executed n*m*l times.  No attempt is made to optimize
!r   the outer loops, executed n*m times.
!r     Examples: product of (n,l) subblock of a into (l,m) subblock of b
!r   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
!r     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
!r     To generate the tranpose of that product, use:
!r   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
! ----------------------------------------------------------------

      implicit none 
! Passed parameters                                                     
      integer nca,nra,ncb,nrb,ncc,nrc,n,m,l 
      double precision a(0:*), b(0:*), c(0:*) 
! Local parameters                                                      
      double precision sum 
      integer i,j,k,nakpi,nbjpk 
!
!#ifdefC CRAY
!      CALL MXMA(A,NRA,NCA,B,NRB,NCB,C,NRC,NCC,N,L,M)
!#else
DO i = n-1, 0, -1
  DO j = m-1, 0, -1
    sum = 0.d0
    nakpi = nra*i
    nbjpk = ncb*j
    DO k = l-1, 0, -1
      sum = sum + a(nakpi)*b(nbjpk)
      nakpi = nakpi + nca
      nbjpk = nbjpk + nrb
    END DO
    c(i*nrc+j*ncc) = sum
  END DO
END DO
!#endif
END SUBROUTINE dmpy
