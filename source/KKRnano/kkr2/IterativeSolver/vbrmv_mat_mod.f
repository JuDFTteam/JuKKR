C     Taken from SPARSKIT
C     modified for double complex and for
C     multiplying several vectors at once


      module vbrmv_mat_mod
      contains

C     kvstc ... dimension: number of matrix block columns + 1

c-----------------------------------------------------------------------
      subroutine vbrmv_mat(nr, ia, ja, a, kvstr, kvstc, x, b, ncols)
c-----------------------------------------------------------------------
      integer nr, ia(nr+1), ja(*), kvstr(nr+1), kvstc(*)
      integer ncols
      double complex  a(*), x(:,:), b(:,:)
c-----------------------------------------------------------------------
c     Sparse matrix-full vector product, in VBR format.
c-----------------------------------------------------------------------
c     On entry:
c--------------
c     nr      = number of block rows in matrix A
c     ia,ja,a,kvstr,kvstc = matrix A in variable block row format
c     x       = multiplier vector in full format
c
c     ncols = number of columns of matrix A
c
c     On return:
c---------------
c     b = product of matrix A times vector x in full format
c
c     Algorithm:
c---------------
c     Perform multiplication by traversing a in order.
c
c-----------------------------------------------------------------------
c-----local variables
      integer i, j, ii, jj, k, istart, istop
      double complex  xjj(ncols)

c     does not seem right...
c     should be:
c     n = kvstr(nr+1)-1
c---------------------------------
c     n = kvstc(nc+1)-1
c     do i = 1, n
c        b(i) = 0.d0
c     enddo
c---------------------------------
      b = 0.d0

      k = 1
c     can parallelise this loop-
      do i = 1, nr
         istart = kvstr(i)
         istop  = kvstr(i+1)-1
         do j = ia(i), ia(i+1)-1
            do jj = kvstc(ja(j)), kvstc(ja(j)+1)-1
               xjj = x(jj,:)
               do ii = istart, istop
                  b(ii,:) = b(ii,:) + xjj*a(k)
                  k = k + 1
               enddo
            enddo
         enddo
      enddo
c---------------------------------
      return
      end subroutine
c-----------------------------------------------------------------------
c----------------------end-of-vbrmv-------------------------------------
c-----------------------------------------------------------------------

      end module
