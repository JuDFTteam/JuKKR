c     @PROCESS HOT=noarraypad:level=1:simd:vector
C     Taken from SPARSKIT
C     modified for double complex and for
C     multiplying several vectors at once


      module vbrmv_mat_mod
      contains

C     kvstc ... dimension: number of matrix block columns + 1

c-----------------------------------------------------------------------
      subroutine vbrmv_mat(nr, ia, ja, ka, a, kvstr, kvstc, x, b, ncols)
c-----------------------------------------------------------------------
      implicit none
      integer nr, ia(nr+1), ja(:), ka(:), kvstr(nr+1), kvstc(*)
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
      integer i, j, k, istart, istop
c     double complex  xjj(ncols)
      integer icols

c     does not seem right...
c     should be:
c     n = kvstr(nr+1)-1
c---------------------------------
c     n = kvstc(nc+1)-1
c     do i = 1, n
c        b(i) = 0.d0
c     enddo
c---------------------------------
      ! TODO: replace hardcoded stuff
cIBM* ALIGN(32, buffer)
      double complex buffer(16*43,16)
      integer nrowbuf, rowbuf, sum_nrowbuf
      integer startrow, leaddim_b, num_rows, leaddim_buffer
      integer jlow, jhigh

      double complex, parameter :: CZERO = (0.0d0, 0.0d0)
      double complex, parameter :: CONE  = (1.0d0, 0.0d0)

      leaddim_b = size(b, 1)
      leaddim_buffer = size(buffer, 1)

      b = 0.d0

c     can parallelise this loop

C$OMP PARALLEL PRIVATE(i, istart, istop,num_rows,rowbuf, 
C$OMP&                  sum_nrowbuf,j,jlow,jhigh,
C$OMP&                  startrow,nrowbuf,icols,buffer,k)

C$OMP DO
      do i = 1, nr
         istart = kvstr(i)
         istop  = kvstr(i+1)-1
         num_rows = istop - istart + 1
         rowbuf = 1
         sum_nrowbuf = 0

         jlow = ia(i)
         jhigh = ia(i+1)-1

         k = ka(jlow)

         do j = jlow, jhigh
             startrow = kvstc(ja(j))
             nrowbuf = kvstc(ja(j)+1) - startrow
cIBM* ASSERT(ITERCNT(16))
             do icols = 1, ncols
                !call ZCOPY(nrowbuf, x(startrow, icols), 1, 
     &          !               buffer(rowbuf,   icols), 1)
                buffer(rowbuf:(rowbuf + nrowbuf -1),      icols) =
     &               x(startrow:(startrow + nrowbuf - 1), icols) 
             enddo

             sum_nrowbuf = sum_nrowbuf + nrowbuf
             rowbuf = rowbuf + nrowbuf
         enddo

         !k = ka(ia(i))

         call ZGEMM('N','N',num_rows,ncols,sum_nrowbuf, 
     &              CONE,a(k),num_rows,  
     &              buffer,leaddim_buffer,  
     &              CZERO,b(istart,1),leaddim_b)

      enddo      
C$OMP END DO

C$OMP END PARALLEL

c---------------------------------
      return
      end subroutine
c-----------------------------------------------------------------------
c----------------------end-of-vbrmv-------------------------------------
c-----------------------------------------------------------------------

      end module
