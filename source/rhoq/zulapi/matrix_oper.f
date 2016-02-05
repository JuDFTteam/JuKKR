      subroutine invert_matrix(a, ainv, n)

      implicit none
      integer, intent(in)        ::  n
      complex,    intent(in)     ::  a(n,n)
      complex,    intent(inout)  ::  ainv(n,n)

      integer                    ::  info, ipiv(n), lwork
      complex, allocatable       ::  work(:)
      ainv = a
      lwork = 64 * n
      allocate(work(lwork))

      call ZGETRF(n, n, ainv, n, ipiv, info)
      if (info.ne.0) then
        print *, info
        stop 'error in matrix multiplication: Zgetrf'
      end if
      call ZGETRI(n, ainv, n, ipiv, work, lwork, info)
      if (info.ne.0) stop 'error in matrix multiplication: Zgetri'

      deallocate(work)
      end subroutine invert_matrix


      subroutine print_matrix(a,n,m)
      implicit none
      integer, intent(in)  ::  n,m
      complex, intent(in)  ::  a(n,m)

      integer  ::  k
      do k = 1, n
         print '(4(f0.2,2x))', a(k,:)
      end do
      end subroutine print_matrix
