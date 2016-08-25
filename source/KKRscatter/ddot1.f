      double precision function ddot1(n,dx,incx,dy,incy)
!- Forms the dot product of two vectors.
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :lenght of dx and dy
!i   dx    :first vector to mutiply
!i   incx  :incrementation for x
!i   dy    :second vector to mutiply
!i   incy  :incrementation for y
!o Outputs:
!o   ddot  :dot product of two vectors
!r Remarks:
!r    Adapted from: jack dongarra, linpack, 3/11/78.
! ----------------------------------------------------------------------
      implicit none
! Passed parameters:
      integer incx,incy,n
      double precision dx(*),dy(*)
! Local parameters:
      double precision dtemp
      integer i,ix,iy,m,mp1
!
      ddot1 = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.ne.1.or.incy.ne.1)then
! ----- code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do i = 1,n
          dtemp = dtemp + dx(ix)*dy(iy)
          ix = ix + incx
          iy = iy + incy
        enddo
        ddot1 = dtemp
      else
! ----- code for both increments equal to 1
        m = mod(n,5)
        if( m .ne. 0 ) then
          do i = 1,m
            dtemp = dtemp + dx(i)*dy(i)
          enddo
          if( n .lt. 5 ) go to 60
        endif
        mp1 = m + 1
        do i = mp1,n,5
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +                   &
     &          dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
        enddo
   60   ddot1 = dtemp
      endif
      end
