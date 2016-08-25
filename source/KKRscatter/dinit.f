      subroutine dinit(array,leng)
!- Initializes double precision array to zero
!----------------------------------------------------------------------
! Inputs:
!   leng  :length of array
! Outputs:
!   array :double precision array set to zero
!----------------------------------------------------------------------
      implicit none
! Passed parameters:
      integer leng
      double precision array(leng)
! Local parameters:
      integer i,m,mp1
!----------------------------------------------------------------------
      m = mod(leng,5)
      if( m .ne. 0 ) then
         do i = 1,m
            array(i) = 0.d0
         enddo
         if( leng .lt. 5 ) return
      endif
      mp1 = m + 1
      do i = mp1,leng,5
         array(i) = 0.d0
         array(i + 1) = 0.d0
         array(i + 2) = 0.d0
         array(i + 3) = 0.d0
         array(i + 4) = 0.d0
      enddo
      end
