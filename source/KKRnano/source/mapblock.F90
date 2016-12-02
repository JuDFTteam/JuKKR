      integer function mapblock(itercurr, iterfirst, iterlast, iterstep, nodefirst, nodelast)
!-----------------------------------------------------------------------
!
!     Maps loop iteration to nodes (once a call).
!     The mapping is done in blocks, such that:
!       1) all nodes work on q consecutive iterations
!          (equally sized blocks)
!       2) the first r nodes get (q+1) iterations and (n-r) nodes get
!          q iterations
!
!     Descrition of input parameters:
!
!       itercurr  : index of current iteration
!       iterfirst : index of first iteration
!       iterlast  : index of last iteration
!       iterstep  : iteration step
!       nodefirst : index of first node
!       nodelast  : index of last node
!
!     where:
!
!       iterfirst <= itercurr <= iterlast
!       nodefirst <= nodelast
!       iterationstep > 0
!
!
!     Example for Intel NX-calls:
!
!      DO i=iterstart, iterend, iterstep
!
!         node = mapblock(i,iterstart,iterend,iterstep,0,numnodes()-1)
!
!         IF (mynode() .EQ. node) THEN
!                  ! the code for an iteration
!         ENDIF
!
!      ENDDO
!
!
!     Algorithm:
!
!              |  number_iterations  |
!     let q := |  -----------------  |
!              |    number_nodes     |
!              ---                 ---
!
!     number_nodes = q * number_nodes + r
!                  = (q + 1) * r + q * (number_nodes - r)
!
!     IF q = number_iterations / number_nodes
!       THEN
!         you have equally sized blocks
!       ELSE
!         you have r blocks of size q+1
!         and (number_nodes - r) blocks of size q
!
!
!                                           Rudolf Berrendorf, July 1992
!                                           last update: February 1994
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: itercurr, iterfirst, iterlast, iterstep, nodefirst, nodelast
      
      integer :: litercurr, literlast, lnodelast, nprime, q, r

!---  normalize ranges
      literlast = max(0, (iterlast - iterfirst + iterstep)/iterstep)
      litercurr = ((itercurr - iterfirst)/iterstep) + 1
      lnodelast = nodelast - nodefirst + 1


      q = literlast/lnodelast ! calculate blocksize

!---  calculate iteration mapping
      if (q*lnodelast == literlast) then

        mapblock = (nodefirst - 1) + ceiling(dble(litercurr)/dble(q)) ! equally sized blocks

      else
        r = literlast - (q*lnodelast) ! unequal blocks
        nprime = (q + 1)*r ! up to nprime blocks of size (q+1)

        if (litercurr <= nprime) then
          mapblock = (nodefirst - 1) + ceiling(dble(litercurr)/dble(q + 1))
        else
          mapblock = (nodefirst - 1) + ceiling(dble(litercurr - nprime)/dble(q) + dble(r))
        endif

      endif

      if ((mapblock < nodefirst) .or. (mapblock > nodelast)) then
        write (6, fmt=*) 'internal error in mapblock'
#ifdef __GFORTRAN__
#define PR(VAR) "  $=",VAR
#else
#define PR(VAR) "  ",#VAR," =",VAR
#endif
        write (6, fmt='(99(3a,i0))') PR(mapblock) ! result
        write (6, fmt='(99(3a,i0))') PR(itercurr),PR(iterfirst),PR(iterlast),PR(iterstep),PR(nodefirst),PR(nodelast) ! arguments
        write (6, fmt='(99(3a,i0))') PR(literlast),PR(litercurr),PR(lnodelast),PR(q),PR(r),PR(nprime) ! internal numbers
#undef  PR
        stop
      endif

      endfunction mapblock
