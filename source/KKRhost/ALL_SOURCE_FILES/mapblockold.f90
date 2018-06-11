    Integer Function mapblockold(itercurr, iterfirst, iterlast, iterstep, &
      nodefirst, nodelast)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------

!     Maps loop iteration to nodes (once a call).
!     The mapping is done in blocks, such that:
!       1) all nodes work on q consecutive iterations
!          (equally sized blocks)
!       2) the first r nodes get (q+1) iterations and (n-r) nodes get
!          q iterations

!     Descrition of input parameters:

!       itercurr  : index of current iteration
!       iterfirst : index of first iteration
!       iterlast  : index of last iteration
!       iterstep  : iteration step
!       nodefirst : index of first node
!       nodelast  : index of last node

!     where:

!       iterfirst <= itercurr <= iterlast
!       nodefirst <= nodelast
!       iterationstep > 0


!     Example for Intel NX-calls:

!      DO i=iterstart, iterend, iterstep

!         node = mapblock(i,iterstart,iterend,iterstep,0,numnodes()-1)

!         IF (mynode() .EQ. node) THEN
!                  ! the code for an iteration
!         ENDIF

!      ENDDO


!     Algorithm:

!              |  number_iterations  |
!     let q := |  -----------------  |
!              |    number_nodes     |
!              ---                 ---

!     number_nodes = q * number_nodes + r
!                  = (q + 1) * r + q * (number_nodes - r)

!     IF q = number_iterations / number_nodes
!       THEN
!         you have equally sized blocks
!       ELSE
!         you have r blocks of size q+1
!         and (number_nodes - r) blocks of size q


!                                           Rudolf Berrendorf, July 1992
!                                           last update: February 1994
!-----------------------------------------------------------------------
!.. Scalar Arguments ..
      Integer :: itercurr, iterfirst, iterlast, iterstep, nodefirst, nodelast
!..
!.. Intrinsic Functions ..

      Intrinsic :: real, max
!..
!.. External Functions ..
      Integer :: ioben
      External :: ioben
!..
!.. Local Scalars ..
      Integer :: litercurr, literlast, lnodelast, nprime, q, r


!---  normalize ranges
      literlast = max(0, (iterlast-iterfirst+iterstep)/iterstep)
      litercurr = ((itercurr-iterfirst)/iterstep) + 1
      lnodelast = nodelast - nodefirst + 1


!---  calculate blocksize
      q = literlast/lnodelast

!---    calculate iteration mapping
      If (q*lnodelast==literlast) Then

!---       equally sized blocks
        mapblock = (nodefirst-1) + ioben(real(litercurr,kind=dp)/real(q,kind= &
          dp))

      Else

!---       unequal blocks
        r = literlast - (q*lnodelast)

!---       up to nprime blocks of size (q+1)
        nprime = (q+1)*r

        If (litercurr<=nprime) Then

          mapblock = (nodefirst-1) + ioben(real(litercurr,kind=dp)/real(q+1, &
            kind=dp))

        Else

          mapblock = (nodefirst-1) + ioben(real(litercurr-nprime,kind=dp)/real &
            (q,kind=dp)+real(r,kind=dp))
        End If

      End If

      If ((mapblock<nodefirst) .Or. (mapblock>nodelast)) Then
        Write (6, Fmt=*) 'internal error in mapblock'
        Write (6, Fmt=*) mapblock, itercurr, iterfirst, iterlast, iterstep, &
          nodefirst, nodelast, literlast, litercurr, lnodelast, q, r, nprime
        Stop
      End If

      Return
    End Function
