INTEGER FUNCTION mapblockold(itercurr,iterfirst,iterlast,iterstep,  &
        nodefirst,nodelast)
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
      INTEGER ITERCURR,ITERFIRST,ITERLAST,ITERSTEP,NODEFIRST,NODELAST
!..
!.. Intrinsic Functions ..

      INTRINSIC DBLE,MAX
!..
!.. External Functions ..
      INTEGER IOBEN
      EXTERNAL IOBEN
!..
!.. Local Scalars ..
      INTEGER LITERCURR,LITERLAST,LNODELAST,NPRIME,Q,R


!---  normalize ranges
literlast = MAX(0, (iterlast-iterfirst+iterstep)/iterstep)
litercurr = ((itercurr-iterfirst)/iterstep) + 1
lnodelast = nodelast - nodefirst + 1


!---  calculate blocksize
q = literlast/lnodelast

!---    calculate iteration mapping
IF (q*lnodelast == literlast) THEN
  
!---       equally sized blocks
  mapblock = (nodefirst-1) + ioben(DBLE(litercurr)/DBLE(q))
  
ELSE
  
!---       unequal blocks
  r = literlast - (q*lnodelast)
  
!---       up to nprime blocks of size (q+1)
  nprime = (q+1)*r
  
  IF (litercurr <= nprime) THEN
    
    mapblock = (nodefirst-1) + ioben(DBLE(litercurr)/DBLE(q+1))
    
  ELSE
    
    mapblock = (nodefirst-1) + ioben(DBLE(litercurr-nprime)/ DBLE(q)+DBLE(r))
  END IF
  
END IF

IF ((mapblock < nodefirst) .OR. (mapblock > nodelast)) THEN
  WRITE (6,FMT=*) 'internal error in mapblock'
  WRITE (6,FMT=*) mapblock,itercurr,iterfirst,iterlast,iterstep,  &
      nodefirst,nodelast,literlast,litercurr,lnodelast,q,r,nprime
  STOP
END IF

RETURN
END FUNCTION mapblockold
