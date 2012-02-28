      INTEGER FUNCTION MAPBLOCK(ITERCURR,ITERFIRST,ITERLAST,ITERSTEP,
     +                          NODEFIRST,NODELAST)
c-----------------------------------------------------------------------
c
c     Maps loop iteration to nodes (once a call).
c     The mapping is done in blocks, such that:
c       1) all nodes work on q consecutive iterations
c          (equally sized blocks)
c       2) the first r nodes get (q+1) iterations and (n-r) nodes get
c          q iterations
c
c     Descrition of input parameters:
c
c       itercurr  : index of current iteration
c       iterfirst : index of first iteration
c       iterlast  : index of last iteration
c       iterstep  : iteration step
c       nodefirst : index of first node
c       nodelast  : index of last node
c
c     where:
c
c       iterfirst <= itercurr <= iterlast
c       nodefirst <= nodelast
c       iterationstep > 0
c
c
c     Example for Intel NX-calls:
c
c      DO i=iterstart, iterend, iterstep
c
c         node = mapblock(i,iterstart,iterend,iterstep,0,numnodes()-1)
c
c         IF (mynode() .EQ. node) THEN
c                  ! the code for an iteration
c         ENDIF
c
c      ENDDO
c
c
c     Algorithm:
c
c              |  number_iterations  |
c     let q := |  -----------------  |
c              |    number_nodes     |
c              ---                 ---
c
c     number_nodes = q * number_nodes + r
c                  = (q + 1) * r + q * (number_nodes - r)
c
c     IF q = number_iterations / number_nodes
c       THEN
c         you have equally sized blocks
c       ELSE
c         you have r blocks of size q+1
c         and (number_nodes - r) blocks of size q
c
c
c                                           Rudolf Berrendorf, July 1992
c                                           last update: February 1994
c-----------------------------------------------------------------------
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER ITERCURR,ITERFIRST,ITERLAST,ITERSTEP,NODEFIRST,NODELAST
C     ..
C     .. Intrinsic Functions ..

      INTRINSIC DBLE,MAX
C     ..
C     .. Local Scalars ..
      INTEGER LITERCURR,LITERLAST,LNODELAST,NPRIME,Q,R
C     ..


C---  normalize ranges
      LITERLAST = MAX(0, (ITERLAST-ITERFIRST+ITERSTEP)/ITERSTEP)
      LITERCURR = ((ITERCURR-ITERFIRST)/ITERSTEP) + 1
      LNODELAST = NODELAST - NODEFIRST + 1


C---  calculate blocksize
      Q = LITERLAST/LNODELAST

C---    calculate iteration mapping
      IF (Q*LNODELAST.EQ.LITERLAST) THEN

C---       equally sized blocks
        MAPBLOCK = (NODEFIRST-1) + CEILING(DBLE(LITERCURR)/DBLE(Q))

      ELSE

C---       unequal blocks
        R = LITERLAST - (Q*LNODELAST)

C---       up to nprime blocks of size (q+1)
        NPRIME = (Q+1)*R

        IF (LITERCURR.LE.NPRIME) THEN

          MAPBLOCK = (NODEFIRST-1) + CEILING(DBLE(LITERCURR)/DBLE(Q+1))

        ELSE

          MAPBLOCK = (NODEFIRST-1) + CEILING(DBLE(LITERCURR-NPRIME)/
     +               DBLE(Q)+DBLE(R))
        END IF

      END IF

      IF ((MAPBLOCK.LT.NODEFIRST) .OR. (MAPBLOCK.GT.NODELAST)) THEN
        WRITE (6,FMT=*) 'internal error in mapblock'
        WRITE (6,FMT=*) MAPBLOCK,ITERCURR,ITERFIRST,ITERLAST,ITERSTEP,
     +    NODEFIRST,NODELAST,LITERLAST,LITERCURR,LNODELAST,Q,R,NPRIME
        STOP
      END IF

      RETURN
      END
