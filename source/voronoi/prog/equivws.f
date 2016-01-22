c***********************************************************************
      LOGICAL FUNCTION EQUIVWS(NFACE1,A1,B1,C1,D1,NFACE2,A2,B2,C2,D2)
c Given two WS- (or Voronoi-) cells, with NFACE1 faces defined by
c A1*x+B1*y+C1*z=D1, for the first, and correspondingly for the second,
c both centered at the origin, this function takes the value .TRUE. if 
c the cells coincide fully. Else, the value .FALSE. is returned.  It is
c naturally assumed that, if the faces councide one-by-one, so will the
c edges. The latter can be searched faster, and thus are used here. 
c Since the coefficients that define the faces are arbitary up to a
c multiplicative constant, they are first normalized to D1=1., D2=1.. 
c (The case D1=0. or D2=0. should not occur, since WS-cell faces never 
c contain the origin. This should be checked before calling this 
c function, so that the check is done only once for each WS-cell).
      implicit none
c Input:
      INTEGER NFACE1,NFACE2
      REAL*8           A1(*),B1(*),C1(*),D1(*),A2(*),B2(*),C2(*),D2(*)
c Inside:
      INTEGER IFACE1,IFACE2
      REAL*8           AA1,BB1,CC1,AA2,BB2,CC2 ! Temporary plane coefficients.
      LOGICAL LTEMP

c---------------------------------------------------------------
c Initialize:
      EQUIVWS = .FALSE.
c---------------------------------------------------------------
      IF (NFACE1.NE.NFACE2) GOTO 100

c---------------------------------------------------------------
c Loop over all faces of the first cell:
      DO IFACE1 = 1,NFACE1
c Normalize to D1=1:
         AA1 = A1(IFACE1)/D1(IFACE1)
         BB1 = B1(IFACE1)/D1(IFACE1)
         CC1 = C1(IFACE1)/D1(IFACE1)

         LTEMP = .FALSE. ! Becomes .TRUE. if a match for IFACE1 is found
         DO IFACE2 = 1,NFACE2

            AA2 = A2(IFACE2)/D1(IFACE2)
            IF (DABS(AA2-AA1).LT.1.D-8) THEN
               BB2 = B2(IFACE2)/D1(IFACE2)
               IF (DABS(BB2-BB1).LT.1.D-8) THEN
                  CC2 = C2(IFACE2)/D1(IFACE2)
                  IF (DABS(CC2-CC1).LT.1.D-8) THEN
                     LTEMP = .TRUE. ! Found a match for IFACE1
                     GOTO 50        ! Jump out of loop
                  ENDIF
               ENDIF
            ENDIF

         ENDDO

 50      CONTINUE               ! Jumped out of loop at GOTO 50 earlier
         IF (.NOT.LTEMP) GOTO 100 ! Did not find any match, for IFACE1

      ENDDO
      EQUIVWS = .TRUE. ! Loops were completed, found mach for all IFACE2
c---------------------------------------------------------------


 100  CONTINUE                  ! From two GOTO 100 statements earlier
      RETURN
      END
