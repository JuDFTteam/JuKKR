c***********************************************************************
      SUBROUTINE WSCLASSES(NFACEMAX,NFACE,NCELL,A3,B3,C3,D3,
     &                     NCLASS,CLASS,CLASSREP)
c Given a number of WS-cells, defined by their faces A3*x+B3*y+C3*z=D3
c in cell-centered coordinates, this subroutine sorts them into
c equivalence classes. Two WS-cells belong to the same class, if they
c are characterized by the same faces. The total number of classes
c found is NCLASS. Each cell, indexed ICELL, belongs to some class,
c namely CLASS(ICELL). Each class, indexed ICLASS, is represented by
c the cell CLASSREP(ICLASS).
c
c Uses logical function EQUIVWS.
      implicit none
c Input:
      INTEGER NFACEMAX   ! Max. number of faces per cell (dimension)
      INTEGER NFACE(*)   ! Number of faces per cell (indexed ICELL)
      INTEGER NCELL      ! Number of cells to be investigated.
      REAL*8 A3(NFACEMAX,*),B3(NFACEMAX,*) !First index is for the face,
      REAL*8 C3(NFACEMAX,*),D3(NFACEMAX,*) !second for the cell.
c Output:
      INTEGER NCLASS      ! Number of different classes found.
      INTEGER CLASS(*)    ! Class to which each cell belongs.
      INTEGER CLASSREP(*) ! Pointing to the cell that represents
c                           ! a certain class.
c Inside:
      INTEGER ICELL,ICELL2,ICLASS
      LOGICAL EQUIVWS       ! Function used

c Initialize:
      NCLASS = 0

c Loop over all cells:
      DO 100 ICELL = 1,NCELL

c Compare current cell with all existing classes:
         DO ICLASS = 1,NCLASS
            ICELL2 = CLASSREP(ICLASS)
            IF (EQUIVWS(
     &           NFACE(ICELL),
     &           A3(1,ICELL),B3(1,ICELL),C3(1,ICELL),D3(1,ICELL),
     &           NFACE(ICELL2),
     &           A3(1,ICELL2),B3(1,ICELL2),C3(1,ICELL2),D3(1,ICELL2)))
     &           THEN
c Current cell belongs to this class:
               CLASS(ICELL) = ICLASS ! Place it there...
               GOTO 100      ! ...and step out of the loop over classes.
            ENDIF
         ENDDO    ! ICLASS = 1,NCLASS

c The loop over existing classes has ended with no success. This means
c that the current cell belongs to a new class. It shall also be its
c representative.
         NCLASS = NCLASS + 1
         CLASS(ICELL) = NCLASS
         CLASSREP(NCLASS) = ICELL

 100  CONTINUE

      RETURN
      END

