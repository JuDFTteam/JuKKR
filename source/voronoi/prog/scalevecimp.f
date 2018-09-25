      SUBROUTINE SCALEVECIMP(
     >        NUMIMP,NKILLATOM,BRAVAIS,LINTERFACE,LCARTESIMP,
     X        RIMPURITY,RKILL,DXIMP,DYIMP,DZIMP)
      implicit none
c#@# KKRtags: VORONOI KKRimp geometry
C     Changes the impurity-atom coordinates from "internal"
C     (Bravais-basis) to cartesian.

c Input:
      INTEGER NUMIMP,NKILLATOM
      REAL*8 BRAVAIS(3,3)
      LOGICAL LINTERFACE
c I/O:
      REAL*8 RIMPURITY(3,*),RKILL(3,*),DXIMP(*),DYIMP(*),DZIMP(*)
c Local:
      LOGICAL LCARTESIMP
      INTEGER IX,NX,IAT,IER
      REAL*8 RAUX(3)
      CHARACTER*256 UIO



      IF (LCARTESIMP) RETURN  ! Do nothing if coordinates are already cartesian

c Transform imp. positions to cartesian coordinates.

      NX = 3
      IF (LINTERFACE) NX = 2     ! In 2d-geometry, only in-plane parts
                                  ! can be in internal coordinates
      RAUX(1:3) = 0.D0
      DO IAT = 1,NUMIMP
         RAUX(1:NX) = RIMPURITY(1:NX,IAT)
         DO IX=1,NX
            RIMPURITY(IX,IAT)=  BRAVAIS(IX,1)*RAUX(1) +
     &                          BRAVAIS(IX,2)*RAUX(2) +
     &                          BRAVAIS(IX,3)*RAUX(3)  
         ENDDO
      ENDDO

      DO IAT = 1,NKILLATOM
         RAUX(1:NX) = RKILL(1:NX,IAT)
         DO IX=1,NX
            RKILL(IX,IAT)=  BRAVAIS(IX,1)*RAUX(1) +
     &                      BRAVAIS(IX,2)*RAUX(2) +
     &                      BRAVAIS(IX,3)*RAUX(3)  
         ENDDO
      ENDDO

      IF (.NOT.LINTERFACE) THEN
         DO IAT = 1,NUMIMP
            RAUX(1) = DXIMP(IAT)
            RAUX(2) = DYIMP(IAT)
            RAUX(3) = DZIMP(IAT)
            DXIMP(IAT) = BRAVAIS(1,1) * RAUX(1) +
     &                   BRAVAIS(1,2) * RAUX(2) +
     &                   BRAVAIS(1,3) * RAUX(3) 
            DYIMP(IAT) = BRAVAIS(2,1) * RAUX(1) +
     &                   BRAVAIS(2,2) * RAUX(2) +
     &                   BRAVAIS(2,3) * RAUX(3) 
            DZIMP(IAT) = BRAVAIS(3,1) * RAUX(1) +
     &                   BRAVAIS(3,2) * RAUX(2) +
     &                   BRAVAIS(3,3) * RAUX(3) 
         ENDDO

      ELSE

         DO IAT = 1,NUMIMP
            RAUX(1) = DXIMP(IAT)
            RAUX(2) = DYIMP(IAT)
            DXIMP(IAT) = BRAVAIS(1,1) * RAUX(1) +
     &                   BRAVAIS(1,2) * RAUX(2) 
            DYIMP(IAT) = BRAVAIS(2,1) * RAUX(1) +
     &                   BRAVAIS(2,2) * RAUX(2) 
         ENDDO

      ENDIF



      RETURN
      END

