      SUBROUTINE EPATHTB(EZ,DF,EFERMI,NPNT,IESEMICORE,IDOSEMICORE,
     &                   EBOTVAL,EMUVAL,TKVAL,NPOLVAL,N1VAL,N2VAL,N3VAL,
     &                   EBOTSEM,EMUSEM,TKSEM,NPOLSEM,N1SEM,N2SEM,N3SEM,
     &                   IEMXD)
C **********************************************************************
C *                                                                    *
C * Generating the energy mesh.                                        *
C *                                                                    *
C * Calls the routine EMESHT once for the valence contour and once for *
C * the semicore contour.                                              *
C * In the semicore range, -NPOLSEM is used to create a rectangular    *
C * contour.                                                           *
C *              ph. mavropoulos, v.popescu Juelich/Munich 2004        *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
      INTEGER IEMXD
      DOUBLE COMPLEX EZ(*),DF(*),EZSEMI(IEMXD),DFSEMI(IEMXD)
      DOUBLE COMPLEX EZVAL(IEMXD),DFVAL(IEMXD)
      DOUBLE PRECISION EBOTSEM,EMUSEM,TKSEM,EBOTVAL,EMUVAL,TKVAL,EFERMI
      INTEGER NPOLSEM,N1SEM,N2SEM,N3SEM
      INTEGER NPOLVAL, N1VAL, N2VAL, N3VAL
      INTEGER IESEMICORE,NPNT,NPNTSEMI,NPNTVAL
      INTEGER IE,JE
      INTEGER IDOSEMICORE

C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,*)
      WRITE (6,'(79(1H=))')
      WRITE (6,'(20X,A)') 'EPATHTB: generates a complex E contour'
      WRITE (6,'(79(1H=))')
      WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      IESEMICORE = 0 
      IF ( IDOSEMICORE.EQ.1 ) THEN
         WRITE(6,99001) 'semi-core contour'
         CALL EMESHT(EZSEMI,DFSEMI,NPNTSEMI,EBOTSEM,EMUSEM,EFERMI,
     &               TKSEM,-NPOLSEM,N1SEM,N2SEM,N3SEM,IEMXD)
         IESEMICORE = NPNTSEMI
         WRITE(6,99001) 'valence contour'
      ENDIF
      CALL EMESHT(EZVAL,DFVAL,NPNTVAL,EBOTVAL,EMUVAL,EFERMI,TKVAL,
     &            NPOLVAL,N1VAL,N2VAL,N3VAL,IEMXD)
C
      NPNT = IESEMICORE + NPNTVAL
C
      DO IE = 1,IESEMICORE
         EZ(IE) = EZSEMI(IE)
         DF(IE) = DFSEMI(IE)
      ENDDO
C
      DO IE = IESEMICORE+1,NPNT
         JE = IE - IESEMICORE
         EZ(IE) = EZVAL(JE)
         DF(IE) = DFVAL(JE)
      ENDDO
C
99001 FORMAT(7X,'* ',A,/,7X,20(1H-),/)
      END
