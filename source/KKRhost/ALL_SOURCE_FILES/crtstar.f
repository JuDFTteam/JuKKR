c 20.07.96 ***************************************************************
      SUBROUTINE CRTSTAR(RATOM,NSHELL,ND,IROT,ISYMINDEX,RROT)
      implicit none
c ************************************************************************
C  THE SYMMETRY OPERATIONS OF THE SYMMETRY GROUP ARE APPLIED TO THE
C  INPUT VECTOR RATOM
c ------------------------------------------------------------------------
      INTEGER IROT,NSHELL
      DOUBLE PRECISION ND(64,3,*),RATOM(3,*),RROT(48,3,*)
      INTEGER ISYMINDEX(*)

      INTEGER I,ID,NS,K,J,ISYM
      LOGICAL TEST
      EXTERNAL TEST
c ------------------------------------------------------------------------

      DO 80 NS = 1,NSHELL
        DO 70 ID = 1,IROT
          ISYM = ISYMINDEX(ID)
          DO 60 I = 1,3
            RROT(ID,I,NS) = ND(ISYM,I,1)*RATOM(1,NS) +
     +                      ND(ISYM,I,2)*RATOM(2,NS) + 
     +                      ND(ISYM,I,3)*RATOM(3,NS)
 60       CONTINUE
 70     CONTINUE
 80   CONTINUE                      ! NS = 1,NSHELL

      IF (TEST('ND      ')) 
     +     WRITE(1337,FMT='((I3,3(/,3f6.2)))') 
     +     (k,((ND(K,I,J),J=1,3), I=1,3),k=1,irot)

      RETURN

      END
