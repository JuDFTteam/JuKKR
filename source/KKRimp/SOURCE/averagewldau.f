      MODULE MOD_AVERAGEWLDAU
      CONTAINS

      SUBROUTINE AVERAGEWLDAU(NATOM,NSPIN,LDAU)
      USE nrtype
      USE type_ldau
      implicit none
c Input:
      INTEGER NATOM,NSPIN
c Output:
! I/O
      TYPE(LDAU_TYPE)       LDAU(:)
c Local:
      INTEGER I1,IS,LMLO,LMHI,M1,MMAX
      REAL*8 SUM

c     Initialize:
      DO I1 = 1,NATOM
         DO IS = 1,NSPIN
            LDAU(I1)%WLDAUAV(IS) = 0.d0
         ENDDO
      ENDDO

c     Average WLDAU for atoms treated with lda+u:
      DO I1 = 1,NATOM

         IF (LDAU(I1)%LOPT.GE.0) THEN
            LMLO = LDAU(I1)%LOPT**2 + 1
            LMHI = (LDAU(I1)%LOPT+1)**2
            MMAX = LMHI - LMLO + 1

            DO IS = 1,NSPIN
               SUM = 0.D0
               DO M1 = 1,MMAX
                  SUM = SUM + LDAU(I1)%WLDAU(M1,M1,IS)
               ENDDO
               LDAU(I1)%WLDAUAV(IS) = SUM / DFLOAT(MMAX)
            ENDDO

         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE AVERAGEWLDAU

      END MODULE MOD_AVERAGEWLDAU
