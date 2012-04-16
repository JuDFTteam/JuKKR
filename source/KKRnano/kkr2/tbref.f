C-----------------------------------------------------------------------
C>    Calculate reference system's T-matrix.
C>    @param     E complex energy
C>    @param     VREF repulsive reference potential field strength
C>    @param     LMAX angular momentum cutoff
C>    @param     RMTREF repulsive reference potential muffin-tin radius
C>    @param     TREFLL reference system T-matrix
C>    @param     DTREFLL energy derivative of reference system T-matrix
C>               only calculated if LLY=1
C>    @param     LLY do Lloyd's formula calculation 0=no/1=yes
      SUBROUTINE TREF(E,VREF,LMAX,RMTREF,TREFLL,DTREFLL,
C                     new input parameter after inc.p removal
     &                LLY)
C
      IMPLICIT NONE

      INTEGER LLY

C     INTEGER LMGF0D
C     PARAMETER (LMGF0D= (LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION RMTREF,VREF
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX A1,B1,DA1,DB1,E
      INTEGER I,L,LM1
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX BESSJW1(0:LMAX+1),BESSJW2(0:LMAX+1),
     +               BESSYW1(0:LMAX+1),BESSYW2(0:LMAX+1),
     +               HANKWS1(0:LMAX+1),HANKWS2(0:LMAX+1)
C     DOUBLE COMPLEX TREFLL(LMGF0D,LMGF0D)
      DOUBLE COMPLEX TREFLL((LMAX+1)**2, (LMAX+1)**2)

      DOUBLE COMPLEX DBESSJW1(0:LMAX+1),DBESSJW2(0:LMAX+1),
     +               DHANKWS1(0:LMAX+1)
C     DOUBLE COMPLEX DTREFLL(LMGF0D,LMGF0D)
      DOUBLE COMPLEX DTREFLL((LMAX+1)**2, (LMAX+1)**2)

C     .. External Subroutines ..
      EXTERNAL BESSEL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Save statement ..
C     SAVE
      INTEGER LMAXD
      INTEGER LMGF0D

C-----------------------------------------------------------------------
C---- t-matrix and derivative of t-matrix of the reference system
C
C     the analytical formula for the derivative of spherical Bessel
C     functions is used:
C
C     d                     l+1
C     --  j (x) = j   (x) - --- j (x)   
C     dx   l       l-1       x   l
C
C     d
C     --  j (x) = - j (x)
C     dx   0         1
C
C     which for x = sqrt(E0)*r leads to
C
C      d          r*r             (l+1)
C     --- j (x) = --- ( j   (x) - ----- j (x) )
C     dE0  l      2 x    l-1        x    l
C
C      d            r*r
C     --- j (x) = - --- j (x)
C     dE0  0        2 x  1
C
C-----------------------------------------------------------------------

      LMAXD = LMAX
      LMGF0D= (LMAXD+1)**2

          A1 = SQRT(E)*RMTREF
          B1 = SQRT(E-VREF)*RMTREF
          CALL BESSEL(BESSJW1,BESSYW1,HANKWS1,A1,LMAXD+1)
          CALL BESSEL(BESSJW2,BESSYW2,HANKWS2,B1,LMAXD+1)

      IF(LLY.EQ.1) THEN

          DBESSJW1(0) = - BESSJW1(1)/A1
          DBESSJW2(0) = - BESSJW2(1)/B1
          DHANKWS1(0) = - HANKWS1(1)/A1

          DO L = 1,LMAX + 1
          DBESSJW1(L) = (BESSJW1(L-1) - (L+1)*BESSJW1(L)/A1)/A1
          DBESSJW2(L) = (BESSJW2(L-1) - (L+1)*BESSJW2(L)/B1)/B1
          DHANKWS1(L) = (HANKWS1(L-1) - (L+1)*HANKWS1(L)/A1)/A1
          END DO

          DO L = 0,LMAX + 1
          DBESSJW1(L) = 0.5D0*DBESSJW1(L)*RMTREF**2 
          DBESSJW2(L) = 0.5D0*DBESSJW2(L)*RMTREF**2 
          DHANKWS1(L) = 0.5D0*DHANKWS1(L)*RMTREF**2 
          END DO

      ENDIF

          DO L = 0,LMAX
            A1 = SQRT(E)*BESSJW1(L+1)*BESSJW2(L) -
     +           SQRT(E-VREF)*BESSJW1(L)*BESSJW2(L+1)

            B1 = SQRT(E)*HANKWS1(L+1)*BESSJW2(L) -
     +           SQRT(E-VREF)*HANKWS1(L)*BESSJW2(L+1)

            DO I = -L,L
              LM1 = L* (L+1) + I + 1
              TREFLL(LM1,LM1) = -1.D0/SQRT(E)*A1/B1
            END DO

          END DO

C
      CALL CINIT(LMGF0D*LMGF0D,DTREFLL)
C
      IF(LLY.EQ.1) THEN
C
          DO L = 0,LMAX
            A1 = SQRT(E)*BESSJW1(L+1)*BESSJW2(L) -
     +           SQRT(E-VREF)*BESSJW1(L)*BESSJW2(L+1)
C
            DA1 = 0.5D0/SQRT(E)*BESSJW1(L+1)*BESSJW2(L) -
     +            0.5D0/SQRT(E-VREF)*BESSJW1(L)*BESSJW2(L+1) +
     +            SQRT(E)*DBESSJW1(L+1)*BESSJW2(L) -
     +            SQRT(E-VREF)*DBESSJW1(L)*BESSJW2(L+1) +
     +            SQRT(E)*BESSJW1(L+1)*DBESSJW2(L) -
     +            SQRT(E-VREF)*BESSJW1(L)*DBESSJW2(L+1)
C
            B1 = SQRT(E)*HANKWS1(L+1)*BESSJW2(L) -
     +           SQRT(E-VREF)*HANKWS1(L)*BESSJW2(L+1)
C
            DB1 = 0.5D0/SQRT(E)*HANKWS1(L+1)*BESSJW2(L) -
     +            0.5D0/SQRT(E-VREF)*HANKWS1(L)*BESSJW2(L+1) +
     +            SQRT(E)*DHANKWS1(L+1)*BESSJW2(L) -
     +            SQRT(E-VREF)*DHANKWS1(L)*BESSJW2(L+1) +
     +            SQRT(E)*HANKWS1(L+1)*DBESSJW2(L) -
     +            SQRT(E-VREF)*HANKWS1(L)*DBESSJW2(L+1)
C
            DO I = -L,L
              LM1 = L* (L+1) + I + 1
              DTREFLL(LM1,LM1) = 0.5D0/SQRT(E)**3*A1/B1
     +                         - 1.D0/SQRT(E)*(DA1/B1-A1*DB1/B1**2)
            END DO

          END DO
        ENDIF
C
      END



      SUBROUTINE GREF(E,ALATC,IEND,NCLS,NAEZ,
     +                CLEB,RCLS,ATOM,CLS,ICLEB,LOFLM,NACLS,
     +                REFPOT,
     +                TREFLL,DTREFLL,GREFN,DGREFN,
     +                IE,
     +                LLY_G0TR,
     +                LMPIC,MYLRANK,LCOMM,LSIZE,
     &                naezd, lmaxd, naclsd, ncleb, nrefd, iemxd, nclsd,
     &                LLY, prod_lmpid_smpid_empid)
C
      IMPLICIT NONE
C
C     .. Parameters ..
      INCLUDE 'mpif.h'
C     INCLUDE 'inc.p'
C     INCLUDE 'inc.cls'
C
      integer  naezd
      integer  lmaxd
      integer  naclsd
      integer  ncleb
      integer  nrefd
      integer  iemxd
      integer  nclsd
      integer  LLY
      integer  prod_lmpid_smpid_empid

C     INTEGER          LMGF0D
C     PARAMETER       (LMGF0D= (LMAXD+1)**2)
C     INTEGER          LMMAXD
C     PARAMETER       (LMMAXD= (LMAXD+1)**2)
C     INTEGER          LM2D
C     PARAMETER       (LM2D= (2*LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALATC
      INTEGER          IE,IEND,NCLS,NAEZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(NCLEB,2),RCLS(3,NACLSD,NCLSD)
      INTEGER          ATOM(NACLSD,NAEZD),CLS(NAEZD),ICLEB(NCLEB,3)
C     INTEGER          LOFLM(LM2D)
      INTEGER          LOFLM((2*LMAXD+1)**2)
      INTEGER          NACLS(NCLSD)
      INTEGER          REFPOT(NAEZD)
      DOUBLE COMPLEX   LLY_G0TR(IEMXD,NCLSD)

C     DOUBLE COMPLEX   TREFLL(LMGF0D,LMGF0D,NREFD),
C    +                 DTREFLL(LMGF0D,LMGF0D,NREFD),
C    +                 DGREFN(LMGF0D,LMGF0D,NACLSD,NCLSD),
C    +                 GREFN(LMGF0D,LMGF0D,NACLSD,NCLSD)

      DOUBLE COMPLEX   TREFLL((LMAXD+1)**2,(LMAXD+1)**2,NREFD),
     +                 DTREFLL((LMAXD+1)**2,(LMAXD+1)**2,NREFD)
      DOUBLE COMPLEX   DGREFN((LMAXD+1)**2,(LMAXD+1)**2,NACLSD,NCLSD),
     +                 GREFN((LMAXD+1)**2,(LMAXD+1)**2,NACLSD,NCLSD)

C     .. L-MPI
      INTEGER          MYLRANK(prod_lmpid_smpid_empid),
     +                 LCOMM(prod_lmpid_smpid_empid),
     +                 LSIZE(prod_lmpid_smpid_empid),
     +                 LMPIC

C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   E
      INTEGER          I1,IC,ICLS,IG,IG1,LM1,LM2
C     ..

C     Local arrays
C     Fortran 90 automatic arrays - medium size

C     DOUBLE COMPLEX   DGINP(NACLSD*LMGF0D,LMGF0D),
C    +                 GINP(NACLSD*LMGF0D,LMGF0D),
C    +                 GBCAST(LMMAXD,LMMAXD,NACLSD)

      DOUBLE COMPLEX   DGINP(NACLSD*(LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX   GINP(NACLSD*(LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX   GBCAST((LMAXD+1)**2,(LMAXD+1)**2,NACLSD)
C     ..
C     .. External Subroutines ..
      EXTERNAL         GLL95

      INTEGER          IERR,MAPBLOCK

C
      EXTERNAL         MPI_BCAST
C
C     .. Save statement ..
C     SAVE
C     ..

      INTEGER          LMGF0D
      INTEGER          LMMAXD

      LMMAXD= (LMAXD+1)**2
      LMGF0D= (LMAXD+1)**2

C attention in this subroutine I3 labels the fixed atom - I1 is a variable !
C
C=====================================================================
        DO ICLS = 1,NCLS
C=====================================================================
C
C NCLS can by no means be larger than NAEZ there distribute as follows
C note that the parallelization of this routine is only active if there
C are more than four non-identical reference clusteres
C
          IF (MYLRANK(LMPIC).EQ.
     +    MAPBLOCK(ICLS,1,NAEZ,1,0,LSIZE(LMPIC)-1).OR.NCLS.LT.5) THEN
C
          I1 = 1
          IC = 0
          DO WHILE (IC.EQ.0 .AND. I1.LE.NAEZ)
            IF (CLS(I1).EQ.ICLS) IC = I1
            I1 = I1 + 1
          END DO
          IF (IC.EQ.0) STOP 'Error in CLS(*) array in tbref'

          CALL GLL95(E,CLEB(1,2),ICLEB,LOFLM,IEND,TREFLL,DTREFLL,
     +               ATOM(1,IC),REFPOT,RCLS(1,1,ICLS),NACLS(ICLS),
     +               ALATC,GINP,DGINP,
     +               LLY_G0TR(IE,ICLS),
     &               lmaxd, naclsd, ncleb, nrefd, LLY )

          DO IG=1,NACLSD
          DO LM2=1,LMGF0D
            IG1 = (IG-1)*LMGF0D + LM2 
            DO LM1=1,LMGF0D
               GREFN(LM2,LM1,IG,ICLS)=GINP(IG1,LM1)
               DGREFN(LM2,LM1,IG,ICLS)=DGINP(IG1,LM1)
            ENDDO
          ENDDO
          ENDDO
C
          ENDIF
C
C=====================================================================
        END DO
C=====================================================================
C
C ok, now MPI_BCAST the results to all processors
        IF (NCLS.GT.4) THEN
C
C=====================================================================
          DO ICLS=1, NCLS
C=====================================================================
C 1st broadcast reference structure constants
C
            DO IG=1, NACLSD
              DO LM2=1, LMMAXD
                DO LM1=1, LMMAXD
                  GBCAST(LM1,LM2,IG) = GREFN(LM1,LM2,IG,ICLS)
                ENDDO
              ENDDO
            ENDDO
C
            CALL MPI_BCAST(GBCAST,LMMAXD*LMMAXD*NACLSD,
     +                   MPI_DOUBLE_COMPLEX,
     +                   MAPBLOCK(ICLS,1,NAEZ,1,0,LSIZE(LMPIC)-1),
     +                   LCOMM(LMPIC),IERR)
C     
            CALL MPI_BARRIER(LCOMM(LMPIC),IERR)
C
            DO IG=1, NACLSD
              DO LM2=1, LMMAXD
                DO LM1=1, LMMAXD
                  GREFN(LM1,LM2,IG,ICLS) = GBCAST(LM1,LM2,IG)
                ENDDO
              ENDDO
            ENDDO
C
C
C 2nd and if Lloyd's formula is going to be applied broadcast
C derivative of reference structure constants
C FIXME: Code duplication: same as above but with GREFN -> DGREFN
C        -check! -extract to separate routine
C
            IF (LLY.EQ.1) THEN
C
            DO IG=1, NACLSD
              DO LM2=1, LMMAXD
                DO LM1=1, LMMAXD
                  GBCAST(LM1,LM2,IG) = DGREFN(LM1,LM2,IG,ICLS)
                ENDDO
              ENDDO
            ENDDO
C
            CALL MPI_BCAST(GBCAST,LMMAXD*LMMAXD*NACLSD,
     +                   MPI_DOUBLE_COMPLEX,
     +                   MAPBLOCK(ICLS,1,NAEZ,1,0,LSIZE(LMPIC)-1),
     +                   LCOMM(LMPIC),IERR)
C     
            CALL MPI_BARRIER(LCOMM(LMPIC),IERR)
C
            DO IG=1, NACLSD
              DO LM2=1, LMMAXD
                DO LM1=1, LMMAXD
                  DGREFN(LM1,LM2,IG,ICLS) = GBCAST(LM1,LM2,IG)
                ENDDO
              ENDDO
            ENDDO
C
            ENDIF
C
C=====================================================================
          ENDDO
C=====================================================================
C
        ENDIF
C
      END
