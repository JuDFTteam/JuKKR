      SUBROUTINE VMADELBLK(CMOM,CMINST,LMAX,NSPIN,NAEZ,VONS,ZAT,R,
     &                     IRCUT,IPAN,VMAD,
     &                     LMPOT,SMAT,CLEB,ICLEB,IEND,
     &                     LMXSPD,NCLEBD,LOFLM,DFAC,I1,
     >                     LMPIC,MYLRANK,
     >                     LGROUP,LCOMM,LSIZE)
C **********************************************************************
C
C     calculate the madelung potentials and add these to the poten-
C     tial v  (in the spin-polarized case for each spin-direction
C     this is the same)
C     it uses the structure dependent matrices AVMAD and BVMAD which
C     are calculated once in the subroutine MADELUNG3D and saved in
C     the DA-file abvmad.unformatted                     ( may 2004)
C     the charge-moments are calculated in the subroutine vintras,
C     therefore vintras has to be called first.
C     the madelung-potential is expanded into spherical harmonics.
C     the lm-term of the potential v of the atom i is given by
C
C      v(r,lm,i) =  (-r)**l * {avmad(i,i2,lm,l'm')*cmom(i2,l'm')
C                                               +bvmad(i,i2,lm)*z(i2)}
C
C     summed over i2 (all atoms) and l'm'
C     (see notes by b.drittler)
C
C                               b.drittler   nov. 1989
C
C     adopted for the case of more atoms on the same site, summation is
C     done over the occupants of that site, the charge is weighted with
C     the appropriate concentration of the occupant  v.popescu feb. 2002
C     
C     impurity-program adopted feb. 2004 (according to n.papanikalou)
C     
C **********************************************************************
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'inc.p'
C     ..
C     .. PARAMETER definitions
      INTEGER LMPOTD
      PARAMETER (LMPOTD=(LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IEND,LMAX,LMXSPD,NCLEBD,LMPOT,NSPIN,NAEZ
      DOUBLE PRECISION VMAD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CMOM(LMPOTD),CMINST(LMPOTD),
     &                 CMOM_SAVE(LMPOTD),CMINST_SAVE(LMPOTD),
     &                 VONS(IRMD,LMPOTD,2),R(IRMD,*),ZAT(*)
      DOUBLE PRECISION SMAT(LMXSPD,*),CLEB(NCLEBD)
      DOUBLE PRECISION DFAC(0:LPOTD,0:LPOTD)
      INTEGER ICLEB(NCLEBD,3),LOFLM(LMXSPD)
      INTEGER IRCUT(0:IPAND,*),IPAN(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC(LMPOTD),PI,FPI
      INTEGER I,L,L1,L2,LM,LM1,LM2,LM3,LMMAX,M,IPOT,I1,I2,
     &        IRS1,ISPIN,ILM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD)
      DOUBLE PRECISION BVMAD(LMPOTD)
C----- MPI ---------------------------------------------------------------
C     .. N-MPI .. 
      INTEGER MYRANK,NROFNODES
      COMMON /MPI/MYRANK,NROFNODES
      INTEGER IERR,MAPBLOCK
C     .. L-MPI
      INTEGER      MYLRANK(LMPID*SMPID*EMPID),
     +             LCOMM(LMPID*SMPID*EMPID),
     +             LGROUP(LMPID*SMPID*EMPID),
     +             LSIZE(LMPID*SMPID*EMPID),
     +             LMPI,LMPIC
C
      EXTERNAL MPI_BCAST
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..................................................................
C
      PI = 4.0D0*ATAN(1.0D0) 
      FPI = 4.0D0*PI
C
      LMMAX = (LMAX+1)*(LMAX+1)
         IRS1 = IRCUT(IPAN(I1),I1)

C===== save information on each proc/atom in order to bcast ==============
      
      DO ILM = 1, LMPOT
        CMOM_SAVE(ILM) = CMOM(ILM)
        CMINST_SAVE(ILM) = CMINST(ILM)
      ENDDO

C===== begin loop over all atoms =========================================

      DO I2 = 1,NAEZ

C-------- bcast information of cmom and cminst begin --------------------

         IF (MYLRANK(LMPIC).EQ.
     +       MAPBLOCK(I2,1,NAEZ,1,0,LSIZE(LMPIC)-1)) THEN
         DO ILM = 1, LMPOT
           CMOM(ILM) = CMOM_SAVE(ILM)
           CMINST(ILM) = CMINST_SAVE(ILM)
         ENDDO
         ENDIF

         CALL MPI_BCAST(CMOM,LMPOTD,MPI_DOUBLE_PRECISION,
     +                  MAPBLOCK(I2,1,NAEZ,1,0,LSIZE(LMPIC)-1),
     +                  LCOMM(LMPIC),IERR)
         CALL MPI_BCAST(CMINST,LMPOTD,MPI_DOUBLE_PRECISION,
     +                  MAPBLOCK(I2,1,NAEZ,1,0,LSIZE(LMPIC)-1),
     +                  LCOMM(LMPIC),IERR)
 
C-------- bcast information of cmom and cminst end ----------------------

C
C --> calculate avmad(lm1,lm2)
C
      DO LM1 = 1,LMPOT
         DO LM2 = 1,LMPOT
            AVMAD(LM1,LM2) = 0.0D0
         END DO
      END DO
      DO I = 1,IEND
         LM1 = ICLEB(I,1)
         LM2 = ICLEB(I,2)
         LM3 = ICLEB(I,3)
         L1 = LOFLM(LM1)
         L2 = LOFLM(LM2)
C
C --> this loop has to be calculated only for l1+l2=l3
C
         AVMAD(LM1,LM2) = AVMAD(LM1,LM2) + 
     &                 2.0D0*DFAC(L1,L2)*SMAT(LM3,I2)*CLEB(I)
      END DO
C
C
C --> calculate bvmad(lm1)
C
      DO LM1 = 1,LMPOT
         BVMAD(LM1) = 0.0D0
      END DO
      IF(NAEZ.GT.1) THEN
C---> lm = 1 component disappears if there is only one host atom
      DO LM1 = 1,LMPOT
         L1 = LOFLM(LM1)
         BVMAD(LM1) = BVMAD(LM1) -
     &             2.0D0*FPI/DBLE(2*L1+1)*SMAT(LM1,I2)
      END DO
      END IF
C
      DO LM = 1,LMPOT
         IF(I2.EQ.1) AC(LM) = 0.0D0
         AC(LM) = AC(LM) + BVMAD(LM)*ZAT(I2)
C
C---> take moments of sphere
C
         DO LM2 = 1,LMMAX
            AC(LM) = AC(LM) + AVMAD(LM,LM2)*CMOM(LM2)
         END DO
         DO LM2 = 1,LMMAX
            AC(LM) = AC(LM) + AVMAD(LM,LM2)*CMINST(LM2)
         END DO
      END DO

      END DO
C===== end loop over all atoms =========================================
C
C
C
C===== save information on each proc/atom in order to bcast ==============
      
      DO ILM = 1, LMPOT
        CMOM(ILM) = CMOM_SAVE(ILM)
        CMINST(ILM) = CMINST_SAVE(ILM)
      ENDDO

C========================================================================
C
C
C
C
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO L = 0,LMAX
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         DO M = -L,L
            LM = L*L + L + M + 1
C     
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF ( LM.EQ.1 ) VMAD = AC(1)/SQRT(4.D0*PI)
C
C---> add to v the intercell-potential
C
C ================================================================= SPIN
            DO ISPIN = 1,NSPIN
C
C---> in the case of l=0 : r(1)**l is not defined
C
               IF ( L.EQ.0 ) VONS(1,1,ISPIN) = VONS(1,1,ISPIN) + AC(LM)
               DO I = 2,IRS1
                  VONS(I,LM,ISPIN) = VONS(I,LM,ISPIN)
     &                             + (-R(I,I1))**L*AC(LM)
               END DO
            END DO
C ================================================================= SPIN
C 
         END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      RETURN
C
      END
