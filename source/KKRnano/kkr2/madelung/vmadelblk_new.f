      SUBROUTINE VMADELBLK_new_com(CMOM,CMINST,LPOT,NSPIN,NAEZ,
     &                     VONS,ZAT,R,
     &                     IRCUT,IPAN,VMAD,
     &                     LMPOT,SMAT,CLEB,ICLEB,IEND,
     &                     LMXSPD,NCLEBD,LOFLM,DFAC,
     >                     MYLRANK,
     >                     communicator,comm_size,
     &                     irmd, ipand)
C **********************************************************************
C
C     calculate the madelung potentials and add these to the poten-
C     tial v  (in the spin-polarized case for each spin-direction
C     this is the same)
C     it uses the structure dependent matrices AVMAD and BVMAD which
C     are calculated once in the subroutine MADELUNG3D
C     ( may 2004)
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

      INTEGER irmd
      INTEGER ipand
C     ..
C     .. PARAMETER definitions
C     INTEGER LMPOTD
C     PARAMETER (LMPOTD=(LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IEND,LPOT,LMXSPD,NCLEBD,LMPOT,NSPIN,NAEZ
      DOUBLE PRECISION VMAD
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION CMOM(LMPOTD)
C     DOUBLE PRECISION CMINST(LMPOTD)
C     DOUBLE PRECISION VONS(IRMD,LMPOTD,2)
      DOUBLE PRECISION CMOM((LPOT+1)**2)
      DOUBLE PRECISION CMINST((LPOT+1)**2)
      DOUBLE PRECISION VONS(IRMD,(LPOT+1)**2,2)

      DOUBLE PRECISION R(IRMD)
      DOUBLE PRECISION ZAT(*)
      DOUBLE PRECISION SMAT(LMXSPD,*)
      DOUBLE PRECISION CLEB(NCLEBD)
      DOUBLE PRECISION DFAC(0:LPOT,0:LPOT)
      INTEGER ICLEB(NCLEBD,3)
      INTEGER LOFLM(LMXSPD)
      INTEGER IRCUT(0:IPAND)
      INTEGER IPAN
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,FPI
      INTEGER I,L,L1,L2,LM,LM1,LM2,LM3,LMMAX,M,I2,
     &        IRS1,ISPIN,ILM
C     ..
C     .. Local Arrays ..
C     Fortran 90 automatic arrays
C     DOUBLE PRECISION AC(LMPOTD)
C     DOUBLE PRECISION CMOM_SAVE(LMPOTD)
C     DOUBLE PRECISION CMINST_SAVE(LMPOTD)
C     DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD)
C     DOUBLE PRECISION BVMAD(LMPOTD)

      DOUBLE PRECISION AC((LPOT+1)**2)
      DOUBLE PRECISION CMOM_SAVE((LPOT+1)**2)
      DOUBLE PRECISION CMINST_SAVE((LPOT+1)**2)
      DOUBLE PRECISION AVMAD((LPOT+1)**2,(LPOT+1)**2)
      DOUBLE PRECISION BVMAD((LPOT+1)**2)

C----- MPI ---------------------------------------------------------------
      INTEGER IERR,MAPBLOCK
C     .. L-MPI
      INTEGER      MYLRANK,
     +             communicator,
     +             comm_size
C
      EXTERNAL MPI_BCAST
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT

      INTEGER LMPOTD
      LMPOTD=(LPOT+1)**2
C     ..................................................................
C
      PI = 4.0D0*ATAN(1.0D0) 
      FPI = 4.0D0*PI
C
      LMMAX = (LPOT+1)*(LPOT+1)
         IRS1 = IRCUT(IPAN)

C===== save information on each proc/atom in order to bcast ==============
      
      DO ILM = 1, LMPOT
        CMOM_SAVE(ILM) = CMOM(ILM)
        CMINST_SAVE(ILM) = CMINST(ILM)
      ENDDO

C===== begin loop over all atoms =========================================

      DO I2 = 1,NAEZ

C-------- bcast information of cmom and cminst begin --------------------

         IF (MYLRANK .EQ.
     +       MAPBLOCK(I2,1,NAEZ,1,0,comm_size-1)) THEN
         DO ILM = 1, LMPOT
           CMOM(ILM) = CMOM_SAVE(ILM)
           CMINST(ILM) = CMINST_SAVE(ILM)
         ENDDO
         ENDIF

         CALL MPI_BCAST(CMOM,LMPOTD,MPI_DOUBLE_PRECISION,
     +                  MAPBLOCK(I2,1,NAEZ,1,0,comm_size-1),
     +                  communicator,IERR)
         CALL MPI_BCAST(CMINST,LMPOTD,MPI_DOUBLE_PRECISION,
     +                  MAPBLOCK(I2,1,NAEZ,1,0,comm_size-1),
     +                  communicator,IERR)
 
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
      DO L = 0,LPOT
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
     &                             + (-R(I))**L*AC(LM)
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
