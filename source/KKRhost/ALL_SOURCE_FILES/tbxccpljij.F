      MODULE MOD_TBXCCPLJIJ

      implicit none

       private
       public TBXCCPLJIJ

      contains


      SUBROUTINE TBXCCPLJIJ(IFTMAT,IELAST,EZ,WEZ,NSPIN,NCPA,NAEZ, NATYP,
     &                      NOQ,ITOQ,IQAT,NSHELL,NATOMIMP,ATOMIMP,RATOM,
     &                      FACTL,NOFGIJ,NQCALC,IQCALC,IJTABCALC,
     &                      IJTABSYM,IJTABSH,ISH,JSH,DSYMLL,IPRINT,
     &                      NATYPD,NAEZD,NSHELD,LMMAXD)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-off diagonal  XC-coupling parameters  J_ij  *
C   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
C   *                                                                  *
C   *  adopted for TB-KKR code from Munich SPR-KKR package Sep 2004    *
C   *                                                                  *
C   *  for mpi-parallel version: moved energy loop from main1b into    *
C   *   here.                                B. Zimmermann, Dez 2015   *
C   *                                                                  *
C   ********************************************************************
C
#ifdef CPP_MPI
      use mpi
#endif
      use mod_types, only: t_tgmat, t_mpi_c_grid
      use mod_mympi, only: myrank, master
      IMPLICIT NONE
C     ..
C     .. Parameters
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER ( CONE  = (1D0,0D0) )
      PARAMETER ( CZERO = (0D0,0D0) )
      INTEGER NIJMAX
      PARAMETER ( NIJMAX = 200 )
C     ..
C     .. Scalar arguments
      INTEGER IELAST,NSPIN,IFTMAT,IPRINT,LMMAXD,NAEZ,NAEZD,NATOMIMP,
     &        NATYP,NATYPD,NCPA,NOFGIJ,NQCALC,NSHELD
      DOUBLE COMPLEX EZ(*), WEZ(*)
C     ..
C     .. Array arguments
      INTEGER ATOMIMP(*),IJTABCALC(*),IJTABSH(*),IJTABSYM(*),
     &        IQAT(*),IQCALC(*),ISH(NSHELD,*),ITOQ(NATYPD,*),
     &        JSH(NSHELD,*),NOQ(*),NSHELL(0:NSHELD)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,*),FACTL(LMMAXD,LMMAXD)
      DOUBLE PRECISION RATOM(3,*)
C     ..
C     .. Local scalars
      INTEGER I1,IA,IFGMAT,IFMCPA,IQ,IREC,ISPIN,ISYM,IT,J1,JA,
     &        JQ,JT,L1,LM1,LM2,LSTR,NS,NSEFF,NSHCALC,NSMAX,NTCALC,
     &        IE, ie_end, ie_num, ie_start, IERR
      DOUBLE COMPLEX CSUM, WGTEMP, XINTEGDTMP
      CHARACTER*8 FMT1
      CHARACTER*22 FMT2
      CHARACTER*8 JFBAS
      CHARACTER*10 JFINTEG,JFNAME
      CHARACTER*13 JFNAM
      CHARACTER*26 JFNAM2
      CHARACTER*80 STRBAR,STRTMP
C     ..
C     .. Local arrays
CF90-------------------------------------------------------------
      INTEGER NIJCALC(:),KIJSH(:,:),JIJDONE(:,:,:)
      DOUBLE COMPLEX JXCIJINT(:,:,:),XINTEGD(:,:,:),
     &               csum_store(:,:,:,:), csum_store2(:,:,:,:)
      ALLOCATABLE NIJCALC,JIJDONE,KIJSH,JXCIJINT,XINTEGD,csum_store,
     &            csum_store2
CF90-------------------------------------------------------------
CF77CF77-------------------------------------------------------------
CF77      INTEGER NTLOC,NSHLOC
CF77      PARAMETER (NTLOC = 5, NSHLOC = 40 )
CF77      INTEGER NIJCALC(NSHLOC),KIJSH(NIJMAX,NSHLOC)
CF77      INTEGER JIJDONE(NTLOC,NTLOC,NSHLOC)
CF77      DOUBLE COMPLEX JXCIJINT(NTLOC,NTLOC,NSHLOC)
CF77CF77-------------------------------------------------------------
CF90CF90-------------------------------------------------------------
CF90      INTEGER NIJCALC(:),KIJSH(:,:),JIJDONE(:,:,:)
CF90      DOUBLE COMPLEX JXCIJINT(:,:,:)
CF90      ALLOCATABLE NIJCALC,JIJDONE,KIJSH,JXCIJINT
CF90CF90-------------------------------------------------------------
      DOUBLE COMPLEX DELTSST(LMMAXD,LMMAXD,NATYP),
     &               DMATTS(LMMAXD,LMMAXD,NATYP,2),
     &               DTILTS(LMMAXD,LMMAXD,NATYP,2),GMIJ(LMMAXD,LMMAXD),
     &               GMJI(LMMAXD,LMMAXD),GS(LMMAXD,LMMAXD,2),
     &               TSST(LMMAXD,LMMAXD,NATYP,2),W1(LMMAXD,LMMAXD),
     &               W2(LMMAXD,LMMAXD),W3(LMMAXD,LMMAXD)
      DOUBLE PRECISION RSH(NSHELD),PI
      INTEGER JTAUX(NATYP)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,SQRT
C     ..
C     .. External Subroutines ..
C     EXTERNAL CINIT,CMATMUL,INITABJIJ,ZCOPY
C     ..
C     .. Save statement
      SAVE IFGMAT,IFMCPA,JIJDONE,JXCIJINT,XINTEGD,KIJSH,NIJCALC
      SAVE NSHCALC,NSMAX
C     ..
C     .. Data statement
      DATA JFBAS/'Jij.atom'/
c,JFINTEG/'integ.atom'/
      DATA PI /3.14159265358979312D0/
C     ..
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C ==>                   -- initialisation step --                    <==
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
c          open(22,STATUS='unknown',FILE='integrand1.dat',
c     &                         FORM='formatted')
c           open(44,STATUS='unknown',FILE='integrand2.dat',
c     &                         FORM='formatted')
c      write(*,*) 'test brahim 2'
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
         WRITE(1337,99000)
#ifdef CPP_MPI
      ie_end = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#else
      ie_end = ielast
#endif
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
c         open(22,STATUS='unknown',FILE='integrand.dat',
c     &                         FORM='formatted')         
         IFGMAT = IFTMAT + 1
         IFMCPA = IFGMAT + 1
         NSMAX = MAX(NAEZ,NATYP)
         NSHCALC = NSHELL(0) - NSMAX
CF77----------------------------------------------------------------
Cccc         IF ( NSHLOC.LT.NSHELD) STOP ' NSHLOC'
Cccc         IF ( NTLOC.LT.NATYP) STOP ' NTLOC'
CF77----------------------------------------------------------------
C
CF90----------------------------------------------------------------
         ALLOCATE (KIJSH(NIJMAX,NSHELL(0)),STAT=LM1)
         IF ( LM1.NE.0 ) THEN
            WRITE(6,99001) 'KIJSH'
            STOP
         END IF
         ALLOCATE (NIJCALC(NSHELL(0)),JIJDONE(NATYP,NATYP,NSHELL(0)),
     &        STAT=LM1)
         IF ( LM1.NE.0 ) THEN
            WRITE(6,99001) 'JIJDONE/NIJCALC'
            STOP
         END IF

c        The array CSUM_STORE could become quite large -> store to MPI-IO-file in future????
c        Only allocate it for MPI-usage. If MPI is not used, two smaller array XINTEGD arrays is needed.
#ifdef CPP_MPI
         ALLOCATE (csum_store(NATYP,NATYP,NSHELL(0),IELAST),STAT=LM1)
         IF ( LM1.NE.0 ) THEN
            WRITE(6,99001) 'csum_store'
            STOP
         END IF
#else
c        ALLOCATE (XINTEGD1(NATYP,NATYP,IELAST),STAT=LM1)
         ALLOCATE (XINTEGD(NATYP,NATYP,NSHELL(0)),STAT=LM1)
         IF ( LM1.NE.0 ) THEN
            WRITE(6,99001) 'XINTEGD'
            STOP
         END IF
#endif
         ALLOCATE (JXCIJINT(NATYP,NATYP,NSHELL(0)),STAT=LM1)
         IF ( LM1.NE.0 ) THEN
            WRITE(6,99001) 'JXCIJINT'
            STOP
         END IF

CF90----------------------------------------------------------------
C..................
#ifdef CPP_MPI
         DO IE=1,IELAST
          DO NS = 1,NSHELL(0)
            DO JT = 1,NATYP
               DO IT = 1,NATYP
                  csum_store(IT,JT,NS,IE)= CZERO
               END DO
            END DO
          END DO
         END DO!IE
#else
         DO NS = 1,NSHELL(0)
            DO JT = 1,NATYP
               DO IT = 1,NATYP
                  XINTEGD(IT,JT,NS)= CZERO
c                  XINTEGD1(IT,JT,IE)= CZERO
               END DO
            END DO
         END DO
#endif
         DO NS = 1,NSHELL(0)
            NIJCALC(NS) = 0
            DO JT = 1,NATYP
               DO IT = 1,NATYP
                  JIJDONE(IT,JT,NS) = 0
                  JXCIJINT(IT,JT,NS) = CZERO
               END DO
            END DO
         END DO
C
         CALL INITABJIJ(IPRINT,NAEZ,NATYP,NATOMIMP,NOFGIJ,NQCALC,
     &                  NSMAX,NSHELL,IQCALC,ATOMIMP,ISH,JSH,IJTABCALC,
     &                  IJTABSH,IJTABSYM,NIJCALC(1),KIJSH(1,1),
     &                  NIJMAX,NSHELL(0),NSHELD)
C --------------------------------------------------------------------
         DO NS = NSMAX + 1,NSHELL(0)
            NSEFF = NS - NSMAX
            DO I1 = 1,NIJCALC(NS)
               IA = ISH(NS,KIJSH(I1,NS))
               JA = JSH(NS,KIJSH(I1,NS))
               LM1 = (IA-1)*NATOMIMP + JA
               IQ = ATOMIMP(IA)
               JQ = ATOMIMP(JA)
C     
               DO L1 = 1,NOQ(JQ)
                  JT = ITOQ(L1,JQ)
                  DO J1 = 1,NOQ(IQ)
                     IT = ITOQ(J1,IQ)
                     JIJDONE(IT,JT,NSEFF) = 1
                  END DO
               END DO
            END DO
         END DO
C ----------------------------------------------------------------------
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
         IF ( IPRINT.GT.0 ) THEN
            DO IT = 1,NATYP
               DO JT = 1,NATYP
                  JTAUX(JT) = 0
                  DO NS = NSMAX + 1,NSHELL(0)
                     NSEFF = NS - NSMAX
                     IF ( JIJDONE(IT,JT,NSEFF).NE.0 ) 
     &                    JTAUX(JT) = JTAUX(JT) + 1
                  END DO
               END DO
Cccc               WRITE (6,99012) IT,(JTAUX(JT),JT=1,MIN(25,NATYP))
            END DO
            WRITE (1337,99013)
         END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C ==>                   INITIALISATION END                           <==
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C
C ==>  read in the single site matrices (TSST) in the LOCAL frame
C
C ==>  set up the Delta t_i matrices DELTSST = TSST(UP) - TSST(DOWN)
C
C ==>  read in projection matrices DMAT and DTIL used to get
C                     ij            ij    _ 
C                    G     =  D  * G    * D
C                     ab       a    CPA    b
C     with a/b the atom of type a/b sitting on site i/j
C     for an atom having occupancy 1, DMAT/DTIL = unit matrix
C
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES
#ifdef CPP_MPI
      ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
      ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

      DO ie_num=1,ie_end
         IE = ie_start+ie_num
#else
      DO IE = 1,IELAST
         ie_num = ie
         ie_end = ielast
#endif
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS SPIN
      DO ISPIN = 1,2
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NATYP
         DO IT = 1,NATYP
c            write(*,*) 'test brahim 3'
            if (t_tgmat%tmat_to_file) then
               IREC = IE + IELAST*(ISPIN-1) + IELAST*2*(IT-1)
               READ (IFTMAT,REC=IREC) W1               
            else
               ie_end = t_mpi_c_grid%ntot2
               IREC = ie_num + ie_end*(ISPIN-1) + ie_end*2*(IT-1)
               W1(:,:) = t_tgmat%tmat(:,:,irec)
            end if
            DO J1 = 1,LMMAXD
               CALL ZCOPY(LMMAXD,W1(1,J1),1,TSST(1,J1,IT,ISPIN),1)
            END DO
C ----------------------------------------------------------------------
            IF ( ISPIN.EQ.2 ) THEN
               DO J1 = 1,LMMAXD
                  DO I1 = 1,LMMAXD
                     DELTSST(I1,J1,IT) = TSST(I1,J1,IT,2)
     &                  - TSST(I1,J1,IT,1)
                  END DO
               END DO
            END IF
C ----------------------------------------------------------------------
            IF ( NCPA.EQ.0 ) THEN
               CALL CINIT(LMMAXD*LMMAXD,DMATTS(1,1,IT,ISPIN))
               DO I1 = 1,LMMAXD
                  DMATTS(I1,I1,IT,ISPIN) = CONE
               END DO
               CALL ZCOPY(LMMAXD*LMMAXD,DMATTS(1,1,IT,ISPIN),1,
     &                    DTILTS(1,1,IT,ISPIN),1)
            ELSE
               READ (IFMCPA,REC=IREC) W1,W2
               DO J1 = 1,LMMAXD
                  DO I1 = 1,LMMAXD
                     DMATTS(I1,J1,IT,ISPIN) = W1(I1,J1)
                     DTILTS(I1,J1,IT,ISPIN) = W2(I1,J1)
                  END DO
               END DO
            END IF
C ----------------------------------------------------------------------
         END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      IF ( IPRINT.GT.1 ) THEN
         WRITE(1337,'(8X,60(1H-),/,10X,
     &        "Single-site and projection matrices read in",
     &        " for IT=1,",I3,/)') NATYP
         DO ISPIN = 1,2
            WRITE(1337,'(8X,60(1H+),/,30X,
     &           " ISPIN = ",I1,/,8X,60(1H+),/)') ISPIN
            DO IT = 1,NATYP
               WRITE(1337,'(12X," IE = ",I2," IT =",I3)') IE,IT
               CALL CMATSTR(' T MAT ',7,TSST(1,1,IT,ISPIN),
     &              LMMAXD,LMMAXD,0,0,0,1D-8,6)
               CALL CMATSTR(' D MAT ',7,DMATTS(1,1,IT,ISPIN),
     &              LMMAXD,LMMAXD,0,0,0,1D-8,6)
               CALL CMATSTR(' D~ MAT',7,DTILTS(1,1,IT,ISPIN),
     &              LMMAXD,LMMAXD,0,0,0,1D-8,6)
               IF ( IT.NE.NATYP) WRITE(1337,'(8X,60(1H-),/)')
            END DO
            WRITE(1337,'(8X,60(1H+),/)')
         END DO
         WRITE(1337,'(8X,60(1H-),/,10X,
     &        "Delta_t = t(it,DN) - t(it,UP) matrices for IT=1,",
     &        I3,/)') NATYP
         DO IT = 1,NATYP
            WRITE(1337,'(12X," IE = ",I2," IT =",I3)') IE,IT
            CALL CMATSTR(' DEL T ',7,DELTSST(1,1,IT),
     &           LMMAXD,LMMAXD,0,0,0,1D-8,6)
            IF ( IT.NE.NATYP) WRITE(1337,'(8X,60(1H-),/)')
         END DO
         WRITE(1337,'(8X,60(1H-),/)')
      END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
C ***************************************************** loop over shells
      DO NS = NSMAX + 1,NSHELL(0)
C
C ======================================================================
C
C ==>  get the off-diagonal Green function matrix Gij(UP) and Gji(DOWN)
C      using the appropiate rotation (similar to ROTGLL routine)
C      step one: read in GS for the representative shell
C ======================================================================
         DO ISPIN = 1,2
          if (t_tgmat%gmat_to_file) then
            IREC = IE + IELAST*(ISPIN-1) + IELAST*2*(NS-1)
            READ (IFGMAT,REC=IREC) W1
          else
            ie_end = t_mpi_c_grid%ntot2
            IREC = ie_num + ie_end*(ISPIN-1) + ie_end*2*(NS-1)
            W1(:,:) = t_tgmat%gmat(:,:,irec) 
          end if
          CALL ZCOPY(LMMAXD*LMMAXD,W1,1,GS(1,1,ISPIN),1)
         END DO
C ======================================================================
C      step two: scan all I,J pairs needed out of this shell 
C                transform with the appropriate symmetry to get Gij
C
C ====================================================== loop over pairs
         NSEFF = NS - NSMAX
         DO L1 = 1,NIJCALC(NS)
C
            IA = ISH(NS,KIJSH(L1,NS))
            JA = JSH(NS,KIJSH(L1,NS))
C
            LM1 = (IA-1)*NATOMIMP + JA
            ISYM = IJTABSYM(LM1)
C ----------------------------------------------------------------------
            DO ISPIN = 1,2
               CALL ZGEMM('C','N',LMMAXD,LMMAXD,LMMAXD,
     +                    CONE,DSYMLL(1,1,ISYM),LMMAXD,
     +                    GS(1,1,ISPIN),LMMAXD,CZERO,W2,LMMAXD)   
C                     
               CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,
     +                    CONE,W2,LMMAXD,DSYMLL(1,1,ISYM),LMMAXD,
     +                    CZERO,W1,LMMAXD)
C     
               IF ( ISPIN.EQ.1 ) THEN
                  CALL ZCOPY(LMMAXD*LMMAXD,W1,1,GMIJ,1)
               ELSE
                  DO LM2 = 1,LMMAXD
                     DO LM1 = 1,LMMAXD
C
C -> use Gji = Gij^T
C
                        GMJI(LM1,LM2) = W1(LM2,LM1)
                     END DO
                  END DO
               END IF
            END DO
C ----------------------------------------------------------------------
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            IF ( IPRINT.GT.2 ) THEN
               WRITE(1337,'(8X,60(1H-),/,10X,
     &              " G_ij(DN) and G_ji(UP) matrices I =",I3," J =",I3,
     &              " IE =",I3,/,8X,60(1H-))') IA,JA,IE
               CALL CMATSTR(' Gij DN',7,GMIJ,LMMAXD,LMMAXD,0,0,0,1D-8,6)
               CALL CMATSTR(' Gji UP',7,GMJI,LMMAXD,LMMAXD,0,0,0,1D-8,6)
            END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C

C ----------------------------------------------------------------------
C
C ==> calculate the exchange coupling constant J_ij via Eq. (19)
C     modified for G instead of tau:
C          J_ij ~ Trace [ (t_i(D)-t_i(U)) * Gij(U)
C                       * (t_j(D)-t_j(U)) * Gji(D)]
C       in case of alloy system: perform projection on atom types
C
C ----------------------------------------------------------------------
            IQ = ATOMIMP(IA)
            JQ = ATOMIMP(JA)
C -------------------------------------------------- loop over occupants
            DO J1 = 1,NOQ(JQ)
               JT = ITOQ(J1,JQ)
               DO I1 = 1,NOQ(IQ)
                  IT = ITOQ(I1,IQ)
C     
C --> Gjq,iq is projected on jt,it ==> Gjt,it
C
                  CALL CMATMUL(LMMAXD,LMMAXD,GMJI,DTILTS(1,1,IT,2),W2)
                  CALL CMATMUL(LMMAXD,LMMAXD,DMATTS(1,1,JT,2),W2,W1)
C
C --> Delta_j * Gjt,it
C
                  CALL CMATMUL(LMMAXD,LMMAXD,DELTSST(1,1,JT),W1,W2)
C     
C --> Giq,jq is projected on it,jt ==> Git,jt * Delta_j * Gjt,it
C
                  CALL CMATMUL(LMMAXD,LMMAXD,GMIJ,DTILTS(1,1,JT,1),W3)
                  CALL CMATMUL(LMMAXD,LMMAXD,DMATTS(1,1,IT,1),W3,W1)
C
C --> Delta_i * Git,jt
C
                  CALL CMATMUL(LMMAXD,LMMAXD,DELTSST(1,1,IT),W1,W3)
C
C --> Delta_i * Git,jt * Delta_j * Gjt,it
C
                  CALL CMATMUL(LMMAXD,LMMAXD,W3,W2,W1)
C
                  CSUM = CZERO
                  DO LM1 = 1,LMMAXD
                     CSUM = CSUM + W1(LM1,LM1)
                  END DO

#ifdef CPP_MPI
!BZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BZ! to ensure the correct order of energy points, the csum=values !!
!BZ! are stored for the MPI version, and the evaluation of JXCIJINT!!
!BZ! and XINTEGD is performed later                                !!
!BZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  csum_store(IT,JT,NSEFF,IE)=CSUM
#else
 
                  JXCIJINT(IT,JT,NSEFF) = JXCIJINT(IT,JT,NSEFF)
     &               - WEZ(IE)*CSUM/DBLE(NSPIN)
                 
                  XINTEGD(IT,JT,NSEFF)= CSUM/(PI*4.D0)

C                  -------> perform substraction instead of addition 
C                           because WGTE ~ -1/pi (WGTE = WEZ(IE)/NSPIN)
                  ! Write out energy-resorved integrand and integral
                  ! Phivos Mavropoulos 24.10.2012
                  FMT2 = '(A,I5.5,A,I5.5,A,I5.5)'
                  WRITE (JFNAM2,FMT2) 'Jij_enrg.',IT,'.',JT,'.',NS
                  IF (IE.EQ.1) THEN
                     OPEN (499,FILE=JFNAM2,STATUS='UNKNOWN')
                     WRITE(499,FMT='(A)')
     &                    '# Energy Re,Im ; j(E) Re,Im; J(E) Re,Im '
                     WRITE(499,FMT='(3(A,I5))')
     &                    '# IT=',IT,' JT=',JT,' SHELL=',NS
                     WRITE(499,FMT='(A,I6)')
     &                    '#ENERGIES: ',IELAST
                  ELSE
                   OPEN (499,FILE=JFNAM2,STATUS='OLD',POSITION='APPEND')
                  ENDIF
                  WRITE (499,FMT='(6E12.4)')
     &            EZ(IE),XINTEGD(IT,JT,NSEFF),JXCIJINT(IT,JT,NSEFF)/4.D0
                  CLOSE(499)
#endif
                  
               END DO !I1
            END DO   !J1, loop over occupants
C ----------------------------------------------------------------------
         END DO   !L1 = 1,NIJCALC(NS)
            
C ======================================================================
        END DO  !loop over shells
C **********************************************************************
c       write(22,*)DIMAG(XINTEGD(1,1,1)),DIMAG(JXCIJINT(1,1,1))
c       write(44,*)DIMAG(XINTEGD(2,2,1)),DIMAG(XINTEGD(2,3,1))
     
      END DO!IE
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES
        

C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC MPI Communication
#ifdef CPP_MPI
        if(myrank==master)then
          !allocate
          ALLOCATE (csum_store2(NATYP,NATYP,NSHELL(0),IELAST),STAT=LM1)
          IF ( LM1.NE.0 ) THEN
             WRITE(6,99001) 'csum_store2'
             STOP
          END IF

          !initialize with zero
          DO IE=1,IELAST
           DO NS = 1,NSHELL(0)
            DO JT = 1,NATYP
               DO IT = 1,NATYP
                  csum_store2(IT,JT,NS,IE)= CZERO
               END DO
            END DO
           END DO
          END DO!IE
        end if!myrank==master

        !perform communication and collect resutls in csum_store2
        LM1 = NATYP*NATYP*NSHELL(0)*IELAST
        call MPI_REDUCE(csum_store,csum_store2,LM1,MPI_DOUBLE_COMPLEX,
     &                  MPI_SUM,master,MPI_COMM_WORLD,IERR)
        if(IERR/=MPI_SUCCESS)then
           write(*,*) 'Problem in MPI_REDUCE(csum_store)'
           stop 'TBXCCPLJIJ'
        end if

        !copy array back to original one
        if(myrank==master)then
          DO IE=1,IELAST
           DO NS = 1,NSHELL(0)
            DO JT = 1,NATYP
               DO IT = 1,NATYP
                  csum_store(IT,JT,NS,IE)= csum_store2(IT,JT,NS,IE)
               END DO
            END DO
           END DO
          END DO!IE

          DEALLOCATE(csum_store2)
        end if!myrank==master



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!~~~~~~~ NEW WRITEOUT: integrand and energy-resolved integral ~~~~~~~~~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      if(myrank==master)then
      DO NS = NSMAX + 1,NSHELL(0)
         NSEFF = NS - NSMAX
         DO L1 = 1,NIJCALC(NS)
            IA = ISH(NS,KIJSH(L1,NS))
            JA = JSH(NS,KIJSH(L1,NS))
            IQ = ATOMIMP(IA)
            JQ = ATOMIMP(JA)
            DO J1 = 1,NOQ(JQ)
               JT = ITOQ(J1,JQ)
               DO I1 = 1,NOQ(IQ)
                  IT = ITOQ(I1,IQ)
C                  -------> perform substraction instead of addition 
C                           because WGTE ~ -1/pi (WGTE = WEZ(IE)/NSPIN)
                  ! Write out energy-resorved integrand and integral
                  ! Phivos Mavropoulos 24.10.2012
                  FMT2 = '(A,I5.5,A,I5.5,A,I5.5)'
                  WRITE (JFNAM2,FMT2) 'Jij_enrg.',IT,'.',JT,'.',NS
                  OPEN (499,FILE=JFNAM2,STATUS='UNKNOWN')
                  WRITE(499,FMT='(A)')
     &                 '# Energy Re,Im ; j(E) Re,Im; J(E) Re,Im '
                  WRITE(499,FMT='(3(A,I5))')
     &                 '# IT=',IT,' JT=',JT,' SHELL=',NS
                  WRITE(499,FMT='(A,I6)')
     &                 '#ENERGIES: ',IELAST

                 DO IE=1,IELAST
                  JXCIJINT(IT,JT,NSEFF) = JXCIJINT(IT,JT,NSEFF)
     &               - WEZ(IE)*csum_store(IT,JT,NSEFF,IE)/DBLE(NSPIN)
                  XINTEGDTMP= 
     &                csum_store(IT,JT,NSEFF,IE)/(PI*4.D0)
                   WRITE (499,FMT='(6E12.4)')
     &            EZ(IE),XINTEGDTMP,JXCIJINT(IT,JT,NSEFF)/4.D0

                 END DO!IE

                 close(499)

               END DO !I1
            END DO   !J1, loop over occupants
         END DO   !L1 = 1,NIJCALC(NS)
        END DO  !NS, loop over shells
      end if!myrank==master
#endif

       if(myrank==master)then
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
       WGTEMP = CONE/4D0
C     -------> here factor 1/pi omitted since it is included in WGTE
      DO NS = 1,NSHCALC
         DO JT = 1,NATYP
            DO IT = 1,NATYP
               JXCIJINT(IT,JT,NS) = WGTEMP*JXCIJINT(IT,JT,NS)
            END DO
         END DO
         NSEFF = NS + NSMAX
         RSH(NS) = 0D0
         DO I1 = 1,3
            RSH(NS) = RSH(NS) + RATOM(I1,NSEFF)*RATOM(I1,NSEFF)
         END DO
         RSH(NS) = SQRT(RSH(NS))
      END DO
C
      LM2 = 0
      NTCALC = 0
      DO IT = 1,NATYP
         L1 = 0
         DO NS = 1,NSHCALC
            LM1 = 0
            DO JT = 1,NATYP
               IF ( JIJDONE(IT,JT,NS).NE.0 ) LM1 = LM1 + 1
            END DO
            LM2 = MAX(LM1,LM2)
            L1 = MAX(L1,LM1)
         END DO
         IF ( L1.GT.0 ) THEN
            NTCALC = NTCALC + 1
            JTAUX(NTCALC) = IT
         END IF
      END DO
      WRITE (STRBAR,'(19(1H-))')
      LSTR = 19
      DO I1 = 1,LM2
         WRITE(STRTMP,'(A,15(1H-))') STRBAR(1:LSTR)
         LSTR = LSTR+15
         STRBAR(1:LSTR)=STRTMP(1:LSTR)
      END DO
      
      WRITE(1337,99002) STRBAR(1:LSTR),STRBAR(1:LSTR)
      DO I1 = 1,NTCALC
         IT = JTAUX(I1)
         WRITE (1337,99003) IT,IQAT(IT)
         L1 = 0
         DO NS = 1,NSHCALC
            LM1 = 0
            DO JT = 1,NATYP
               IF ( JIJDONE(IT,JT,NS).NE.0 ) LM1 = LM1 + 1
            END DO
            IF ( LM1.NE.0 ) THEN
               LM2 = 0
               IF ( L1.EQ.0 ) THEN
                  WRITE(1337,99005) RSH(NS)
                  L1 = 1
               ELSE
                  WRITE (1337,99006) RSH(NS)
               END IF
               DO JT = 1,NATYP
                  IF ( JIJDONE(IT,JT,NS).NE.0 ) THEN
                     LM2 = LM2 + 1
                     IF ( LM2.EQ.1 ) THEN
                        WRITE (1337,99007) IQAT(JT),
     &                       DIMAG(JXCIJINT(IT,JT,NS))*1D3,JT
                        write(1337,*) ns+nsmax,' shell'
                     ELSE
                        WRITE (1337,99008) 
     &                       DIMAG(JXCIJINT(IT,JT,NS))*1D3,JT
                        write(1337,*) ns+nsmax,' shell'
                     END IF
                     IF ( LM2.EQ.LM1 ) WRITE (1337,*)
                  END IF
               END DO
            END IF
         END DO
         WRITE (1337,99004) STRBAR(1:LSTR)
      END DO !I1 = 1,NTCALC
      WRITE (1337,*)
C ----------------------------------------------------------------------
C --> prepare output files 
C
      DO I1 = 1,NTCALC
         IT = JTAUX(I1)
         DO NS = 1,NSHCALC
            L1 = 1
            DO JT = 1,NATYP
               IF ( JIJDONE(IT,JT,NS).NE.0 ) THEN
                  JIJDONE(IT,L1,NS) = JT
                  JXCIJINT(IT,L1,NS) = JXCIJINT(IT,JT,NS)
                  L1 = L1 + 1
               END IF
            END DO
            DO JT = L1,NATYP
               JIJDONE(IT,JT,NS) = 0
            END DO
         END DO
      END DO
C
      DO I1 = 1,NTCALC
         IT = JTAUX(I1)
!         FMT1 = '(A,I1)'
!         IF ( IT.GE.10 ) THEN
!            FMT1 = '(A,I2)'
!            IF ( IT.GE.100 ) FMT1 = '(A,I3)'
!         END IF
         FMT1 = '(A,I5.5)'
         WRITE (JFNAM,FMT1) JFBAS,IT

c         write(JFNAME,FMT1)JFINTEG, IT
         OPEN (49,FILE=JFNAM)
         WRITE (49,99009) IT,IQAT(IT)
c         open(22,FILE=JFNAME)
c          write(22,99014)IT,IQAT(IT)
         DO L1 = 1,NATYP
            LM1 = 0
            DO NS = 1,NSHCALC
               IF ( JIJDONE(IT,L1,NS).NE.0 ) THEN
                  LM1 = LM1 + 1
                  WRITE (49,99010) RSH(NS),
     &                 DIMAG(JXCIJINT(IT,L1,NS)),
     &                 JIJDONE(IT,L1,NS),NS+NSMAX ! fivos added NS+NSMAX
c              do IE=1,IELAST
c                XINTEGD(IT,JT,IE)= XINTEGD(IT,JT,IE)*(1.D0/(PI*4.D0))
                
c                write(22,99015)XINTEGD(IT,JT,IE),JIJDONE(IT,L1,NS)
c             end do 
               END IF
            END DO
c            end do
            IF ( (LM1.NE.0) .AND. (L1.NE.NATYP) ) WRITE (49,*) '&'
c           IF ( (LM1.NE.0) .AND. (L1.NE.NATYP) ) WRITE (22,*) '&'
         END DO !l1=1,natyp
         CLOSE (49)
c         Close(22)
      END DO !i1=1,ntcalc
      WRITE (1337,99011) 
C ----------------------------------------------------------------------
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      endif!myrank==master
C
CF90----------------------------------------------------------------
      DEALLOCATE (KIJSH,NIJCALC,JIJDONE,JXCIJINT,STAT=LM1)
CF90----------------------------------------------------------------
99000 FORMAT(79(1H=),/,10X,
     &     "TBXCCPLJIJ : Off-diagonal exchange coupling",
     &     " constants J_ij",/,79(1H=),/)
99001 FORMAT(6X,"ERROR: Cannot allocate array(s) :",A,/)
99002 FORMAT(4X,A,4(1H-),/,5X," IT/IQ ",3X,"R_IQ,JQ",2X,"JQ ",
     &       " (J_IT,JT  JT)",/,15X," [ ALAT ] ",4X,"[ mRy ]",/,
     &       4X,A,4(1H-))
99003 FORMAT(5X,I3,1X,I3,$)
99004 FORMAT(4X,A,4(1H-))
99005 FORMAT(F10.6,$)
99006 FORMAT(12X,F10.6,$)
99007 FORMAT(I4,F12.8,I3,$)
99008 FORMAT(F12.8,I3,$)
99009 FORMAT("# off-diagonal exchange coupling constants ",/,
     &     "# for atom IT = ",I3," on site IQ = ",I3,/,
     &     "# R_IQ,JQ      J_IT,JT       JT",/,
     &     "# ( ALAT )       ( Ry )",/,"#      ")
99010 FORMAT(F12.8,2X,E14.8,2X,I3,I5)
99011 FORMAT(6X,"Output written into the files Jij.atomX",/,
     &       6X,"  X = atom index in the input file")
99012       FORMAT (10X,I3,3X,25(I3))
99013       FORMAT (/,8X,60('-'),/)

99014       FORMAT("# Integrand  ",/,
     &     "# for atom IT = ",I3," on site IQ = ",I3,/,
     &     "#      J_IT,JT       JT",/,
     &     "#        ( Ry )",/,"#      ")    
99015 format(F12.8)      
              END SUBROUTINE TBXCCPLJIJ

        END MODULE MOD_TBXCCPLJIJ
