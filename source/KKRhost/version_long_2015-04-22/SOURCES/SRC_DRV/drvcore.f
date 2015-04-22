C*==drvcore.f    processed by SPAG 6.05Rc at 11:35 on 10 May 2004
      SUBROUTINE DRVCORE(IPRINT,ITPRT,LCORE,NCORE,CSCL,VTIN,BTIN,RIN,
     &                   A,B,DRDIIN,R2DRDIIN,ZAT,JWS,ISHIFT,RHOC,
     &                   ECOREREL,NKCORE,KAPCORE,ECORE,LMAXD,IRMD)
C   ********************************************************************
C   *                                                                  *
C   * driving routine to call relativistic < CORE > routine            *
C   * counterpart of < COREL > of the non- or scalar-relativistic mode *
C   *                                                                  *
C   * ATTENTION: all the variables connected with hyperfine fields are *
C   *            OFF                                                   *
C   *                                                                  *
C   * The non-relativistic variables NCORE and LCORE are used as read  *
C   * in from the potential file; they are not modified and are again  *
C   * written out for the next iteration ( routine <RITES>)            *
C   *                                                                  *
C   * Relativistic variables passed outside this routine:              *
C   *    NKCORE(1..NCORE)       = number of KAPPA values for a given   *
C   *                             (n,l) core state                     *
C   *    KAPCORE(1..NCORE,1..NCORE) = the (maximum 2) values of KAPPA  *
C   *
C   *    ECOREREL(1..NCORE,1..NCORE) = for a given (n,l) state the core*
C   *                              energies corresponding first/second *
C   *                              KAPPA value, AVERAGED over \mu's    *
C   *                              These values are written out to the *
C   *                              potential file (routine <RITES>),   *
C   *                              but the read in (routine <STARTB1>) *
C   *                              updates the ECORE array             *
C   *     Please note that ALL the core energies (also \mu-resolved)   *
C   *     are output by <CORE> routine but not passed out of here      *
C   *                                                                  *
C   *    ECORE(1..NCORE,1..2) = this (non- or scalar relativistic)     *
C   *                           variable is updated here to be used in *
C   *                           calculating/printing out the spin-     *
C   *                           resolved energies (see <ESPCB> )       *
C   *                           A SUMMATION is done here:              *
C   *        ECORE(L,1/2) = SUM_{\kappa=-L-1,L} E(\kappa,\mu)          *
C   *                       /(2*L+1)                                   *
C   *                           with negative \mu's for E(*,1) and     *
C   *                           positive \mu's for E(*,2)              *
C   *                                                                  *
C   *                           v.popescu July/2002                    *
C   ********************************************************************
      IMPLICIT NONE
C
C PARAMETER definitions
C
      INTEGER NTMAX,NMMAX,NCSTMAX,NMEMAX,NLMAX,NKMMAX
      PARAMETER (NTMAX=1,NMMAX=1,NCSTMAX=6,NMEMAX=5,NLMAX=5,
     &           NKMMAX=2*NLMAX**2) ! NLMAX should be >= LCOREMAX + 1
      INTEGER NRMAX
      PARAMETER (NRMAX = 750)
      DOUBLE PRECISION DZERO
      PARAMETER (DZERO=0.0D0)
C
C Dummy arguments
C
      DOUBLE PRECISION A,B
      INTEGER IPRINT,IRMD,ITPRT,LMAXD,NCORE,ISHIFT
C
C --> obs: here, in contrast to DRVRHO, one works with local
C     arrays, since they have to be set up as far as NRMAX (see CORE)
C
      DOUBLE PRECISION VTIN(IRMD),BTIN(IRMD)
      DOUBLE PRECISION DRDIIN(IRMD),R2DRDIIN(IRMD)
      DOUBLE PRECISION RIN(IRMD),CSCL(LMAXD+1)
C
      DOUBLE PRECISION ECORE(20,2),ECOREREL(20*2)
      INTEGER KAPCORE(20*2),LCORE(20,2),NKCORE(20)
      INTEGER ZAT(NTMAX),JWS(NMMAX)
      DOUBLE PRECISION RHOC(IRMD,2)
C
C Local variables
C
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),
     &       EA,ECOR(NCSTMAX),ECORTAB(120,NTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),QDIA(NKMMAX),
     &       QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX),
     &       RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX),
     &       SDIA(NKMMAX),SMDIA(NKMMAX),SMOFF(NKMMAX),SOFF(NKMMAX),
     &       SZCOR(NCSTMAX)
      DOUBLE PRECISION VT(NRMAX,NTMAX),BT(NRMAX,NTMAX),CTL(NTMAX,NLMAX)
      DOUBLE PRECISION DRDI(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX)
      DOUBLE PRECISION R(NRMAX,NMMAX)
      INTEGER I,ICALL,IKMCOR(NCSTMAX,2),IMT(NTMAX),IP,ISMQHFI,
     &        IT,ITXRAY,IZERO(NCSTMAX),J,KAPCOR(NCSTMAX),
     &        LCXRAY(NTMAX),MM05COR(NCSTMAX),NCORT(NTMAX),NCXRAY(NTMAX),
     &        NKPCOR(NCSTMAX),NT,NUCLEUS
      INTEGER LCOREMAX
C
      SAVE IMT,ITXRAY,NT,NUCLEUS,QDIA,QMDIA,QMOFF,QOFF,SDIA,SMDIA,
     &     SMOFF,SOFF
C
      DATA NCXRAY/NTMAX*0/,LCXRAY/NTMAX*0/,ISMQHFI/0/
C
      DATA ICALL/0/
C
      ICALL = ICALL + 1
C
C=======================================================================
C       initialise relativistic and dummy variables and SAVE them
C=======================================================================
      IF ( ICALL.EQ.1 ) THEN
C
         IF ( LMAXD.GT.NLMAX-1 ) THEN
            WRITE (6,*) ' LMAXD = ',LMAXD,' > NLMAX-1 = ',NLMAX - 1
            STOP ' Increase NLMAX in < DRVCORE > '
         END IF
C
         IF ( IRMD.GT.NRMAX ) THEN
            WRITE (6,*) ' IRMD = ',IRMD,' > NRMAX = ',NRMAX
            WRITE (6,*) ' Increase NRMAX in < sprkkr_rmesh.dim > '
            STOP ' In < DRVCORE > '
         END IF
C
         ITXRAY = 0
C
         DO IT = 1,NTMAX
            IMT(IT) = 1
         END DO
C     
         NT = 1
         NUCLEUS = 0
C     
         DO IT = 1,NKMMAX
            SDIA(IT) = DZERO
            SMDIA(IT) = DZERO
            SOFF(IT) = DZERO
            SMOFF(IT) = DZERO
C     
            QDIA(IT) = DZERO
            QMDIA(IT) = DZERO
            QOFF(IT) = DZERO
            QMOFF(IT) = DZERO
         END DO
C     
      END IF                    ! ICALL.EQ.1
C=======================================================================
C
C --> fill up CTL array for the case of core states with higher L values
C     than those used in the valence band
C
      LCOREMAX = 0
      DO IT = 1,NCORE
         J = LCORE(IT,1)
         LCOREMAX = MAX(LCOREMAX,J)
      END DO
      IF ( LCOREMAX.GT.NLMAX-1 ) THEN
         WRITE (6,*) ' LCOREMAX = ',LCOREMAX,' > NLMAX-1 = ',NLMAX - 1
         STOP ' Increase NLMAX in < DRVCORE > '
      END IF
      DO J = 1,LMAXD+1
         CTL(1,J) = CSCL(J)
      END DO
      IF ( LCOREMAX.GT.0 ) THEN
         DO J = LCOREMAX+1,NLMAX
            CTL(1,J) = CSCL(LCOREMAX)
         END DO
      END IF
C
      CALL DCOPY(JWS(1),VTIN,1,VT(1,1),1)
      CALL DCOPY(JWS(1),BTIN,1,BT(1,1),1)
      CALL DCOPY(JWS(1),RIN,1,R(1,1),1)
      CALL DCOPY(JWS(1),DRDIIN,1,DRDI(1,1),1)
      CALL DCOPY(JWS(1),R2DRDIIN,1,R2DRDI(1,1),1)
C
      DO J = JWS(1) + 1,NRMAX
         EA = DEXP(A*DBLE(J+ISHIFT-1)) ! corrected from (J-1) 07.05.2004
         R(J,1) = B*(EA-1D0)
         DRDI(J,1) = A*B*EA
         R2DRDI(J,1) = R(J,1)*R(J,1)*DRDI(J,1)
         VT(J,1) = 0D0
         BT(J,1) = 0D0
      END DO
C
      NCORT(1) = 0  ! no. of core electrons = no. of diff. core states
C
      DO IT = 1,NCORE
         NCORT(1) = NCORT(1) + 2*(2*LCORE(IT,1)+1)
      END DO
C
      CALL CORE(IPRINT,ITPRT,NT,NCORT,CTL,VT,BT,ZAT,NUCLEUS,R,R2DRDI,
     &          DRDI,JWS,IMT,RHOCHR,RHOSPN,ECORTAB,GCOR,FCOR,ECOR,SZCOR,
     &          KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,NCXRAY,LCXRAY,ITXRAY,
     &          BCOR,BCORS,SDIA,SMDIA,SOFF,SMOFF,QDIA,QOFF,QMDIA,QMOFF,
     &          NKMMAX,NMEMAX,ISMQHFI,NTMAX,NRMAX,NMMAX,NCSTMAX,NLMAX)
C
      CALL RINIT(2*IRMD,RHOC(1,1))
C
      DO I = 1,JWS(1)
         IP = I + ISHIFT
         RHOC(IP,2) = (RHOCHR(I,1)+RHOSPN(I,1))*0.5D0*(R(I,1)**2)
         RHOC(IP,1) = (RHOCHR(I,1)-RHOSPN(I,1))*0.5D0*(R(I,1)**2)
      END DO
C
      CALL SUMECORE(NCORE,LCORE(1,1),ECORTAB(1,1),NKCORE,ECOREREL,ECORE,
     &              KAPCORE)
C
      END
C*==sumecore.f    processed by SPAG 6.05Rc at 11:35 on 10 May 2004
C
      SUBROUTINE SUMECORE(NCORE,LCORE,ECORTAB,NKCORE,ECOREREL,ECORE,
     &                    KAPCORE)
      IMPLICIT NONE
C
C Dummy arguments
C
      INTEGER NCORE
      DOUBLE PRECISION ECORE(20,2),ECOREREL(20*2),ECORTAB(*)
      INTEGER KAPCORE(20*2),LCORE(*),NKCORE(*)
C
C Local variables
C
      DOUBLE PRECISION DBLE
      INTEGER I,IC,ICREL,JREL,KFG(4),L,LMP1,LMXC,LP1,MUEM05,NC,NMAX,NN,
     &        NSOL,WGT(2)
      DOUBLE PRECISION MJ
      INTRINSIC ABS
C
C --> find the principal quantum numbers
C
      DO IC = 1,4
         KFG(IC) = 0
      END DO
      DO IC = 1,NCORE
         IF ( LCORE(IC).EQ.0 ) KFG(1) = KFG(1) + 1
         IF ( LCORE(IC).EQ.1 ) KFG(2) = KFG(2) + 1
         IF ( LCORE(IC).EQ.2 ) KFG(3) = KFG(3) + 1
         IF ( LCORE(IC).EQ.3 ) KFG(4) = KFG(4) + 1
      END DO
C
      IF ( KFG(2).NE.0 ) KFG(2) = KFG(2) + 1
      IF ( KFG(3).NE.0 ) KFG(3) = KFG(3) + 2
      IF ( KFG(4).NE.0 ) KFG(4) = KFG(4) + 3
C
      LMXC = 0
      IF ( KFG(2).NE.0 ) LMXC = 1
      IF ( KFG(3).NE.0 ) LMXC = 2
      IF ( KFG(4).NE.0 ) LMXC = 3
C
      LMP1 = LMXC + 1
      NC = 0
      DO LP1 = 1,LMP1
         L = LP1 - 1
         NMAX = KFG(LP1)
         DO NN = LP1,NMAX
            NC = NC + 1
            ECOREREL(NC) = 0.0D0
            ECOREREL(NC+20) = 0.0D0
C
C  ECOREREL(NC..NC+20) = 1st/2nd value of \kappa for current l
C
C
C --> icrel is pointing a core-state (n,l) = (nn,l) from the
C     relativistic-routine sequence 1s,2s,2p,3s,3p,... to the
C     ECOREREL array  1s,2s,3s,2p,3p,...
C
C     (nn,l) -> nn*(nn-1)*(2*nn-1)/3 + 2*l**2
C
            ICREL = NN*(NN-1)*(2*NN-1)/3 + 2*L**2
            JREL = 0
            WGT(1) = 0
            WGT(2) = 0
            DO MUEM05 = -L - 1, + L
               MJ = MUEM05 + 0.5D0
               IF ( ABS(MJ).GT.L ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
               DO I = 1,NSOL
                  WGT(I) = WGT(I) + 1
                  JREL = JREL + 1
                  ECOREREL((I-1)*20+NC) = ECOREREL((I-1)*20+NC)
     &               + ECORTAB(ICREL+JREL)
               END DO
            END DO
            NKCORE(NC) = 1
            IF ( L.NE.0 ) NKCORE(NC) = 2
            KAPCORE(NC) = -L - 1
            KAPCORE(NC+20) = L
C
            DO I = 1,NKCORE(NC)
               ECOREREL((I-1)*20+NC) = ECOREREL((I-1)*20+NC)
     &                                 /DBLE(WGT(I))
            END DO
C
C --> update the array ECORE(1..NCORE,UP/DOWN) as
C
C         ECORE(L,SIGMA) = 1/(2*L+1) * 
C              SUM_KAPPA(L) SUM_(MUE,SIGN(MUE)=SIGN(SIGMA)) ECORTAB(KAP,MUE)
C     i.e., states with negative MUE are added to SPIN DOWN, those with
C     positive MUE to SPIN UP.
C
C     ECORE is used later in calculating the total energy
C
            ECORE(NC,1) = 0.0D0
            ECORE(NC,2) = 0.0D0
            DO I = 1,2*L + 1
               ECORE(NC,1) = ECORE(NC,1) + ECORTAB(ICREL+I)
               ECORE(NC,2) = ECORE(NC,2) + ECORTAB(ICREL+2*L+1+I)
            END DO
            DO I=1,2
               ECORE(NC,I) = ECORE(NC,I)/DBLE(2*L+1)
            END DO
         END DO
      END DO
      END
