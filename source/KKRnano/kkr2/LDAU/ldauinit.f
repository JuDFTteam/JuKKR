C=======================================================================
C initialization of LDA+U calculation operating
C=======================================================================
C
      SUBROUTINE LDAUINIT(
     >                    I1,ITER,NSRA,NLDAU,LLDAU,ULDAU,JLDAU,EREFLDAU,
     >                    VISP,NSPIN,
     >                    R,DRDI,ZAT,IPAN,IRCUT,
     <                    PHILDAU,UMLDAU,WMLDAU,
C                         new input parameters after inc.p removal
     &                    lmax, irmd, ipand)
C
      IMPLICIT NONE

      INTEGER lmax
      INTEGER irmd
      INTEGER ipand
C
C     INTEGER             LMAXD1
C     PARAMETER          (LMAXD1 = LMAXD + 1)
C     INTEGER             MMAXD
C     PARAMETER          (MMAXD=2*LMAXD+1)
C
C global arrays ..
C     DOUBLE COMPLEX     PHILDAU(IRMD,LMAXD1)
C     DOUBLE PRECISION   UMLDAU(MMAXD,MMAXD,MMAXD,MMAXD,LMAXD1),
C    +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1),
C    +                   ULDAU(LMAXD1),
C    +                   JLDAU(LMAXD1),
C    +                   VISP(IRMD,2),
C    +                   R(IRMD),
C    +                   DRDI(IRMD)
C     INTEGER            LLDAU(LMAXD1),
C    +                   IRCUT(0:IPAND)

      DOUBLE COMPLEX     PHILDAU(IRMD,LMAX + 1)
      DOUBLE PRECISION   UMLDAU(2*LMAX+1, 2*LMAX+1,
     &                          2*LMAX+1, 2*LMAX+1, LMAX + 1)

      DOUBLE PRECISION   WMLDAU(2*LMAX+1,2*LMAX+1,LMAX + 1,NSPIN)
      DOUBLE PRECISION   ULDAU(LMAX + 1)
      DOUBLE PRECISION   JLDAU(LMAX + 1)
      DOUBLE PRECISION   VISP(IRMD,2)
      DOUBLE PRECISION   R(IRMD)
      DOUBLE PRECISION   DRDI(IRMD)
      INTEGER            LLDAU(LMAX + 1)
      INTEGER            IRCUT(0:IPAND)
C ..
C global scalars
      DOUBLE PRECISION   ZAT,EREFLDAU
      INTEGER            I1,ITER,NLDAU,NSPIN,NSRA,IPAN
C ..
C local arrays ..

C     Fortran 90 - automatic arrays
C     DOUBLE PRECISION   AA(MMAXD,MMAXD,MMAXD,MMAXD,0:2*LMAX) ! large?
      DOUBLE PRECISION   AA(2*LMAX+1,2*LMAX+1,
     &                      2*LMAX+1,2*LMAX+1,0:2*LMAX)

      DOUBLE PRECISION   FACT(0:100),
     +                   RPW(IRMD,2*LMAX+1),
     +                   TG(IRMD),TL(IRMD),SG(IRMD),SL(IRMD),
     +                   WINT(IRMD),W2(IRMD),WGTFCLMB(0:2*LMAX+1),
     +                   FCLMB(0:2*LMAX+1)
C ..
C local scalars
      DOUBLE PRECISION   PI,RLOP,SUMFCLMB,WIG3J,SCL,SUM,G12,G34,
     +                   CGCRAC,GAUNTC
      INTEGER            LRECLDAU,
     +                   IU,ILDAU,
     +                   IRC1,LFMAX,IR,L1,LF,LL,KK,
     +                   M1,M2,M3,M4,
     +                   IM1,IM2,IM3,IM4
      LOGICAL            WINFO
C ..
C
      INTEGER             LMAXD1
      INTEGER             MMAXD
      INTEGER             NSPIND
      INTEGER             LMAXD

      LMAXD1 = LMAX + 1
      MMAXD  = 2*LMAX+1
      NSPIND = NSPIN
      LMAXD = LMAX

      PI = 4.D0*ATAN(1.0D0)
      FACT(0) = 1.0D0
      DO IU = 1,100
         FACT(IU) = FACT(IU-1)*DBLE(IU) ! factorial
      END DO
C
      WINFO = .FALSE.                   ! set to .true. if information on initilization 
C                                         should be written to output
C
C
C read information on LDAU in any case ..
C
      LRECLDAU = 4*(1+LMAXD1)                   ! NLDAU & LLDAU
     +         + 8*2*LMAXD1                     ! ULDAU & JLDAU
     +         + 8*MMAXD*MMAXD*NSPIND*LMAXD1    ! WMLDAU
c
      OPEN (65,ACCESS='direct',RECL=LRECLDAU,FILE='wldau.unf',
     +      FORM='unformatted')
      READ (65,REC=I1) NLDAU,LLDAU,ULDAU,JLDAU,WMLDAU
      CLOSE(65)
C
C ..
C
C only proceed if orbitals exist where LDAU should be applied to
C
      IF (NLDAU.GE.1.AND.ITER.EQ.1) THEN
        WRITE(6,*) 'LDA+U: set up basis'
C      IF (NLDAU.GE.1) THEN
C
C loop over orbitals ..
C
        DO ILDAU = 1,NLDAU
C
C -> Calculate test functions Phi. Phi is already normalised to
C    int phi**2 dr =1, thus it also contains factor r.
C
          CALL LDAUPHI(LLDAU(ILDAU),VISP,IPAN,IRCUT,R,DRDI,ZAT,
     &                 EREFLDAU,PHILDAU(1,ILDAU),NSPIN,NSRA,
     &                 NLDAU,LLDAU,
     &                 lmaxd, irmd, ipand)
C
        ENDDO
C
        IRC1 = IRCUT(IPAN)
C
C -> define r**l in array rpw:
C
        DO ILDAU = 1,NLDAU
          LFMAX = 2*LLDAU(ILDAU)
          DO IR = 2,IRMD
            RPW(IR,1) = R(IR)
            DO L1 = 2,2*LMAXD + 1
              RPW(IR,L1) = R(IR)*RPW(IR,L1-1)
            ENDDO
          ENDDO
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C 
C 1.  Calculate slater integrals FCLMB
C     (here using only the large component of sra wavefunction.
C     Whoever wants can also do the small component.)
C
C 1a. Calculate slater integrals
C
          RLOP = DBLE(LLDAU(ILDAU))
          SUMFCLMB = 0.D0
C ----------------------------------------------------------------------
          DO LF = 2,LFMAX,2
            TL(1) = 0.0D0
            TG(1) = 0.0D0
C     
C Note on integrand:
C Integrals are up to IRC1 because we integrate in sphere, 
C without thetas.
C In case of cell integration, from IRCUT(1)+1 to IRC1 a convolution
C with thetas and gaunts is needed: 
C Int dr R_l(r)**2 Sum_L' Gaunt_{lm,lm,l'm'}*thetas_{l'm'}(r).
C But then the result is m-dependent. Here, the easy way is used!
C     
            DO IR = 2,IRC1
              WINT(IR) = 
     +        DREAL( DCONJG(PHILDAU(IR,ILDAU)) * PHILDAU(IR,ILDAU) )
              W2(IR) = 2.D0 * DRDI(IR) * WINT(IR)
              TL(IR) = W2(IR) * RPW(IR,LF)
              TG(IR) = W2(IR) / RPW(IR,LF+1)
            ENDDO
C
            CALL SOUTK(TL,SL,IPAN,IRCUT(0))
            CALL SOUTK(TG,SG,IPAN,IRCUT(0))
C
            SL(1) = 0.0D0
            DO IR = 2,IRC1
              SL(IR) = SL(IR)/RPW(IR,LF+1)
     &               + (SG(IRC1)-SG(IR))*RPW(IR,LF)
            ENDDO
C
            SG(1) = 0.0D0
C
C See Note on integrand above.
C
            DO IR = 2,IRC1
              SG(IR) = WINT(IR) * SL(IR)
            ENDDO
C
            CALL SIMPK(SG,FCLMB(LF),IPAN,IRCUT(0),DRDI(1))
C
            WIG3J = (-1)**NINT(RLOP) * (1D0/SQRT(2D0*RLOP+1D0))
     &            * CGCRAC(FACT,RLOP,DBLE(LF),RLOP,0D0,0D0,0D0)
C
            WGTFCLMB(LF) = ((2*RLOP+1)/(2*RLOP))*WIG3J**2
            SUMFCLMB     = SUMFCLMB + WGTFCLMB(LF)*FCLMB(LF)
C
          ENDDO
C ----------------------------------------------------------------------
C
C 1b.   Normalise slater integrals FCLMB
C
          SCL = JLDAU(ILDAU)/SUMFCLMB
          FCLMB(0) = ULDAU(ILDAU)
C
          DO LF = 2,LFMAX,2
            FCLMB(LF) = SCL*FCLMB(LF)
          ENDDO
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C ======================================================================
C
C 2.   Calculate coefficient matrix AA = a(m1,m2,m3,m4)
C
          CALL RINIT(MMAXD*MMAXD*MMAXD*MMAXD*(2*LMAXD+1),AA)
          LL = LLDAU(ILDAU)
          DO LF = 0,LFMAX,2
            DO M3 = -LL,LL
              IM3 = LL + M3 + 1
              DO M2 = -LL,LL
                IM2 = LL + M2 + 1
                DO M1 = -LL,LL
                  IM1 = LL + M1 + 1
                  M4 = M1 - M2 + M3
                  IF ( -LL.LE.M4 .AND. M4.LE.LL ) THEN
                    IM4 = LL + M4 + 1
                    SUM = 0.D0
C
                    DO KK = -LF,LF
                      G12 = GAUNTC(FACT,LL,M1,LF,KK,LL,M2)
                      G34 = GAUNTC(FACT,LL,M3,LF,-KK,LL,M4)
                      SUM = SUM + G12*G34*(-1)**ABS(KK)
                    ENDDO
C
                    AA(IM1,IM2,IM3,IM4,LF)
     &                     = SUM*4.D0*PI/(2.D0*DBLE(LF)+1.D0)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C ======================================================================
C
C UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C 3.  Calculate ULDAU
C
          CALL RINIT(MMAXD*MMAXD*MMAXD*MMAXD,UMLDAU(1,1,1,1,ILDAU))
C
          DO LF = 0,LFMAX,2
            DO IM4 = 1,2*LL + 1
              DO IM3 = 1,2*LL + 1
                DO IM2 = 1,2*LL + 1
                  DO IM1 = 1,2*LL + 1
                    UMLDAU(IM1,IM2,IM3,IM4,ILDAU)
     &                         = UMLDAU(IM1,IM2,IM3,IM4,ILDAU)
     &                         + AA(IM1,IM2,IM3,IM4,LF)*FCLMB(LF)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C
C UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C
          IF (WINFO) THEN
C
          WRITE (6,'(/,8X,A,I3,/)') 'ATOM: ',I1
          WRITE (6,'(12X,A,F8.4,A)') 'LDA+U reference energy :',
     &        EREFLDAU,' Ry'
          WRITE (6,'(12X,A,2F8.4,A,/)') 'Ueff and Jeff = ',ULDAU(ILDAU),
     &        JLDAU(ILDAU),' Ry'
          WRITE (6,'(12X,"Scaling factor for F^n :",F10.6,/)') SCL
          WRITE (6,'(12X,"  n   F^n calculated   F^n scaled ")')
          DO LF = 2,LFMAX,2
            WRITE (6,'(12X,I3,2(2X,F12.8," Ry"))') 
     &             LF,FCLMB(LF)/SCL,FCLMB(LF)
          ENDDO
C
          ENDIF
C
        ENDDO
C
      ENDIF
C
C
      END
C
C
C
C
C
C
C
C
C*==cgcrac.f    processed by SPAG 6.05Rc at 15:27 on  7 Mar 2003
C
      FUNCTION CGCRAC(FACT,J1,J2,J3,M1,M2,M3)
C   ********************************************************************
C   *                                                                  *
C   *     CLEBSCH GORDAN COEFFICIENTS FOR ARBITRARY                    *
C   *     QUANTUM NUMBERS  J1,J2 ...                                   *
C   *     ACCORDING TO THE FORMULA OF   RACAH                          *
C   *     SEE: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM          *
C   *          EQUATION (3.19)                                         *
C   *          EDMONDS EQ. (3.6.11) PAGE 45                            *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C Dummy arguments
C
      DOUBLE PRECISION J1,J2,J3,M1,M2,M3
      DOUBLE PRECISION CGCRAC
      DOUBLE PRECISION FACT(0:100)
C
C Local variables
C
      DOUBLE PRECISION DSQRT
      INTEGER J,N,N1,N2,N3,N4,N5,NBOT,NTOP
      INTEGER NINT
      DOUBLE PRECISION RFACT
      DOUBLE PRECISION S,SUM,VF,X,Y
C
C INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
      RFACT(X) = FACT(NINT(X))
C
C
      CGCRAC = 0.0D0
      IF ( ABS(M3-(M1+M2)).GT.1.0D-6 ) RETURN
      IF ( ABS(J1-J2).GT.J3 ) RETURN
      IF ( (J1+J2).LT.J3 ) RETURN
      IF ( ABS(M1).GT.(J1+1.0D-6) ) RETURN
      IF ( ABS(M2).GT.(J2+1.0D-6) ) RETURN
      IF ( ABS(M3).GT.(J3+1.0D-6) ) RETURN
C
      DO J = ABS(NINT(2*(J1-J2))),NINT(2*(J1+J2)),2
         IF ( J.EQ.NINT(2*J3) ) GOTO 100
      END DO
      RETURN
C
C
 100  CONTINUE
      X = (2.0D0*J3+1.0D0)*RFACT(J1+J2-J3)*RFACT(J1-J2+J3)
     &    *RFACT(-J1+J2+J3)*RFACT(J1+M1)*RFACT(J1-M1)*RFACT(J2+M2)
     &    *RFACT(J2-M2)*RFACT(J3+M3)*RFACT(J3-M3)
C
      Y = RFACT(J1+J2+J3+1)
C
      VF = DSQRT(X/Y)
C
C
      N1 = NINT(J1+J2-J3)
      N2 = NINT(J1-M1)
      N3 = NINT(J2+M2)
      N4 = NINT(J3-J2+M1)
      N5 = NINT(J3-J1-M2)
      NTOP = MIN(N1,N2,N3)
      NBOT = MAX(0,-N4,-N5)
C
      N = NBOT + 1
      IF ( N.EQ.(2*(N/2)) ) THEN
         S = +1.0D0
      ELSE
         S = -1.0D0
      END IF
      SUM = 0.0D0
C
      DO N = NBOT,NTOP
         S = -S
         Y = FACT(N)*FACT(N1-N)*FACT(N2-N)*FACT(N3-N)*FACT(N4+N)
     &       *FACT(N5+N)
         SUM = SUM + (S/Y)
      END DO
      CGCRAC = VF*SUM
      END
C
C
C*==gauntc.f    processed by SPAG 6.05Rc at 15:27 on  7 Mar 2003
      FUNCTION GAUNTC(FACT,L1,M1,L2,M2,L3,M3)
C   ********************************************************************
C   *                                                                  *
C   *     GAUNT COEFFICIENTS for complex spherical harmonics  Y[l,m]   *
C   *                                                                  *
C   *            G = INT dr^  Y[l1,m1]* Y[l2,m2] Y[l3,m3]              *
C   *                                                                  *
C   * see: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM  Eq. (4.34)  *
C   *                                                                  *
C   * 26/01/95  HE                                                     *
C   ********************************************************************
C
      IMPLICIT NONE
C
C PARAMETER definitions
C
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793238462643D0)
C
C Dummy arguments
C
      INTEGER L1,L2,L3,M1,M2,M3
      DOUBLE PRECISION FACT(0:100)
      DOUBLE PRECISION GAUNTC
C
C Local variables
C
      DOUBLE PRECISION CGCRAC
      DOUBLE PRECISION DBLE
      DOUBLE PRECISION G
C
      IF ( (L1.LT.0) .OR. (L2.LT.0) .OR. (L3.LT.0) ) THEN
         G = 0.0D0
      ELSE
         G = (DBLE(2*L2+1)*DBLE(2*L3+1)/(4.0D0*PI*DBLE(2*L1+1)))
     &       **0.5D0*CGCRAC(FACT,DBLE(L3),DBLE(L2),DBLE(L1),DBLE(M3),
     &       DBLE(M2),DBLE(M1))*CGCRAC(FACT,DBLE(L3),DBLE(L2),DBLE(L1),
     &       0.0D0,0.0D0,0.0D0)
      END IF
      GAUNTC = G
      END
C
