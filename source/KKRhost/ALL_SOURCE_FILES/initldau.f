C*==initldau.f    processed by SPAG 6.05Rc at 15:27 on  7 Mar 2003
      SUBROUTINE INITLDAU(NSRA,NTLDAU,ITLDAU,LOPT,UEFF,JEFF,EREFLDAU,
     &                    VISP,NSPIN,R,DRDI,Z,IPAN,IRCUT,PHI,ULDAU)
C
C   *******************************************************************
C   *  Calculates ULDAU                                               *
C   *  fivos will add some comments later                             *
C   *                                                                 *
C   *  Munich, March 2003, h.ebert, v.popescu and ph.mavropoulos      *
C   *                                                                 *
C   *  It is already later, but no comments have been added yet.      *
C   *  But wait up, they're coming...    Munich, Feb.2004, Phivos     *
C   *                                                                 *
C   *******************************************************************
C
      IMPLICIT NONE
      INCLUDE 'inc.p'
C
C PARAMETER definitions
C
      INTEGER NPOTD,MMAXD
      PARAMETER (NPOTD=(2*KREL+(1-KREL)*NSPIND)*NATYPD)
      PARAMETER (MMAXD=2*LMAXD+1)
C
C Dummy arguments
C
      INTEGER NTLDAU,NSPIN,NSRA
      DOUBLE PRECISION DRDI(IRMD,NATYPD),R(IRMD,NATYPD),
     &                 VISP(IRMD,NPOTD),Z(NATYPD)
!      &                ,ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD)
      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE COMPLEX PHI(IRMD,NATYPD)
      INTEGER ITLDAU(NATYPD),LOPT(NATYPD)
      INTEGER IPAN(NATYPD),IRCUT(0:IPAND,NATYPD)
C
C Local variables
C
      DOUBLE PRECISION AA(MMAXD,MMAXD,MMAXD,MMAXD,0:2*LMAXD),
     &                 EREFLDAU(NATYPD),FACT(0:100),FCLMB(0:2*LMAXD+1),
     &                 G12,G34,JEFF(NATYPD),PI,RLOP,
     &                 RPW(IRMD,2*LMAXD+1),SCL,SG(IRMD),SL(IRMD),
     &                 SUM,SUMFCLMB,TG(IRMD),TL(IRMD),UEFF(NATYPD),
     &                 WGTFCLMB(0:2*LMAXD+1),WIG3J,WINT(IRMD),W2(IRMD)
      DOUBLE PRECISION CGCRAC,GAUNTC
      DOUBLE PRECISION ATAN,DBLE
      INTEGER I1,IM1,IM2,IM3,IM4,IPAN1,IR,IRC1,IT,KK,
     &        L1,LF,LFMAX,LL,M1,M2,M3,M4
      INTEGER NINT
C     ..


      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

      PI = 4.D0*ATAN(1.0D0)
      FACT(0) = 1.0D0
      DO I1 = 1,100
         FACT(I1) = FACT(I1-1)*DBLE(I1) ! factorial
      END DO
      WRITE (1337,'(/,79(1H=),/,22X,A,/,79(1H=))') 
     &                            'LDA+U:  INITIALISE Coulomb matrix U'
C
C -> Calculate test functions Phi. Phi is already normalised to
C    int phi**2 dr =1, thus it also contains factor r.
C
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                        Loop over atoms
C                                         which need LDA+U ( LOPT >= 0 )
      DO IT = 1,NTLDAU
         I1 = ITLDAU(IT)
         IF ( LOPT(I1)+1.EQ.0 ) STOP ' this atom should be LDA+U'
         CALL PHICALC(I1,LOPT(I1),VISP,IPAN,IRCUT,R,DRDI,Z,
     &                EREFLDAU(I1),PHI(1,I1),NSPIN,NSRA)
      ENDDO
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      WRITE (1337,'(6X,43(1H-),/,6X,A,/,6X,43(1H-))')
     &     'Slater integrals F^n'
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      DO IT = 1,NTLDAU
         I1 = ITLDAU(IT)
         IPAN1 = IPAN(I1)
         IRC1 = IRCUT(IPAN1,I1)
C
C -> define r**l in array rpw:
C
         LFMAX = 2*LOPT(I1)
         DO IR = 2,IRMD
            RPW(IR,1) = R(IR,I1)
            DO L1 = 2,2*LMAXD + 1
               RPW(IR,L1) = R(IR,I1)*RPW(IR,L1-1)
            END DO
         END DO
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C 
C 1.  Calculate slater integrals FCLMB
C     (here using only the large component of sra wavefunction.
C     Whoever wants can also do the small component.)
C
C 1a. Calculate slater integrals
C
         RLOP = DBLE(LOPT(I1))
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
               WINT(IR) = DREAL( DCONJG(PHI(IR,I1)) * PHI(IR,I1) )
               W2(IR) = 2.D0 * DRDI(IR,I1) * WINT(IR)
               TL(IR) = W2(IR) * RPW(IR,LF)
               TG(IR) = W2(IR) / RPW(IR,LF+1)
            END DO
C
            CALL SOUTK(TL,SL,IPAN(I1),IRCUT(0,I1))
            CALL SOUTK(TG,SG,IPAN(I1),IRCUT(0,I1))
C
            SL(1) = 0.0D0
            DO IR = 2,IRC1
               SL(IR) = SL(IR)/RPW(IR,LF+1)
     &                   + (SG(IRC1)-SG(IR))*RPW(IR,LF)
            END DO
C
            SG(1) = 0.0D0
C
C See Note on integrand above.
C
            DO IR = 2,IRC1
               SG(IR) = WINT(IR) * SL(IR)
            END DO
C
            CALL SIMPK(SG,FCLMB(LF),IPAN1,IRCUT(0,I1),DRDI(1,I1))
C
            WIG3J = (-1)**NINT(RLOP) * (1D0/SQRT(2D0*RLOP+1D0))
     &                 * CGCRAC(FACT,RLOP,DBLE(LF),RLOP,0D0,0D0,0D0)
C
            WGTFCLMB(LF) = ((2*RLOP+1)/(2*RLOP))*WIG3J**2
            SUMFCLMB = SUMFCLMB + WGTFCLMB(LF)*FCLMB(LF)
         END DO
C ----------------------------------------------------------------------
C
C 1b.   Normalise slater integrals FCLMB
C
         SCL = JEFF(I1)/SUMFCLMB
         FCLMB(0) = UEFF(I1)
C
         DO LF = 2,LFMAX,2
            FCLMB(LF) = SCL*FCLMB(LF)
         END DO
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C ======================================================================
C
C 2.   Calculate coefficient matrix AA = a(m1,m2,m3,m4)
C     
         CALL RINIT(MMAXD*MMAXD*MMAXD*MMAXD*(2*LMAXD+1),AA)
         LL = LOPT(I1)
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
                        END DO
C
                        AA(IM1,IM2,IM3,IM4,LF)
     &                       = SUM*4.D0*PI/(2.D0*DBLE(LF)+1.D0)
                     END IF
                  END DO
               END DO
            END DO
         END DO
C ======================================================================
C
C UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C 3.  Calculate ULDAU
C
         CALL RINIT(MMAXD*MMAXD*MMAXD*MMAXD,ULDAU(1,1,1,1,I1))
C
         DO LF = 0,LFMAX,2
            DO IM4 = 1,2*LL + 1
               DO IM3 = 1,2*LL + 1
                  DO IM2 = 1,2*LL + 1
                     DO IM1 = 1,2*LL + 1
                        ULDAU(IM1,IM2,IM3,IM4,I1)
     &                       = ULDAU(IM1,IM2,IM3,IM4,I1)
     &                       + AA(IM1,IM2,IM3,IM4,LF)*FCLMB(LF)
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
         WRITE (1337,'(/,8X,A,I3,/)') 'ATOM: ',I1
         WRITE (1337,'(12X,A,F8.4,A)') 'LDA+U reference energy :',
     &        EREFLDAU(I1),' Ry'
         WRITE (1337,'(12X,A,2F8.4,A,/)') 'Ueff and Jeff = ',UEFF(I1),
     &        JEFF(I1),' Ry'
         WRITE (1337,'(12X,"Scaling factor for F^n :",F10.6,/)') SCL
         WRITE (1337,'(12X,"  n   F^n calculated   F^n scaled ")')
         DO LF = 2,LFMAX,2
            WRITE (1337,'(12X,I3,2(2X,F12.8," Ry"))') 
     &           LF,FCLMB(LF)/SCL,FCLMB(LF)
         END DO
         IF ( IT.LT.NTLDAU ) WRITE(1337,'(8X,58(1H~))')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
      END DO
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      WRITE (1337,'(/,6X,60(1H-),/,6X,A,/,6X,60(1H-))')
     &     'Coulomb matrix U(m1,m1,m3,m3)'
      DO IT = 1,NTLDAU
         I1 = ITLDAU(IT)
         LL = LOPT(I1)
         LL = MIN(3,LL)
         WRITE(1337,'(/,8X,"ATOM :",I3,/)') I1
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
         DO IM1 = 1,2*LL + 1
            WRITE (1337,99001) (ULDAU(IM1,IM1,IM3,IM3,I1),IM3=1,2*LL+1)
         END DO
         IF ( IT.LT.NTLDAU ) WRITE(1337,'(8X,58(1H~))')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
      END DO
      WRITE (1337,'(/,6X,60(1H-),/)')
99001 FORMAT(6X,7F10.6)
      END
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
