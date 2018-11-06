C*==initldau.f    processed by SPAG 6.05Rc at 15:27 on  7 Mar 2003
      MODULE mod_initldau
      CONTAINS

!-------------------------------------------------------------------------------
!> Summary: Calcualtes ULDAU 
!> Author: Who wrote this subroutine
!> Category: LDA+u, initialization, kkrimp
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------
!> @note  NOTES 
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!-------------------------------------------------------------------------------



      SUBROUTINE INITLDAU(LMAXD,NATYP,NSPIN,VPOT,ZATOM,KVREL,CELL,LDAU)
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
C
C PARAMETER definitions
C
#ifdef CPP_MPI
      use mpi
#endif
      USE NRTYPE
      USE TYPE_LDAU 
      USE TYPE_CELL
      USE mod_simpk
      USE mod_soutk
      USE MOD_RWLDAUPOT
      USE MOD_PHICALC
      use mod_version_info
      IMPLICIT NONE

      TYPE(LDAU_TYPE),ALLOCATABLE    :: LDAU(:)  ! lda+u variables, intended dimension: (natyp) ! lda+u
      TYPE(CELL_TYPE)                :: CELL(:)
      INTEGER NATYP,MMAX,LMAXD,IRMD
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION DZERO
      PARAMETER (DZERO=0.0D0)
C
C Dummy arguments
C
      INTEGER NSPIN,KVREL
!      DOUBLE PRECISION DRDI(:,:),R(:,:),   ,RWS(:),RMT(:)! DRDI(IRMD,NATYP),R(IRMD,NATYP)
      DOUBLE PRECISION VPOT(:,:,:,:),                    ! (IR,LM,ISPIN,IATOM)
     &     ZATOM(:)  ! (NATYP)
!      INTEGER IPAN(:),IRCUT(:,:)! IPAN(NATYP),IRCUT(0:IPAND,NATYP)
C
C Local variables
C
      DOUBLE PRECISION,ALLOCATABLE :: AA(:,:,:,:,:),
     &     FCLMB(:),RPW(:,:),RS(:,:),S(:),
     &     SG(:),SL(:),TG(:),TL(:),
     &     WGTFCLMB(:),WINT(:),W2(:),VSPHER(:,:)
      DOUBLE PRECISION FACT(0:100),
     &                 G12,G34,RFPI,RLOP,SCL,SUM,SUMFCLMB,WIG3J
!      DOUBLE PRECISION CGCRAC,GAUNTC1
      DOUBLE PRECISION DBLE,DSQRT
      INTEGER I1,IM1,IM2,IM3,IM4,IR,IRC1,IRS1,IRUNLDAU,KK
     &        ,L1,LF,LFMAX,LL,M1,M2,M3,M4,NSRA,IAT,ISPIN,IPOT,IPAN1
      INTEGER NINT,IERROR


!      EXTERNAL CGCRAC,GAUNTC1
C
      MMAX = 2*LMAXD + 1
      IRMD = CELL(1)%NRMAXD

C
      OPEN (67,FILE='ldaupot',FORM='FORMATTED') 
      call version_check_header(67)
      READ (67,*) IRUNLDAU
      CLOSE (67)
C
C  Calculates U matrix, or reads from disk. Remove this IF statement to have
C  U calculated anew.
C
      IF ( IRUNLDAU.NE.0 ) THEN ! read from disk and exit
         CALL RWLDAUPOT(.FALSE.,4,NATYP,NSPIN,LMAXD,IRMD,IRUNLDAU,LDAU)
         RETURN
      ENDIF

! Allocate arrays
      ALLOCATE( AA(MMAX,MMAX,MMAX,MMAX,0:2*LMAXD) )
      ALLOCATE( FCLMB(0:2*LMAXD+1), RPW(IRMD,2*LMAXD+1) ) 
      ALLOCATE( RS(IRMD,0:LMAXD), S(0:LMAXD), SG(IRMD), SL(IRMD) ) 
      ALLOCATE( TG(IRMD), TL(IRMD), WGTFCLMB(0:2*LMAXD+1) ) 
      ALLOCATE( WINT(IRMD), W2(IRMD), VSPHER(IRMD,NATYP*NSPIN) ) 

! Copy spherical part of potential
      DO IAT = 1,NATYP
         DO ISPIN = 1,NSPIN
            IPOT = (IAT-1)*NSPIN + ISPIN
            VSPHER(:,IPOT) = VPOT(:,1,ISPIN,IAT)
         ENDDO
      ENDDO


C Read from ldaupot file IRUNLDAU,LOPT,UEFF,JEFF,EREFLDAU. 
! Arrays are allocated in this routine
      CALL RWLDAUPOT(.FALSE.,1,NATYP,NSPIN,LMAXD,IRMD,IRUNLDAU,LDAU)

#ifdef CPP_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
#endif

! Initialize arrays, calculate matrix ULDAU

      DO I1 = 1,NATYP
          IF ( LDAU(I1)%LOPT.GE.0 ) THEN
             LDAU(I1)%WLDAU(:,:,:) = 0.D0
             LDAU(I1)%ULDAU(:,:,:,:) = 0.D0
             LDAU(I1)%PHI(:) = CZERO
             DO IR = 1,IRMD
                LDAU(I1)%CUTOFF(IR) = 1.D0 /                            ! Fermi-like function
     &        ( 1.D0 + EXP( (CELL(I1)%RMESH(IR)-CELL(I1)%RMT)/5.D-03 ) )
             ENDDO
          ENDIF
       ENDDO

C Calculate test functions Phi. Phi is already normalised to
C int phi**2 dr =1, thus it also contains factor r.
      DO I1 = 1,NATYP
          IF ( LDAU(I1)%LOPT.GE.0 ) THEN
c Specify nsra
              NSRA = 1                  ! Schroedinger
              IF (KVREL.GE.1) NSRA = 2  ! Scalar-relativistic 
              CALL PHICALC(LDAU(I1)%LOPT,LDAU(I1)%PHI,VSPHER
     &             ,CELL(I1)%NPAN,CELL(I1)%NRCUT,CELL(I1)%RMESH
     &             ,CELL(I1)%DRMESHDI,ZATOM(I1)
     &             ,LDAU(I1)%EREFLDAU,0
     &             ,LDAU(I1)%WLDAUAV,LDAU(I1)%CUTOFF
     &             ,I1,NSPIN,NSRA,LMAXD,IRMD)
          ENDIF
      ENDDO
C
!     PI = 4.D0*DATAN(1.0D0)
      RFPI = 2.D0*DSQRT(PI)
C
      FACT(0) = 1.0D0
      DO I1 = 1,100
         FACT(I1) = FACT(I1-1)*DBLE(I1) ! factorial
      END DO
      WRITE (6,*)
      WRITE (6,'(79(1H=))')
      WRITE (6,'(A)') '           LDA+U  INITIALISE Coulomb matrix U'
      WRITE (6,*)
      WRITE (6,'(79(1H=))')
C
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                       Loop over atoms...which need LDA+U ( LOPT >= 0 )
C
      DO I1 = 1,NATYP
C
C                                     ...
C
         IF ( LDAU(I1)%LOPT.GE.0 ) THEN

            LL = LDAU(I1)%LOPT

C
            IPAN1 = CELL(I1)%NPAN
            IRS1 = CELL(I1)%NRCUT(1)
            IRC1 = CELL(I1)%NRCUT(IPAN1)
C
C     define r**l in array rpw:
C
            LFMAX = 2*LL
            DO IR = 2,IRMD
               RPW(IR,1) = CELL(I1)%RMESH(IR)
               DO L1 = 2,2*LMAXD + 1
                  RPW(IR,L1) = CELL(I1)%RMESH(IR)*RPW(IR,L1-1)
               END DO
            END DO
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C 
C 1.  Calculate slater integrals FCLMB
C     (here using only the large component of sra wavefunction.
C     Whoever wants can also do the small component.)
C
C 1a. Calculate slater integrals
            RLOP = DBLE(LL)
            SUMFCLMB = 0.D0
            DO LF = 2,LFMAX,2
               TL(1) = 0.0D0
               TG(1) = 0.0D0
C Note on integrand:
C Integrals are up to IRC1 because we integrate in sphere, without thetas.
C In case of cell integration, from IRS1+1 to IRC1 a convolution with
C thetas and gaunts is needed: 
C Int dr R_l(r)**2 Sum_L' Gaunt_{lm,lm,l'm'}*thetas_{l'm'}(r).
C But then the result is m-dependent. Here, the easy way is used!
               DO IR = 2,IRC1
                   WINT(IR) = DREAL( 
     &                     DCONJG(LDAU(I1)%PHI(IR)) * LDAU(I1)%PHI(IR) )
                   W2(IR) = 2.D0 * CELL(I1)%DRMESHDI(IR) * WINT(IR)
                   TL(IR) = W2(IR) * RPW(IR,LF)
                   TG(IR) = W2(IR) / RPW(IR,LF+1)
               END DO

C
               CALL SOUTK(TL,SL,CELL(I1)%NPAN,CELL(I1)%NRCUT(0))
               CALL SOUTK(TG,SG,CELL(I1)%NPAN,CELL(I1)%NRCUT(0))
C
               SL(1) = 0.0D0
               DO IR = 2,IRC1
                  SL(IR) = SL(IR)/RPW(IR,LF+1) + (SG(IRC1)-SG(IR))
     &                     *RPW(IR,LF)
               END DO
C
               SG(1) = 0.0D0
C See Note on integrand above.
               DO IR = 2,IRC1
                  SG(IR) = WINT(IR) * SL(IR)
               END DO
C
               CALL SIMPK(SG,FCLMB(LF),IPAN1,CELL(I1)%NRCUT(0)
     &                                ,CELL(I1)%DRMESHDI,IPAN1)
C
               WIG3J = (-1)**NINT(RLOP)*(1D0/SQRT(2D0*RLOP+1D0))
     &                 *CGCRAC(FACT,RLOP,DBLE(LF),RLOP,0D0,0D0,0D0)
C
               WGTFCLMB(LF) = ((2*RLOP+1)/(2*RLOP))*WIG3J**2
               SUMFCLMB = SUMFCLMB + WGTFCLMB(LF)*FCLMB(LF)
            END DO
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C
C 1b.   Normalise slater integrals FCLMB
C
            IF (DABS(SUMFCLMB).GT.1.D-10) THEN
               SCL = LDAU(I1)%JEFF/SUMFCLMB
            ELSE
               SCL = 1.D-10
            ENDIF
            FCLMB(0) = LDAU(I1)%UEFF
C
            DO LF = 2,LFMAX,2
               FCLMB(LF) = SCL*FCLMB(LF)
            END DO
C
C
C =====================================================================
C
C 2.   Calculate coefficient matrix AA = a(m1,m2,m3,m4)
C
            AA(:,:,:,:,:) = 0.D0
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
                           DO KK = -LF,LF
                              G12 = GAUNTC1(FACT,LL,M1,LF,KK,LL,M2)
                              G34 = GAUNTC1(FACT,LL,M3,LF,-KK,LL,M4)
                              SUM = SUM + G12*G34*(-1)**ABS(KK)
                           END DO
                           AA(IM1,IM2,IM3,IM4,LF)
     &                        = SUM*4.D0*PI/(2.D0*DBLE(LF)+1.D0)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
C
C =====================================================================
C
C
C UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C 3.  Calculate ULDAU
C
            LDAU(I1)%ULDAU(:,:,:,:) = 0.D0
            DO LF = 0,LFMAX,2
               DO IM4 = 1,2*LL + 1
                  DO IM3 = 1,2*LL + 1
                     DO IM2 = 1,2*LL + 1
                        DO IM1 = 1,2*LL + 1
                           LDAU(I1)%ULDAU(IM1,IM2,IM3,IM4)
     &                   = LDAU(I1)%ULDAU(IM1,IM2,IM3,IM4)
     &                        + AA(IM1,IM2,IM3,IM4,LF)*FCLMB(LF)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
C
C UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C
            WRITE (6,*)
            WRITE (6,'(A,i3)') ' ATOM ',I1
            WRITE (6,*)
            WRITE (6,'(4X,A,f8.4,A)') 'LDA+U reference energy ',
     &                                LDAU(I1)%EREFLDAU,' Ry'
            WRITE (6,'(4X,A,2f8.4,A)') 'Ueff and Jeff = ',LDAU(I1)%UEFF,
     &                                 LDAU(I1)%JEFF,' Ry'
            WRITE (6,*)
            WRITE (6,'(4X,"Scaling factor for F^n :",f10.6)') SCL
            WRITE (6,*)
            WRITE (6,'(4X,"  n   F^n calculated   F^n scaled ")')
            DO LF = 2,LFMAX,2
               WRITE (6,'(4X,i3,2(2X,f12.8," Ry"))') LF,FCLMB(LF)/SCL,
     &                FCLMB(LF)
            END DO
C
            WRITE (6,*)
            WRITE (6,'(" Coulomb matrix U(m1,m1,m3,m3) for atom ",i3)') 
     &            I1
            WRITE (6,*)
            DO IM1 = 1,2*LL + 1
               WRITE (6,99001) (LDAU(I1)%ULDAU(IM1,IM1,IM3,IM3),
     &                          IM3=1,2*LL+1)
            END DO
            WRITE (6,'(79(1H=))')
C


         END IF                 ! (LDAU(I1)%LOPT.GE.0)
      END DO                    !  I1 = 1,NATYP

C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

      CALL RWLDAUPOT(.TRUE.,4,NATYP,NSPIN,LMAXD,IRMD,IRUNLDAU,LDAU)
C
99001 FORMAT (5F12.8)
C
! Deallocate arrays
      DEALLOCATE( AA, FCLMB, RPW, RS, S, SG, SL, TG, TL ) 
      DEALLOCATE( WGTFCLMB, WINT, W2, VSPHER ) 

      RETURN
      END SUBROUTINE INITLDAU

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
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
      REAL*8 J1,J2,J3,M1,M2,M3
      REAL*8 CGCRAC
      REAL*8 FACT(0:100)
C
C Local variables
C
      DOUBLE PRECISION DSQRT
      INTEGER J,N,N1,N2,N3,N4,N5,NBOT,NTOP
      INTEGER NINT
      REAL*8 RFACT
      REAL*8 S,SUM,VF,X,Y
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
      END FUNCTION CGCRAC
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
C*==gauntc.f    processed by SPAG 6.05Rc at 15:27 on  7 Mar 2003
      FUNCTION GAUNTC1(FACT,L1,M1,L2,M2,L3,M3)
C   ********************************************************************
C   *                                                                  *
C   *     GAUNT COEFFICIENTS for complex spherical harmonics  Y[l,m]   *
C   *                                                                  *
C   *            G = INT d^3r  Y[l1,m1]* Y[l2,m2] Y[l3,m3]              *
C   *                                                                  *
C   * see: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM  Eq. (4.34)  *
C   *                                                                  *
C   * 26/01/95  HE                                                     *
C   ********************************************************************
C
      use nrtype, only: pi
      IMPLICIT NONE
C
C PARAMETER definitions
C
!      REAL*8 PI
!      PARAMETER (PI=3.141592653589793238462643D0)
C
C Dummy arguments
C
      INTEGER L1,L2,L3,M1,M2,M3
      REAL*8 FACT(0:100)
      REAL*8 GAUNTC1
C
C Local variables
C
!      REAL*8 CGCRAC
      DOUBLE PRECISION DBLE
      REAL*8 G
!      EXTERNAL CGCRAC
C
      IF ( (L1.LT.0) .OR. (L2.LT.0) .OR. (L3.LT.0) ) THEN
         G = 0.0D0
      ELSE
         G = (DBLE(2*L2+1)*DBLE(2*L3+1)/(4.0D0*PI*DBLE(2*L1+1)))
     &       **0.5D0*CGCRAC(FACT,DBLE(L3),DBLE(L2),DBLE(L1),DBLE(M3),
     &       DBLE(M2),DBLE(M1))*CGCRAC(FACT,DBLE(L3),DBLE(L2),DBLE(L1),
     &       0.0D0,0.0D0,0.0D0)
      END IF
      GAUNTC1 = G
      END FUNCTION GAUNTC1


      END MODULE mod_initldau 
