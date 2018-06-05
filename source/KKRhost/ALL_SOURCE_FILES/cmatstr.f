C*==cmatstr.f    processed by SPAG 6.05Rc at 15:50 on 12 Oct 2002
      SUBROUTINE CMATSTR(STR,LSTR,A,N,M,MLIN,MCOL,IJQ,TOLP,K_FMT_FIL)
C   ********************************************************************
C   *                                                                  *
C   *   writes structure of COMPLEX   NxN   matrix   A                 *
C   *                                                                  *
C   *   M           is the actual array - size used for   A            *
C   *   MLIN/COL    MODE for line and column indexing                  *
C   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
C   *   TOL         tolerance for difference                           *
C   *   IJQ         if IJQ > 1000    pick  IQ-JQ-block matrix          *
C   *               assuming  IJQ = IQ*1000 + JQ                       *
C   *               else: no IQ-JQ-indexing                            *
C   *   K_FMT_FIL   output channel                                     *
C   *               a negative sign suppresses table at the end        *
C   *                                                                  *
C   *   any changes should be done in RMATSTR as well !!!!!!!!!!!!!!!  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE COMPLEX CI
      PARAMETER (CI=(0.0D0,1.0D0))
C
C Dummy arguments
C
      INTEGER IJQ,K_FMT_FIL,LSTR,M,MCOL,MLIN,N
      CHARACTER (len=LSTR) :: STR
      DOUBLE PRECISION TOLP
      DOUBLE COMPLEX A(M,M)
C
C Local variables
C
      DOUBLE COMPLEX B(N,N),CA,CB,ARG,DTAB(0:N*N)
      CHARACTER CHAR
      LOGICAL SAME,SMALL
      CHARACTER (len=1) :: CTAB(0:N*N),VZ(-1:+1)
      DOUBLE PRECISION DBLE
      CHARACTER (len=150) ::FMT1,FMT2,FMT3,FMT4
      INTEGER I,I1,IC0,ID,IL,ILSEP(20),IPT(218),IQ,ISL,IW(M),J,
     &        J0,JP,JQ,K,L3,LF,MM,N1,N2,N3,NC,ND,NFIL,NK,NM,NM1,NM2,NM3,
     &        NNON0,NSL
      INTEGER ICHAR,ISIGN,NINT
      DOUBLE PRECISION TOL
C
C*** End of declarations rewritten by SPAG
C
      DATA VZ/'-',' ',' '/
C
      SMALL(ARG) = ABS(ARG*TOL).LT.1.0D0
C
      SAME(CA,CB) = SMALL(1.0D0-CA/CB)
C
      NFIL = ABS(K_FMT_FIL)
C
      TOL = 1.0D0/TOLP
C
C----------------------------------------------- set block indices IQ JQ
C
      IF ( IJQ.GT.1000 ) THEN
         IQ = IJQ/1000
         JQ = IJQ - IQ*1000
         IF ( IQ*N.GT.M .OR. IQ*N.GT.M ) THEN
            WRITE (1337,99002) IJQ,IQ,JQ,IQ*N,JQ*N,N,M
            RETURN
         END IF
      ELSE
         IQ = 1
         JQ = 1
      END IF
C
C----------------------------------------------------- copy matrix block
C
      J0 = N*(JQ-1)
      DO J = 1,N
         I1 = N*(IQ-1)+1
         JP = J0 + J
         CALL ZCOPY(N,A(I1,JP),1,B(1,J),1)
      END DO
C
C------------------------------------------------ set up character table
C
      NC = 0
      DO I = 1,26
         NC = NC + 1
         IPT(NC) = 62 + I
      END DO
      DO I = 1,8
         NC = NC + 1
         IPT(NC) = 96 + I
      END DO
      DO I = 10,26
         NC = NC + 1
         IPT(NC) = 96 + I
      END DO
      DO I = 191,218
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 35,38
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 40,42
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 91,93
         NC = NC + 1
         IPT(NC) = I
      END DO
C
C---------------------------------------------------------------- header
      IC0 = ICHAR('0')
      N3 = N/100
      N2 = N/10 - N3*10
      N1 = N - N2*10 - N3*100
C
      IF ( N.LE.18 ) THEN
         FMT1 = '(8X,I3,''|'','
         FMT2 = '( 9X,''--|'','
         FMT3 = '( 9X,'' #|'','
         FMT4 = '( 9X,''  |'','
      ELSE
         FMT1 = '(   I4,''|'','
         FMT2 = '( 2X,''--|'','
         FMT3 = '( 2X,'' #|'','
         FMT4 = '( 2X,''  |'','
      END IF
C
      LF = 11
      L3 = 11
      IF ( MCOL.EQ.0 ) THEN
         FMT1 = FMT1(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'( 2A1),''|'',I3)'
         FMT2 = FMT2(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'(''--''),''|'',I3)'
         FMT3 = FMT3(1:LF)//'60(2X,I2))'
         FMT4 = FMT4(1:LF)//'60(I2,2X))'
         LF = 21
      ELSE
         IF ( MCOL.EQ.1 ) THEN
            NK = NINT(SQRT(DBLE(N)))
         ELSE IF ( MCOL.EQ.2 ) THEN
            NK = NINT(SQRT(DBLE(N/2)))
         ELSE IF ( MCOL.EQ.3 ) THEN
            NK = 2*NINT(SQRT(DBLE(N/2))) - 1
         END IF
         DO K = 1,NK
            IF ( MCOL.LE.2 ) THEN
               NM = 2*K - 1
            ELSE
               NM = 2*((K+1)/2)
            END IF
            NM2 = NM/10
            NM1 = NM - NM2*10
            NM3 = NM/2
            FMT1 = FMT1(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'( 2A1),''|'','
            FMT2 = FMT2(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'(''--''),''|'','
C
            IF ( MCOL.LE.2 ) THEN
               DO MM = 1,NM
                  IF ( MOD(MM,2).EQ.MOD(K,2) ) THEN
                     FMT3 = FMT3(1:L3)//'2X,'
                     FMT4 = FMT4(1:L3)//'I2,'
                  ELSE
                     FMT3 = FMT3(1:L3)//'I2,'
                     FMT4 = FMT4(1:L3)//'2X,'
                  END IF
                  L3 = L3 + 3
               END DO
               FMT3 = FMT3(1:L3)//'''|'','
               FMT4 = FMT4(1:L3)//'''|'','
               L3 = L3 + 4
            ELSE
               FMT3 = FMT3(1:LF)//CHAR(IC0+NM3)//'(2X,I2),''|'','
               FMT4 = FMT4(1:LF)//CHAR(IC0+NM3)//'(I2,2X),''|'','
               L3 = L3 + 13
            END IF
            LF = LF + 13
         END DO
         IF ( MCOL.EQ.2 ) THEN
            FMT1 = FMT1(1:LF)//FMT1(12:LF)
            FMT2 = FMT2(1:LF)//FMT2(12:LF)
            FMT3 = FMT3(1:L3)//FMT3(12:L3)
            FMT4 = FMT4(1:L3)//FMT4(12:L3)
            LF = 2*LF - 11
            L3 = 2*L3 - 11
         END IF
         FMT1 = FMT1(1:LF)//'I3)'
         FMT2 = FMT2(1:LF)//'I3)'
         FMT3 = FMT3(1:L3)//'I3)'
         FMT4 = FMT4(1:L3)//'I3)'
      END IF
      IF ( MLIN.EQ.0 ) THEN
         NSL = 1
         ILSEP(1) = N
      ELSE IF ( MLIN.EQ.1 ) THEN
         NSL = NINT(SQRT(DBLE(N)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
      ELSE IF ( MLIN.EQ.2 ) THEN
         NSL = NINT(SQRT(DBLE(N/2)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
         DO IL = 1,NSL
            ILSEP(NSL+IL) = ILSEP(NSL) + IL**2
         END DO
         NSL = 2*NSL
      ELSE IF ( MLIN.EQ.3 ) THEN
         NSL = 2*NINT(SQRT(DBLE(N/2))) - 1
         ILSEP(1) = 2
         DO K = 2,NSL
            ILSEP(K) = ILSEP(K-1) + 2*((K+1)/2)
         END DO
      END IF
C
C
      WRITE (NFIL,99001) STR(1:LSTR)
      IF ( IJQ.GT.1000 ) WRITE (NFIL,99003) IQ,JQ
      WRITE (NFIL,FMT3) (I,I=2,N,2)
      WRITE (NFIL,FMT4) (I,I=1,N,2)
      WRITE (NFIL,FMT=FMT2)
C------------------------------------------------------------ header end
      NNON0 = 0
      ND = 0
      CTAB(0) = ' '
      DTAB(0) = 9999D0
C
      DO I = 1,N
         DO J = 1,N
            IF ( .NOT.SMALL(B(I,J)) ) THEN
               NNON0 = NNON0 + 1
               DO ID = 1,ND
                  IF ( SAME(B(I,J),+DTAB(ID)) ) THEN
                     IW(J) = +ID
                     GOTO 50
                  END IF
                  IF ( SAME(B(I,J),-DTAB(ID)) ) THEN
                     IW(J) = -ID
                     GOTO 50
                  END IF
               END DO
C----------------------------------------------------------- new element
               ND = ND + 1
               IW(J) = ND
               DTAB(ND) = B(I,J)
               IF ( ABS(DTAB(ND)-1.0D0)*TOL.LT.1.0D0 ) THEN
                  CTAB(ND) = '1'
               ELSE IF ( ABS(DTAB(ND)+1.0D0)*TOL.LT.1.0D0 ) THEN
                  DTAB(ND) = +1.0D0
                  CTAB(ND) = '1'
                  IW(J) = -ND
               ELSE IF ( ABS(DTAB(ND)-CI)*TOL.LT.1.0D0 ) THEN
                  CTAB(ND) = 'i'
               ELSE IF ( ABS(DTAB(ND)+CI)*TOL.LT.1.0D0 ) THEN
                  DTAB(ND) = +CI
                  CTAB(ND) = 'i'
                  IW(J) = -ND
               ELSE
                  CTAB(ND) = CHAR(IPT(1+MOD((ND+1),NC)))
               END IF
            ELSE
               IW(J) = 0
            END IF
 50      END DO
C------------------------------------------------------------ write line
         WRITE (NFIL,FMT=FMT1) I,
     &                         (VZ(ISIGN(1,IW(J))),CTAB(ABS(IW(J))),J=1,
     &                         N),I
C
         DO ISL = 1,NSL
            IF ( I.EQ.ILSEP(ISL) ) WRITE (NFIL,FMT=FMT2)
         END DO
      END DO
C
C------------------------------------------------------------------ foot
C
      WRITE (NFIL,FMT4) (I,I=1,N,2)
      WRITE (NFIL,FMT3) (I,I=2,N,2)
C
      IF ( K_FMT_FIL.GT.0 ) THEN
         WRITE (NFIL,99004) (ID,CTAB(ID),DTAB(ID),ID=1,ND)
         WRITE (NFIL,99005) NNON0,TOLP,N*N - NNON0,TOLP
      ELSE
         WRITE (NFIL,*) ' '
      END IF
C
99001 FORMAT (/,8X,A,/)
99002 FORMAT (/,1X,79('*'),/,10X,'inconsistent call of <CMATSTR>',/,10X,
     &        'argument IJQ =',I8,'  implies IQ=',I3,'   JQ=',I3,/,10X,
     &        'IQ*N=',I6,' > M   or   JQ*N=',I6,' > M   for N =',I4,
     &        ' M=',I4,/,1X,79('*'),/)
99003 FORMAT (8X,'IQ-JQ-block  for  IQ = ',I3,'   JQ = ',I3,/)
99004 FORMAT (/,8X,'symbols used:',/,(8X,I3,3X,A1,2X,2F20.12))
99005 FORMAT (/,8X,I5,' elements   >',1PE9.1,/,
     &          8X,I5,' elements   <',1PE9.1,/)
      END
