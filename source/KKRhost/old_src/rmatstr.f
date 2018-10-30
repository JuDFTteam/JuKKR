      SUBROUTINE RMATSTR(STR,LSTR,A,N,M,MLIN,MCOL,TOLP,NFIL)
C   ********************************************************************
C   *                                                                  *
C   *   writes structure of REAL      NxN   matrix   A                 *
C   *                                                                  *
C   *   M           is the actual array - size used for   A            *
C   *   MLIN/COL    MODE for line and comlun indexing                  *
C   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
C   *   TOL         tolerance for difference                           *
C   *                                                                  *
C   *                                                         03/01/96 *
C   ********************************************************************
      IMPLICIT NONE
 
      INTEGER, PARAMETER :: NDIFMAX =  250
      INTEGER, intent(in) :: lstr, m, n, mlin, mcol, nfil
      double precision, intent(in) :: tolp

      CHARACTER (len=LSTR) :: STR
      CHARACTER (len=150) :: FMT1,FMT2,FMT3,FMT4
      CHARACTER (len=1) :: CTAB(0:NDIFMAX), VZ(-1:+1)
      INTEGER    IW(M), ILSEP(20)
      INTEGER :: ICAUC , ICZUC, NATOZ, ICALC, ICILC, K, NK, NM, NM2, 
     &           NM1, NM3, IC0, N3, N2, N1, LF, L3, MM, NSL, IL, I,
     &           NNON0, ND, NALF, J, ID, ICND, ISL
      double precision :: tol
      DOUBLE PRECISION     A(M,M), CNUM,CA,CB,DTAB(0:NDIFMAX)  
      LOGICAL  CSMALL, CSAME

      SAVE VZ
      DATA VZ / '-', ' ', ' ' /
 
      CSMALL(CNUM) =   ABS(CNUM*TOL) .LT. 1.0D0

      CSAME(CA,CB) =   CSMALL(1.0D0 - CA/CB)

      TOL = 1.0D0 / TOLP

      ICAUC = ICHAR('A')
      ICZUC = ICHAR('Z')
      NATOZ = ICZUC - ICAUC + 1
      ICALC = ICHAR('a')
      ICILC = ICHAR('i')
        
c---------------------------------------------------------------- header
      IC0 = ICHAR('0')
      N3=N/100
      N2=N/10 - N3*10
      N1=N    - N2*10 - N3*100 

      FMT1='(8X,I3,''|'','
      FMT2='( 9X,''--|'','
      FMT3='( 9X,'' #|'','
      FMT4='( 9X,''  |'','

      LF  = 11
      L3  = 11
      IF( MCOL .EQ. 0 ) THEN
         FMT1=FMT1(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)
     &                  //CHAR(IC0+N1)//'( 2A1),''|'',I3)'
         FMT2=FMT2(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)
     &                  //CHAR(IC0+N1)//'(''--''),''|'',I3)'
         FMT3=FMT3(1:LF)//'60(2X,I2))'
         FMT4=FMT4(1:LF)//'60(I2,2X))'
         LF = 21
      ELSE 
         IF( MCOL .EQ. 1 ) THEN
            NK = NINT(SQRT( DBLE(N) )) 
         ELSE IF( MCOL .EQ. 2 ) THEN         
            NK = NINT(SQRT( DBLE(N/2) )) 
         ELSE IF( MCOL .EQ. 3 ) THEN         
            NK = 2*NINT(SQRT( DBLE(N/2) )) - 1
         END IF
         DO K=1,NK
            IF( MCOL .LE. 2 ) THEN
               NM = 2*K-1
            ELSE
               NM = 2*((K+1)/2)
            END IF
            NM2=NM/10
            NM1=NM - NM2*10 
            NM3=NM/2 
            FMT1=FMT1(1:LF)//CHAR(IC0+NM2)
     &                     //CHAR(IC0+NM1)//'( 2A1),''|'','
            FMT2=FMT2(1:LF)//CHAR(IC0+NM2)
     &                     //CHAR(IC0+NM1)//'(''--''),''|'','
            
            IF( MCOL .LE. 2 ) THEN
               DO MM=1,NM
                  IF( MOD(MM,2) .EQ. MOD(K,2) ) THEN
                     FMT3=FMT3(1:L3)//'2X,'
                     FMT4=FMT4(1:L3)//'I2,'
                  ELSE
                     FMT3=FMT3(1:L3)//'I2,'
                     FMT4=FMT4(1:L3)//'2X,'
                  END IF
                  L3=L3+3
               END DO
               FMT3=FMT3(1:L3)//'''|'','
               FMT4=FMT4(1:L3)//'''|'','
               L3=L3+4
            ELSE
               FMT3=FMT3(1:LF)//CHAR(IC0+NM3)//'(2X,I2),''|'','
               FMT4=FMT4(1:LF)//CHAR(IC0+NM3)//'(I2,2X),''|'','
               L3=L3+13
            END IF
            LF  = LF + 13
         END DO
         IF( MCOL .EQ. 2 ) THEN
            FMT1=FMT1(1:LF)//FMT1(12:LF)
            FMT2=FMT2(1:LF)//FMT2(12:LF)
            FMT3=FMT3(1:L3)//FMT4(12:L3)
            FMT4=FMT4(1:L3)//FMT3(12:L3)
            LF  = 2*LF - 11
         END IF
         FMT1=FMT1(1:LF)//'I3)'
         FMT2=FMT2(1:LF)//'I3)'
         FMT3=FMT3(1:L3)//'I3)'
         FMT4=FMT4(1:L3)//'I3)'
      END IF
      IF( MLIN .EQ. 0 ) THEN
         NSL = 1
         ILSEP(1) = N
      ELSE IF( MLIN .EQ. 1 ) THEN
         NSL = NINT(SQRT( DBLE(N) ))
         DO IL=1,NSL
            ILSEP(IL) = IL**2
         END DO
      ELSE IF( MLIN .EQ. 2 ) THEN
         NSL = NINT(SQRT( DBLE(N/2) ))
         DO IL=1,NSL
            ILSEP(IL) = IL**2
         END DO
         DO IL=1,NSL
            ILSEP(NSL+IL) = ILSEP(NSL) + IL**2
         END DO
         NSL = 2*NSL
      ELSE IF( MLIN .EQ. 3 ) THEN
         NSL = 2*NINT(SQRT( DBLE(N/2) )) - 1
         ILSEP(1) = 2
         DO K=2,NSL
            ILSEP(K) = ILSEP(K-1) + 2*((K+1)/2)
         END DO
      END IF


      WRITE(NFIL,9000) STR(1:LSTR)
      WRITE(NFIL,FMT3) (I,I=2,N,2)
      WRITE(NFIL,FMT4) (I,I=1,N,2)
      WRITE(NFIL,FMT=FMT2)
c------------------------------------------------------------ header end
      NNON0   = 0
      ND      = 0
      NALF    = 0
      CTAB(0) = ' '
      DTAB(0) = 9999D0

      DO 10 I=1,N
         DO 20 J=1,N
            IF( .NOT. CSMALL( A(I,J) ) ) THEN
               NNON0 = NNON0 + 1
               DO 30 ID=1,ND    
                  IF( CSAME(A(I,J),+DTAB(ID)) ) THEN
                     IW(J) = + ID
                     GOTO 40
                  END IF
                  IF( CSAME(A(I,J),-DTAB(ID)) ) THEN
                     IW(J) = - ID
                     GOTO 40
                  END IF
   30          CONTINUE     
c----------------------------------------------------------- new element
               ND = ND + 1
               IF( ND .GT. NDIFMAX ) THEN
                  WRITE(NFIL,'('' trouble in <RMATSTR> !!!!'',/,
     &              '' ND > array size NDIFMAX='',I3)') NDIFMAX
                  STOP
               END IF
               IW(J)    = ND
               DTAB(ND) = A(I,J)
               IF( ABS(DTAB(ND)-1.0D0)*TOL .LT. 1.0D0 ) THEN
                  CTAB(ND) = '1'
               ELSE IF( ABS(DTAB(ND)+1.0D0)*TOL .LT. 1.0D0 ) THEN
                  DTAB(ND) = +1.0D0
                  CTAB(ND) = '1'
                  IW(J) = - ND
               ELSE
                  NALF = NALF + 1     
                  IF( NALF .LE. NATOZ ) THEN
                     CTAB(ND) = CHAR( ICAUC + NALF - 1 )
                  ELSE 
                     ICND = ICALC + NALF-NATOZ - 1
                     IF( ICND .LT. ICILC ) THEN
                        CTAB(ND) = CHAR(ICND)
                     ELSE
                        CTAB(ND) = CHAR(ICND+1)
                     END IF
                  END IF
               END IF
   40          CONTINUE    
            ELSE
               IW(J) = 0
            END IF
   20    CONTINUE
c------------------------------------------------------------ write line
         WRITE(NFIL,FMT=FMT1)
     &        I,(VZ(ISIGN(1,IW(J))),CTAB(ABS(IW(J))),J=1,N),I

            
         DO ISL=1,NSL
            IF( I .EQ. ILSEP(ISL) ) WRITE(NFIL,FMT=FMT2)
         END DO
   10 CONTINUE

c------------------------------------------------------------------ foot
      
      WRITE(NFIL,FMT4) (I,I=1,N,2)
      WRITE(NFIL,FMT3) (I,I=2,N,2)

      WRITE(NFIL,9030) (ID,CTAB(ID),DTAB(ID),ID=1,ND)
      WRITE(NFIL,9040)  NNON0,N*N-NNON0
      
 9000 FORMAT(/,8X,A,/)
 9030 FORMAT(/,8X,'symbols used:',/,(8X,I3,3X,A1,2X, F20.12) )
 9040 FORMAT(/,8X,'elements <> 0:',I4,/,
     &         8X,'elements  = 0:',I4)
      RETURN
      END
