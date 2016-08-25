c ************************************************************************
      SUBROUTINE SP2(T,IN,NUMT,N,DIM,XX,INX,NUMX,NAUX,IO,IOPT)
c ************************************************************************
c
c     XX = T**(-1)
c     T      nonzero blocks of matrix T to invert
c     IN     index array IN(I,J) = position of block (I,J)
c            in linear array T(*)
c     NUMT   number of nonzero blocks in T
c     N      block dimension of matrix T
c     DIM    dimension of blocks in T
c     XX     nonzero blocks of T**(-1)
c     INX    index array INX(I,J) for XX (see IN)
c     NUMX   number of nonzero blocks in XX
c     NAUX   block dimension of linear array T (.ge. NUMT)
c     IO     option for output of interim results for test purposes 
c            (io.ge.1)
c     IOPT   returns number of required fields in array T for
c            left-right decomposition (dimension NAUX)
c
c
cc ************************************************************************
      include 'inc.p'
      PARAMETER (LMMAXD = (LMAXD+1)**2,
     +           NDIM   = LMMAXD*1) ! NSPBLOCK = 1
      DOUBLE COMPLEX 
     +     T(NDIM,NDIM,*),
     +     XX(NDIM,NDIM,*)
      INTEGER IO,IOPT,DIM,N,NAUX,NUMT,NUMX
      INTEGER IN(NLSPD,*),INX(NLSPD,*)
C
C     .. LOCALS
      DOUBLE COMPLEX 
     +     DD(NDIM,NDIM,NLSPD),
     +     RR(NDIM,NDIM,NLSPD),
     +     U(NDIM,NDIM),
     +     UNIT(NDIM,NDIM),
     +     U1(NDIM,NDIM),
     +     X1(NDIM,NDIM,NLSPD)
      INTEGER INR(NLSPD),INXX(NLSPD),IPVT(NDIM)

      INTEGER I,IEND,IND,INFO,J,K,KN,LM,NUMR
      DOUBLE COMPLEX FAC
      
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER(CONE=(1.0D0,0.0D0),CZERO=(0.0D0,0.0D0))

      EXTERNAL CINIT
c ------------------------------------------------------------------------
      IF (DIM.NE.NDIM) STOP 'DIM.NE.NDIM IN ROUTINE SP2'
C
c
c ---> initialize unitary matrix UNIT
c
      CALL CINIT(DIM*DIM,UNIT)
      DO 5 LM = 1,DIM
        UNIT(LM,LM) = CONE
 5    END DO
c
      IND = NUMT + 1
      IF (IO.GT.0) write(6,*) 'NUMT',IND-1
c
c ---> LR decomposition of T = L * R,
c      L and R are stored in T
c
      DO 10 J = 1,N-1
        IF (IN(J,J).EQ.0) THEN
          STOP '1 - SP'
        END IF
c
c --->  U = T(1,1,IN(J,J))**(-1)
c
c        CALL MINV(T(1,1,IN(J,J)),U,U1,DIM)
        CALL ZCOPY(DIM*DIM,UNIT,1,U,1)
        CALL ZCOPY(DIM*DIM,T(1,1,IN(J,J)),1,U1,1)
        CALL ZGETRF(DIM,DIM,U1,DIM,IPVT,INFO)
        CALL ZGETRS('N',DIM,DIM,U1,DIM,IPVT,U,DIM,INFO)
        DO 20 I = J+1,N
          IF (IN(I,J).NE.0) THEN
c
c --->      T(1,1,IN(I,J)) = T(1,1,IN(I,J)) * U
c
c            CALL MMUL(T(1,1,IN(I,J)),U,U1,DIM)
            CALL ZGEMM('N','N',DIM,DIM,DIM,CONE,
     +           T(1,1,IN(I,J)),DIM,U,DIM,CZERO,U1,DIM)
            CALL ZCOPY(DIM*DIM,U1,1,T(1,1,IN(I,J)),1)
            DO 30 K = J+1,N
              IF (IN(J,K).NE.0) THEN
                FAC = CONE
                IF (IN(I,K).EQ.0) THEN
                  IN(I,K) = IND
                  IF (IND.GT.NAUX) THEN
                       WRITE(6,*) IND,NAUX 
                       STOP 'INCREASE NAUX IN ROUTINE SP1 (LR)'
                  END IF
                  IND = IND+1
                  FAC = CZERO
                END IF
                CALL ZGEMM('N','N',DIM,DIM,DIM,
     +               -CONE,T(1,1,IN(I,J)),DIM,T(1,1,IN(J,K)),DIM,
     +               FAC,T(1,1,IN(I,K)),DIM)
              END IF
 30         END DO
          END IF
 20     END DO
 10   END DO

      IOPT = IND-1
      IF (IO.GT.0) write(6,*) 'NUM_LR',IND-1

      IF (IO.GT.1) THEN
        write(6,*) 'L,R :'
        DO 33 I=1,NLSPD
          DO 34 LMI=1,NDIM
            write(6,FMT=9000) 
     +           ((DREAL(T(LMI,LMJ,IN(I,J))),LMJ=1,NDIM),J=1,NLSPD)
 34       CONTINUE
 33     CONTINUE
        DO 35 I=1,NLSPD
          write(6,FMT=9010) (IN(I,J),J=1,NLSPD)
 35     CONTINUE
      END IF

      DO 100 K = 1,N
c
c --->  DD(K) = T(K,K)**(-1)
c
c        CALL MINV(T(1,1,IN(K,K)),DD(1,1,K),U1,DIM)
        CALL ZCOPY(DIM*DIM,UNIT,1,DD(1,1,K),1)
        CALL ZCOPY(DIM*DIM,T(1,1,IN(K,K)),1,U1,1)
        CALL ZGETRF(DIM,DIM,U1,DIM,IPVT,INFO)
        CALL ZGETRS('N',DIM,DIM,U1,DIM,IPVT,DD(1,1,K),DIM,INFO)
 100  END DO

      DO 40 K = 1,N
c
c ---> forward eleminitation of  1 = L * RR (inverse of L)
c
C     initialize arrays RR and INR
        CALL CINIT(DIM*DIM*NLSPD,RR)
        DO 31 I=1,NLSPD
          INR(I) = 0
 31     END DO
        IND = 1
c
        INR(K) = IND
        IND = IND+1
        CALL ZCOPY(DIM*DIM,UNIT,1,RR(1,1,K),1)
        DO 50 I = K+1,N
          DO 60 J = K,I-1
            IF (IN(I,J)*INR(J).NE.0) THEN
              FAC = CONE
              IF (INR(I).EQ.0) THEN
                INR(I) = IND
                IND = IND+1
                FAC = CZERO
              END IF
              CALL ZGEMM('N','N',DIM,DIM,DIM,
     +             -CONE,T(1,1,IN(I,J)),DIM,RR(1,1,J),DIM,
     +             FAC,RR(1,1,I),DIM)
            END IF
 60       END DO
 50     END DO

c        IF (IO.GT.0) write(6,*) 'NUM_RR(',K,')',IND-1

        IF (IO.GT.1) THEN
          write(6,*) 'RR :'
          DO 90 I=1,NLSPD
            DO 91 LMI=1,NDIM
              write(6,FMT=9000) (DREAL(RR(LMI,LMJ,I)),LMJ=1,DIM)
 91         CONTINUE
 90       CONTINUE
          DO 94 I=1,NLSPD
            write(6,FMT=9010) INR(I)
 94       CONTINUE
        END IF
        
c     
c --->   backward elimination
c
        CALL CINIT(DIM*DIM*NLSPD,X1)
        DO 32 I=1,NLSPD
          INXX(I) = 0
 32     END DO
        IND = 1
c
c --->  X1(N) = T(N,N)**(-1) * RR(N)
c
        IF (INR(N).NE.0) THEN
          INXX(N) = IND
          IND = IND + 1 
          CALL ZGEMM('N','N',DIM,DIM,DIM,
     +         CONE,DD(1,1,N),DIM,RR(1,1,N),DIM,
     +         CZERO,X1(1,1,N),DIM)
        END IF
c
c --->  lower bound for backward elemination due to index array INX()
c
        IEND = N
        DO 150 KN = N,1,-1
          IF (INX(KN,K).NE.0) IEND = KN
 150    END DO

        DO 120 I = N-1,IEND,-1
c
c --->    X1(I) = RR(I)
c
          IF (INR(I).NE.0) THEN
            INXX(I) = IND
            IND = IND+1
            CALL ZCOPY(DIM*DIM,RR(1,1,I),1,X1(1,1,I),1)
          END IF
          DO 130 J = I+1,N
            IF (IN(I,J)*INXX(J).NE.0) THEN
              FAC = CONE
              IF (INXX(I).EQ.0) THEN
                INXX(I) = IND
                IND = IND+1
                FAC = CZERO
              END IF
              CALL ZGEMM('N','N',DIM,DIM,DIM,
     +             -CONE,T(1,1,IN(I,J)),DIM,X1(1,1,J),DIM,
     +             FAC,X1(1,1,I),DIM)
            END IF
 130      END DO                    ! J = I+1,N
          IF (INXX(I).NE.0) THEN
c
c --->      X1(1,1,I) = DD(1,1,I)*X1(1,1,I)
c
c            CALL MMUL(DD(1,1,I),X1(1,1,I),U(1,1),DIM)
            CALL ZGEMM('N','N',DIM,DIM,DIM,CONE,DD(1,1,I),DIM,
     +           X1(1,1,I),DIM,CZERO,U(1,1),DIM)
            CALL ZCOPY(DIM*DIM,U(1,1),1,X1(1,1,I),1)
          END IF
 120    END DO                      ! I = N-1,K,-1

        DO 140 KN = 1,N
c
c --->    copy blocks labeled in INX() to array XX
c
          IF (INX(KN,K)*INXX(KN).NE.0) 
     +         CALL ZCOPY(DIM*DIM,X1(1,1,KN),1,XX(1,1,INX(KN,K)),1)
 140    END DO                      ! KN = 1,N

 40   END DO                        ! K = 1,N


      IF (IO.GT.1) THEN
        write(6,*) 'XX :'
        DO 190 I=1,NLSPD
          DO 191 LMI=1,NDIM
            write(6,FMT=9000) 
     +        ((DREAL(XX(LMI,LMJ,INX(I,J))),LMJ=1,DIM),J=1,NLSPD)
 191      CONTINUE
 190    CONTINUE
        DO 194 I=1,NLSPD
          write(6,FMT=9010) (INX(I,J),J=1,NLSPD)
 194    CONTINUE
      END IF

      RETURN

 9000 FORMAT( 1p,  12d10.2)
 9010 FORMAT(   12i10)

      END
