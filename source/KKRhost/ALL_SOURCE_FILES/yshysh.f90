!-------------------------------------------------------------------------------
subroutine YSHYSH(X,Y,Z,R,YREALY)
   INTEGER LMAX
   PARAMETER (LMAX=LMAXD)
   INTEGER LMAX2
   PARAMETER (LMAX2=2*LMAX)
   INTEGER LMAX2P,LMXP
   PARAMETER (LMAX2P=LMAX2+1,LMXP=LMAX2P* (LMAX2P+1)/2)
   !     ..
   !     .. Scalar Arguments ..
   DOUBLE PRECISION R,X,Y,Z
   !     ..
   !     .. Array Arguments ..
   DOUBLE PRECISION YREALY(*)
   !     ..
   !     .. Local Scalars ..
   DOUBLE PRECISION A,ARG,B,C,COSPHI,COSTHE,EAT,P,PSQ,RSSQ,SAVE,SINPHI,SINTHE,TAVE,TENT,W,XA,YA
   INTEGER I,IA,IB,IC,ISTOPZ,ISUZY,J,JC,K,KOV2,L,LAVE,LSUZY,LTWOQP,M,MAVE,MP1,MSUZY,N
   !     ..
   !     .. Local Arrays ..
   DOUBLE PRECISION COSMPH(LMAX2P),FACTOR(50),PLM(LMXP),SINMPH(LMAX2P)
   !     ..
   !     .. Intrinsic Functions ..
   INTRINSIC DABS,DSIGN,DSQRT
   !     ..
   FACTOR(1) = 1.0D00
   do I = 2,50
      XA = I - 1
      FACTOR(I) = XA*FACTOR(I-1)
   enddo
   PSQ = X*X + Y*Y
   RSSQ = PSQ + Z*Z
   IF (RSSQ-1.0D-10) 20,20,30
   20 SINTHE = 0.0D00
      COSTHE = 1.0D00
      COSPHI = 1.0D00
      SINPHI = 0.0D00
      R = 0.0D00
   GO TO 60

   30 IF (PSQ-1.0D-10) 40,40,50
      40 R = DSQRT(RSSQ)
      SINTHE = 0.0D00
      COSTHE = DSIGN(1.0D00,Z)
      SINPHI = 0.0D00
      COSPHI = 1.0D00
      GO TO 60

   50 R = DSQRT(RSSQ)
      P = DSQRT(PSQ)
      SINTHE = P/R
      COSTHE = Z/R
      SINPHI = Y/P
      COSPHI = X/P
      !
   60 XA = DABS(COSTHE)
      YA = DABS(SINTHE)
      !      write(6,*) costhe,sinthe
      !
   IF (XA-1.0D-08) 70,70,150
   70 L = 0
      J = 0
      TAVE = 1.0D00
      LSUZY = 1
   80 M = 0
      MSUZY = 1
   90 J = J + 1
      ISUZY = LSUZY*MSUZY
   IF (ISUZY.GT.0) GO TO 100
   PLM(J) = 0.0D00
   GO TO 110

   100 K = L + M
      KOV2 = K/2
      IA = K + 1
      IB = KOV2 + 1
      JC = KOV2 - M
      IC = JC + 1
      PLM(J) = (((-1)**JC)*FACTOR(IA))/ (TAVE*FACTOR(IB)*FACTOR(IC))
   110 IF (M-L) 120,130,130
   120 M = M + 1
      MSUZY = -MSUZY
   GO TO 90

   130 IF (L-LMAX2) 140,370,370
   140 L = L + 1
      LSUZY = -LSUZY
      TAVE = 2.0D00*TAVE
   GO TO 80
      !
   150 IF (XA-0.99999999D00) 250,160,160
   160 PLM(1) = 1.0D00
      PLM(2) = COSTHE
      L = 2
      J = 2
   170 J = J + L
      A = L
      LTWOQP = L + L
      B = LTWOQP - 1
      C = L - 1
      K = J - L
      M = J - LTWOQP + 1
      PLM(J) = (B*COSTHE*PLM(K)-C*PLM(M))/A
   IF (L-LMAX2) 180,190,190
   180 L = L + 1
      GO TO 170
      !
      190 L = 1
      LAVE = 1
      200 M = 1
      LAVE = LAVE + L
      210 J = LAVE + M
      PLM(J) = 0.0D00
      IF (M-L) 220,230,230
      220 M = M + 1
      GO TO 210

      230 IF (L-LMAX2) 240,370,370
      240 L = L + 1
      GO TO 200
      C
      250 TENT = (2.0D00*COSTHE)/YA
      PLM(1) = 1.0D00
      PLM(2) = COSTHE
      PLM(3) = YA
      PLM(5) = 3.0D00*YA*COSTHE
      L = 2
      J = 2
      260 J = J + L
      A = L
      LTWOQP = L + L
      B = LTWOQP - 1
      C = L - 1
      K = J - L
      M = J - LTWOQP + 1
      PLM(J) = (B*COSTHE*PLM(K)-C*PLM(M))/A
      IF (L-LMAX2) 270,280,280
      270 L = L + 1
      GO TO 260

      280 L = 3
      J = 5
      290 J = J + L
      A = L - 1
      LTWOQP = L + L
      B = LTWOQP - 1
      C = L
      K = J - L
      M = J - LTWOQP + 1
      PLM(J) = (B*COSTHE*PLM(K)-C*PLM(M))/A
      IF (L-LMAX2) 300,310,310
      300 L = L + 1
      GO TO 290

      310 L = 2
      LAVE = 3
      320 LAVE = LAVE + L
      M = 1
      330 J = LAVE + M
      K = J - 1
      N = K - 1
      EAT = M
      A = TENT*EAT
      B = (M+L)* (L-M+1)
      PLM(J) = A*PLM(K) - B*PLM(N)
      IF (M+1-L) 340,350,350
      340 M = M + 1
      GO TO 330

      350 IF (L-LMAX2) 360,370,370
      360 L = L + 1
      GO TO 320
      C
      370 SINMPH(1) = 0.0D00
      COSMPH(1) = 1.0D00
      ISTOPZ = LMAX2 + 1
      DO 380 I = 2,ISTOPZ
         J = I - 1
         SINMPH(I) = SINPHI*COSMPH(J) + COSPHI*SINMPH(J)
         COSMPH(I) = COSPHI*COSMPH(J) - SINPHI*SINMPH(J)
         380 CONTINUE
         C
         L = 0
         390 M = 0
         LAVE = L* (L+1) + 1
         MAVE = ((L* (L+1))/2) + 1
         SAVE = 2*L + 1
         400 IF (M.NE.0) GO TO 410
         ARG = SAVE/12.5663706144D00
         W = DSQRT(ARG)
         YREALY(LAVE) = W*PLM(MAVE)
         GO TO 420

         410 IA = L - M + 1
         IB = L + M + 1
         ARG = (SAVE*FACTOR(IA))/ (6.28318530718D00*FACTOR(IB))
         MP1 = M + 1
         W = DSQRT(ARG)
         I = LAVE + M
         J = MAVE + M
         YREALY(I) = W*PLM(J)*COSMPH(MP1)
         I = LAVE - M
         YREALY(I) = W*PLM(J)*SINMPH(MP1)
         420 IF (M.GE.L) GO TO 430
         M = M + 1
         GO TO 400

         430 IF (L.GE.LMAX2) GO TO 440
         L = L + 1
         GO TO 390
         C
         440 RETURN
         C
         C
         C
      end subroutine YSHYSH
