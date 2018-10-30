      SUBROUTINE ROTBRILL(DSYMLL,LMAX,NSYMAT,RSYMAT,ISYMINDEX,
     &                    LPOT,YR,WTYR,RIJ,IJEND)
c-----------------------------------------------------------------------
      IMPLICIT NONE
C     .. Parameters ..
      include 'inc.p'
c      INTEGER LMX,LPOTD
c      PARAMETER (LMX=4,LPOTD=6)
       
      INTEGER LMMAXD,NSYMAXD
      PARAMETER (LMMAXD=(LMAXD+1)**2,NSYMAXD=48)
      INTEGER LMPOTD,N
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER BRLSYM,LMAX,NSYMAT,LPOT,IJEND
C     ..
C     .. Array Arguments ..
      INTEGER ISYMINDEX(NSYMAXD)
      DOUBLE PRECISION C(0:LMAXD,-LMAXD:LMAXD,-LMAXD:LMAXD,NSYMD),
     +                 RSYMAT(64,3,3)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMD)
      DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*)
C     ..
C     .. Local Scalars ..
      LOGICAL TEST
      DOUBLE PRECISION ROTR,ROTR1,ROTR2,ROTR3
      INTEGER IJ,IL,K,L1,L2,LM1,LM2,M1,M2,M3,M4,LMMAX,KSYM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION YRROT(LMMAXD,IJD)
C     ..
C     .. External Subroutines ..
      EXTERNAL YMY,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
c      COMMON /RSPHERE/YR,WTYR,RIJ,X,W,IJEND
C     ..
C     .. Save statement ..
c      SAVE /RSPHERE/
C     ..


      lmmax = (lmax+1)**2
      if (NSYMAT.gt.NSYMAXD) then
         write(6,*)'SUB ROTBRILL ERROR STOP LOOK DIMENSIONS'
         STOP
      end if 
      DO 70 K = 1,NSYMAT

          KSYM=ISYMINDEX(K)

          do lm1=1,lmmax
             do lm2=1,lmmax
                dsymll(lm1,lm2,k) = 0.d0
             end do
          end do
c
c---> generate the rotated spherical harmonics on the angular mesh
c      (mesh points are generated in deck sphere)
c
          DO 80 IJ = 1,IJEND
c          write(6,*) "IJ",IJ
            ROTR1 = RSYMAT(KSYM,1,1)*RIJ(IJ,1) +
     &              RSYMAT(KSYM,1,2)*RIJ(IJ,2) +
     &              RSYMAT(KSYM,1,3)*RIJ(IJ,3)
            ROTR2 = RSYMAT(KSYM,2,1)*RIJ(IJ,1) + 
     &              RSYMAT(KSYM,2,2)*RIJ(IJ,2) +
     &              RSYMAT(KSYM,2,3)*RIJ(IJ,3)
            ROTR3 = RSYMAT(KSYM,3,1)*RIJ(IJ,1) + 
     &              RSYMAT(KSYM,3,2)*RIJ(IJ,2) +
     &              RSYMAT(KSYM,3,3)*RIJ(IJ,3)
c
c            write(6,*) rsymat(k,1:3,1:3)
c           write(6,*) rij(ij,1),rij(ij,2),rij(ij,3)
c             write(6,*) rotr1,rotr2,rotr3
              CALL YMY(ROTR1,ROTR2,ROTR3,ROTR,YRROT(1,IJ),LMAX)
   80     CONTINUE
c
c
c---> now integrate
c
c         write(6,*) "after YMY"
          DO 90 IL = 0,LMAX
c             write(6,*) "IL",IL
              DO 100 M1 = -IL,IL
                  DO 110 M2 = -IL,IL
                      C(IL,M1,M2,K) = (0.D0,0.D0)
                      LM1 = IL* (IL+1) + M1 + 1
                      LM2 = IL* (IL+1) + M2 + 1
                      DO 120 IJ = 1,IJEND
                          C(IL,M1,M2,K) = C(IL,M1,M2,K) +
     +                                    WTYR(IJ,LM1)*YRROT(LM2,IJ)
  120                 CONTINUE
                      IF (ABS(C(IL,M1,M2,K)).LE.1.0D-10) C(IL,M1,M2,
     +                    K) = 0.0D0
                      dsymll(lm1,lm2,k) = C(IL,M1,M2,K)
c                      write(6,*) lm1,lm2,k,dsymll(lm1,lm2,k)
  110             CONTINUE
  100         CONTINUE
   90     CONTINUE

 70    CONTINUE

      IF(TEST('flow    '))  
     +        write(6,*) '<<<<< ROTBRILL ended'

      END








