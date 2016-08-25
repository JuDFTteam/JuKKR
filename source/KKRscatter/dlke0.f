c 04.10.95 ***************************************************************
      SUBROUTINE DLKE0(GLLKE,GSPARSE,INGSP,NSPBLOCK,ALAT,NAEZ,CLS,EQINV,
     &                 NACLS,RR,EZOA,ATOM,BZKP,IE,KAOEZ,RCLS,GINP)
c
c
c                                      update for sparse matrix 14.6.2001
c   simple to turm of GLLKE array in case of big calculation
c   gsparse can be turned of from inc.p directly when not used
c   set   nauxspd to 1 to turn it off!          
c ************************************************************************
      implicit none
C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER ALM,CLM
      PARAMETER (ALM=LMAXSQ*NDIMGK,CLM=LMAXSQ*NACLSD)
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
      DOUBLE COMPLEX CONE,CONEM
      PARAMETER (CONE= (1.0D0,0.0D0), CONEM= (-1.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IE,NAEZ,NDIM,NSPBLOCK
C     ..
C     .. Array Arguments ..
      INTEGER, intent(in) :: 
     +     ATOM(NACLSD,NAEZD),
     +     CLS(NAEZ),
     +     EZOA(NACLSD,NAEZ),
     +     EQINV(NAEZ),
     +     KAOEZ(NAEZ+NEMBD),
     +     NACLS(*),INGSP(NLSPD,NLSPD)
c
      DOUBLE COMPLEX, intent(in)  :: 
     +     GINP(LMAXSQ*NACLSD,LMAXSQ,NCLSD)
      DOUBLE COMPLEX, intent(out) :: 
     +     GLLKE(ALM,ALM),
     &     GSPARSE(NSPBLOCK*LMAXSQ,NSPBLOCK*LMAXSQ,*)
c
      DOUBLE PRECISION 
     +     BZKP(6),
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IC,IM,INV,J,JM,M,N,JN,II,DI,JJ,DJ,LM1,LM2
      DOUBLE COMPLEX FAC
      LOGICAL OPT
c     ..
c     .. Local Arrays ..
      DOUBLE COMPLEX GLLKE1(ALM,LMAXSQ)
      DOUBLE PRECISION KP(6)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,DLKE1,OPT,ZAXPY
C     ..
C     .. Save statement ..
      SAVE
C     ..
c      write(6,*) '>>> DLKE0 : Fourier-transforms the ',
c     +           'GF of reference system'
c ------------------------------------------------------------------------

c      WRITE(6,*) "IN DLKE0"
      CALL CINIT(ALM*ALM,GLLKE(1,1))
      GLLKE=0d0
      GLLKE1=0d0
c      WRITE(6,*) "(BZKP(J),J=1,3)"
c      WRITE(6,*) (BZKP(J),J=1,3)
      NDIM = NSPBLOCK*LMAXSQ
      CALL CINIT(NDIM*NDIM*NAUXSPD,GSPARSE(1,1,1))
      DO 20 I=1,NAEZ

        FAC = CONE
c
c -->     cluster around i is invers symmetric 
c         to cluster around eqinv(i)
c
        IF ( I.NE.EQINV(I) ) FAC = CONEM

        KP(1) = BZKP(1)
        KP(2) = BZKP(2)
        KP(3) = BZKP(3)
        IF (OPT('COMPLEX ')) THEN
          KP(4) = BZKP(4)
          KP(5) = BZKP(5)
          KP(6) = BZKP(6)
        END IF

        IC  = CLS(KAOEZ(I))

c        write(6,*) 'fac :',fac
c        write(6,*) 'in dlke0 checing cluster ic :',ic,'for atom',i
c       WRITE(6,*) "ALAT",ALAT
c        WRITE(6,*) "NACLS",NACLS(I)
c       WRITE(6,*) "EZOA ",EZOA(1,I)
c       WRITE(6,*) "ATOM ",ATOM(1,I)
c       WRITE(6,*) "KP",(KP(J),J=1,6)
c       WRITE(6,*) "IE",IE            
c       WRITE(6,*) "IC",IC            
c       WRITE(6,*) "FAC",FAC          

c       WRITE(203,*) "IC",IC            
c       do LM1=1,LMAXSQ
c         do LM2=1,LMAXSQ*NACLSD
c           IF (LM2 .LE. LMAXSQ) THEN
c              write(203,"((2I5),(4e17.9))") LM1,LM2,
c    +             GINP(LM2,LM1,IC),GINP(LM1,LM2,IC)
c           ELSE
c              write(203,"((2I5),(4e17.9))") LM1,LM2,
c    +             GINP(LM2,LM1,IC)!,GINP(LM1,LM2,IC)
c           END IF
c         END DO
c       END DO
c        WRITE(6,*) "before DLKE1"

        CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,
     +       EZOA(1,I),ATOM(1,I),KP,IE,IC,FAC,GINP(1,1,IC),
     +       RCLS(1,1,IC))   ! Changed on 22.03.2000 ! forgoten I 3.11.2000
                               ! Correction for fourier trans.

c        WRITE(6,*) "after DLKE1"
c       do LM1=1,ALM
c         do LM2=1,LMAXSQ
c           write(205,"((2I5),(2e17.9))") LM1,LM2,GLLKE1(LM1,LM2)
c         end do
c       end do
c       write(205,*) " " 

c       WRITE(6,*) "after DLKE1"
c       do LM1=1,ALM
c         do LM2=1,LMAXSQ
c           write(205,"((2I5),(2e17.9))") LM1,LM2,GLLKE1(LM1,LM2)
c         end do
c       end do
c       write(205,*) " " 

        IF (OPT('sparse   ')) THEN
           
           II = (I-1)/NSPBLOCK + 1
           DI = I - (II-1)*NSPBLOCK
           
           DO J=1, NAEZ
              JJ = (J-1)/NSPBLOCK + 1
              DJ = J - (JJ-1)*NSPBLOCK
              IF (INGSP(JJ,II).NE.0) 
     &             CALL GLLCOPYA(GLLKE1,GSPARSE(1,1,INGSP(JJ,II)),
     &             NSPBLOCK,J,DJ,DI)
           END DO
           
        ELSE                    ! NOT sparse 
      
c           write(206,*) "I",I
           DO 140 M=1,LMAXSQ
              IM=(I-1)*LMAXSQ+M
              DO 150 JN=1,LMAXSQ*NAEZ
                 GLLKE(JN,IM) = GLLKE(JN,IM)+ GLLKE1(JN,M)
c                 write(206,"((3I5),(2e17.9))") I,M,JN,GLLKE(JN,IM)
 150          CONTINUE
 140       CONTINUE
           
        END IF                  ! end if sparse
        
 20   CONTINUE                  ! I=1,NAEZ
      
c     DO IM=1,ALM
c       DO JN=1,ALM
c         WRITE(207,"((2I5),(4e17.9))") IM,JN,GLLKE(JN,IM) 
c       END DO
c     END DO
c     ------------------------------------------------------------------------
      IF (OPT('symG(k) ')) THEN
         IF (OPT('sparse   ')) THEN
            write(6,*) ' No symmetrization for sparse !!'
            STOP
         END IF 
c     
c     -->   symmetrization
c     
         DO 90 I=1,NAEZ
            
            FAC = CONE
c
c -->     cluster around i is invers symmetric 
c         to cluster around inv = eqinv(i)
c
          IF (I.NE.EQINV(I)) FAC = CONEM
            
          KP(1) = -BZKP(1)
          KP(2) = -BZKP(2)
          KP(3) = -BZKP(3)
          IF (OPT('COMPLEX ')) THEN
            KP(4) = -BZKP(4)
            KP(5) = -BZKP(5)
            KP(6) = -BZKP(6)
          END IF
          
          IC  = CLS(KAOEZ(I))
          
          CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,
     +         EZOA(1,I),ATOM(1,I),KP,IE,IC,FAC,GINP(1,1,IC),
     +         RCLS(1,1,IC))

          DO 120 J=1,NAEZ
            DO 110 M=1,LMAXSQ
              IM=(I-1)*LMAXSQ+M
              DO 100 N=1,LMAXSQ
                JN=(J-1)*LMAXSQ+N
                GLLKE(IM,JN) = (GLLKE(IM,JN)+ GLLKE1(JN,M))/2.0D0
 100          CONTINUE
 110        CONTINUE
 120      CONTINUE
          
 90     CONTINUE                    ! I=1,NAEZ

      END IF                        ! (OPT('symG(k) '))
c ------------------------------------------------------------------------

      RETURN

      END                           ! DLKE0










