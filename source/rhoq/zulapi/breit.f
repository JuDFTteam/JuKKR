c 20.09.95 ***************************************************************
      SUBROUTINE BREIT(BREITI,NONSRA,IR,IS,NRE,L,FOURPI,GAUSS,R2RHO1,
     +                 R2RHO2,RSQ,WGT,R,DRDI,G,F)
c ************************************************************************
c
c this function calculates the breit-integral g(r)*f(r)/r/r for spin up
c in relativistic case and the difference of fermi-contact-terms in non-
c relativistic case
c al is a weight depending on l which comes from integral of spherical-
c harmonics (a*a+b*b)*ylm*ylm l is real l ( l=0,1,2,3....)
c
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
c      INTEGER IRMD
c      PARAMETER (IRMD=424)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION BREITI,FOURPI,GAUSS,R2RHO1,R2RHO2,RSQ
      INTEGER IR,IS,L,NONSRA,NRE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(*),F(IRMD,*),G(IRMD,*),R(*),WGT(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AL,C,RL,RR,SUM1,SUM2,ZERO
      INTEGER IDD,IH,IRP1,IUU,JR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD,REAL
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA ZERO,C/0.0D0,274.072D0/
C     ..
c
      IF (IS.NE.2) THEN
        BREITI = ZERO
c changed from NONSRA.EQ.0   p.zahn, 20.4.98
      ELSE IF (NONSRA.NE.1) THEN
c
        RL = REAL(L)
        AL = (3.0D0*RL*RL-RL-1.0D0)/ (4.0D0*RL*RL-1)
        IDD = 1
        IUU = IS
        SUM1 = ZERO
        SUM2 = ZERO
c
c simpson rule for integration of g*f/r/r
c
        IRP1 = IR + 1
        DO 10 JR = IRP1,NRE
          RR = R(JR)*R(JR)
          RR = DRDI(JR)/RR
          IH = 1 + MOD(IRP1+1+JR,2)
          SUM1 = SUM1 + REAL(IH)*RR*G(JR,IDD)*F(JR,IDD)
          SUM2 = SUM2 + REAL(IH)*RR*G(JR,IUU)*F(JR,IUU)
c           if(ip.eq.2.and.l.eq.0.and.ir.eq.8)write(6,60)jr,rr,sum1,sum2
   10   CONTINUE
c
c  60 format(1x,i5,1p3d14.8)
c
        RR = DRDI(IR)/R(IR)/R(IR)
        SUM1 = WGT(IDD)* (SUM1+SUM1+RR*G(IR,IDD)*F(IR,IDD))
        SUM2 = WGT(IUU)* (SUM2+SUM2+RR*G(IR,IUU)*F(IR,IUU))
        BREITI = -2.0D0*C*GAUSS/FOURPI*AL* (SUM2-SUM1)/3.0D0

      ELSE
        BREITI = GAUSS* (R2RHO2-R2RHO1)/RSQ
      END IF

      END
