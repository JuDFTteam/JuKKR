      SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
     +                  DRDI,IRWS,NATREF)
      use mod_types, only: t_inc
      IMPLICIT NONE
c>>>>>BEWARE!!! RM commented away!!! -->Dipole Tensor is useless      
c     SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
c    +                  RM,NSHELL,DRDI,IRWS,NATREF)
c-----------------------------------------------------------------------
c     calculates the force on nucleus m
c     from a given non spherical charge density at the nucleus site r
c     with core correction(exchange contribution)
 
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER LMAX,NATREF,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMD,*),
     +       RHOC(IRMD,*),
c     $     RM(3,*),
     $     V(IRMD,LMPOTD,*)
c      INTEGER IRWS(*),NSHELL(*)
       INTEGER IRWS(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DV,DVOL,FAC,PI,RWS,TRP,VINT1,VOL
      INTEGER I,IATYP,IPER,IPOT,IREP,IRWS1,ISPIN,J,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION F(3,NATYPD),FLMH(-1:1,NATYPD),
     +       FLMXC(-1:1,NATYPD),P(NATYPD),V1(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3
C     ..
C     .. Save statement ..
      SAVE PI
C     ..
c
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DSQRT
C     ..
      PI = 4.D0*ATAN(1.D0)
      FAC = DSQRT((4.0D0*PI)/3.0D0)
      TRP = 0.0D0
      IF (LMAX.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
c
      if(t_inc%i_write>0) WRITE (1337,FMT=9200)
      if(t_inc%i_write>0) WRITE (1337,FMT=9100)
      if(t_inc%i_write>0) WRITE (1337,FMT=9200)
c
      IREP = 1
      DO 10 IATYP = NSTART,NEND
c
         IPER = IATYP - NATREF
         P(IPER) = 0.0D0
         if(t_inc%i_write>0) WRITE (1337,FMT=9400) IPER
c
         IRWS1 = IRWS(IATYP)
         RWS = R(IRWS1,IATYP)
         VOL = 0.25*ALAT**3
c
c---> determine the right potential numbers
c
 
         DO 20 M = -1,1
            LM = 2 + M + 1
c
            DO 30 I = 1,IRWS1
               V1(I) = 0.0d0
   30       CONTINUE
c
            DO 40 ISPIN = 1,NSPIN
c
               IPOT = NSPIN* (IATYP-1) + ISPIN
c
c
               DV = (-3.0D0*V(1,LM,IPOT)-10.0D0*V(2,LM,IPOT)+
     +             18.0D0*V(3,LM,IPOT)-6.0D0*V(4,LM,IPOT)+V(5,LM,IPOT))/
     +              (12.0D0*DRDI(2,IATYP))
c
               V1(2) = RHOC(2,IPOT)* (2.0D0*V(2,LM,IPOT)/R(2,IATYP)+DV)/
     +                 (4.0D0*PI) + V1(2)
c
               DO 50 I = 3,IRWS1 - 2
c
                  DV = (V(I-2,LM,IPOT)-V(I+2,LM,IPOT)+
     +                 8.0D0* (V(I+1,LM,IPOT)-V(I-1,LM,IPOT)))/
     +                 (12.0D0*DRDI(I,IATYP))
c
                  V1(I) = RHOC(I,IPOT)* (2.0D0*V(I,LM,IPOT)/R(I,IATYP)+
     +                    DV)/ (4.0D0*PI) + V1(I)
   50          CONTINUE
c
               DV = (-V(IRWS1-4,LM,IPOT)+6.0D0*V(IRWS1-3,LM,IPOT)-
     +              18.0D0*V(IRWS1-2,LM,IPOT)+10.0D0*V(IRWS1-1,LM,IPOT)+
     +             3.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1-1,IATYP))
               V1(IRWS1-1) = RHOC(IRWS1-1,IPOT)*
     +                       (2.0D0*V(IRWS1-1,LM,IPOT)/R(IRWS1-1,IATYP)+
     +                       DV)/ (4.0D0*PI) + V1(IRWS1-1)
c
               DV = (3.0D0*V(IRWS1-4,LM,IPOT)-16.0D0*V(IRWS1-3,LM,IPOT)+
     +              36.0D0*V(IRWS1-2,LM,IPOT)-48.0D0*V(IRWS1-1,LM,IPOT)+
     +              25.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1,IATYP))
c
               V1(IRWS1) = RHOC(IRWS1,IPOT)*
     +                     (2.0D0*V(IRWS1,LM,IPOT)/R(IRWS1,IATYP)+DV)/
     +                     (4.0D0*PI) + V1(IRWS1)
   40       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,DRDI(1,IATYP))
c
            FLMH(M,IATYP) = FLM(M,IATYP) - FLMC(M,IATYP)
            FLMXC(M,IATYP) = -FAC*VINT1 - FLMC(M,IATYP)
            FLM(M,IATYP) = FLM(M,IATYP) + FLMXC(M,IATYP)
c
 
   20    CONTINUE
c
         if(t_inc%i_write>0) then
         WRITE (1337,FMT=9600) FLMH(1,IATYP),FLMC(1,IATYP),
     +     FLMXC(1,IATYP),FLM(1,IATYP)
         WRITE (1337,FMT=9601) FLMH(-1,IATYP),FLMC(-1,IATYP),
     +     FLMXC(-1,IATYP),FLM(-1,IATYP)
         WRITE (1337,FMT=9602) FLMH(0,IATYP),FLMC(0,IATYP),
     +     FLMXC(0,IATYP),FLM(0,IATYP)
         endif
c
         F(1,IATYP) = FLM(1,IATYP)
         F(2,IATYP) = FLM(-1,IATYP)
         F(3,IATYP) = FLM(0,IATYP)
c
c
c         DO 60 J = 1,3
c            P(IPER) = P(IPER) + RM(J,IREP)*NSHELL(IPER)*F(J,IATYP)*ALAT
c   60    CONTINUE
c         TRP = TRP + P(IPER)
c
c         IREP = IREP + NSHELL(IPER)
c
c        write (6,*) '-->Tensor is useless'
c        WRITE (6,FMT=9700) P(IPER)
c
   10 CONTINUE
c
c     DVOL = TRP/ (3.0D0*VOL)
c
      if(t_inc%i_write>0) then
        WRITE (1337,FMT=9200)
c       WRITE (6,FMT=9101)
c       WRITE (6,FMT=9200)
c       WRITE (6,FMT=9800) DVOL
c       WRITE (6,FMT=9200)
        WRITE (1337,FMT=9102)
        WRITE (1337,FMT=9200)
      endif
c
 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least ')
 9101 FORMAT (1x,33 ('-'),' volume change ',33 ('-'),/,34x,
     +       ' in units Ry/(a(Bohr)**3 ')
 9102 FORMAT (1x,81 ('-'))
 9100 FORMAT (1x,33 ('-'),' force on the nucleus ',33 ('-'),/,34x,
     +       ' in units Ry/(a(Bohr) ')
 9200 FORMAT (1x,'>')
 9400 FORMAT (3x,i5,'th shell')
 9600 FORMAT (7x,'fhx=',e12.6,2x,'fcx=',e12.6,2x,'fxcx=',e12.6,2x,'fx=',
     +       e12.6,' Ry/(a(Bohr))')
 9601 FORMAT (7x,'fhy=',e12.6,2x,'fcy=',e12.6,2x,'fxcy=',e12.6,2x,'fy=',
     +       e12.6,' Ry/(a(Bohr))')
 9602 FORMAT (7x,'fhz=',e12.6,2x,'fcz=',e12.6,2x,'fxcz=',e12.6,2x,'fz=',
     +       e12.6,' Ry/(a(Bohr))')
 9700 FORMAT (10x,'contribution to the trace of the dipol force tensor:'
     +       ,3x,e12.6,' Ry')
 9800 FORMAT (7x,' volume change dvol/vol=',2x,e12.6,' Ry/(a(Bohr))**3',
     +       /,7x,'( notice: has to be divided',
     +       ' by the bulk modulus of the host)')
 
      END
