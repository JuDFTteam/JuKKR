      MODULE MOD_FORCXC
      CONTAINS
      SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NATOM,VPOT,DENSITY,CELL,ALAT,LMPOTD,IRMD,INS)
      USE TYPE_DENSITY
      USE TYPE_CELL
      USE MOD_SIMP3     
      IMPLICIT NONE
! c>>>>>BEWARE!!! RM commented away!!! -->Dipole Tensor is useless      
! c     SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
! c    +                  RM,NSHELL,DRDI,IRWS,NATREF)
! c-----------------------------------------------------------------------
! c     calculates the force on nucleus m
! c     from a given non spherical charge density at the nucleus site r
! c     with core correction(exchange contribution)
!  
! c-----------------------------------------------------------------------
! C     .. Parameters ..
!       include 'inc.p'
      INTEGER LMPOTD,IRMD,INS
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
! C     ..
! C     .. Scalar Arguments ..
      TYPE(CELL_TYPE) :: CELL(NATOM)
      TYPE(DENSITY_TYPE) :: DENSITY(NATOM)
      DOUBLE PRECISION ALAT
      INTEGER LMAX,NATOM,NSPIN,NSTART
! C     ..
! C     .. Array Arguments ..
      DOUBLE PRECISION FLM(-1:1,NATOM),FLMC(-1:1,NATOM) !,R(IRMD,*),
! !      +       RHOC(IRMD,*),
! c     $     RM(3,*),
      DOUBLE PRECISION VPOT(IRMD,LMPOTD,NSPIN,NATOM)
! c      INTEGER IRWS(*),NSHELL(*)
!        INTEGER IRWS(*)
! C     ..
! C     .. Local Scalars ..
      DOUBLE PRECISION DV,DVOL,FAC,RWS,TRP,VINT1,VOL
      INTEGER I,IPER, ISPIN,IATOM,IATOMP,IRWS1,J,LM,M
! C     ..
! C     .. Local Arrays ..
      DOUBLE PRECISION F(3,NATOM),FLMH(-1:1,NATOM), &
             FLMXC(-1:1,NATOM),V1(IRMD)
! C     ..
! C     .. External Subroutines ..
!       EXTERNAL SIMP3
! C     ..
! C     .. Save statement ..
!       SAVE PI
! C     ..
! c
! C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DSQRT
! C     ..
!       PI = 4.D0*ATAN(1.D0)
      FAC = DSQRT((4.0D0*PI)/3.0D0)
      TRP = 0.0D0
      IF (LMAX.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
! c
      WRITE (6,FMT=9200)
      WRITE (6,FMT=9100)
      WRITE (6,FMT=9200)
! c
!       IREP = 1
      DO 10 IATOM = 1,NATOM
! c
!          IPER = IATYP - NATREF
!          P(IPER) = 0.0D0
! !          WRITE (6,FMT=9400) IPER
! c
         IF (INS==1) THEN
           IRWS1 = CELL(IATOM)%NRCORE !IRWS(IATOM)
         ELSE
           IRWS1 = CELL(IATOM)%NRMAX !IRWS(IATOM)
         END IF
!          IRWS1 = CELL(IATOM)%NRMAX !IRWS(IATYP)
         RWS =   CELL(IATOM)%RMESH(IRWS1)
         VOL = 0.25*ALAT**3
! c
! c---> determine the right potential numbers
! c
 
         DO 20 M = -1,1
            LM = 2 + M + 1
! c
            DO 30 I = 1,IRWS1
               V1(I) = 0.0d0
   30       END DO
! c
            DO 40 ISPIN = 1,NSPIN
! c
!                ISPIN,IATOM = NSPIN* (IATYP-1) + ISPIN
! c
! c
               DV = (-3.0D0*VPOT(1,LM,ISPIN,IATOM)-10.0D0*VPOT(2,LM,ISPIN,IATOM)+ &
                   18.0D0*VPOT(3,LM,ISPIN,IATOM)-6.0D0*VPOT(4,LM,ISPIN,IATOM)+VPOT(5,LM,ISPIN,IATOM))/&
                    (12.0D0*CELL(IATOM)%DRMESHDI(2))
! c
               V1(2) = DENSITY(IATOM)%RHOC(2,ISPIN)* (2.0D0*VPOT(2,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(2)+DV)/&
                       (4.0D0*PI) + V1(2)
! c
               DO 50 I = 3,IRWS1 - 2
! c
                  DV = (VPOT(I-2,LM,ISPIN,IATOM)-VPOT(I+2,LM,ISPIN,IATOM)+&
                       8.0D0* (VPOT(I+1,LM,ISPIN,IATOM)-VPOT(I-1,LM,ISPIN,IATOM)))/&
                       (12.0D0*CELL(IATOM)%DRMESHDI(I))
! c
                  V1(I) = DENSITY(IATOM)%RHOC(I,ISPIN)* (2.0D0*VPOT(I,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(I)+&
                          DV)/ (4.0D0*PI) + V1(I)
   50          END DO
! c
               DV = (-VPOT(IRWS1-4,LM,ISPIN,IATOM)+6.0D0*VPOT(IRWS1-3,LM,ISPIN,IATOM)- &
                    18.0D0*VPOT(IRWS1-2,LM,ISPIN,IATOM)+10.0D0*VPOT(IRWS1-1,LM,ISPIN,IATOM)+ &
                   3.0D0*VPOT(IRWS1,LM,ISPIN,IATOM))/ (12.0D0*CELL(IATOM)%DRMESHDI(IRWS1-1))
               V1(IRWS1-1) = DENSITY(IATOM)%RHOC(IRWS1-1,ISPIN)* &
                             (2.0D0*VPOT(IRWS1-1,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(IRWS1-1)+ &
                             DV)/ (4.0D0*PI) + V1(IRWS1-1)
! c
               DV = (3.0D0*VPOT(IRWS1-4,LM,ISPIN,IATOM)-16.0D0*VPOT(IRWS1-3,LM,ISPIN,IATOM)+ &
                    36.0D0*VPOT(IRWS1-2,LM,ISPIN,IATOM)-48.0D0*VPOT(IRWS1-1,LM,ISPIN,IATOM)+ &
                    25.0D0*VPOT(IRWS1,LM,ISPIN,IATOM))/ (12.0D0*CELL(IATOM)%DRMESHDI(IRWS1))
! c
               V1(IRWS1) = DENSITY(IATOM)%RHOC(IRWS1,ISPIN)* &
                           (2.0D0*VPOT(IRWS1,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(IRWS1)+DV)/ &
                           (4.0D0*PI) + V1(IRWS1)
   40       END DO
! c
! c---> integrate with simpson subroutine
! c
            CALL SIMP3(V1,VINT1,1,IRWS1,CELL(IATOM)%DRMESHDI(1))
! c
            FLMH(M,IATOM) = FLM(M,IATOM) - FLMC(M,IATOM)
            FLMXC(M,IATOM) = -FAC*VINT1 - FLMC(M,IATOM)
            FLM(M,IATOM) = FLM(M,IATOM) + FLMXC(M,IATOM)
! c
 
   20    END DO
! c
         WRITE (6,FMT=9600) FLMH(1,IATOM),FLMC(1,IATOM),FLMXC(1,IATOM), &
           FLM(1,IATOM)
         WRITE (6,FMT=9601) FLMH(-1,IATOM),FLMC(-1,IATOM), &
           FLMXC(-1,IATOM),FLM(-1,IATOM)
         WRITE (6,FMT=9602) FLMH(0,IATOM),FLMC(0,IATOM),FLMXC(0,IATOM), &
           FLM(0,IATOM)
! c
         F(1,IATOM) = FLM(1,IATOM)
         F(2,IATOM) = FLM(-1,IATOM)
         F(3,IATOM) = FLM(0,IATOM)
! c
! c
! c         DO 60 J = 1,3
! c            P(IPER) = P(IPER) + RM(J,IREP)*NSHELL(IPER)*F(J,IATYP)*ALAT
! c   60    END DO
! c         TRP = TRP + P(IPER)
! c
! c         IREP = IREP + NSHELL(IPER)
! c
! c        write (6,*) '-->Tensor is useless'
! c        WRITE (6,FMT=9700) P(IPER)
! c
   10 END DO
! c
! c     DVOL = TRP/ (3.0D0*VOL)
! c
      WRITE (6,FMT=9200)
! c     WRITE (6,FMT=9101)
! c     WRITE (6,FMT=9200)
! c     WRITE (6,FMT=9800) DVOL
! c     WRITE (6,FMT=9200)
      WRITE (6,FMT=9102)
      WRITE (6,FMT=9200)
! c
 9000 FORMAT (13x,'error stop in subroutine force :', &
             ' the charge density has to contain non spherical',&
             ' contributions up to l=1 at least ')
 9101 FORMAT (1x,33 ('-'),' volume change ',33 ('-'),/,34x, &
             ' in units Ry/(a(Bohr)**3 ')
 9102 FORMAT (1x,81 ('-'))
 9100 FORMAT (1x,33 ('-'),' force on the nucleus ',33 ('-'),/,34x,&
             ' in units Ry/(a(Bohr) ')
 9200 FORMAT (1x,'>')
 9400 FORMAT (3x,i5,'th shell')
 9600 FORMAT (7x,'fhx=',e12.6,2x,'fcx=',e12.6,2x,'fxcx=',e12.6,2x,'fx=',&
             e12.6,' Ry/(a(Bohr))')
 9601 FORMAT (7x,'fhy=',e12.6,2x,'fcy=',e12.6,2x,'fxcy=',e12.6,2x,'fy=', &
             e12.6,' Ry/(a(Bohr))')
 9602 FORMAT (7x,'fhz=',e12.6,2x,'fcz=',e12.6,2x,'fxcz=',e12.6,2x,'fz=',&
             e12.6,' Ry/(a(Bohr))')
 9700 FORMAT (10x,'contribution to the trace of the dipol force tensor:'&
             ,3x,e12.6,' Ry')
 9800 FORMAT (7x,' volume change dvol/vol=',2x,e12.6,' Ry/(a(Bohr))**3',&
             /,7x,'( notice: has to be divided',&
             ' by the bulk modulus of the host)')
 
      END SUBROUTINE
      END MODULE MOD_FORCXC
