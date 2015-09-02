C             New call for main2
C             call EPOTINB_NEW(EPOTIN,NSPIN,RHO2NS,VISP,R(:,I1),DRDI(:,I1), &
C             IRMIN(I1),IRWS(I1),LPOT,VINS,IRCUT(:,I1),IPAN(I1),ZAT(I1), &
C             irmd, irnsd, ipand)

c 13.10.95 ***************************************************************
      SUBROUTINE EPOTINB_NEW(EPOTIN,NSPIN,RHO2NS,VM2Z,R,DRDI,
     &                   IRMIN,IRWS,LPOT,VINS,IRCUT,IPAN,Z,
C                        new input parameters after inc.p removal
     &                   irmd, irnsd, ipand)
c ************************************************************************
c
c     attention : energy zero ---> electro static zero
c
c                 since input potential and single particle energies
c                 are using muffin tin zero as zero the energy shift
c                 is cancelled in the kinetic energy contribution !
c
c
c     calculate the energy of the input potential
c     the energy for the representive atom i is given by
c
c                               rws
c       epotin(i) = - sqrt(4 pi) {  dr' vm2z(r',i)*rho2ns(r',1,i)
c                                0
c
c     in case of non spherical input potential one has to add
c
c                 rirt
c            {  -  {  dr' vins(r',lm,i)rho2ns(r',lm,i)   }
c                 rmin
c                                        (summed over lm)
c
c     remember : the non spherical part of the input potential is
c                different from zero only between r(irmin) and r(irt)
c
c             (see notes by b.drittler)
c
c     attention: vm2z is the spherically averaged input potential ,
c                vins contains the non spherical contribution of the
c                potential and rho2ns(...,1) is the  real charge density
c                times r**2. vins and rho2ns are expanded into spherical
c                harmonics . (see deck rholm or rhons)
c
c     remember :  in case of shape corrections  the contribution of
c                 the nuclear potential - 2*Z/r has to be explicitly
c                 taken into account between muffin tin sphere and
c                 circum scribed sphere .
c                 only within the muffin tin sphere this term is
c                 analytically cancelled wtih the contribution of
c                 the coulomb potential - see deck ecoulom
c
c
c                 modified for non spherical potential and shape correc-
c                  tions
c
c                               b.drittler   oct. 1989
c-----------------------------------------------------------------------
      use Quadrature_mod, only: simpson
      IMPLICIT NONE

      INTEGER irmd
      INTEGER irnsd
      INTEGER ipand

C     .. Parameters ..
C     ..
C      INTEGER LMPOTD
C      PARAMETER (LMPOTD= (LPOTD+1)**2)
C      INTEGER IRMIND
C      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NSPIN
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION DRDI(IRMD,*),EPOTIN,R(IRMD,*),
C    +                 RHO2NS(IRMD,LMPOTD,2),
C    +                 VINS(IRMIND:IRMD,LMPOTD,2),VM2Z(IRMD,2),Z(*)

      DOUBLE PRECISION DRDI(IRMD),EPOTIN,R(IRMD),
     &                 RHO2NS(IRMD,(LPOT+1)**2,2),
     &                 VINS((IRMD-IRNSD):IRMD,(LPOT+1)**2,2),
     &                 VM2Z(IRMD,2)
      DOUBLE PRECISION Z

      INTEGER IPAN,IRCUT(0:IPAND),IRMIN,IRWS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,R2RHOD,R2RHOU,RFPI,TEMP,ZZOR
      INTEGER I,IC,IPAN1,IPOTD,IPOTU,IRC1,IRMIN1,IRS1,L1,LM,M1
C     ..
C     .. Local Arrays ..(automatic Fortran arrays)
      DOUBLE PRECISION ER(IRMD)
      INTEGER IRCUTM(0:IPAND)

      PI = 4.0D0*ATAN(1.0D0)
      RFPI = SQRT(4.0D0*PI)

        IPAN1 = IPAN
        IRC1 = IRCUT(IPAN1)

        IF (IPAN1.GT.1) THEN
          IRS1 = IRCUT(1)
        ELSE
          IRS1 = IRWS
        END IF
c
        IF (NSPIN.EQ.1) THEN
          IPOTU = 1
          IPOTD = 1
        ELSE
          IPOTU = 1
          IPOTD = 2
        END IF

        DO 10 I = 1,IRS1
c
c---> calculate charge density times input potential
c
          R2RHOU = (RHO2NS(I,1,1)-RHO2NS(I,1,NSPIN))/2.0D0
          R2RHOD = (RHO2NS(I,1,1)+RHO2NS(I,1,NSPIN))/2.0D0
          ER(I) = - (R2RHOU*VM2Z(I,IPOTU)+R2RHOD*VM2Z(I,IPOTD))*RFPI
   10   CONTINUE                    ! I = 1,IRS1
c
c--->  remember the form of vm2z between mt sphere and rirc
c
        IF (IPAN1.GT.1) THEN
          DO 20 I = IRS1 + 1,IRC1
            R2RHOU = (RHO2NS(I,1,1)-RHO2NS(I,1,NSPIN))/2.0D0
            R2RHOD = (RHO2NS(I,1,1)+RHO2NS(I,1,NSPIN))/2.0D0
            ZZOR = 2.0D0*Z/R(I)
            ER(I) = - (R2RHOU* (VM2Z(I,IPOTU)-ZZOR)+
     +              R2RHOD* (VM2Z(I,IPOTD)-ZZOR))*RFPI
   20     CONTINUE
        END IF
c
c--->   now integrate er to get epotin
c
        IF (IPAN1.GT.1) THEN
c         CALL SIMPK(ER,TEMP,IPAN,IRCUT,DRDI)
          TEMP = simpson(ER, IPAN, IRCUT, DRDI)
        ELSE
c         CALL SIMP3(ER,TEMP,1,IRS1,DRDI)
          TEMP = simpson(ER, 1, IRS1, DRDI)
        END IF

        EPOTIN = TEMP
c
c--->   add non spher. contribution in case of non spher. input potential
c
        IRMIN1 = IRMIN

C STRANGE, TODO: CHECK E.R.
        IF (IRMIN1.LE.IRS1) THEN

          IRCUTM(0) = IRMIN1 - 1
          DO 40 IC = 1,IPAN1
            IRCUTM(IC) = IRCUT(IC)
   40     CONTINUE
c
          DO 80 L1 = 1,LPOT

            DO 50 I = 1,IRMD
              ER(I) = 0.0D0
   50       CONTINUE

            DO 70 M1 = -L1,L1
              LM = L1* (L1+1) + M1 + 1
              DO 60 I = IRMIN1,IRC1
c
c---> calculate charge density times potential
c
                R2RHOU = (RHO2NS(I,LM,1)-
     +                   RHO2NS(I,LM,NSPIN))/2.0D0
                R2RHOD = (RHO2NS(I,LM,1)+
     +                   RHO2NS(I,LM,NSPIN))/2.0D0
                ER(I) = ER(I) - R2RHOU*VINS(I,LM,IPOTU) -
     +                  R2RHOD*VINS(I,LM,IPOTD)
   60         CONTINUE
   70       CONTINUE
c           CALL SIMPK(ER,TEMP,IPAN1,IRCUTM,DRDI)
            TEMP = simpson(ER, IPAN1, IRCUTM, DRDI)
c
            EPOTIN = EPOTIN + TEMP
!ART            ENS(L1,IATYP) = TEMP
   80     CONTINUE

        END IF                    ! (IRMIN1.LE.IRS1)

      RETURN
      END
