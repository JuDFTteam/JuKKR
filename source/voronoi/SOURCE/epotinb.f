      MODULE MOD_EPOTINB
        CONTAINS

c 13.10.95 ***************************************************************
      SUBROUTINE EPOTINB(EPOTIN,NSPIN,NATOM,VM2Z,INS,
     +                   LMAXATOM,ZATOM,CELL,DENSITY,
     +                   IPAND, IRMD,LPOTD)
      USE TYPE_DENSITY
      USE TYPE_CELL
      USE MOD_SIMPK
      USE MOD_SIMP3
      IMPLICIT NONE
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
C     .. Parameters ..
!       include 'inc.p'
C     ..
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
!       INTEGER IRMIND
!       PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      INTEGER INS,LPOT,NATOM,NSPIN,IPAND,IRMD,LPOTD
      INTEGER LMAXATOM(NATOM)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EPOTIN(NATOM),ZATOM(NATOM),
     +                VM2Z(IRMD,(LPOTD+1)**2,NSPIN,NATOM) !,VM2Z(IRMD,*)
      TYPE(CELL_TYPE)        :: CELL(NATOM)
      TYPE(DENSITY_TYPE)     :: DENSITY(NATOM)

!       DOUBLE PRECISION DRDI(IRMD,*),EPOTIN(*),R(IRMD,*),
!      +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),
!      +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*)
!       INTEGER IPAN(*),IRCUT(0:IPAND,*),IRMIN(*),IRWS(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R2RHOD,R2RHOU,TEMP,ZZOR,RFPI
      INTEGER I,IATOM,IC,IPAN1,IPOTD,IPOTU,IRC1,IRMIN1,IRS1,L1,LM,M1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ENS(0:LPOTD,NATOM),ER(IRMD)
      INTEGER IRCUTM(0:IPAND)
C     ..
C     .. External Subroutines ..
!       EXTERNAL SIMP3,SIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
!       PI = 4.0D0*ATAN(1.0D0)
      RFPI = SQRT(4.0D0*PI)

      DO 90 IATOM = 1,NATOM

        IPAN1 = CELL(IATOM)%NPAN
        IRC1 = CELL(IATOM)%NRCUT(IPAN1)

        IF (IPAN1.GT.1) THEN
          IRS1 = CELL(IATOM)%NRCUT(1)
        ELSE
          IRS1 = CELL(IATOM)%NRMAX !IRS1 = IRWS(IATOM)
        END IF
c
!         IF (NSPIN.EQ.1) THEN
!           IPOTU = IATOM
!           IPOTD = IATOM
!         ELSE
!           IPOTU = 2*IATOM - 1
!           IPOTD = 2*IATOM
!         END IF

        DO 10 I = 1,IRS1
c
c---> calculate charge density times input potential
c
          R2RHOU = (DENSITY(IATOM)%RHO2NS(I,1,1)
     +              -DENSITY(IATOM)%RHO2NS(I,1,NSPIN))/2.0D0
          R2RHOD = (DENSITY(IATOM)%RHO2NS(I,1,1)
     +              +DENSITY(IATOM)%RHO2NS(I,1,NSPIN))/2.0D0

!          write(*,*) VM2Z(I,1,1,IATOM),R2RHOD

          ER(I) = - (R2RHOU*VM2Z(I,1,1,IATOM)+
     +               R2RHOD*VM2Z(I,1,NSPIN,IATOM))*RFPI
!           write(*,*) er(i)

   10   CONTINUE                    ! I = 1,IRS1
c
c--->  remember the form of vm2z between mt sphere and rirc
c
        IF (IPAN1.GT.1) THEN
          DO 20 I = IRS1 + 1,IRC1
            R2RHOU = (DENSITY(IATOM)%RHO2NS(I,1,1)
     +                -DENSITY(IATOM)%RHO2NS(I,1,NSPIN))/2.0D0
            R2RHOD = (DENSITY(IATOM)%RHO2NS(I,1,1)
     +               +DENSITY(IATOM)%RHO2NS(I,1,NSPIN))/2.0D0
            ZZOR = 2.0D0*ZATOM(IATOM)/CELL(IATOM)%RMESH(I)
            ER(I) = - (R2RHOU* (VM2Z(I,1,1,IATOM)-ZZOR)+
     +              R2RHOD* (VM2Z(I,1,NSPIN,IATOM)-ZZOR))*RFPI
   20     CONTINUE
        END IF
c
c--->   now integrate er to get epotin
c
        IF (IPAN1.GT.1) THEN
          CALL SIMPK(ER,TEMP,CELL(IATOM)%NPAN,CELL(IATOM)%NRCUT,
     +               CELL(IATOM)%DRMESHDI,IPAND)
        ELSE
          CALL SIMP3(ER,TEMP,1,IRS1,CELL(IATOM)%DRMESHDI(1))
        END IF
c
        EPOTIN(IATOM) = TEMP
        ENS(0,IATOM) = TEMP
c
c--->   add non spher. contribution in case of non spher. input potential
c
        LPOT=2*LMAXATOM(IATOM)
        DO 30 L1 = 1,LPOT
          ENS(L1,IATOM) = 0.0D0
   30   CONTINUE
c
         
!         write(*,*) 'EPOTIN',iatom,EPOTIN(IATOM)


        IF (INS.NE.0) THEN

          IRMIN1 = CELL(IATOM)%NRMIN_NS
          IF (IRMIN1.LE.IRS1) THEN

            IRCUTM(0) = IRMIN1 - 1
            DO 40 IC = 1,IPAN1
              IRCUTM(IC) = CELL(IATOM)%NRCUT(IC)
   40       CONTINUE
c
            DO 80 L1 = 1,LPOT

              DO 50 I = 1,IRMD
                ER(I) = 0.0D0
   50         CONTINUE

              DO 70 M1 = -L1,L1
                LM = L1* (L1+1) + M1 + 1
                DO 60 I = IRMIN1,IRC1
c
c---> calculate charge density times potential
c
                  R2RHOU = (DENSITY(IATOM)%RHO2NS(I,LM,1)-
     +                     DENSITY(IATOM)%RHO2NS(I,LM,NSPIN))/2.0D0
                  R2RHOD = (DENSITY(IATOM)%RHO2NS(I,LM,1)+
     +                     DENSITY(IATOM)%RHO2NS(I,LM,NSPIN))/2.0D0
                  ER(I) = ER(I) - R2RHOU*VM2Z(I,LM,1,IATOM) -
     +                    R2RHOD*VM2Z(I,LM,NSPIN,IATOM)
   60           CONTINUE
   70         CONTINUE
              CALL SIMPK(ER,TEMP,IPAN1,IRCUTM,
     +                   CELL(IATOM)%DRMESHDI,IPAND)
c
              EPOTIN(IATOM) = EPOTIN(IATOM) + TEMP
              ENS(L1,IATOM) = TEMP
   80       CONTINUE

          END IF                    ! (IRMIN1.LE.IRS1)

        END IF                      ! (INS.NE.0)

   90 CONTINUE                      ! IATOM = 1,NATYP

      END SUBROUTINE
      END MODULE MOD_EPOTINB
