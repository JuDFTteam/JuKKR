  !-------------------------------------------------------------------------------
  !> Summary: Calculates force on nucleus with core contribution (Coulomb contribution)
  !> Author: 
  !>
  !> Calculates the force on nucleus m
  !> from a given non spherical charge density at the nucleus site r
  !> with core correction (coulomb contribution)
  !-------------------------------------------------------------------------------
      MODULE MOD_FORCE
      CONTAINS
  !-------------------------------------------------------------------------------
  !> Summary: Calculates force on nucleus with core contribution (Coulomb contribution)
  !> Author: 
  !> Category: KKRimp, physical-observables
  !> Deprecated: False
  !>
  !> Calculates the force on nucleus m
  !> from a given non spherical charge density at the nucleus site r
  !> with core correction (coulomb contribution)
  !-------------------------------------------------------------------------------
      SUBROUTINE FORCE(FLM,FLMC,LMAX,NSPIN,NATOM,VPOT,DENSITY,CELL, &
                       IRMD,LMPOTD,INS)
      USE MOD_SIMP3
      USE TYPE_DENSITY
      USE TYPE_CELL
      IMPLICIT NONE
!C     .. Parameters ..
!       include 'inc.p'
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
!C     ..
!C     .. Scalar Arguments ..
      INTEGER :: IRMD,LMPOTD,INS
      INTEGER LMAX,NATOM,NSPIN
      TYPE(CELL_TYPE) :: CELL(NATOM)
      TYPE(DENSITY_TYPE) :: DENSITY(NATOM)
!C     ..
!C     .. Array Arguments ..
!       DOUBLE PRECISION DRDI(IRMD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMD,*),
!      +       RHOC(IRMD,*),VPOT(IRMD,LMPOTD,*)
      DOUBLE PRECISION FLM(-1:1,NATOM),FLMC(-1:1,NATOM)
      DOUBLE PRECISION VPOT(IRMD,LMPOTD,NSPIN,NATOM)
!       INTEGER IRWS(*)
!C     ..
!C     .. Local Scalars ..
      DOUBLE PRECISION DV,FAC,RWS,VINT1
      INTEGER I,IATOM,IRWS1,ISPIN,LM,M
!C     ..
!C     .. Local Arrays ..
      DOUBLE PRECISION FLMH(-1:1,NATOM),V1(IRMD)
!C     ..
!C     .. External Subroutines ..
!       EXTERNAL SIMP3
!C     ..
!C     .. Save statement ..
!       SAVE PI
!C     ..
!c
!C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DSQRT
!C     ..
!       PI = 4.D0*ATAN(1.D0)
      FAC = DSQRT((4.0D0*PI)/3.0D0)
      IF (LMAX.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
!c
!c---> loop over rep. atoms
!c
      DO IATOM = 1,NATOM
!c
!c
         IF (INS==1) THEN
           IRWS1 = CELL(IATOM)%NRCORE !IRWS(IATOM)
         ELSE
           IRWS1 = CELL(IATOM)%NRMAX !IRWS(IATOM)
         END IF

!          IRWS1 = CELL(IATOM)%NRMAX !IRWS(IATOM)
         RWS =   CELL(IATOM)%RMESH(IRWS1)
!c
!c
 
         DO M = -1,1
            LM = 2 + M + 1
!c
!c---> initialize v1
!c
            DO  I = 1,IRWS1
               V1(I) = 0.0D0
            END DO
!c
            DO  ISPIN = 1,NSPIN
!c
!c---> determine the right potential numbers
!c
!                IPOT = NSPIN* (IATOM-1) + ISPIN
!c
!c---> determine the derivative of the potential using a 5-point formular
!c
               DV = (-3.0D0*VPOT(1,LM,ISPIN,IATOM)-10.0D0*VPOT(2,LM,ISPIN,IATOM)+ &
                  18.0D0*VPOT(3,LM,ISPIN,IATOM)-6.0D0*VPOT(4,LM,ISPIN,IATOM)+VPOT(5,LM,ISPIN,IATOM))/ &
                    (12.0D0*CELL(IATOM)%DRMESHDI(2))
!c
               V1(2) = DENSITY(IATOM)%RHOC(2,ISPIN)* (2.0D0*VPOT(2,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(2)+DV)/&
                       (4.0D0*PI) + V1(2)
!c
               DO  I = 3,IRWS1 - 2
!c
                  DV = (VPOT(I-2,LM,ISPIN,IATOM)-VPOT(I+2,LM,ISPIN,IATOM)+ &
                       8.0D0* (VPOT(I+1,LM,ISPIN,IATOM)-VPOT(I-1,LM,ISPIN,IATOM)))/ &
                       (12.0D0*CELL(IATOM)%DRMESHDI(I))
!c
                  V1(I) = DENSITY(IATOM)%RHOC(I,ISPIN)* (2.0D0*VPOT(I,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(I)+ &
                          DV)/ (4.0D0*PI) + V1(I)
               END DO !I
!c
               DV = (-VPOT(IRWS1-4,LM,ISPIN,IATOM)+6.0D0*VPOT(IRWS1-3,LM,ISPIN,IATOM)- &
                    18.0D0*VPOT(IRWS1-2,LM,ISPIN,IATOM)+10.0D0*VPOT(IRWS1-1,LM,ISPIN,IATOM)+ &
                   3.0D0*VPOT(IRWS1,LM,ISPIN,IATOM))/ (12.0D0*CELL(IATOM)%DRMESHDI(IRWS1-1))
               V1(IRWS1-1) = DENSITY(IATOM)%RHOC(IRWS1-1,ISPIN)* &
                             (2.0D0*VPOT(IRWS1-1,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(IRWS1-1)+ &
                             DV)/ (4.0D0*PI) + V1(IRWS1-1)
!c
               DV = (3.0D0*VPOT(IRWS1-4,LM,ISPIN,IATOM)-16.0D0*VPOT(IRWS1-3,LM,ISPIN,IATOM)+ &
                    36.0D0*VPOT(IRWS1-2,LM,ISPIN,IATOM)-48.0D0*VPOT(IRWS1-1,LM,ISPIN,IATOM)+ &
                    25.0D0*VPOT(IRWS1,LM,ISPIN,IATOM))/ (12.0D0*CELL(IATOM)%DRMESHDI(IRWS1))
!c
               V1(IRWS1) = DENSITY(IATOM)%RHOC(IRWS1,ISPIN)* &
                           (2.0D0*VPOT(IRWS1,LM,ISPIN,IATOM)/CELL(IATOM)%RMESH(IRWS1)+DV)/ &
                           (4.0D0*PI) + V1(IRWS1)
            END DO !ISPIN
!c
!c---> integrate with simpson subroutine
!c
            CALL SIMP3(V1,VINT1,1,IRWS1,CELL(IATOM)%DRMESHDI(1))
!c
            FLMH(M,IATOM) = FAC*FLM(M,IATOM)
            FLMC(M,IATOM) = -FAC*VINT1
            FLM(M,IATOM) = FLMH(M,IATOM) + FLMC(M,IATOM)
!c
 
         END DO !M
!c
!c
      END DO !NATOM
!c
 9000 FORMAT (13x,'error stop in subroutine force :', &
             ' the charge density has to contain non spherical', &
             ' contributions up to l=1 at least ')
 
      END SUBROUTINE
      END MODULE MOD_FORCE
