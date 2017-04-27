      MODULE MOD_FORCEH
      CONTAINS
      SUBROUTINE FORCEH(CMOM,FLMH,LMAX,LMPOTD,NSPIN,NATOM,
     +                  CELL,DENSITY,VPOT,
     +                  IRMD, ZATOM,INS) !,NCORE,IRWS,Z)
      USE TYPE_DENSITY
      USE TYPE_CELL
      USE MOD_SIMP3

        IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m with hellmann - feynman theorem
c     from a given non spherical charge density at the nucleus site r
c
 
c-----------------------------------------------------------------------
C     .. Parameters ..
!       include 'inc.p'
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMPOTD,NATOM,IRMD,INS
      INTEGER LMAX,NEND,NSPIN,NSTART
      TYPE(DENSITY_TYPE) :: DENSITY(NATOM)
      TYPE(CELL_TYPE)    :: CELL(NATOM)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CMOM(LMPOTD,NATOM),FLMH(-1:1,NATOM),
!      +       R(IRMD,*),
     +       R2RHO(IRMD,LMPOTD,NATOM,NSPIN),
     +       VPOT(IRMD,LMPOTD,NSPIN,NATOM),ZATOM(NATOM)
!       INTEGER IRWS(NATOM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RWS,VINT1
      INTEGER I,IATOM,IPOT,IRWS1,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FLM(-1:1,2),V1(IRMD)
C     ..
C     .. External Subroutines ..
!       EXTERNAL SIMP3
C     ..
C     .. Save statement ..
!       SAVE PI
C     ..
c
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
!       PI = 4.D0*ATAN(1.D0)
      IF (LMAX.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
c
c---> loop over the rep. atoms
c
      DO 10 IATOM = 1,NATOM
c
c---> reading the right Wigner-S. radius
c
         IF (INS==1) THEN
           IRWS1 = CELL(IATOM)%NRCORE !IRWS(IATOM)
         ELSE
           IRWS1 = CELL(IATOM)%NRMAX !IRWS(IATOM)
         END IF
!          IRWS1 = CELL(NATOM)%NRMAX !IRWS(IATOM)
         RWS = CELL(NATOM)%RMESH(IRWS1)
c
c---> determine the right potential numbers
c
!          IPOT = NSPIN* (IATOM-1) + 1
 
         DO 20 M = -1,1
            LM = 2 + M + 1
c
            V1(1) = 0.0D0
            DO 30 I = 2,IRWS1
               V1(I) = DENSITY(IATOM)%RHO2NS(I,LM,1)* 
     +         (CELL(IATOM)%RMESH(I)** (-2.0D0))
   30       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,CELL(IATOM)%DRMESHDI)
c
            FLM(M,1) = 2.0D0*VINT1
c
c---> use coulomb potential to determine extra atomic contribution
c
            FLM(M,2) = VPOT(IRWS1,LM,1,IATOM)* (3.0D0/ 
     +                 (4.0D0*PI*RWS)) -
     +                 2.0D0*CMOM(LM,IATOM)/ (RWS**3)
c
c---> total Hellman-Feynman force
c
            FLMH(M,IATOM) = (FLM(M,1)+FLM(M,2))*ZATOM(IATOM)
   20    CONTINUE
   10 CONTINUE
c
c
 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least ')
 
      END SUBROUTINE
      END MODULE MOD_FORCEH
