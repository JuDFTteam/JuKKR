      SUBROUTINE RHOTOTB(NAEZ,IATYP,NSPIN,RHO2NS,RHOC,
     +                   Z,DRDI,
     +                   IRCUT,LPOT,NFU,LLMSP,THETAS,ICELL,IPAN,
     +                   CATOM)
      implicit none
c ************************************************************************
c     add core and valence density expanded in spherical harmonics
c         ( convention see subroutine rholm )
c     in the paramagnetic case (nspin=1) the core valence charge times
c         r**2 is add to the valence charge density times r**2
c         then only rho2ns(irmd,lmxtsq,natypd,1) is used .
c     in the spin-polarized case (nspin=2) the spin-splitted core
c         charge density times r**2 is converted into core charge
c         density times r**2 and core spin density times r**2 .
c         then these parts are added to corresponding parts of
c         the valence densities times r**2 , that are rho2ns(...,1)
c         which contains the charge density  and rho2ns(...,2) which
c         contains in that case the spin density .
c             (see notes by b.drittler)
c
c     attention : the core density is spherically averaged and multi-
c                 plied by 4 pi. therefore the core density is only
c                 added to l=0 part .
c
c                               b.drittler   nov. 1989
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NSPIN,NAEZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),RHO2NS(IRMD,LMPOTD,2),
     +                 RHOC(IRMD,2),THETAS(IRID,NFUND,*),Z(*)
      DOUBLE PRECISION CATOM(NSPIND)

      INTEGER IPAN(*),IRCUT(0:IPAND,*),LLMSP(*),NFU(*)

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RFPI
      INTEGER I,I1,IATYP,ICELL,IFUN,IPAN1,IRC1,IRS1,ISPIN,
     +        LM,LMPOT
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RHO(IRMD)
c
      LOGICAL OPT
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMPK,OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
      RFPI = SQRT(16.0D0*ATAN(1.0D0))
      LMPOT = (LPOT+1)**2

C ======================================================================
c
       IPAN1 = IPAN(IATYP)
       IRS1 = IRCUT(1,IATYP)
       IRC1 = IRCUT(IPAN1,IATYP)

C-----------------------------------------------------------------------
       IF(NSPIN.EQ.2) THEN
       DO I = 2,IRS1
         RHO2NS(I,1,1) = RHO2NS(I,1,1) + (RHOC(I,1)+RHOC(I,2))/RFPI
         RHO2NS(I,1,2) = RHO2NS(I,1,2) + (RHOC(I,2)-RHOC(I,1))/RFPI

       END DO
       ELSE
       DO I = 2,IRS1
         RHO2NS(I,1,1) = RHO2NS(I,1,1) + RHOC(I,1)/RFPI
       END DO
       END IF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
c
c--->   calculate  charge and moment of the atom
c
       DO ISPIN = 1,NSPIN
c
c--->       convolute charge density with shape function to get the
c           charge in the exact cell
c
         DO I = 1,IRS1
           RHO(I) = RHO2NS(I,1,ISPIN)*RFPI
         END DO
C
         DO I = IRS1 + 1,IRC1
           RHO(I) = 0.0D0
         END DO
C
         DO IFUN = 1,NFU(ICELL)
           LM = LLMSP(IFUN)
           IF (LM.LE.LMPOT) THEN
             DO I = IRS1 + 1,IRC1
               RHO(I) = RHO(I) + RHO2NS(I,LM,ISPIN)*
     +                  THETAS(I-IRS1,IFUN,ICELL)
             END DO 
           END IF
         END DO
c
c--->       integrate over circumscribed sphere
c
         CALL SIMPK(RHO,CATOM(ISPIN),IPAN1,
     +                            IRCUT(0,IATYP),DRDI(1,IATYP))

       END DO                      ! ISPIN = 1,NSPIN
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      RETURN

      END
