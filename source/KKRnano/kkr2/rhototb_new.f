C           new call in main2.f
C           call RHOTOTB_NEW(NSPIN,RHO2NS,RHOCAT, &
C                        DRDI(:,I1),IRCUT(:,I1), &
C                        LPOT,NFU(ICELL),LLMSP(1,ICELL),THETAS(:,:,ICELL),IPAN(I1), &
C                        CATOM, &
C                        irmd, irid, ipand, nfund)


C>    @param[out] CATOM  CATOM(1) charge, CATOM(2) magn. moment
C>    @param[in,out] RHO2NS is modified on output! - core charge added
      SUBROUTINE RHOTOTB_NEW(NSPIN,RHO2NS,RHOC,
     +                   DRDI,
     +                   IRCUT,LPOT,NFU,LLMSP,THETAS,IPAN,
     +                   CATOM,
C                        new input parameters after inc.p removal
     &                   irmd, irid, ipand, nfund)
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

      INTEGER irmd
      INTEGER irid
      INTEGER ipand
      INTEGER nfund

C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAX+1)**2
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NSPIN

      DOUBLE PRECISION DRDI(IRMD)
      DOUBLE PRECISION RHO2NS(IRMD,(LPOT+1)**2,2)
      DOUBLE PRECISION RHOC(IRMD,2)
      DOUBLE PRECISION THETAS(IRID,NFUND)
      DOUBLE PRECISION CATOM(NSPIN)

      INTEGER IPAN,IRCUT(0:IPAND),LLMSP(NFU),NFU

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RFPI
      INTEGER I,IFUN,IPAN1,IRC1,IRS1,ISPIN,
     +        LM,LMPOT
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RHO(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT

C     ..
      RFPI = SQRT(16.0D0*ATAN(1.0D0))
      LMPOT = (LPOT+1)**2

C ======================================================================
c
       IPAN1 = IPAN
       IRS1 = IRCUT(1)
       IRC1 = IRCUT(IPAN1)

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
         DO IFUN = 1,NFU
           LM = LLMSP(IFUN)
           IF (LM.LE.LMPOT) THEN
             DO I = IRS1 + 1,IRC1
               RHO(I) = RHO(I) + RHO2NS(I,LM,ISPIN)*
     +                  THETAS(I-IRS1,IFUN)
             END DO 
           END IF
         END DO
c
c--->       integrate over circumscribed sphere
c
         CALL SIMPK(RHO,CATOM(ISPIN),IPAN1,
     +                            IRCUT, DRDI)

       END DO                      ! ISPIN = 1,NSPIN
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      END
