c 13.10.95 ***************************************************************
      SUBROUTINE RHOTOTB(IPF,NATYP,NAEZ,NSPIN,RHO2NS,RHOC,RHOORB,
     +                   Z,DRDI,IRWS,
     +                   IRCUT,LPOT,NFU,LLMSP,THETAS,NTCELL,KSHAPE,IPAN,
     +                   CHRGNT,ITC,NSHELL,NOQ,CONC,KAOEZ,CATOM)
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
c     total orbital moment within the WS sphere is also calculated
c     in the relativistic case; orbital density is normalised in the
c     same way as the charge density. v.popescu march 2002
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CHRGNT
      INTEGER ITC,IPF,KSHAPE,LPOT,NATYP,NSPIN,NAEZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 RHOC(IRMD,*),THETAS(IRID,NFUND,*),Z(*)
      DOUBLE PRECISION CATOM(NATYPD,2*KREL+(1-KREL)*NSPIND)
      DOUBLE PRECISION CONC(NATYPD)
C----------------------------------- orbital density and moment
      DOUBLE PRECISION  RHOORB(IRMD*KREL + (1-KREL),NATYPD)
      DOUBLE PRECISION  OMOM(NATYPD)
C----------------------------------------------------------

      INTEGER IPAN(*),IRCUT(0:IPAND,*),IRWS(*),LLMSP(NATYPD,*),NFU(*),
     +        NSHELL(0:NSHELD),NTCELL(*),KAOEZ(NATYPD,NAEZD+NEMBD),
     +        NOQ(NAEZD)

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DIFF,FACTOR,RFPI,SUM,TOTSMOM,TOTOMOM,SUMO
      INTEGER I,I1,IATYP,ICELL,IFUN,IPAN1,IPOTD,IPOTU,IRC1,IRS1,ISPIN,
     +        LM,LMPOT,IQEZ,IOEZ
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RHO(IRMD)
      DOUBLE PRECISION CSITE(NAEZD,2*KREL+(1-KREL)*NSPIND)
      DOUBLE PRECISION MUOSITE(KREL*NAEZD+(1-KREL))
c
c
      LOGICAL OPT
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3,SIMPK,OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
      RFPI = SQRT(16.0D0*ATAN(1.0D0))
      LMPOT = (LPOT+1)**2

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c---> loop over atomic sites
c
      DO IQEZ = 1, NAEZ 
C ======================================================================
c
c---> loop over atoms located on IQEZ
c
        DO ISPIN=1,NSPIN
          CSITE(IQEZ,ISPIN) = 0.0D0
        END DO

        IF (KREL.EQ.1) MUOSITE(IQEZ) = 0.0D0

        DO IOEZ =1, NOQ(IQEZ)

          IATYP = KAOEZ(IOEZ,IQEZ)
c
c--->   determine the right potential numbers for rhoc
c
          IF (NSPIN.EQ.2) THEN
            IPOTD = 2*IATYP - 1
            IPOTU = 2*IATYP
            FACTOR = 1.0D0
          ELSE
            IPOTD = IATYP
            IPOTU = IATYP
            FACTOR = 0.5D0
          END IF

          IF (KSHAPE.NE.0) THEN
            IPAN1 = IPAN(IATYP)
            IRS1 = IRCUT(1,IATYP)
            IRC1 = IRCUT(IPAN1,IATYP)
          ELSE
            IRS1 = IRWS(IATYP)
            IRC1 = IRS1
          END IF
c

C-----------------------------------------------------------------------
          DO I = 2,IRS1
c
c--->     convert core density
c
            SUM = (RHOC(I,IPOTD)+RHOC(I,IPOTU))*FACTOR/RFPI
            DIFF = (RHOC(I,IPOTU)-RHOC(I,IPOTD))/RFPI
c
c--->     add this to the lm=1 component of rho2ns
c
            RHO2NS(I,1,IATYP,1) = RHO2NS(I,1,IATYP,1) + SUM
            RHO2NS(I,1,IATYP,NSPIN) = RHO2NS(I,1,IATYP,NSPIN) + DIFF

          END DO
          
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
c
c--->   calculate  charge and moment of the atom
c
          DO ISPIN = 1,NSPIN
c
            IF (KSHAPE.EQ.0) THEN
c
c--->       integrate over wigner seitz sphere - no shape correction
c
              CALL SIMP3(RHO2NS(1,1,IATYP,ISPIN),
     &                   SUM,1,IRS1,DRDI(1,IATYP))
c
c--->       the result has to be multiplied by sqrt(4 pi)
c           (4 pi for integration over angle and 1/sqrt(4 pi) for
c           the spherical harmonic y(l=0))
c
              SUM = SUM*RFPI
            ELSE                      ! (KSHAPE.EQ.0)
c
c--->       convolute charge density with shape function to get the
c           charge in the exact cell - if kshape .gt. 0
c
              ICELL = NTCELL(IATYP)
  
              DO I = 1,IRS1
                RHO(I) = RHO2NS(I,1,IATYP,ISPIN)*RFPI
              END DO
C
              DO I = IRS1 + 1,IRC1
                RHO(I) = 0.0D0
              END DO
C
              DO IFUN = 1,NFU(ICELL)
                LM = LLMSP(ICELL,IFUN)
                IF (LM.LE.LMPOT) THEN
                  DO I = IRS1 + 1,IRC1
                    RHO(I) = RHO(I) + RHO2NS(I,LM,IATYP,ISPIN)*
     +                       THETAS(I-IRS1,IFUN,ICELL)
                  END DO 
                END IF
              END DO
c
c--->       integrate over circumscribed sphere
c
              CALL SIMPK(RHO,SUM,IPAN1,IRCUT(0,IATYP),DRDI(1,IATYP))

            END IF                    ! (KSHAPE.EQ.0)

            CATOM(IATYP,ISPIN) = SUM
            CSITE(IQEZ,ISPIN) = CSITE(IQEZ,ISPIN) 
     &                          + CATOM(IATYP,ISPIN) *CONC(IATYP)

            IF (ISPIN.NE.1) THEN
c
c ---> calculate orbital moment (ASA) and add it to the total
c
               IF ((KREL.EQ.1) .AND. (KSHAPE.EQ.0)) THEN
                  CALL SIMP3(RHOORB(1,IATYP),SUMO,1,IRS1,DRDI(1,IATYP))
                  SUMO = SUMO*RFPI
                  OMOM(IATYP) = SUMO
                  MUOSITE(IQEZ) = MUOSITE(IQEZ) 
     &                            + OMOM(IATYP) * CONC(IATYP)
               END IF

               IF (KSHAPE.NE.0) THEN
                  WRITE (IPF,FMT=9010) SUM
               ELSE
                  WRITE (IPF,FMT=9050) SUM
                  IF (KREL.EQ.1) THEN
                     WRITE (IPF,FMT=9051) OMOM(IATYP)
                     WRITE (IPF,FMT=9052) SUM+OMOM(IATYP)
                  END IF
               END IF
            ELSE                      ! (ISPIN.NE.1)

              IF (KSHAPE.NE.0) THEN
                WRITE (IPF,FMT=9000) IATYP,SUM
              ELSE
                WRITE (IPF,FMT=9040) IATYP,SUM
              END IF

            END IF                    ! (ISPIN.NE.1)

          END DO                      ! ISPIN = 1,NSPIN
C
C-----------------------------------------------------------------------
          IF (IOEZ.NE.NOQ(IQEZ)) WRITE(6,'(2X,77(1H-))')
        END DO
C                                     ! IOEZ = 1, NOQ(IQEZ)
C ======================================================================

        IF (NOQ(IQEZ).GT.1) THEN
          WRITE(6,'(2X,77(1H=))')
          WRITE(IPF,FMT=9071) IQEZ,CSITE(IQEZ,1)
          IF (NSPIN.EQ.2) THEN
              WRITE(IPF,FMT=9072) CSITE(IQEZ,NSPIN)
              IF (KREL.EQ.1) THEN
                  WRITE(IPF,FMT=9073) MUOSITE(IQEZ)
                  WRITE(IPF,FMT=9074) CSITE(IQEZ,NSPIN)+MUOSITE(IQEZ)
              END IF
          END IF
          IF (IQEZ.NE.NAEZ) WRITE(6,'(2X,77(1H=))')
        ELSE
          IF (IQEZ.NE.NAEZ) WRITE(6,'(2X,77(1H=))')
        END IF

      END DO
      WRITE(6,*)
C                                     ! IQEZ = 1, NAEZ
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      CHRGNT = 0.0D0
      DO I1 = 1,NATYP
        CHRGNT = CHRGNT 
     &           + DBLE(NSHELL(I1))*(CATOM(I1,1) - Z(I1))*CONC(I1)
      END DO

      WRITE(6,'(79(1H+))')
      WRITE (IPF,FMT=9020) ITC,CHRGNT

      IF (NSPIN.EQ.2) THEN
        TOTSMOM = 0.0D0
        IF (KREL.EQ.1) TOTOMOM = 0.0D0
        DO I1 = 1,NATYP
          TOTSMOM = TOTSMOM 
     &           + DBLE(NSHELL(I1))*CATOM(I1,NSPIN)*CONC(I1)
          IF (KREL.EQ.1) TOTOMOM = TOTOMOM 
     &         + DBLE(NSHELL(I1))*OMOM(I1)*CONC(I1)
        END DO

        IF (KREL.EQ.0) THEN
          WRITE (IPF,FMT=9030) TOTSMOM
        ELSE
          WRITE (IPF,FMT=9030) TOTSMOM+TOTOMOM
          WRITE (IPF,FMT=9031) TOTSMOM
          WRITE (IPF,FMT=9032) TOTOMOM
        END IF
      END IF
      WRITE (IPF,*)

      RETURN

 9000 FORMAT ('  Atom ',I4,' charge in wigner seitz cell =',f10.6)
 9010 FORMAT (7X,'spin moment in wigner seitz cell =',f10.6)
 9040 FORMAT ('  Atom ',I4,' charge in wigner seitz sphere =',f10.6)
 9050 FORMAT (7X,'spin moment in wigner seitz sphere =',f10.6)
 9051 FORMAT (7X,'orb. moment in wigner seitz sphere =',f10.6)
 9052 FORMAT (7X,'total magnetic moment in WS sphere =',f10.6)
 9020 FORMAT ('      ITERATION',I4,
     &     ' charge neutrality in unit cell = ',f12.6)
 9030 FORMAT ('                   ',
     &     ' TOTAL mag. moment in unit cell = ',f12.6)
 9031 FORMAT ('                   ',
     &     '           spin magnetic moment = ',f12.6)
 9032 FORMAT ('                   ',
     &     '        orbital magnetic moment = ',f12.6)
 9071 FORMAT ('      Site ',i3,' total charge =',f10.6)
 9072 FORMAT ('         ',' total spin moment =',f10.6)
 9073 FORMAT ('         ',' total orb. moment =',f10.6)
 9074 FORMAT ('      total magnetic moment =',f10.6)

      END
