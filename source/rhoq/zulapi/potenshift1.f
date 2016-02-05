      SUBROUTINE POTENSHIFT1(VISP,VINS,NATYP,NSPIN,
     +     IRCUT,IRC,IRMIN,NTCELL,IMAXSH,ILM,IFUNM,LMSP,LMPOT,GSH,
     +     THETAS,RMESH,KSHAPE,VSHIFT)
      implicit none
c Adds a constant (=VSHIFT) to the potentials of atoms
c
c Parameters:
      include 'inc.p'
      INTEGER NPOTD,LMPOTD,LMXSPD,IRMIND
      PARAMETER (NPOTD=NSPIND*NATYPD,LMPOTD= (LPOTD+1)**2
     &     ,LMXSPD= (2*LPOTD+1)**2,IRMIND=IRMD-IRNSD)
c Input
      INTEGER KSHAPE,LMPOT,NATYP,NSPIN
      INTEGER IRCUT(0:IPAND,NATYPD),IRC(NATYPD),NTCELL(NATYPD)
     &     ,IMAXSH(0:LMPOTD),ILM(NGSHD,3),IFUNM(NATYPD,LMXSPD)
     &     ,LMSP(NATYPD,LMXSPD),IRMIN(NATYPD)
      REAL*8 GSH(NGSHD),THETAS(IRID,NFUND,NCELLD),RFPI,PI,FPI
     &     ,RMESH(IRMD,NATYPD)
      REAL*8 THESME(IRID,NFUND,NCELLD)
      REAL*8 VSHIFT


c Input/Output:
      REAL*8 VISP(IRMD,NPOTD),VINS(IRMIND:IRMD,LMPOTD,NSPOTD)

c Inside
      REAL*8 PSHIFTLMR(IRMD,LMPOTD),PSHIFTR(IRMD)
      INTEGER ISPIN,IH,IPOT,IR,LM,IMT1,IRC1,IRMIN1
      INTEGER IER
      CHARACTER*80 UIO


      THESME=0d0
      PI=4.0D0*ATAN(1.0D0)
      FPI=4.0D0*PI
      RFPI=SQRT(FPI)

      DO IH = 6,16

         IMT1 = IRCUT(1,IH)
         IRC1 = IRC(IH)
         IRMIN1 = IRMIN(IH)

         DO ISPIN = 1,NSPIN

            WRITE (6,*) 'SHIFTING OF THE POTENTIALS OF ATOM',IH,
     &           ' BY', VSHIFT, 'RY.'
            IPOT = NSPIN * (IH-1) + ISPIN

            CALL RINIT(IRMD*LMPOTD,PSHIFTLMR)
            CALL RINIT(IRMD,PSHIFTR)
            DO IR = 1,IRC1
               PSHIFTLMR(IR,1) = VSHIFT
            ENDDO

            IF (KSHAPE.EQ.0) THEN ! ASA

               DO IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               END DO

            ELSE                ! Full-potential


               CALL CONVOL(IMT1,IRC1,NTCELL(IH),
     &              IMAXSH(LMPOT),ILM,IFUNM,LMPOT,GSH,
     &              THETAS,THESME,0.d0,RFPI,
     &              RMESH(1,IH),PSHIFTLMR,PSHIFTR,LMSP)


               DO IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               ENDDO


               DO LM = 2,LMPOT
                  DO IR = IRMIN1,IRC1
                VINS(IR,LM,IPOT)=VINS(IR,LM,IPOT)+PSHIFTLMR(IR,LM)*RFPI
                  ENDDO
               ENDDO

            END IF              ! (KSHAPE.EQ.0)


         END DO

      END DO
c      DO IR=1,IRMD
c       DO LM=1,NPOTD
c        WRITE(55,'(1e17.9)') VISP(IR,LM)
c       ENDDO
c      ENDDO
c      DO IR=IRMIND,IRMD
c       DO LM=1,LMPOTD
c        DO IH=1,NSPOTD
c         WRITE(56,'(1e17.9)') VINS(IR,LM,IH)
c        ENDDO
c       ENDDO
c      ENDDO 
      END
