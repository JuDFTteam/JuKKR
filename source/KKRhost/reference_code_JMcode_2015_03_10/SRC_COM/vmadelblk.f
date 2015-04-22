C*==vmadelblk.f    processed by SPAG 6.05Rc at 13:42 on  1 Feb 2002
      SUBROUTINE VMADELBLK(CMOM,CMINST,LMAX,NSPIN,NAEZ,NATYP,V,ZAT,R,
     &                     IRWS,IRCUT,IPAN,KSHAPE,NOQ,KAOEZ,IQAT,
     &                     CONC,CATOM,ICC,HOSTIMP,IELAST,ATOMIMP,
     &                     NATOMIMP,ALAT,VBC,VINTERS)
C **********************************************************************
C
C     calculate the madelung potentials and add these to the poten-
C     tial v  (in the spin-polarized case for each spin-direction
C     this is the same)
C     it uses the structure dependent matrices AVMAD and BVMAD which
C     are calculated once in the subroutine MADELUNG3D and saved in
C     the DA-file abvmad.unformatted                     ( may 2004)
C     the charge-moments are calculated in the subroutine vintras,
C     therefore vintras has to be called first.
C     the madelung-potential is expanded into spherical harmonics.
C     the lm-term of the potential v of the atom i is given by
C
C      v(r,lm,i) =  (-r)**l * {avmad(i,i2,lm,l'm')*cmom(i2,l'm')
C                                               +bvmad(i,i2,lm)*z(i2)}
C
C     summed over i2 (all atoms) and l'm'
C     (see notes by b.drittler)
C
C                               b.drittler   nov. 1989
C
C     adopted for the case of more atoms on the same site, summation is
C     done over the occupants of that site, the charge is weighted with
C     the appropriate concentration of the occupant  v.popescu feb. 2002
C     
C     impurity-program adopted feb. 2004 (according to n.papanikalou)
C     
C **********************************************************************
      IMPLICIT NONE
      INCLUDE 'inc.p'
C     ..
C     .. PARAMETER definitions
      INTEGER LMPOTD
      PARAMETER (LMPOTD=(LPOTD+1)**2)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
      INTEGER ATOMIMP(NATOMIMPD),NATOMIMP

      DOUBLE PRECISION VBC(2)

C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX,NSPIN,NAEZ,NATYP,KSHAPE,ICC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CMOM(LMPOTD,*),CMINST(LMPOTD,*),
     &                 V(IRMD,LMPOTD,*),R(IRMD,*),ZAT(*)
      DOUBLE PRECISION CONC(NATYPD),CATOM(NATYPD)
      INTEGER IRWS(*),IRCUT(0:IPAND,*),IPAN(*),
     &        NOQ(NAEZD),KAOEZ(NATYPD,NAEZD+NEMBD),
     &        IQAT(NATYPD),HOSTIMP(0:NATYPD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,PI
      INTEGER LRECABMAD,IREC
      INTEGER I,L,LM,LM2,LMMAX,M,IO1,IO2,IPOT,IQ1,IQ2,
     &        IRS1,ISPIN,IT1,IT2,IATOM,NOQVAL
      DOUBLE PRECISION ALAT
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
      DOUBLE PRECISION VINTERS(LMPOTD,NAEZD)
      LOGICAL OPT
      INTEGER IELAST,IE,I1,LRECTMT
      DOUBLE COMPLEX TMAT0(LMMAXD,LMMAXD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..................................................................
C
      WRITE(6,FMT=99001)
      WRITE(6,FMT=99002)
C
      PI = 4.0D0*ATAN(1.0D0) 
      LRECABMAD = WLENGTH*2*LMPOTD*LMPOTD + WLENGTH*2*LMPOTD
      OPEN (69,ACCESS='direct',RECL=LRECABMAD,FILE='abvmad.unformatted',
     +     FORM='unformatted')
C
      LMMAX = (LMAX+1)*(LMAX+1)
C
      IF (ICC.NE.0) THEN
         DO IQ1=1,NAEZD
            DO LM=1,LMPOTD
               VINTERS(LM,IQ1) = 0D0
            END DO
         END DO
      END IF
C
C ************************************** loop over all types in unit cell

!     DO IT1 = 1,NATYP         ! comment out bauer 2/7/2012

      DO IQ1 = 1,NAEZ            ! added bauer 2/7/2012
      NOQVAL=NOQ(IQ1)            ! added bauer 2/7/2012
      IF (NOQVAL<1) NOQVAL=1     ! added bauer 2/7/2012
      DO IO1 = 1,NOQVAL          ! added bauer 2/7/2012
         IT1 = KAOEZ(IO1,IQ1)    ! added bauer 2/7/2012

C ====================================== take a site occupied by atom IT1
!        IQ1 = IQAT(1,IT1) ! comment out bauer 2/7/2012
C
         IF (IT1/=-1) THEN                ! added bauer 2/7/2012
          IF ( KSHAPE.NE.0 ) THEN
              IRS1 = IRCUT(IPAN(IT1),IT1)
          ELSE
              IRS1 = IRWS(IT1)
          END IF
         END IF                           ! added bauer 2/7/2012
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,LMAX
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO M = -L,L
               LM = L*L + L + M + 1
C     
               AC = 0.0D0
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               IF ( NAEZ.EQ.1 ) THEN
                  IREC = IQ1 + NAEZ*(IQ1-1) 
                  READ(69,REC=IREC) AVMAD,BVMAD
C ============================== loop over all occupants of site IQ2=IQ1
                  DO IO2 = 1,NOQ(IQ1)
                     IT2 = KAOEZ(IO2,IQ1)
C
C---> lm = 1 component disappears if there is only one host atom
C---> take moments of sphere
C
                     DO LM2 = 2,LMMAX
                        AC = AC + AVMAD(LM,LM2)*CMOM(LM2,IT2)
     &                       *CONC(IT2)
                     END DO
C
C---> add contribution of interstial in case of shapes
C
                     IF ( KSHAPE.NE.0 ) THEN
                        DO LM2 = 2,LMMAX
                           AC = AC + AVMAD(LM,LM2)*CMINST(LM2,IT2)
     &                          *CONC(IT2)
                        END DO
                     END IF
                  END DO
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               ELSE
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ loop over all sites
                  DO IQ2 = 1,NAEZ
                     IREC = IQ2 + NAEZ*(IQ1-1) 
                     READ(69,REC=IREC) AVMAD,BVMAD
C =================================== loop over all occupants of site IQ2
                     DO IO2 = 1,NOQ(IQ2)
C
                        IT2 = KAOEZ(IO2,IQ2)
                        AC = AC + BVMAD(LM)*ZAT(IT2)*CONC(IT2)
C
C---> take moments of sphere
C
                        DO LM2 = 1,LMMAX
                           AC = AC + AVMAD(LM,LM2)*CMOM(LM2,IT2)
     &                          *CONC(IT2)
                        END DO
C
C---> add contribution of interstial in case of shapes
C
                        IF ( KSHAPE.NE.0 ) THEN
                           DO LM2 = 1,LMMAX
                              AC = AC + AVMAD(LM,LM2)*
     &                             CMINST(LM2,IT2)*CONC(IT2)
                           END DO
                        END IF
                     END DO    ! IO2 = 1, NOQ(IQ2)
C ======================================================================
                  END DO       ! IQ2 = 1, NAEZ
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               END IF          ! NAEZ.GT.1
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
               IF ( LM.EQ.1 ) 
     &              WRITE (6,FMT=99003) IT1,(CATOM(IT1)-ZAT(IT1)),
     &                    (AC/SQRT(4.D0*PI))
C
C---> add to v the intercell-potential
C
C ================================================================= SPIN
               DO ISPIN = 1,NSPIN
C
C---> determine the right potential number
C
                  IPOT = NSPIN*(IT1-1) + ISPIN
C
C---> in the case of l=0 : r(1)**l is not defined
C
                  IF (IT1/=-1) THEN                ! added bauer 2/7/2012
                    IF ( L.EQ.0 ) V(1,1,IPOT) = V(1,1,IPOT) + AC
                    DO I = 2,IRS1
                       V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IT1))**L*AC
                    END DO
                  END IF
               END DO                              ! added bauer 2/7/2012
C ================================================================= SPIN
               IF (ICC.NE.0 .or. OPT('KKRFLEX ')) THEN
                  LM = L*L + L + M + 1
                  write(*,*) 'ac',iq1,lm,ac
                  VINTERS(LM,IQ1) = AC
               END IF
C 
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      END DO
      END DO
C *********************************************************************
      CLOSE(69)
C *********************************************************************
      WRITE(*,*) 'ICC in VMADELBLK',ICC
      WRITE(6,'(25X,30(1H-),/)')
      WRITE(6,'(79(1H=))')
C
      IF (ICC.EQ.0 .and. OPT('KKRFLEX ')==.false.) RETURN
C *********************************************************************
C
C Now Prepare output for Impurity calculation 
C
C *********************************************************************
      OPEN (91,FILE='intercell_ref',STATUS='unknown',FORM='formatted')
      WRITE(6,*) 
      WRITE(6,*) '                     ',
     &     'Writing intercell potential for impurity'
      WRITE(6,'(/,20X,55(1H-))')
      WRITE(6,99004) HOSTIMP(0),LMMAX
      WRITE(6,'(20X,55(1H-),/,35X,"  i host lm  Vint")') 
      DO I=1,HOSTIMP(0)
         WRITE(6,*)
         LM = 1
         WRITE(6,'(35X,I4,I4,I3,1X,F10.6)') I, HOSTIMP(I),
     &        LM,VINTERS(LM,HOSTIMP(I))
         DO LM=2,9
            WRITE (6,'(43X,I3,1X,F10.6)') LM,VINTERS(LM,HOSTIMP(I))
         END DO
         WRITE(6,'(20X,55(1H-))')
      END DO
      WRITE(6,'(79(1H=),/)')
C         
      WRITE(91,99005) HOSTIMP(0),LMMAX
      DO I=1,HOSTIMP(0)
         WRITE(91,99006) (VINTERS(LM,HOSTIMP(I)),LM=1,LMMAX)
      END DO
      CLOSE(91)

!       write(*,*) 'KKRFLEX WRITEOUT'
!       writE(*,*) OPT('KKRFLEX ')
!       IF ( OPT('KKRFLEX ') ) THEN
!         OPEN (6699,FILE='kkrflex_tmat',STATUS='unknown')
!         write(6699,*) '#',NATOMIMP,NSPIN,IELAST,LMMAXD
!         OPEN (69,ACCESS='direct',RECL=16*LMMAXD*LMMAXD,
!      &             FILE='tmat',  FORM='unformatted')
! 
!         DO IATOM = 1,NATOMIMP
!           I1=ATOMIMP(IATOM)
!           DO ISPIN=1,NSPIN
!             DO IE=1,IELAST
!               IF (I1<=NATYP) THEN 
!                   IREC = IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) 
!                   write(*,*) irec
! !                 READ (69,REC=IREC) TMAT
! !                 IREC = IE + IELAST* (ISPIN-1) + IELAST*NSPIN* (IATOM-1) 
!                   READ (69,REC=IREC) TMAT0
!               ELSE
!                  TMAT0=(0.0D0,0.0D0)
!               END IF 
! !                   write(6699,*) IREC,IE,ISPIN,IATOM,I1
!                   write(6699,'(4I,50000E)') IATOM,ISPIN,IE,0,TMAT0
!             END DO
!           END DO !IE=1,IELAST
!         END DO !ISPIN=1,NSPIN
!         CLOSE(69)
!         CLOSE(6699)
! 
!         OPEN (91,FILE='kkrflex_intercell_ref',STATUS='unknown')
!         WRITE(91,*) '# Intercell potential of each atom'
!         WRITE(91,*) '# '
!         WRITE(91,*) '# NATOMIMP',NATOMIMP
!         WRITE(91,*) '# LMMAX',LMMAX
!         WRITE(91,*) '# KSHAPE',KSHAPE
!         WRITE(91,*) '# NATOMIMP, LMMAX, ALAT VBC(1), VBC(2)'
! !         WRITE(91,*) HOSTIMP(0),LMMAX
!         WRITE(91,'(2I,10F)') NATOMIMP,LMMAX,ALAT,VBC(1),VBC(2)
!         DO IATOM = 1,NATOMIMP
!           I=ATOMIMP(IATOM)
! !         DO I=1,HOSTIMP(0)
!          write(*,*) 'ac2',I,HOSTIMP(I),LM,
!      +                (VINTERS(LM,HOSTIMP(I)),LM=1,LMMAX)
!           WRITE(91,'(5000F)') (VINTERS(LM,HOSTIMP(I)),LM=1,LMMAX)
!         END DO
!         CLOSE(91)
! 
!         OPEN (91,FILE='kkrflex_intercell_cmoms',STATUS='unknown')
!         WRITE(91,*) '# Charge moments of each atom in the unit cell'
!         WRITE(91,*) '# Values given are CMOM + CMINST'
!         WRITE(91,*) '# Firts colums is the core charge other'
!         WRITE(91,*) '# other colums the charge moments'
!         WRITE(91,*) '# NATOMIMP',NATOMIMP
!         WRITE(91,*) '# LMMAX',LMMAX
!         WRITE(91,*) '# KSHAPE',KSHAPE
! 
!         DO IATOM = 1,NATOMIMP
!           I=ATOMIMP(IATOM)
! !         DO IATOM=1,HOSTIMP(0)
! !           I=HOSTIMP(IATOM)
!           WRITE(*,*) 'NOQ',I,NOQ(I)
!           IF (NOQ(I)/=1 .and. NOQ(I)/=0) 
!      +      stop '[vmadelblk] VIRATOMS: NOQ/=1'
!           IF (NOQ(I)==0) then
!             WRITE(91,'(5000F)') 0.0D0,(0.0D0,LM=1,LMMAX)
!           ELSE
!             IF ( KSHAPE.NE.0 ) THEN
!               WRITE(91,'(5000F)') ZAT(KAOEZ(1,I)), 
!      +                          ((CMOM(LM,KAOEZ(1,I)) 
!      +             + CMINST(LM,KAOEZ(1,I)))*CONC(KAOEZ(1,I)),LM=1,LMMAX)
!             ELSE
!               WRITE(91,'(5000F)') ZAT(KAOEZ(1,I)), 
!      +            (CMOM(LM,KAOEZ(1,I))*CONC(KAOEZ(1,I)),LM=1,LMMAX)
!             END IF 
!           END IF
!         END DO
!         CLOSE(91)
! 
! 
!       END IF

      RETURN
C
99001 FORMAT (79(1H=),/,18X,' MADELUNG POTENTIALS ',
     &        '(spherically averaged) ')
99002 FORMAT (/,25X,' ATOM ','  Delta_Q  ','     VMAD',/,25X,30(1H-))
99003 FORMAT (25X,I4,2X,F10.6,1X,F12.6)
99004 FORMAT (22X,I4,' host atoms, LMPOT = ',I2,' output up to LM = 9')
99005 FORMAT (3I6)
99006 FORMAT (4D20.10)
      END
