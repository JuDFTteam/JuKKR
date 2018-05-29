      SUBROUTINE MADELUNG2D(LPOT,YRG,WG,NAEZ,ALAT,VOL,
     &                      BRAVAIS,RECBV,RBASIS,RMAX,GMAX,
     &                      NLBASIS,NLEFT,ZPERLEFT,TLEFT,
     &                      NRBASIS,NRIGHT,ZPERIGHT,TRIGHT,
     &                      LMXSPD,LASSLD,LPOTD,LMPOTD,
     &                      NMAXD,ISHLD,NEMBD1,WLENGTH)
C **********************************************************************
C *                                                                    *
C * This subroutine calculates the Madelung potential coefficients     *
C * in the 2D case and stores them in the DA-file abmad.unformatted    *
C * For each layer in the slab, the summation is split into three      *
C * parts (see also VINTERFACE):                                       *
C * within the slab, over the NLEFT*NLBASIS left host sites and over   *
C * the NRIGHT*NRBASIS right host sites, the last two steps only in    *
C * case of decimation run                                             *
C *                                                                    *
C * all positions must be scaled with ALAT to get them correct         *
C * (done in EWALD2D)                                                  *
C *                                                                    *
C * The record index is:                                               *
C *   (IQ1-1)*NAEZ + IQ2                 for (IQ1,IQ2) within the slab *
C *   NAEZ*NAEZ + (IQ1-1)*NLEFT*NLBASIS  for (IQ1,(IL,IBL)), IQ1 in    *
C *                  + (IL-1)*NLEFT+IBL  slab, (IL,IBL) in the left    *
C *   NAEZ*NAEZ + NAEZ*NLEFT*NLBASIS                                   *
C *             + (IQ1-1)*NRIGHT*NRBASIS for (IQ1,(IR,IBR)), IQ1 in    *
C *             + (IR-1)*NRIGHT+IBR      slab, (IR,IBR) in the right   *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NAEZ,WLENGTH
      INTEGER NLBASIS,NLEFT,NRBASIS,NRIGHT
      INTEGER LASSLD,LPOTD,LMPOTD,LMXSPD,NMAXD,ISHLD,NEMBD1
      DOUBLE PRECISION ALAT,VOL,RMAX,GMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION YRG(LASSLD,0:LASSLD,0:LASSLD),WG(LASSLD)
      DOUBLE PRECISION RBASIS(3,*)
      DOUBLE PRECISION ZPERIGHT(3),ZPERLEFT(3)
      DOUBLE PRECISION TLEFT(3,NEMBD1),TRIGHT(3,NEMBD1)
C     ..
C     .. Local Scalars ..
      INTEGER IQ1,IQ2,IEND,NCLEBD,IPRINT
      INTEGER I,IB,IH,ILEFT,IRIGHT
      INTEGER LRECAMAD,IREC,NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL
      INTEGER NGMAX,NRMAX,NSHLG,NSHLR
      LOGICAL OPT
C     ..
C     .. Local Arrays ..
C     .. Attention: LMXSPD*LMPOTD appears as NCLEB1 in other routines
      DOUBLE PRECISION CLEB(LMXSPD*LMPOTD)
      DOUBLE PRECISION BM(LMPOTD),VEC2(3),SUM(LMXSPD)
      DOUBLE PRECISION GN2(2,NMAXD),RM2(2,NMAXD)
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD)
      INTEGER NSG(ISHLD),NSR(ISHLD)
      INTEGER ICLEB(LMXSPD*LMPOTD,3)
C     ..
C     .. External Functions/Subroutines
      EXTERNAL EWALD2D,OPT,MADELGAUNT,MADELCOEF
C ......................................................................
      IPRINT = 0
      NCLEBD = LMXSPD*LMPOTD
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(79(1H=))')
      WRITE (1337,'(18X,A)') 
     &                  'MADELUNG2D: setting 2D Madelung coefficients'
      WRITE (1337,'(79(1H=))')
      WRITE (1337,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ======================================================================
      CALL LATTICE2D(ALAT,BRAVAIS,RECBV,NGMAX,NRMAX,NSHLG,NSHLR,NSG,NSR,
     +               GN2,RM2,RMAX,GMAX,IPRINT,NMAXD,ISHLD)
C ======================================================================
C
      LRECAMAD = WLENGTH*2*LMPOTD*LMPOTD
      OPEN (69,ACCESS='direct',RECL=LRECAMAD,FILE='avmad.unformatted',
     +     FORM='unformatted')
C
C --> calculate the gaunt coefs
C
      CALL MADELGAUNT(LPOT,YRG,WG,CLEB,ICLEB,IEND,LASSLD,NCLEBD)
C
C --> calculate the madelung coefficients to be used for VMAD
C
C **********************************************************************
C ********************************************** loop over atoms in slab
      DO IQ1 = 1,NAEZ
C
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                     1.  Summation in all layers in the slab
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C ++++++++++++++++++++++++++++++++ loop over all other sites in the slab
         DO IQ2 = 1,NAEZ
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
            IF ( IQ1.EQ.1 .AND. IQ2.EQ.1 ) THEN
               WRITE (1337,'(5X,2A,/)')
     &              '< EWALD2D > : calculating 2D-lattice sums ',
     &              'inside the slab'
               IF ( IPRINT.GE.2 ) WRITE(1337,99001)
            END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C
C make ewald sumation in plane and inverse space 
C sum if rz<>0 (out of plane)
C
!             WRITE(99,*) 'Layer pair:',IQ1,IQ2
            CALL EWALD2D(LPOT,ALAT,RBASIS(1,IQ1),RBASIS(1,IQ2),IQ1,IQ2,
     &                   RM2,NRMAX,NSHLR,NSR,GN2,
     &                   NGMAX,NSHLG,NSG,SUM,VOL,LASSLD,LMXSPD)
!             WRITE(99,*) 'SUM: ',SUM
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
            IF ( IPRINT.GE.2 ) THEN
               WRITE(1337,99002) IQ1,IQ2,SUM(1)
               IF ( IQ2.EQ.NAEZ .AND. IQ1.NE.NAEZ ) 
     &              WRITE(1337,'(20X,20(1H-))')
            END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
            CALL MADELCOEF(.TRUE.,LPOT,AVMAD,BM,SUM,CLEB,ICLEB,IEND,
     &                     LPOTD,LMPOTD,LMXSPD,NCLEBD)
C
            IREC = IQ2 + NAEZ*(IQ1-1) 
            WRITE(69,REC=IREC) AVMAD
         END DO                 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END DO
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      IF ( IPRINT.GE.2 ) WRITE(1337,'(18X,22(1H-),/)')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C ********************************************** loop over atoms in slab

C ######################################################################
      IF ( OPT('DECIMATE') ) THEN

         NLEFTOFF = NAEZ * NAEZ                        ! record offsets
         NRIGHTOFF = NLEFTOFF + NAEZ * NLEFT * NLBASIS ! left and right
         NLEFTALL = NLEFT * NLBASIS
         NRIGHTALL = NRIGHT * NRBASIS

C ********************************************** loop over atoms in slab
         DO IQ1 = 1,NAEZ
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                     2.  Summation in the LEFT bulk side
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            ILEFT = 0
C ++++++++++++++++++++++++++++++++ loop over all sites in the left host
            DO IH = 1,NLEFT
               DO IB = 1,NLBASIS
                  DO I = 1,3
                     VEC2(I) = (TLEFT(I,IB)+(IH-1)*ZPERLEFT(I))
                  END DO
                  ILEFT = ILEFT + 1
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
                  IF ( IQ1.EQ.1 .AND. ILEFT.EQ.1 ) THEN
                     WRITE (1337,'(5X,2A,/)')
     &                    '< EWALD2D > : calculating 2D-lattice sums ',
     &                    'slab - left host'
                     IF ( IPRINT.GE.2 ) WRITE(1337,99001)
                  END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C
C-->  make ewald sumation for m= 0 l<5 rz=0 (in plane) and
C     Inverse space sum if rz<>0 (out of plane)
C
                  CALL EWALD2D(LPOT,ALAT,RBASIS(1,IQ1),VEC2,IQ1,IH,
     &                         RM2,NRMAX,NSHLR,NSR,GN2,
     &                         NGMAX,NSHLG,NSG,SUM,VOL,LASSLD,LMXSPD)
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
                  IF ( IPRINT.GE.2 ) THEN
                     WRITE(1337,99002) IQ1,ILEFT,SUM(1)
                     IF ( ILEFT.EQ.NLEFTALL .AND. IQ1.NE.NAEZ ) 
     &                    WRITE(1337,'(20X,20(1H-))')
                  END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
                  CALL MADELCOEF(.TRUE.,LPOT,AVMAD,BM,SUM,CLEB,ICLEB,
     &                           IEND,LPOTD,LMPOTD,LMXSPD,NCLEBD)
C
                  IREC = ILEFT + NLEFTALL*(IQ1-1) + NLEFTOFF
                  WRITE(69,REC=IREC) AVMAD
               END DO           ! ib loop in left host basis
            END DO              ! ih loop in layers to get convergence
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF ( ILEFT.NE.NLEFTALL) THEN
               WRITE(6,*) ' < MADELUNG2D > : index error ',
     &           'ILEFT <> NLEFT*NLBASIS'
               STOP
            END IF
         END DO                 ! ILAY1 loop
C **********************************************************************
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      IF ( IPRINT.GE.2 ) WRITE(1337,'(18X,22(1H-),/)')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ********************************************** loop over atoms in slab
         DO IQ1 = 1,NAEZ
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                     3.  Summation in the RIGHT bulk side
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C ++++++++++++++++++++++++++++++++ loop over all sites in the right host
            IRIGHT = 0
            DO IH = 1,NRIGHT
               DO IB = 1,NRBASIS
                  DO I = 1,3
                     VEC2(I) = (TRIGHT(I,IB)+(IH-1)*ZPERIGHT(I))
                  END DO
                  IRIGHT = IRIGHT + 1
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
                  IF ( IQ1.EQ.1 .AND. IRIGHT.EQ.1 ) THEN
                     WRITE (1337,'(5X,2A,/)')
     &                    '< EWALD2D > : calculating 2D-lattice sums ',
     &                    'slab - right host'
                     IF ( IPRINT.GE.2 ) WRITE(1337,99001)
                  END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C-->  make ewald sumation (in plane) and
C     Inverse space sum if rz<>0 (out of plane)
C     
                  CALL EWALD2D(LPOT,ALAT,RBASIS(1,IQ1),VEC2,IQ1,IH,
     &                         RM2,NRMAX,NSHLR,NSR,GN2,
     &                         NGMAX,NSHLG,NSG,SUM,VOL,LASSLD,LMXSPD)
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
                  IF ( IPRINT.GE.2 ) THEN
                     WRITE(1337,99002) IQ1,IRIGHT,SUM(1)
                     IF ( IRIGHT.EQ.NRIGHTALL .AND. IQ1.NE.NAEZ ) 
     &                    WRITE(1337,'(20X,20(1H-))')
                  END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
                  CALL MADELCOEF(.TRUE.,LPOT,AVMAD,BM,SUM,CLEB,ICLEB,
     &                           IEND,LPOTD,LMPOTD,LMXSPD,NCLEBD)
C
                  IREC = IRIGHT + NRIGHTALL*(IQ1-1) + NRIGHTOFF
                  WRITE(69,REC=IREC) AVMAD
               END DO           ! ib loop in right host basis
            END DO              ! ih loop in layers to get convergence
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF ( IRIGHT.NE.NRIGHTALL ) THEN
               WRITE(6,*) ' < MADELUNG2D > : index error ',
     &                'IRIGHT <> NRIGHT*NRBASIS'
               STOP
            END IF
         END DO                 ! ILAY1 loop
C **********************************************************************
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      IF ( IPRINT.GE.2 ) WRITE(1337,'(18X,22(1H-),/)')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      END IF
C ######################################################################
      CLOSE(69)
C
      IF ( IPRINT.LT.1 ) RETURN
C ======================================================================
C
      CALL MADEL2OUT(IPRINT,NAEZ,LRECAMAD,LMPOTD,
     &               NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL)
C
99001 FORMAT (8X,'2D Lattice sum (LMXSP = 1)',/,
     &        18X,'  IQ1  IQ2  SUM',/,18X,23(1H-))
99002 FORMAT (18X,2I5,D12.4)
      END
