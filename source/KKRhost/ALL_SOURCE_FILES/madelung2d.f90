SUBROUTINE madelung2d(lpot,yrg,wg,naez,alat,vol,  &
        bravais,recbv,rbasis,rmax,gmax,  &
        nlbasis,nleft,zperleft,tleft,  &
        nrbasis,nright,zperight,tright,  &
        lmxspd,lassld,lpotd,lmpotd,  &
        nmaxd,ishld,nembd1,wlength)
! **********************************************************************
! *                                                                    *
! * This subroutine calculates the Madelung potential coefficients     *
! * in the 2D case and stores them in the DA-file abmad.unformatted    *
! * For each layer in the slab, the summation is split into three      *
! * parts (see also VINTERFACE):                                       *
! * within the slab, over the NLEFT*NLBASIS left host sites and over   *
! * the NRIGHT*NRBASIS right host sites, the last two steps only in    *
! * case of decimation run                                             *
! *                                                                    *
! * all positions must be scaled with ALAT to get them correct         *
! * (done in EWALD2D)                                                  *
! *                                                                    *
! * The record index is:                                               *
! *   (IQ1-1)*NAEZ + IQ2                 for (IQ1,IQ2) within the slab *
! *   NAEZ*NAEZ + (IQ1-1)*NLEFT*NLBASIS  for (IQ1,(IL,IBL)), IQ1 in    *
! *                  + (IL-1)*NLEFT+IBL  slab, (IL,IBL) in the left    *
! *   NAEZ*NAEZ + NAEZ*NLEFT*NLBASIS                                   *
! *             + (IQ1-1)*NRIGHT*NRBASIS for (IQ1,(IR,IBR)), IQ1 in    *
! *             + (IR-1)*NRIGHT+IBR      slab, (IR,IBR) in the right   *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
!..
!.. Scalar Arguments ..
      INTEGER LPOT,NAEZ,WLENGTH
      INTEGER NLBASIS,NLEFT,NRBASIS,NRIGHT
      INTEGER LASSLD,LPOTD,LMPOTD,LMXSPD,NMAXD,ISHLD,NEMBD1
      DOUBLE PRECISION ALAT,VOL,RMAX,GMAX
!..
!.. Array Arguments ..
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION YRG(LASSLD,0:LASSLD,0:LASSLD),WG(LASSLD)
      DOUBLE PRECISION RBASIS(3,*)
      DOUBLE PRECISION ZPERIGHT(3),ZPERLEFT(3)
      DOUBLE PRECISION TLEFT(3,NEMBD1),TRIGHT(3,NEMBD1)
!..
!.. Local Scalars ..
      INTEGER IQ1,IQ2,IEND,NCLEBD,IPRINT
      INTEGER I,IB,IH,ILEFT,IRIGHT
      INTEGER LRECAMAD,IREC,NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL
      INTEGER NGMAX,NRMAX,NSHLG,NSHLR
      LOGICAL OPT
!..
!.. Local Arrays ..
!.. Attention: LMXSPD*LMPOTD appears as NCLEB1 in other routines
      DOUBLE PRECISION CLEB(LMXSPD*LMPOTD)
      DOUBLE PRECISION BM(LMPOTD),VEC2(3),SUM(LMXSPD)
      DOUBLE PRECISION GN2(2,NMAXD),RM2(2,NMAXD)
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD)
      INTEGER NSG(ISHLD),NSR(ISHLD)
      INTEGER ICLEB(LMXSPD*LMPOTD,3)
!..
!.. External Functions/Subroutines
! ......................................................................
iprint = 0
nclebd = lmxspd*lmpotd

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,'(79(1H=))')
WRITE (1337,'(18X,A)') 'MADELUNG2D: setting 2D Madelung coefficients'
WRITE (1337,'(79(1H=))')
WRITE (1337,*)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ======================================================================
CALL lattice2d(alat,bravais,recbv,ngmax,nrmax,nshlg,nshlr,nsg,nsr,  &
    gn2,rm2,rmax,gmax,iprint,nmaxd,ishld)
! ======================================================================

lrecamad = wlength*2*lmpotd*lmpotd
OPEN (69,ACCESS='direct',RECL=lrecamad,FILE='avmad.unformatted',  &
    FORM='unformatted')

! --> calculate the gaunt coefs

CALL madelgaunt(lpot,yrg,wg,cleb,icleb,iend,lassld,nclebd)

! --> calculate the madelung coefficients to be used for VMAD

! **********************************************************************
! ********************************************** loop over atoms in slab
DO iq1 = 1,naez
  
!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                     1.  Summation in all layers in the slab
!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! ++++++++++++++++++++++++++++++++ loop over all other sites in the slab
  DO iq2 = 1,naez
    
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    IF ( iq1 == 1 .AND. iq2 == 1 ) THEN
      WRITE (1337,'(5X,2A,/)')  &
          '< EWALD2D > : calculating 2D-lattice sums ', 'inside the slab'
      IF ( iprint >= 2 ) WRITE(1337,99001)
    END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    
    
! make ewald sumation in plane and inverse space
! sum if rz<>0 (out of plane)
    
!             WRITE(99,*) 'Layer pair:',IQ1,IQ2
    CALL ewald2d(lpot,alat,rbasis(1,iq1),rbasis(1,iq2),iq1,iq2,  &
        rm2,nrmax,nshlr,nsr,gn2, ngmax,nshlg,nsg,sum,vol,lassld,lmxspd)
!             WRITE(99,*) 'SUM: ',SUM
    
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    IF ( iprint >= 2 ) THEN
      WRITE(1337,99002) iq1,iq2,sum(1)
      IF ( iq2 == naez .AND. iq1 /= naez ) WRITE(1337,'(20X,20(1H-))')
    END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    
    CALL madelcoef(.true.,lpot,avmad,bm,sum,cleb,icleb,iend,  &
        lpotd,lmpotd,lmxspd,nclebd)
    
    irec = iq2 + naez*(iq1-1)
    WRITE(69,REC=irec) avmad
  END DO
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END DO
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
IF ( iprint >= 2 ) WRITE(1337,'(18X,22(1H-),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
! ********************************************** loop over atoms in slab

! ######################################################################
IF ( opt('DECIMATE') ) THEN
  
  nleftoff = naez * naez                        ! record offsets
  nrightoff = nleftoff + naez * nleft * nlbasis ! left and right
  nleftall = nleft * nlbasis
  nrightall = nright * nrbasis
  
! ********************************************** loop over atoms in slab
  DO iq1 = 1,naez
!     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                     2.  Summation in the LEFT bulk side
!     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ileft = 0
! ++++++++++++++++++++++++++++++++ loop over all sites in the left host
    DO ih = 1,nleft
      DO ib = 1,nlbasis
        DO i = 1,3
          vec2(i) = (tleft(i,ib)+(ih-1)*zperleft(i))
        END DO
        ileft = ileft + 1
        
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        IF ( iq1 == 1 .AND. ileft == 1 ) THEN
          WRITE (1337,'(5X,2A,/)')  &
              '< EWALD2D > : calculating 2D-lattice sums ', 'slab - left host'
          IF ( iprint >= 2 ) WRITE(1337,99001)
        END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        
        
!-->  make ewald sumation for m= 0 l<5 rz=0 (in plane) and
!     Inverse space sum if rz<>0 (out of plane)
        
        CALL ewald2d(lpot,alat,rbasis(1,iq1),vec2,iq1,ih,  &
            rm2,nrmax,nshlr,nsr,gn2, ngmax,nshlg,nsg,sum,vol,lassld,lmxspd)
        
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        IF ( iprint >= 2 ) THEN
          WRITE(1337,99002) iq1,ileft,sum(1)
          IF ( ileft == nleftall .AND. iq1 /= naez )  &
              WRITE(1337,'(20X,20(1H-))')
        END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        
        CALL madelcoef(.true.,lpot,avmad,bm,sum,cleb,icleb,  &
            iend,lpotd,lmpotd,lmxspd,nclebd)
        
        irec = ileft + nleftall*(iq1-1) + nleftoff
        WRITE(69,REC=irec) avmad
      END DO           ! ib loop in left host basis
    END DO              ! ih loop in layers to get convergence
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( ileft /= nleftall) THEN
      WRITE(6,*) ' < MADELUNG2D > : index error ', 'ILEFT <> NLEFT*NLBASIS'
      STOP
    END IF
  END DO                 ! ILAY1 loop
! **********************************************************************
  
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  IF ( iprint >= 2 ) WRITE(1337,'(18X,22(1H-),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  
! ********************************************** loop over atoms in slab
  DO iq1 = 1,naez
!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                     3.  Summation in the RIGHT bulk side
!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! ++++++++++++++++++++++++++++++++ loop over all sites in the right host
    iright = 0
    DO ih = 1,nright
      DO ib = 1,nrbasis
        DO i = 1,3
          vec2(i) = (tright(i,ib)+(ih-1)*zperight(i))
        END DO
        iright = iright + 1
        
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        IF ( iq1 == 1 .AND. iright == 1 ) THEN
          WRITE (1337,'(5X,2A,/)')  &
              '< EWALD2D > : calculating 2D-lattice sums ', 'slab - right host'
          IF ( iprint >= 2 ) WRITE(1337,99001)
        END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        
!-->  make ewald sumation (in plane) and
!     Inverse space sum if rz<>0 (out of plane)
        
        CALL ewald2d(lpot,alat,rbasis(1,iq1),vec2,iq1,ih,  &
            rm2,nrmax,nshlr,nsr,gn2, ngmax,nshlg,nsg,sum,vol,lassld,lmxspd)
        
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        IF ( iprint >= 2 ) THEN
          WRITE(1337,99002) iq1,iright,sum(1)
          IF ( iright == nrightall .AND. iq1 /= naez )  &
              WRITE(1337,'(20X,20(1H-))')
        END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        
        CALL madelcoef(.true.,lpot,avmad,bm,sum,cleb,icleb,  &
            iend,lpotd,lmpotd,lmxspd,nclebd)
        
        irec = iright + nrightall*(iq1-1) + nrightoff
        WRITE(69,REC=irec) avmad
      END DO           ! ib loop in right host basis
    END DO              ! ih loop in layers to get convergence
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( iright /= nrightall ) THEN
      WRITE(6,*) ' < MADELUNG2D > : index error ', 'IRIGHT <> NRIGHT*NRBASIS'
      STOP
    END IF
  END DO                 ! ILAY1 loop
! **********************************************************************
  
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  IF ( iprint >= 2 ) WRITE(1337,'(18X,22(1H-),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  
END IF
! ######################################################################
CLOSE(69)

IF ( iprint < 1 ) RETURN
! ======================================================================

CALL madel2out(iprint,naez,lrecamad,lmpotd,  &
    nleftoff,nrightoff,nleftall,nrightall)

99001 FORMAT (8X,'2D Lattice sum (LMXSP = 1)',/,  &
    18X,'  IQ1  IQ2  SUM',/,18X,23(1H-))
99002 FORMAT (18X,2I5,d12.4)
END SUBROUTINE madelung2d
