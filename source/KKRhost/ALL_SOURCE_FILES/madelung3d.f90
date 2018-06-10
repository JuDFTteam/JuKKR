SUBROUTINE madelung3d(lpot,yrg,wg,naez,alat,volume0,  &
    bravais,recbv,rbasis,rmax,gmax, naezd,lmxspd,lassld,lpotd,lmpotd,  &
    nmaxd,ishld,nembd,wlength)
! **********************************************************************
! *                                                                    *
! * This subroutine calculates the Madelung potential coefficients     *
! * in the 3D case and stores them in the DA-file abvmad.unformatted   *
! * The record index is simply (IQ1-1)*NAEZ + IQ2 for record (IQ1,IQ2) *
! *                                                                    *
! **********************************************************************

      IMPLICIT NONE
!..
!.. Scalar Arguments ..
      INTEGER LPOT,NAEZ,WLENGTH
      INTEGER NAEZD,LMXSPD,LASSLD,LPOTD,LMPOTD,NMAXD,ISHLD,NEMBD
      DOUBLE PRECISION ALAT,VOLUME0,RMAX,GMAX
!..
!.. Array Arguments ..
      DOUBLE PRECISION YRG(LASSLD,0:LASSLD,0:LASSLD),WG(LASSLD)
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION RBASIS(3,NAEZD+NEMBD)
!..
!.. Local Scalars ..
      INTEGER IEND,IPRINT,IQ1,IQ2,NCLEBD
      INTEGER NGMAX,NRMAX,NSHLG,NSHLR
      INTEGER LRECABMAD,IREC
!..
!.. Local Arrays ..
!.. Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
      DOUBLE PRECISION CLEB(LMXSPD*LMPOTD)
      DOUBLE PRECISION MADELSMAT(LMXSPD,NAEZD,NAEZD)
      DOUBLE PRECISION SMAT1(6,6),SMAT2(6,6)
      DOUBLE PRECISION GN(3,NMAXD),RM(3,NMAXD)
      INTEGER ICLEB(LMXSPD*LMPOTD,3)
      INTEGER NSG(ISHLD),NSR(ISHLD)
!..
!.. External subroutines
      EXTERNAL LATTICE3D,STRMAT,MADELGAUNT,MADELCOEF
! ......................................................................
iprint = 0
nclebd = lmxspd*lmpotd

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,'(79(1H=))')
WRITE (1337,'(18X,A)') 'MADELUNG3D: setting bulk Madelung coefficients'
WRITE (1337,'(79(1H=))')
WRITE (1337,*)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ======================================================================
CALL lattice3d(alat,bravais,recbv,ngmax,nrmax,nshlg,nshlr,nsg,nsr,  &
    gn,rm,rmax,gmax,iprint,nmaxd,ishld)

CALL strmat(alat,lpot,naez,ngmax,nrmax,nsg,nsr,nshlg,nshlr,gn,rm,  &
    rbasis,madelsmat,volume0,iprint,lassld,lmxspd,naezd)
! ======================================================================

lrecabmad = wlength*2*lmpotd*lmpotd + wlength*2*lmpotd
OPEN (69,ACCESS='direct',RECL=lrecabmad,FILE='abvmad.unformatted',  &
    FORM='unformatted')

! --> calculate the gaunt coefficients

CALL madelgaunt(lpot,yrg,wg,cleb,icleb,iend,lassld,nclebd)

! --> calculate the madelung coefficients to be used for VMAD
!     call MADELCOEF with first arg. .FALSE. = 3D case

DO iq1 = 1,naez
  DO iq2 = 1,naez
    CALL madelcoef(.false.,lpot,avmad,bvmad,madelsmat(1,iq1,iq2)  &
        ,cleb,icleb,iend,lpotd,lmpotd,lmxspd,nclebd)
    
    irec = iq2 + naez*(iq1-1)
    WRITE (69,REC=irec) avmad,bvmad
!-----------------------------------------------------------------------
    IF ( (iq1 <= 6) .AND. (iq2 <= 6) ) THEN
      smat1(iq1,iq2) = avmad(1,1)
      smat2(iq1,iq2) = bvmad(1)
    END IF
!-----------------------------------------------------------------------
  END DO
END DO
CLOSE (69)

IF ( iprint < 1 ) RETURN
! ======================================================================

CALL madel3out(iprint,naez,lrecabmad,smat1,smat2,lmpotd)

END SUBROUTINE madelung3d
