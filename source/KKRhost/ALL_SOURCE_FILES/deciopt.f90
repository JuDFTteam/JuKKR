SUBROUTINE deciopt(alat,ins,krel,kvrel,kmrot,nspin,naez,lmmax,  &
    bravais,tk,npol,npnt1,npnt2,npnt3, ez,ielast,kaoez,  &
    lefttinvll,righttinvll,vacflag,nlbasis,nrbasis,  &
    cmomhost,vref,rmtref,nref,refpot,  &
    lmaxd,lmgf0d,lmmaxd,lm2d,nembd1,iemxd,nspind, lmpotd,natypd,irmd,ipand)
! **********************************************************************
! *                                                                    *
! * This routine treats the DECIMATION case setting up the single-site *
! * (Delta t)^(-1) matrices and the charge moments of the host(s).     *
! *                                                                    *
! * This is realised in two ways:                                      *
! *      - either reading in the matrices (and moments - if SCFSTEPS   *
! *        is greater than 1) as written out in a previous (bulk) run  *
! *        -- DECIFILES token points to the files containing the nece- *
! *        ssary information                                           *
! *      - or reading in the self-consistent potential for each host   *
! *        and effectively calculating the matrices; the potential     *
! *        must have the specific format set in < OUTPOTHOST > routine *
! *        and the DECIPOTS token should point to the corresponding    *
! *        potential file(s)                                           *
! *                                                                    *
! * Notes:                                                             *
! *        - DECIFILES token is considered by default and sought first *
! *                          is requiring the same energy mesh for the *
! *                          system as for the host                    *
! *        - DECIPOTS token is not restrictive in this sense           *
! *                         however, it does not calculate charge mo-  *
! *                         ments -- does not work in SCF mode         *
! *                         is not dealing with CPA host so far        *
! *                         is not dealing with NON-SPHERICAL poten-   *
! *                         tials so far                               *
! *                                                                    *
! *                                     v.popescu - munich, Dec 04     *
! *                                                                    *
! **********************************************************************

      IMPLICIT NONE
!..
!.. Scalar arguments
      INTEGER LMMAXD,NEMBD1,IEMXD,NSPIND,LMPOTD,NATYPD,IPAND,IRMD,LMAXD
      INTEGER LM2D,NREF,LMGF0D
      INTEGER INS,KREL,KMROT,NSPIN,NAEZ,LMMAX,NPOL,NPNT1,NPNT2,NPNT3
      INTEGER IELAST,NLBASIS,NRBASIS,KVREL
      DOUBLE PRECISION ALAT,TK
!..
!.. Array arguments
      INTEGER KAOEZ(NATYPD,*),REFPOT(NEMBD1)
      DOUBLE PRECISION BRAVAIS(3,3),CMOMHOST(LMPOTD,*)
      DOUBLE PRECISION VREF(*),RMTREF(*)
      DOUBLE COMPLEX LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD), &
                     RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD)
      DOUBLE COMPLEX EZ(IEMXD)
      LOGICAL VACFLAG(2)
!..
!.. Local scalars
      INTEGER IERROR,IL,IE,ISPIN,NSPINSO           ! ruess: for tmat new solver
      DOUBLE COMPLEX CFCTOR
      CHARACTER*40 FILELEFT,FILERIGHT
      CHARACTER*256 UIO ! NCOLIO=256
!..                                  ! ruess: for NEWSOSOL running option
!.. External Functions ..
      LOGICAL OPT
      EXTERNAL OPT

! ======================================================================
WRITE (1337,'(79(1H=))')
WRITE (1337,'(15X,A,/,79(1H=),/)')  &
    'DECIOPT: reading left/right host decimation files'
il = 1
ierror = 0
CALL ioinput('DECIFILES       ',uio,il,7,ierror)
! :::::::::::::::::::::::::::::::::::::::::::::::: decifiles (tmatrices)
IF ( ierror == 0 ) THEN
  READ (UNIT=uio,FMT='(A40)') fileleft
  CALL ioinput('DECIFILES       ',uio,il+1,7,ierror)
  READ (UNIT=uio,FMT='(A40)') fileright
! ----------------------------------------------------------------------
  
! --> first call to read the header ( IE = 0 )
  
  ie = 0
  CALL decimaread(ez,tk,npnt1,npnt2,npnt3,npol,nspin,  &
      lefttinvll(1,1,1,1,1),righttinvll(1,1,1,1,1),vacflag,  &
      ie,nlbasis,nrbasis,naez,kaoez,kmrot,  &
      ins,nspin,lmmax,ielast,fileleft,fileright, krel,natypd,lmmaxd,nembd1)
  
! --> get the left and right host Delta_t matrices
  
  cfctor = alat/(8.d0*ATAN(1.0D0)) ! = ALAT/(2*PI)
  nspinso = nspin
  IF (opt('NEWSOSOL')) nspinso = 1 ! ruess: only combined l-s index for newsolver
  DO ispin = 1,nspinso
    DO ie = 1,ielast
      CALL decimaread(ez,tk,npnt1,npnt2,npnt3,npol,ispin,  &
          lefttinvll(1,1,1,ispin,ie), righttinvll(1,1,1,ispin,ie),vacflag,  &
          ie,nlbasis,nrbasis,naez,kaoez,kmrot,ins,  &
          nspin,lmmax,ielast,fileleft,fileright, krel,natypd,lmmaxd,nembd1)
      
! --> host matrices have been written out in true units
!     they are used in p.u. units (see kloopz) --> convert them here
      
      CALL zscal(lmmaxd*lmmaxd*nembd1,cfctor, lefttinvll(1,1,1,ispin,ie),1)
      CALL zscal(lmmaxd*lmmaxd*nembd1,cfctor, righttinvll(1,1,1,ispin,ie),1)
    END DO
  END DO
  
! --> get the left and right host charge moments
!     ( not needed in single-iteration mode calculations )
  
!fivos        IF ( SCFSTEPS.GT.1 )
  CALL cmomsread(nlbasis,nrbasis,naez, cmomhost,vacflag,kaoez,  &
      natypd,nembd1,lmpotd)
  CLOSE (37)
  CLOSE (38)
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE
! :::::::::::::::::::::::::::::::::::::::::::::::: decipots (calc tmats)
  ierror = 0
  CALL ioinput('DECIPOTS        ',uio,il,7,ierror)
  IF ( ierror /= 0 ) THEN
    WRITE(6,99010)
    STOP
  END IF
  READ (UNIT=uio,FMT='(A40)') fileleft
  CALL ioinput('DECIPOTS        ',uio,il+1,7,ierror)
  READ (UNIT=uio,FMT='(A40)') fileright
  CALL decitset(alat,bravais,ez,ielast, nlbasis,nrbasis,fileleft,fileright,  &
      ins,kvrel,krel,nspin,kmrot, vref,rmtref,nref,refpot,  &
      lefttinvll,righttinvll,vacflag, nembd1,iemxd,irmd,ipand,  &
      lmaxd,lmgf0d,lmmaxd,lm2d,nspind)
END IF
! ======================================================================

99010 FORMAT (/,6X,'ERROR : Missing decimation files (t-mat or pot)',/,  &
    14X,'Please use one of the tokens DECIFILES/DECIPOTS',/)
END SUBROUTINE deciopt
