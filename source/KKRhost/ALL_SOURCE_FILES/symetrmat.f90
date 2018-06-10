SUBROUTINE symetrmat(nsym,cpref,dsymll,symunitary,matq,  &
        iqs,matsym,lmmaxd)
! **********************************************************************
! *                                                                    *
! *  Symmetrising the t/G matrix (or their inverses):                  *
! *                                                                    *
! *                       nsym                                         *
! *     MATSYM = CPREF *  SUM  [ DLL(i) * MATQ(iqs(i)) * DLL(i)^T ]    *
! *                      i = 1                                         *
! *                                                                    *
! *  IQS - set outside the routine - has either the same value         *
! *        regardless of i (e.g. in case of the single-site matrices)  *
! *        or is taking on the value of i => MATSYM is a cummulative   *
! *        sum over MATQ[1...nsym] (e.g. in case of the BZ-integration *
! *        of G)                                                       *
! *                                                                    *
! * CPREF = 1/NSYM  or 1/VBZ                                           *
! *                                                                    *
! *                                    v.popescu, munich nov. 2004     *
! **********************************************************************
      IMPLICIT NONE
!..
!.. Parameters ..
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER ( CZERO = (0D0,0D0), CONE = (1D0,0D0) )
!..
!.. Arguments ..
      INTEGER LMMAXD,NSYM
      INTEGER IQS(*)
      DOUBLE COMPLEX CPREF
      LOGICAL SYMUNITARY(*)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,*),MATSYM(LMMAXD,LMMAXD), &
                     MATQ(LMMAXD,LMMAXD,*)
!..
!.. Locals ..
      INTEGER ISYM,L1,L2
      CHARACTER*1 CNT
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD),TS(LMMAXD,LMMAXD)
!..
      EXTERNAL ZCOPY,ZGEMM
!..

CALL zcopy(lmmaxd*lmmaxd,matq(1,1,iqs(1)),1,ts,1)

! ----------------------------------------------------------------------
DO isym = 2,nsym
  
! --> look if the symmetry operation is unitary / anti-unitary
  
  cnt = 'N'
  IF ( .NOT.symunitary(isym) ) cnt = 'T'
  
  CALL zgemm('N',cnt,lmmaxd,lmmaxd,lmmaxd,cone,dsymll(1,1,isym),  &
      lmmaxd,matq(1,1,iqs(isym)),lmmaxd,czero,w1,lmmaxd)
  
  CALL zgemm('N','C',lmmaxd,lmmaxd,lmmaxd,cone,w1,  &
      lmmaxd,dsymll(1,1,isym),lmmaxd,cone,ts,lmmaxd)
END DO
! ----------------------------------------------------------------------

DO l2 = 1,lmmaxd
  DO l1 = 1,lmmaxd
    matsym(l1,l2) = cpref * ts(l1,l2)
  END DO
END DO

END SUBROUTINE symetrmat
