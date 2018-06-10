SUBROUTINE gijdmat(tauq,tsst,mssq,dmat,dtil,cfctorinv,iprint,  &
        ie,it,krel,lmmaxd)
! **********************************************************************
! * Subroutine to get the projection matrices                          *
! *                                                                    *
! *                ii     i       i    (-1)                            *
! *    D  = [ 1 + G   * (t    -  t  ) ]                                *
! *     a          CPA    CPA     a                                    *
! *                                                                    *
! *    _            i       i      ii   (-1)                           *
! *    D  = [ 1 + (t    -  t  ) * G    ]                               *
! *     a           CPA     a      CPA                                 *
! *                                                                    *
! * Used is made of the same routine GETDMAT as in the case of TAU     *
! * projection. Note, however, that there the matrices are defined for *
! * example as                                                         *
! *                                                                    *
! *                ij     i       i      (-1)                          *
! *    D  = [ 1 + tau * (m    -  m    ) ]                              *
! *     a          CPA    a       CPA                                  *
! *                                                                    *
! * with m = (t)**(-1)                                                 *
! *                                                                    *
! *                                      v.popescu Oct. 2004           *
! **********************************************************************

      IMPLICIT NONE
!..
!.. Arguments ..
INTEGER IPRINT,LMMAXD
INTEGER IE,IT,KREL
DOUBLE COMPLEX CFCTORINV
DOUBLE COMPLEX TAUQ(LMMAXD,LMMAXD),TSST(LMMAXD,LMMAXD), &
               MSSQ(LMMAXD,LMMAXD)
DOUBLE COMPLEX DMAT(LMMAXD,LMMAXD),DTIL(LMMAXD,LMMAXD)
!..
!.. Locals ..
INTEGER IK,INFO,I1,I2
INTEGER IPVT(LMMAXD)
DOUBLE COMPLEX CONE,CZERO
DOUBLE COMPLEX GLL(LMMAXD,LMMAXD),TSSQ(LMMAXD,LMMAXD), &
               TPG(LMMAXD,LMMAXD),XC(LMMAXD,LMMAXD)
CHARACTER*18 BANNER
!..
!.. Externals
      EXTERNAL CMATSTR,GETDMAT,ZCOPY,ZGEMM,ZGETRF,ZGETRI,ZSCAL
!..
!.. Data
      DATA CONE / (1D0,0D0) /
      DATA CZERO / (0D0,0D0) /
! --> get G(CPA) using the same procedure as for GMATLL in < KLOOPZ >
!           G(CPA) = -MSSQ - MSSQ * TAUQ * MSSQ

DO i2 = 1,lmmaxd
  DO i1 = 1,lmmaxd
    tpg(i1,i2) = -cone*tauq(i1,i2)
  END DO
END DO

CALL zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,mssq,  &
    lmmaxd,tpg,lmmaxd,czero,xc,lmmaxd)

CALL zcopy(lmmaxd*lmmaxd,mssq,1,gll,1)

CALL zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,-cone,xc,lmmaxd,  &
    mssq,lmmaxd,-cone,gll,lmmaxd)

CALL zscal(lmmaxd*lmmaxd,cfctorinv,gll,1)

! --> invert MSSQ to get TSSQ

DO i2 = 1,lmmaxd
  DO i1 = 1,lmmaxd
    tssq(i1,i2) = mssq(i1,i2) * cfctorinv
  END DO
END DO
CALL zgetrf(lmmaxd,lmmaxd,tssq,lmmaxd,ipvt,info)
CALL zgetri(lmmaxd,tssq,lmmaxd,ipvt,xc,lmmaxd*lmmaxd,info)

! --> get projection matrices from G(CPA)
!     call getdmat with t,c for the local arg. c,t

CALL getdmat(gll,dmat,dtil,xc,lmmaxd,tsst,tssq,lmmaxd)

IF ( iprint <= 1 ) RETURN

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
ik = 2*krel + 1
WRITE(banner,'("DMAT IE=",I3," IT=",I3)') ie,it
CALL cmatstr(banner,18,dmat,lmmaxd,lmmaxd,ik,ik,0,1D-10,6)

WRITE(banner,'("DTIL IE=",I3," IT=",I3)') ie,it
CALL cmatstr(banner,18,dtil,lmmaxd,lmmaxd,ik,ik,0,1D-10,6)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

END SUBROUTINE gijdmat
