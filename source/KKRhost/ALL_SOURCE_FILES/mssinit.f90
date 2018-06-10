!**********************************************************************
SUBROUTINE mssinit(ncpa,icpastart,tsst,msst,mssq,trefll,drotq,  &
    refpot,iqat,itoq,noq,conc,  &
    kmrot,natyp,naez,lmmaxd) ! nrefd was taken out of calling list 1.2.2012

      use mod_mympi, only: myrank, master
      IMPLICIT NONE
      include 'inc.p' ! Included  1.2.2012

!.. Parameters
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO = (0.0D0,0.0D0) )
      PARAMETER (CONE = (1.0D0,0.0D0) )

!.. Dummy arguments
INTEGER KMROT, NATYP, NAEZ, LMMAXD, NCPA, ICPASTART
INTEGER IQAT(NATYPD),ITOQ(NATYPD,NAEZD)
INTEGER REFPOT(NAEZD),NOQ(NAEZD)
DOUBLE PRECISION CONC(NATYPD)
DOUBLE COMPLEX  &
     TSST(LMMAXD,LMMAXD,NATYPD), &
     TREFLL(LMMAXD,LMMAXD,NREFD)
DOUBLE COMPLEX MSST(LMMAXD,LMMAXD,NATYPD)
DOUBLE COMPLEX MSSQ(LMMAXD,LMMAXD,NAEZD)
DOUBLE COMPLEX DROTQ(LMMAXD,LMMAXD,NAEZD)

!.. Local variables
      INTEGER IT,IQ,RF,J,IO,INFO,LM1,LM2,LP,LD,LMP,LMD
      INTEGER IPVT(LMMAXD)
      DOUBLE COMPLEX ZC
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD)
      DOUBLE COMPLEX W2(LMMAXD,LMMAXD)
!..
!.. External Subroutines ..
      LOGICAL TEST,OPT
      EXTERNAL TEST,OPT
! ======================================================================

!   --> set up the Delta_t^-1 matrix (MSST) in the LOCAL frame




DO it = 1,natyp
  
  DO j = 1,lmmaxd
    CALL zcopy(lmmaxd,tsst(1,j,it),1,msst(1,j,it),1)
  END DO
  
  iq = iqat(it)
  rf = refpot(iq)
  
  
  IF (kmrot /= 0) THEN
    CALL rotate(trefll(1,1,rf),'G->L',w1,lmmaxd, drotq(1,1,iq),lmmaxd)
  ELSE
    DO j=1,lmmaxd
      CALL zcopy(lmmaxd,trefll(1,j,rf),1,w1(1,j),1)
    END DO
  END IF
  
  
! ---> determine Delta_t = t(sys) - tmat(ref) = TSST - TREFLL
!      in local frame
  
  DO lm1 = 1,lmmaxd
    DO lm2 = 1,lmmaxd
      
      msst(lm2,lm1,it) = ( msst(lm2,lm1,it) - w1(lm2,lm1) )
      
    END DO
  END DO
  
!  --> inversion
  
  IF ( .NOT. opt('VIRATOMS') ) THEN
    IF (.NOT.test('testgmat')) THEN
      CALL zgetrf(lmmaxd,lmmaxd,msst(1,1,it),lmmaxd,ipvt,info)
!         CALL ZGETRI(LMMAXD,MSST(1,1,IT),LMMAXD,IPVT,W1,
!     &               LMMAXD*LMMAXD,INFO)
      DO lm1=1,lmmaxd
        DO lm2=1,lmmaxd
          IF (lm1 == lm2) THEN
            w2(lm1,lm2)=cone
          ELSE
            w2(lm1,lm2)=czero
          END IF
        END DO
      END DO
      CALL zgetrs('N',lmmaxd,lmmaxd,msst(1,1,it),lmmaxd,ipvt,w2, lmmaxd,info)
      DO lm1=1,lmmaxd
        DO lm2=1,lmmaxd
          msst(lm1,lm2,it)=w2(lm1,lm2)
        END DO
      END DO
    END IF !( TEST('testgmat') ) THEN
  END IF !( OPT('VIRATOMS') ) THEN
  
END DO


! ======================================================================

! ---> determine tmat(sys) - tmat(ref) = MSSQ - TREFLL

!      because TSST is calculated in the LOCAL frame, if KMROT<>0
!      it needs to be rotated prior to set the Delta matrix

! ---> set up the effective (on-site) Delta_t- and Delta_m-matrices
!      using the Average T-matrix Approximation

!   ICPASTART=1:
!       m(IQ) = t(ATA) = SUM(it)  c(it) * t(it)
!       m(IQ) = t(ATA) - t_ref = Delta_t(ATA)
!       m(IQ) = (Delta_t(ATA))^(-1)

!   ICPASTART=2:
!       m(IQ) = (Delta_t(ATA))^(-1) for l = 2
!       m(IQ) = SUM(it) c(it)*m(it) with m(it)=t(it)^(-1)

!       mssq(IQ)  refer to the GLOBAL frame
!       tsst(IT),msst(IT)  refer to the LOCAL  frame


! ----------------------------------------------------------------------

CALL cinit(lmmaxd*lmmaxd*naez,mssq)

DO iq = 1,naez
  
  DO io = 1,noq(iq)
    it = itoq(io,iq)
!         write(*,*) 'test fivos mssinit IO,IQ,IT',IO,IQ,IT ! test fivos
    zc = conc(it)
    DO j = 1,lmmaxd
      CALL zaxpy(lmmaxd,zc,tsst(1,j,it),1,mssq(1,j,iq),1)
    END DO
  END DO
  
! ---> rotate MSSQ from the LOCAL to the GLOBAL frame if necessary
  
  IF ( kmrot /= 0 ) THEN
    CALL rotate(mssq(1,1,iq),'L->G',w1,lmmaxd, drotq(1,1,iq),lmmaxd)
    DO j=1,lmmaxd
      CALL zcopy(lmmaxd,w1(1,j),1,mssq(1,j,iq),1)
    END DO
  END IF
  
  
  
  
! ---> determine Delta_t = t(sys) - tmat(ref) = TSSQ - TREFLL
!      in the GLOBAL frame
  
  rf = refpot(iq)
!           write(*,*) 'RF',RF
  
  DO lm1 = 1,lmmaxd
    DO lm2 = 1,lmmaxd
      
      mssq(lm2,lm1,iq) = ( mssq(lm2,lm1,iq) - trefll(lm2,lm1,rf) )
      
    END DO
  END DO
  
  IF ( test('tmat    ') ) THEN
    WRITE(1337,*) 'IQ,IT,RF',iq,it,rf
    WRITE (1337,*) 'DELTA_TMATLL (',iq,' )'
    CALL cmatstr(' ',1,mssq(1,1,iq),lmmaxd,lmmaxd, 2*krel+1,2*krel+1,0,1D-8,6)
    WRITE (1337,*)
  END IF
  
  
END DO
! ----------------------------------------------------------------------
!    store the Delta_t matrix
! ----------------------------------------------------------------------
IF (opt('FERMIOUT') .AND. myrank==master) THEN             ! fswrt
  WRITE(6801,'(A)') 'TMATLL(ie):'                          ! fswrt
  DO iq=1,naez                                             ! fswrt
    DO lm2=1,lmmaxd                                        ! fswrt
      DO lm1=1,lmmaxd                                      ! fswrt
        WRITE(6801,'(2ES25.16)') mssq(lm1,lm2,iq)          ! fswrt
      END DO                                               ! fswrt
    END DO                                                 ! fswrt
  END DO                                                   ! fswrt
END IF                                                     ! fswrt
! ----------------------------------------------------------------------

!    MSSQ is now the Delta_t matrix in the GLOBAL frame
!    below, we determine (Delta_t)^(-1) in the GLOBAL frame

! ----------------------------------------------------------------------

!=======================================================================

! ---> loop over all atoms in unit cell, get Delta_t^(-1) = MSSQ

DO iq = 1,naez
  
! ---> inversion
  
  
  IF ( .NOT. opt('VIRATOMS') ) THEN
    IF (.NOT.test('testgmat')) THEN
      CALL zgetrf(lmmaxd,lmmaxd,mssq(1,1,iq),lmmaxd,ipvt,info)
!        CALL ZGETRI(LMMAXD,MSSQ(1,1,IQ),LMMAXD,IPVT,W1,
!     &              LMMAXD*LMMAXD,INFO)
      DO lm1=1,lmmaxd
        DO lm2=1,lmmaxd
          IF (lm1 == lm2) THEN
            w2(lm1,lm2)=cone
          ELSE
            w2(lm1,lm2)=czero
          END IF
        END DO
      END DO
      CALL zgetrs('N',lmmaxd,lmmaxd,mssq(1,1,iq),lmmaxd,ipvt,w2, lmmaxd,info)
      DO lm1=1,lmmaxd
        DO lm2=1,lmmaxd
          mssq(lm1,lm2,iq)=w2(lm1,lm2)
        END DO
      END DO
    END IF !( .not. TEST('testgmat') ) THEN
  END IF !( .not. OPT('VIRATOMS') ) THEN
  
  
END DO                    ! IQ = 1,NAEZ
!            stop


!============================================================IQ = 1,NAEZ
IF ((ncpa /= 0).AND.(icpastart == 2)) THEN
!----------------------------------------------------------------------
! s-, p-, and f-terms:    m(ata) = sum(q) c(q) * m(q)         >>>  AKAI
  
  lp=1
  ld=2
  lmp=(krel+1)*(lp+1)**2
  lmd=(krel+1)*(ld+1)**2 + 1
  DO iq = 1,naez
! ------------------------------------------ s,p blocks
    DO lm1 = 1, lmp
      DO lm2 = 1, lmp
        mssq(lm1,lm2,iq) = czero
        
        DO io=1,noq(iq)
          it=itoq(io,iq)
          mssq(lm1,lm2,iq) = mssq(lm1,lm2,iq) + conc(it) * msst(lm1,lm2,it)
        END DO
        
      END DO
    END DO
    
! ------------------------------------------ f block
    DO lm1 = lmd,lmmaxd
      DO lm2 = lmd,lmmaxd
        mssq(lm1,lm2,iq) = czero
        
        DO io=1,noq(iq)
          it=itoq(io,iq)
          mssq(lm1,lm2,iq) = mssq(lm1,lm2,iq) + conc(it) * msst(lm1,lm2,it)
        END DO
        
      END DO
    END DO
    
  END DO                 ! IQ=1,NAEZ
  
END IF                    ! ICPASTART.EQ.2

RETURN
END SUBROUTINE mssinit
