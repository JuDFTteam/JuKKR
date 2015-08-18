C**********************************************************************
      SUBROUTINE MSSINIT(NCPA,ICPASTART,TSST,MSST,MSSQ,TREFLL,DROTQ,
     &                   REFPOT,IQAT,ITOQ,NOQ,CONC,
     &                   KMROT,NATYP,NAEZ,LMMAXD) ! nrefd was taken out of calling list 1.2.2012
C
      IMPLICIT NONE
      include 'inc.p' ! Included  1.2.2012
C
C     .. Parameters
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO = (0.0D0,0.0D0) )
      PARAMETER (CONE = (1.0D0,0.0D0) )
C
C     .. Dummy arguments
      INTEGER KMROT, NATYP, NAEZ, LMMAXD, NCPA, ICPASTART
      INTEGER IQAT(NATYPD),ITOQ(NATYPD,NAEZD)
      INTEGER REFPOT(NAEZD),NOQ(NAEZD)
      DOUBLE PRECISION CONC(NATYPD)
      DOUBLE COMPLEX 
     +     TSST(LMMAXD,LMMAXD,NATYPD),
     +     TREFLL(LMMAXD,LMMAXD,NREFD)
      DOUBLE COMPLEX MSST(LMMAXD,LMMAXD,NATYPD)
      DOUBLE COMPLEX MSSQ(LMMAXD,LMMAXD,NAEZD)
      DOUBLE COMPLEX DROTQ(LMMAXD,LMMAXD,NAEZD)
C
C     .. Local variables
      INTEGER IT,IQ,RF,J,IO,INFO,LM1,LM2,LP,LD,LMP,LMD
      INTEGER IPVT(LMMAXD)
      DOUBLE COMPLEX ZC
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD)
      DOUBLE COMPLEX W2(LMMAXD,LMMAXD)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST,OPT
      EXTERNAL TEST,OPT
C ======================================================================
C
C   --> set up the Delta_t^-1 matrix (MSST) in the LOCAL frame
C
!           write(*,*) 'natyp',natyp



      DO IT = 1,NATYP
C     
         DO J = 1,LMMAXD
            CALL ZCOPY(LMMAXD,TSST(1,J,IT),1,MSST(1,J,IT),1)
         END DO

         IQ = IQAT(IT)
         RF = REFPOT(IQ)


         IF (KMROT.NE.0) THEN
           CALL ROTATE(TREFLL(1,1,RF),'G->L',W1,LMMAXD,
     &                 DROTQ(1,1,IQ),LMMAXD)
         ELSE
           DO J=1,LMMAXD
             CALL ZCOPY(LMMAXD,TREFLL(1,J,RF),1,W1(1,J),1)
           END DO
         END IF

C
C ---> determine Delta_t = t(sys) - tmat(ref) = TSST - TREFLL 
C      in local frame
C
         DO LM1 = 1,LMMAXD
           DO LM2 = 1,LMMAXD
C
             MSST(LM2,LM1,IT) = 
     &                        ( MSST(LM2,LM1,IT) - W1(LM2,LM1) )
C
           END DO
         END DO
C
C  --> inversion
C
      IF ( .not. OPT('VIRATOMS') ) THEN
        IF (.not.TEST('testgmat')) THEN
         CALL ZGETRF(LMMAXD,LMMAXD,MSST(1,1,IT),LMMAXD,IPVT,INFO)
c         CALL ZGETRI(LMMAXD,MSST(1,1,IT),LMMAXD,IPVT,W1,
c     &               LMMAXD*LMMAXD,INFO)
        DO LM1=1,LMMAXD
         DO LM2=1,LMMAXD
          IF (LM1.EQ.LM2) THEN
           W2(LM1,LM2)=CONE
          ELSE
           W2(LM1,LM2)=CZERO
          ENDIF
         ENDDO         
        ENDDO         
        CALL ZGETRS('N',LMMAXD,LMMAXD,MSST(1,1,IT),LMMAXD,IPVT,W2,
     &              LMMAXD,INFO)
        DO LM1=1,LMMAXD
         DO LM2=1,LMMAXD
          MSST(LM1,LM2,IT)=W2(LM1,LM2)
         ENDDO
        ENDDO
      END IF !( TEST('testgmat') ) THEN
      END IF !( OPT('VIRATOMS') ) THEN
C    
      END DO
      

C ======================================================================
C
C ---> determine tmat(sys) - tmat(ref) = MSSQ - TREFLL
C           
C      because TSST is calculated in the LOCAL frame, if KMROT<>0
C      it needs to be rotated prior to set the Delta matrix
C
C ---> set up the effective (on-site) Delta_t- and Delta_m-matrices 
C      using the Average T-matrix Approximation
C
C   ICPASTART=1:     
C       m(IQ) = t(ATA) = SUM(it)  c(it) * t(it)
C       m(IQ) = t(ATA) - t_ref = Delta_t(ATA)
C       m(IQ) = (Delta_t(ATA))^(-1)
C
C   ICPASTART=2:
C       m(IQ) = (Delta_t(ATA))^(-1) for l = 2
C       m(IQ) = SUM(it) c(it)*m(it) with m(it)=t(it)^(-1)
C
C       mssq(IQ)  refer to the GLOBAL frame
C       tsst(IT),msst(IT)  refer to the LOCAL  frame
C
C
C ----------------------------------------------------------------------
C
      CALL CINIT(LMMAXD*LMMAXD*NAEZ,MSSQ)

      DO IQ = 1,NAEZ
C
        DO IO = 1,NOQ(IQ)
          IT = ITOQ(IO,IQ)
!         write(*,*) 'test fivos mssinit IO,IQ,IT',IO,IQ,IT ! test fivos
          ZC = CONC(IT)
          DO J = 1,LMMAXD
            CALL ZAXPY(LMMAXD,ZC,TSST(1,J,IT),1,MSSQ(1,J,IQ),1)
          END DO
        END DO
C
C ---> rotate MSSQ from the LOCAL to the GLOBAL frame if necessary
C
        IF ( KMROT.NE.0 ) THEN
          CALL ROTATE(MSSQ(1,1,IQ),'L->G',W1,LMMAXD,
     &                DROTQ(1,1,IQ),LMMAXD)
          DO J=1,LMMAXD
            CALL ZCOPY(LMMAXD,W1(1,J),1,MSSQ(1,J,IQ),1)
          END DO
        END IF



C
C ---> determine Delta_t = t(sys) - tmat(ref) = TSSQ - TREFLL 
C      in the GLOBAL frame
C
        RF = REFPOT(IQ)       
!           write(*,*) 'RF',RF

        DO LM1 = 1,LMMAXD
          DO LM2 = 1,LMMAXD
C
            MSSQ(LM2,LM1,IQ) = 
     &                       ( MSSQ(LM2,LM1,IQ) - TREFLL(LM2,LM1,RF) )
C
          END DO
        END DO
C
        IF ( TEST('tmat    ') ) THEN
          WRITE(*,*) 'IQ,IT,RF',IQ,IT,RF
          WRITE (6,*) 'DELTA_TMATLL (',IQ,' )'
          CALL CMATSTR(' ',1,MSSQ(1,1,IQ),LMMAXD,LMMAXD,
     &                 2*KREL+1,2*KREL+1,0,1d-8,6)
          WRITE (6,*)
        END IF
C
!            open(23234,file='test1')
!            do lm1=1,LMMAXD
!              write(23234,'(50000F)') TREFLL(lm1,:,RF)
!            end do !lm1=1,alm
!            close(23234)
!            write(*,*) NCPA
!            write(*,*) 'stop'
!               stop


      END DO
C ----------------------------------------------------------------------
C
C    MSSQ is now the Delta_t matrix in the GLOBAL frame
C    below, we determine (Delta_t)^(-1) in the GLOBAL frame
C
C ----------------------------------------------------------------------

C=======================================================================
C
C ---> loop over all atoms in unit cell, get Delta_t^(-1) = MSSQ
C
      DO IQ = 1,NAEZ
C
C ---> inversion 
C

!            open(23234,file='test1')
!            do lm1=1,LMMAXD
!              write(23234,'(50000F)') MSSQ(lm1,:,1)
!            end do !lm1=1,alm
!            close(23234)
!            write(*,*) NCPA
!            write(*,*) 'stop'
! !            stop

      IF ( .not. OPT('VIRATOMS') ) THEN
        IF (.not.TEST('testgmat')) THEN
        CALL ZGETRF(LMMAXD,LMMAXD,MSSQ(1,1,IQ),LMMAXD,IPVT,INFO)
c        CALL ZGETRI(LMMAXD,MSSQ(1,1,IQ),LMMAXD,IPVT,W1,
c     &              LMMAXD*LMMAXD,INFO)
        DO LM1=1,LMMAXD
         DO LM2=1,LMMAXD
          IF (LM1.EQ.LM2) THEN
           W2(LM1,LM2)=CONE
          ELSE
           W2(LM1,LM2)=CZERO
          ENDIF
         ENDDO         
        ENDDO         
        CALL ZGETRS('N',LMMAXD,LMMAXD,MSSQ(1,1,IQ),LMMAXD,IPVT,W2,
     &              LMMAXD,INFO)
        DO LM1=1,LMMAXD
         DO LM2=1,LMMAXD
          MSSQ(LM1,LM2,IQ)=W2(LM1,LM2)
         ENDDO
        ENDDO
       END IF !( .not. TEST('testgmat') ) THEN
      END IF !( .not. OPT('VIRATOMS') ) THEN

!            open(23234,file='test2')
!            do lm1=1,LMMAXD
!              write(23234,'(50000F)') MSSQ(lm1,:,IQ)
!            end do !lm1=1,alm
!            close(23234)
!            write(*,*) NCPA
!            write(*,*) 'stop'
!            stop

C
      END DO                    ! IQ = 1,NAEZ
!            stop

C
C============================================================IQ = 1,NAEZ
      IF ((NCPA.NE.0).AND.(ICPASTART.EQ.2)) THEN
C----------------------------------------------------------------------
C s-, p-, and f-terms:    m(ata) = sum(q) c(q) * m(q)         >>>  AKAI
C
         LP=1
         LD=2
         LMP=(KREL+1)*(LP+1)**2
         LMD=(KREL+1)*(LD+1)**2 + 1
         DO IQ = 1,NAEZ
C ------------------------------------------ s,p blocks
            DO LM1 = 1, LMP
               DO LM2 = 1, LMP
                  MSSQ(LM1,LM2,IQ) = CZERO
C
                  DO IO=1,NOQ(IQ)
                     IT=ITOQ(IO,IQ)
                     MSSQ(LM1,LM2,IQ) = MSSQ(LM1,LM2,IQ) +
     &                    CONC(IT) * MSST(LM1,LM2,IT)
                  END DO
C
               END DO
            END DO
C
C ------------------------------------------ f block
            DO LM1 = LMD,LMMAXD
               DO LM2 = LMD,LMMAXD
                  MSSQ(LM1,LM2,IQ) = CZERO
C
                  DO IO=1,NOQ(IQ)
                     IT=ITOQ(IO,IQ)
                     MSSQ(LM1,LM2,IQ) = MSSQ(LM1,LM2,IQ) +
     &                    CONC(IT) * MSST(LM1,LM2,IT)
                  END DO
C
               END DO
            END DO
C
         END DO                 ! IQ=1,NAEZ
C
      END IF                    ! ICPASTART.EQ.2
C
      RETURN
      END
