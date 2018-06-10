SUBROUTINE decimaread(ez,tk,nptp1,nptp2,nptp3,npol,ispin,  &
    lefttinvll,righttinvll,vacflag,ienergy,  &
    nlbasis,nrbasis,naez,kaoez,kmrot,ins,nspin,  &
    lmmax,ielast,fileleft,fileright, krel,natypd,lmmaxd,nembd1)!,KORBIT)
! **********************************************************************
! *                                                                    *
! * This subroutine reads in the t-matrices of the left                *
! * and right host for the decimation method.                          *
! *                                                                    *
! * The t-matrices are writen out in kloopz1  (option 'deci-out')      *
! *                                                                    *
! * The host files contain the CMOMHOST data neeed to create the       *
! * interatomic potential in subroutine < vinterface >                 *
! * This is going to be read in < cmomsread >, the decimation files    *
! * are for this reason not rewinded.                                  *
! *                                                                    *
! * In case of 'vacuum' setting on one of the sides, no energy points  *
! * are read in and the VACLAG is set to TRUE.                         *
! *                                                                    *
! * A call of this routine with IENERGY = 0 means that we are not in   *
! * the energy loop - only the header of each decimation file read in  *
! *                                                                    *
! * IENERGY <> 0 reads in the matrices/energy at IENERGY -> returned   *
! *                                                                    *
! **********************************************************************
      use mod_version_info
      IMPLICIT NONE
!..
      INTEGER KREL,NATYPD,NEMBD1,LMMAXD!,KORBIT
!..
!.. Scalar arguments
      INTEGER NPTP1,NPTP2,NPTP3,NPOL,ISPIN,IENERGY
      INTEGER NLBASIS,NRBASIS,NAEZ,KMROT,INS,NSPIN,LMMAX,IELAST
      DOUBLE PRECISION TK
      CHARACTER*40 FILELEFT,FILERIGHT  
!..
!.. Array arguments
      INTEGER KAOEZ(NATYPD,*)
      DOUBLE COMPLEX EZ(*)
      DOUBLE COMPLEX LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1), &
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1)
      LOGICAL VACFLAG(2)
!..
!.. Local variables
      INTEGER IOS,NPT1L,NPT2L,NPT3L,NPOLL,INSL
      INTEGER NSPINL,NAEZL,LMMAXL,KRELL,KMROTL
      INTEGER IEL,IH,LM1,LM2,I,IDUM,IH1,IHOST,ILHOST,NATHOST
      DOUBLE PRECISION ALATL,TEMPL,E2L,QMT,QMP
      DOUBLE COMPLEX EL
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD)
      CHARACTER*80 BANER1
      CHARACTER*40 FILEHOST
      CHARACTER*5 CHHOST(2),STR5
      DOUBLE PRECISION BRAVAISL(3,3),RBASISL(3,NEMBD1)
!..
!.. External Functions
      INTEGER LNGSTRING
      EXTERNAL LNGSTRING
!..
!.. External Subroutines
      EXTERNAL ZCOPY
!..
!.. Data statements
      DATA CHHOST/'LEFT ','RIGHT'/
!     ..
! ========================================================== IENERGY = 0
IF ( ienergy < 1 ) THEN
  WRITE (1337,'(5X,A,/,8X,65(1H-))') 'Reading in host Delta_t matrices'
  vacflag(1) = .false.
  vacflag(2) = .false.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
  DO ihost = 1,2
    filehost = fileleft
    nathost = nlbasis
    IF ( ihost == 2 ) THEN
      filehost = fileright
      nathost = nrbasis
    END IF
    ilhost = lngstring(filehost,40)
    
    WRITE (1337,'(8X,A5," side host: ",$)') chhost(ihost)
! ----------------------------------------------------------------------
    IF ( filehost(1:7) == 'vacuum' ) THEN
      WRITE (1337,'(A)') 'VACUUM will be used'
      vacflag(ihost) = .true.
! ----------------------------------------------------------------------
    ELSE
      WRITE (1337,'(A,/)') filehost(1:ilhost)
      OPEN (36+ihost,FILE=filehost,STATUS='OLD',IOSTAT=ios)
      CALL version_check_header(36+ihost)
! ......................................................................
      IF ( ios > 0 ) THEN
        WRITE (6,'(/,5X,2A)') 'ERROR: Can not open host file '  &
            ,filehost(1:ilhost)
        WRITE (6,'(12X,A,A)') 'vacuum    entry should be used to set'  &
            ,' a vacuum-host on this side'
        STOP '       < DECIMAREAD > '
      END IF
! ......................................................................
      DO ih = 1,3
        READ (36+ihost,99007) baner1
      END DO
      READ (36+ihost,99009) alatl,nspinl,naezl,lmmaxl,insl, krell,kmrotl
! ......................................................................
      IF ( (krell /= krel) .OR. (kmrotl /= kmrot) .OR.  &
            (insl /= ins) .OR. (nspinl /= nspin) .OR.  &
            (lmmaxl /= lmmax) .OR. (naezl /= nathost) ) THEN
        WRITE (6,'(/,5X,2A)') 'ERROR: ',  &
            'host not compatible with your input/sytem'
        WRITE (6,'(14X,6(A6),/,8X,42(1H-))') '  KREL',  &
            ' KMROT','   INS',' NSPIN',' LMMAX',' BASIS'
        WRITE (6,'(8X,A6,6I6)') 'syst: ',krel,kmrot,ins,nspin, lmmax,nathost
        WRITE (6,'(8X,A6,6I6,/)') 'host: ',krell,kmrotl,insl,  &
            nspinl,lmmaxl,naezl
        STOP '       < DECIMAREAD > '
      END IF
! ......................................................................
      WRITE (1337,99001) alatl,nspinl,naezl,lmmaxl,insl,krell, kmrotl
      
      READ (36+ihost,99007) baner1
      READ (36+ihost,99003) bravaisl
      WRITE (1337,99002) bravaisl
      READ (36+ihost,99007) baner1
      ih = lngstring(baner1,80)
      WRITE (1337,99008) baner1(1:ih)
! ......................................................................
      IF ( krel == 0 ) THEN
        DO ih = 1,naezl
          READ (36+ihost,99003) (rbasisl(i,ih),i=1,3)
          WRITE (1337,99003) (rbasisl(i,ih),i=1,3)
        END DO
      ELSE
        DO ih = 1,naezl
          READ (36+ihost,99014) (rbasisl(i,ih),i=1,3),qmt,qmp
          WRITE (1337,99004) (rbasisl(i,ih),i=1,3),qmt,qmp
        END DO
      END IF
! ......................................................................
      READ (36+ihost,99010) e2l,templ
      READ (36+ihost,99011) npt1l,npt2l,npt3l,npoll
! ......................................................................
      IF ( (npt1l /= nptp1) .OR. (npt2l /= nptp2) .OR.  &
            (npt3l /= nptp3) .OR. (npoll /= npol) .OR.  &
            (ABS(templ-tk) > 1.d-4) ) THEN
        WRITE (6,'(/,5X,2A)') 'ERROR: Host energy mesh',  &
            ' not compatible with your input'
        WRITE (6,'(14X,5(A6),/,8X,40(1H-))') '   NE1',  &
            '   NE2','   NE3','  NPOL','  TEMP'
        WRITE (6,'(8X,A6,4I6,2X,F10.6)') 'input:',nptp1,nptp2, nptp3,npol,tk
        WRITE (6,'(8X,A6,4I6,2X,F10.6,/)') ' host:',npt1l,  &
            npt2l,npt3l,npoll,templ
        STOP '       < DECIMAREAD > '
      END IF
! ......................................................................
      WRITE (1337,99005) e2l,templ
      WRITE (1337,99006) npt1l,npt2l,npt3l,npoll
    END IF
! ----------------------------------------------------------------------
    WRITE (1337,'(8X,65(1H-))')
  END DO
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
!                            Headers read in
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  RETURN
END IF
! ========================================================== IENERGY = 0


! ======================================================================
!                  READ IN Delta_t matrices for IENERGY > 0
! ======================================================================
IF ( (ispin /= 2) .OR. (nspin /= 1) ) THEN
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
  DO ihost = 1,2
    nathost = nlbasis
    IF ( ihost == 2 ) nathost = nrbasis
! ----------------------------------------------------------------------
! --> read in Delta_t if it is not a vacuum host
!     attention: new format Dec. 2004
    
    IF ( .NOT.vacflag(ihost) ) THEN
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ HOST ATOM LOOP
      DO ih = 1,nathost
        CALL cinit(lmmaxd*lmmaxd,w1)
        
        9910             CONTINUE
        READ (36+ihost,'(A5)',END=300) str5
        IF ( str5 /= '*****' ) GO TO 9910
        
        READ (36+ihost,99012) iel,el,idum
! ......................................................................
        IF ( (ABS(el-ez(ienergy)) > 1.d-4) .OR. (iel /= ienergy) ) THEN
          WRITE (6,'(/,5X,2A)') 'ERROR: Host energy mesh',  &
              ' not compatible with your input'
          WRITE (6,'(14X,2(A6),/,8X,43(1H-))') '    IE', ' E(IE)'
          WRITE (6,'(8X,A6,I6,1X,2(1X,F10.6))') 'input:', ienergy,ez(ienergy)
          WRITE (6,'(8X,A6,I6,1X,2(1X,F10.6),/)') ' host:', iel,el
          STOP '       < DECIMAREAD > '
        END IF
! ......................................................................
        IF ( idum /= ih ) THEN
          WRITE (6,'(/,5X,2A,/)') 'ERROR reading host file',  &
              ' basis indexing wrong'
          STOP '       < DECIMAREAD > '
        END IF
! ......................................................................
        
! --> map the t-matrix corectly
        
        ih1 = kaoez(1,naez+(ihost-1)*nlbasis+ih)
        
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        9920             CONTINUE
        READ(36+ihost,99013) lm1,lm2,el
        IF ( (lm1+lm2) /= 0 ) THEN
          w1(lm1,lm2) = el
          IF ( (lm1+lm2) < 2*lmmaxd ) GO TO 9920
        END IF
        
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        IF ( ihost == 1 ) THEN
          DO lm1 = 1,lmmaxd  !(KREL+KORBIT+1)*LMMAX
!                         CALL ZCOPY((KREL+KORBIT+1)*LMMAX,W1(1,LM1),1,
            CALL zcopy(lmmaxd,w1(1,lm1),1, lefttinvll(1,lm1,ih1),1)
          END DO
        ELSE
          DO lm1 = 1,lmmaxd  !(KREL+KORBIT+1)*LMMAX
!                         CALL ZCOPY((KREL+KORBIT+1)*LMMAX,W1(1,LM1),1,
            CALL zcopy(lmmaxd,w1(1,lm1),1, righttinvll(1,lm1,ih1),1)
          END DO
        END IF
      END DO
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    END IF              ! not vacuum
! ----------------------------------------------------------------------
  END DO
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
END IF                    ! ispin.ne.2 .or. .nspin.ne.1
! ======================================================================

IF ( ienergy == ielast ) WRITE (1337,*)
RETURN
300  CONTINUE
STOP '        Error reading hostfile'

99001 FORMAT (10X,'ALAT=',f9.6,' NSPIN=',i2,'  NAEZ=',i3,' LMMAX=',i3,  &
    ' INS=',i1,' KREL=',i1,' KMROT=',i1)
99002 FORMAT (10X,'BRAVAIS '/10X,3F8.4/10X,3F8.4/10X,3F8.4)
99003 FORMAT (10X,3F8.4)
99004 FORMAT (10X,3F8.4,2F9.4)
99005 FORMAT (10X,'EF=',f10.6,' TEMP=',f10.4,' Kelvin')
99006 FORMAT (10X,'N1=',i3,' N2=',i3,' N3=',i3,' NPOL=',i3)

99007 FORMAT (a80)
99008 FORMAT (10X,a)
99009 FORMAT (5X,f9.6,7X,i2,7X,i3,7X,i3,5X,i1,6X,i1,7X,i1)
99010 FORMAT (3X,f10.6,6X,f10.4)
99011 FORMAT (3X,i3,4X,i3,4X,i3,6X,i3)
99012 FORMAT (7X,i5,2D16.8,6X,i3)
99013 FORMAT (2I5,1P,2D22.14)
99014 FORMAT (3F8.4,2F9.4)

END SUBROUTINE decimaread
