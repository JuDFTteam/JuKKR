      MODULE MOD_RWLDAUPOT
      CONTAINS

      SUBROUTINE RWLDAUPOT(LWRITE,IPOS,NATYP,NSPIN,LMAXD,IRMD,IRUNLDAU,LDAU)
! Reads or writes the LDA+U arrays to the disk.
!
! LWRITE = .TRUE.  : write out
! LWRITE = .FALSE. : read in
!
! IPOS = 1,2,3,4 : up to which position to read or write, in the order:
!
! IRUNLDAU,LOPT,UEFF,JEFF,EREFLDAU,NATYP
! ----- 1 -----
! WLDAU
! ----- 2 -----
! ULDAU
! ----- 3 -----
! PHI
! ----- 4 -----
      USE TYPE_LDAU 
      implicit none
! Input:
      INTEGER                        :: IPOS,NATYP,NSPIN,LMAXD,IRMD
      LOGICAL                        :: LWRITE
! Input/Output:
      TYPE(LDAU_TYPE),ALLOCATABLE    :: LDAU(:)  ! lda+u variables, intended dimension: (natyp) ! lda+u
      INTEGER                        :: IRUNLDAU
! Inside
      CHARACTER*20                   :: TEXT,TEXT1
      INTEGER                        :: IR,M1,M2,M3,M4,I2,IS,I1,MMAXD,LMMAXD,NSPIND
 
      MMAXD = 2*LMAXD + 1
      LMMAXD = (LMAXD + 1)**2
      NSPIND = 2

      OPEN (67,FILE='ldaupot',FORM='FORMATTED')
      REWIND (67)

      IF (LWRITE) THEN

          WRITE (67,*) IRUNLDAU,                                          &
              (LDAU(I1)%LOPT,I1=1,NATYP),(LDAU(I1)%UEFF,I1=1,NATYP),      &
              (LDAU(I1)%JEFF,I1=1,NATYP),(LDAU(I1)%EREFLDAU,I1=1,NATYP)   

          IF (IPOS.EQ.1) RETURN

          WRITE (67,9010) 'WLDAU'
          DO I1 = 1,NATYP
             IF (LDAU(I1)%LOPT.GE.0) THEN
             WRITE (67,9020) 'ATOM ', I1
             WRITE (67,9000) (((LDAU(I1)%WLDAU(M1,M2,IS),M1=1,MMAXD),M2=1,MMAXD),IS=1,NSPIN)
             ENDIF
          ENDDO

          IF (IPOS.EQ.2) RETURN

          WRITE (67,9010) 'ULDAU'
          DO I1 = 1,NATYP
             IF (LDAU(I1)%LOPT.GE.0) THEN
                WRITE (67,9020) 'ATOM ', I1
                WRITE (67,9000) ((((LDAU(I1)%ULDAU(M1,M2,M3,M4),M1=1,MMAXD),M2=1,MMAXD),M3=1,MMAXD),M4=1,MMAXD)
             ENDIF
          ENDDO

          IF (IPOS.EQ.3) RETURN

          WRITE (67,9010) 'PHI  '
          DO I1 = 1,NATYP
             IF (LDAU(I1)%LOPT.GE.0) THEN
                WRITE (67,9020) 'ATOM ', I1
                WRITE (67,9005) (LDAU(I1)%PHI(IR),IR=1,IRMD)
             ENDIF
          ENDDO

      ELSE  ! Read from file
          READ (67,*) IRUNLDAU,                                         &
              (LDAU(I1)%LOPT,I1=1,NATYP),(LDAU(I1)%UEFF,I1=1,NATYP),    &
              (LDAU(I1)%JEFF,I1=1,NATYP),(LDAU(I1)%EREFLDAU,I1=1,NATYP) 

! Allocate everything that is not allocated but should be:
          DO I1 = 1,NATYP
             IF (LDAU(I1)%LOPT.GE.0) THEN
                IF (.NOT.ALLOCATED(LDAU(I1)%WLDAU)) ALLOCATE( LDAU(I1)%WLDAU(MMAXD,MMAXD,NSPIN) )
                IF (.NOT.ALLOCATED(LDAU(I1)%ULDAU)) ALLOCATE( LDAU(I1)%ULDAU(MMAXD,MMAXD,MMAXD,MMAXD) )
                IF (.NOT.ALLOCATED(LDAU(I1)%PHI))   ALLOCATE( LDAU(I1)%PHI(IRMD) )
                IF (.NOT.ALLOCATED(LDAU(I1)%CUTOFF)) ALLOCATE( LDAU(I1)%CUTOFF(IRMD) )
                IF (.NOT.ALLOCATED(LDAU(I1)%DENMATC)) ALLOCATE( LDAU(I1)%DENMATC(MMAXD,MMAXD,NSPIN) )
             ENDIF
          ENDDO

          IF (IPOS.EQ.1) RETURN

          READ (67,*) TEXT
          DO I1 = 1,NATYP
             IF (LDAU(I1)%LOPT.GE.0) THEN
                READ (67,*) TEXT,TEXT1
                READ (67,*) (((LDAU(I1)%WLDAU(M1,M2,IS),M1=1,MMAXD),M2=1,MMAXD),IS=1,NSPIN)
             ENDIF
          ENDDO

          IF (IPOS.EQ.2) RETURN

          READ (67,*) TEXT
          DO I1 = 1,NATYP
             IF (LDAU(I1)%LOPT.GE.0) THEN
                READ (67,*) TEXT,TEXT1
                READ (67,*) ((((LDAU(I1)%ULDAU(M1,M2,M3,M4),M1=1,MMAXD),M2=1,MMAXD),M3=1,MMAXD),M4=1,MMAXD)
             ENDIF
          ENDDO

          IF (IPOS.EQ.3) RETURN

          READ (67,*) TEXT
          DO I1 = 1,NATYP
             IF (LDAU(I1)%LOPT.GE.0) THEN
                READ (67,*) TEXT,TEXT1
                READ (67,9005) (LDAU(I1)%PHI(IR),IR=1,IRMD)
             ENDIF
          ENDDO
      ENDIF

      CLOSE (67)

 9000 FORMAT (4D22.14)
 9005 FORMAT (4E20.12)
 9010 FORMAT (A5)
 9020 FORMAT (A5,I5)

      RETURN
      END SUBROUTINE RWLDAUPOT


      END MODULE MOD_RWLDAUPOT
