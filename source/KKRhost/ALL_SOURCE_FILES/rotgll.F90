SUBROUTINE rotgll(gmatll,natomimp,ijtabsym,ijtabsh,  &
        dsymll,symunitary,igf,rc,crel,rrel,  &
        krel,lmmaxd,irec)
! **********************************************************************
! *                                                                    *
! *   it calculates all the elements of the Green Function of          *
! *   the impurity cluster using the GF calculated for the             *
! *   representative pairs.                                            *
! *   the representative pair and the symmetry operation D are given   *
! *   through the arrays IJTABSH and IJTABSYM set up in < SHELLGEN2K > *
! *                                                                    *
! *     _ _                                                            *
! *     n n'                      n n'                                 *
! *     m m'               T      m m'                                 *
! *    G    (E) = SUM    D     * G    (E) * D                          *
! *     L L'      L1 L2   L L1    L1L2      L2 L'                      *
! *                                                                    *
! *                   _                 _                              *
! *              n    n            n'   n'                             *
! *   where   D R  = R     and  D R  = R                               *
! *              m    m            m    m                              *
! *                                                                    *
! **********************************************************************
use mod_mympi, only: myrank, master

IMPLICIT NONE
!     ..
!     .. Parameter definitions
DOUBLE COMPLEX czero,cone
PARAMETER (czero= (0.0D0,0.0D0),cone= (1.d0,0.d0))
!     ..
!     .. Scalar arguments
INTEGER :: ngclus,lmmaxd,irec
INTEGER :: igf,krel,natomimp
!     ..
!     .. Array arguments
INTEGER :: ijtabsym(*),ijtabsh(*)
DOUBLE COMPLEX gmatll(lmmaxd,lmmaxd,*), dsymll(lmmaxd,lmmaxd,*)
DOUBLE COMPLEX crel(lmmaxd,lmmaxd),rc(lmmaxd,lmmaxd), rrel(lmmaxd,lmmaxd)
LOGICAL :: symunitary(*)
!     ..
!     .. Local arrays
DOUBLE COMPLEX gll(:,:,:,:),tpg(:,:)
COMPLEX*8 gclust(:)
allocatable gll,tpg,gclust
!     ..
!     .. Local scalars
INTEGER :: ilin,iq,icall,ish,isym,jq
INTEGER :: lm1,lm2,nlin, ilm, jlm
CHARACTER (LEN=1) :: cnt
CHARACTER (LEN=4) :: str4i,str4j
CHARACTER (LEN=18) :: str18
!     ..
!     .. External Subroutines
EXTERNAL changerep,cmatstr,zgemm,opt
!     ..
!     .. External Functions
LOGICAL :: test,opt
EXTERNAL test
!     ..
!     .. Data statement
DATA icall / 1 /
!     ..
!     .. Save statement
SAVE icall
!     ..
allocate (gll(lmmaxd,lmmaxd,natomimp,natomimp), tpg(lmmaxd,lmmaxd),stat=lm1)
IF ( lm1 /= 0 ) THEN
  WRITE(6,99001) ' GLL/TPG'
  STOP '           < ROTGLL > '
endif
99001 FORMAT(6X,"ERROR: failed to allocate array(s) :",a,/)
! **********************************************************************
IF ( icall == 1) THEN
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  IF(myrank==master) WRITE (1337,'(79("="))')
  IF(myrank==master) WRITE (1337,'(6X,2A)')  &
      'ROTGLL : Expand GF for all pairs by rotation',  &
      ' and write out (all E-points)'
  IF(myrank==master) WRITE (1337,'(79("="))')
  IF(myrank==master) WRITE (1337,*)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
endif
!***********************************************************************
DO iq = 1,natomimp
  DO jq = 1,natomimp
!-----------------------------------------------------------------------
    ilin = (iq-1)*natomimp + jq
    ish = ijtabsh(ilin)
    isym = ijtabsym(ilin)
!-----------------------------------------------------------------------
!    for REL CASE look if it is a unitary / ANTI - unitary rotation
!-----------------------------------------------------------------------
    cnt = 'N'
    IF ( .NOT.symunitary(isym) ) cnt = 'T'
    
    CALL zgemm('C',cnt,lmmaxd,lmmaxd,lmmaxd,cone,  &
        dsymll(1,1,isym),lmmaxd,gmatll(1,1,ish), lmmaxd,czero,tpg,lmmaxd)
    
    CALL zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,  &
        tpg,lmmaxd,dsymll(1,1,isym), lmmaxd,czero,gll(1,1,iq,jq),lmmaxd)
!-----------------------------------------------------------------------
  END DO
END DO
!***********************************************************************

!     visualise Gij

IF ( test('Gmatij  ') ) THEN
  WRITE (1337,'(/,4X,70("+"),/,4X,A,I4)')  &
      'cluster G_ij matrices for i,j = 1,',natomimp
  
  DO iq = 1,natomimp
    WRITE(1337,'(/,8X,66("="))')
    DO jq = 1,natomimp
      
      WRITE(str4i,'(I4)') iq
      WRITE(str4j,'(I4)') jq
      str18 = '   i ='//str4i(1:4)//' j ='//str4j(1:4)
      IF (krel == 0) THEN
        CALL cmatstr(str18,18,gll(1,1,iq,jq),lmmaxd,lmmaxd, 0,0,0,1D-8,6)
      ELSE
        CALL changerep(gll(1,1,iq,jq),'REL>RLM',tpg,lmmaxd,  &
            lmmaxd,rc,crel,rrel,str18,18)
      endif
      IF (jq < natomimp) WRITE(1337,'(/,9X,65("-"))')
    END DO
    IF (iq == natomimp) WRITE(1337,'(/,8X,66("="),/)')
  END DO
  WRITE(1337,'(4X,70("+"))')
endif

!***********************************************************************

! --> output of GF

ngclus=lmmaxd*natomimp
IF ( igf /= 0 ) THEN
  icall = icall + 1
  nlin = 0
  
  allocate (gclust(ngclus*ngclus),stat=lm1)
  IF ( lm1 /= 0 ) THEN
    WRITE(6,99001) ' GCLUST'
    STOP '           < ROTGLL > '
  endif
  DO  jq=1,natomimp
    DO lm2=1,lmmaxd
      jlm = (jq-1)*lmmaxd+lm2
      DO  iq=1,natomimp
        DO lm1=1,lmmaxd
          nlin = nlin + 1
          IF ( nlin > ngclus*ngclus )  &
              STOP "<ROTGLL>: NLIN.GT.(NATOMIMP*LMMAXD)**2"
          gclust(nlin) = gll(lm1,lm2,iq,jq)
!test
!                    WRITE(214321,'(4i,2E)') LM1,LM2,IQ,JQ,GCLUST(NLIN)
! writeout of green_host for WRTGREEN option
          IF(opt('WRTGREEN') .AND. myrank==master) THEN
            ilm=(iq-1)*lmmaxd+lm1
            WRITE(58,'((2I5),(2e17.9))') jlm, ilm, gll(lm1,lm2,iq,jq)
          endif
        END DO
      END DO
    END DO
  END DO
  
  
  
  
  IF ( ( opt('KKRFLEX ') ) ) THEN
#ifdef CPP_MPI
    irec = irec
#ELSE
    irec = icall
#ENDIF
  WRITE(888,REC=irec) gclust
  IF ( ( opt('GPLAIN  ') ) ) THEN
    WRITE(8888,'(50000E25.16)') gclust
  endif
endif

!==== the following write-out has been disabled, because it was assumed to be     !no-green
!====  obsolete with the implementation of the MPI-communicated arrays. If I am   !no-green
!====  wrong and the write-out is needed in subsequent parts, construct a         !no-green
!====  test-option around it so that it is only written out in this case.         !no-green
!        IF ( .not. OPT('KKRFLEX ') ) THEN                                        !no-green
!          WRITE(88,REC=ICALL) GCLUST                                             !no-green
!        endif                                                                   !no-green

deallocate (gclust)
endif !IGF/=0
!***********************************************************************
deallocate (gll,tpg)
RETURN
END                       ! SUBROUTINE ROTGLL
