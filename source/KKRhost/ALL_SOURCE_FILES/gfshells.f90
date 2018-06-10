SUBROUTINE gfshells(icc,natomimp,nsh1,nsh2,  &
        ijtabsym,ijtabsh,ijtabcalc,  &
        iofgij,jofgij,nofgij,ish,jsh,  &
        nshell,naez,natyp,noq,rbasis,bravais,  &
        ifilimp,ratom,rclsimp,  &
        nsymat,isymindex,rsymat,  &
        kaoez,atomimp,  &
        rotname,hostimp,lmaxd,lmmaxd,  &
        naezd,natypd,natomimpd,nembd,nsheld)
! **********************************************************************
! *                                                                    *
! * This subroutine constructs mainly the index arrays                 *
! * NSHELL, NSH1, NSH2 -- NSHELL(0) number of different GF blocks that *
! * have to be calculated, NSH1(I),NSH2(I) the sites connected for     *
! * the block I, I = 1,NSHELL(0)                                       *
! *                                                                    *
! **********************************************************************
      use mod_types, only: t_imp
      IMPLICIT NONE
      INTEGER  LMAXD,LMMAXD,NAEZD,NATYPD,NATOMIMPD,NEMBD,NSHELD
!..
!.. Scalar arguments
      INTEGER ICC,NAEZ,NATOMIMP,NATYP,NSYMAT,NOFGIJ
      CHARACTER*40 IFILIMP
!..
!.. Array arguments
      CHARACTER*10 ROTNAME(64)
!..
      INTEGER ATOMIMP(NATOMIMPD),HOSTIMP(0:NATYPD)
      INTEGER ISYMINDEX(*),KAOEZ(NATYPD,NAEZD+NEMBD)
      INTEGER NOQ(NAEZD),NSH1(*),NSH2(*),NSHELL(0:NSHELD)
      INTEGER ISH(NSHELD,*),JSH(NSHELD,*)
      INTEGER IJTABSYM(*),IJTABSH(*),IJTABCALC(*),IOFGIJ(*),JOFGIJ(*)
!..
      DOUBLE PRECISION BRAVAIS(3,3),RATOM(3,NSHELD)
      DOUBLE PRECISION RBASIS(3,*),RCLSIMP(3,NATOMIMPD)
      DOUBLE PRECISION RSYMAT(64,3,3)
!.. 
!.. Local scalars
      INTEGER NB,I,J,POS,II,IO,NS,IN,NDIM,NSIZE,IHOST,IERR
      CHARACTER*9 STR9
      LOGICAL LSURF,OPT
!..
!.. External subroutines
      EXTERNAL IMPCHECK,IMPCOEFS,SHELLGEN2K,OPT

WRITE (1337,99000)

nsize = natomimpd*lmmaxd

! **********************************************************************

! --> construction of ratom, nsh1 and nsh2 for a self-consistent
!     calculation

IF ( .NOT. opt('VIRATOMS') ) THEN
  nshell(0) = natyp
ELSE
  nshell(0) = naez
END IF !( .not. OPT('VIRATOMS') ) THEN

IF ( nshell(0) > nsheld ) THEN
  WRITE(6,99001) 'NSHELD',nshell(0)
  STOP
END IF

DO i=1,nshell(0)
  ratom(1,i) = 0.0D0
  ratom(2,i) = 0.0D0
  ratom(3,i) = 0.0D0
  nshell(i) = 0
  
  DO j=1,naez
    DO io=1,noq(j)
      IF (kaoez(io,j) == i) THEN
        nshell(i) = nshell(i) + 1
        IF (nshell(i) == 1) THEN
          nsh1(i) = j
          nsh2(i) = j
        END IF
      END IF
    END DO
  END DO
  IF ( opt('VIRATOMS') ) THEN
    nshell(i)=1
    nsh1(i) = i
    nsh2(i) = i
  END IF
  
  
  
  IF (nshell(i) == 0) THEN
    WRITE(6,99002)
    STOP
  END IF
END DO

IF ( icc == 0 ) THEN
  WRITE(1337,99003) nshell(0)
  RETURN
END IF

!      end of simple SCF-calculation part.
! **********************************************************************

!heck if we are in surface mode

lsurf = .false.
IF ( bravais(1,3) == 0D0 .AND. bravais(2,3) == 0D0 .AND.  &
    bravais(3,3) == 0D0 ) lsurf = .true.
ndim = 3
IF (lsurf) ndim = 2

! **********************************************************************
!      NATOMIMP=0   ! BUG: This initialization breaks the shell generation for
!                   ! ICC=-1, which is set by option XCPL.  B. Zimmermann

IF (icc < 0) THEN
  
! --->  ICC.LT.1 all shells are (should be) prepared
  
  WRITE(1337,99011) natomimp
ELSE
  
! --> read-in the cluster coordinates from an external file
  
  REWIND 25
  READ (25,FMT=*) natomimp
  
  IF (natomimp > natomimpd ) THEN
    WRITE(6,99001) 'NATOMIMPD',natomimp
    STOP
  END IF
  WRITE(1337,99004) ifilimp,natomimp
  
  DO i=1,natomimp
    READ (25,FMT=*) (rclsimp(j,i),j=1,3),atomimp(i)
    atomimp(i) = atomimp(i) + icc - 1
  END DO
  
  IF (opt('GREENIMP') .OR. opt('OPERATOR')) THEN
    ihost=0
    DO  i=1,natypd
      DO j=1,natomimp
        IF (atomimp(j) == i) THEN
          ihost=ihost+1
          hostimp(ihost)=atomimp(j)
          CYCLE
        END IF
      END DO
    END DO
    
! save stuff to t_imp for later use
    t_imp%ihost = ihost
    t_imp%natomimp = natomimp
    allocate(t_imp%hostimp(ihost), stat=ierr)
    IF(ierr/=0) STOP 'Error allocating t_imp%HOSTIMP'
    t_imp%hostimp(1:ihost) = hostimp(1:ihost)
    allocate(t_imp%atomimp(natomimp), stat=ierr)
    IF(ierr/=0) STOP 'Error allocating t_imp%ATOMIMP'
    t_imp%atomimp(1:natomimp) = atomimp(1:natomimp)
    
  END IF!GREENIMP
  
END IF ! ICC>=0
! **********************************************************************

CALL impcheck(atomimp,natomimp,naez,rclsimp,rbasis,bravais,ndim)

! **********************************************************************
IF ( icc > 0 ) THEN
  WRITE(1337,99005)
  
! --> set up the number of all (I,J)-pairs to be looked for,
!     avoid considering again the diagonal elements
  
  nofgij = 0
  DO i = 1,natomimp
    nb = (i-1)*natomimp
    DO j = 1,natomimp
      ijtabcalc(nb+j) = 0
    END DO
    IF ( atomimp(i) >= 0 ) THEN
      DO j = 1,natomimp
        IF ( ( atomimp(j) >= 0 ).AND.( i /= j ) ) THEN
          nofgij = nofgij + 1
          IF ( nofgij > natomimpd*natomimpd ) THEN
            WRITE(6,99001) 'NATOMIMPD',nofgij/natomimp
            STOP
          END IF
          iofgij(nofgij) = i
          jofgij(nofgij) = j
          ijtabcalc(nb+j) = 1
        END IF
      END DO
    END IF
  END DO
END IF
! **********************************************************************

CALL shellgen2k(icc,natomimp,rclsimp(1,1),atomimp(1), nofgij,iofgij,jofgij,  &
    nsymat,rsymat,isymindex,rotname, nshell,ratom(1,1),nsh1,nsh2,ish,jsh,  &
    ijtabsym,ijtabsh,ijtabcalc,2,nsheld)

! **********************************************************************

! --> now write out the impurity.coefs file for the impurity calculation
!                                                         n.papanikolaou

IF ( icc > 0 .OR. opt('KKRFLEX '))  &
    CALL impcoefs(natomimp,naez,atomimp,rclsimp,nshell,  &
    nsh1,nsh2,ratom,nsymat,isymindex,rotname, hostimp,natypd,lmaxd,nsheld,nsize)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE(1337,99003) nshell(0)
WRITE(1337,99006)
nb = MAX(natyp,naez)
DO ns = 1,nshell(0)
  IF ( ns == nb+1 ) WRITE(1337,99012)
  IF ( ns <= nb ) THEN
    CALL setpairstr(nsh1(ns),nsh2(ns),str9)
    WRITE(1337,99007) ns,nsh1(ns),nsh2(ns),(ratom(ii,ns),ii=1,3),  &
        SQRT(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3,ns)**2), str9
  ELSE
    WRITE(1337,99008) ns,nsh1(ns),nsh2(ns),(ratom(ii,ns),ii=1,3),  &
        SQRT(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3,ns)**2)
    io = MIN(2,nshell(ns))
    DO i = 1,io
      CALL setpairstr(ish(ns,i),jsh(ns,i),str9)
      WRITE(1337,'(A9,$)') str9
    END DO
    WRITE(1337,*)
    pos = (nshell(ns)+1)/2
    DO i = 2,pos
      io = (i-1)*2
      in = MIN(2,nshell(ns)-io)
      WRITE(1337,99009)
      DO j = 1,in
        CALL setpairstr(ish(ns,io+j),jsh(ns,io+j),str9)
        WRITE(1337,'(A9,$)') str9
      END DO
      WRITE(1337,*)
    END DO
  END IF
END DO
WRITE(1337,'(6X,72(1H-))')
nb = 0
DO ns=1,nshell(0)
  nb = nb + nshell(ns)
END DO
WRITE(1337,99010) nb
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ----------------------------------------------------------------------
99000 FORMAT(5X,"< GFSHELLS > : setting up indices of the GF blocks",/)
99001 FORMAT(6X,"Dimension ERROR: please increase the global parameter",  &
    /,6X,a," to a value >=",i5,/)
99002 FORMAT(6X,"ERROR: there are some inconsistencies in your input",/,  &
    13X,"not all atoms defined by NATYP have been found",/)
99003 FORMAT(8X,"Different shells for GF calculation : ",i3,/)
99004 FORMAT(8X,"Reading in cluster impurity sites from file",/,  &
    12X,"file name        : ",a,/,12X,"atoms in cluster : ",i3)
99005 FORMAT(8X, "Preparing indexing for impurity GF",/,11X,  &
    "- unsymmetrised GF is written out (v. 20.09.2001)",/,11X,  &
    "- files that will be created: impurity.coefs",/,41X,  &
    "intercell_ref",/,41X,"green",/)
99006 FORMAT(6X,72(1H-),/,6X,"shell|"," IQ ",  &
    " JQ"," | ",10X,"vec R_IJ ",11X,"R_IJ   | equiv. pairs",/, 6X,72(1H-))
99007 FORMAT(5X,i5," |",i3,1X,i3," | ",3F9.4,f9.5,1X,"|",a9)
99008 FORMAT(5X,i5," |",i3,1X,i3," | ",3F9.4,f9.5,1X,"|",$)
99009 FORMAT(5X,5X," |",7X," | ",27X,9X,1X,"|",$)
99010 FORMAT(8X,"Number of block elements to be calculated : ",i3,/)
99011 FORMAT(8X,"Setting pairs for task-defined cluster sites ",  &
    "and connections",/, 12X,"atoms in cluster : ",i3)
99012 FORMAT(6X,72(1H:),/,  &
    22X,"(impurity) cluster related data/indexing",/, 6X,72(1H:))

END SUBROUTINE gfshells
! **********************************************************************

SUBROUTINE setpairstr(i,j,str9)
      IMPLICIT NONE
      CHARACTER*9 STR9,STRD
      INTEGER I,J,L,LSTR
      CHARACTER*20 FMT1
!     ..
fmt1 = '("(",I'
fmt1 = fmt1(1:6)//'1'
lstr = 4
IF ( i >= 10 ) THEN
  fmt1 = fmt1(1:6)//'2'
  lstr = lstr + 1
  IF (i >= 100 ) THEN
    fmt1 = fmt1(1:6)//'3'
    lstr = lstr + 1
  END IF
END IF
fmt1 = fmt1(1:7)//',",",I'
fmt1 = fmt1(1:13)//'1'
lstr = lstr + 1
IF ( j >= 10 ) THEN
  fmt1 = fmt1(1:13)//'2'
  lstr = lstr + 1
  IF ( j >= 100 ) THEN
    fmt1 = fmt1(1:13)//'3'
    lstr = lstr + 1
  END IF
END IF
fmt1 = fmt1(1:14)//',")")'
WRITE(strd,fmt1) i,j
DO l = 1,9-lstr
  str9(l:l) = ' '
END DO
str9 = str9(1:9-lstr)//strd(1:lstr)
END SUBROUTINE setpairstr
! **********************************************************************
