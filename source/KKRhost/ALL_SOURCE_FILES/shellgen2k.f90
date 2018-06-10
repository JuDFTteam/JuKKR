! 23.2.2000/ 27.9.2004 *************************************************
SUBROUTINE shellgen2k(icc,natom,rcls,atom,nofgij,iofgij,jofgij,  &
    nrot,rsymat,isymindex,rotname, nshell,ratom,nsh1,nsh2,ish,jsh,  &
    ijtabsym,ijtabsh,ijtabcalc, iprint,nsheld)
! **********************************************************************
! *    Determines the number of different atomic pairs in a cluster by *
! * symmetry considerations, assigning a "shell" pointer (used to set  *
! * up the GF matrix), to each representative pair.                    *
! *                                                                    *
! * NATOM       number of atoms in the cluster                         *
! * RCLS(3,*)   atom positions in the cluster                          *
! * ATOM(*)     corresponding site index in the unit cell              *
! * NROT        actual number of symmetry operations                   *
! * RSYMAT      symmetry operation matrices                            *
! * ISYMINDEX   symmetry operation pointer                             *
! * NSHELD      dimension parameter ( max number of different shells)  *
! * IJTABCALC   flag to calculate the pair (I,J) - 1/0 for YES/NO      *
! *             (e.g. for impurity calc IJTABCALC(I,J) = 1 - delta_ij) *
! * NOFGIJ      total number of ij pairs (equals number of non-zero    *
! *             IJTABCALC elements                                     *
! * IOFGIJ      cluster indices i for pair ij                          *
! * JOFGIJ                      j for pair ij                          *
! *                                                                    *
! * NSHELL(0)   number of different shells (ij pairs)                  *
! * NSHELL(NS)  number of equivalent pairs in shell NS                 *
! * NSH1(NS),                                                          *
! * NSH2(NS)    site indices i,j of shell (representative pair) NS     *
! * ISH/JSH     cluster indices i,j of all NSHELL(NS) equivalent pairs *
!               described by shell NS                                  *
! * IJTABSH     the index of the representative shell NS for G_ij      *
! * IJTABSYM    the index of the symmetry operation which brings G(NS) *
! *             into G_ij                                              *
! * RATOM(3,NS) diference vector R_i(NS) - R_j(NS)                     *
! *                                                                    *
! **********************************************************************
IMPLICIT NONE
!..
!.. Parameters
INTEGER NSHELL0
PARAMETER(NSHELL0 = 10000)
!..
!.. Scalar arguments
INTEGER ICC,NOFGIJ,NATOM,NROT,IPRINT,NSHELD
!..
!.. Array arguments
INTEGER ATOM(*),ISYMINDEX(*),IJTABSYM(*),IJTABSH(*),IJTABCALC(*)
INTEGER NSHELL(0:NSHELD),NSH1(*),NSH2(*)
INTEGER ISH(NSHELD,*),JSH(NSHELD,*)
INTEGER IOFGIJ(*),JOFGIJ(*)
DOUBLE PRECISION RCLS(3,*),RSYMAT(64,3,*)
DOUBLE PRECISION RATOM(3,*)
CHARACTER*10 ROTNAME(*)
!..
!.. Local scalars
INTEGER AI,AJ,I,J,K,NS,NSNEW,NSGEN,ID,ISYM,II,IJ,IGIJ
DOUBLE PRECISION R1,SMALL
LOGICAL LFOUND
!..
!.. Local arrays
DOUBLE PRECISION RI(3),RJ(3)
INTEGER NSH1I(:),NSH2I(:),NSHELLI(:)
DOUBLE PRECISION RATOMI(:,:)
ALLOCATABLE NSH1I,NSH2I,NSHELLI,RATOMI
!..
!.. Data statements
DATA SMALL /  1.0D-10/
!..

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99001)
IF ( iprint > 1 ) CALL printijtab(natom,ijtabcalc)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

IF ( nsheld >= nshell0 ) THEN
  WRITE(6,99000) 'local','NSHELL0',nsheld
  STOP
END IF

IF ( nofgij <= 0 ) THEN
  WRITE(6,'(6X,"WARNING: no off-diagonal Gij elements found.",  &
      " ICC set to 0",/, 6X,"         maybe you should check your input?",/)')
  icc = 0 ! Bauer Long 2011-10-11
  RETURN
END IF
allocate(nsh1i(nshell0),nsh2i(nshell0), nshelli(nshell0),stat=ns)
IF ( ns /= 0 ) STOP '   < shellgen2k > allocate NSHELLI arrays'
allocate(ratomi(3,nshell0),stat=ns)
IF ( ns /= 0 ) STOP '   < shellgen2k > allocate RATOMI array'
! ======================================================================

! --> initialise number of shells found for this cluster, setup the
!     working arrays NSH1I,NSH2I,NSHELLI,RATOMI and set the number of
!     new found shells (NSNEW) to zero

DO i = 1,nshell(0)
  nsh1i(i) = nsh1(i)
  nsh2i(i) = nsh2(i)
  nshelli(i) = nshell(i)
  DO j = 1,3
    ratomi(j,i) = ratom(j,i)
  END DO
END DO
nsnew=0

! **********************************************************************
!                                         loop over I,J-pairs in cluster
DO igij = 1,nofgij
  
! --> search for a symmetric equivalent pair of atoms, LFOUND takes
!     on the value false/true if this equivalent pair is found
  
  i = iofgij(igij)
  j = jofgij(igij)
  ai = atom(i)
  aj = atom(j)
  
  lfound = .false.
  nsgen = nshell(0) + nsnew
  
! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
  DO id = 1,nrot
    isym = isymindex(id)
! ----------------------------------------------------------------------
    DO ii = 1,3
      ri(ii) = rsymat(isym,ii,1)*rcls(1,i) + rsymat(isym,ii,2)*rcls(2,i) +  &
          rsymat(isym,ii,3)*rcls(3,i)
      
      rj(ii) = rsymat(isym,ii,1)*rcls(1,j) + rsymat(isym,ii,2)*rcls(2,j) +  &
          rsymat(isym,ii,3)*rcls(3,j)
    END DO
! ----------------------------------------------------------------------
    
! --> search for an equivalent pair within the already generated
!     shells (1..NSHELL(0)+NSNEW)
    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ns = 0
    DO WHILE ( ( .NOT.lfound ).AND.( ns < nsgen ) )
      ns = ns + 1
! ----------------------------------------------------------------------
!              IF ( ( AI.EQ.NSH1I(NS) .AND. AJ.EQ.NSH2I(NS) ).OR.
!    &              ( AI.EQ.NSH2I(NS) .AND. AJ.EQ.NSH1I(NS) )  ) THEN
! Commented out by Phivos Mavropoulos 31 Oct 2008. The problem is that if (I,J) and (J,I)
! are assigned to the same shell, then G(I,J) should be transposed to obtain G(J,I).
! However, this transposition is not performed in account in kkr1b (subr. tbxccpljij).
! There, only the real-space rotations (DSYMLL) are performed to generate each pair GF from the
! representative pair, but the transposition is forgotten. Thus there are two ways to resolve this:
! Either flag the pairs to be transposed, which is is a little faster but complicated
! to program, or do not consider the (I,J) and (J,I) pairs as belonging to the same shell,
! which is done now:
      IF ( ( ai == nsh1i(ns) .AND. aj == nsh2i(ns) ) ) THEN
        
        r1 = (ri(1)-rj(1)+ratomi(1,ns))**2  +  &
            (ri(2)-rj(2)+ratomi(2,ns))**2  + (ri(3)-rj(3)+ratomi(3,ns))**2
        
        IF ( r1 < small ) THEN
          lfound = .true.
          nshelli(ns) = nshelli(ns) + 1
          IF ( ns <= nshell(0) ) WRITE(1337,99002) ai,(rcls(ii,i),ii=1,3),  &
              aj,(rcls(ii,j),ii=1,3),ns
          ish(ns,nshelli(ns)) = i
          jsh(ns,nshelli(ns)) = j
        END IF
        
      END IF
! ----------------------------------------------------------------------
    END DO              ! NS = 1..NSGEN while .NOT.LFOUND
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END DO
! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
  
! --> if the rotation and the representative pair (shell) that
!     identify a pair of atoms was found LFOUND=.TRUE. and
!     the search for a different pair of atoms starts; otherwise
!     the pair (I,J) requires a new shell
  
  IF ( .NOT.lfound ) THEN
    nsnew = nsnew + 1
    IF ( nsnew+nshell(0) > nshell0 ) THEN
      WRITE(6,99000) 'local','NSHELL0',nsnew+nshell(0)
      STOP
    END IF
    IF ( nsnew+nshell(0) > nsheld ) THEN
      WRITE(6,99000) 'global','NSHELD',nsnew+nshell(0)
      STOP
    END IF
    
    nsh1i(nshell(0)+nsnew) = ai
    nsh2i(nshell(0)+nsnew) = aj
    nshelli(nshell(0)+nsnew) = 1
    ish(nshell(0)+nsnew,1) = i
    jsh(nshell(0)+nsnew,1) = j
    DO ii=1,3
      ratomi(ii,nshell(0)+nsnew) = rcls(ii,j)-rcls(ii,i)
    END DO
  END IF
  
END DO
! **********************************************************************

! --> test number of shells

IF ( nsnew+nshell(0) > nsheld ) THEN
  WRITE(6,99000) 'global','NSHELD',nsnew+nshell(0)
  STOP
END IF

! --> update the argument arrays

DO i = 1,nshell(0) + nsnew
  nsh1(i) = nsh1i(i)
  nsh2(i) = nsh2i(i)
  nshell(i) = nshelli(i)
  DO j = 1,3
    ratom(j,i) = ratomi(j,i)
  END DO
END DO

nshell(0) = nshell(0) + nsnew
deallocate(nsh1i,nsh2i,nshelli,ratomi,stat=ns)
IF ( ns /= 0 ) STOP '   < shellgen2k > deallocate arrays'

! **********************************************************************

! --> scan once again the shells to find the corresponding symmetry
!     index bringing GS(1..NSHELL(0)) to Gij.
!     Setup the tables IJTABSH  assigning (I,J) --> NS
!                      IJTABSYM assigning (I,J) --> ISYM
!     G_ij = D^\dagger(ISYM) * G(NS) * D(ISYM)

! **********************************************************************
DO i = 1,natom
  ai = (i-1)*natom
  DO j = 1,natom
    ij = ai + j
    ijtabsh(ij) = 0
    ijtabsym(ij) = 0
  END DO
END DO
! **********************************************************************
DO i = 1,natom
  ai = atom(i)
  DO j = 1,natom
    aj = atom(j)
!=======================================================================
    DO ii = 1,nshell(0)
!-----------------------------------------------------------------------
      DO id = 1,nrot
        isym = isymindex(id)
        
        DO k=1,3
          ri(k) = rsymat(isym,k,1)*ratom(1,ii) +  &
              rsymat(isym,k,2)*ratom(2,ii) + rsymat(isym,k,3)*ratom(3,ii)
        END DO
        
        IF ( (ai == nsh1(ii) .AND. aj == nsh2(ii)) .OR.  &
              (ai == nsh2(ii) .AND. aj == nsh1(ii)) ) THEN
          
          r1 = (rcls(1,j)-rcls(1,i)-ri(1))**2  +  &
              (rcls(2,j)-rcls(2,i)-ri(2))**2  + (rcls(3,j)-rcls(3,i)-ri(3))**2
          
          IF (r1 < small) THEN
            ij = (i-1)*natom + j
            ijtabsh(ij) = ii
            ijtabsym(ij) = id
            GO TO 10
          END IF
        END IF
      END DO
!-----------------------------------------------------------------------
      10            CONTINUE
    END DO
!=======================================================================
  END DO
END DO
!***********************************************************************
IF ( iprint <= 0 ) RETURN

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99003) 'assigned shells and symmetries'
DO i = 1,natom
  ai = (i-1)*natom + j
  DO j = 1,natom
    ij = ai + j
    IF ( ijtabcalc(ij)>0 ) WRITE(1337,99004) i,j,ijtabsh(ij),ijtabsym(ij),  &
        rotname(ijtabsym(ij))
  END DO
END DO
WRITE (1337,99005)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

99000 FORMAT(6X,"Dimension ERROR: please increase the ",a," parameter",  &
    /,6X,a," to a value >=",i5,/)
99001 FORMAT (9X,'< SHELLGEN2K > : assigning representative pairs',  &
    ' (shells) ',/,26X,'for the off-diagonal elements Gij',/)
99002 FORMAT(9X,'INFO: For the atomic sites   I=',i3,' :',3F10.6,/,  &
    9X,29X,'J=',i3,' :',3F10.6,/,  &
    9X,6X,'an already generated equivalent shell (',i3, ') was found',/)
99003 FORMAT(13X,30(1H-),/,13X,a,/,13X,30(1H-),/,13X,  &
    " I ",1X," J "," | ","shell",4X,"isym",/,13X,30(1H-))
99004 FORMAT(13X,i3,1X,i3," | ",1X,i4,4X,i2,2X,a10)
99005 FORMAT(13X,30(1H-),/)
END                           ! SUBROUTINE SHELLGEN

! **********************************************************************

SUBROUTINE printijtab(natom,ijtab)
IMPLICIT NONE
!     ..
INTEGER :: natom
INTEGER :: ijtab(*)
!     ..
INTEGER :: i,j,ij
INTEGER :: lgmax
!     ..
lgmax = 59
WRITE(1337,99000) '  searched for pairs marked with 1 in the table below'
DO j = 1,MIN(natom+3,lgmax)
  WRITE(1337,'(1H-,$)')
END DO
WRITE(1337,*)
DO i = 1,natom
  WRITE(1337,'(14X,I3," | ",$)') i
  ij = (i-1)*natom
  DO j = 1,natom
    WRITE(1337,'(I1,$)') ijtab(ij+j)
  END DO
  WRITE(1337,*)
END DO
WRITE(1337,'(13X,6(1H-),$)')
DO j = 1,MIN(natom+3,lgmax)
  WRITE(1337,'(1H-,$)')
END DO
WRITE(1337,'(/)')
!     ...........................................
99000 FORMAT(13X,65(1H-),/,18X,a,/,13X,65(1H-),/,13X,  &
    "   J |",/,13X,"I    | 1..NATCLUS",/,13X,6(1H-),$)
END SUBROUTINE printijtab
