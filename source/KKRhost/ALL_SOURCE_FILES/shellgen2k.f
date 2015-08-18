C 23.2.2000/ 27.9.2004 *************************************************
      SUBROUTINE SHELLGEN2K(ICC,NATOM,RCLS,ATOM,NOFGIJ,IOFGIJ,JOFGIJ,
     &                      NROT,RSYMAT,ISYMINDEX,ROTNAME,
     &                      NSHELL,RATOM,NSH1,NSH2,ISH,JSH,
     &                      IJTABSYM,IJTABSH,IJTABCALC,
     &                      IPRINT,NSHELD)
C **********************************************************************
C *    Determines the number of different atomic pairs in a cluster by *
C * symmetry considerations, assigning a "shell" pointer (used to set  *
C * up the GF matrix), to each representative pair.                    *
C *                                                                    *
C * NATOM       number of atoms in the cluster                         *
C * RCLS(3,*)   atom positions in the cluster                          *
C * ATOM(*)     corresponding site index in the unit cell              *
C * NROT        actual number of symmetry operations                   *
C * RSYMAT      symmetry operation matrices                            *
C * ISYMINDEX   symmetry operation pointer                             *
C * NSHELD      dimension parameter ( max number of different shells)  *
C * IJTABCALC   flag to calculate the pair (I,J) - 1/0 for YES/NO      *
C *             (e.g. for impurity calc IJTABCALC(I,J) = 1 - delta_ij) *
C * NOFGIJ      total number of ij pairs (equals number of non-zero    *
C *             IJTABCALC elements                                     *
C * IOFGIJ      cluster indices i for pair ij                          *
C * JOFGIJ                      j for pair ij                          *
C *                                                                    *
C * NSHELL(0)   number of different shells (ij pairs)                  *
C * NSHELL(NS)  number of equivalent pairs in shell NS                 *
C * NSH1(NS),                                                          *
C * NSH2(NS)    site indices i,j of shell (representative pair) NS     *
C * ISH/JSH     cluster indices i,j of all NSHELL(NS) equivalent pairs *
C               described by shell NS                                  *
C * IJTABSH     the index of the representative shell NS for G_ij      *
C * IJTABSYM    the index of the symmetry operation which brings G(NS) *
C *             into G_ij                                              *
C * RATOM(3,NS) diference vector R_i(NS) - R_j(NS)                     *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Parameters
      INTEGER NSHELL0
      PARAMETER(NSHELL0 = 10000)
C     ..
C     .. Scalar arguments
      INTEGER ICC,NOFGIJ,NATOM,NROT,IPRINT,NSHELD
C     ..
C     .. Array arguments
      INTEGER ATOM(*),ISYMINDEX(*),IJTABSYM(*),IJTABSH(*),IJTABCALC(*)
      INTEGER NSHELL(0:NSHELD),NSH1(*),NSH2(*)
      INTEGER ISH(NSHELD,*),JSH(NSHELD,*)
      INTEGER IOFGIJ(*),JOFGIJ(*)
      DOUBLE PRECISION RCLS(3,*),RSYMAT(64,3,*)
      DOUBLE PRECISION RATOM(3,*)
      CHARACTER*10 ROTNAME(*)
C     ..
C     .. Local scalars
      INTEGER AI,AJ,I,J,K,NS,NSNEW,NSGEN,ID,ISYM,II,IJ,IGIJ
      DOUBLE PRECISION R1,SMALL
      LOGICAL LFOUND
C     ..
C     .. Local arrays
      DOUBLE PRECISION RI(3),RJ(3)
      INTEGER NSH1I(:),NSH2I(:),NSHELLI(:)
      DOUBLE PRECISION RATOMI(:,:)
      ALLOCATABLE NSH1I,NSH2I,NSHELLI,RATOMI
CF77CF77-------------------------------------------------------------
CF77      INTEGER NSH1I(NSHELL0),NSH2I(NSHELL0),NSHELLI(NSHELL0)
CF77      DOUBLE PRECISION RATOMI(3,NSHELL0)
CF77CF77-------------------------------------------------------------
CF90CF90-------------------------------------------------------------
CF90      INTEGER NSH1I(:),NSH2I(:),NSHELLI(:)
CF90      DOUBLE PRECISION RATOMI(:,:)
CF90      ALLOCATABLE NSH1I,NSH2I,NSHELLI,RATOMI
CF90CF90-------------------------------------------------------------
C     ..
C     .. Data statements
      DATA SMALL /  1.0D-10/
C     .. 
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99001)
      IF ( IPRINT.GT.1 ) CALL PRINTIJTAB(NATOM,IJTABCALC)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      IF ( NSHELD.GE.NSHELL0 ) THEN 
         WRITE(6,99000) 'local','NSHELL0',NSHELD
         STOP
      END IF
C
      IF ( NOFGIJ.LE.0 ) THEN
         WRITE(6,'(6X,"WARNING: no off-diagonal Gij elements found.",
     &                " ICC set to 0",/,
     &        6X,"         maybe you should check your input?",/)')
          ICC = 0 ! Bauer Long 2011-10-11
         RETURN
      END IF
CF90--------------------------------------------------------------------
      ALLOCATE(NSH1I(NSHELL0),NSH2I(NSHELL0),
     &         NSHELLI(NSHELL0),STAT=NS)
      IF ( NS.NE.0 ) STOP '   < shellgen2k > allocate NSHELLI arrays'
      ALLOCATE(RATOMI(3,NSHELL0),STAT=NS)
      IF ( NS.NE.0 ) STOP '   < shellgen2k > allocate RATOMI array'
CF90--------------------------------------------------------------------
C ======================================================================
C
C --> initialise number of shells found for this cluster, setup the
C     working arrays NSH1I,NSH2I,NSHELLI,RATOMI and set the number of
C     new found shells (NSNEW) to zero
C
      DO I = 1,NSHELL(0)
         NSH1I(I) = NSH1(I)
         NSH2I(I) = NSH2(I)
         NSHELLI(I) = NSHELL(I)
         DO J = 1,3
            RATOMI(J,I) = RATOM(J,I)
         END DO
      END DO
      NSNEW=0      
C
C **********************************************************************
C                                         loop over I,J-pairs in cluster
      DO IGIJ = 1,NOFGIJ
C
C --> search for a symmetric equivalent pair of atoms, LFOUND takes
C     on the value false/true if this equivalent pair is found
C
         I = IOFGIJ(IGIJ)
         J = JOFGIJ(IGIJ)
         AI = ATOM(I)
         AJ = ATOM(J)
C
         LFOUND = .FALSE.
         NSGEN = NSHELL(0) + NSNEW
C
C RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
         DO ID = 1,NROT
            ISYM = ISYMINDEX(ID)
C ----------------------------------------------------------------------
            DO II = 1,3
               RI(II) = RSYMAT(ISYM,II,1)*RCLS(1,I) +
     &                  RSYMAT(ISYM,II,2)*RCLS(2,I) + 
     &                  RSYMAT(ISYM,II,3)*RCLS(3,I)
C
               RJ(II) = RSYMAT(ISYM,II,1)*RCLS(1,J) +
     &                  RSYMAT(ISYM,II,2)*RCLS(2,J) + 
     &                  RSYMAT(ISYM,II,3)*RCLS(3,J)
            END DO
C ----------------------------------------------------------------------
C
C --> search for an equivalent pair within the already generated 
C     shells (1..NSHELL(0)+NSNEW)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            NS = 0
            DO WHILE ( ( .NOT.LFOUND ).AND.( NS.LT.NSGEN ) )
               NS = NS + 1
C ----------------------------------------------------------------------
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
               IF ( ( AI.EQ.NSH1I(NS) .AND. AJ.EQ.NSH2I(NS) ) ) THEN 
C
                  R1 = (RI(1)-RJ(1)+RATOMI(1,NS))**2  +
     &                 (RI(2)-RJ(2)+RATOMI(2,NS))**2  +
     &                 (RI(3)-RJ(3)+RATOMI(3,NS))**2
C
                  IF ( R1.LT.SMALL ) THEN
                     LFOUND = .TRUE.
                     NSHELLI(NS) = NSHELLI(NS) + 1
                     IF ( NS.LE.NSHELL(0) ) WRITE(6,99002) 
     &                    AI,(RCLS(II,I),II=1,3),
     &                    AJ,(RCLS(II,J),II=1,3),NS
                     ISH(NS,NSHELLI(NS)) = I
                     JSH(NS,NSHELLI(NS)) = J
                  END IF     
C
               END IF
C ----------------------------------------------------------------------
            END DO              ! NS = 1..NSGEN while .NOT.LFOUND
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         END DO
C RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
C --> if the rotation and the representative pair (shell) that
C     identify a pair of atoms was found LFOUND=.TRUE. and 
C     the search for a different pair of atoms starts; otherwise
C     the pair (I,J) requires a new shell
C
         IF ( .NOT.LFOUND ) THEN
            NSNEW = NSNEW + 1
            IF ( NSNEW+NSHELL(0).GT.NSHELL0 ) THEN
               WRITE(6,99000) 'local','NSHELL0',NSNEW+NSHELL(0)
               STOP 
            END IF
            IF ( NSNEW+NSHELL(0).GT.NSHELD ) THEN
               WRITE(6,99000) 'global','NSHELD',NSNEW+NSHELL(0)
               STOP 
            END IF
C
            NSH1I(NSHELL(0)+NSNEW) = AI
            NSH2I(NSHELL(0)+NSNEW) = AJ
            NSHELLI(NSHELL(0)+NSNEW) = 1
            ISH(NSHELL(0)+NSNEW,1) = I
            JSH(NSHELL(0)+NSNEW,1) = J
            DO II=1,3
               RATOMI(II,NSHELL(0)+NSNEW) = RCLS(II,J)-RCLS(II,I)
            END DO
         END IF
C
      END DO
C **********************************************************************
C
C --> test number of shells
C
      IF ( NSNEW+NSHELL(0).GT.NSHELD ) THEN
         WRITE(6,99000) 'global','NSHELD',NSNEW+NSHELL(0)
         STOP 
      END IF
C
C --> update the argument arrays
C 
      DO I = 1,NSHELL(0) + NSNEW
         NSH1(I) = NSH1I(I) 
         NSH2(I) = NSH2I(I) 
         NSHELL(I) = NSHELLI(I)
         DO J = 1,3
            RATOM(J,I) = RATOMI(J,I)
         END DO
      END DO
C
      NSHELL(0) = NSHELL(0) + NSNEW
CF90--------------------------------------------------------------------
      DEALLOCATE(NSH1I,NSH2I,NSHELLI,RATOMI,STAT=NS)
      IF ( NS.NE.0 ) STOP '   < shellgen2k > deallocate arrays'
CF90--------------------------------------------------------------------
C
C **********************************************************************
C
C --> scan once again the shells to find the corresponding symmetry
C     index bringing GS(1..NSHELL(0)) to Gij. 
C     Setup the tables IJTABSH  assigning (I,J) --> NS
C                      IJTABSYM assigning (I,J) --> ISYM
C     G_ij = D^\dagger(ISYM) * G(NS) * D(ISYM)
C     
C **********************************************************************
      DO I = 1,NATOM
         AI = (I-1)*NATOM
         DO J = 1,NATOM
            IJ = AI + J
            IJTABSH(IJ) = 0
            IJTABSYM(IJ) = 0
         END DO
      END DO
C **********************************************************************
      DO I = 1,NATOM
         AI = ATOM(I)
         DO J = 1,NATOM
            AJ = ATOM(J)
C=======================================================================
            DO II = 1,NSHELL(0)
C-----------------------------------------------------------------------
               DO ID = 1,NROT
                  ISYM = ISYMINDEX(ID)
C
                  DO K=1,3
                     RI(K) = RSYMAT(ISYM,K,1)*RATOM(1,II) +
     &                       RSYMAT(ISYM,K,2)*RATOM(2,II) + 
     &                       RSYMAT(ISYM,K,3)*RATOM(3,II)
                  ENDDO
C               
                  IF ( (AI.EQ.NSH1(II) .AND.
     +                  AJ.EQ.NSH2(II)) .OR.
     +                 (AI.EQ.NSH2(II) .AND. 
     +                  AJ.EQ.NSH1(II)) ) THEN                  
C
                     R1 = (RCLS(1,J)-RCLS(1,I)-RI(1))**2  +
     +                    (RCLS(2,J)-RCLS(2,I)-RI(2))**2  +
     +                    (RCLS(3,J)-RCLS(3,I)-RI(3))**2
C                  
                     IF (R1.LT.SMALL) THEN
                        IJ = (I-1)*NATOM + J 
                        IJTABSH(IJ) = II
                        IJTABSYM(IJ) = ID
                        GOTO 10
                     ENDIF
                  ENDIF
               ENDDO
C-----------------------------------------------------------------------
 10            CONTINUE
            END DO
C=======================================================================
         END DO
      END DO
C***********************************************************************
      IF ( IPRINT.LE.0 ) RETURN
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99003) 'assigned shells and symmetries'
      DO I = 1,NATOM
         AI = (I-1)*NATOM + J
         DO J = 1,NATOM
            IJ = AI + J
            IF ( IJTABCALC(IJ).NE.0 ) 
     &           WRITE(6,99004) I,J,IJTABSH(IJ),IJTABSYM(IJ),
     &                          ROTNAME(IJTABSYM(IJ))
         END DO
      END DO
      WRITE (6,99005)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
99000 FORMAT(6X,"Dimension ERROR: please increase the ",A," parameter",
     &     /,6X,A," to a value >=",I5,/)
99001 FORMAT (9X,'< SHELLGEN2K > : assigning representative pairs',
     &        ' (shells) ',/,26X,'for the off-diagonal elements Gij',/)
99002 FORMAT(9X,'INFO: For the atomic sites   I=',I3,' :',3F10.6,/,
     &       9X,29X,'J=',I3,' :',3F10.6,/,
     &       9X,6X,'an already generated equivalent shell (',I3,
     &             ') was found',/)
99003 FORMAT(13X,30(1H-),/,13X,A,/,13X,30(1H-),/,13X,
     &     " I ",1X," J "," | ","shell",4X,"isym",/,13X,30(1H-))
99004 FORMAT(13X,I3,1X,I3," | ",1X,I4,4X,I2,2X,A10)
99005 FORMAT(13X,30(1H-),/)
      END                           ! SUBROUTINE SHELLGEN
C
C **********************************************************************
C
      SUBROUTINE PRINTIJTAB(NATOM,IJTAB)
      IMPLICIT NONE
C     ..
      INTEGER NATOM
      INTEGER IJTAB(*)
C     ..
      INTEGER I,J,IJ
      INTEGER LGMAX
C     ..
      LGMAX = 59
      WRITE(6,99000) 
     &     '  searched for pairs marked with 1 in the table below'
      DO J = 1,MIN(NATOM+3,LGMAX)
         WRITE(6,'(1H-,$)')
      END DO
      WRITE(6,*)
      DO I = 1,NATOM
         WRITE(6,'(14X,I3," | ",$)') I
         IJ = (I-1)*NATOM
         DO J = 1,NATOM
            WRITE(6,'(I1,$)') IJTAB(IJ+J)
         END DO
         WRITE(6,*)
      END DO
      WRITE(6,'(13X,6(1H-),$)')
      DO J = 1,MIN(NATOM+3,LGMAX)
         WRITE(6,'(1H-,$)')
      END DO
      WRITE(6,'(/)')
C     ...........................................         
99000 FORMAT(13X,65(1H-),/,18X,A,/,13X,65(1H-),/,13X,
     &     "   J |",/,13X,"I    | 1..NATCLUS",/,13X,6(1H-),$)
      END 
