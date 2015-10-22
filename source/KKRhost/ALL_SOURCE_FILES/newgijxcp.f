C*==gijxcpl.f    processed by SPAG 6.05Rc at 15:40 on 18 Oct 2004
       SUBROUTINE GIJXCPL(IDO,NAEZ,RBASIS,BRAVAIS,LINTERFACE,
     &                   NIQCALC,IQCALC,NATOMIMP,RCLSIMP,ATOMIMP,
     &                   IJTABCALC,NATOMIMPD)
C **********************************************************************
C *                                                                    *
C * In case of tasks requiring Gij blocks calculation, set variables:  *
C *                                                                    *
C * NATOMIMP, RCLSIMP(3,1..NATOMIMP), ATOMIMP(1..NATOMIMP)             *
C * IJTABCALC flag to which pair is needed: I,J --> (I-1)*NATOMIMP + J *
C *           indexing refers to the generated cluster                 *
C * NIQCALC   number of sites in the first unit cell contained in the  *
C *           cluster                                                  *
C * IQCALC()  correspondence index to the 1..NAEZ sites                *
C * IDO takes on the value 1 or 0 if setting up process was OK or not  *
C *                                                                    *
C * EXCHANGE COUPLING CONSTANTS calculation case                       *
C *                                                                    *
C **********************************************************************
       IMPLICIT NONE
C ..
C ..  Arguments
       INTEGER IDO,NAEZ,NATOMIMP,NIQCALC,NATOMIMPD
       INTEGER ATOMIMP(*),IJTABCALC(*),IQCALC(*)
       DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),RCLSIMP(3,*)
       LOGICAL LINTERFACE
C ..
C ..  Locals
       INTEGER I,I1,I2,I3,IBR(3),IEQVEC,IQ,IQS,J,JQ,MM,NBR(3),NDIM,
     &        NN,NOUT,JQS
C.......................................................................
C     uniquely identify a vector R_j - r_i = (r_j + T_n) - r_i by a 
C     quintuple integer in IVECI2J(5,*):
C     index  meaning
C        1    unit-cell site of atom i (always i - first unit cell)
C        2    unit-cell site of atom j 
C        3..5 the three coeficients of the translation vector
C             T_n = n1*a1 + n2*a2 + n3*a3
C     NVECI2J(I) gives the number of identical vectors Rj (initially 1)
C                in the array IVECI2J(2/3/4/5,I)
C     IREF(1..NVECI2J(I),I) points to the NVECI2J(I) identical vectors 
C     in the array IVECI2J (one site is kept only once)
C.......................................................................
CF90----------------------------------------------------------------
       INTEGER NB3MAX
       INTEGER NJQCALC
       INTEGER IVECI2J(:,:),NVECI2J(:),IREF(:,:),JQCALC(:)
       ALLOCATABLE IVECI2J,NVECI2J,IREF,JQCALC
CF90----------------------------------------------------------------
CF77CF77----------------------------------------------------------------
CF77      INTEGER NBMAX,NB3MAX
CF77      PARAMETER ( NBMAX = 4, NB3MAX = NBMAX**3 )
CF77      INTEGER IVECI2J(5,NAEZ*NB3MAX),NVECI2J(NAEZ*NB3MAX)
CF77      INTEGER IREF(NAEZ*NB3MAX,NAEZ*NB3MAX)
CF77CF77----------------------------------------------------------------
CF90CF90----------------------------------------------------------------
CF90      INTEGER NB3MAX
CF90      INTEGER IVECI2J(:,:),NVECI2J(:),IREF(:,:)
CF90      ALLOCATABLE IVECI2J,NVECI2J,IREF
CF90CF90----------------------------------------------------------------
       DOUBLE PRECISION CLURAD,CLURADSQ,DQ(3),DR(3),DRSQ,TOL,TOLSQ
       DOUBLE PRECISION CLURADXY,CLURADXYSQ,DRXYSQ
       LOGICAL LSPHER
       CHARACTER*256 UIO ! NCOLIO=256
C     ..
C     .. Externals
       EXTERNAL GETCLUSNXYZ,IOINPUT
C     ..
       IDO = 0
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
       WRITE (1337,99000)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
       TOL = 1.0D-4
       TOLSQ = TOL*TOL
       NDIM = 3
       IF ( LINTERFACE ) NDIM = 2
C
       IQ = 0
       CLURAD = 0D0
       LSPHER = .TRUE.
       CALL IOINPUT('JIJRAD          ',UIO,0,7,IQ)
       IF ( IQ.EQ.0 ) THEN
          READ (UNIT=UIO,FMT=*) CLURAD
          IF ( CLURAD.LT.0D0 ) CLURAD = 0D0
          CLURADXY = CLURAD
          IF ( CLURAD.GT.0D0 ) THEN
             IQ = 0
             CALL IOINPUT('JIJRADXY        ',UIO,0,7,IQ)
             IF ( IQ.EQ.0 ) THEN
                READ (UNIT=UIO,FMT=*) CLURADXY
                IF ( CLURADXY.LE.0D0 ) CLURADXY = CLURAD
             END IF
          END IF
          LSPHER = ( ABS(CLURAD-CLURADXY).LT.TOL )
       ELSE
          WRITE (1337,99001)
       END IF
C
       DO I = 1,3
          NBR(I) = 0
       END DO
       CLURADXYSQ = MAX(CLURAD,CLURADXY)
       CALL GETCLUSNXYZ(CLURADXYSQ,BRAVAIS,NDIM,CLURADSQ,NBR)
C
       CLURADSQ = (CLURAD+TOL)*(CLURAD+TOL)
       CLURADXYSQ = (CLURADXY+TOL)*(CLURADXY+TOL)
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
       IF ( CLURAD.GT.0D0 ) THEN
          IF ( LSPHER ) THEN
             WRITE (1337,99002) CLURAD,((2*NBR(I)+1),I=1,3)
          ELSE
             WRITE (1337,99003) CLURADXY,CLURAD,((2*NBR(I)+1),I=1,3)
          END IF
       ELSE
          WRITE (1337,99004)
       END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C --> set the reference I and J sites for J_IJ within the unit cell 
C     tokens JIJSITEI and JIJSITEJ allows one to specify which sites 
C     in the unit cell should be considered in the calculation 
C     (default ALL)
C     Usage:     JIJSITEx  NUMBER_OF_SITES SITE_INDICES_1_TO_NOS
C     with x=I/J
C     Examples:  assume a unit cell with 4 atomic sites 
C           JIJSITEI= 1 1
C        will take only 1 'I' site having the index 1, as sites 'J'
C        all the others in the u.c.
C           JIJSITEI= 2 1 3  JIJSITEJ= 1 4
C        will take 2 sites 'I', 1st and 3rd, and 1 site 'J', the 4th
C        hence the calculated pairs will be (1,4),(3,4)
C
CF90----------------------------------------------------------------
       ALLOCATE(JQCALC(NAEZ),STAT=IQ)
       IF ( IQ.NE.0 ) STOP '    Allocate JQCALC'
CF90----------------------------------------------------------------
       NIQCALC = NAEZ
       NJQCALC = NIQCALC
       DO IQ = 1,NAEZ
          IQCALC(IQ) = IQ
          JQCALC(IQ) = IQ
       END DO
       NN = 0
       CALL IOINPUT('JIJSITEI        ',UIO,0,7,NN)
       IF ( NN.EQ.0 ) THEN
          READ (UNIT=UIO,FMT=*) NIQCALC,(IQCALC(IQ),IQ=1,NIQCALC)
          IF ( NIQCALC.LE.0 ) RETURN
          DO IQ = 1,NIQCALC
             IF (( IQCALC(IQ).LE.0 ).OR.( IQCALC(IQ).GT.NAEZ )) RETURN
          END DO
       END IF
       CALL IOINPUT('JIJSITEJ        ',UIO,0,7,NN)
       IF ( NN.EQ.0 ) THEN
          READ (UNIT=UIO,FMT=*) NJQCALC,(JQCALC(IQ),IQ=1,NJQCALC)
          IF ( NJQCALC.LE.0 ) RETURN
          DO IQ = 1,NJQCALC
             IF (( JQCALC(IQ).LE.0 ).OR.( JQCALC(IQ).GT.NAEZ )) RETURN
          END DO
       END IF
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
       IF ( NIQCALC.EQ.NAEZ ) THEN
          WRITE (1337,99005) NAEZ
       ELSE
          WRITE (1337,99006) NIQCALC,NAEZ
       END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ======================================================================
C
C --> determine the size of the arrays IVECI2J,NVECI2J,IREF
C
C ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
       NN = 0
       DO IQS = 1,NIQCALC
          NN = NN + 1
          IQ = IQCALC(IQS)
C **************************************** (selected) sites in unit cell
          DO JQS = 1,NJQCALC
             JQ = JQCALC(JQS)
             DO I = 1,3
                DQ(I) = RBASIS(I,JQ) - RBASIS(I,IQ)
             END DO
C -------------------------------------------------- translation vectors
             DO I1 = -NBR(1),NBR(1)
                DO I2 = -NBR(2),NBR(2)
                   DO I3 = -NBR(3),NBR(3)
                      IBR(1) = I1
                      IBR(2) = I2
                      IBR(3) = I3
                      DO I = 1,3
                          DR(I) = DQ(I)
                          DO J = 1,3
                            DR(I) = DR(I) + DBLE(IBR(J))*BRAVAIS(I,J)
                         END DO
                      END DO
                      DRXYSQ = 0D0
                      DO I = 1,2
                         DRXYSQ = DRXYSQ + DR(I)*DR(I)
                      END DO
                     DRSQ = DR(3)*DR(3)
C
                      IF ( (DRXYSQ.LE.CLURADXYSQ) .AND. 
     &                    (DRSQ.LE.CLURADSQ) )
     &                    NN = NN + 1
                   END DO
                END DO
             END DO
C ----------------------------------------------------------------------
          END DO
C **********************************************************************
       END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CF77----------------------------------------------------------------
Cccc      IF ( NN.GT.NAEZ*NB3MAX ) THEN
Cccc         WRITE (6,99007) 'local','NBMAX',NBMAX + 1
Cccc         STOP
Cccc      END IF
CF77----------------------------------------------------------------
C
CF90----------------------------------------------------------------
       NB3MAX = NN 
       ALLOCATE(IVECI2J(5,NB3MAX),NVECI2J(NB3MAX),
     &     IREF(NB3MAX,NB3MAX),STAT = IQ)
       IF ( IQ.NE.0 ) STOP '    Allocate IVECI2J/NVECI2J/IREF'
CF90----------------------------------------------------------------
       
C
C --> set the first NAEZ vectors (inside the unit-cell at [0,0,0])
C
       NN = 0
       DO IQ = 1,NIQCALC
          NN = NN + 1
          NVECI2J(NN) = 1
          IREF(1,NN) = NN
          IVECI2J(1,NN) = IQCALC(IQ)
          IVECI2J(2,NN) = IQCALC(IQ)
          DO J = 3,5
             IVECI2J(J,NN) = 0
          END DO
       END DO
C ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
       DO IQS = 1,NIQCALC
          IQ = IQCALC(IQS)
C *************************************** (selected)  sites in unit cell
          DO JQS = 1,NJQCALC
             JQ = JQCALC(JQS)
C
C --> set up vector Rij = R_j - r_i = (r_j+T_n) - r_i = (r_j - r_i) + T_n
C                   DR  =                               DQ          + T_n
C
             DO I = 1,3
                DQ(I) = RBASIS(I,JQ) - RBASIS(I,IQ)
             END DO
C ================================================== translation vectors
             DO I1 = -NBR(1),NBR(1)
                DO I2 = -NBR(2),NBR(2)
                   DO I3 = -NBR(3),NBR(3)
                      IBR(1) = I1
                      IBR(2) = I2
                      IBR(3) = I3
                      DO I = 1,3
                         DR(I) = DQ(I)
                         DO J = 1,3
                            DR(I) = DR(I) + DBLE(IBR(J))*BRAVAIS(I,J)
                         END DO
                      END DO
C
C --> calculate Rij(xy)**2 -- DRXYSQ
C     and       Rij(z)**2  -- DRSQ
C
                      DRXYSQ = 0D0
                      DO I = 1,2
                         DRXYSQ = DRXYSQ + DR(I)*DR(I)
                      END DO
                      DRSQ = DR(3)*DR(3)
                      IF ( LSPHER ) DRSQ = DRSQ + DRXYSQ
C
C --> TOL <= Rij(xy)**2 <= CLURADXY**2 and
C     TOL <= Rij(z)**2  <= CLURADZ**2  --> keep the vector Rij by its
C     beginning and end points and the translation indices IBR(1..3)
C
C ------------------------------------------- TOL <= Rij**2 <= CLURAD**2
                      IF ( (DRXYSQ.LE.CLURADXYSQ) .AND. 
     &                    (DRSQ.LE.CLURADSQ) ) THEN
                         NN = NN + 1
                         NVECI2J(NN) = 1
                         IREF(1,NN) = NN
                         IVECI2J(1,NN) = IQ
                         IVECI2J(2,NN) = JQ
                         DO J = 3,5
                            IVECI2J(J,NN) = IBR(J-2)
                         END DO
                      END IF
C ----------------------------------------------------------------------
                   END DO
                END DO
             END DO
C ======================================================================
          END DO
C **********************************************************************
       END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C --> array IVECI2J contains now the positions of all NIQCALC clusters
C     next step eliminates common positions but keeps the information
C     about the connections I,J in IREF(*,*) to be used in setting up
C     IJTABCALC array
C     NOUT is the number of eliminated (repeated) sites
C
       NOUT = 0
C **********************************************************************
      DO I = 1,NN
C ======================================================================
          IF ( NVECI2J(I).EQ.1 ) THEN
C
C --> check vector IVECI2J(I) only if NVECI2J(I).EQ.1
C     finding a J for which the R_j is the same, increment NVECI2J(I)
C     and set NVEVI2(J) to 0
C     same R_j means same site j, same translation vector (indices 2..5)
C
C ----------------------------------------------------------------------
             DO J = I + 1,NN
C
                IF ( NVECI2J(J).EQ.1 ) THEN
                   IEQVEC = 0
                   DO IQ = 2,5
                      IEQVEC = IEQVEC + ABS(IVECI2J(IQ,J)-IVECI2J(IQ,I))
                   END DO
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                   IF ( IEQVEC.EQ.0 ) THEN
                      NVECI2J(J) = 0
                      NVECI2J(I) = NVECI2J(I) + 1
                      IREF(NVECI2J(I),I) = J
                      NOUT = NOUT + 1
                   END IF
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                END IF
C
             END DO
C ----------------------------------------------------------------------
          END IF
C ======================================================================
       END DO
C **********************************************************************
C
C --> get now the actual NATOMIMP cluster positions R_j to be scanned
C     R_j is obtained from IVECI2J
C
       NATOMIMP = 0
C **********************************************************************
       DO I = 1,NN
C ======================================================================
          IF ( NVECI2J(I).NE.0 ) THEN
C
            DO J = 1,3
                DR(J) = RBASIS(J,IVECI2J(2,I))
                DO IQ = 3,5
                   DR(J) = DR(J) + IVECI2J(IQ,I)*BRAVAIS(J,IQ-2)
                END DO
             END DO
C
             NATOMIMP = NATOMIMP + 1
             IF ( NATOMIMP.GT.NATOMIMPD ) THEN
                WRITE (6,99007) 'global','NATOMIMPD',NATOMIMP
                STOP
             END IF
             DO J = 1,3
                RCLSIMP(J,NATOMIMP) = DR(J)
             END DO
C
             ATOMIMP(NATOMIMP) = IVECI2J(2,I)
             NVECI2J(NATOMIMP) = NVECI2J(I)
             DO J = 1,NVECI2J(NATOMIMP)
                IREF(J,NATOMIMP) = IREF(J,I)
             END DO
          END IF
C ======================================================================
       END DO
C **********************************************************************
C --> crosscheck -- if something went wrong return with IDO=0
C
       IF ( (NN-NOUT.NE.NATOMIMP) .OR. (NATOMIMP.LE.0) ) GOTO 100
       IF ( (NAEZ.EQ.1) .AND. (NATOMIMP.EQ.NAEZ) ) GOTO 100
C
C **********************************************************************
C --> set table IJTABCALC(I,J)
c     IJTABCALC(I,J) = 1 for each I = 1,NIQCALC
C                    and for each J which was previously obtained as
C                    connected to I ( IVECI2J(1,IREF)=I )
C
       DO IQ = 1,2*NATOMIMP
          IJTABCALC(IQ) = 0
       END DO
C ======================================================================
       WRITE (1337,99008)
       DO IQ = 1,NIQCALC
          NN = (IQ-1)*NATOMIMP
          NOUT = 0
          DO JQ = 1,NATOMIMP
C ----------------------------------------------------------------------
             IF ( JQ.NE.IQ ) THEN
                DO I = 1,NVECI2J(JQ)
                   IF ( IVECI2J(1,IREF(I,JQ)).EQ.ATOMIMP(IQ))THEN
c                  IF ( IVECI2J(1,IREF(I,JQ)).EQ.1.AND.
c     +               IVECI2J(2,IREF(I,JQ)).EQ.1 ) THEN
                      IJTABCALC(NN+JQ) = 1
                      NOUT = NOUT + 1
c                  END IF
                   END IF
                END DO
             END IF
C ----------------------------------------------------------------------
          END DO
          WRITE (1337,99009) IQ,NOUT
       END DO
       WRITE (1337,99010)
C ======================================================================
       IDO = 1
  100  CONTINUE
CF90----------------------------------------------------------------
       DEALLOCATE(IVECI2J,NVECI2J,IREF,STAT = IQ)
       IF ( IQ.NE.0 ) STOP '    Deallocate IVECI2J/NVECI2J/IREF'
       DEALLOCATE(JQCALC,STAT=IQ)
       IF ( IQ.NE.0 ) STOP '    Deallocate JQCALC'
CF90----------------------------------------------------------------
C ..
99000 FORMAT (5X,'< GIJXCPL > : Exchange coupling constants calculation'
     &        ,/)
99001   FORMAT (6X,
     &       'WARNING: Calculation range JIJRAD missing from your input'
     &       ,/,6X,'         Default value JIJRAD = 0.0 will be assumed'
     &       ,/)
99002 FORMAT (6X,'Range of calculating Jij around each atom',/,
     &        6X,'      spherical cluster of radius :',F7.4,' (ALAT)'//,
     &        6X,'Sites j sought within a parallelipiped',/,6X,
     &        'of size  (',I3,'*a_1) X (',I3,'*a_2) X (',I3,'*a_3)')
99003 FORMAT (6X,'Range of calculating Jij around each atom',/,
     &        6X,'    cylindrical cluster of radius :',F7.4,' (ALAT)'/,
     &        6X,'                           height :',F7.4,' (ALAT)'//,
     &        6X,'Sites j sought within a parallelipiped',/,6X,
     &        'of size  (',I3,'*a_1) X (',I3,'*a_2) X (',I3,'*a_3)')
99004  FORMAT (6X,'Calculations restricted within the unit cell')
99005  FORMAT (6X,' - all of the',I3,' atoms of the u.c. ',
     &     'will be taken into account',/)
99006 FORMAT (6X,' - only',I3,' atom(s) (out of',I3,') in the u.c.',
     &     'will be calculated',/)
99007 FORMAT (6X,'Dimension ERROR: please increase the ',A,' parameter',
     &        /,6X,A,' to a value >=',I5,/)
99008 FORMAT (8X,'Jij connections set',/,10X,20('-'),/,11X,' I ',3X,
     &        'no. of J''s',/,10X,20('-'))
99009 FORMAT (10X,I4,8X,I6)
99010 FORMAT (10X,20('-'),/)
      END
