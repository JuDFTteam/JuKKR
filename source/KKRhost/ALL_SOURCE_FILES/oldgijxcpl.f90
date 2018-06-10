SUBROUTINE gijxcpl(ido,naez,rbasis,bravais,linterface,  &
    niqcalc,iqcalc,natomimp,rclsimp,atomimp, ijtabcalc,natomimpd)
! **********************************************************************
! *                                                                    *
! * In case of tasks requiring Gij blocks calculation, set variables:  *
! *                                                                    *
! * NATOMIMP, RCLSIMP(3,1..NATOMIMP), ATOMIMP(1..NATOMIMP)             *
! * IJTABCALC flag to which pair is needed: I,J --> (I-1)*NATOMIMP + J *
! *           indexing refers to the generated cluster                 *
! * NIQCALC   number of sites in the first unit cell contained in the  *
! *           cluster                                                  *
! * IQCALC()  correspondence index to the 1..NAEZ sites                *
! * IDO takes on the value 1 or 0 if setting up process was OK or not  *
! *                                                                    *
! * EXCHANGE COUPLING CONSTANTS calculation case                       *
! *                                                                    *
! **********************************************************************
IMPLICIT NONE

! Arguments
INTEGER IDO,NAEZ,NATOMIMP,NIQCALC,NATOMIMPD
INTEGER ATOMIMP(*),IJTABCALC(*),IQCALC(*)
DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),RCLSIMP(3,*)
LOGICAL LINTERFACE

! Locals
INTEGER I,I1,I2,I3,IBR(3),IEQVEC,IQ,IQS,J,JQ,MM,NBR(3),NDIM, &
        NN,NOUT
!.......................................................................
!     uniquely identify a vector R_j - r_i = (r_j + T_n) - r_i by a
!     quintuple integer in IVECI2J(5,*):
!     index  meaning
!        1    unit-cell site of atom i (always i - first unit cell)
!        2    unit-cell site of atom j
!        3..5 the three coeficients of the translation vector
!             T_n = n1*a1 + n2*a2 + n3*a3
!     NVECI2J(I) gives the number of identical vectors Rj (initially 1)
!                in the array IVECI2J(2/3/4/5,I)
!     IREF(1..NVECI2J(I),I) points to the NVECI2J(I) identical vectors
!     in the array IVECI2J (one site is kept only once)
!.......................................................................
      INTEGER NB3MAX
      INTEGER IVECI2J(:,:),NVECI2J(:),IREF(:,:)
      ALLOCATABLE IVECI2J,NVECI2J,IREF
      DOUBLE PRECISION CLURAD,CLURADSQ,DQ(3),DR(3),DRSQ,TOL,TOLSQ
      DOUBLE PRECISION CLURADXY,CLURADXYSQ,DRXYSQ
      LOGICAL LSPHER
      CHARACTER*256 UIO !NCOLIO=256
!..
!.. Externals
      EXTERNAL GETCLUSNXYZ,IOINPUT
!     ..
ido = 0

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (6,99000)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

tol = 1.0D-4
tolsq = tol*tol
ndim = 3
IF ( linterface ) ndim = 2

iq = 0
clurad = 0D0
lspher = .true.
CALL ioinput('JIJRAD          ',uio,0,7,iq)
IF ( iq == 0 ) THEN
  READ (UNIT=uio,FMT=*) clurad
  IF ( clurad < 0D0 ) clurad = 0D0
  cluradxy = clurad
  IF ( clurad > 0D0 ) THEN
    iq = 0
    CALL ioinput('JIJRADXY        ',uio,0,7,iq)
    IF ( iq == 0 ) THEN
      READ (UNIT=uio,FMT=*) cluradxy
      IF ( cluradxy <= 0D0 ) cluradxy = clurad
    END IF
  END IF
  lspher = ( ABS(clurad-cluradxy) < tol )
ELSE
  WRITE (6,99001)
END IF

DO i = 1,3
  nbr(i) = 0
END DO
cluradxysq = MAX(clurad,cluradxy)
CALL getclusnxyz(cluradxysq,bravais,ndim,cluradsq,nbr)

cluradsq = (clurad+tol)*(clurad+tol)
cluradxysq = (cluradxy+tol)*(cluradxy+tol)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
IF ( clurad > 0D0 ) THEN
  IF ( lspher ) THEN
    WRITE (6,99002) clurad,((2*nbr(i)+1),i=1,3)
  ELSE
    WRITE (6,99003) cluradxy,clurad,((2*nbr(i)+1),i=1,3)
  END IF
ELSE
  WRITE (6,99004)
END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! --> set the reference I sites for J_IJ within the unit cell
!     token JIJSITES allows one to specify which sites in the
!     unit cell should be considered in the calculation (default ALL)

niqcalc = naez
DO iq = 1,naez
  iqcalc(iq) = iq
END DO
nn = 0
CALL ioinput('JIJSITES        ',uio,0,7,nn)
IF ( nn == 0 ) THEN
  READ (UNIT=uio,FMT=*) niqcalc,(iqcalc(iq),iq=1,niqcalc)
  IF ( niqcalc <= 0 ) RETURN
  DO iq = 1,niqcalc
    IF (( iqcalc(iq) <= 0 ).OR.( iqcalc(iq) > naez )) RETURN
  END DO
END IF

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
IF ( niqcalc == naez ) THEN
  WRITE (6,99005) naez
ELSE
  WRITE (6,99006) niqcalc,naez
END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ======================================================================

! --> determine the size of the arrays IVECI2J,NVECI2J,IREF

! ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
nn = 0
DO iqs = 1,niqcalc
  nn = nn + 1
  iq = iqcalc(iqs)
! *********************************************** all sites in unit cell
  DO jq = 1,1            !NAEZ  change Fivos
    DO i = 1,3
      dq(i) = rbasis(i,jq) - rbasis(i,iq)
    END DO
! -------------------------------------------------- translation vectors
    DO i1 = -nbr(1),nbr(1)
      DO i2 = -nbr(2),nbr(2)
        DO i3 = -nbr(3),nbr(3)
          ibr(1) = i1
          ibr(2) = i2
          ibr(3) = i3
          DO i = 1,3
            dr(i) = dq(i)
            DO j = 1,3
              dr(i) = dr(i) + DBLE(ibr(j))*bravais(i,j)
            END DO
          END DO
          drxysq = 0D0
          DO i = 1,2
            drxysq = drxysq + dr(i)*dr(i)
          END DO
          drsq = dr(3)*dr(3)
          
          IF ( (drxysq <= cluradxysq) .AND. (drsq <= cluradsq) )  &
              nn = nn + 1
        END DO
      END DO
    END DO
! ----------------------------------------------------------------------
  END DO
! **********************************************************************
END DO
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

nb3max = nn
allocate(iveci2j(5,nb3max),nveci2j(nb3max), iref(nb3max,nb3max),stat = iq)
IF ( iq /= 0 ) STOP '    Allocate IVECI2J/NVECI2J/IREF'


! --> set the first NAEZ vectors (inside the unit-cell at [0,0,0])

nn = 0
DO iq = 1,niqcalc
  nn = nn + 1
  nveci2j(nn) = 1
  iref(1,nn) = nn
  iveci2j(1,nn) = iqcalc(iq)
  iveci2j(2,nn) = iqcalc(iq)
  DO j = 3,5
    iveci2j(j,nn) = 0
  END DO
END DO
! ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
DO iqs = 1,niqcalc
  iq = iqcalc(iqs)
! *********************************************** all sites in unit cell
  DO jq = 1,1           !NAEZ  change Fivos
    
! --> set up vector Rij = R_j - r_i = (r_j+T_n) - r_i = (r_j - r_i) + T_n
!                   DR  =                               DQ          + T_n
    
    DO i = 1,3
      dq(i) = rbasis(i,jq) - rbasis(i,iq)
    END DO
! ================================================== translation vectors
    DO i1 = -nbr(1),nbr(1)
      DO i2 = -nbr(2),nbr(2)
        DO i3 = -nbr(3),nbr(3)
          ibr(1) = i1
          ibr(2) = i2
          ibr(3) = i3
          DO i = 1,3
            dr(i) = dq(i)
            DO j = 1,3
              dr(i) = dr(i) + DBLE(ibr(j))*bravais(i,j)
            END DO
          END DO
          
! --> calculate Rij(xy)**2 -- DRXYSQ
!     and       Rij(z)**2  -- DRSQ
          
          drxysq = 0D0
          DO i = 1,2
            drxysq = drxysq + dr(i)*dr(i)
          END DO
          drsq = dr(3)*dr(3)
          IF ( lspher ) drsq = drsq + drxysq
          
! --> TOL <= Rij(xy)**2 <= CLURADXY**2 and
!     TOL <= Rij(z)**2  <= CLURADZ**2  --> keep the vector Rij by its
!     beginning and end points and the translation indices IBR(1..3)
          
! ------------------------------------------- TOL <= Rij**2 <= CLURAD**2
          IF ( (drxysq <= cluradxysq) .AND. (drsq <= cluradsq) ) THEN
            nn = nn + 1
            nveci2j(nn) = 1
            iref(1,nn) = nn
            iveci2j(1,nn) = iq
            iveci2j(2,nn) = jq
            DO j = 3,5
              iveci2j(j,nn) = ibr(j-2)
            END DO
          END IF
! ----------------------------------------------------------------------
        END DO
      END DO
    END DO
! ======================================================================
  END DO
! **********************************************************************
END DO
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! --> array IVECI2J contains now the positions of all NIQCALC clusters
!     next step eliminates common positions but keeps the information
!     about the connections I,J in IREF(*,*) to be used in setting up
!     IJTABCALC array
!     NOUT is the number of eliminated (repeated) sites

nout = 0
! **********************************************************************
DO i = 1,nn
! ======================================================================
  IF ( nveci2j(i) == 1 ) THEN
    
! --> check vector IVECI2J(I) only if NVECI2J(I).EQ.1
!     finding a J for which the R_j is the same, increment NVECI2J(I)
!     and set NVEVI2(J) to 0
!     same R_j means same site j, same translation vector (indices 2..5)
    
! ----------------------------------------------------------------------
    DO j = i + 1,nn
      
      IF ( nveci2j(j) == 1 ) THEN
        ieqvec = 0
        DO iq = 2,5
          ieqvec = ieqvec + ABS(iveci2j(iq,j)-iveci2j(iq,i))
        END DO
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        IF ( ieqvec == 0 ) THEN
          nveci2j(j) = 0
          nveci2j(i) = nveci2j(i) + 1
          iref(nveci2j(i),i) = j
          nout = nout + 1
        END IF
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      END IF
      
    END DO
! ----------------------------------------------------------------------
  END IF
! ======================================================================
END DO
! **********************************************************************

! --> get now the actual NATOMIMP cluster positions R_j to be scanned
!     R_j is obtained from IVECI2J

natomimp = 0
! **********************************************************************
DO i = 1,nn
! ======================================================================
  IF ( nveci2j(i) /= 0 ) THEN
    
    DO j = 1,3
      dr(j) = rbasis(j,iveci2j(2,i))
      DO iq = 3,5
        dr(j) = dr(j) + iveci2j(iq,i)*bravais(j,iq-2)
      END DO
    END DO
    
    natomimp = natomimp + 1
    IF ( natomimp > natomimpd ) THEN
      WRITE (6,99007) 'global','NATOMIMPD',natomimp
      STOP
    END IF
    DO j = 1,3
      rclsimp(j,natomimp) = dr(j)
    END DO
    
    atomimp(natomimp) = iveci2j(2,i)
    nveci2j(natomimp) = nveci2j(i)
    DO j = 1,nveci2j(natomimp)
      iref(j,natomimp) = iref(j,i)
    END DO
  END IF
! ======================================================================
END DO
! **********************************************************************
! --> crosscheck -- if something went wrong return with IDO=0

IF ( (nn-nout /= natomimp) .OR. (natomimp <= 0) ) GO TO 100
IF ( (naez == 1) .AND. (natomimp == naez) ) GO TO 100

! **********************************************************************
! --> set table IJTABCALC(I,J)
!     IJTABCALC(I,J) = 1 for each I = 1,NIQCALC
!                    and for each J which was previously obtained as
!                    connected to I ( IVECI2J(1,IREF)=I )

DO iq = 1,2*natomimp
  ijtabcalc(iq) = 0
END DO
! ======================================================================
WRITE (6,99008)
DO iq = 1,niqcalc
  nn = (iq-1)*natomimp
  nout = 0
  DO jq = 1,natomimp
! ----------------------------------------------------------------------
    IF ( jq /= iq ) THEN
      DO i = 1,nveci2j(jq)
        IF ( iveci2j(1,iref(i,jq)) == atomimp(iq))THEN
!                  IF ( IVECI2J(1,IREF(I,JQ)).EQ.1.AND.
!     +               IVECI2J(2,IREF(I,JQ)).EQ.1 ) THEN
          ijtabcalc(nn+jq) = 1
          nout = nout + 1
!                  END IF
        END IF
      END DO
    END IF
! ----------------------------------------------------------------------
  END DO
  WRITE (6,99009) iq,nout
END DO
WRITE (6,99010)
! ======================================================================
ido = 1
100  CONTINUE
deallocate(iveci2j,nveci2j,iref,stat = iq)
IF ( iq /= 0 ) STOP '    Deallocate IVECI2J/NVECI2J/IREF'
! ..
99000 FORMAT (5X,'< GIJXCPL > : Exchange coupling constants calculation' ,/)
99001 FORMAT (6X,  &
    'WARNING: Calculation range JIJRAD missing from your input'  &
    ,/,6X,'         Default value JIJRAD = 0.0 will be assumed' ,/)
99002 FORMAT (6X,'Range of calculating Jij around each atom',/,  &
    6X,'      spherical cluster of radius :',f7.4,' (ALAT)'//,  &
    6X,'Sites j sought within a parallelipiped',/,6X,  &
    'of size  (',i3,'*a_1) X (',i3,'*a_2) X (',i3,'*a_3)')
99003 FORMAT (6X,'Range of calculating Jij around each atom',/,  &
    6X,'    cylindrical cluster of radius :',f7.4,' (ALAT)'/,  &
    6X,'                           height :',f7.4,' (ALAT)'//,  &
    6X,'Sites j sought within a parallelipiped',/,6X,  &
    'of size  (',i3,'*a_1) X (',i3,'*a_2) X (',i3,'*a_3)')
99004 FORMAT (6X,'Calculations restricted within the unit cell')
99005 FORMAT (6X,' - all of the',i3,' atoms of the u.c. ',  &
    'will be taken into account',/)
99006 FORMAT (6X,' - only',i3,' atom(s) (out of',i3,') in the u.c.',  &
    'will be calculated',/)
99007 FORMAT (6X,'Dimension ERROR: please increase the ',a,' parameter',  &
    /,6X,a,' to a value >=',i5,/)
99008 FORMAT (8X,'Jij connections set',/,10X,20('-'),/,11X,' I ',3X,  &
    'no. of J''S',/,10X,20('-'))
99009 FORMAT (10X,i4,8X,i6)
99010 FORMAT (10X,20('-'),/)
END SUBROUTINE gijxcpl
