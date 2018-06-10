! ************************************************************************
SUBROUTINE clsgen_tb(naez,nemb,nvirt,rr,nr,rbasis,kaoez,zat,  &
    cls,ncls,nacls,atom,ezoa, nlbasis,nrbasis,nleft,nright,zperleft,zperight,  &
    tleft,tright,rmtref,rmtrefat,vref,  &
    irefpot,nrefpot,rcls,rcut,rcutxy,l2dim,alat,  &
    naezd,natyp,nembd,nrd,naclsd,nclsd,nrefd)
! ************************************************************************
! This subroutine is used to create the clusters around each atom
! (Based on clsgen99.f). Also the reference potential height and radius is set
! (vref and rmtref).
!
! STRATEGY :
! Calculate the cluster of each atom by the lattice
! parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
! compare the positions with the previous clusters to see if there is
! a difference. If not keep only previous clusters and make indexing if
! a new cluster is found then check dimensions and continue for the new
! atom.
!
!
!
      use mod_version_info
implicit none
!.. arguments
INTEGER NAEZD,NATYP,NEMBD,NRD,NACLSD,NCLSD,NREFD
DOUBLE PRECISION       ALAT         ! lattice constant A
DOUBLE PRECISION       RCUT,RCUTXY
DOUBLE PRECISION  RBASIS(3,naez+nembd),     &        ! pos. of basis atoms in EZ
     RCLS(3,NACLSD,ncls),&        ! real space position of atom in cluster
     RR(3,0:NR),     &        ! set of lattice vectors
     ZAT(natyp),          &        ! nucleus charge
     RMTREF(nrefd)
!     RWS(*),
!     BBOX(3)                  ! bounding box for povray plots

INTEGER NAEZ,         &           ! number of atoms in EZ
     NEMB,         &           ! number of embedding postions
     NCLS,         &           ! number of diff. clusters
     NR,           &           ! number of lattice vectors RR
     NLR,          &           ! =NEMB in decimation, =0 in slab or bulk
     NVIRT,        &           ! Number of virtual atoms
     NPRINC                   ! Calculated number of layers in a principal layer

INTEGER CLS(naez+nembd),                   & ! type of cluster around atom
     KAOEZ(NATYP,NAEZ+NEMBD),& ! type of atom at position in EZ
     NACLS(natyp),                 & ! number of atoms in cluster
     ATOM(NACLSD,naez+nembd),           & ! index to atom in elem/cell at site in cluster
     EZOA (NACLSD,naez+nembd)             ! index to bravais lattice  at site in cluster

!.. locals
INTEGER ILAY,N1,IR,ISITE,JSITE,IAT1, &
     NA,NUMBER,MAXNUMBER, & !IX,
     POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,I1,ICLUSTER,NCLSALL
INTEGER IATOM(NACLSD),IEZOA(NACLSD), &
     ISORT(NACLSD),ICOUPLMAT(NAEZ,NAEZ), &
     IREP(NCLS) ! representative atom of cluster (inverse of CLS)
INTEGER IREFPOT(NAEZ+NEMBD),NREFPOT
DOUBLE PRECISION RMTREFAT(NAEZ+NEMBD),RMTREF1(NAEZ+NEMBD)
DOUBLE PRECISION VREFAT(NAEZ+NEMBD),VREF1(NAEZ+NEMBD), &
                 VREF(NREFD)


DOUBLE PRECISION R2,EPSSHL,TOL,TOL2,DISTMIN,&
     RCLS1(3,NACLSD),&
     RG(3,NACLSD),TMP(3),RSORT(NACLSD)
INTEGER NLBASIS,NRBASIS,                    &
        NLEFT,NRIGHT           
DOUBLE PRECISION ZPERLEFT(3),ZPERIGHT(3),            &
        TLEFT(3,nlbasis),TRIGHT(3,nrbasis)
DOUBLE PRECISION       RCUT2,RCUTXY2,RXY2,DIST

LOGICAL LFOUND
LOGICAL  L2DIM,CLUSTCOMP_TB

EXTERNAL DSORT,CLUSTCOMP_TB
INTRINSIC MIN,SQRT

DATA     EPSSHL   / 1.D-4 /
DATA     TOL   / 1.D-7 /
DATA     TOL2   / 1.D-7 /

! ------------------------------------------------------------------------
WRITE(1337,*) '>>> CLSGEN_TB: generation of cluster coordinates'
! This is generating the clusters which have a distance smaller
! than RCUT and RCUTXY in plane .
! The cluster atoms are ordered with radius and then z>y>x
! The ordering allows an easy comparison of clusters
! The principal layer for each layer (atom in unit cell) is
! calculated also for each cluster and the maximum number
! is returned. Some dimension tests are also done
WRITE(1337,*) 'RCUT = ',rcut,' RCUTXY = ',rcutxy
IF (ABS(rcutxy - rcut) < 1.d-4) THEN
  WRITE(1337,*) 'Spherical Clusters are created'
!          LSPHER=.TRUE.
END IF
OPEN(8,FILE='clusters',STATUS='unknown')
CALL version_print_header(8)
WRITE(8,9005) naez
WRITE(8,9030) alat
WRITE(8,9010) (zat(kaoez(1,iat)),iat=1,naez-nvirt)
WRITE(8,9020) (kaoez(1,iat),iat=1,naez)

rcutxy2 = rcutxy**2
rcut2 = rcut**2
nlr = 0
IF (l2dim) nlr = nemb
vrefat(:) = 8.d0      ! Set to 8 Rydbergs
vref1(:) = 8.d0      ! Set to 8 Rydbergs
vref(:) = 8.d0      ! Set to 8 Rydbergs

!===============================================================
! Check the if the dimension NACLSD is enough before starting to
! assign clusters. Write out necessary dimension.
! Also find touching MT radius.
DO  jatom = 1,naez + nlr ! loop in all sites incl. left/right host in decimation case
  
  maxnumber = 0
  NUMBER = 0             ! counter for atoms in cluster
  distmin = 1.d100       ! Large initial value for RMT**2
  DO na = 1,naez         ! loop in all sites in unit cell
    IF (kaoez(1,na) /= -1) THEN ! Exclude virtual atoms from clusters (except clust. center)
      DO ir = 0,nr     ! loop in all bravais vectors
        tmp(1:3) = rr(1:3,ir)+rbasis(1:3,na)-rbasis(1:3,jatom)
        rxy2 =  tmp(1)**2+tmp(2)**2
        r2   =  tmp(3)**2 + tmp(1)**2+tmp(2)**2
        IF ( (rxy2 <= rcutxy2).AND.(r2 <= rcut2) )  THEN
          NUMBER = NUMBER + 1
        END IF
        IF (r2 > tol2) distmin = MIN(distmin,r2)
      END DO           ! IR loop in bravais
    END IF               ! (KAOEZ(1,NA).NE.-1)
  END DO                 ! NA loop in NAEZ
  
!
!        In the case of 2 dimensional case loop in the atoms outside.
!
  IF (l2dim) THEN
!        Somehow messy (lionel messi?)
    DO ir = 0,nr
      DO ilay = nleft,1,-1  ! loop in some layers on left side
        DO i1 = nlbasis,1,-1 ! loop in representative atoms on left side
          tmp(1:3) = rr(1:3,ir) + tleft(1:3,i1) +  &
              (ilay-1)*zperleft(1:3) - rbasis(1:3,jatom)
          rxy2 =  tmp(1)**2 + tmp(2)**2
          r2  =   tmp(3)**2 + rxy2
          
          IF ((rxy2 <= rcutxy2).AND.(r2 <= rcut2)) THEN
            NUMBER = NUMBER + 1
          END IF
          IF (r2 > tol2) distmin = MIN(distmin,r2)
        END DO
      END DO
      
      DO ilay = 1,nright
        DO i1 = 1,nrbasis
          tmp(1:3) = rr(1:3,ir)+ tright(1:3,i1) +  &
              (ilay-1)*zperight(1:3) - rbasis(1:3,jatom)
          rxy2 = tmp(1)**2 + tmp(2)**2
          r2  =  tmp(3)**2 + rxy2
          IF ((rxy2 <= rcutxy2).AND.(r2 <= rcut2)) THEN
            NUMBER = NUMBER + 1
          END IF
          IF (r2 > tol2) distmin = MIN(distmin,r2)
        END DO
      END DO
!
      
    END DO              ! loop in all bravais lattices
  END IF                 ! L2DIM Interface calculation
  
  distmin = DSQRT(distmin)/2.d0 ! Touching RMT
! Define MT-radius of TB-ref. potentials if not read-in from the input
  IF (rmtrefat(jatom) < 0.d0) & ! I.e., if not read in  &
      rmtrefat(jatom) = INT( distmin*alat*99.d0 )/100.d0 ! round up the decimals
  IF (kaoez(1,jatom) == -1) THEN ! Virtual atom
    rmtrefat(jatom) = 1.d-20
    vrefat(jatom) = 0.d0
  END IF
  
!
  maxnumber = MAX(maxnumber,NUMBER) ! Find largest cluster size
  
  WRITE(1337,*) 'clsgen_tb: cluster size of site:',jatom,':', NUMBER
  WRITE(1337,*) 'clsgen_tb: Touching RMT of site:',jatom,':', distmin
  
  
END DO

IF (NUMBER > naclsd) THEN
  WRITE(6,*) '(a) Increase the parameter NACLSD ',  &
      'to a value greater equal ',maxnumber,'.'
  STOP 'clsgen_tb: Dimension error (a).'
END IF

!===============================================================

! Find different types of ref. potential according to the rmtref and vref
irefpot(1) = 1
rmtref1(1) = rmtrefat(1)
vref1(1) = vrefat(1)
nrefpot = 1
DO isite = 2,naez + nemb
  lfound = .false.
  DO jsite = 1,isite - 1
    IF ( ABS(rmtrefat(isite) - rmtrefat(jsite)) +  &
          ABS(vrefat(isite)   - vrefat(jsite)      ) <= tol) THEN
      irefpot(isite) = irefpot(jsite)
      lfound = .true.
    END IF
  END DO
  IF (.NOT.lfound) THEN
    nrefpot = nrefpot + 1
    irefpot(isite) = nrefpot
! RMTREFAT goes over all sites, RMTREF1 only over all inequivalent ref. potentials.
    rmtref1(nrefpot) = rmtrefat(isite)
    vref1(nrefpot) = vrefat(isite)
  END IF
END DO

IF (nrefpot > nrefd) THEN
  WRITE(*,*) 'clsgen_tb: NREFPOT.GT.NREFD:',nrefpot,nrefd
  STOP 'clsgen_tb: NREFPOT.GT.NREFD'
END IF
! Now that the dimension is known, copy to array RMTREF
DO i1 = 1,nrefpot
  IF (rmtref(i1) < 0.d0) rmtref(i1) = rmtref1(i1)
END DO
vref(1:nrefpot) = vref1(1:nrefpot)



icluster = 1
rcutxy2 = (rcutxy+epsshl)*(rcutxy+epsshl)
rcut2   = (rcut+epsshl)*(rcut+epsshl)



!===============================================================

DO  jatom = 1,naez + nlr   ! loop in all atoms or layers
  
  cls(jatom) = 0
  
  NUMBER = 0                      ! counter for sites in cluster
  DO na = 1,naez                  ! loop in all sites
    IF (kaoez(1,na) /= -1) THEN  ! proceed only if the neighbour is not virtual atom
      DO ir = 0,nr              ! loop in all bravais vectors
        tmp(1:3) = rr(1:3,ir)+rbasis(1:3,na)-rbasis(1:3,jatom)
        rxy2 =  tmp(1)**2+tmp(2)**2
        r2   =  tmp(3)**2 + tmp(1)**2+tmp(2)**2
        
        IF ( (rxy2 <= rcutxy2).AND.(r2 <= rcut2) )  THEN
          NUMBER = NUMBER + 1
          atom(NUMBER,jatom) = na ! store the atom in elem cell
          ezoa(NUMBER,jatom) = ir ! store the bravais vector
          rcls1(1:3,NUMBER) = tmp(1:3)
        END IF
      END DO           ! IR loop in bravais
    END IF
  END DO                 ! NA loop in NAEZ
  
!
!        In the case of 2 dimensional case loop in the atoms outside.
!
  IF (l2dim) THEN
!        Somehow messy (eh? lionel?)
!        Index ATOM gives the kind of atom
!
    DO ir = 0,nr
      DO ilay = nleft,1,-1  ! loop in some layers on left side
        DO i1 = nlbasis,1,-1 ! loop in representative atoms on left side
          tmp(1:3) = rr(1:3,ir) + tleft(1:3,i1) +  &
              (ilay-1)*zperleft(1:3) - rbasis(1:3,jatom)
          rxy2 =  tmp(1)**2 + tmp(2)**2
          r2  =   tmp(3)**2 + rxy2
          
          IF ((rxy2 <= rcutxy2).AND.(r2 <= rcut2)) THEN
            NUMBER = NUMBER + 1
            atom(NUMBER,jatom) = -naez - i1 ! negative values are used in dlke1.f
            ezoa(NUMBER,jatom) = ir ! ILAY,I1 are negative
            rcls1(1:3,NUMBER) = tmp(1:3)
          END IF
        END DO
      END DO
!
!
      DO ilay = 1,nright
        DO i1 = 1,nrbasis
          tmp(1:3) = rr(1:3,ir)+ tright(1:3,i1) +  &
              (ilay-1)*zperight(1:3) - rbasis(1:3,jatom)
          rxy2 = tmp(1)**2 + tmp(2)**2
          r2  =  tmp(3)**2 + rxy2
          IF ((rxy2 <= rcutxy2).AND.(r2 <= rcut2)) THEN
            NUMBER = NUMBER + 1
            atom(NUMBER,jatom) = -naez - nlbasis - i1
            ezoa(NUMBER,jatom) = ir
            rcls1(1:3,NUMBER) = tmp(1:3)
          END IF
        END DO
      END DO
!
      
    END DO              ! loop in all bravais lattices
  END IF                 ! L2DIM Interface calculation
!
!     Now the atom JATOM has its cluster.
!     Sort the atoms of the cluster in increasing order.
!     First by distance, then by z, then by y, then by x.
!
  IF (NUMBER > naclsd) THEN ! should not hit here, this was checked earlier
    WRITE(6,*) '(b) Increase the parameter NACLSD ',  &
        'to a value greater equal ',NUMBER,'.'
    STOP 'clsgen_tb: Dimension error (b).'
  END IF
  
  DO ia=1,NUMBER
    rsort(ia) = DSQRT(rcls1(1,ia)**2+ rcls1(2,ia)**2+  &
        rcls1(3,ia)**2)
    
    rsort(ia) = 1.d9*rsort(ia)+ 1.d6*rcls1(3,ia)+  &
        1.d3*rcls1(2,ia)+ 1.d0*rcls1(1,ia)
  END DO
!
  CALL dsort(rsort,isort,NUMBER,pos)
!     Rearange exchange ia with ib
!     MAP temporarily to another array
  DO ia=1,NUMBER
    rg(1:3,ia) = rcls1(1:3,ia)
    iatom(ia) = atom(ia,jatom)
    iezoa(ia) = ezoa(ia,jatom)
  END DO
! Now use correct order
  DO ia =1,NUMBER
    ib = isort(ia)
    rcls1(1:3,ia) = rg(1:3,ib)
    atom(ia,jatom) = iatom(ib)
    ezoa(ia,jatom) = iezoa(ib)
  END DO
!
!     Now the clusters have a unique sorting and can be compared with
!     each other Check if ICLUSTER was found previously
!
  DO icu = 1,icluster-1
    n1 = nacls(icu)
    iat1 = irep(icu)
    IF ( clustcomp_tb( rcls,irefpot,atom,iat1,  &
        icu,n1,rcls1,NUMBER,jatom,naclsd) ) cls(jatom) = icu
! return true if found before
  END DO
  
  IF (cls(jatom) == 0) THEN ! no equivalent found, add new cluster
    nclsall = icluster ! incl. embedded atoms of left/right host
    IF (jatom <= naez) ncls = icluster  ! Excl. embedded atoms of left/right host
    IF (nclsall > nclsd) THEN
      WRITE(6,*) '(c) Increase the parameter NCLSD ',  &
          '  to a value greater equal ',icluster,' .'
      STOP 'clsgen_tb: Dimension error (c).'
    END IF
    cls(jatom) = icluster
    nacls(icluster) = NUMBER
    irep(icluster) = jatom ! cluster-class is represented by the cluster around jatom
    DO in = 1,NUMBER
      rcls(1:3,in,icluster) = rcls1(1:3,in)
!               WRITE(6,800) JATOM,ATOM(IN,JATOM),EZOA(IN,JATOM),
!     &                  (RCLS1(IX,IN),IX=1,3),
!     &              SQRT(RCLS1(1,IN)**2+RCLS1(2,IN)**2+RCLS1(3,IN)**2)
! 800           FORMAT(3I5,4F8.4)
    END DO
    icluster = icluster + 1
  END IF
! ******************************************************
  enddo                     ! JATOM = 1,NAEZ + NEMB

! Now all clusters of all atoms are found.


WRITE(1337,*) 'Clusters from clsgen_tb:'
DO jatom = 1,naez
  WRITE(1337,8000) jatom,irefpot(jatom),rmtrefat(jatom),  &
      vrefat(jatom),cls(jatom),nacls(cls(jatom))
END DO
IF (l2dim) THEN
  WRITE(1337,*) 'Clusters from clsgen_tb in outer region, left:'
  DO ia = 1,nlbasis
    jatom = naez + ia
    WRITE(1337,8000) jatom,irefpot(jatom),rmtrefat(jatom),  &
        vrefat(jatom),cls(jatom),nacls(cls(jatom))
  END DO
  
  WRITE(1337,*) 'Clusters from clsgen_tb in outer region, right:'
  DO ia = 1,nrbasis
    jatom = naez + nlbasis + ia
    WRITE(1337,8000) jatom,irefpot(jatom),rmtrefat(jatom),  &
        vrefat(jatom),cls(jatom),nacls(cls(jatom))
  END DO
END IF

! Write out clusters in file
OPEN(8,FILE='clusters',STATUS='UNKNOWN')
WRITE(8,9005) naez
WRITE(8,9030) alat
WRITE(8,9010) (zat(kaoez(1,i1)),i1=1,naez)
WRITE(8,9020) (kaoez(1,i1),i1=1,naez)
DO jatom = 1,naez
  ic = cls(jatom)
  NUMBER = nacls(ic)
  WRITE(8,FMT=1030) NUMBER
  WRITE(8,FMT=1030) jatom,ic
  DO i1=1,NUMBER
    dist = SQRT( rcls(1,i1,ic)**2+rcls(2,i1,ic)**2+rcls(3,i1,ic)**2)
    WRITE(8,1041) (rcls(ii,i1,ic),ii=1,3),atom(i1,jatom),  &
        zat(ABS(atom(i1,jatom))),dist
  END DO
  
END DO

! Write out the coupling matrix
WRITE(1337,*) 'Coupling matrix:'
DO jatom = 1,naez
  DO iat = 1,naez
    icouplmat(jatom,iat) = 0
    DO i1=1,NUMBER
      IF (atom(i1,jatom) == iat) THEN
        icouplmat(jatom,iat) = 1
      END IF
    END DO
  END DO
  WRITE(1337,9060) jatom,(icouplmat(jatom,iat),iat=1,naez)
END DO

IF (l2dim) THEN
! Calculate number of layers in principal layer
  nprinc = 1
  DO jatom = 1,naez  ! loop over rows
    DO iat = 1,jatom - 1 ! loop over columns before the diagonal
      IF (icouplmat(jatom,iat) == 1) nprinc = MAX(nprinc,jatom-iat)
    END DO
    DO iat = jatom + 1,naez ! loop over columns after the diagonal
      IF (icouplmat(jatom,iat) == 1) nprinc = MAX(nprinc,iat-jatom)
    END DO
  END DO
  WRITE(1337,*)  &
      'CLSGEN_TB: Number of layers in a principal layer: NPRINC=', nprinc
END IF



! ------------------------------------------------------------------------
WRITE(1337,*) ' Sub clsgen_tb  exiting <<<<<<<<<<<<<'
! ------------------------------------------------------------------------

1000   FORMAT(' cluster around atom     ',10I4/,  &
    ' number atoms in cluster ',10I4)
1001   FORMAT(i4,2I5,3F8.2,i6,4F7.2)
1002   FORMAT(' cocls : naez =',i3,' (x,y,z)= ',3F10.4)
1010   FORMAT(12X,i6,3F10.4)
1020   FORMAT('  Nr  naez kaoez     x       y       z',  &
    '  ezoa  RR(1)  RR(2)  RR(3)      R')
1030   FORMAT(3I8)
1040   FORMAT('Cu  ',3D24.16,2I8,f18.12)
1041   FORMAT(3E27.19,i5,f5.1,e17.9)
1050   FORMAT(3F12.7,'  scaling factor')
1060   FORMAT(i4,3F12.7,'  center',/,(i4,3F12.7))
1070   FORMAT(i4,3F12.7,'  center of gravity')
1080   FORMAT('contains ',i4,'  atoms.')
8000   FORMAT('CLSGEN_TB: Atom',i5,' Refpot',i5,' Rmtref',f10.7,  &
    ' Vref',f10.7,' TB-cluster',i5,' Sites',i5)
9005   FORMAT(i4)
9010   FORMAT(('# ZAT   ',20F6.2))
9020   FORMAT(('# KAOEZ ',20I4))
9030   FORMAT(f12.7,6X,'ALAT')
9040   FORMAT('> cluster ',i4,' at atom ',i4, ' of type ',i4,'.')
9060   FORMAT(i4,1X,500I1)
9080   FORMAT(i6,3F15.6)
9090   FORMAT('The Number of layers is each Principal Layer = ', i5)


RETURN
END SUBROUTINE clsgen_tb


