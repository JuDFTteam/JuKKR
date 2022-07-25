PROGRAM EMPTYCELLS
implicit none
! This program finds optimal positions for empty spheres in 3D-periodic lattices
! by means of a Monte-Carlo optimization.
! On input it requires the Bravais vectors, number and positions of basis sites 
! of the lattice (NFIX, RFIX), and wished number of empty spheres (NVAR) that it
! positions on output at RVAR.

! Created by Phivos Mavropoulos, 2013-2014.

INTEGER NFIX ! Number of fixed atoms
INTEGER NVAR ! Number of moving atoms with variable positions to be optimized
INTEGER NSUPER ! Supercell extends from -NSUPER to +NSUPER in 3 dimensions
INTEGER NEXPANDVAR,NEXPANDFIX ! Number of moving and fixed atoms in supercell
REAL*8 BRAVAIS(3,3),BRINV(3,3),RECBV(3,3) ! Bravais vectors, inverse Bravais matrix, Bravais of rec. space
LOGICAL CARTESIAN
REAL*8,ALLOCATABLE :: RFIX(:,:),RVAR(:,:) ! Positions of fixed and moving atoms in primitive cell
REAL*8,ALLOCATABLE :: REXPANDFIX(:,:),REXPANDVAR(:,:) ! Positions of fixed and moving atoms in supercell
REAL*8,ALLOCATABLE :: DISTFIX(:),DISTVAR(:) ! Distance of an atom to all other atoms
REAL*8,ALLOCATABLE :: RMTFIX(:),RMTVAR(:) ! Muffin tin radii of moving and fixed atoms
REAL*8,ALLOCATABLE :: RHFIX(:),RHVAR(:) ! Hard-sphere radii of moving and fixed atoms
INTEGER,ALLOCATABLE :: TYPEFIX(:),TYPEVAR(:) ! Type of atom in extended supercell corresponding to prim. cell
REAL*8 ROLD(3),RNEW(3) ! Atom position before and after move
REAL*8 SHIFTMAX,SHIFTSTART,TEMPR,TSTART  ! Max. shift for moving the atom with respect to Bravais vector, temperature
REAL*8 SDIV,TDIV
REAL*8 EEXTI,EINTI,EOLD,EDIFF,XMETR,EXPEDIFF,ACCEPTANCE,EINTRN,EEXTRN,ETOT,EFLUCT,ENEW,EINTRN1,EEXTRN1,ETOT1 ! energies
REAL*8 EINTIOLD,EINTINEW,EINTIDIFF,EEXTIOLD,EEXTINEW,EEXTIDIFF ! energies
REAL*8 VOLMT,VOLUME,VOLFIL ! MT-occupied volume, total volume, and volume-filling fraction
REAL*8 HSVNUM ! Number of hard-sphere violations
REAL*8 RTMP(3) ! Temporary
INTEGER MCSTEPS,ISEED,ITEMPR,ISTEP,NTEMPRMC,METHOD ! MC-parameters
INTEGER N0,NN0
LOGICAL LACCEPT,LREAD
INTEGER IFIX,IVAR,IX,IB,IVAR1,COUNT ! Running indices
CHARACTER*1 HASH
CHARACTER*200 UIO  ! For IoInput routine
INTEGER IER        ! For IoInput routine

REAL*8 XRAND,PI
EXTERNAL XRAND
PI = 4.D0 * DATAN(1.D0)

! Read in and initialize
OPEN(21,FILE='inputcard_constructed')
CALL READPARA(BRAVAIS,CARTESIAN,NFIX,NVAR,ISEED,NSUPER,METHOD,TSTART,NTEMPRMC,TDIV,SHIFTSTART,SDIV,MCSTEPS)
! Calculate inverse of Bravais matrix BRINV(3,3)
CALL INVERTBR(BRAVAIS,BRINV,RECBV)
ALLOCATE( RFIX(3,NFIX) , RVAR(3,NVAR) )
CALL READPOS(NFIX,NVAR,BRAVAIS,BRINV,CARTESIAN,ISEED,RFIX,RVAR)
CLOSE(21)

! Rationalize basis to bring atoms in the primitive unit cell
CALL RATIONALBASIS(BRAVAIS,BRINV,NFIX,RFIX)
! For the moving atoms it is not necessary because they are
! positioned in the primitive cell by construction.
! CALL RATIONALBASIS(BRAVAIS,BRINV,NVAR,RVAR) 


! Allocate fixed positions and fixed supercell
NEXPANDFIX = (2*NSUPER+1)**3*NFIX
ALLOCATE( REXPANDFIX(3,NEXPANDFIX) , TYPEFIX(NEXPANDFIX) , DISTFIX(NFIX) , RMTFIX(NFIX) , RHFIX(NFIX) )
CALL EXPANDUC3D(BRAVAIS,NFIX,RFIX,NSUPER,1,NFIX,REXPANDFIX,TYPEFIX)

! Allocate and create moving positions and supercell
NEXPANDVAR = (2*NSUPER+1)**3*NVAR
ALLOCATE( REXPANDVAR(3,NEXPANDVAR) , TYPEVAR(NEXPANDVAR) , DISTVAR(NVAR) , RMTVAR(NVAR) , RHVAR(NVAR))
CALL EXPANDUC3D(BRAVAIS,NVAR,RVAR,NSUPER,1,NVAR,REXPANDVAR,TYPEVAR)


! Find initial volume filling (before moving atoms are introduced)
N0 = 0
NN0 = 0
RHFIX(:) = 0.D0
RHVAR(:) = 0.D0

CALL VOLUMEFILH(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,N0,RVAR,NN0,REXPANDVAR,TYPEVAR, &  ! >
                     RMTVAR,RMTFIX,RHVAR,RHFIX,VOLMT,VOLUME,VOLFIL)                                ! <
WRITE(*,*) 'Initial volume filling (without empty spheres):',VOLFIL 


! Define hard-sphere radius of fixed atoms
IF (METHOD.EQ.-1) RHFIX(1:NFIX) = RMTFIX(1:NFIX) ! equal to their MT sphere
! Define hard-sphere radius of moving atoms
RHVAR(1:NVAR) = 0.D0 ! equal to 0


! Initialize total energy
CALL ETOTAL(NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR,RHFIX,RHVAR,HSVNUM,EINTRN,EEXTRN,ETOT)
WRITE(*,FMT='(A22,E16.8)') 'Initial total energy=',ETOT
WRITE(*,*) 'Initial number of hard-sphere violations:',HSVNUM


IF (METHOD.EQ.0) THEN
   WRITE(*,*) 'Maximize volume filling, METHOD=',METHOD
ELSE IF (METHOD.EQ.1) THEN
   WRITE(*,*) 'Minimize energy, METHOD=',METHOD
ELSE IF (METHOD.EQ.-1) THEN
   WRITE(*,*) 'Minimize energy with hard spheres for fixed atoms, METHOD=',METHOD
ELSE
   WRITE(*,*) 'METHOD should be -1,0, or 1 but on input METHOD=',METHOD
   STOP 'Invalid METHOD'
ENDIF


WRITE(*,*) 'Number of MC steps=',MCSTEPS
WRITE(*,*) 'Starting MC, number of temperatures:',NTEMPRMC

OPEN(68,FILE='pos_var_intemediate.dat')
COUNT = 0
TEMPR = TSTART
SHIFTMAX = SHIFTSTART
! Main temperature loop
DO ITEMPR = 1,NTEMPRMC



IF (METHOD.EQ.0) THEN
   CALL METROPOLIS2( &
        ISEED,MCSTEPS,TEMPR,SHIFTMAX,BRAVAIS,BRINV,NSUPER,NFIX,NEXPANDFIX,RFIX,REXPANDFIX,TYPEFIX,NVAR,NEXPANDVAR, & !>
        RVAR,REXPANDVAR,TYPEVAR,ETOT, & ! X
        EFLUCT,ACCEPTANCE)  ! <
ELSE
   CALL METROPOLIS( &
        ISEED,MCSTEPS,TEMPR,SHIFTMAX,BRAVAIS,BRINV,NSUPER,NFIX,NEXPANDFIX,RFIX,REXPANDFIX,TYPEFIX,NVAR,NEXPANDVAR,RHFIX,RHVAR, & !>
        RVAR,REXPANDVAR,TYPEVAR,EINTRN,EEXTRN,ETOT, & ! X
        EFLUCT,ACCEPTANCE)  ! <
ENDIF

   CALL VOLUMEFILH(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR, &  ! >
                     RMTVAR,RMTFIX,RHVAR,RHFIX,VOLMT,VOLUME,VOLFIL)                                                    ! <

   CALL ETOTAL(NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR,RHFIX,RHVAR,HSVNUM,EINTRN1,EEXTRN1,ETOT1)



   WRITE(*,FMT='(i6,a4,e10.3,a11,e10.3,a7,e9.2,a10,e12.4,a10,e10.3,a13,e16.8)') &
        ITEMPR,'T=',TEMPR,'Shiftmax=',SHIFTMAX,'Acc.=',ACCEPTANCE/DFLOAT(NVAR*MCSTEPS),'ENERGY=',ETOT,'Efluct=',EFLUCT,'Vol.fill.=',VOLFIL
   IF (HSVNUM.NE.0.D0) WRITE(*,*) 'Number of hard-sphere violations:',HSVNUM


!   TEMPR = DABS(EFLUCT)/NVAR/TDIV
   TEMPR = TEMPR / TDIV
   SHIFTMAX = SHIFTMAX / SDIV
   IF (ACCEPTANCE.EQ.0.D0) THEN
      COUNT = COUNT + 1
      SHIFTMAX = SHIFTMAX / 10.D0
   ENDIF
   IF (COUNT.GT.0.AND.COUNT.LE.10.AND.ACCEPTANCE.NE.0.D0) COUNT = 0
   IF (ACCEPTANCE.EQ.0.D0.AND.COUNT.GT.10) THEN ! save and anneal
      WRITE(68,*) 'New configuration'
      WRITE(68,FMT='(2E16.8,A16)') VOLFIL,ETOT, 'volfil, etot'
      WRITE(68,FMT='(3E16.8)') RVAR
      TEMPR = TSTART
      SHIFTMAX = SHIFTSTART
      COUNT = 0
   ENDIF


ENDDO ! ITEMPR = 1,NTEMPRMC
! End of main temperature loop

CLOSE(68)

! Evaluate final total energy
WRITE(*,fmt='(a66,3e16.8)') 'Final total energy from differences (internal, external, sum):   ',EINTRN,EEXTRN,ETOT
CALL ETOTAL(NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR,RHFIX,RHVAR,HSVNUM,EINTRN,EEXTRN,ETOT)
WRITE(*,fmt='(a66,3e16.8)') 'Final total energy, direct calculation (internal, external, sum):',EINTRN,EEXTRN,ETOT
WRITE(*,*) 'Final number of hard-sphere violations:',HSVNUM

! Evaluate volume filling fraction
CALL VOLUMEFILH(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR, &  ! >
                     RMTVAR,RMTFIX,RHVAR,RHFIX,VOLMT,VOLUME,VOLFIL)                                                    ! <
WRITE(*,*) 'Unit cell volume: ', VOLUME
WRITE(*,*) 'Muffin-tin volume:', VOLMT
WRITE(*,*) 'Volume filling fraction',VOLFIL

! Write out results in files
OPEN(67,FILE='positions_var.dat')
WRITE(67,*) '# ',NVAR
WRITE(67,*) '# x y z rmt, rmt**2 cartesian coordinates'
DO IVAR = 1,NVAR
   WRITE(67,FMT='(5E16.8)') (RVAR(IX,IVAR),IX=1,3),RMTVAR(IVAR),RMTVAR(IVAR)**2
ENDDO
CLOSE(67)
OPEN(67,FILE='positions_fix.dat')
WRITE(67,*) '# ',NFIX
WRITE(67,*) '# x y z rmt, rmt**2 cartesian coordinates'
DO IFIX = 1,NFIX
   WRITE(67,FMT='(5E16.8)') (RFIX(IX,IFIX),IX=1,3),RMTFIX(IFIX),RMTFIX(IFIX)**2
ENDDO
CLOSE(67)
WRITE(*,fmt='(a84)') 'End positions of empty-spheres in cartesian coords written in file positions_var.dat'
WRITE(*,*) 'Copy this to positions_var_old.dat if you wish to restart from these.'
! Write out resuls to be copy-pasted for Voronoi input
OPEN(67,FILE='input_voronoi_cartesian.txt')
WRITE(67,*) '#Copy-paste this into the inputcard for the Voronoi program.'
WRITE(67,*) '#First come the fixed coords, then the empty-sphere coords.'
WRITE(67,FMT='(A7)') 'BRAVAIS'
WRITE(67,FMT='(2X,3E16.8)') BRAVAIS(1:3,1:3)
WRITE(67,FMT='(A6,I7)')'NAEZ= ', NFIX + NVAR
WRITE(67,FMT='(A18)')'CARTESIAN= .TRUE. '
WRITE(67,FMT='(A76)')'RBASIS                                             <MTWAL> (rmt, alat units)'
DO IFIX = 1,NFIX
   WRITE(67,FMT='(3E16.8,3X,E16.8)') (RFIX(IX,IFIX),IX=1,3),RMTFIX(IFIX)
ENDDO
DO IVAR = 1,NVAR
   WRITE(67,FMT='(3E16.8,3X,E16.8)') (RVAR(IX,IVAR),IX=1,3),RMTVAR(IVAR)
ENDDO
CLOSE(67)


! Write out results in internal coordinates in files
OPEN(67,FILE='positions_var_intcoord.dat')
WRITE(67,*) '# ',NVAR
WRITE(67,*) '# x y z internal coordinates; rmt , rmt**2'
DO IVAR = 1,NVAR
   RTMP(1:3) = RVAR(1:3,IVAR)
   CALL CART2INTERN(BRINV,RTMP)
   WRITE(67,FMT='(5E16.8)') (RTMP(IX),IX=1,3),RMTVAR(IVAR),RMTVAR(IVAR)**2
ENDDO
CLOSE(67)
OPEN(67,FILE='positions_fix_intcoord.dat')
WRITE(67,*) '# ',NFIX
WRITE(67,*) '# x y z internal coordinates; rmt , rmt**2'

DO IFIX = 1,NFIX
   RTMP(1:3) = RFIX(1:3,IFIX)
   CALL CART2INTERN(BRINV,RTMP)
   WRITE(67,FMT='(5E16.8)') (RTMP(IX),IX=1,3),RMTFIX(IFIX),RMTFIX(IFIX)**2
ENDDO
CLOSE(67)

! Write out resuls to be copy-pasted for Voronoi input
OPEN(67,FILE='input_voronoi_internal.txt')
WRITE(67,*) '#Copy-paste this into the inputcard for the Voronoi program.'
WRITE(67,*) '#First come the fixed coords, then the empty-sphere coords.'
WRITE(67,FMT='(A7)') 'BRAVAIS'
WRITE(67,FMT='(2X,3E16.8)') BRAVAIS(1:3,1:3)
WRITE(67,FMT='(A6,I7)')'NAEZ= ',NFIX + NVAR
WRITE(67,FMT='(A18)')'CARTESIAN= .FALSE.'
WRITE(67,FMT='(A76)')'RBASIS                                             <MTWAL> (rmt, alat units)'
DO IFIX = 1,NFIX
   RTMP(1:3) = RFIX(1:3,IFIX)
   CALL CART2INTERN(BRINV,RTMP)
   WRITE(67,FMT='(3E16.8,3X,E16.8)') (RTMP(IX),IX=1,3),RMTFIX(IFIX)
ENDDO
DO IVAR = 1,NVAR
   RTMP(1:3) = RVAR(1:3,IVAR)
   CALL CART2INTERN(BRINV,RTMP)
   WRITE(67,FMT='(3E16.8,3X,E16.8)') (RTMP(IX),IX=1,3),RMTVAR(IVAR)
ENDDO
CLOSE(67)



OPEN(67,FILE='MTradii.dat')
WRITE(67,FMT='(A30,3E16.8)') '# Volume: MT, tot, filling:',VOLMT,VOLUME,VOLFIL
WRITE(67,FMT='(A20,I6)') '# Fixed atoms:',NFIX
WRITE(67,FMT='(A30)') '# RMT and volume per atom:'
DO IFIX = 1,NFIX
   WRITE(67,FMT='(I6,2E16.8)') IFIX,RMTFIX(IFIX), 4.D0 * PI * RMTFIX(IFIX)**3 /3.D0
ENDDO
WRITE(67,FMT='(A20,I6)') '# Moving atoms:',NVAR
DO IVAR = 1,NVAR
   WRITE(67,FMT='(I6,2E16.8)') IVAR,RMTVAR(IVAR), 4.D0 * PI * RMTVAR(IVAR)**3 /3.D0
ENDDO
CLOSE(67)


! Finito la musica
STOP 'TELOS KALO OLA KALA'
END PROGRAM EMPTYCELLS

!========================================================================
SUBROUTINE DISTUC3D(RPOS,NEXPAND,REXPAND,ATYPE,NBASIS,DISTANCE)
implicit none
! Given a position RPOS(3) this subroutine returns the distance to the positions of all atoms
! in the unit cell under supercell conditions (i.e. for every particular atom-type JTYPE it compares
! the distance to all atoms of JTYPE in the supercell and chooses the smallest).

! Input:
REAL*8 RPOS(3)      
REAL*8 REXPAND(3,*) ! REXPAND(3,NEXPAND) : Positions of atoms in supercell
INTEGER NEXPAND     ! Nubmer of atoms in supercell
INTEGER NBASIS      ! Number of basis atoms (atom types)
INTEGER ATYPE(*)

! Output:
REAL*8 DISTANCE(*)  ! DISTANCE(NBASIS) distance to all atoms in primitive cell
!INTEGER NEIGHB(NBASIS) ! Neighbour atom with this distance
!REAL*8 RNEIGHB(3,NBASIS)

! Internal
REAL*8 RAT(3) ! Position of probe atom
REAL*8 DPR    ! Distance to probe atom
INTEGER ITYPE,IAT,ICELL,IX

DISTANCE(1:NBASIS) = 1.D100
DO IAT = 1,NEXPAND
   ITYPE = ATYPE(IAT)
   RAT(1:3) = REXPAND(1:3,IAT)
   DPR = (RAT(1) - RPOS(1))**2 + (RAT(2) - RPOS(2))**2 + (RAT(3) - RPOS(3))**2
   IF (DPR.LT.DISTANCE(ITYPE)) THEN 
      DISTANCE(ITYPE) = DPR
!      INEIGHB(ITYPE) = IAT
   ENDIF
ENDDO

DO ITYPE = 1,NBASIS
   DISTANCE(ITYPE) = DSQRT(DISTANCE(ITYPE))
!   IAT = NEIGHB(ITYPE)
!   DO IX = 1,3
!      RNEIGHB(IX,ITYPE) = REXPAND(IX,IAT) - RPOS(IX) ! Vector pointing to the neighbor
!   ENDDO
ENDDO

RETURN
END SUBROUTINE DISTUC3D

!========================================================================
SUBROUTINE EXPANDUC3D(BRAVAIS,NBASIS,RBASIS,NSUPER,ISTART,IEND,REXPAND,ATYPE)
implicit none
! Given a unit cell in 3D with Bravais vectors and basis atoms, this subroutine
! expands the cell to form a NxNxN super-cell by adding Bravais vectors to the basis atoms.

! Input:
REAL*8 BRAVAIS(3,3) ! Bravais (1,2,3 ; x,y,z)
INTEGER NBASIS      ! Number of basis atoms
REAL*8 RBASIS(3,*)  ! Position of basis atoms (3,NBASIS)
INTEGER NSUPER  ! How large is the supercell (2*NSUPER+1)**3, usually NSUPER=1 is enough
INTEGER ISTART,IEND ! Start-atom and end-atom, only expand for these.

! Output
REAL*8 REXPAND(3,*)    ! Position of original + shifted atoms (3,(2*NSUPER+1)**3 * NBASIS)
                       ! First come all atoms of cell (-1,-1,-1), then of cell (0,-1,-1), then (+1,-1,-1)
                       ! in a linear array including all triplets (i,j,k) with -NSUPER<i,j,k<+NSUPER
INTEGER ATYPE(*)       ! Type of atom in expanded cell corresponding to the
                       ! original atom type in the basis.

! Internal:
INTEGER IX,IBAS,IB1,IB2,IB3,INEW
INTEGER IB1C,IB2C,IB3C,ICELL,NLIN
REAL*8 RSHIFT(3)



NLIN = 2*NSUPER + 1
DO IB3 = -NSUPER,NSUPER
   IB3C = IB3 + NSUPER  ! IBC1,2,3 are for counter purposes, counting 0...2*NSUPER
   DO IB2 = -NSUPER,NSUPER
      IB2C = IB2 + NSUPER 
      DO IB1 = -NSUPER,NSUPER
         IB1C = IB1 + NSUPER 

         ICELL = NBASIS * (NLIN**2 * IB3C + NLIN * IB2C + IB1C) ! ICELL counts 0,NBASIS,2*NBASIS,3*NBASIS,...
         RSHIFT(1) = IB1*BRAVAIS(1,1) + IB2*BRAVAIS(1,2) + IB3*BRAVAIS(1,3)
         RSHIFT(2) = IB1*BRAVAIS(2,1) + IB2*BRAVAIS(2,2) + IB3*BRAVAIS(2,3)
         RSHIFT(3) = IB1*BRAVAIS(3,1) + IB2*BRAVAIS(3,2) + IB3*BRAVAIS(3,3)


         DO IBAS = ISTART,IEND
            INEW = IBAS + ICELL
            ATYPE(INEW) = IBAS
            REXPAND(1,INEW) = RBASIS(1,IBAS) + RSHIFT(1)
            REXPAND(2,INEW) = RBASIS(2,IBAS) + RSHIFT(2)
            REXPAND(3,INEW) = RBASIS(3,IBAS) + RSHIFT(3)
         ENDDO

      ENDDO
   ENDDO
ENDDO


RETURN
END SUBROUTINE EXPANDUC3D

!========================================================================

SUBROUTINE RANDPOS(BRAVAIS,BRINV,ISEED,SHIFTMAX,ROLD,RNEW)
implicit none
! Pick a position in the unit cell by a random shift of ROLD
! Input
REAL*8 BRAVAIS(3,3),BRINV(3,3)
REAL*8 ROLD(3) 
REAL*8 SHIFTMAX
INTEGER ISEED

! Output
REAL*8 RNEW(3)

! Internal
REAL*8 RSHIFTI(3),ROLDI(3),RNEWI(3)  ! Ending "I" means internal coords.
INTEGER IX,IX1,IB
REAL*8 XRAND
EXTERNAL XRAND

! Shift in internal (Bravais) coordinates
RSHIFTI(1) = SHIFTMAX * (XRAND(ISEED) - 0.5D0)
RSHIFTI(2) = SHIFTMAX * (XRAND(ISEED) - 0.5D0)
RSHIFTI(3) = SHIFTMAX * (XRAND(ISEED) - 0.5D0)

! ROLD in internal coordinates:
ROLDI(:) = ROLD(:)
CALL CART2INTERN(BRINV,ROLDI)

!Shift to RNEW in internal coords:
DO IB = 1,3
   RNEWI(IB) = ROLDI(IB) + RSHIFTI(IB)
   RNEWI(IB) = RNEWI(IB) - DFLOAT(FLOOR(RNEWI(IB))) ! shift back to primitive cell (make integer_part=0)
   IF (RNEWI(IB).LT.0.D0.OR.RNEWI(IB).GE.1.D0) STOP 'error in RANDPOS'
!   IF (RNEWI(IB).GE.1.D0) RNEWI(IB) = RNEWI(IB) - INT(RNEWI(IB))
!   IF (RNEWI(IB).LT.0.D0) RNEWI(IB) = RNEWI(IB) - INT(RNEWI(IB)) + 1.D0
ENDDO


! Back to Cartesian coords.
RNEW(:) = RNEWI(:)

CALL CART2INTERN(BRAVAIS,RNEW)


RETURN
END SUBROUTINE RANDPOS

!========================================================================
      subroutine randatom(natom,iseed,iat)
! Choose randomly one atom 
      implicit none
      integer natom ! Input, Number of atoms
      integer iseed ! Input, Random number seed
      integer iat   ! Output, chosen atom
      real*8 xrand  ! Random number function
      external xrand

      iat = int(natom * xrand(iseed)) + 1
! Integer number generation seems to work at least up to 10**6
! (xrand seems reliable to at least 6th digit) if mersene twister is used.

      end subroutine randatom
!========================================================================

SUBROUTINE EEXT(DISTFIX,NFIX,RHFIX,RHVAR1,EEXTI,LHSV)
implicit none
! Interaction energy of a moving atom to all fixed atoms
! Input
REAL*8 DISTFIX(*),RHFIX(*),RHVAR1 ! Distance and hard-sphere radii
INTEGER NFIX
! Output
REAL*8 EEXTI
LOGICAL LHSV ! Hard-sphere violation
! Internal
INTEGER IFIX
REAL*8 DD,EIJ,DHARD


EEXTI = 0.D0
LHSV = .FALSE.
DO IFIX = 1,NFIX
   DHARD = RHFIX(IFIX) + RHVAR1 ! sum of hard-sphere radii
   DD = DISTFIX(IFIX)
   IF (DD.LT.DHARD) LHSV = .TRUE.
   CALL EPAIR(DD,DHARD,EIJ)
   EEXTI = EEXTI + EIJ
ENDDO

RETURN
END SUBROUTINE EEXT

!========================================================================

SUBROUTINE EINT(DISTVAR,NVAR,IVAR,RHVAR,EINTI,LHSV)
implicit none
! Interaction energy of a moving atom to all moving atoms except itself
! Input
REAL*8 DISTVAR(*),RHVAR(*)
INTEGER NVAR,IVAR
! Output
REAL*8 EINTI ! energy 
LOGICAL LHSV ! Hard-sphere violation
! Internal
INTEGER IVAR1,IVAR2
REAL*8 DD,EIJ,DHARD,RHVAR1


RHVAR1 = RHVAR(IVAR)
EINTI = 0.D0
LHSV = .FALSE.
DO IVAR2 = 1,NVAR
   IF (IVAR2.NE.IVAR) THEN
      DHARD = RHVAR(IVAR2) + RHVAR1 ! sum of hard-sphere radii
      DD = DISTVAR(IVAR2)
      IF (DD.LT.DHARD) LHSV = .TRUE.
      CALL EPAIR(DD,DHARD,EIJ)
      EINTI = EINTI + EIJ
   ENDIF
ENDDO

RETURN
END SUBROUTINE EINT

!========================================================================

SUBROUTINE EPAIR(DIST,DHARD,EIJ)
implicit none
! Pair interaction energy
! Input:
REAL*8 DIST,DHARD ! Distance, sum of hard-sphere radii
! Output
REAL*8 EIJ,FIJ  ! Pair energy,force

IF (DIST.LT.DHARD) THEN
   EIJ = 1.D4
ELSE
   EIJ = 1./(DIST-DHARD+1.D-4)
ENDIF

END SUBROUTINE EPAIR
!========================================================================

SUBROUTINE FPAIR(DIST,FIJ)
implicit none
! Pair interaction energy
! Input:
REAL*8 DIST ! Distance
! Output
REAL*8 EIJ,FIJ  ! Pair energy,force

IF (DIST.LT.1.D-3) THEN
   FIJ = -1.D20
ELSE
   EIJ = -DIST
ENDIF

END SUBROUTINE FPAIR

!========================================================================

SUBROUTINE ETOTAL(NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR,RHFIX,RHVAR,HSVNUM,EINTRN,EEXTRN,ETOT)
implicit none
! Total energy
! Input:
REAL*8 RFIX(3,*),RVAR(3,*),REXPANDFIX(3,*),REXPANDVAR(3,*),RHFIX(*),RHVAR(*)
REAL*8 HSVNUM ! Total of hard-sphere violations
INTEGER NFIX,NVAR,NEXPANDFIX,NEXPANDVAR
INTEGER TYPEFIX,TYPEVAR(*)
! Output:
REAL*8 ETOT,EINTRN,EEXTRN ! Total, internal and external energy of the system of moving atoms
! Internal 
INTEGER IVAR,IFIX,IVAR1,IVAR2
REAL*8 DD,EIJ,DHARD,RHVAR1
REAL*8 RVEC(3)
REAL*8,ALLOCATABLE :: DISTFIX(:),DISTVAR(:)

ALLOCATE(DISTFIX(NFIX),DISTVAR(NVAR))

ETOT = 0.D0
EINTRN = 0.D0
EEXTRN = 0.D0
HSVNUM = 0.D0

DO IVAR = 1,NVAR
   RVEC(:) = RVAR(1:3,IVAR)
   ! First the interaction with fixed atoms
   CALL DISTUC3D(RVEC,NEXPANDFIX,REXPANDFIX,TYPEFIX,NFIX,DISTFIX)
   RHVAR1 = RHVAR(IVAR)
   DO IFIX = 1,NFIX
      DHARD = RHFIX(IFIX) + RHVAR1 ! Sum of hard-sphere radii
      DD = DISTFIX(IFIX)
      IF (DD.LT.DHARD) HSVNUM = HSVNUM + 1.D0
      CALL EPAIR(DD,DHARD,EIJ)
      EEXTRN = EEXTRN + EIJ
   ENDDO
   ! Now with moving atoms
   CALL DISTUC3D(RVEC,NEXPANDVAR,REXPANDVAR,TYPEVAR,NVAR,DISTVAR)
   DO IVAR2 = IVAR + 1,NVAR
      DHARD = RHVAR(IVAR) + RHVAR1 ! Sum of hard-sphere radii
      DD = DISTVAR(IVAR2)
      IF (DD.LT.DHARD) HSVNUM = HSVNUM + 1.D0
      CALL EPAIR(DD,DHARD,EIJ)
      EINTRN = EINTRN + EIJ
   ENDDO
ENDDO

ETOT = EINTRN + EEXTRN


RETURN
END SUBROUTINE ETOTAL
!========================================================================

SUBROUTINE VOLUMEFIL(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR, &  ! >
                     RMTVAR,RMTFIX,VOLMT,VOLUME,VOLFIL)                                                          ! <
implicit none
! Volume filling
! Input:
REAL*8 BRAVAIS(3,3)
REAL*8 RFIX(3,*),RVAR(3,*),REXPANDFIX(3,*),REXPANDVAR(3,*)
INTEGER NFIX,NVAR,NEXPANDFIX,NEXPANDVAR
INTEGER TYPEFIX(*),TYPEVAR(*)
! Output:
REAL*8 VOLMT,VOLUME,VOLFIL ! Volume of MT-spheres, of unit cell, and volume filling fraction
REAL*8 RMTVAR(*),RMTFIX(*) ! Muffin tin radii of moving and fixed atoms
! Internal 
INTEGER IVAR,IFIX,IVAR1,IAT,IAT1,NTOT,NEXPANDTOT
REAL*8 DD
REAL*8 RVEC(3)
REAL*8,ALLOCATABLE :: RMT(:),RAT(:,:),REXPAND(:,:)
REAL*8 PI

PI = 4.D0 * DATAN(1.D0)

NTOT = NVAR + NFIX
NEXPANDTOT = NEXPANDVAR + NEXPANDFIX
ALLOCATE(RMT(NTOT),RAT(3,NTOT),REXPAND(3,NEXPANDTOT))

!Pass all atoms in one array
DO IAT = 1,NVAR
   RAT(1:3,IAT) = RVAR(1:3,IAT)
ENDDO
DO IFIX = 1,NFIX
   IAT = NVAR + IFIX
   RAT(1:3,IAT) = RFIX(1:3,IFIX)
ENDDO

DO IAT = 1,NEXPANDVAR
   REXPAND(1:3,IAT) = REXPANDVAR(1:3,IAT)
ENDDO
DO IFIX = 1,NEXPANDFIX
   IAT = NEXPANDVAR + IFIX
   REXPAND(1:3,IAT) = REXPANDFIX(1:3,IFIX)
ENDDO

RMT(1:NTOT) = 1.D100

DO IAT = 1,NTOT
   DO IAT1 = 1,NEXPANDTOT
      DD = (RAT(1,IAT)-REXPAND(1,IAT1))**2+(RAT(2,IAT)-REXPAND(2,IAT1))**2+(RAT(3,IAT)-REXPAND(3,IAT1))**2
      IF (DD.LT.RMT(IAT).AND.DD.GT.0.D0) RMT(IAT) = DD
   ENDDO
   RMT(IAT) = 0.5D0 * DSQRT(RMT(IAT))
ENDDO



VOLMT = 0.D0
DO IAT = 1,NTOT
   VOLMT = VOLMT + RMT(IAT)**3
ENDDO
VOLMT = 4.D0 * PI * VOLMT / 3.D0

VOLUME =  BRAVAIS(1,1) * (BRAVAIS(2,2)*BRAVAIS(3,3)-BRAVAIS(2,3)*BRAVAIS(3,2)) &
        + BRAVAIS(2,1) * (BRAVAIS(3,2)*BRAVAIS(1,3)-BRAVAIS(1,2)*BRAVAIS(3,3)) &
        + BRAVAIS(3,1) * (BRAVAIS(1,2)*BRAVAIS(2,3)-BRAVAIS(2,2)*BRAVAIS(1,3)) 



VOLUME = DABS(VOLUME)
VOLFIL = VOLMT / VOLUME

DO IAT = 1,NVAR
   RMTVAR(IAT) = RMT(IAT)
ENDDO
DO IAT = NVAR + 1,NTOT
   IFIX = IAT - NVAR
   RMTFIX(IFIX) = RMT(IAT)
ENDDO

!WRITE(*,*) 'Unit Cell Volume :',VOLUME
!WRITE(*,*) 'Muffin Tin Volume:',VOLMT
!WRITE(*,*) 'Volume fill fraction:',VOLFIL


RETURN
END SUBROUTINE VOLUMEFIL
!========================================================================

SUBROUTINE VOLUMEFILH_OLD(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR, &  ! >
                     RMTVAR,RMTFIX,RHVAR,RHFIX,VOLMT,VOLUME,VOLFIL)                                                          ! <
implicit none
! Volume filling
! Input:
REAL*8 BRAVAIS(3,3)
REAL*8 RFIX(3,*),RVAR(3,*),REXPANDFIX(3,*),REXPANDVAR(3,*),RHVAR(*),RHFIX(*)
INTEGER NFIX,NVAR,NEXPANDFIX,NEXPANDVAR
INTEGER TYPEFIX(*),TYPEVAR(*)
! Output:
REAL*8 VOLMT,VOLUME,VOLFIL ! Volume of MT-spheres, of unit cell, and volume filling fraction
REAL*8 RMTVAR(*),RMTFIX(*) ! Muffin tin radii of moving and fixed atoms
! Internal 
INTEGER IVAR,IFIX,IVAR1,IAT,IAT1,NTOT,NEXPANDTOT
REAL*8 DD
REAL*8 RVEC(3)
REAL*8,ALLOCATABLE :: RMT(:),RAT(:,:),REXPAND(:,:),RH(:)
INTEGER,ALLOCATABLE :: ATYPE(:)
REAL*8 PI

PI = 4.D0 * DATAN(1.D0)

NTOT = NVAR + NFIX
NEXPANDTOT = NEXPANDVAR + NEXPANDFIX
ALLOCATE(RMT(NTOT),RAT(3,NTOT),REXPAND(3,NEXPANDTOT),ATYPE(NEXPANDTOT),RH(NTOT))

! Pass all atoms in one array, first moving atoms, then fixed atoms
DO IAT = 1,NVAR
   RAT(1:3,IAT) = RVAR(1:3,IAT)   ! Atom positions
   RH(IAT) = RHVAR(IAT)           ! Hard sphere radii
ENDDO
DO IFIX = 1,NFIX
   IAT = NVAR + IFIX
   RAT(1:3,IAT) = RFIX(1:3,IFIX)  ! Atom positions
   RH(IAT) = RHFIX(IFIX)          ! Hard sphere radii
ENDDO

DO IAT = 1,NEXPANDVAR
   REXPAND(1:3,IAT) = REXPANDVAR(1:3,IAT)
   ATYPE(IAT) = TYPEVAR(IAT)
ENDDO
DO IFIX = 1,NEXPANDFIX
   IAT = NEXPANDVAR + IFIX
   REXPAND(1:3,IAT) = REXPANDFIX(1:3,IFIX)
   ATYPE(IAT) = NVAR + TYPEFIX(IFIX)
ENDDO


RMT(1:NTOT) = 1.D100

DO IAT = 1,NTOT
   DO IAT1 = 1,NEXPANDTOT
      DD = (RAT(1,IAT)-REXPAND(1,IAT1))**2+(RAT(2,IAT)-REXPAND(2,IAT1))**2+(RAT(3,IAT)-REXPAND(3,IAT1))**2
      DD = DSQRT(DD)
      DD = MIN(0.5D0*DD,DD-RH(ATYPE(IAT1)))
      IF (DD.LT.RMT(IAT).AND.DD.GT.1.D-10) RMT(IAT) = DD
   ENDDO
   ! In case RMT became smaller than hard-sphere radius, return it to hard-sphere radius
   RMT(IAT) = MAX(RMT(IAT),RH(IAT))
ENDDO




VOLMT = 0.D0
DO IAT = 1,NTOT
   VOLMT = VOLMT + RMT(IAT)**3
ENDDO
VOLMT = 4.D0 * PI * VOLMT / 3.D0

VOLUME =  BRAVAIS(1,1) * (BRAVAIS(2,2)*BRAVAIS(3,3)-BRAVAIS(2,3)*BRAVAIS(3,2)) &
        + BRAVAIS(2,1) * (BRAVAIS(3,2)*BRAVAIS(1,3)-BRAVAIS(1,2)*BRAVAIS(3,3)) &
        + BRAVAIS(3,1) * (BRAVAIS(1,2)*BRAVAIS(2,3)-BRAVAIS(2,2)*BRAVAIS(1,3)) 



VOLUME = DABS(VOLUME)
VOLFIL = VOLMT / VOLUME

DO IAT = 1,NVAR
   RMTVAR(IAT) = RMT(IAT)
ENDDO
DO IAT = NVAR + 1,NTOT
   IFIX = IAT - NVAR
   RMTFIX(IFIX) = RMT(IAT)
ENDDO


!WRITE(*,*) 'Unit Cell Volume :',VOLUME
!WRITE(*,*) 'Muffin Tin Volume:',VOLMT
!WRITE(*,*) 'Volume fill fraction:',VOLFIL

DEALLOCATE(RMT,RAT,REXPAND,ATYPE,RH)


RETURN
END SUBROUTINE VOLUMEFILH_OLD

!========================================================================
    FUNCTION XRAND(INIT)
      USE STDTYPES
      USE MTPRNG
      IMPLICIT NONE
      REAL*8 XRAND,URAND
      INTEGER INIT
      LOGICAL LFIRST
      TYPE(MTPRNG_STATE) :: STATE
      SAVE STATE
      DATA LFIRST/.TRUE./
      SAVE LFIRST


! Initialize in first call:
      IF (LFIRST) THEN
         CALL MTPRNG_INIT(INIT,STATE)
         LFIRST = .FALSE.
      ENDIF

      XRAND =  MTPRNG_RAND_REAL2(STATE)

    END FUNCTION XRAND


!========================================================================
SUBROUTINE READPARA(BRAVAIS,CARTESIAN,NFIX,NVAR,ISEED,NSUPER,METHOD,TEMPR,NTEMPRMC,TDIV,SHIFTMAX,SDIV,MCSTEPS)
implicit none

INTEGER NFIX ! Number of fixed atoms
INTEGER NVAR ! Number of moving atoms with variable positions to be optimized
INTEGER NSUPER ! Supercell extends from -NSUPER to +NSUPER in 3 dimensions
INTEGER NEXPANDVAR,NEXPANDFIX ! Number of moving and fixed atoms in supercell
REAL*8 BRAVAIS(3,3)
REAL*8 SHIFTMAX,TEMPR  ! Max. shift for moving the atom with respect to Bravais vector, temperature
REAL*8 SDIV,TDIV
INTEGER MCSTEPS,ISEED,ITEMPR,ISTEP,IVAR1,NTEMPRMC,METHOD
LOGICAL CARTESIAN ! True if positions are in cartesian units; false if they are in Bravais units.
INTEGER N0,NN0
LOGICAL LACCEPT,LREAD
INTEGER IFIX,IVAR,IX,IB
CHARACTER*1 HASH
CHARACTER*200 UIO  ! For IoInput routine
INTEGER IER        ! For IoInput routine


! Read in Bravais vectors
CALL IoInput('BRAVAIS   ',UIO,1,7,IER)
IF (IER.NE.0) STOP 'Key BRAVAIS not found in inputcard'
DO IB = 1,3
   CALL IoInput('BRAVAIS   ',UIO,IB,7,IER)
   READ (UNIT=UIO,FMT=*) (BRAVAIS(IX,IB), IX=1,3)
ENDDO
WRITE(*,FMT='(A24)') 'BRAVAIS vectors read in.'
WRITE(*,FMT='(3E16.8)') ((BRAVAIS(IX,IB),IX=1,3),IB=1,3)
WRITE(21,FMT='(A10)') 'BRAVAIS   '
WRITE(21,FMT='(3E16.8)') ((BRAVAIS(IX,IB),IX=1,3),IB=1,3)


! Read number of fixed positions
CALL IoInput('NAEZ      ',UIO,1,7,IER)
IF (IER.NE.0) STOP 'Key NAEZ not found in inputcard'
READ (UNIT=UIO,FMT=*) NFIX
WRITE(*,*) '--> Number of fixed basis atoms NAEZ=',NFIX
WRITE(21,*)    'Number of fixed basis atoms NAEZ=',NFIX

! Read if positions are in cartesian coords
CARTESIAN = .TRUE.
CALL IoInput('CARTESIAN ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default cartesian coords'
ELSE
   READ (UNIT=UIO,FMT=*) CARTESIAN
ENDIF
WRITE(*,*) '--> Cartesian coords: CARTESIAN=',CARTESIAN
WRITE(21,*)    'Cartesian coords: CARTESIAN=',CARTESIAN


! Read number of empty spheres
INQUIRE(FILE='positions_var_old.dat',EXIST=LREAD)
IF (LREAD) THEN
   WRITE(*,*) 'Reading number of empty-spheres from file positions_var_old.dat'
   OPEN(10,FILE='positions_var_old.dat')
   READ(10,*) HASH,NVAR
   CLOSE(10)
   WRITE(*,*) '--> No. of empty spheres NVAR=',NVAR
ELSE
   WRITE(*,*) 'Reading number of empty-sphere positions from inputcard'
   CALL IoInput('NEMPTYSPH ',UIO,1,7,IER)
   IF (IER.NE.0) STOP 'Key NEMPTYSPH for number of empty positions not found in inputcard'
   READ (UNIT=UIO,FMT=*) NVAR
   WRITE(*,*) '--> No. of empty spheres NEMPTYSPH=',NVAR
   WRITE(21,*)    'No. of empty spheres NEMPTYSPH=',NVAR
ENDIF

! Initialize random number seed
ISEED = 1
CALL IoInput('ISEED     ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default random number seed'
ELSE
   READ (UNIT=UIO,FMT=*) ISEED
ENDIF
WRITE(*,*) '--> Random number seed ISEED=',ISEED
WRITE(21,*)    'Random number seed ISEED=',ISEED

! Initialize supercell size
NSUPER = 1
CALL IoInput('NSUPER    ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default supercell size'
ELSE
   READ (UNIT=UIO,FMT=*) NSUPER
ENDIF
WRITE(*,*) '--> Supercell size NSUPER=',NSUPER
WRITE(21,*)    'Supercell size NSUPER=',NSUPER

! Initialize method (0=volume-fill ; 1=energy)
METHOD = -1
CALL IoInput('METHOD    ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default method'
ELSE
   READ (UNIT=UIO,FMT=*) METHOD
ENDIF
WRITE(*,*) '-->METHOD=',METHOD
WRITE(21,*)   'METHOD=',METHOD

! Initialize temperature
TEMPR = 1.D0
CALL IoInput('TMC       ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default starting temperature'
ELSE
   READ (UNIT=UIO,FMT=*) TEMPR
ENDIF
WRITE(*,*) '-->Starting temperature TMC=',TEMPR
WRITE(21,*)   'Starting temperature TMC=',TEMPR


! Read nuber of temperatures
NTEMPRMC = 20
CALL IoInput('NTEMPRMC  ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default number of temperatures'
ELSE
   READ (UNIT=UIO,FMT=*) NTEMPRMC
ENDIF
WRITE(*,*) '-->Number of temperature steps NTEMPRMC=',NTEMPRMC
WRITE(21,*)   'Number of temperature steps NTEMPRMC=',NTEMPRMC

! Read temperature divisor
TDIV = 2.D0
CALL IoInput('TDIV      ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default temperature divisor'
ELSE
   READ (UNIT=UIO,FMT=*) TDIV
ENDIF
WRITE(*,*) '-->Temperature divisor TDIV=',TDIV
WRITE(21,*)   'Temperature divisor TDIV=',TDIV

! Initialize max. shift
SHIFTMAX = 1.D0
CALL IoInput('SHIFTMAX  ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default max. shift'
ELSE
   READ (UNIT=UIO,FMT=*) SHIFTMAX
ENDIF
WRITE(*,*) '-->Starting max. shift SHIFTMAX=',SHIFTMAX
WRITE(21,*)   'Starting max. shift SHIFTMAX=',SHIFTMAX

! Read shift divisor
SDIV = 1.3D0
CALL IoInput('SDIV      ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default shift divisor'
ELSE
   READ (UNIT=UIO,FMT=*) SDIV
ENDIF
WRITE(*,*) '-->Shift divisor SDIV=',SDIV
WRITE(21,*)   'Shift divisor SDIV=',SDIV

! Read MC steps
MCSTEPS = 1000
CALL IoInput('MCSTEPS   ',UIO,1,7,IER)
IF (IER.NE.0) THEN 
   WRITE(*,*) 'Using default number of MC steps'
ELSE
   READ (UNIT=UIO,FMT=*) MCSTEPS
ENDIF
WRITE(*,*) '-->Number of MC steps MCSTEPS=',MCSTEPS
WRITE(21,*)   'Number of MC steps MCSTEPS=',MCSTEPS

RETURN
END SUBROUTINE READPARA

!========================================================================

SUBROUTINE READPOS(NFIX,NVAR,BRAVAIS,BRINV,CARTESIAN,ISEED,RFIX,RVAR)
implicit none
! Input:
INTEGER NFIX,NVAR,ISEED
REAL*8 BRAVAIS(3,3),BRINV(3,3)
LOGICAL CARTESIAN
! Output:
REAL*8 RFIX(3,NFIX),RVAR(3,NVAR)

!Internal
INTEGER IFIX,IVAR,IX,IB
REAL*8 SHIFT,ROLD(3),RNEW(3),RTMP(3)
CHARACTER*1 HASH
LOGICAL LREAD
CHARACTER*200 UIO  ! For IoInput routine
INTEGER IER        ! For IoInput routine

!WRITE(*,*) NFIX,NVAR,iseed!,ALLOCATED(RFIX)
!WRITE(*,*) 'size rfix',size(rfix,1),size(rfix,2)
!WRITE(*,*) 'size rvar',size(rvar,1),size(rvar,2)


! Read fixed positions
DO IFIX = 1,NFIX
   CALL IoInput('RBASIS    ',UIO,IFIX,7,IER)
   IF (IER.NE.0) STOP 'Key RBASIS not found in inputcard'
   READ (UNIT=UIO,FMT=*) (RFIX(IX,IFIX), IX=1,3)
ENDDO
WRITE(*,*) 'Fixed positions read in.'
WRITE(*,FMT='(A24)') 'RBASIS  vectors read in.'
WRITE(*,FMT='(3E16.8)') ((RFIX(IX,IFIX),IX=1,3),IFIX=1,NFIX)
WRITE(21,FMT='(A10)') 'RBASIS    '
WRITE(21,FMT='(3E16.8)') ((RFIX(IX,IFIX),IX=1,3),IFIX=1,NFIX)


! Note fivos: this was done after creating the empty sphere pos. before.
IF (.NOT.CARTESIAN) THEN ! Transform all positions to cartesian coord. system
   DO IFIX = 1,NFIX
      CALL CART2INTERN(BRAVAIS,RFIX(1,IFIX)) 
   ENDDO

!   DO IVAR = 1,NVAR
!      CALL CART2INTERN(BRAVAIS,RVAR(1,IVAR)) ! Note fivos: this was not commented out before
!   ENDDO
ENDIF



! Initialize moving positions
INQUIRE(FILE='positions_var_old.dat',EXIST=LREAD)
IF (LREAD) THEN
   WRITE(*,*) 'Reading initial empty-sphere positions from file positions_var_old.dat'
   OPEN(10,FILE='positions_var_old.dat')
   READ(10,*) HASH
   READ(10,*) HASH
   DO IVAR = 1,NVAR
      READ(10,*) (RVAR(IX,IVAR),IX=1,3)
   ENDDO
   CLOSE(10)
ELSE
   WRITE(*,*) 'Constructing initial empty-sphere positions randomly'
   ROLD(1:3) = 0.D0 ! Origin
   DO IVAR = 1,NVAR
      SHIFT = 1.D0
      CALL RANDPOS(BRAVAIS,BRINV,ISEED,SHIFT,ROLD,RNEW)
      RVAR(1:3,IVAR) = RNEW(1:3)
   ENDDO
   WRITE(21,*) 'Create file positions_var_old.dat if you wish to start from pre-given positions'
   WRITE(*,*)  'Create file positions_var_old.dat if you wish to start from pre-given positions'
ENDIF

OPEN(67,FILE='init_pos.dat')
WRITE(67,*) '# INITIAL POSITIONS OF EMPTY SPHERES:'
DO IVAR = 1,NVAR
   WRITE(67,FMT='(3E16.8)') (RVAR(IX,IVAR),IX=1,3)
ENDDO
CLOSE(67)



RETURN
END SUBROUTINE READPOS

!========================================================================

SUBROUTINE INVERTBR(BRAVAIS,BRINV,RECBV)
implicit none
! Invert bravais matrix, store it in BRINV
! Find bravais matrix in reciprocal space, sore in RECBV! Input
REAL*8 BRAVAIS(3,3)
! Output
REAL*8 BRINV(3,3),RECBV(3,3) 
! Internal 
REAL*8 UNIT(3,3)
INTEGER INFO,IPIV(3)
INTEGER I1,I2,I,J
REAL*8 SUM(3,3)
REAL*8 WORK(3)


!UNIT(:,:) = 0.D0
!UNIT(1,1) = 1.D0
!UNIT(2,2) = 1.D0
!UNIT(3,3) = 1.D0

BRINV(:,:) = BRAVAIS(:,:)

CALL DGETRF( 3, 3, BRINV, 3, IPIV, INFO )
CALL DGETRI( 3, BRINV, 3, IPIV, WORK, 3, INFO )
IF (INFO.NE.0) THEN
   WRITE(*,*) 'INFO FROM DGESV IN INVERTBR',INFO
   STOP 'ERROR'
ENDIF

CALL CROSPR(BRAVAIS(1,2),BRAVAIS(1,3),RECBV(1,1))
CALL CROSPR(BRAVAIS(1,3),BRAVAIS(1,1),RECBV(1,2))
CALL CROSPR(BRAVAIS(1,1),BRAVAIS(1,2),RECBV(1,3))

! BRINV and RECBV should differ only by 2*pi.

!SUM(:,:)=0.D0
!DO I=1,3
!DO J=1,3
!DO I1 = 1,3
!SUM(I,J) = SUM(I,J) +BRAVAIS(I,I1)*BRINV(I1,J)
!ENDDO
!ENDDO
!ENDDO


RETURN
END SUBROUTINE INVERTBR

!========================================================================

SUBROUTINE CART2INTERN(BRMAT,RPOS)
implicit none
! Transforms a vector RPOS between cartesian and internal coordinates
! Returns in the same array
! Cartesian --> Internal : BRMAT should be the Inverse Bravais Matrix BRINV(3,3)
! Internal --> Cartesian : BRMAT should be the Bravais Matrix BRAVAIS(3,3)
! Input
REAL*8 BRMAT(3,3) ! Either Bravais or inverse of Bravais matrix
! Input/Output
REAL*8 RPOS(3)
! Internal
REAL*8 RTMP(3)
INTEGER IX,IB

RTMP(:) = 0.D0
DO IB = 1,3
   DO IX = 1,3
      RTMP(IX) = RTMP(IX) + BRMAT(IX,IB) * RPOS(IB)  
   ENDDO
ENDDO
RPOS(:) = RTMP(:)

RETURN
END SUBROUTINE CART2INTERN



!========================================================================

SUBROUTINE METROPOLIS( &
ISEED,MCSTEPS,TEMPR,SHIFTMAX,BRAVAIS,BRINV,NSUPER,NFIX,NEXPANDFIX,RFIX,REXPANDFIX,TYPEFIX,NVAR,NEXPANDVAR,RHFIX,RHVAR, & !>
RVAR,REXPANDVAR,TYPEVAR,EINTRN,EEXTRN,ETOT, & ! X
EFLUCT,ACCEPTANCE)  ! <
implicit none
! Input
INTEGER ISEED,MCSTEPS,NSUPER
INTEGER NFIX,NVAR,NEXPANDFIX,NEXPANDVAR
REAL*8 SHIFTMAX,TEMPR
REAL*8 BRAVAIS(3,3),BRINV(3,3)
REAL*8 RFIX(3,NFIX),REXPANDFIX(3,NEXPANDFIX)
REAL*8 RHFIX(NFIX),RHVAR(NVAR) ! Hard-sphere radii
INTEGER TYPEFIX(NEXPANDFIX),TYPEVAR(NEXPANDVAR)
! Input/Output
REAL*8 RVAR(3,NVAR),REXPANDVAR(3,NEXPANDVAR)
REAL*8 EINTRN,EEXTRN,ETOT
! Output
REAL*8 ACCEPTANCE,EFLUCT
! Inside 
INTEGER ISTEP,IVAR,IVAR1
REAL*8 EOLD,EINTIOLD,EEXTIOLD,ENEW,EINTINEW,EEXTINEW,EDIFF,EXPEDIFF,EINTIDIFF,EEXTIDIFF,ETOT0,EAV,ESQAV,EFLUCTSQ
REAL*8 ETOT1,EINTRN1,EEXTRN1
REAL*8 XMETR,XRAND
REAL*8 ROLD(3),RNEW(3)
REAL*8,ALLOCATABLE :: DISTFIX(:),DISTVAR(:)
LOGICAL LACCEPT,LHSV,LHSVINT,LHSVEXT,LHSVDUM
EXTERNAL XRAND

ALLOCATE(  DISTFIX(NFIX), DISTVAR(NVAR) )


ETOT0 = 0.D0
EAV = 0.D0
ESQAV = 0.D0
ACCEPTANCE = 0.D0

DO ISTEP = 1,MCSTEPS

   DO IVAR1 = 1,NVAR ! Loop over all atoms
      ! Pick one of the moving atoms at random
      CALL RANDATOM(NVAR,ISEED,IVAR)
      ROLD(1:3) = RVAR(1:3,IVAR)
      ! Find distance of unshifted atom to all fixed and moving atoms
      CALL DISTUC3D(ROLD,NEXPANDFIX,REXPANDFIX,TYPEFIX,NFIX,DISTFIX)
      CALL DISTUC3D(ROLD,NEXPANDVAR,REXPANDVAR,TYPEVAR,NVAR,DISTVAR)
      ! Find energy of shifted atom to all fixed and moving atoms
      CALL EINT(DISTVAR,NVAR,IVAR,RHVAR,EINTIOLD,LHSVINT)
      CALL EEXT(DISTFIX,NFIX,RHFIX,RHVAR(IVAR),EEXTIOLD,LHSVEXT)  
      EOLD = EINTIOLD + EEXTIOLD

      ! Choose new position at random in unit cell; shift within +/- (1/2)fraction
      CALL RANDPOS(BRAVAIS,BRINV,ISEED,SHIFTMAX,ROLD,RNEW)

      ! Find distance of shifted atom to all fixed and moving atoms
      CALL DISTUC3D(RNEW,NEXPANDFIX,REXPANDFIX,TYPEFIX,NFIX,DISTFIX)
      CALL DISTUC3D(RNEW,NEXPANDVAR,REXPANDVAR,TYPEVAR,NVAR,DISTVAR)
      ! Find energy of shifted atom to all fixed and moving atoms
      CALL EINT(DISTVAR,NVAR,IVAR,RHVAR,EINTINEW,LHSVINT)
      CALL EEXT(DISTFIX,NFIX,RHFIX,RHVAR(IVAR),EEXTINEW,LHSVEXT)

      LHSV = LHSVINT.AND.LHSVEXT   ! True if there is a "hard sphere violation"

      ENEW = EINTINEW + EEXTINEW
      EDIFF = ENEW - EOLD
      EINTIDIFF = EINTINEW - EINTIOLD
      EEXTIDIFF = EEXTINEW - EEXTIOLD



      LACCEPT = .FALSE.
      IF (.NOT.LHSV) THEN ! If there is no hard-sphere violation then apply metropolis
         IF (EDIFF.LT.0.D0) THEN
            LACCEPT = .TRUE.
         ELSE
            IF (TEMPR.EQ.0.D0) THEN
               LACCEPT = .FALSE.
            ELSE
               XMETR = XRAND(ISEED)
               EXPEDIFF = DEXP(-EDIFF/TEMPR)
               IF (EXPEDIFF.GT.XMETR) LACCEPT = .TRUE.
            ENDIF
         ENDIF
      ENDIF

      IF (LACCEPT) THEN
         ! Place new position into array and recreate the supercell
         RVAR(1:3,IVAR) = RNEW(1:3)
         CALL EXPANDUC3D(BRAVAIS,NVAR,RVAR,NSUPER,IVAR,IVAR,REXPANDVAR,TYPEVAR)
         ACCEPTANCE = ACCEPTANCE + 1.D0
         ETOT = ETOT + EDIFF
         EINTRN = EINTRN + EINTIDIFF
         EEXTRN = EEXTRN + EEXTIDIFF
         ETOT0 = ETOT0 + EDIFF
    ENDIF

      EAV = EAV + ETOT0
      ESQAV = ESQAV + ETOT0**2


   ENDDO

ENDDO ! ISTEP = 1,MCSTEPS

EAV = EAV / DFLOAT(MCSTEPS*NVAR)
ESQAV = ESQAV / DFLOAT(MCSTEPS*NVAR)
EFLUCTSQ = ESQAV - EAV**2

IF (EFLUCTSQ.GE.0.D0) THEN
   EFLUCT = DSQRT(EFLUCTSQ)
ELSE
   EFLUCT = -DSQRT(DABS(EFLUCTSQ))
ENDIF


DEALLOCATE( DISTFIX,DISTVAR )

END SUBROUTINE METROPOLIS

!========================================================================

SUBROUTINE METROPOLIS2( &
ISEED,MCSTEPS,TEMPR,SHIFTMAX,BRAVAIS,BRINV,NSUPER,NFIX,NEXPANDFIX,RFIX,REXPANDFIX,TYPEFIX,NVAR,NEXPANDVAR, & !>
RVAR,REXPANDVAR,TYPEVAR,ETOT, & ! X
EFLUCT,ACCEPTANCE)  ! <
implicit none
! Input
INTEGER ISEED,MCSTEPS,NSUPER
INTEGER NFIX,NVAR,NEXPANDFIX,NEXPANDVAR
REAL*8 SHIFTMAX,TEMPR
REAL*8 BRAVAIS(3,3),BRINV(3,3)
REAL*8 RFIX(3,NFIX),REXPANDFIX(3,NEXPANDFIX)
INTEGER TYPEFIX(NEXPANDFIX),TYPEVAR(NEXPANDVAR)
! Input/Output
REAL*8 RVAR(3,NVAR),REXPANDVAR(3,NEXPANDVAR)
REAL*8 EINTRN,EEXTRN,ETOT
! Output
REAL*8 ACCEPTANCE,EFLUCT
! Inside 
INTEGER ISTEP,IVAR,IVAR1
REAL*8 EOLD,EINTIOLD,EEXTIOLD,ENEW,EINTINEW,EEXTINEW,EDIFF,EXPEDIFF,EINTIDIFF,EEXTIDIFF,ETOT0,EAV,ESQAV,EFLUCTSQ
REAL*8 VOLFIL
REAL*8 XMETR,XRAND
REAL*8 ROLD(3),RNEW(3)
REAL*8,ALLOCATABLE :: DISTFIX(:),DISTVAR(:)
LOGICAL LACCEPT
EXTERNAL XRAND


ETOT0 = 0.D0
EAV = 0.D0
ESQAV = 0.D0
ACCEPTANCE = 0.D0

! Evaluate volume filling fraction
CALL VOLUMEFIL2(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR, &  ! >
     VOLFIL)                                                    ! <
EOLD = -VOLFIL


DO ISTEP = 1,MCSTEPS

   DO IVAR1 = 1,NVAR ! Loop over all atoms
      ! Pick one of the moving atoms at random
      CALL RANDATOM(NVAR,ISEED,IVAR)
      ROLD(1:3) = RVAR(1:3,IVAR)
      ! Choose new position at random in unit cell; shift within +/- (1/2)fraction
      CALL RANDPOS(BRAVAIS,BRINV,ISEED,SHIFTMAX,ROLD,RNEW)
      RVAR(1:3,IVAR) = RNEW(1:3)
      CALL EXPANDUC3D(BRAVAIS,NVAR,RVAR,NSUPER,IVAR,IVAR,REXPANDVAR,TYPEVAR)
      CALL VOLUMEFIL2(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR, &  ! >
                     VOLFIL)                                                    ! <


      ENEW = -VOLFIL
      EDIFF = ENEW - EOLD


      LACCEPT = .FALSE.
      IF (EDIFF.LT.0) THEN
         LACCEPT = .TRUE.
      ELSE
         IF (TEMPR.EQ.0.D0) THEN
            LACCEPT = .FALSE.
         ELSE
            XMETR = XRAND(ISEED)
            EXPEDIFF = DEXP(-EDIFF/TEMPR)
            IF (EXPEDIFF.GT.XMETR) LACCEPT = .TRUE.
         ENDIF
      ENDIF

      IF (LACCEPT) THEN
         ACCEPTANCE = ACCEPTANCE + 1.D0
         EOLD = ENEW
         ETOT0 = ETOT0 + EDIFF
      ELSE ! Restore old atom
         ! Place old position into array and recreate the supercell
         RVAR(1:3,IVAR) = ROLD(1:3)
         CALL EXPANDUC3D(BRAVAIS,NVAR,RVAR,NSUPER,IVAR,IVAR,REXPANDVAR,TYPEVAR)
      ENDIF

      EAV = EAV + ETOT0
      ESQAV = ESQAV + ETOT0**2

   ENDDO

ENDDO ! ISTEP = 1,MCSTEPS


EAV = EAV / DFLOAT(MCSTEPS*NVAR)
ESQAV = ESQAV / DFLOAT(MCSTEPS*NVAR)
EFLUCTSQ = ESQAV - EAV**2

IF (EFLUCTSQ.GE.0.D0) THEN
   EFLUCT = DSQRT(EFLUCTSQ)
ELSE
   EFLUCT = -DSQRT(DABS(EFLUCTSQ))
ENDIF


END SUBROUTINE METROPOLIS2
!========================================================================

SUBROUTINE VOLUMEFIL2(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR, &  ! >
                     VOLFIL)                                                          ! <
implicit none
! Volume filling
! Input:
REAL*8 BRAVAIS(3,3)
REAL*8 RFIX(3,*),RVAR(3,*),REXPANDFIX(3,*),REXPANDVAR(3,*)
INTEGER NFIX,NVAR,NEXPANDFIX,NEXPANDVAR
! Output:
REAL*8 VOLMT,VOLUME,VOLFIL ! Volume of MT-spheres, of unit cell, and volume filling fraction
! Internal 
INTEGER IVAR,IFIX,IVAR1,IAT,IAT1,NTOT,NEXPANDTOT
REAL*8 DD
REAL*8 RVEC(3)
REAL*8,ALLOCATABLE :: RMT(:),RAT(:,:),REXPAND(:,:)
REAL*8 PI

PI = 4.D0 * DATAN(1.D0)

NTOT = NVAR + NFIX
NEXPANDTOT = NEXPANDVAR + NEXPANDFIX
ALLOCATE(RMT(NTOT),RAT(3,NTOT),REXPAND(3,NEXPANDTOT))

!Pass all atoms in one array
DO IAT = 1,NVAR
   RAT(1:3,IAT) = RVAR(1:3,IAT)
ENDDO
DO IFIX = 1,NFIX
   IAT = NVAR + IFIX
   RAT(1:3,IAT) = RFIX(1:3,IFIX)
ENDDO

DO IAT = 1,NEXPANDVAR
   REXPAND(1:3,IAT) = REXPANDVAR(1:3,IAT)
ENDDO
DO IFIX = 1,NEXPANDFIX
   IAT = NEXPANDVAR + IFIX
   REXPAND(1:3,IAT) = REXPANDFIX(1:3,IFIX)
ENDDO

RMT(1:NTOT) = 1.D100

DO IAT = 1,NTOT
   DO IAT1 = 1,NEXPANDTOT
      DD = (RAT(1,IAT)-REXPAND(1,IAT1))**2+(RAT(2,IAT)-REXPAND(2,IAT1))**2+(RAT(3,IAT)-REXPAND(3,IAT1))**2
      IF (DD.LT.RMT(IAT).AND.DD.GT.0.D0) RMT(IAT) = DD
   ENDDO
   RMT(IAT) = 0.5D0 * DSQRT(RMT(IAT))
ENDDO

VOLMT = 0.D0
DO IAT = 1,NTOT
   VOLMT = VOLMT + RMT(IAT)**3
ENDDO
VOLMT = 4.D0 * PI * VOLMT / 3.D0

VOLUME =  BRAVAIS(1,1) * (BRAVAIS(2,2)*BRAVAIS(3,3)-BRAVAIS(2,3)*BRAVAIS(3,2)) &
        + BRAVAIS(2,1) * (BRAVAIS(3,2)*BRAVAIS(1,3)-BRAVAIS(1,2)*BRAVAIS(3,3)) &
        + BRAVAIS(3,1) * (BRAVAIS(1,2)*BRAVAIS(2,3)-BRAVAIS(2,2)*BRAVAIS(1,3)) 



VOLUME = DABS(VOLUME)
VOLFIL = VOLMT / VOLUME

DEALLOCATE(RMT,RAT,REXPAND)

!WRITE(*,*) 'Unit Cell Volume :',VOLUME
!WRITE(*,*) 'Muffin Tin Volume:',VOLMT
!WRITE(*,*) 'Volume fill fraction:',VOLFIL


RETURN
END SUBROUTINE VOLUMEFIL2

!========================================================================

SUBROUTINE CROSPR(X,Y,Z)
implicit none
! ------------------------------------------------------------------------
!     CROSP COMPUTES THE CROSS PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------
REAL*8           X(*), Y(*), Z(*)
Z(1)=X(2)*Y(3)-X(3)*Y(2)
Z(2)=X(3)*Y(1)-X(1)*Y(3)
Z(3)=X(1)*Y(2)-X(2)*Y(1)
RETURN
END SUBROUTINE CROSPR

!========================================================================

SUBROUTINE RATIONALBASIS(BRAVAIS,BRINV,NFIX,RFIX)
implicit none
! Bring basis atoms into the primitive cell by Bravais translation
!Input:
REAL*8 BRAVAIS(3,3),BRINV(3,3)
INTEGER NFIX 
! Input/output
REAL*8 RFIX(3,NFIX)
!Internal
INTEGER IFIX,IB(3)
REAL*8 RVEC(3)

DO IFIX = 1,NFIX
   RVEC(1:3) = RFIX(1:3,IFIX)
   CALL CART2INTERN(BRINV,RVEC) ! convert to internal coords
   IB(1:3) = FLOOR(RVEC(1:3))   ! Inreger part
   RVEC(1:3) = RVEC(1:3) - DFLOAT(IB(1:3)) ! Shift to primitive cell
   CALL CART2INTERN(BRAVAIS,RVEC) ! Back to cartesian
   RFIX(1:3,IFIX) = RVEC(1:3)
ENDDO

RETURN
END SUBROUTINE RATIONALBASIS

!========================================================================

SUBROUTINE VOLUMEFILH(BRAVAIS,NFIX,RFIX,NEXPANDFIX,REXPANDFIX,TYPEFIX,NVAR,RVAR,NEXPANDVAR,REXPANDVAR,TYPEVAR, &  ! >
                     RMTVAR,RMTFIX,RHVAR,RHFIX,VOLMT,VOLUME,VOLFIL)                                                          ! <
implicit none
! Inflate spheres at given positions so that they are touching at the end
! Input:
REAL*8 BRAVAIS(3,3)
REAL*8 RFIX(3,*),RVAR(3,*),REXPANDFIX(3,*),REXPANDVAR(3,*),RHVAR(*),RHFIX(*)
INTEGER NFIX,NVAR,NEXPANDFIX,NEXPANDVAR
INTEGER TYPEFIX(*),TYPEVAR(*)
! Output:
REAL*8 VOLMT,VOLUME,VOLFIL ! Volume of MT-spheres, of unit cell, and volume filling fraction
REAL*8 RMTVAR(*),RMTFIX(*) ! Muffin tin radii of moving and fixed atoms
! Internal 
INTEGER IVAR,IFIX,IVAR1,IAT,IAT1,NTOT,NEXPANDTOT
REAL*8 DD
REAL*8 RVEC(3)
REAL*8,ALLOCATABLE :: RMT(:),RAT(:,:),REXPAND(:,:),RH(:),RMTNEW(:)
INTEGER,ALLOCATABLE :: ATYPE(:)
LOGICAL,ALLOCATABLE:: LINFL(:)
LOGICAL LTOUCH
REAL*8 TOLTOUCH ! Tolerance for touching mt spheres
REAL*8 DSTEP ! step size while inflating
INTEGER ISTEP,NSTEP ! Inflation steps
INTEGER NINFL ! How many spheres are inflated
DATA TOLTOUCH /1.D-5/
DATA DSTEP /1.D-5/
DATA NSTEP /10000/
REAL*8 PI

PI = 4.D0 * DATAN(1.D0)

NTOT = NVAR + NFIX
NEXPANDTOT = NEXPANDVAR + NEXPANDFIX
ALLOCATE(RMT(NTOT),RAT(3,NTOT),REXPAND(3,NEXPANDTOT),ATYPE(NEXPANDTOT),RH(NTOT),LINFL(NTOT),RMTNEW(NTOT))

! Pass all atoms in one array, first moving atoms, then fixed atoms
DO IAT = 1,NVAR
   RAT(1:3,IAT) = RVAR(1:3,IAT)   ! Atom positions
   RH(IAT) = RHVAR(IAT)           ! Hard sphere radii
ENDDO
DO IFIX = 1,NFIX
   IAT = NVAR + IFIX
   RAT(1:3,IAT) = RFIX(1:3,IFIX)  ! Atom positions
   RH(IAT) = RHFIX(IFIX)          ! Hard sphere radii
ENDDO

DO IAT = 1,NEXPANDVAR
   REXPAND(1:3,IAT) = REXPANDVAR(1:3,IAT)
   ATYPE(IAT) = TYPEVAR(IAT)
ENDDO
DO IFIX = 1,NEXPANDFIX
   IAT = NEXPANDVAR + IFIX
   REXPAND(1:3,IAT) = REXPANDFIX(1:3,IFIX)
   ATYPE(IAT) = NVAR + TYPEFIX(IFIX)
ENDDO


RMT(1:NTOT) = 1.D100

! Find initial MT radii by bisection of distance to neighbours.
! Take care not to enter hard-sphere radius (RH) of neighbours.
DO IAT = 1,NTOT ! loop over atoms in primitive cell
   DO IAT1 = 1,NEXPANDTOT ! Loop over atoms in expanded cell (supercell)
      DD = (RAT(1,IAT)-REXPAND(1,IAT1))**2+(RAT(2,IAT)-REXPAND(2,IAT1))**2+(RAT(3,IAT)-REXPAND(3,IAT1))**2
      DD = DSQRT(DD)
      DD = MIN(0.5D0*DD,DD-RH(ATYPE(IAT1)))
      IF (DD.LT.RMT(IAT).AND.DD.GT.1.D-10) RMT(IAT) = DD
   ENDDO
   ! In case RMT became smaller than hard-sphere radius, return it to hard-sphere radius
   RMT(IAT) = MAX(RMT(IAT),RH(IAT))
ENDDO

! Find which atoms can be further inflated
NINFL = 0
DO IAT = 1,NTOT
   LINFL(IAT) = .TRUE. ! LINFL is true if the atom can be further inflated
   DO IAT1 = 1,NEXPANDTOT ! Loop over atoms in expanded cell (supercell)
      DD = (RAT(1,IAT)-REXPAND(1,IAT1))**2+(RAT(2,IAT)-REXPAND(2,IAT1))**2+(RAT(3,IAT)-REXPAND(3,IAT1))**2 
      DD = DSQRT(DD)
      IF (DD.LT.1.D-10) EXIT ! take care that atom is not compared to itself
      LTOUCH = ( RMT(IAT)+RMT(ATYPE(IAT1))+TOLTOUCH.GT.DD ) 
      LINFL(IAT) = LINFL(IAT).AND.LTOUCH  ! LINFL becomes false if the MT spheres touch
      IF (.NOT.LINFL(IAT)) EXIT ! Break loop if at least one neighbor is touching
   ENDDO
   IF (LINFL(IAT)) NINFL = NINFL + 1
ENDDO
WRITE(*,*) 'VOLUMEFILH: Found inflatable spheres NINFL=',NINFL

! Perform futher inflation 
RMTNEW(1:NTOT) = RMT(1:NTOT)
DO ISTEP = 1,NSTEP ! Loop over inflation steps

   ! Try to inflate all atoms simultanteously, not one-by-one, for symmetry reasons.
   DO IAT = 1,NTOT
      IF (LINFL(IAT)) THEN  ! Only inflate atoms that do not touch their neighbours
         RMTNEW(IAT) = RMT(IAT) + DSTEP
      ENDIF
   ENDDO

   ! Now atoms have new trial RMT, test if they can be accepted.
   NINFL = 0
   DO IAT = 1,NTOT
      IF (LINFL(IAT)) THEN
         DO IAT1 = 1,NEXPANDTOT
            DD = (RAT(1,IAT)-REXPAND(1,IAT1))**2+(RAT(2,IAT)-REXPAND(2,IAT1))**2+(RAT(3,IAT)-REXPAND(3,IAT1))**2  
            DD = DSQRT(DD)
            IF (DD.LT.1.D-10) EXIT ! take care that atom is not compared to itself
            LTOUCH = ( RMTNEW(IAT)+RMTNEW(ATYPE(IAT1))+TOLTOUCH.GT.DD ) 
            LINFL(IAT) = LINFL(IAT).AND.LTOUCH  ! LINFL becomes false if the MT spheres touch
            IF (.NOT.LINFL(IAT)) EXIT ! Break loop if at least one neighbor is touching
         ENDDO
         IF (LINFL(IAT)) RMT(IAT) = RMTNEW(IAT)  ! Accept inflation if radius did not hit neighbour
      ENDIF
      IF (LINFL(IAT)) NINFL = NINFL + 1
   ENDDO

ENDDO
WRITE(*,*) 'VOLUMEFILH: Remaining inflatable spheres NINFL=',NINFL


VOLMT = 0.D0
DO IAT = 1,NTOT
   VOLMT = VOLMT + RMT(IAT)**3
ENDDO
VOLMT = 4.D0 * PI * VOLMT / 3.D0

VOLUME =  BRAVAIS(1,1) * (BRAVAIS(2,2)*BRAVAIS(3,3)-BRAVAIS(2,3)*BRAVAIS(3,2)) &
        + BRAVAIS(2,1) * (BRAVAIS(3,2)*BRAVAIS(1,3)-BRAVAIS(1,2)*BRAVAIS(3,3)) &
        + BRAVAIS(3,1) * (BRAVAIS(1,2)*BRAVAIS(2,3)-BRAVAIS(2,2)*BRAVAIS(1,3)) 



VOLUME = DABS(VOLUME)
VOLFIL = VOLMT / VOLUME

DO IAT = 1,NVAR
   RMTVAR(IAT) = RMT(IAT)
ENDDO
DO IAT = NVAR + 1,NTOT
   IFIX = IAT - NVAR
   RMTFIX(IFIX) = RMT(IAT)
ENDDO


!WRITE(*,*) 'Unit Cell Volume :',VOLUME
!WRITE(*,*) 'Muffin Tin Volume:',VOLMT
!WRITE(*,*) 'Volume fill fraction:',VOLFIL

DEALLOCATE(RMT,RAT,REXPAND,ATYPE,RH)


RETURN
END SUBROUTINE VOLUMEFILH

