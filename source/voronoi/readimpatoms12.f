       SUBROUTINE READIMPATOMS12(
     >     ALATC,LCARTESIAN,
     <     NUMIMP,RIMPURITY,NKILLATOM,RKILL,DXIMP,DYIMP,DZIMP,
     <     RMTIMP,WEIGHT,ZIMP,LCARTESIMP)   
!
! Read in atomic positions for preparation for the impurity program
! from the "inputcard", if option IMPURITY is used.
! Contents should look like this:
!
! IMPINFO            (keyword)
! N      (=Number of sites to read in (at correct positions))
! 1 x(1)   y(1)    z(1)   rmtcore(1)  weight(1)  Z(1)   ! Large Z means atomic number
! 2 x(2)   y(2)    z(2)   rmtcore(2)  weight(2)  Z(2)  
! 3 x(3)   y(3)    z(3)   rmtcore(3)  weight(3)  Z(3)
!  .............
! N Z(N) x(N)   y(N)    z(N)   killatom(N)   weight(N)
! M    (= Number of unshifted positions among the above positions )
! index(1) xold(1)   yold(1)    zold(1)
!  .............
! index(N) xold(N)   yold(N)    zold(N)
! K (= Number of positions to be killed)
! xkill(1)  ykill(1) zkill(1)
!  .......................
! xkill(N)  ykill(N) zkill(N)
!
! (Put M=0 if there are no atoms to be shifted,
!  K=0 if there are no atoms to be killed.)
! 
! The index(i) shows which atom among the initial read-in atoms
! is to be shifted to new position (xnew,ynew,znew).
!
! The index named killatom(i) is >0 if the corresponding host-atom 
! should be excluded in the impurity calculation. Such atoms are
! excluded from the atom clusters for the Voronoi cells.
!
      implicit none
c#@# KKRtags: VORONOI input-output KKRimp
      INCLUDE 'inc.geometry'
! Input:
      REAL*8 ALATC  ! Lattice parameter
      LOGICAL LCARTESIAN
! Output:
      INTEGER NUMIMP,NKILLATOM  ! Number of impurity atoms to keep and killed atoms
      REAL*8 RIMPURITY(3,NIMPD), RKILL(3,NIMPD) ! Corresponding coordinates
      REAL*8 DXIMP(NIMPD),DYIMP(NIMPD),DZIMP(NIMPD) ! Shift of atom to new position
      REAL*8 WEIGHT(NIMPD) ! Weight for the Voronoi construction
      REAL*8 ZIMP(NIMPD)   ! Impurity atomic number
      REAL*8 RMTIMP(NIMPD)
      LOGICAL LCARTESIMP ! Imp. potitions in cartesian (true) or internal (false) coords.
! Local:
      REAL*8  R0(3,NIMPD),R1(3,NIMPD),ZREAD(NIMPD)
      REAL*8  RMTREAD(NIMPD),WREAD(NIMPD)
      INTEGER KILLATOM(NIMPD)  ! (0/1) Host atoms to be removed in impurity calculation

      INTEGER INDEX,IREAD,NREAD1,NREAD2,IAT,IX,ILINE,IER
      LOGICAL LSHIFT
      CHARACTER*256 UIO

! Array R0 containd coordinates of unshifted positions,
! array R1 containd coordinates of shifted positions,
! DX,DY,DZ are the shifting vectors R1-R0.
      WRITE(*,*) 'Entering READIMPATOMS12'

      LCARTESIMP = LCARTESIAN
      CALL IoInput('CARTESIMP       ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) LCARTESIMP
      WRITE(*,*) 'readimpatoms: CARTESIMP=',LCARTESIMP


      ILINE = 1
      CALL IoInput('IMPINFO         ',UIO,ILINE,7,IER)
      ILINE = ILINE + 1
      READ (UNIT=UIO,FMT=*) NUMIMP
     
      IF (NUMIMP.GT.NIMPD) STOP 'READIMPATOMS12: NUMIMP.GT.NIMPD'

      WRITE(6,*) 'Reading impurity atoms from inputcard'
      WRITE(6,*) 
     & 'INDEX      X          Y          Z       RMT    WEIGHT    Z'

! Read in impurity-atom positions from file. Weight should be wished MT radius
      DO IAT = 1,NUMIMP

         CALL IoInput('IMPINFO         ',UIO,ILINE,7,IER)
         ILINE = ILINE + 1
         READ(UNIT=UIO,FMT=*)  INDEX,(RIMPURITY(IX,IAT),IX=1,3),
     &                         RMTIMP(IAT),WEIGHT(IAT),ZIMP(IAT)

         WRITE(6,1015)         INDEX,(RIMPURITY(IX,IAT),IX=1,3),
     &                         RMTIMP(IAT),WEIGHT(IAT),ZIMP(IAT)
      END DO 
! Read in unshifted positions of only the atoms that are to be shifted.
! Mapping to previous positions is made by index.
      R0(:,:) = RIMPURITY(:,:)
      WRITE(6,*) 'Unshifted positions:'
      CALL IoInput('IMPINFO         ',UIO,ILINE,7,IER)
      ILINE = ILINE + 1
      READ(UNIT=UIO,FMT=*) NREAD2
      IF (NREAD2.GT.NUMIMP) STOP 'READIMPATOMS12: NREAD2.GT.NUMIMP'

      DO IREAD = 1,NREAD2
         CALL IoInput('IMPINFO         ',UIO,ILINE,7,IER)
         ILINE = ILINE + 1
         READ(UNIT=UIO,FMT=*) INDEX,(R0(IX,INDEX),IX=1,3)

         WRITE(6,1010) INDEX,(R0(IX,INDEX),IX=1,3)
      END DO 

! Read in killed positions

      WRITE(6,*) 'Killed positions:'

      CALL IoInput('IMPINFO         ',UIO,ILINE,7,IER)
      ILINE = ILINE + 1
      READ(UNIT=UIO,FMT=*) NKILLATOM
      IF (NKILLATOM.GT.NIMPD) STOP 'READIMPATOMS12: NKILLATOM.GT.NIMPD'

      DO IREAD = 1,NKILLATOM
         CALL IoInput('IMPINFO         ',UIO,ILINE,7,IER)
         ILINE = ILINE + 1
         READ (UNIT=UIO,FMT=*) 
     &        RKILL(1,IREAD),RKILL(2,IREAD),RKILL(3,IREAD)

         WRITE(6,1010) IREAD,(R0(IX,INDEX),IX=1,3)
      ENDDO

! Re-define impurity weights: should become square of wished MT radius
! in units of latt. constant. If the command <MTWAU> was given in the
! inputcard, then it is assumed that the host atom weights are in 
! atomic units, and the same is true for the impurity weights.
      CALL IoInput('<MTWAU>         ',UIO,1,7,IER)
      IF (IER.EQ.0)  WEIGHT(1:NUMIMP) = WEIGHT(1:NUMIMP)/ALATC
      WEIGHT(1:NUMIMP) = WEIGHT(1:NUMIMP)**2

! The Voronoi cell will be centered at R0, but the shape
! will be expanded around RIMPURITY = R0 + DX,Y,Z.

      WRITE(6,*) 'Impurity atoms read in from inputcard'

! Now put unshifted positions are in array RIMPURITY, extra shift in
! arrays DXIMP,DYIMP,DZIMP

      DXIMP(1:NUMIMP) = RIMPURITY(1,1:NUMIMP) - R0(1,1:NUMIMP)
      DYIMP(1:NUMIMP) = RIMPURITY(2,1:NUMIMP) - R0(2,1:NUMIMP)
      DZIMP(1:NUMIMP) = RIMPURITY(3,1:NUMIMP) - R0(3,1:NUMIMP)
      RIMPURITY(:,:) = R0(:,:)


      WRITE(*,*) 'READIMPATOMS12:'
      WRITE(*,*) 'Found ',NUMIMP,' impurities and ',NKILLATOM,
     &           ' sites to be killed.'

      WRITE(*,*) 'Unshifted impurity positions and extra shift:'
      DO IAT = 1,NUMIMP
         WRITE(*,1015) IAT,(RIMPURITY(IX,IAT),IX=1,3),
     &                 DXIMP(IAT),DYIMP(IAT),DZIMP(IAT)
      ENDDO
      IF (NKILLATOM.GT.0) WRITE(*,*) 'Killed-atom positions:'
      DO IAT = 1,NKILLATOM
         WRITE(*,1015) IAT,(RKILL(IX,IAT),IX=1,3)
      ENDDO

      
      WRITE(*,*) 'Exiting READIMPATOMS12'
        
 1000 FORMAT(A93)
 1005 FORMAT(A66)
 1010 FORMAT(I5,3F12.8)
 1015 FORMAT(I5,6F12.8)
       END





