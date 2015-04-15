      program magdip
c
c     Index:
c
c     1) Purpose of the program
c     2) Input data and their physical meaning
c     3) Physical meaning of another variables
c     4) Compilation and execution of the program
c
c
c     1) Purpose of the program
c
c     This program makes three types of calculations: 
c
c     a) Calculation of the magnetostatic dipolar energy per unit cell of a lattice 
c     of parallel magnetic dipoles for ndir directions of the magnetic dipoles and
c     the magnetostatic dipolar anisotropy energy, i.e. the difference between the 
c     magnetostatic dipolar energies of these ndir directions (first - second), 
c     ..., (first - ndir)
c
c     b) Calculation of the magnetostatic dipolar field per unit cell at npoints 
c     specific points r, due to a lattice of parallel or non-parallel magnetic 
c     dipoles. These points must not coincide with the positions of the basis 
c     atoms in the unit cell. The next type of calculation deals with these points
c
c     c) Calculation of the magnetostatic dipolar energy per unit cell of a lattice 
c     of parallel or non-parallel magnetic dipoles and the corresponding magnetic
c     dipolar field at the positions of the basis atoms in the unit cell
c
c
c     The isotropic Fermi contact contribution is not included in the magnetostatic
c     dipolar energies and magnetic dipolar fields. In any case, this contribution
c     must not be included in the calculations of type b)
c
c
c     2) Input data and their physical meaning
c
c     Line 1:
c     natom: total number of atoms in the unit cell
c     magmod(i) (i=1,natom): the modulus of the magnetic moment of the
c                             atom i in units of the Bohr magneton
c
c     Line 2:
c     ndir: total number of directions. A direction means here the set of 
c           Euler angles that defined the direction of the magnetic moments 
c           of a lattice of parallel magnetic dipoles.
c
c     Lines 3 to 2+ndir:
c     alphadeg(id), betadeg(id), gammadeg(id) (id=1,ndir): Euler angles in
c                        degrees of the id direction, not in radians, to
c                        rotate the frame of reference from the one of
c                        the unit cell, (x,y,z), to that of the magnetic
c                        moment of the atoms. When option=0, the magnetic
c                        moments of all the atoms are oriented in the same
c                        direction defined by these Euler angles.
c                        Examples: 
c                        Euler angles of a moment || x axis:  0,90,0
c                          "     "    "  "   "    || y axis: 90,90,0
c                          "     "    "  "   "    || z axis:  0, 0,0
c     
c     The lines 3+ndir to 7+ndir+natom are used by the madelung subroutine
c
c     Line 3+ndir:
c     bound: precision of the matrix elements smat. It is also approximately
c            the precision of the anisotropy energy in eV
c
c     Line 4+ndir:
c     alat: lattice parameter in bohrs
c     rmax: maximum radius in units of alat of the shells of atoms in the
c           real space included in the calculation of the magnetostatic dipolar 
c           energy
c     gmax: maximum radius in units of 1/alat of the shells of atoms in the
c           reciprocal space included in the calculation of the magnetostatic
c           dipolar energy
c
c     Lines 5+ndir to 7+ndir:
c     br(1,i), br(2,i), br(3,i) (i=1,3): coordinates in a cartesian reference
c              frame of the i primitive vector of the lattice, in units of alat.
c              br(a,i), a = 1,2,3 --> x,y,z
c
c     Lines 8+ndir to 7+ndir+natom:
c     qi(1,i), qi(2,i), qi(3,i) (i=1,natom): coordinates in a cartesian reference
c              frame of the i basis vector of the lattice, in units of alat.
c              qi(a,i), a = 1,2,3 --> x,y,z
c
c     Line 8+ndir+natom:
c     option: type of calculation. option=0,1,2 means type a), b) and c),
c             respectively
c
c     The following lines only exist and must be read when option>0
c
c     Lines 9+ndir+natom to 8+ndir+2*natom:
c     magmom(i,id) (i=1,natom, id=1,3): cartesian coordinates of the unitary 
c                                       vector of the magnetic moment of 
c                                       atom i. id=1,2,3 --> x,y,z
c
c     Lines 9+ndir+2*natom to 8+ndir+3*natom:
c     alphadegmom(i), betadegmom(i), gammadegmom(i) (i=1,natom): Euler angles
c                        in degrees, not in radians, to rotate the frame of 
c                        reference from the one of the unit cell, (x,y,z), to 
c                        that of the magnetic moment of every particular atom i 
c                        of the unit cell. This is the most general case: the 
c                        magnetic moments of all the atoms of the unit cell are 
c                        not oriented in the same direction, and the orientation 
c                        is defined by these Euler angles
c
c     The following lines only exist and must be read when option=1      
c
c     Line 9+ndir+3*natom:
c     npoints: number of points where the magnetic dipolar field will be calculated
c
c     Lines 10+ndir+3*natom to 9+ndir+3*natom+npoints:
c     r(1,i), r(2,i), r(3,i) (i=1,natom): the cartesian coordinates of the point i 
c                                         in units of the lattice parameter alat. 
c                                         r(a,i), a = 1,2,3 --> x,y,z
c
c
c     3) Physical meaning of another variables
c
c     sum(i,j): multiplied by magmod(i)*magmod(j)*c, this is the magnetostatic dipolar
c               energy coming from the interaction between a lattice of magnetic
c               moments of the atoms of type i with only one magnetic moment of
c               the atom of type j in the unit cell, parallel to the magnetic moments
c               of the atoms of type i
c
c     magdipenergy(id): Magnetostatic dipolar energy of the lattice of magnetic
c                       dipoles for the id direction in Rydbergs
c
c     magdipanienergy: Magnetostatic dipolar anisotropy energy of the lattice of magnetic
c                      dipoles in Rydbergs. It is the difference between the mag.
c                      dip. energy for the first and another direction
c
c     magfield(i,j,m): Component m of the magnetic dipolar field per unit cell at the 
c                      point where is the atom of type j, due to the atoms of type i.
c                      m=1,2,3 --> x,y,z
c
c     magfieldtotal(j,m): Component m of the total magnetic dipolar field per unit 
c                         cell at the point where is the atom of type j. 
c                         m=1,2,3 --> x,y,z
c
c     natomd: maximum number of atoms per unit cell plus maximum number of points where
c             the magnetic dipolar field will be calculated
c
c     nmaxpoints: maximum number of points where the magnetic dipolar field will be 
c                 calculated
c
c     mmax: maximum magnetic quantum number
c
c     nmaxdir: maximum number of directions of the magnetic moments when option=0
c
c
c     4) Compilation and execution of the program
c
c     To compile the program, type f77 -i8 -r8 -o magdip.exe magdip.f, where magdip.exe 
c     and magdip.f are the executable and the source file, respectively.
c
c     To execute the program, type magdip.exe < input > output
c
c
c     Ivan Cabria      November 15th 2000
c     Forschungszentrum Juelich
c     D-52425 Juelich
c     Germany
c
c
c     Declaration of parameters
c
      implicit none

      integer natypmax,natomd,lmaxd,mmax,nkmmax,nmaxdir,nmaxpoints
      parameter ( natypmax=5, natomd=80, lmaxd=16, mmax=5, nkmmax=9)
      parameter ( nmaxdir=4, nmaxpoints = 10 )

      real pi,ry,muB,eV,clightryd,c,d,tesla
      parameter ( pi = 3.141592653589793238462643d0 )
c                      eV/Rydberg
      parameter ( ry = 13.6058d0 )
c                      Joule/eV
      parameter ( eV = 1.60219d-19 )
c
c     In atomic units: e = hbar = m = 1
c     In atomic Rydbergs units: e = sqrt(2.d0), hbar = 1, 2m = 1
c
c     1 Borh magneton = e*hbar/(2m) = 1.d0/2.d0   in atomic units 
c                                   = sqrt(2.d0)  in atomic Rydbergs units
c                                   = 9.27408d-24 Joule/Tesla
c 
c                       Joule/Tesla
      parameter ( muB = 9.27408d-24 )
c                         Tesla/Rydberg
      parameter ( tesla = ry*eV*sqrt(2.d0)/muB )
c
c     Speed of light in atomic Rydbergs units
      parameter ( clightryd = 2.d0*137.036d0 )
c
c     A constant used in the program
      parameter ( c = -dsqrt(16.d0*pi/5.d0)/(clightryd**2) )
c
c     Declaration of variables
c
      integer i,j,m,natom,ndir,npoints,id,option,lmax

      real    smat(natomd,natomd,(2*lmaxd+1)**2),br(3,3),qi(3,natomd)
      real    bound,alat,rmax,gmax
      real    alphadeg(nmaxdir),betadeg(nmaxdir),gammadeg(nmaxdir)
      real    alphadegmom(natomd),betadegmom(natomd),gammadegmom(natomd)
      real    magmod(natomd),fact(0:100),r(3,nmaxpoints)
      real    magdipenergy(nmaxdir),magdipanienergy
      real    magfield(natomd,nmaxpoints,3),magfieldtotal(nmaxpoints,3)
      real    magmom(natomd,3),sum(natomd,natomd),sumtotal

      complex rot(nkmmax,nkmmax)
c
c     Intrinsic functions
c
      intrinsic sqrt
c
c     Initialization of arrays
c
      call rinit(natomd,magmod)
      call rinit(natomd*3,magmom)
      call rinit(nmaxpoints*3,magfieldtotal)
      call rinit(natomd*nmaxpoints*3,magfield)
      call rinit(natomd*natomd*(2*lmaxd+1)**2,smat)
      call rinit(3*3,br)
      call rinit(3*natomd,qi)
      call rinit(3*nmaxpoints,r)
      call rinit(nmaxdir,magdipenergy)
      call rinit(natomd,alphadegmom)
      call rinit(natomd,betadegmom)
      call rinit(natomd,gammadegmom)
      call rinit(nmaxdir,alphadeg)
      call rinit(nmaxdir,betadeg)
      call rinit(nmaxdir,gammadeg)
      call rinit(natom*natomd,sum)

      call cinit(nkmmax*nkmmax,rot)
c
c     Reading input data
c
      read(5,*) natom,(magmod(i), i=1,natom)
      read(5,*) ndir
      do i = 1,ndir
         read(5,*) alphadeg(i),betadeg(i),gammadeg(i)
         print*,alphadeg(i),betadeg(i),gammadeg(i)
      end do

      read(5,*) bound
      read(5,*) alat,rmax,gmax
      do i = 1,3
         read(5,*) br(1,i),br(2,i),br(3,i)
      end do
      do i = 1,natom
         read(5,*) qi(1,i),qi(2,i),qi(3,i)
      end do

      read(5,*) option
      if ((option.ne.0).and.(option.ne.1).and.(option.ne.2)) then
         write(6,100)
         stop
      end if
      if (option.gt.0) then
         do i = 1,natom
            read(5,*) (magmom(i,id), id=1,3)
            if ( (magmom(i,1).eq.0.d0).and.(magmom(i,2).eq.0.d0).and.
     &         (magmom(i,3).eq.0.d0) ) then
               write(6,110) i
               stop
            end if
         end do
         do i = 1,natom
            read(5,*) alphadegmom(i),betadegmom(i),gammadegmom(i)
         end do
      end if
      if (option.eq.1) then
         read(5,*) npoints
         if (npoints.le.0) then
            write(6,120)
            stop
         end if
         do j = 1,npoints
            read(5,*) (r(i,j),i=1,3)
            do i = 1,natom
               if ((qi(1,i).eq.r(1,j)).and.
     &            (qi(2,i).eq.r(2,j)).and.
     &            (qi(3,i).eq.r(3,j))) then
                  write(6,130) j
                  stop
               end if
               do m = 1,3
                  if ((br(1,m)+qi(1,i).eq.r(1,j)).and.
     &               (br(2,m)+qi(2,i).eq.r(2,j)).and.
     &               (br(3,m)+qi(3,i).eq.r(3,j))) then
                     write(6,130) j
                     stop
                  end if
               end do
            end do
         end do
      end if
c
c     End of reading input data
c
      if (option.eq.1) then
         do j = 1,npoints
            do i = 1,3
               qi(i,natom+j) = r(i,j)
            end do
            magmod(natom+j) = 1.d0
         end do
      else
         npoints = 0
      end if
c
c     Calculation of the matrix elements smat(i,j,lm) for l = 2.
c     To calculate the magnetostatic dipolar energy it is necessary to calculate
c     the matrix elements smat(i,j,k+4), which correspond to l = 2 according
c     to the following correspondence:
c     k = 1, 2, 3, 4, 5 --> magnetic quantum number m = -2, -1, 0, 1, 2
c     The indexes i and j refer to the basis atoms
c
      call madelung(smat,2,2,alat,natom,rmax,gmax,br,qi,bound,option,
     &              npoints)
c
c     Symmetry properties of the matrix elements smat(i,j,lm) for l = 2
c
      do i = 1,natom+npoints
         do j = i,natom+npoints
            do m = 5,9
               smat(j,i,m) = smat(i,j,m)
            end do
         end do
      end do

      fact(0) = 1.d0
      do i=1,100
         fact(i) = fact(i-1)*i
      end do

      if (option.ne.1) write(6,140)
c
c     Calculation of the magnetostatic dipolar energy per unit cell of a lattice 
c     of parallel magnetic dipoles for ndir directions of the magnetic dipoles 
c     and the magnetostatic dipolar anisotropy energy, i.e. the difference between the 
c     magnetostatic dipolar energies of these ndir directions (first - second), 
c     ..., (first - ndir)
c
      if (option.eq.0) then
         do id = 1,ndir
c
c     Calculation of the matrix elements rot(i,j) for l = 2 and for every 
c     direction id. To calculate the magnetostatic dipolar energy it is necessary 
c     to calculate only the matrix elements rot(i,3), which are equal to 
c     D2m0(alphadeg,betadeg,gammadeg) according to the following correspondence:
c     i = 1, 2, 3, 4, 5 --> magnetic quantum number m = -2, -1, 0, 1, 2
c
            call rotmat(3,3,0,alphadeg(id),betadeg(id),gammadeg(id),rot,
     &                 fact,nkmmax)

            write(6,'(/)')
            do m = 1,5
               write(6,150) m,rot(m,3)
            end do
c
c     We have taken into account that the madelung subroutine deals
c     with REAL spherical harmonics and that the rotmat subroutine
c     deals with COMPLEX spherical harmonics
c
            sumtotal=0.d0
            do i = 1,natom
               do j = 1,natom
                  sum(i,j)= -sqrt(2.d0)*smat(i,j,5)*imag(rot(5,3))
     &                      +sqrt(2.d0)*smat(i,j,6)*imag(rot(4,3))
     &                                 +smat(i,j,7)*rot(3,3)
     &                      -sqrt(2.d0)*smat(i,j,8)*real(rot(4,3)) 
     &                      +sqrt(2.d0)*smat(i,j,9)*real(rot(5,3))
                  sumtotal = sumtotal + magmod(i)*magmod(j)*sum(i,j)
               end do
            end do
c
c     Attention: in the following lines magdipenergy is in Rydbergs
c
            magdipenergy(id) = c*sumtotal 

            write(6,160) id,alphadeg(id),betadeg(id),gammadeg(id),
     &                   magdipenergy(id)/2d0,magdipenergy(id),
     &                   magdipenergy(id)*ry
         end do
c
c     End of do loop over directions
c
c
c     Attention: in the following lines magdipanienergy is in Rydbergs
c
         do id = 2,ndir
            magdipanienergy=magdipenergy(1)-magdipenergy(id)
            write(6,170) 1,id,magdipanienergy/2d0,magdipanienergy,
     &                   magdipanienergy*ry
         end do
      end if
c
c     Calculation of the magnetostatic dipolar field per unit cell at npoints 
c     specific points r, due to a lattice of parallel or non-parallel magnetic 
c     dipoles. These points must not coincide with the positions of the basis 
c     atoms in the unit cell. The next type of calculation deals with these points
c
      if (option.eq.1) then
c
c     In the unit cell there are, in general, several basis atoms. The 
c     magnetic moments of the atom i of the basis form a lattice of magnetic 
c     moments. The magnetic field due to this lattice has only a component 
c     and it is parallel to the direction of the magnetic moments. The other 
c     two components perpendicular to the mentioned one are null. The program
c     calculates the component parallel to the magnetic moments and later 
c     multiplies this quantity by the unitary vector of the corresponding 
c     magnetic moments in cartesian coordinates
c
         do j = 1,npoints
            do i = 1,natom
               call rotmat(3,3,0,alphadegmom(i),betadegmom(i),
     &                    gammadegmom(i),rot,fact,nkmmax)
               sum(i,natom+j) = 
     &                      -sqrt(2.d0)*smat(i,natom+j,5)*imag(rot(5,3))
     &                      +sqrt(2.d0)*smat(i,natom+j,6)*imag(rot(4,3))
     &                                 +smat(i,natom+j,7)*rot(3,3)
     &                      -sqrt(2.d0)*smat(i,natom+j,8)*real(rot(4,3)) 
     &                      +sqrt(2.d0)*smat(i,natom+j,9)*real(rot(5,3))
c
c     Attention: in the following lines magfield and magfieldtotal are in 
c     atomic Rydbergs units. muB = e*hbar/(2m) = sqrt(2.d0) in atomic 
c     Rydbergs units
c
               do m = 1,3
                  magfield(i,j,m) = 
     &            -magmod(i)*magmom(i,m)*c*sum(i,natom+j)/sqrt(2.d0)
                  magfieldtotal(j,m) = magfieldtotal(j,m) 
     &                                 + magfield(i,j,m)
               end do
            end do   ! do i = 1,natom
c
c     Writing magnetic dipolar fields in atomic Rydbergs units
c
            write(6,180) j,(r(i,j),i=1,3)
            do i = 1,natom
               write(6,190) i,(magfield(i,j,m),m=1,3)
            end do
            write(6,200) (magfieldtotal(j,m),m=1,3)
c
c     Writing magnetic dipolar fields in Teslas
c
            write(6,210) j,(r(i,j),i=1,3)
            do i = 1,natom
               write(6,220) i,(tesla*magfield(i,j,m),m=1,3)
            end do
            write(6,230) (tesla*magfieldtotal(j,m),m=1,3)
         end do
      end if
c
c     Calculation of the magnetostatic dipolar energy per unit cell of a lattice 
c     of parallel or non-parallel magnetic dipoles and the corresponding magnetic
c     dipolar field at the positions of the basis atoms in the unit cell
c
      if (option.eq.2) then
         do j = 1,natom
            do i = 1,natom
               call rotmat(3,3,0,alphadegmom(i),betadegmom(i),
     &                    gammadegmom(i),rot,fact,nkmmax)
               sum(i,j) = -sqrt(2.d0)*smat(i,j,5)*imag(rot(5,3))
     &                    +sqrt(2.d0)*smat(i,j,6)*imag(rot(4,3))
     &                               +smat(i,j,7)*rot(3,3)
     &                    -sqrt(2.d0)*smat(i,j,8)*real(rot(4,3)) 
     &                    +sqrt(2.d0)*smat(i,j,9)*real(rot(5,3))
c
c     Attention: in the following lines magfield and magfieldtotal are in 
c     atomic Rydbergs units. muB = e*hbar/(2m) = sqrt(2.d0) in atomic 
c     Rydbergs units
c
               do m = 1,3
                  magfield(i,j,m) = 
     &            -magmod(i)*magmom(i,m)*c*sum(i,j)/sqrt(2.d0)
                  magfieldtotal(j,m) = magfieldtotal(j,m) 
     &                                 + magfield(i,j,m)
               end do
            end do   ! do i = 1,natom
         end do   ! do j = 1,natom

         magdipenergy(1) = 0.d0
         do i = 1,natom
            do m = 1,3
               magdipenergy(1) = magdipenergy(1) 
     &         - magfieldtotal(i,m)*magmom(i,m)*magmod(i)*sqrt(2.d0)
            end do
         end do
c
c     Writing magnetic dipolar fields in atomic Rydbergs units
c
         do j = 1,natom
            write(6,180) j,(qi(i,j)/alat,i=1,3)
            do i = 1,natom
               write(6,190) i,(magfield(i,j,m),m=1,3)
            end do
            write(6,200) (magfieldtotal(j,m),m=1,3)
c
c     Writing magnetic dipolar fields in Teslas
c
            write(6,210) j,(qi(i,j)/alat,i=1,3)
            do i = 1,natom
               write(6,220) i,(tesla*magfield(i,j,m),m=1,3)
            end do
            write(6,230) (tesla*magfieldtotal(j,m),m=1,3)
         end do
c
c     Writing magnetostatic dipolar energy for the present magnetic configuration
c
         write(6,240) magdipenergy(1)/2d0,magdipenergy(1),
     &                magdipenergy(1)*ry
      end if

 100  format('Option is not equal to 0, 1 or 2',/,'Program stops')
 110  format('Magnetic moment of atom',i2,' is equal to zero',/,
     & 'Program stops')
 120  format('Number of points must be bigger or equal to 1',/,
     & 'Program stops')
 130  format('Point',i2,'is equal to the position of one basis atom',
     & ' in the unit cell',/,'Program stops')
 140  format(//,'The isotropic Fermi contact contribution is not',
     & ' included in',/,'the magnetostatic dipolar energies and',
     & ' magnetic dipolar fields')
 150  format('rot(',i1,') = ',2f12.8)
 160  format(/,'Magnetostatic dipolar energy per unit cell for direction
     & ',i1,': ',/,'Euler angles alpha, beta and gamma in degrees: ',
     & 3f9.5,/,f12.8,' hartrees',/,f12.8,' rydbergs',/,f12.8,' eV')
 170  format(//,'Magnetostatic dipolar anisotropy energy per unit cell'
     & ,/,'Direction',i2,' - Direction',i2,': ',
     & /,f12.8,' hartrees',/,f12.8,' rydbergs',/,f12.8,' eV')
 180  format(//,'Magnetic dipolar fields in atomic Rydbergs units',
     & ' at the point',i2,/,'in lattice parameter units:',3f7.3,/)
 190  format('Due to the atoms of type ',i2,': ',3f12.8)
 200  format(/,'Total magnetic dipolar field:',3f12.8)
 210  format(//,'Magnetic dipolar fields in Teslas at the point',i2,/,
     & 'in lattice parameter units:',3f7.3,/)
 220  format('Due to the atoms of type ',i2,': ',3f12.5)
 230  format(/,'Total magnetic dipolar field:',3f12.5)
 240  format(/,'Magnetostatic dipolar energy per unit cell for the '
     & 'present magnetic configuration: ',/,f12.8,' hartrees',/,f12.8,
     & ' rydbergs',/,f12.8,' eV')
      end

      SUBROUTINE ROTMAT( NK1, NK2,IREL, ALFDEG,BETDEG,GAMDEG,
     &                   ROT, FACT, NKMMAX )
C   ********************************************************************
C   *                                                                  *
C   *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
C   *           ( ALFDEG, BETDEG, GAMDEG )                             *
C   *                                                                  *
C   *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
C   *            EQS. (4.8), (4.12) AND (4.13)                         *
C   *                                                                  *
C   *   12/11/96  HE  deal with beta = 0                               *
C   ********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMPLEX*16 CI, CZ
      PARAMETER ( CI = (0.0D0,1.0D0), CZ = (0.0D0,0.0D0) )
      REAL*8 PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )   
C
      REAL*8     NUM, MSB05, MSB05SQ, MSB05PW, J,M1,M2, RFAC,X
      REAL*8     FACT(0:100)

      INTEGER    S, SLOW, SHIGH, OFF, NK1, NK2
      COMPLEX*16 EMIM2A, EMIM1G, ROT(NKMMAX,NKMMAX) 
C                       
C INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
      RFAC(X) = FACT( NINT(X) )
C
      DO 20 I2=1,NKMMAX 
      DO 20 I1=1,NKMMAX 
20    ROT(I1,I2) = CZ
C
       CB05   =   DCOS( BETDEG*0.5D0*PI/180.0D0 )
       CB05SQ =   CB05 *  CB05
      MSB05   = - DSIN( BETDEG*0.5D0*PI/180.0D0 )
      MSB05SQ =  MSB05 * MSB05
C     
      OFF = 0
      DO 100 K=NK1,NK2 
      IF( IREL .LT. 2 ) THEN
         L = K - 1
         J = L
      ELSE
         L = K/2
         IF( L*2 .EQ. K ) THEN
            J = L - 0.5D0
         ELSE 
            J = L + 0.5D0
         END IF         
      END IF

      NMUE = NINT( 2*J + 1 )
C
         DO 90 IM2 = 1, NMUE
         M2 = - J + (IM2-1.0D0)
         EMIM2A = CDEXP( -CI*M2*ALFDEG*PI/180.0D0 )
C
            DO 80 IM1 = 1, NMUE
            M1 = - J + (IM1-1.0D0)
            EMIM1G = CDEXP( -CI*M1*GAMDEG*PI/180.0D0 )
C
            IF( ABS(BETDEG) .LT. 1D-8 ) THEN
               IF( IM1 .EQ. IM2 ) THEN 
                  SUM = 1.0D0
               ELSE
                  SUM = 0.0D0
               END IF
            ELSE
               SLOW   = MAX(          0, NINT(M1-M2) )
               SHIGH  = MIN( NINT(J-M2), NINT( J+M1) )
                CB05PW =  CB05**NINT(2*J+M1-M2-2*SLOW    +2)
               MSB05PW = MSB05**NINT(    M2-M1+2*SLOW    -2)
               DOM = (-1.0D0)**(SLOW-1) *
     &          DSQRT( RFAC(J+M1)*RFAC(J-M1)*RFAC(J+M2)*RFAC(J-M2) )
               SUM = 0.0D0
C
               DO S=SLOW,SHIGH
                  DOM = -DOM
                  NUM =    FACT(S) * RFAC(J-M2-S)
     &                             * RFAC(J+M1-S) * RFAC(M2-M1+S)
                   CB05PW =  CB05PW /  CB05SQ
                  MSB05PW = MSB05PW * MSB05SQ
                  SUM = SUM + (DOM/NUM) * CB05PW * MSB05PW
               END DO
            END IF
C           
         ROT(OFF+IM2,OFF+IM1) = EMIM1G * SUM * EMIM2A
80         continue
C
90       CONTINUE

      OFF = OFF + NMUE
100   CONTINUE
C
      RETURN
      END

      SUBROUTINE MADELUNG(SMAT,linit,lmax,alat,natom,rmax,gmax,br,qi,
     +                    bound,option,npoints)
      PARAMETER (NMAXD=22000,NATOMD=80,ISHLD=3000,LMAXD=2)
      parameter (lmmaxd = (lmaxd+1)**2)
      real gn(4,nmaxd),qi(3,natomd),rm(4,nmaxd),ylm(lmmaxd),
     +     smat(natomd,natomd, (2*lmaxd+1)**2),
     +     amat(natomd,natomd,lmmaxd,lmmaxd),
     +     bmat(natomd,natomd,lmmaxd),br(3,3)
      integer nsg(ishld),nsr(ishld),linit,lmax,natom,option,npoints
      real alat,rmax,gmax,bound

cb ic
c      write(6,*) ' lmax : ',lmax
c      lmmax = (lmax + 1 )**2
ce ic
      call gaunt2
      call latvec(alat,natom+npoints,ngmax,nrmax,nshlg,nshlr,nsg,nsr,1,
     +            gn,rm,qi,vol,rmax,gmax,br)
      call strmat(alat,linit,lmax,natom,ngmax,nrmax,nsg,nsr,
     +            nshlg,nshlr,gn,rm,qi,smat,vol,bound,option,npoints)
      end
*DECK erfcex
      real function erfcex(z)
c-----------------------------------------------------------------------
c
c     calculates complementary errorfunction times sqrt(pi)
c      times exp(z*z)  by continued fractions
c
c-----------------------------------------------------------------------
c     .. scalar arguments ..
      real z
c     ..
c     .. local scalars ..
      real bound,erf1,exzz,f,fa,q,ratio,sqrtpi,term,u,ua,v,x,xa,y,z2,zz
c     ..
c     .. intrinsic functions ..
      intrinsic abs,atan,exp,sqrt
c     ..
c     .. data statements ..
      data bound/3.e-11/
c     ..
      sqrtpi = sqrt(4.0*atan(1.0))
      zz = z*z
      exzz = exp(zz)
c
c---> choose algorithm
c
      if (z.lt.1.5) then
 
         z2 = 2.0*zz
         erf1 = z
         ratio = 1.0
         term = z
   10    continue
         ratio = ratio + 2.0
         term = term*z2/ratio
         erf1 = erf1 + term
         if (term.gt.bound) go to 10
         erfcex = sqrtpi*exzz - 2.0*erf1
 
      else
c
c---> continued fraction expansion : abramowitz p. 298, eq. (7.1.14)
c
         u = 1.0
         v = 0.0
         x = z
         y = 1.0
         q = 0.5
         f = (u+v*q)/ (x+y*q)
   20    continue
         ua = u
         u = u*z + v*q
         v = ua
         xa = x
         x = x*z + y*q
         y = xa
         q = q + 0.5
         fa = f
         f = (u+v*q)/ (x+y*q)
         if (abs(fa-f).gt.bound*f) go to 20
         erfcex = f
      end if
 
      end
*DECK gamfc
      subroutine gamfc(alpha,glh,lmax,r)
c----------------------------------------------------------------------
c
c      calculation of convergence function
c
c       glh = i(alpha,l)/r**(l+1)*sqrt(pi)
c
c      with
c            alpha = r times the splitting paramter lamda
c      and
c            i(x,l) = erfc(x) + exp(-x*x)/sqrt(pi) *
c
c                                sum ( 2**i * x**(2i-1) / (2i-1)!! )
c                              1..i..l
c
c-----------------------------------------------------------------------
c     .. scalar arguments ..
      real alpha,r
      integer lmax
c     ..
c     .. array arguments ..
      real glh(0:lmax)
c     ..
c     .. local scalars ..
      real arg,facl,fex
      integer l
c     ..
c     .. external functions ..
      real erfcex
      external erfcex
c     ..
c     .. intrinsic functions ..
      intrinsic exp,real
c     ..
      arg = alpha*alpha
      glh(0) = erfcex(alpha)
      facl = 2.0*alpha
c
c---> recursion
c
      do 10 l = 1,lmax
         glh(l) = glh(l-1) + facl
         facl = facl*arg/ (real(l)+0.5)
   10 continue
      fex = 1.0/exp(arg)
      do 20 l = 0,lmax
         fex = fex/r
         glh(l) = glh(l)*fex
   20 continue
 
      end
*DECK gaunt2
      subroutine gaunt2
c-----------------------------------------------------------------------
c     sets up values needed for gaunt1
c        m. weinert  january 1982
c
c     changed for calculating with real spherical harmonics
c                                           b.drittler  july 1987
c-----------------------------------------------------------------------
c     .. parameters ..
      integer natomd,lmaxd,lpotd
*CALL I1
      PARAMETER (NATOMD=80,LMAXD=2,LPOTD=2)
      integer n,lassld
      parameter (n=4*lmaxd,lassld=n)
c     ..
c     .. arrays in common ..
      real w(n),x(n),yr(n,0:lassld,0:lassld)
c     ..
c     .. local scalars ..
      real a,cd,cth,fac,fpi,rf,sth,t
      integer k,l,lomax,m,nn
c     ..
c     .. local arrays ..
      real p(0:lassld+1,0:lassld)
c     ..
c     .. external subroutines ..
      external grule
c     ..
c     .. intrinsic functions ..
      intrinsic atan,sqrt
c     ..
c     .. common blocks ..
      common /assleg/w,x,yr
c     ..
c     .. save statement ..
      save
c     ..
      fpi = 16.*atan(1.0)
      rf = fpi** (1.0/3.0)
      lomax = lassld
c
c--->    obtain gauss-legendre points and weights
c
      nn = 2*n
      call grule(nn,x,w)
c
c--->    generate associated legendre functions for m.ge.0
c
      do 10 k = 1,n
         cth = x(k)
         sth = sqrt(1.0-cth*cth)
         fac = 1.0
c
c--->    loop over m values
c
         do 20 m = 0,lomax
            fac = - (2*m-1)*fac
            p(m,m) = fac
            p(m+1,m) = (2*m+1)*cth*fac
c
c--->    recurse upward in l
c
            do 30 l = m + 2,lomax
               p(l,m) = ((2*l-1)*cth*p(l-1,m)- (l+m-1)*p(l-2,m))/ (l-m)
   30       continue
            fac = fac*sth
   20    continue
c
c--->    multiply in the normalization factors
c
         do 40 l = 0,lomax
            a = rf*sqrt((2*l+1)/fpi)
            cd = 1
            yr(k,l,0) = a*p(l,0)
            do 50 m = 1,l
               t = (l+1-m)* (l+m)
               cd = cd/t
               yr(k,l,m) = a*sqrt(2.0*cd)*p(l,m)
   50       continue
   40    continue
   10 continue
      end
*DECK grule
      subroutine grule(n,x,w)
c
c***********************************************************************
c
c     determines the (n+1)/2 nonnegative points x(i) and
c     the corresponding weights w(i) of the n-point
c     gauss-legendre integration rule, normalized to the
c     interval [-1,1]. the x(i) appear in descending order.
c
c     this routine is from 'methods of numerical integration',
c     p.j. davis and p. rabinowitz, page 369.
c
c***********************************************************************
c
 
c     .. scalar arguments ..
      integer n
c     ..
c     .. array arguments ..
      real w(n),x(n)
c     ..
c     .. local scalars ..
      real d1,d2pn,d3pn,d4pn,den,dp,dpn,e1,fx,h,p,pi,pk,pkm1,pkp1,t,t1,
     +     u,v,x0
      integer i,it,k,m
c     ..
c     .. intrinsic functions ..
      intrinsic atan,cos
c     ..
      pi = 4.*atan(1.)
      m = (n+1)/2
      e1 = n* (n+1)
      do 10 i = 1,m
         t = (4*i-1)*pi/ (4*n+2)
         x0 = (1.- (1.-1./n)/ (8.*n*n))*cos(t)
c
c--->    iterate on the value  (m.w. jan. 1982)
c
         do 20 it = 1,2
            pkm1 = 1.
            pk = x0
            do 30 k = 2,n
               t1 = x0*pk
               pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1
               pkm1 = pk
               pk = pkp1
   30       continue
            den = 1. - x0*x0
            d1 = n* (pkm1-x0*pk)
            dpn = d1/den
            d2pn = (2.*x0*dpn-e1*pk)/den
            d3pn = (4.*x0*d2pn+ (2.-e1)*dpn)/den
            d4pn = (6.*x0*d3pn+ (6.-e1)*d2pn)/den
            u = pk/dpn
            v = d2pn/dpn
            h = -u* (1.+.5*u* (v+u* (v*v-u*d3pn/ (3.*dpn))))
            p = pk + h* (dpn+.5*h* (d2pn+h/3.* (d3pn+.25*h*d4pn)))
            dp = dpn + h* (d2pn+.5*h* (d3pn+h*d4pn/3.))
            h = h - p/dp
            x0 = x0 + h
   20    continue
         x(i) = x0
         fx = d1 - h*e1* (pk+.5*h* (dpn+h/3.* (d2pn+.25*h* (d3pn+
     +        .2*h*d4pn))))
         w(i) = 2.* (1.-x(i)*x(i))/ (fx*fx)
   10 continue
      if (m+m.gt.n) x(m) = 0.
      end
*DECK latvec
      subroutine latvec(alat,natom,ngmax,nrmax,nshlg,nshlr,nsg,nsr,
     +                  iprint,gn,rm,qi,vol,rmax,gmax,br)
c-----------------------------------------------------------------------
c
c     generate lattice vectors of direct and reciprocal space from
c      basic translation vectors br
c
c-----------------------------------------------------------------------
 
c     .. parameters ..
      integer nmaxd,natomd,ishld
      PARAMETER (NMAXD=22000,NATOMD=80,ISHLD=3000)
c     ..
c     .. scalar arguments ..
      real alat,vol
      integer iprint,natom,ngmax,nrmax,nshlg,nshlr
c     ..
c     .. array arguments ..
      real gn(4,nmaxd),qi(3,natomd),rm(4,nmaxd)
      integer nsg(ishld),nsr(ishld)
c     ..
c     .. local scalars ..
      real a,absgm,absrm,ag,ar,b,c,da,db,gmax,gx,gy,gz,pi,rmax,rx,ry,rz,
     +     vmin
      integer i,i1,i2,k,l,m,n,n1,ng,nr,nsh,nshl,numg,numgh,numr,numrh
c     ..
c     .. local arrays ..
      real absg(3),absr(3),bg(3,3),br(3,3),cj(4,nmaxd)
c     ..
c     .. intrinsic functions ..
      intrinsic abs,atan,ifix,max,mod,real,sqrt
c     ..
      pi = 4.0*atan(1.0)
c
c---> read lattice constant , no. of atoms per unit cell and cutoffs
c
!      read (5,*) alat,natom,rmax,gmax
!c     read (5,fmt=9000) alat,natom,rmax,gmax
      write (6,fmt=9010) alat,rmax,gmax
      rmax = rmax*alat
      gmax = gmax/alat
c
c---> basic trans. vectors and basis vectors
c
      write (6,fmt=9020)
      do 10 i = 1,3
!         read (5,*) br(1,i),br(2,i),br(3,i)
!c        read (5,fmt=9040) br(1,i),br(2,i),br(3,i)
         write (6,fmt=9040) br(1,i),br(2,i),br(3,i)
         br(1,i) = br(1,i)*alat
         br(2,i) = br(2,i)*alat
         br(3,i) = br(3,i)*alat
   10 continue
c
      write (6,fmt=9030) natom
      if (natom.gt.natomd) then
         stop ' 0 - latvec '
 
      else
c
         do 20 i = 1,natom
!            read (5,*) qi(1,i),qi(2,i),qi(3,i)
!c           read (5,fmt=9040) qi(1,i),qi(2,i),qi(3,i)
            write (6,fmt=9040) qi(1,i),qi(2,i),qi(3,i)
            qi(1,i) = qi(1,i)*alat
            qi(2,i) = qi(2,i)*alat
            qi(3,i) = qi(3,i)*alat
   20    continue
c
c---> generate primitive vectors bg of reciprocal space
c
         do 30 i = 1,3
            i1 = 1 + mod(i,3)
            i2 = 1 + mod(i1,3)
c
c---> cross product
c
            bg(1,i) = br(2,i1)*br(3,i2) - br(2,i2)*br(3,i1)
            bg(2,i) = br(3,i1)*br(1,i2) - br(3,i2)*br(1,i1)
            bg(3,i) = br(1,i1)*br(2,i2) - br(1,i2)*br(2,i1)
   30    continue
c
         vol = abs(br(1,1)*bg(1,1)+br(2,1)*bg(2,1)+br(3,1)*bg(3,1))
c
         write (6,fmt=9050)
         do 40 i = 1,3
            bg(1,i) = bg(1,i)/vol*2.0*pi
            bg(2,i) = bg(2,i)/vol*2.0*pi
            bg(3,i) = bg(3,i)/vol*2.0*pi
            write (6,fmt=9040) bg(1,i),bg(2,i),bg(3,i)
   40    continue
c
c---> estimate no. of lattice vectors
c
         do 50 i = 1,3
            absr(i) = sqrt(br(1,i)**2+br(2,i)**2+br(3,i)**2)
            absg(i) = sqrt(bg(1,i)**2+bg(2,i)**2+bg(3,i)**2)
   50    continue
         absrm = max(absr(1),absr(2),absr(3))
         absgm = max(absg(1),absg(2),absg(3))
         absrm = 2.0*pi/absrm
         absgm = 2.0*pi/absgm
         numr = 2* (ifix(rmax/absgm)+1) + 1
         numg = 2* (ifix(gmax/absrm)+1) + 1
         numrh = numr/2 + 1
         numgh = numg/2 + 1
         write (6,fmt=9110) numr,numg
c
c---> generate lattice vectors of real space
c
         if (iprint.gt.0) write (6,fmt=9060)
         nr = 0
         do 60 l = 1,numr
            a = real(l-numrh)
            do 70 m = 1,numr
               b = real(m-numrh)
               do 80 n = 1,numr
                  c = real(n-numrh)
                  rx = a*br(1,1) + b*br(1,2) + c*br(1,3)
                  ry = a*br(2,1) + b*br(2,2) + c*br(2,3)
                  rz = a*br(3,1) + b*br(3,2) + c*br(3,3)
                  ar = sqrt(rx*rx+ry*ry+rz*rz)
                  if (ar.le.rmax) then
                     nr = nr + 1
                     cj(1,nr) = rx
                     cj(2,nr) = ry
                     cj(3,nr) = rz
                     cj(4,nr) = ar
                  end if
 
   80          continue
 
   70       continue
 
   60    continue
c
         if (nr.gt.nmaxd) then
            write(6,*) 'nr = ',nr,' nmaxd = ',nmaxd
            stop ' 1 - latvec '
 
         else
 
            nrmax = nr
c
c---> sort vectors in order of increasing absolute value
c
            da = 1.e-06
            nsh = 0
            nshl = -1
            do 90 k = 1,nr
               vmin = rmax + 1.0
               do 100 n = 1,nr
                  if (cj(4,n)-vmin.lt.0) then
                     vmin = cj(4,n)
                     n1 = n
                  end if
 
  100          continue
 
 
               nshl = nshl + 1
               rm(1,k) = cj(1,n1)
               rm(2,k) = cj(2,n1)
               rm(3,k) = cj(3,n1)
               rm(4,k) = cj(4,n1)
               db = vmin
               if (db.gt.da+1.e-06) then
 
                  nsh = nsh + 1
                  if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
                  nsr(nsh) = nshl
                  if (iprint.gt.0) write (6,fmt=9070) k,rm(1,k),rm(2,k),
     +                rm(3,k),db
                  nshl = 0
                  da = db
 
               else if (iprint.gt.0) then
                  write (6,fmt=9070) k,rm(1,k),rm(2,k),rm(3,k),db
               end if
 
               cj(4,n1) = rmax + 1.0
   90       continue
            nsh = nsh + 1
            nshl = nshl + 1
            nsr(nsh) = nshl
            if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
            if (nsh.gt.ishld) then
               write(6,*) 'nsh = ',nsh,' ishld = ',ishld
               stop ' 2 - latvec '
 
            else
 
               nshlr = nsh
c
c---> generate lattice vectors of reciprocal space
c
               if (iprint.gt.0) write (6,fmt=9090)
               ng = 0
               do 110 l = 1,numg
                  a = real(l-numgh)
                  do 120 m = 1,numg
                     b = real(m-numgh)
                     do 130 n = 1,numg
                        c = real(n-numgh)
                        gx = a*bg(1,1) + b*bg(1,2) + c*bg(1,3)
                        gy = a*bg(2,1) + b*bg(2,2) + c*bg(2,3)
                        gz = a*bg(3,1) + b*bg(3,2) + c*bg(3,3)
                        ag = sqrt(gx*gx+gy*gy+gz*gz)
                        if (ag.le.gmax) then
                           ng = ng + 1
                           cj(1,ng) = gx
                           cj(2,ng) = gy
                           cj(3,ng) = gz
                           cj(4,ng) = ag
                        end if
 
  130                continue
 
  120             continue
 
  110          continue
c
               if (ng.gt.nmaxd) then
                  write(6,*) 'ng = ',ng,' nmaxd = ',nmaxd
                  stop ' 3 - latvec '
 
               else
 
                  ngmax = ng
c
c---> sort vectors in order of increasing abs. value
c
                  da = 1.e-06
                  nsh = 0
                  nshl = -1
                  do 140 k = 1,ng
                     vmin = gmax + 1.0
                     do 150 n = 1,ng
                        if (cj(4,n)-vmin.lt.0) then
                           vmin = cj(4,n)
                           n1 = n
                        end if
 
  150                continue
 
 
                     nshl = nshl + 1
                     gn(1,k) = cj(1,n1)
                     gn(2,k) = cj(2,n1)
                     gn(3,k) = cj(3,n1)
                     gn(4,k) = cj(4,n1)
                     db = vmin
                     if (db.gt.da+1.e-07) then
 
                        nsh = nsh + 1
                        if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
                        if (iprint.gt.0) write (6,fmt=9070) k,gn(1,k),
     +                      gn(2,k),gn(3,k),db
                        nsg(nsh) = nshl
                        nshl = 0
                        da = db
 
                     else if (iprint.gt.0) then
                        write (6,fmt=9070) k,gn(1,k),gn(2,k),gn(3,k),db
                     end if
 
                     cj(4,n1) = gmax + 1.0
  140             continue
                  nsh = nsh + 1
                  nshl = nshl + 1
                  nsg(nsh) = nshl
                  if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
                  if (nsh.gt.ishld) then
                     write(6,*) 'nsh = ',nsh,' ishld = ',ishld
                     stop ' 4 - latvec '
 
                  else
 
                     nshlg = nsh
 
                     write (6,fmt=9100) nrmax,ngmax,vol
 
 
 
                  end if
 
               end if
 
            end if
 
         end if
 
      end if
 
 
 9000 format (f10.5,i5,2f10.5)
 9010 format (/,1x,' lattice constant : ',f10.5,/,7x,' rmax : ',f10.5,
     +       ' gmax : ',f10.5)
 9020 format (/,1x,' primitive vectors of direct lattice ',/)
 9030 format (/,1x,i3,' basis vectors per unit cell ',/)
 9040 format (3f10.5)
 9050 format (/,1x,' primitive vectors of reciprocal lattice ',/)
 9060 format (/,1x,' real space lattice vectors ',/)
 9070 format (5x,i5,4f10.5)
 9080 format (1x,' shell no. ',i3,' contains ',i3,' points ',/)
 9090 format (/,1x,' reciprocal space lattice vectors ',/)
 9100 format (/,1x,' no. of lattice vectors   : ',i4,/,1x,
     +       ' no. of rec. lat. vectors : ',i4,/,1x,
     +       ' volume of the unit cell  : ',f10.5)
 9110 format (/,5x,' numr ',i3,' numg ',i3)
 
      end
*DECK strmat
      subroutine strmat(alat,linit,lmax,natom,ngmax,nrmax,nsg,nsr,
     +                  nshlg,nshlr,gn,rm,qi,smat,vol,bound,option,
     +                  npoints)
c-----------------------------------------------------------------------
c
c     calculation of lattice sums for l.ge.linit and l.le.lmax :
c
c                      ylm( q(i) - q(j) - rm )
c           sum      ===========================
c                    | q(i) - q(j) - rm |**(l+1)
c
c            - summed over all lattice vectors rm  -
c
c     ylm       : real spherical harmonic to given l,m
c     q(i),q(j) : basis vectors of the unit cell
c
c     in the case of i = j is rm = 0 omitted .
c
c     the ewald method is used to perform the lattice summations
c     the splitting parameter lamda is set equal sqrt(pi)/alat
c     (alat is the lattice constant) .
c
c     if the contribution of the last shell of the direct and the
c     reciprocal lattice is greater than bound a message is written
c
c                                     b.drittler may 1989
c
c-----------------------------------------------------------------------
c     .. parameters ..
*CALL I1
      integer natomd,lmaxd,lpotd
      PARAMETER (NATOMD=80,LMAXD=2,LPOTD=2)
      integer l2maxd,l2mmxd
      parameter (l2maxd=2*lpotd,l2mmxd= (l2maxd+1)**2)
c     ..
c     .. scalar arguments ..
      real alat,vol
      integer lpot,natom,ngmax,nrmax,nshlg,nshlr
c     ..
c     .. array arguments ..
      real gn(4,ngmax),qi(3,natomd),rm(4,nrmax),
     +     smat(natomd,natomd,l2mmxd)
      integer nsg(nshlg),nsr(nshlr)
c     ..
c     .. local scalars ..
      complex bfac,cfac,ci
      real alpha,beta,bound,dq1,dq2,dq3,dqdotg,expbsq,fpi,g1,g2,g3,ga,
     +     lamda,pi,r,r1,r2,r3,rfac,s,sgm
      integer i,i1,i2,it,l,lm,lmax,lmmax,m,nge,ngs,nre,nrs,nstart
cb ic
      integer linit,lminit,option,npoints,init2
ce ic
c     ..
c     .. local arrays ..
      complex stest(l2mmxd)
      real g(0:l2maxd),ylm(l2mmxd)
c     ..
c     .. external subroutines ..
      external gamfc,ymy
c     ..
c     .. intrinsic functions ..
      intrinsic abs,aimag,atan,exp,real,sqrt
c     ..
c     .. save statement ..
cb ic
c      save ci,bound
      save ci
ce ic
c     ..
c     .. data statements ..
cb ic
c      data ci/ (0.0,1.0)/,bound/1.0e-8/
      data ci/ (0.0,1.0)/
ce ic
c     ..
cb ic
c      lmax = 2*lpot
      if (linit.eq.lmax) then
         lminit = 1+(lmax*lmax)
      else
         lminit = (linit+1)*(linit+1)
      end if
ce ic
      lmmax = (lmax+1)* (lmax+1)
      pi = 4.0*atan(1.0)
      fpi = 4.0*pi
c
c---> choose proper splitting parameter
c
      lamda = sqrt(pi)/alat
c
c---> loop over atoms per unit cell
c
cb ic
      write(6,'(/)')
ce ic
      do 10 i1 = 1,natom
cb ic
c        do 20 i2 = 1,natom
         if ((option.eq.0).or.(option.eq.2)) init2 = i1
         if (option.eq.1) init2 = natom + 1
         do 20 i2 = i1,natom+npoints
ce ic
c
            dq1 = qi(1,i1) - qi(1,i2)
            dq2 = qi(2,i1) - qi(2,i2)
            dq3 = qi(3,i1) - qi(3,i2)
c
            stest(1) = -sqrt(fpi)/vol/ (4.0*lamda*lamda)
            do 30 lm = 2,lmmax
               stest(lm) = 0.0
   30       continue
c
c---> exclude the origine and add correction if i1.eq.i2
c
            if (i1.eq.i2) then
               stest(1) = stest(1) - lamda/pi
               nstart = 2
 
            else
               nstart = 1
 
            end if
c
c---> loop first over n-1 shells of real and recipro. lattice - then
c      add the contribution of the last shells to see convergence
c
            do 40 it = 1,2
               if (it.eq.1) then
                  nrs = nstart
                  ngs = 2
                  nre = nrmax - nsr(nshlr)
                  nge = ngmax - nsg(nshlg)
 
               else
                  nrs = nre + 1
                  ngs = nge + 1
                  nre = nrmax
                  nge = ngmax
               end if
c
c---> sum over real lattice
c
               do 50 i = nrs,nre
                  r1 = dq1 - rm(1,i)
                  r2 = dq2 - rm(2,i)
                  r3 = dq3 - rm(3,i)
                  call ymy(r1,r2,r3,r,ylm,lmax)
c
                  alpha = lamda*r
c
                  call gamfc(alpha,g,lmax,r)
c
                  do 60 l = 0,lmax
c
                     rfac = g(l)/sqrt(pi)
                     do 70 m = -l,l
                        lm = l* (l+1) + m + 1
                        stest(lm) = stest(lm) + ylm(lm)*rfac
   70                continue
   60             continue
 
   50          continue
c
c---> sum over reciprocal lattice
c
               do 80 i = ngs,nge
                  g1 = gn(1,i)
                  g2 = gn(2,i)
                  g3 = gn(3,i)
                  call ymy(g1,g2,g3,ga,ylm,lmax)
c
                  beta = ga/lamda
                  expbsq = exp(beta*beta/4.0)
                  dqdotg = dq1*g1 + dq2*g2 + dq3*g3
c
                  bfac = fpi*exp(ci*dqdotg)/ (ga*ga*expbsq*vol)
                  do 90 l = 0,lmax
c
                     do 100 m = -l,l
                        lm = l* (l+1) + m + 1
                        stest(lm) = stest(lm) + ylm(lm)*bfac
  100                continue
                     bfac = bfac*ga/ci/real(2*l+1)
   90             continue
   80          continue
c
               if (it.eq.1) then
cb ic
c                  do 110 lm = 1,lmmax
                  do 110 lm = lminit,lmmax
ce ic
                     if (abs(aimag(stest(lm))).gt.bound) then
                        go to 120
 
                     else
                        smat(i1,i2,lm) = real(stest(lm))
                        stest(lm) = 0.0
                     end if
 
  110             continue
 
               else
c
c---> test convergence
c
cb ic
c                  do 130 lm = 1,lmmax
                  do 130 lm = lminit,lmmax
ce ic
                     s = real(stest(lm))
                     smat(i1,i2,lm) = smat(i1,i2,lm) + s
cb ic
c                     if (abs(s).gt.bound) write (6,fmt=9000) i1,i2,lm,
c     +                   abs(s)
                      if (abs(s).gt.bound) write (6,fmt=9000) i1,i2,lm,
     +                   abs(s),bound
ce ic
 130             continue
               end if
 
   40       continue
cb ic
c            do 140 lm = 1,lmmax
            do 140 lm = lminit,lmmax
ce ic
               if (abs(smat(i1,i2,lm)).gt.bound) write (6,200) i1,i2,lm,
     +            smat(i1,i2,lm)
  140       continue
   20    continue
   10 continue
      return
 
  120 stop ' imaginary contribution to real lattice sum '
 
cb ic
 200  format ('smat(',i3,',',i3,',',i3,') =',f12.8)
c 9000 format (1x,' convergence of smat(',i2,i2,i4,') : ',e10.5,
c     +       ' is less than 1.0e-8 - use more lattice vectors ')
 9000 format (1x,' convergence of smat(',i2,i2,i4,') : ',e10.5,
     +       ' is more than ',e7.1,' - use more lattice vectors ')
ce ic
 
      end
*DECK ymy
      subroutine ymy(v1,v2,v3,r,ylm,lmax)
c-----------------------------------------------------------------------
c    this subroutine calculates real spherical harmonics with the
c     normalization : <y|y> =1
c    returns also r = length of vector v
c
c     generate the complex spherical harmonics for the vector v
c     using a stable upward recursion in l.  (see notes
c     by m. weinert.)
c                                  m.weinert  1982
c
c     converted to real spherical harmonics .
c                                  b.drittler 1987
c-----------------------------------------------------------------------
c
c     .. parameters ..
      integer lmaxd
      PARAMETER (LMAXD=4)
c     ..
c     .. scalar arguments ..
      real r,v1,v2,v3
      integer lmax
c     ..
c     .. array arguments ..
      real ylm((lmax+1)**2)
c     ..
c     .. local scalars ..
      real a,cd,cph,cth,fac,fpi,pi,rtwo,sgm,snull,sph,sth,t,xy,xyz
      integer i,l,m
c     ..
c     .. local arrays ..
      real c(0:lmaxd),p(0:lmaxd,0:lmaxd),s(0:lmaxd)
c     ..
c     .. intrinsic functions ..
      intrinsic sqrt
c     ..
c     .. save statement ..
      save
c     ..
c     .. data statements ..
      data snull/1.0e-20/
c     ..
      pi = 4.0*atan(1.0)
      fpi = 4.e0*pi
      rtwo = sqrt(2.0)
c
      if (lmax.gt.lmaxd) then
         stop 'ylm3'
 
      else
 
c
c--->    calculate sin and cos of theta and phi
c
         xy = v1**2 + v2**2
         xyz = xy + v3**2
c
         r = sqrt(xyz)
         if (xyz.le.0.0) then
            stop 'ylm=0'
 
         else
 
            if (xy.gt.snull*xyz) then
               xy = sqrt(xy)
               xyz = sqrt(xyz)
               cth = v3/xyz
               sth = xy/xyz
               cph = v1/xy
               sph = v2/xy
 
            else
 
               sth = 0.0
               cth = 1.0
               if (v3.lt.0) cth = -1.0
               cph = 1.0
               sph = 0.0
            end if
c
c--->    generate associated legendre functions for m.ge.0
c        loop over m values
c
            fac = 1.0
            do 10 m = 0,lmax - 1
               fac = - (2*m-1)*fac
               p(m,m) = fac
               p(m+1,m) = (2*m+1)*cth*fac
c
c--->    recurse upward in l
c
               do 20 l = m + 2,lmax
                  p(l,m) = ((2*l-1)*cth*p(l-1,m)- (l+m-1)*p(l-2,m))/
     +                     (l-m)
   20          continue
               fac = fac*sth
   10       continue
            p(lmax,lmax) = - (2*lmax-1)*fac
c
c--->    determine sin and cos of phi
c
            s(0) = 0.0
            s(1) = sph
            c(0) = 1.0
            c(1) = cph
            do 30 m = 2,lmax
               s(m) = 2*cph*s(m-1) - s(m-2)
               c(m) = 2*cph*c(m-1) - c(m-2)
   30       continue
c
c--->    multiply in the normalization factors
c
            i = 0
            do 40 l = 0,lmax
               i = i + l + 1
               a = sqrt((2*l+1)/fpi)
               cd = 1
               ylm(i) = a*p(l,0)
               sgm = -rtwo
               do 50 m = 1,l
                  t = (l+1-m)* (l+m)
                  cd = cd/t
                  t = a*sqrt(cd)
                  ylm(i+m) = sgm*t*p(l,m)*c(m)
                  ylm(i-m) = sgm*t*p(l,m)*s(m)
                  sgm = -sgm
   50          continue
               i = i + l
   40       continue
 
         end if
 
      end if
 
      end

      subroutine rinit(n,a)
c
c     Initialization of the first n values of the real array a with zero
c
      integer n
      real    a(n)

      do 10 i = 1,n
         a(n) = 0.d0
 10   continue
      end

      subroutine cinit(n,a)
c
c     Initialization of the first n values of the complex array a with zero
c
      integer n
      complex a(n)

      do 10 i = 1,n
         a(n) = (0.d0,0.d0)
 10   continue
      end





