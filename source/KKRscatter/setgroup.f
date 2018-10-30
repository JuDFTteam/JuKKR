      subroutine setgroup(bravais,recbv,rbasis,rfctor,nbasis,
     &                     rsymat,rotname,isymindex,nsymat)
c **********************************************************
c This subroutine set the rotation matrices It can be used
c for test purposes when the symmetries have to be set by hand...
c 
c input:  bravais(i,j)    true bravais lattice vectors
c                         i = x,y,z ; j = A, B, C (a.u.)
c         recbv(i,j)      reciprocal basis vectors 
c         rbasis          coordinates of basis atoms
c         nbasis          number of basis atoms
c         rfctor          alat/4/pi
c         rsymat          all 64 rotation matrices.
c         rotname         names for the rotation matrices
c output: nsymat          number of rotations that restore the lattice.
c         ISYMINDEX       index for the symmeties found
c         
c The array ISYMINDEX holds the numbers of the symmetry operations
c that are stored in array RSYMAT
c **********************************************************
      implicit none
      include 'inc.p'
      integer NSYMAXD,NMATD
      parameter (NSYMAXD=48)
      integer nbasis,nsymat
      integer isymindex(NSYMAXD)
      double precision BRAVAIS(3,3),RBASIS(3,NAEZD+NEMBD)
      double precision RSYMAT(64,3,3),recbv(3,3)
      double precision rfctor
c
c Local variables
c
      double precision r(3,4),rotrbas(3,naezd+nembd)
      double precision alat,bravais1(3,3)
      integer i,j,isym,nsym,i0,ia,ns
      Character*10 ROTNAME(64),symop(1)
c     Character*10 ROTNAME(64),symop(8)
      CHARACTER*10 CHAR(64)
      logical llatbas,latvec,LBULK



c      DATA SYMOP/'E         ','C2z       ','IC2a      ','IC2b      ',
c     &           'IC2x      ','IC2y      ','C4z       ','C4z-1     '/
c      DATA NSYM/8/

       DATA SYMOP/'E         '/
       DATA NSYM/1/

c
c for more or different symmery op's change the dimensions of symop
c
c     -------------------------------------------------------------
      alat = RFCTOR*8.D0*DATAN(1.0D0)
c     the ISYMINDEX array has the numbers of the symmetries found
c     
c     
      WRITE(6,*) ' Entered in sub SETGROUP '
      write(6,*) ' Using user defined symmetry operations '
      write(6,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'      
      NSYMAT = NSYM 
      DO NS=1,NSYM
            ISYM = 0
            do i0=1,64
               if (rotname(i0).eq.symop(ns)) THEN 
                  ISYM = I0
                  ISYMINDEX(NS) = I0
                  !write(6,*) rotname(i0)
               end if
            end do

            if (ISYM.eq.0) then
               WRITE (6,*) ' No symmetry found ERROR'
               STOP
            end if
c     
c     the symmetry is isym rotate
c      
      END DO ! ns =1,nsym

      write(6,*) 'Information from SetGroup'
      write(6,1020) NSYMAT
      write(6,*) 
      do i=1,nsymat
         I0 = ISYMINDEX(I)
         CHAR(I) =  ROTNAME(I0) 
      end do
c      write(6,1010) (CHAR(I),I=1,NSYMAT)
      write(6,*) '----------- * setgroup ends here * ---------------'
      write(6,*) 
 1010 FORMAT(5(A10,2X))
 1020 FORMAT(' Symmetries given for this lattice: ',I5)
c     
      END



























