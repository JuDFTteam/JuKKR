!-------------------------------------------------------------------------------
! SUBROUTINE: CLSGEN99
!> @brief This subroutine is used to create the clusters around each atom
!> where repulsive potentials will be positioned.
!> @details Calculate the cluster of each atom by the lattice parameters avaliable.
!> Sort the atoms in a unique way : big r, big z, big y
!> compare the positions with the previous clusters to see if there is
!> a difference. If not keep only previous clusters and make indexing if
!> a new cluster is found then check dimensions and continue for the new
!> atom.
!> @note Small bug in assigning clusters removed 29/04/2003 v.popescu
!> IATCLUS(NCLSD) is pointing to the first atomic site associated with a tb-cluster
!> @note JC: This routine seems to be deprecated and not used in the actual code
!-------------------------------------------------------------------------------
subroutine CLSGEN99(NAEZ,RR,NR,RBASIS,KAOEZ,Z,CLS,NACLS,REFPOT,ATOM,EZOA,  &
      NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,TLEFT,TRIGHT,RCLS,    &
      RCUT,RCUTXY,L2DIM,ALAT,NATYP,NEMB,NPRINCD,NACLSD,NCLSD)

   use mod_version_info

   implicit none
   ! .. Arguments ..
   integer, intent(in) :: NR        !< Number of real space vectors rr
   integer, intent(in) :: NEMB      !< Number of 'embedding' positions
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NCLSD     !< Maximum number of different TB-clusters
   integer, intent(in) :: NLEFT     !< Number of repeated basis for left host to get converged electrostatic potentials
   integer, intent(in) :: NRIGHT    !< Number of repeated basis for right host to get converged electrostatic potentials
   integer, intent(in) :: NACLSD    !< Maximum number of atoms in a TB-cluster
   integer, intent(in) :: NPRINCD   !< Number of principle layers, set to a number >= NRPINC in output of main0
   integer, intent(in) :: NLBASIS   !< Number of basis layers of left host (repeated units)
   integer, intent(in) :: NRBASIS   !< Number of basis layers of right host (repeated units)
   double precision, intent(in) :: ALAT !< Lattice constant in a.u.
   double precision, intent(in) :: RCUT   !< Parameter for the screening cluster along the z-direction
   double precision, intent(in) :: RCUTXY !< Parameter for the screening cluster along the x-y plane
   integer, dimension(NAEZ+NEMB), intent(in) :: REFPOT !< Ref. pot. card  at position
   double precision, dimension(NATYP), intent(in) :: ZAT       !< Nuclear charge
   double precision, dimension(3,0:NR), intent(in) :: RR       !< Set of real space vectors (in a.u.)
   double precision, dimension(3,NAEZ+NEMB), intent(in) :: RBASIS   !< Position of atoms in the unit cell in units of bravais vectors
   double precision, dimension(3,NACLSD,NCLSD), intent(in) :: RCLS  !< Real space position of atom in cluster
   !
   integer, dimension(NAEZ+NEMB), intent(inout) :: CLS   !< Cluster around atomic sites
   integer, dimension(NCLSD), intent(inout)     :: NACLS !< Number of atoms in cluster
   integer, dimension(NACLSD,NAEZ+NEMB), intent(in)   :: ATOM  !< Atom at site in cluster
   integer, dimension(NACLSD,NAEZ+NEMB), intent(in)   :: EZOA  !< EZ of atom at site in cluster
   integer, dimension(NATYP,NAEZ+NEMB), intent(in)    :: KAOEZ !< Kind of atom at site in elem. cell
   ! .. Local variables
   !
   integer :: I,N1,INUM,ISUM,NA,NUMBER,N,NPRIN,ITEST1,ITEST
   integer :: POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,I0,I1,ICLUSTER
   integer, dimension(NACLSD) :: IATOM
   integer, dimension(NACLSD) :: IEZOA
   integer, dimension(NACLSD) :: ISORT
   integer, dimension(NCLSD)  :: IATCLS
   integer, dimension(NAEZ,NAEZ) :: ICOUPLMAT
   !
   double precision :: R,R2,EPSSHL
   double precision :: RCUT2,RCUTXY2,RXY2
   double precision, dimension(3)      :: TMP
   double precision, dimension(NACLSD) :: RSORT
   double precision, dimension(3)      :: ZPERLEFT
   double precision, dimension(3)      :: ZPERIGHT
   double precision, dimension(3,NACLSD) :: RG
   double precision, dimension(3,NACLSD) :: RCLS1
   double precision, dimension(3,NEMB+1) :: TLEFT
   double precision, dimension(3,NEMB+1) :: TRIGHT
   !
   logical :: L2DIM,CLUSTCOMP
   !
   logical :: TEST,LSPHER
   external :: DSORT,CLUSTCOMP
   intrinsic :: MIN,SQRT
   !
   data EPSSHL   / 1.0D-4 /
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This is generating the clusters which have a distance smaller
   ! than RCUT and RCUTXY in plane .
   ! The cluster atoms are ordered with radious and then z>y>x
   ! The ordering allows an easy comparison of clusters
   ! The principal layer for each layer (atom in unit cell) is
   ! calculated also for each cluster and the maximum number
   ! is returned. Some dimension tests are also done
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !---------------------------------------------------------------------------
   ! OUTPUT
   !---------------------------------------------------------------------------
   write (1337,'(79(1H=))')
   write (1337,'(16X,A)') 'CLSGEN99: generation of TB-clusters coordinates'
   write (1337,'(79(1H=))')
   write (1337,*)
   !---------------------------------------------------------------------------
   ! OUTPUT
   !---------------------------------------------------------------------------
   LSPHER = .FALSE.
   write(1337,*) 'RCUT = ',rcut,' RCUTXY = ',rcutxy
   if (ABS(rcutxy - rcut).LT.1.D-4) then
      write(1337,*) 'Spherical Clusters are created'
      LSPHER = .TRUE.
   end if
   !----------------------------------------------------------------------------
   if (TEST('clusters')) then
      open(8,FILE='clusters',STATUS='UNKNOWN')
      call version_print_header(8)
      write(8,9005) NAEZ
      write(8,9030) ALAT
      write(8,9010) (Z(KAOEZ(1,I)),I=1,NAEZ)
      write(8,9020) (KAOEZ(1,I),I=1,NAEZ)
   end if
   !----------------------------------------------------------------------------
   ICLUSTER = 1
   do N = 1,NCLSD
      IATCLS(N) = 0
      NACLS(N) = 0
   end do
   call RINIT(3*NACLSD*NCLSD,RCLS)
   !
   RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
   RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)

   do 2 JATOM = 1,NAEZ       ! loop in all atoms or layers
      CLS(JATOM) = 0
      NUMBER = 0             ! counter for atoms in cluster
      do NA = 1,NAEZ  ! loop in all atoms
         do N=0,NR    ! loop in all bravais vectors
            do I=1,3
               TMP(I) = RR(I,N)+RBASIS(I,NA)-RBASIS(I,JATOM)
            end do
            RXY2 =  TMP(1)**2+TMP(2)**2
            R2   =  TMP(3)**2
            if (LSPHER) R2 = R2 + RXY2

            if ( (RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2) )  then
               NUMBER = NUMBER + 1
               if (NUMBER.GT.NACLSD) then
                  write (6,*)' ERROR: Dimension NACLSD in inc.cls too small', &
                  NUMBER, NACLSD
                  stop '   < CLSGEN99 >'
               end if
               !
               ATOM(NUMBER,JATOM) = NA ! store the atom in elem cell
               EZOA(NUMBER,JATOM) = N ! store the bravais vector
               do I=1,3
                  RCLS1(I,NUMBER) = TMP(I)
               end do
            end if
         end do              ! N loop in bravais
      end do                 ! NA loop in NAEZ

      !-------------------------------------------------------------------------
      ! In the case of 2 dimensional case loop in the atoms
      ! outside.
      !-------------------------------------------------------------------------
      IF (L2DIM) THEN
         !----------------------------------------------------------------------
         ! Somehow meshy
         ! ATOM gives the kind of atom
         !----------------------------------------------------------------------
         do N=0,NR
            do I=NLEFT,1,-1  ! loop in some layers on left side
               do I1=NLBASIS,1,-1 ! loop in representative atoms on left side
                  do I0=1,3
                     TMP(I0) = RR(I0,N) + TLEFT(I0,i1) + (I-1)*ZPERLEFT(I0)   &
                              - RBASIS(I0,JATOM)
                  end do
                  RXY2 =  TMP(1)**2+TMP(2)**2
                  R2   =  TMP(3)**2
                  if (LSPHER) R2 = R2 + RXY2

                  if ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) then

                     NUMBER = NUMBER + 1
                     if (NUMBER.GT.NACLSD) then
                        write (6,*)'ERROR: Dimension NACLSD in inc.cls too small',  &
                           NUMBER, NACLSD
                        stop '   < CLSGEN99 >'
                     end if
                     !
                     ATOM(NUMBER,JATOM) = -NAEZ-I1 ! negative values are used in dlke1.f
                     EZOA(NUMBER,JATOM) = N ! I,I1 are negative
                     do I0=1,3
                        RCLS1(I0,NUMBER) = TMP(I0)
                     end do
                  end if
               end do
            end do
            !
            !
            do I=1,NRIGHT
               do I1=1,NRBASIS
                  do I0=1,3
                     TMP(I0) = RR(I0,N)+ TRIGHT(I0,i1) + (I-1)*ZPERIGHT(I0)
                     &                             - RBASIS(I0,JATOM)
                  end do
                  RXY2 =  TMP(1)**2+TMP(2)**2
                  R2  =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
                  if ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) then
                     NUMBER = NUMBER + 1
                     if (NUMBER.GT.NACLSD) then
                        write (6,*)' ERROR: Dimension NACLSD in inc.cls too small', &
                           NUMBER, NACLSD
                        stop '   < CLSGEN99 >'
                     end if
                     ATOM(NUMBER,JATOM) = -NAEZ-NLBASIS-I1
                     EZOA(NUMBER,JATOM) = N
                     do I0=1,3
                        RCLS1(I0,NUMBER) = TMP(I0)
                     end do
                  end if
               end do
            end do
         end do                    ! loop in all bravais lattices
      end if                 ! l2dim interface calculation
      !-------------------------------------------------------------------------
      ! Now the atom JATOM Has it's cluster first
      ! sort the atoms of the cluster in increasing order. First by distance
      ! Then by z then by y
      !-------------------------------------------------------------------------
      do ia=1,number
         rsort(ia) = SQRT(RCLS1(1,ia)**2+RCLS1(2,ia)**2+RCLS1(3,ia)**2)
         rsort(ia) = 100000000.d0*rsort(ia)+ &
                     10000.d0*RCLS1(3,ia)+   &
                     10.d0*RCLS1(2,ia)+      &
                     0.1d0*RCLS1(1,ia)
      end do
      !
      call DSORT(RSORT,ISORT,NUMBER,POS)
      ! Rearange exchange ia with ib
      ! MAP temporarily to another array
      do IA=1,NUMBER
         do I=1,3
            RG(I,IA)    = RCLS1(I,IA)
         end do
         IATOM(IA) = ATOM(IA,JATOM)
         IEZOA(IA) = EZOA(IA,JATOM)
      end do
      ! Now use correct order
      do IA =1,NUMBER
         IB = ISORT(IA)
         do I=1,3
            RCLS1(I,IA) = RG(I,IB)
         end do
         ATOM(IA,JATOM) = IATOM(IB)
         EZOA(IA,JATOM) = IEZOA(IB)
      end do
      !-------------------------------------------------------------------------
      ! Now the clusters have a unique sorting and can be compared with
      ! each other Check if ICLUSTER was found previously
      !-------------------------------------------------------------------------
      do ICU = 1,ICLUSTER-1

         N1 = NACLS(ICU)
         ! return true if found before
         if (CLUSTCOMP(RCLS,REFPOT,ATOM,IATCLS(ICU),ICU,N1,RCLS1,NUMBER,JATOM,NACLSD))then
             CLS(JATOM) = ICU
         endif
      end do
      if (CLS(JATOM).EQ.0) then
         if (ICLUSTER.GT.NCLSD) then
            write(6,*) 'Please, increase the parameter NCLSD in', &
            ' inc.cls to a value greater equal ',ICLUSTER,' .'
            stop 'Dimension error.'
         end if
         CLS(JATOM) = ICLUSTER
         NACLS(ICLUSTER) = NUMBER
         IATCLS(ICLUSTER) = JATOM
         do IN = 1,NUMBER
            do II=1,3
               RCLS(II,IN,ICLUSTER) = RCLS1(II,IN)
            end do
            write(1337,800) jatom,atom(in,jatom),ezoa(in,jatom),  &
            (rcls1(i,in),i=1,3),sqrt(rcls1(1,in)**2+rcls1(2,in)**2+rcls1(3,in)**2)
         end do
         ICLUSTER = ICLUSTER + 1
      end if
      !-------------------------------------------------------------------------
      write(1337,*) 'Atom ',JATOM,' has cluster ', CLS(JATOM),'with ',NUMBER,' sites'
   2 continue                      ! JATOM = 1,NAEZ
   !----------------------------------------------------------------------------
   ! Now all clusters of all atoms are found print out
   ! and test the results...
   !----------------------------------------------------------------------------
   do 22 JATOM = 1,NAEZ
      !-------------------------------------------------------------------------
      IC = CLS(JATOM)
      NUMBER = NACLS(IC)
      !-------------------------------------------------------------------------
      if (TEST('clusters')) then
         write(8,FMT=1030) NUMBER
         write(8,FMT=1030) JATOM,IC
         do I=1,NUMBER
            R = SQRT(CLS(1,I,IC)**2+RCLS(2,I,IC)**2+RCLS(3,I,IC)**2)
            write(8,1041) (RCLS(II,I,IC),II=1,3),ATOM(I,JATOM),Z(ABS(ATOM(I,JATOM))),R
         end do
      end if
      !-------------------------------------------------------------------------
      !  Now print out the coupling matrix
      !-------------------------------------------------------------------------
      do IAT = 1,NAEZ
         ICOUPLMAT(JATOM,IAT) = 0
         do I=1,NUMBER
            if (ATOM(I,JATOM).EQ.IAT) then
               ICOUPLMAT(JATOM,IAT) = 1
            end if
         end do
      end do
      write(1337,9060) JATOM,(ICOUPLMAT(JATOM,IAT),IAT=1,NAEZ)

   22 end do ! Do loop in JATOM (second test loop)
   !----------------------------------------------------------------------------
   ! Now testing the shape of the dyson equation
   !----------------------------------------------------------------------------
   ITEST1 = ICOUPLMAT(NAEZ,1)+ICOUPLMAT(1,NAEZ)
   IF (ITEST1.NE.0) THEN
      WRITE (1337,*) ' This is not a banded matrix '
   END IF
   nprin = 0
   do i=1,naez-1
      itest = icouplmat(1,i)
      itest1 = icouplmat(1,i+1)
      if (itest.eq.1.and.itest1.eq.0) then
         nprin = i-1
      end if
   end do
   if (nprin.eq.0) nprin = 1
   write(1337,*) '***********************************************'
   write(1337,*) '********** TESTING THE COUPLING MATRIX ********'
   write(1337,*) '***********************************************'
   write(1337,9090) NPRIN
   if (NPRIN.NE.NPRINCD) then
      write(1337,*) 'Please change NPRINCD in your inc.p file'
      write(1337,*) 'from ',NPRINCD, 'to ',NPRIN
      write(1337,*) ' ******** RESULTS COULD BE WRONG ******** '
   end if
   ! Now check if you can divide the matrix correctly
   if (MOD(NAEZ,NPRIN).NE.0) then
      write(1337,*) ' Your matrix cannot be divided in '
      write(1337,*) ' Principal layers. Use a number of layers '
      write(1337,*) ' which is multiple of ',NPRIN
   end if
   !> @note \f$ NL + 2*NPR*NL - 4 \sum_{n=1^{NPR}}  {n} \f$
   isum = 0
   do i=1,nprin
      isum = isum + i
   end do
   inum = NAEZ + 2*NPRIN*NAEZ - 2*ISUM
   !
   !
   isum = 0
   do i=1,naez
      do i0=1,naez
         isum = isum + icouplmat(i,i0)
      end do
   end do

   if (ISUM.EQ.INUM) then
      write(1337,*) ' Your matrix is BAND DIAGONAL'
   else
      write(1337,*) ' Your matrix is *NOT* BAND DIAGONAL',ISUM,INUM
   end if
   write(1337,*)
   write(1337,*) ' Sub clsgen99  exiting <<<<<<<<<<<<<'
   !----------------------------------------------------------------------------
   800    format(3I5,4F8.4)
   1030   format(3I8)
   ! 1041   formaT(3F15.8,I6,F7.1,F12.6)
   1041   format(3E27.19,I5,F5.1,E17.9)
   9005   format(I4)
   9010   format(('# Z     ',20F4.0))
   9020   format(('# KAOEZ ',20I4))
   9030   format(F12.7,6x,'ALAT')
   9060   format(I4,1X,200I1)
   9090   format('The Number of layers in each Principal Layer = ',I5)
   return
end subroutine CLSGEN99
