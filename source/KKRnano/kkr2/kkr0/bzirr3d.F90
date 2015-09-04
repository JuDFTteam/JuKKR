      subroutine bzirr3d(nkp, nkxyz, kpoibz, kp, recbv, bravais, wtkp, volbz, rsymat, nsymat, isymindex, irr, iprint)
!===========================================================================
!info 
!info   find irreducible BZ and create mesh in it.
!info           original version Arthur Ernst, Daresbury, 27/03/98
!                                fcc lattice bug corrected
!                                modified 26/5/99, 2/06/99
! Modified on 20.01.2000 To use the 2d inversion!
! Modified Dec 2001 - Apr 2002 to deal with relativistic case, bugs
!          removed (full BZ integration gives now same results as 
!          symetrised calculation)   HE/VP, Munich
!
!===========================================================================
!  Input:
!         nkxyz : original k-mesh net in the 3 directions of the reciprocal
!                 lattice vectors (not xyz directions).
!         recbv : reciprocal unit cell 
!                 (normalized: recbv(i,j)=G(i,j)*alat/(2*pi), i=x,y,z j=1,2,3)
!                 if G is the "normal" rec. lattice vector.
!      bravais  : direct unit cell in a.u.
!      rsymat   : symmetry rotation matrices of real lattice
!      nsymat   : number of rotation matrices of real lattice
!         irr   : if .true. then irreducible BZ else take all BZ
!
!  Output:
!         nkp   : number of points in irreducible BZ
!         kp    : k-point mesh
!         wtkp  : weights for k-points
!        volbz  : volume of the BZ
!  Inside:
!         lbk   : Flag showing if the mesh point jx,jy,jz has already been 
!                 taken. Could also be a logical variable.
!==========================================================================
      implicit none
      integer, intent(inout) :: nkxyz(3)
      integer, intent(in) :: kpoibz, nsymat, iprint
      integer, intent(in) :: isymindex(*)
      double precision, intent(in) :: recbv(3,3), bravais(3,3), rsymat(64,3,3) ! todo: transpose into (3,3,*)
      logical, intent(in) :: irr !< use irreducible BZ only
      integer, intent(out) :: nkp
      double precision, intent(out) :: kp(3,*), wtkp(*), volbz
      
      double precision, external :: ddet33, ddot
      
      double precision bginv(3,3), bgmat(3,3), bgp(3,3), bv(3), cf(3), u(3,3,48), gq(3,3), v1
      integer nbgp(3,3,48), mkxyz(3), ind1(3), ind2(3)
      integer i,j,jx,jy,jz,nsym,iws,iwt,is,nk,n,ndim,isym
      integer(kind=1), allocatable :: ibk(:,:,:)
      logical surface

! --> check if we are in surface mode
      surface = .false.
      if (bravais(1,3) == 0.d0 .and. bravais(2,3) == 0.d0 .and. bravais(3,3) == 0.d0) surface = .true. 

      mkxyz(1:3) = nkxyz(1:3)
      
      ndim = 3
      if (surface) then
        mkxyz(3) = 1
        nkxyz(3) = 1 ! this is the only operation why nkxyz is not intent(in)
        ndim = 2
      endif ! surface

      allocate(ibk(0:mkxyz(1)-1,0:mkxyz(2)-1,0:mkxyz(3)-1))
      
!-------------------
!
! create small unit cell for the integration in the reciprocal space (gq), steps along the basis vectors
      do j = 1, 3
        gq(1:3,j) = recbv(1:3,j)/dble(mkxyz(j))
      enddo ! j
      
      do j = 1,3
        do i = 1, 3
          bgmat(i,j) = ddot(3,gq(1,i),1,gq(1,j),1) ! ==? dot_product(gq(1:3,i), gq(1:3,j)) 
        enddo ! i
      enddo ! j

      call rinvgj(bginv, bgmat, 3, ndim)
      
      if (irr) then
        nsym = nsymat
        
        if (nsym > 48) stop 'DIM problem in bzirr3d.'
        if (surface .and. nsym > 12) stop 'DIM problem in bzirr3d, surf mode.'

        do n = 1, nsym
          isym = isymindex(n)
          u(1:3,1:3,n) = rsymat(isym,1:3,1:3)
        enddo ! n
      else  ! irreducible
        nsym = 1
        u(1:3,1:3,1) = dble(reshape([1,0,0,0,1,0,0,0,1], [3,3])) ! unit matrix
      endif ! irreducible

!-----------------------------------------------------------------------
!          rotate the 3 step vectors   gq  --  it should be possible
!          to express the rotated vectors  b(j) in terms of the old ones
!     b(i) = sum(j) gq(j) * n(j,i)  with integer coefficients n(j,i)
!
!==========================================================================
      do is = 1, nsym
        if (iprint > 2) write(*,fmt="(5x,'rotated GQ  for IROT=',i3)") is
    
        do i = 1, 3
          call dgemv('n',3,3,1.d0,u(1,1,is),3,gq(1,i),1,0.d0,bgp(1,i),1)
          do j = 1, 3
            bv(j) = ddot(3,gq(1,j),1,bgp(1,i),1) ! ==? dot_product(gq(1:3,j), bgp(1:3,i))
          enddo ! j
          call dgemv('n',3,3,1d0,bginv,3,bv,1,0.d0,cf,1)
          do j = 1, 3
            if (abs(nint(cf(j)) - cf(j)) > 1d-8) write (*,fmt="(5x,2i3,3f7.3)") i, j, cf(j)
            nbgp(j,i,is) = nint(cf(j))
          enddo !
          if (iprint > 2) write(*,fmt="(5x,i3,3f7.3,2x,3f7.3,2x,3f7.3,2x,3i3)") i, bgp(1:3,i), bv, cf, nbgp(1:3,i,is)
        enddo ! i
      enddo ! is
!========================================================================

      nk = product(mkxyz(1:3))

      ibk(:,:,:) = 0 ! since the legacy code has x as outermost of the three loops, we index ibk in (3,2,1)-order

      nkp = 0
      iws = 0

      do jx = 0, mkxyz(1)-1
        ind1(1) = jx
        do jy = 0, mkxyz(2)-1
          ind1(2) = jy
          do jz = 0, mkxyz(3)-1
            ind1(3) = jz
            
            if (ibk(ind1(3),ind1(2),ind1(1)) == 0) then
              ! a new irreducible k-point is found
              nkp = nkp + 1
              iwt = 0
              
!========================================================================
              do isym = 1, nsym ! loop over symmetries

!               rotate k-vector  nkp  and transform into parallelepiped
!          rot k = sum(i) m(i) rot gq(i)
!                = sum(i) m(i) sum(j) n(i,j) gq(j)
!                = sum(j) [sum(i) m(i) n(i,j)] gq(j)

                do j = 1, 3
                  is = 0
                  do i = 1, 3
                    is = is + ind1(i)*nbgp(j,i,isym)
                  enddo !
                  is = modulo(is, mkxyz(j))
!                 if (is < 0) is = is + mkxyz(j) ! this should never happen using modulo instead of mod in the line above
                  ind2(j) = is
                enddo ! j

                if (ibk(ind2(3),ind2(2),ind2(1)) == 0) then
                    ibk(ind2(3),ind2(2),ind2(1)) = isym ! mark this k-point to be covered by symmetry number isym
                    iwt = iwt + 1 ! increase the weight of this irreducible k-point
                endif
!               mesh point in the unit cell found.
              enddo ! isym
!========================================================================
              
              if (nkp > kpoibz) then
                write(6,*) ' No. of k-points ',nkp,' > KPOIBZ '
                stop ' <BZIRR3D> increase parameter KPOIBZ'
              endif
              
              do i = 1, 3
                kp(i,nkp) = 0.d0
                do j = 1, 3
                  kp(i,nkp) = kp(i,nkp) + gq(i,j)*dble(ind1(j))
                enddo ! j
              enddo ! i
              wtkp(nkp) = dble(iwt)/dble(nk)
              
            endif ! not lbk
            iws = iws + iwt ! add op to global weight sum
            
          enddo ! jz
        enddo ! jy
      enddo ! jx

      deallocate(ibk)
      
      if (surface) then
        v1 = abs(recbv(1,1)*recbv(2,2) - recbv(1,2)*recbv(2,1))
      else
        v1 = ddet33(recbv) 
      endif

      volbz = 0.d0
      do i = 1, nkp
        volbz = volbz + wtkp(i)*v1
        wtkp(i) = wtkp(i)*v1/dble(nsym)
      enddo ! i
      
      endsubroutine bzirr3d
