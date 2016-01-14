module BrillouinZone_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: bzkint0, readKpointsFile
  
  integer, parameter :: maxmshd = 8
  
  contains

  
  !------------------------------------------------------------------------
  ! Read k-mesh file
  subroutine readKpointsFile(maxmesh, nofks, bzkp, volcub, volbz)
    integer, intent(in) :: maxmesh
    integer, intent(out) :: nofks(:)
    double precision, intent(out) :: bzkp(:,:,:), volcub(:,:), volbz(:)

    integer, parameter :: fu = 52 ! file unit
    integer :: i, l, ios, maxm
    character(len=*), parameter :: defname='kpoints', altname='new.kpoints'

    open (fu, file=altname, form='formatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then ! default to reading the old kpoint file
      open (fu, file=defname, form='formatted', status='old', action='read', iostat=ios)
      if (ios /= 0) then ! failed to read the old kpoint file
        write(*,*) 'ERROR: file "',defname,'" not found, nor file "'-altname-'".'
        die_here('file "'-defname-'" not found!')
      endif
    else ! file new.kpoints exists - use those kpoints
      write(*,*) 'WARNING: rejecting file "',defname,'" - using file "',altname,'" instead.'
      warn(6, 'rejecting file "'-defname-'" - use file "'-altname-'" instead.')
    endif

    rewind (fu) ! rewind the files to the beginning (probably legacy code)
    
    read (fu, fmt=*) maxm ! read the maximum number of k-meshes given in that file
    do l = 1, min(maxmesh, maxm)
      read (fu, fmt='(i8,f15.10)') nofks(l), volbz(l)
      do i = 1, nofks(l)
        read (fu, fmt=*) bzkp(1:3,i,l), volcub(i,l)
      enddo ! i
    enddo ! l
    close (fu)
    
  endsubroutine ! readKpointsFile
  
  
  subroutine bzkint0(naez, rbasis, bravais, recbv, nsymat, isymindex, &
                     dsymll, intervxyz, ielast, ez, iemxd, kmesh, maxmesh, lmax, krel, ekmd, nowrite, kpms)
    use Symmetry_mod, only: pointgrp, findgroup, symtaumat
    use BrillouinZoneMesh_mod, only: BrillouinZoneMesh!, create, load, store, destroy

    integer, parameter :: nsymaxd=48

    integer, intent(out) :: ekmd
    integer, intent(in)  :: naez, krel, lmax, iemxd
    integer, intent(out) :: nsymat, maxmesh
    integer, intent(in)  :: intervxyz(3), ielast
    double complex, intent(out) :: dsymll((lmax+1)**2,(lmax+1)**2,nsymaxd) ! (lmmaxd,lmmaxd,nsymaxd)
    double complex, intent(in) :: ez(iemxd)
    double precision, intent(in) :: bravais(3,3), rbasis(3,naez), recbv(3,3)
    integer, intent(out) :: isymindex(nsymaxd)
    integer, intent(out) :: kmesh(iemxd)
    logical, intent(in) :: nowrite
    type(BrillouinZoneMesh), intent(out) :: kpms(:)

!!  logical, external :: test
#define test(STRING) .false.
    
    logical :: lirr
    double precision :: rsymat(3,3,64)
    character(len=10) :: rotname(64)
    integer :: lmmaxd, iprint

    lmmaxd = (lmax+1)**2

    write(6,'(79(1h=),/,15x,a)') 'BZKINT0: finding symmetry, setting BZ integration'
    write(6,'(79(1h=),/)')

    call pointgrp(rsymat, rotname) ! set up the default symmetries and their names
    call findgroup(bravais, recbv, rbasis, naez, rsymat, rotname, isymindex, nsymat, naez)

    lirr = .true.
    iprint = 0 ; if (test('TAUSTRUC')) iprint = 2

    ! test: full Brillouin zone integration
    if (test('fullBZ  ')) then
      nsymat = 1 ! limit the number of applied symmetries to the 1st one (which is always unity)
      lirr = .false.
      write(6,'(8x,2a,/)') 'Test option < fullBZ > : overriding NSYMAT,', ' generate full BZ k-mesh'
    endif ! full Brillouin zone

    ! generate Brillouin zone k-meshes
    call bzkmesh(intervxyz, maxmesh, lirr, bravais, recbv, nsymat, rsymat, isymindex, ielast, ez, kmesh, iprint, iemxd, ekmd, nowrite, kpms)

    call symtaumat(rotname, rsymat, dsymll, nsymat, isymindex, naez, lmmaxd, naez, lmax+1, krel, iprint, nsymaxd)

    ! now dsymll hold nsymat symmetrization matrices
  endsubroutine ! bzkint0

      
  subroutine bzkmesh(nbxyz, maxmesh, lirr, bravais, recbv, nsymat, rsymat, isymindex, ielast, ez, kmesh, iprint, iemxd, ekmd, nowrite, kpms)
    use BrillouinZoneMesh_mod, only: BrillouinZoneMesh, create, load, store, destroy
    integer, intent(in) :: iemxd, nbxyz(3), nsymat, iprint, ielast
    logical, intent(in) :: lirr, nowrite
    double precision, intent(in) :: bravais(3,3), recbv(3,3), rsymat(64,3,3)
    integer, intent(in) :: isymindex(:) !< dim(nsymaxd)
    double complex, intent(in) :: ez(iemxd) !< energy contour points

    integer, intent(out) :: maxmesh, ekmd
    integer, intent(out) :: kmesh(iemxd) !< mapping of k-meshes and energy contour points

    integer :: ie, ik, im, n, nb(3), ekmin, nxyz(3), nofks(maxmshd), newnofks(maxmshd), fu
    logical :: newkp, oldkp
    double precision :: wxyz(0:3), volBz, newvolBz
    double precision, allocatable :: kwxyz(:,:) !< dim(0:3,product(nbxyz(1:3)))
    type(BrillouinZoneMesh), intent(out) :: kpms(:)

    ! set number of different k-meshes 
    maxmesh = 1
    if (test('fix mesh')) then
      kmesh(1:ielast) = 1
      warn(6, "kmesh fixed to"+nbxyz)
    else
      do ie = 1, ielast
        n = int(1.001d0 + log(dimag(ez(ie))/dimag(ez(ielast)))/log(2.d0))
        kmesh(ie) = max(1, n)
        maxmesh = max(maxmesh, kmesh(ie))
!       write(*, '(a,f16.6,f12.6,9(a,i0))') 'renormalized energy:',ez(ie)/dimag(ez(ielast)),'  n=',n,'  kmesh(',ie,')=',kmesh(ie),'  max=',maxmesh
      enddo ! ie
      kmesh(1:2) = maxmesh ! why is this done?
      warn(6, "kmesh(1:2) is set to mesh #"-maxmesh)
    endif ! variable mesh

    if (maxmesh > maxmshd) die_here("Dimension ERROR: Please increase hard limit MAXMSHD to"+maxmesh)

    write(6,'(79(1h=))')
    write(6, fmt="(12x,' BZKMESH : creating k-mesh,',' write to file kpoints')")
    write(6,'(79(1h=))')
    write(6,*)
    write(6, fmt="(8x,'number of different k-meshes :',i2)") maxmesh
    write(6, fmt="(8x,'the direct lattice',i3,' symmetries will be used')") nsymat
    write(6, fmt="(/,8x,35(1h-),/,8x,'k-mesh NofKs N kx N ky N kz vol BZ',/,8x,35(1h-))")

    if (nowrite) then
      fu = 0
      inquire(file='kpoints', exist=oldkp)
      if (oldkp) then
        warn(6, "a file 'kpoints' exists and is not overwritten.")
      else
        warn(6, "file 'kpoints' is not created.")
      endif
    else
      fu = 54
      open(unit=fu, file='kpoints', action='write') ! create or overwrite existing file with the same name
      write(unit=fu, fmt='(9(a,i0))') ' ',maxmesh,'   =NumberOfMeshes'
    endif

    nb(1:3) = nbxyz(1:3)
    allocate(kwxyz(0:3,product(nbxyz(1:3))))
    do im = 1, maxmesh
      nb = max(1, nb)

      nxyz(1:3) = nb(1:3)
      
      call bzirr3d(nxyz, recbv, bravais, nofks(im), volBz, kwxyz, rsymat, nsymat, isymindex, lirr, iprint)
      
      write(6, fmt="(8x,2i6,3i5,f8.4)") im, nofks(im), nxyz(1:3), volBz
      
      call create(kpms(im), kwxyz(:,1:nofks(im)), volBz)

!     if (fu > 0) write(unit=fu, fmt='(i8,f15.10,/,(3f12.8,d20.10))') nofks(im), volBz, (kwxyz(1:3,ik), kwxyz(0,ik), ik=1,nofks(im))
      if (fu > 0) call store(kpms(im), fu)

      if (test('k-net   ')) write(6, fmt="(3f12.5,f15.8)") (kwxyz(1:3,ik), kwxyz(0,ik), ik=1,nofks(im)) ! output of k-mesh

      nb(1:3) = nb(1:3)/1.4 ! use less k-points on the next level grid
    enddo ! im ! loop over different meshes
    if (fu > 0) close(unit=fu, iostat=ie)
    deallocate(kwxyz, stat=ie)
    write(6, fmt="(8x,35(1h-),/)")

    ! check dimensions of ekmd precond. array ! fix: check regardless of iguessd

    write(6,'(79(1h=))')
    write(6,'(12x,a)') 'BZKMESH: checking dimensions of precond. arrays ...'

    ekmin = sum(nofks(kmesh(1:ielast)))
    ! set ekmd to the minimum value needed
    ekmd = ekmin
    write(6,*) '           EKMIN=',ekmin,'  EKMD=',ekmd

    ! check the format of the file called new.kpoints 
    inquire(file='new.kpoints', exist=newkp)
    if (newkp) then

      open(53, file='new.kpoints', form='formatted', action='read', status='old')
      rewind(53)
      do im = 1, maxmesh
        read(53, fmt='(i8,f15.10)') newnofks(im), newvolBz
        do ik = 1, newnofks(im)
          read(53, fmt=*) wxyz(1:3), wxyz(0) ! read dummies just to see if the format is good
        enddo ! ik
      enddo ! im
      close(53)

      ekmin = sum(newnofks(kmesh(1:ielast)))
      ! set ekmd to the minimum value needed
      ekmd = ekmin
      write(6,*) '      new: EKMIN=',ekmin,'  EKMD=',ekmd
      write(6,'(79(1h=))')
      write(6,*)
    endif ! newkp

  endsubroutine ! bzkmesh

  
  
  
  subroutine bzirr3d(nkxyz, recbv, bravais, nkp, volBz, kwxyz, rsymat, nsymat, isymindex, irr, iprint)
    use VectorMath_mod, only: ddet33
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
!         kwxyz : k-point mesh with weights in the 0-th component
!        volBz  : volume of the BZ
!  Inside:
!         lbk   : Flag showing if the mesh point jx,jy,jz has already been 
!                 taken. Could also be a logical variable.
!==========================================================================
    integer, intent(inout) :: nkxyz(3) !< would be intent(out) without the surface options
    integer, intent(in) :: nsymat, iprint
    integer, intent(in) :: isymindex(:) !< dim(nsymaxd)
    double precision, intent(in) :: recbv(3,3), bravais(3,3), rsymat(3,3,64) ! todo: transpose into (3,3,*)
    logical, intent(in) :: irr !< use irreducible BZ only
    integer, intent(out) :: nkp
    double precision, intent(out) :: volBz, kwxyz(0:,1:) ! (0:3,1:max...)
    
    external :: dgemv ! BLAS
    double precision :: bginv(3,3), bgmat(3,3), bgp(3,3), bv(3), cf(3), u(3,3,48), gq(3,3), v1
    integer :: nbgp(3,3,48), mkxyz(3), ind1(3), ind2(3)
    integer :: i, j, jx, jy, jz, nsym, iws, iwt, is, nk, n, ndim, isym
    integer(kind=1), allocatable :: ibk(:,:,:)
    logical :: surface

! --> check if we are in surface mode
    surface = all(bravais(1:3,3) == 0.d0)

    mkxyz(1:3) = max(1, nkxyz(1:3))
    
    ndim = 3
    if (surface) then
      mkxyz(3) = 1
      nkxyz(3) = 1 ! this is the only operation why nkxyz is not intent(in)
      ndim = 2
    endif ! surface

    
!-------------------
!
! create small unit cell for the integration in the reciprocal space (gq), steps along the basis vectors
    do j = 1, 3
      gq(1:3,j) = recbv(1:3,j)/dble(mkxyz(j))
    enddo ! j

    do j = 1, 3
      do i = 1, 3
        bgmat(i,j) = dot_product(gq(1:3,i), gq(1:3,j)) 
      enddo ! i
    enddo ! j

    call rinvgj(bginv, bgmat, 3, ndim)

    if (irr) then
      nsym = nsymat
      
      if (nsym > 48) stop 'DIM problem in bzirr3d.'
      if (surface .and. nsym > 12) stop 'DIM problem in bzirr3d, surf mode.'

      do n = 1, nsym
        isym = isymindex(n)
        u(1:3,1:3,n) = rsymat(1:3,1:3,isym)
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
      if (iprint > 2) write(*, fmt="(5x,'rotated GQ  for IROT=',i3)") is
  
      do i = 1, 3
        call dgemv('n',3,3,1.d0,u(1,1,is),3,gq(1,i),1,0.d0,bgp(1,i),1)
        do j = 1, 3
          bv(j) = dot_product(gq(1:3,j), bgp(1:3,i))
        enddo ! j
        call dgemv('n',3,3,1.d0,bginv,3,bv,1,0.d0,cf,1)
        do j = 1, 3
          if (abs(nint(cf(j)) - cf(j)) > 1.d-8) write(*, fmt="(5x,2i3,3f7.3)") i, j, cf(j)
          nbgp(j,i,is) = nint(cf(j))
        enddo !
        if (iprint > 2) write(*, fmt="(5x,i3,3f7.3,2x,3f7.3,2x,3f7.3,2x,3i3)") i, bgp(1:3,i), bv, cf, nbgp(1:3,i,is)
      enddo ! i
    enddo ! is
!========================================================================

    nk = product(mkxyz(1:3))

    allocate(ibk(0:mkxyz(3)-1,0:mkxyz(2)-1,0:mkxyz(1)-1))
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
                is = ind1(1)*nbgp(j,1,isym) + ind1(2)*nbgp(j,2,isym) + ind1(3)*nbgp(j,3,isym)
                ind2(j) = modulo(is, mkxyz(j))
              enddo ! j

              if (ibk(ind2(3),ind2(2),ind2(1)) == 0) then
                  ibk(ind2(3),ind2(2),ind2(1)) = isym ! mark this k-point to be covered by symmetry number isym
                  iwt = iwt + 1 ! increase the weight of this irreducible k-point
              endif
!             mesh point in the unit cell found.
            enddo ! isym
!========================================================================
            
            kwxyz(1:3,nkp) = gq(1:3,1)*ind1(1) + gq(1:3,2)*ind1(2) + gq(1:3,3)*ind1(3)
            kwxyz(0,nkp) = dble(iwt) ! un-normalized weight
            
          endif ! not lbk
          iws = iws + iwt ! add op to global weight sum
          
        enddo ! jz
      enddo ! jy
    enddo ! jx

    kwxyz(:,nkp+1:) = 0.d0
    
    deallocate(ibk, stat=i)
    
    if (surface) then
      v1 = abs(recbv(1,1)*recbv(2,2) - recbv(1,2)*recbv(2,1))
    else  ! surface
      v1 = ddet33(recbv) 
    endif ! surface

    volBz = 0.d0
    do i = 1, nkp
      volBz = volBz + kwxyz(0,i)
      kwxyz(0,i) = kwxyz(0,i)*v1/dble(nsym*nk) ! normalize
    enddo ! i
    volBz = volBz*v1/dble(nk)
 
  endsubroutine ! bzirr3d
  
  subroutine rinvgj(ainv,a,M,n)
!   ********************************************************************
!   *                      ainv = a**(-1)                              *
!   *  invert a using the Gauss-Jordan - algorithm                     *
!   *  the 1- matrix is not set up and use is made of its structure    *
!   *                    double precision version                      *
!   ********************************************************************
    integer, intent(in) :: M, n
    double precision, intent(inout) ::  a(M,M)
    double precision, intent(out) :: ainv(M,M)

    integer :: ic, l, ll
    double precision :: t, t1

    ainv(1,1) = 0.d0
    do ic = 1, n ! scan columns
      t1 = 1.d0/a(ic,ic) ! make a(ic,ic) = 1
      do l = ic+1, n
        a(ic,l) = a(ic,l)*t1
      enddo ! l

      do l = 1, ic-1
        ainv(ic,l) = ainv(ic,l)*t1
      enddo ! l
      ainv(ic,ic) = t1
                                
      do ll = 1, n ! make a(ll,ic) = 0 for ll /= ic
        if (ll /= ic) then
          t = a(ll,ic)
          do l = ic+1, n
            a(ll,l) = a(ll,l) - a(ic,l)*t
          enddo ! l
          do l = 1, ic-1
            ainv(ll,l) = ainv(ll,l) - ainv(ic,l)*t
          enddo ! l
          ainv(ll,ic) = -t1*t
        endif ! ll /= ic
      enddo ! ll
    enddo ! ic

  endsubroutine ! rinvgj
  
endmodule ! BrillouinZone_mod