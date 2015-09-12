module Symmetry_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: pointgrp, findgroup, symtaumat

  contains

  subroutine findgroup(bravais, recbv, rbasis, nbasis, rotmat, rotname, isymindex, nsym, naezd)
    use VectorMath_mod, only: ddet33
! **********************************************************
! This subroutine finds the rotation matrices that leave the
! real lattice unchanged. 
! input:  bravais(i,j)    true bravais lattice vectors
!                         i = x,y,z ; j = A, B, C (a.u.)
!         recbv(i,j)      reciprocal basis vectors 
!         rbasis          coordinates of basis atoms
!         nbasis          number of basis atoms
!         rotmat          all 64 rotation matrices.
!         rotname         names for the rotation matrices
! output: nsym            number of rotations that restore the lattice.
!         ISYMINDEX       index for the symmeties found
!
! This sub makes all 64 rotations in the basis vectors and bravais
! vectors and checks if the new rotated vectror belongs in the 
! lattice. The proper rotation must bring all vectors to a lattice
! vector. Information about the rotations found is printed in the end.         
! The array ISYMINDEX holds the numbers of the symmetry operations
! that are stored in array RSYMAT
! **********************************************************
    integer, parameter   :: nsymaxd=48
    integer, intent(out) :: nsym
    integer, intent(in)  :: nbasis, naezd
    integer, intent(out) :: isymindex(nsymaxd)
    double precision, intent(in) :: bravais(3,3), rbasis(3,naezd), recbv(3,3)
    double precision, intent(in) :: rotmat(3,3,64)
    character(len=*), intent(in) :: rotname(64)

!!  logical, external :: test
#define test(STRING) .false.
    
    double precision :: r(3,4), rotrbas(3,naezd), bravais1(3,3)
    integer :: i, j, isym, i0, ia
!   double precision :: mdotmp, mvecq(3,naezd), mvecqp(3,naezd)
    double precision :: mrotr(3,3), symdet, summdotmp
!   double precision :: stet
    character(len=10) :: name(64)
    logical :: llatbas, lbulk

    write(6, fmt="(5x,'< FINDGROUP > : Finding symmetry operations',/)")

    nsym = 0
    bravais1(1:3,1:3) = bravais(1:3,1:3)
    
!     check for surface mode. if so, set bravais1(3,3) very large, so
!     that only the in-plane symmetries are found. not checked, be careful of z--> -z!
    lbulk = .true.
!     now check the bravais vectors if they have a z component 
    if (all(bravais(1:3,3) == 0.d0)) lbulk = .false.
!     
    do isym = 1, 64
!
!--------------------------------- store rotation matrix
!
      mrotr(1:3,1:3) = rotmat(1:3,1:3,isym)
      summdotmp = 0.d0
      symdet = ddet33(mrotr)
!
!     rotate bravais lattice vectors
!     
!     in the case of slab/interface geometry look only for 
!     symmetry opperations that preserve the z axis..
!     
      if (lbulk .or. rotmat(3,3,isym) == 1.d0) then 
!     do rotation only in case bulk or if slab and z axis is restored..  
          
        do i = 1, 3            ! loop on bravais vectors
          do j = 1, 3         ! loop on coordinates
            r(j,i) = rotmat(j,1,isym)*bravais1(1,i) + rotmat(j,2,isym)*bravais1(2,i) + rotmat(j,3,isym)*bravais1(3,i)
          enddo ! j
        enddo ! i
!     
!     rotate the basis atoms p and take rsymat.p - p then
!     find if r = (rsymat.bravais + rsymat.p - p) belongs to the
!     lattice. this is done by function latvec by checking
!     if r.q = integer (q reciprocal lattice vector)
!     
        llatbas = .true.
        do ia = 1, nbasis      ! loop on basis atoms
          do j = 1, 3         ! loop on coordinates
            rotrbas(j,ia) = rotmat(j,1,isym)*rbasis(1,ia) + rotmat(j,2,isym)*rbasis(2,ia) + rotmat(j,3,isym)*rbasis(3,ia)
            rotrbas(j,ia) = rotrbas(j,ia) - rbasis(j,ia)
            r(j,4) = rotrbas(j,ia) 
          enddo ! j
!     do j=1,4
!     do i=1,3
!     write(6,*) 'rrr',i,j,r(i,j)
!     enddo
!     enddo
!     write(6,*) 'latvec',latvec(4,recbv,r)
          if (.not. latvec(4, recbv, r)) llatbas = .false.
          if (test('Oh-symm ') .and. isym <= 48)  llatbas = .true.
          if (test('Td-symm ') .and. isym <= 12)  llatbas = .true.
          if (test('Td-symm ') .and. isym >= 37 .and. isym <= 48) llatbas = .true.
        enddo               ! ia=1,nbasis
!     
!     if llatbas=.true. the rotation does not change the lattice 
!     
!     write(6,*) 'llatbas',llatbas
        if (llatbas) then
          nsym = nsym + 1
          isymindex(nsym) = isym
        endif ! llatbas
      endif ! (lbulk .or. (rotmat(3,3,isym) == 1))
    enddo ! isym = 1, nmatd
!     nsym symmetries were found
!     the isymindex array has the numbers of the symmetries found
 
    write(6, fmt="(8x,60(1h-))")
    if (lbulk) then
      write(6, fmt="(8X,'3D symmetries',$)") 
    else  ! bulk
      write(6, fmt="(8X,'surface symmetries',$)") 
    endif ! bulk
    write(6, fmt="(' found for this lattice: ',i2,/,8x,60(1h-))") nsym
    do i = 1, nsym
      i0 = isymindex(i)
      name(i) = rotname(i0) 
    enddo ! i
    do i = 1, nsym/5 + 1
      isym = min(5, nsym - (i-1)*5)
      write(6, fmt="(8x,5(a10,2x))") (name(j),j=(i-1)*5+1,(i-1)*5+isym)
    enddo ! i
    write(6, fmt="(8x,60(1h-),/)")
    
  endsubroutine findgroup
  
  
  
!*==symtaumat.f    processed by SPAG 6.05Rc at 15:50 on 10 Dec 2002
  subroutine symtaumat(rotname, rotmat, drot, nsym, isymindex, nqmax, nkmmax, nq, nl, krel, iprint, nsymaxd)
    use VectorMath_mod, only: ddet33
    use Constants_mod, only: pi
!   ********************************************************************
!   *                                                                  *
!   *  Find the symmetry matrices DROT that act on t, tau, ....        *
!   *  KREL=0: for real spherical harmonics                            *
!   *  KREL=1: for relativistic represntation                          *
!   *                                                                  *
!   *  The NSYM allowed symmetry operations are indicated by ISYMINDEX *
!   *  in the table rotmat. 
!   *                                                                  *
!   *  The routine determines first the Euler angles correponding      *
!   *  to a symmetry operation. Reflections are decomposed into        *
!   *  inversion + rotation for this reason.                           *
!   *                                                                  *
!   ********************************************************************
    integer, intent(in) :: iprint, krel, nkmmax, nl, nq, nqmax, nsym, nsymaxd
    double complex, intent(out) :: drot(nkmmax,nkmmax,48)
    integer, intent(in) :: isymindex(nsymaxd)
    double precision, intent(in) :: rotmat(3,3,64)
    character(len=*), intent(in) :: rotname(64)
    
    double complex, parameter :: ci=(0.d0,1.d0), c1=(1.d0,0.d0), c0=(0.d0,0.d0)
    double precision :: a, b, fact(0:100),rj, rmj, sk, symeulang(3,48), tet(1:3), co2, si2
    logical :: equal
    double complex :: dinv(nkmmax,nkmmax), dtim(nkmmax,nkmmax), rc(nkmmax,nkmmax), w1(nkmmax,nkmmax), w2(nkmmax,nkmmax)
    integer :: i,i1,i2,ind0q(nqmax),invflag(48),iq,irel,ireleff,isym, itop,j,k,l,loop,m,n,nk,nkeff,nkm,nlm,nok,ns
    double precision :: rmat(3,3), w

    equal(a,b) = (abs(a-b) < 1d-7) ! inline function?

    write(6, fmt="(5x,'<SYMTAUMAT> : rotation matrices acting on t/G/tau',//,8x,57(1h-),/,8x,'ISYM            INV          Euler angles      Unitarity',/,8x,57(1h-))")

    irel = krel*3
    nk = (1-krel)*nl + krel*(2*nl-1)
    nkm = (1+krel)*nl**2

    fact(0) = 1.d0
    do i = 1, 100
      fact(i) = fact(i-1)*dble(i)
    enddo ! i

    ind0q(1) = 0
    do iq = 2, nq
      ind0q(iq) = ind0q(iq-1) + nkm
    enddo ! iq

! ----------------------------------------------------------------------
!    rc  transforms from  real to  complex (l,m,s) - representation
!                 |lc> = sum[lr] |lr> * rc(lr,lc)
! ----------------------------------------------------------------------
    if (krel == 0) then
      nlm = nkm

      call zcopy(nkmmax*nkmmax,c0,0,rc,1)

      w = 1.d0/sqrt(2.d0)

      do l = 0, (nl-1)
        do m = -l, l
          i = l*(l+1) + m + 1
          j = l*(l+1) - m + 1

          if (m < 0) then
            rc(i,i) = -ci*w
            rc(j,i) = w
          endif
          if (m == 0) rc(i,i) = c1
          if (m > 0) then
            rc(i,i) = w*(-1.d0)**m
            rc(j,i) = ci*w*(-1.d0)**m
          endif
        enddo ! m
      enddo ! l
    endif
!
!=======================================================================
!     the routine determines first the euler angles correponding
!     to a symmetry operation. reflections are decomposed into
!     inversion + rotation for this reason.
!=======================================================================
!
    do isym = 1,nsym
      rmat(1:3,1:3) = rotmat(1:3,1:3,isymindex(isym))

      invflag(isym) = 0
      if (ddet33(rmat) < 0.d0) then
        call dscal(9,-1.d0,rmat,1)
        invflag(isym) = 1
      endif

!----------------------------------------------------------------------
      co2 = rmat(3,3)
      tet(2) = acos(co2)
      loop = 0
555   continue
      if (loop == 1) tet(2) = -tet(2)
      si2 = sin(tet(2))

      if (equal(co2, 1.d0)) then
        tet(1) = acos(rmat(1,1))
        if (.not. equal(rmat(1,2), sin(tet(1)))) then
          tet(1) = -tet(1)
          if (.not. equal(rmat(1,2), sin(tet(1)))) write(*,*) '>>>>>>>>>>>>>>> STRAGE 1'
        endif
        tet(2:3) = 0.d0
      elseif (equal(co2, -1.d0)) then
        tet(1) = acos(-rmat(1,1))
        if (.not. equal(rmat(1,2), -sin(tet(1)))) then
          tet(1) = -tet(1)
          if (.not. equal(rmat(1,2), -sin(tet(1)))) write(*,*) '>>>>>>>>>>>>>>> STRAGE 2'
        endif
        tet(2:3) = [pi, 0.d0]
      else
        tet(1) = acos(rmat(3,1)/si2)
        if (.not. equal(rmat(3,2), si2*sin(tet(1)))) then
          tet(1) = -tet(1)
          if (.not. equal(rmat(3,2), si2*sin(tet(1)))) write(*,*) '>>>>>>>>>>>>>>> STRAGE 3'
        endif

        tet(3) = acos(-rmat(1,3)/si2)
        if (.not. equal(rmat(2,3), si2*sin(tet(3)))) then
          tet(3) = -tet(3)
          if (.not. equal(rmat(2,3), si2*sin(tet(3)))) write(*,*) '>>>>>>>>>>>>>>> STRAGE 4'
        endif
      endif

      nok = 0
      do i1 = 1, 3
        do i2 = 1, 3
          if (checkrmat(rmat, cos(tet(1)), sin(tet(1)), cos(tet(2)), sin(tet(2)), cos(tet(3)), sin(tet(3)), i1, i2)) then
            nok = nok + 1
          elseif (loop < 1) then
            loop = loop + 1
            goto 555
          endif
        enddo ! i2
      enddo ! i1

      symeulang(1:3,isym) = tet(1:3)*(180.d0/pi)
      if (nok /= 9) write(*, fmt="(50('>'),' trouble in <SYMTAUMAT>',i3,f10.5)") nok
      write(*, fmt="(8x,i2,3x,a,i3,3f10.5)") isym, rotname(isymindex(isym)), invflag(isym), symeulang(1:3,isym)
    enddo ! isym
    write(6,'(8x,57(1h-),/)')
!
!-----------------------------------------------------------------------
!                    initialize all rotation matrices
!-----------------------------------------------------------------------
    call zcopy(nkmmax*nkmmax*nsym,c0,0,drot,1)
!-----------------------------------------------------------------------
!                       create rotation matrices
!-----------------------------------------------------------------------
!
    if (irel <= 2) then
      ireleff = 0
      nkeff = nl
    else
      ireleff = 3
      nkeff = nk
    endif
!
    do isym = 1,nsym
!
      call calcrotmat(nkeff,ireleff,symeulang(1,isym), symeulang(2,isym),symeulang(3,isym), drot(1,1,isym),fact,nkmmax)
!
    enddo
!-----------------------------------------------------------------------
!                     create matrix for inversion
!-----------------------------------------------------------------------
    call zcopy(nkmmax*nkmmax,c0,0,dinv,1)

    i = 0
    if (irel > 2) then
      ns = 2
    else
      ns = 1
    endif
    do l = 0, (nl-1)
      do m = 1, ns*(2*l+1)
        i = i + 1
        dinv(i,i) = (-1.d0)**l
      enddo ! m
    enddo ! l
    itop = i

!-----------------------------------------------------------------------
!                         include inversion
!-----------------------------------------------------------------------
    do isym = 1,nsym
      if (invflag(isym) /= 0) then

        call zgemm('n','n',nkm,nkm,nkm,c1,drot(1,1,isym),nkmmax,dinv,nkmmax,c0,w2,nkmmax)

        do j = 1, nkm
          call zcopy(nkm,w2(1,j),1,drot(1,j,isym),1)
        enddo ! j
      endif
    enddo ! isym
!
!-----------------------------------------------------------------------
!            add second spin-diagonal block for  irel=2
!            spin off-diagonal blocks have been initialized before
!-----------------------------------------------------------------------
    if (irel == 2) then
      nlm = nkm/2
      if (itop /= nlm) call errortrap('SYMTAUMAT',11,1)
      do isym = 1, nsym
        do j = 1, nlm
          call zcopy(nlm,drot(1,j,isym),1,drot(nlm+1,nlm+j,isym),1)
        enddo ! j
      enddo ! isym
    endif
!-----------------------------------------------------------------------
!            transform to real spherical representation for  krel=0
!-----------------------------------------------------------------------
    n = nkm
    m = nkmmax
    if (krel == 0) then
      do isym = 1, nsym
        call zgemm('n','n',n,n,n,c1,rc,m,drot(1,1,isym),m,c0,w1,m)
        call zgemm('n','c',n,n,n,c1,w1,m,rc,m,c0,drot(1,1,isym),m)
      enddo ! isym
    endif
!-----------------------------------------------------------------------
!                     create matrix for time reversal
!-----------------------------------------------------------------------
    if (irel > 1) then
      call zcopy(nkmmax*nkmmax,c0,0,dtim,1)

      i = 0
      do k = 1, nk
        l = k/2
        if (l*2 == k) then
          sk = -1d0
        else
          sk = +1d0
        endif
        rj = l + sk*0.5d0
!           do rmj = -rj, + rj ! real do loop iterator transformed
!                              ! to integer do loop iterator e.rabel
        do j = 0, (2*l + nint(sk))
          rmj = -rj + dble(j)
          i1 = nint(2*l*(rj + 0.5d0)+rj+rmj+1)
          i2 = nint(2*l*(rj + 0.5d0)+rj-rmj+1)
          dtim(i1,i2) = sk*(-1)**nint(rmj + 0.5d0)
        enddo ! j
      enddo ! k

    endif
!=======================================================================
!            set up of transformation matrices completed
!=======================================================================
!
!
!-----------------------------------------------------------------------
! for testing
!
!ccc      write(6,*) ' number of symmetries : ', nsym
!ccc      
!ccc      do isym = 1, nsym
!ccc         write(6,*) ' isym = ',isym
!ccc         call cmatstr('DROT',4,drot(1,1,isym),nkm,nkmmax,krel*3,krel*3,0,1d-12,6)
!ccc         write(6,*)
!ccc      enddo ! isym
!
!-----------------------------------------------------------------------
!
    if (iprint == 0) return

!=======================================================================
!       find the structure of the site-diagonal tau - matrices  tauq
!=======================================================================

    call taustruct(drot, nsym, nkm, nq, nqmax, nkmmax, iprint)

  endsubroutine symtaumat
  
  subroutine taustruct(drot, nsym, nkm, nq, nqmax, nkmmax, iprint)
!   ********************************************************************
!   *   find the structure of the site-diagonal tau - matrices  tauq   *
!   ********************************************************************
    integer, intent(in) :: iprint, nkm, nkmmax, nq, nqmax, nsym
    double complex, intent(in) :: drot(nkmmax,nkmmax,48)
    
    double complex, parameter :: c0=(0.d0,0.d0)
    double precision :: abst,x
    integer :: i,i0,imweight,iq,isym,iw,iwr,j,k,l,lin,nelmt, nkmq(nqmax), nkmtop,nlin,non0(nqmax)
    double complex :: st(nkmmax,nkmmax), tauk(nkmmax,nkmmax,nqmax)

    nkmq(1:nq) = nkm

    imweight = 0
    nelmt = 0
    nlin = 0
    iw = 6

    do iq = 1,nq
      non0(iq) = 0
      nkmtop = nkmq(iq)

      if (iprint > 0) write(6, fmt="(//,' ===========================================================',/,'   structure of  TAU-matrix   INT <i|t(k)|j>     IQ=',I3,/,' ===========================================================',/)") iq
      do i = 1, nkmtop
        do j = 1, nkmtop
          st(i,j) = 0.d0

          call zcopy(nkmmax*nkmmax*nqmax,c0,0,tauk,1)

          do isym = 1,nsym
            i0 = iq

            do l = 1, nkmtop
              do k = 1, nkmtop
                tauk(k,l,i0) = tauk(k,l,i0) + drot(i,k,isym)*dconjg(drot(j,l,isym))
              enddo ! k
            enddo ! l
            
          enddo ! isym

          lin = 0
          iwr = 0
          do k = 1, nkmq(iq)
            do l = 1, nkmq(iq)
              abst = abs(tauk(k,l,iq))
              st(i,j) = st(i,j) + abst
              if (abst > 1d-8) then
                if (dimag(tauk(k,l,iq)) > 1d-5) then
                  if (iprint > 0) write(*,*) ' Im(Weight) > 1D-5 ',i,j,k,l
                  imweight = 1
                endif
                x = dreal(tauk(k,l,iq))/dble(nsym)

                if (iprint > 1) then
                  if (iwr == 0) then
                    iwr = 1
                    write(iw, fmt="('     TAUQ(',I2,',',I2,',',I2,') =  ','   ',F10.4,' * <',I3,'|T(K)|',I3,'>')") i,j,iq,x,k + (iq-1)*nkm, l + (iq-1)*nkm
                  else
                    write(iw, fmt="(23X,' + ',F10.4,' * <',I3,'|T(K)|',I3,'>')") x,k + (iq-1)*nkm, l + (iq-1)*nkm
                  endif
                endif
                lin = lin + 1
              endif
            enddo ! l
          enddo ! k

          if (lin > 0) then
            nlin = nlin + lin
            nelmt = nelmt + 1
            non0(iq) = non0(iq) + 1
          endif

          if (abs(st(i,j)) > 1d-5) st(i,j) = 2

        enddo ! j
      enddo ! i
    enddo ! iq

    write(6, fmt="(/,5X,'non-0 TAU-elements          ',I5,'   Q:',80I4)") nelmt, non0(1:nq)
    write(6, fmt="(5X,'terms to sum up             ',I5,/)") nlin

    if (imweight /= 0) then
      write(*, fmt="(/,5X,50('#'),/,5X,'WARNING: complex TAU weights found',/,5X,'this may occur for rotated magnetic moments',/,5X,'relevant only for tetrahedron BZ-integration',/,5X,50('#'),/)")
      warn(6, "complex TAU weights found, this may occur for rotated magnetic moments (relevant only for tetrahedron BZ-integration).")
    endif ! imweight nonzero

  endsubroutine taustruct
 
  
  
  subroutine pointgrp(rotmat, rotname)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after 
! J.F. Cornwell (Group Theory??) second edition 
! Appendix D, p 324-325
! 
! *********************************************    
    double precision, intent(out) :: rotmat(3,3,64)
    character(len=*), intent(out) :: rotname(64)
    integer :: is
    double precision, parameter :: RTHREE = SQRT(3.d0)*.5d0, HALF=0.5d0, ONE=1.d0

    rotmat(:,:,:) =  0.d0 ! set all matrices to zero

    rotmat(1,1,1) =  ONE
    rotmat(2,2,1) =  ONE
    rotmat(3,3,1) =  ONE
    rotname(1) = 'E'

    rotmat(1,2,2) =  ONE
    rotmat(2,3,2) = -ONE
    rotmat(3,1,2) = -ONE
    rotname(2) = 'C3alfa'          

    rotmat(1,2,3) = -ONE
    rotmat(2,3,3) = -ONE
    rotmat(3,1,3) =  ONE
    rotname(3) = 'C3beta '

    rotmat(1,2,4) = -ONE
    rotmat(2,3,4) =  ONE
    rotmat(3,1,4) = -ONE
    rotname(4) = 'C3gamma'

    rotmat(1,2,5) = ONE
    rotmat(2,3,5) = ONE
    rotmat(3,1,5) = ONE
    rotname(5) = 'C3delta '

    rotmat(1,3,6) = -ONE
    rotmat(2,1,6) =  ONE
    rotmat(3,2,6) = -ONE
    rotname(6) = 'C3alfa-1'

    rotmat(1,3,7) =  ONE
    rotmat(2,1,7) = -ONE
    rotmat(3,2,7) = -ONE
    rotname(7) = 'C3beta-1 '

    rotmat(1,3,8) = -ONE
    rotmat(2,1,8) = -ONE
    rotmat(3,2,8) =  ONE
    rotname(8) = 'C3gamma-1'

    rotmat(1,3,9) =  ONE
    rotmat(2,1,9) =  ONE
    rotmat(3,2,9) =  ONE
    rotname(9) = 'C3delta-1'

    rotmat(1,1,10) =  ONE
    rotmat(2,2,10) = -ONE
    rotmat(3,3,10) = -ONE
    rotname(10) = 'C2x'
          
    rotmat(1,1,11) = -ONE
    rotmat(2,2,11) =  ONE
    rotmat(3,3,11) = -ONE
    rotname(11) = 'C2y'

    rotmat(1,1,12) = -ONE
    rotmat(2,2,12) = -ONE
    rotmat(3,3,12) =  ONE
    rotname(12) = 'C2z'

    rotmat(1,1,13) =  ONE
    rotmat(2,3,13) =  ONE
    rotmat(3,2,13) = -ONE
    rotname(13) = 'C4x'
          
    rotmat(1,3,14) = -ONE
    rotmat(2,2,14) =  ONE
    rotmat(3,1,14) =  ONE
    rotname(14) = 'C4y '

    rotmat(1,2,15) =  ONE
    rotmat(2,1,15) = -ONE
    rotmat(3,3,15) =  ONE
    rotname(15) = 'C4z'
          
    rotmat(1,1,16) =  ONE
    rotmat(2,3,16) = -ONE
    rotmat(3,2,16) =  ONE
    rotname(16) = 'C4x-1 '

    rotmat(1,3,17) =  ONE
    rotmat(2,2,17) =  ONE
    rotmat(3,1,17) = -ONE
    rotname(17) = 'C4y-1'

    rotmat(1,2,18) = -ONE
    rotmat(2,1,18) =  ONE
    rotmat(3,3,18) =  ONE
    rotname(18) = 'C4z-1'
          
    rotmat(1,2,19) =  ONE
    rotmat(2,1,19) =  ONE
    rotmat(3,3,19) = -ONE
    rotname(19) = 'C2a'

    rotmat(1,2,20) = -ONE
    rotmat(2,1,20) = -ONE
    rotmat(3,3,20) = -ONE
    rotname(20) = 'C2b'

    rotmat(1,3,21) =  ONE
    rotmat(2,2,21) = -ONE
    rotmat(3,1,21) =  ONE
    rotname(21) = 'C2c'

    rotmat(1,3,22) = -ONE
    rotmat(2,2,22) = -ONE
    rotmat(3,1,22) = -ONE
    rotname(22) = 'C2d'

    rotmat(1,1,23) = -ONE
    rotmat(2,3,23) =  ONE
    rotmat(3,2,23) =  ONE
    rotname(23) = 'C2e'

    rotmat(1,1,24) = -ONE
    rotmat(2,3,24) = -ONE
    rotmat(3,2,24) = -ONE
    rotname(24) = 'C2f'
    
    do is = 1, 24
      rotmat(1:3,1:3,is+24) = -rotmat(1:3,1:3,is)
      rotname(is+24) = 'I'//rotname(is)
    enddo ! i1      

!      
!*********************************************
! Trigonal and hexagonal groups
!*********************************************
!
    rotmat(1,1,49) = -HALF
    rotmat(1,2,49) =  RTHREE
    rotmat(2,1,49) = -RTHREE
    rotmat(2,2,49) = -HALF
    rotmat(3,3,49) =  ONE
    rotname(49) = 'C3z'  

    rotmat(1,1,50) = -HALF
    rotmat(1,2,50) = -RTHREE
    rotmat(2,1,50) =  RTHREE
    rotmat(2,2,50) = -HALF
    rotmat(3,3,50) =  ONE
    rotname(50) = 'C3z-1'

    rotmat(1,1,51) =  HALF
    rotmat(1,2,51) =  RTHREE
    rotmat(2,1,51) = -RTHREE
    rotmat(2,2,51) =  HALF
    rotmat(3,3,51) =  ONE
    rotname(51) = 'C6z'

    rotmat(1,1,52) =  HALF
    rotmat(1,2,52) = -RTHREE
    rotmat(2,1,52) =  RTHREE
    rotmat(2,2,52) =  HALF
    rotmat(3,3,52) =  ONE
    rotname(52) = 'C6z-1'

    rotmat(1,1,53) = -HALF
    rotmat(1,2,53) =  RTHREE
    rotmat(2,1,53) =  RTHREE
    rotmat(2,2,53) =  HALF
    rotmat(3,3,53) = -ONE
    rotname(53) = 'C2A'    

    rotmat(1,1,54) = -HALF
    rotmat(1,2,54) = -RTHREE
    rotmat(2,1,54) = -RTHREE
    rotmat(2,2,54) =  HALF
    rotmat(3,3,54) = -ONE
    rotname(54) = 'C2B'

    rotmat(1,1,55) =  HALF
    rotmat(1,2,55) = -RTHREE
    rotmat(2,1,55) = -RTHREE
    rotmat(2,2,55) = -HALF
    rotmat(3,3,55) = -ONE
    rotname(55) = 'C2C'

    rotmat(1,1,56) =  HALF
    rotmat(1,2,56) =  RTHREE
    rotmat(2,1,56) =  RTHREE
    rotmat(2,2,56) = -HALF
    rotmat(3,3,56) = -ONE
    rotname(56) = 'C2D'
    
    do is = 1, 8
      rotmat(1:3,1:3,56+is) = -rotmat(1:3,1:3,48+is)
      rotname(56+is) = 'I'//rotname(48+is) 
    enddo ! is
      
  endsubroutine pointgrp
  
  
  
  subroutine calcrotmat(nk, irel, alfdeg, betdeg, gamdeg, rot, fact, nkmmax)
    use Constants_mod, only: pi
!   ********************************************************************
!   *                                                                  *
!   *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
!   *           (ALFDEG, BETDEG, GAMDEG)                               *
!   *                                                                  *
!   *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
!   *            EQS. (4.8), (4.12) AND (4.13)                         *
!   *                                                                  *
!   *   for IREL=0,1   NK == NL           non-relativistic (l,m_l)     *
!   *       IREL=3     NK == odd          relativistic (kappa,mue)     *
!   *                                                                  *
!   *   12/11/96  HE  deal with beta = 0                               *
!   ********************************************************************
    integer, intent(in) :: nk, irel, nkmmax
    double precision, intent(in) :: alfdeg, betdeg, gamdeg
    
    double complex, parameter :: ci = (0.d0,1.d0), c0 = (0.d0,0.d0)
!   double precision, parameter :: pi = 3.141592653589793238462643d0   

    double precision :: num, msb05, msb05sq, msb05pw, j,m1,m2, rfac,x
    double precision :: fact(0:100)

    integer :: s, slow, shigh, off
    double complex :: emim2a, emim1g, rot(nkmmax,nkmmax) 
    integer :: k, l, nmue, im2, im1
    double precision :: cb05, cb05sq, sm, cb05pw, dom
!                       
! inline function    factorial for real argument
    rfac(x) = fact(nint(x))

    if(irel == 2) call errortrap('calcrotmat', 12, 1) 
    if(irel == 3 .and. mod(nk, 2) == 0) call errortrap('calcrotmat', 13, 1) 

    rot(:,:) = c0

    cb05    =  dcos(betdeg*0.5d0*pi/180.d0)
    cb05sq  =  cb05*cb05
    msb05   = -dsin(betdeg*0.5d0*pi/180.d0)
    msb05sq =  msb05*msb05
    
    off = 0
    do k = 1, nk 
      if (irel < 2) then
        l = k - 1
        j = l
      else
        l = k/2
        if (l*2 == k) then
          j = l - 0.5d0
        else
          j = l + 0.5d0
        endif         
      endif

      nmue = nint(2*j + 1)

      do im2 = 1, nmue
        m2 = - j + (im2 - 1.d0)
        emim2a = cdexp(-ci*m2*alfdeg*pi/180.d0)

          do im1 = 1, nmue
            m1 = - j + (im1 - 1.d0)
            emim1g = cdexp(-ci*m1*gamdeg*pi/180.d0)

            if (dabs(betdeg) < 1d-8) then
              if (im1 == im2) then 
                sm = 1.d0
              else
                sm = 0.d0
              endif
            else
              slow   = max(           0, nint(m1 - m2))
              shigh  = min(nint(j - m2), nint( j + m1))
              cb05pw =  cb05** nint(2*j + m1 - m2 - 2*slow + 2)
              msb05pw = msb05**nint(      m2 - m1 + 2*slow - 2)
              dom = (-1.d0)**(slow-1) * dsqrt(rfac(j+m1)*rfac(j-m1)*rfac(j+m2)*rfac(j-m2))
              sm = 0.d0

              do s = slow, shigh
                dom = -dom
                num = fact(s) * rfac(j-m2-s) * rfac(j+m1-s) * rfac(m2-m1+s)
                cb05pw  = cb05pw  /  cb05sq
                msb05pw = msb05pw * msb05sq
                sm = sm + (dom/num) * cb05pw * msb05pw
              enddo ! s
            endif

            rot(off+im2,off+im1) = emim1g * sm * emim2a
          enddo ! im1
        enddo ! im2

      off = off + nmue
    enddo ! k

  endsubroutine calcrotmat

  
  subroutine errortrap(routine, k, istop)
    character(len=*), intent(in) :: routine
    integer, intent(in) :: k, istop
    
    integer, parameter :: kmax=14
    character(len=60) :: text
    character(len=*), parameter :: t(kmax) = [ &
    'program KKRSCF called for TASK <> SCF    check input file   ', &
    'TASK = SCF in input   but not program KKRSCF called         ', &
    'ITEST should be <=4                                         ', &
    'DATASET not initialized           >>>> check input file     ', &
    'POTFIL  not initialized           >>>> check input file     ', &
    'SCFSTART for   PROGRAM <> KKRSCF  >>>>  call  kkrscf instead', &
    'SUM OF CONC <> 1                                            ', &
    'all SOCTL should be >= 0          >>>> check input file     ', &
    'all SOCTL should be < 0           >>>> check input file     ', &
    'use BZINT = POINTS    for SP-SREL case and for spin spirals ', &
    'I < > NLM   for IREL = 2                                    ', &
    'routine called for IREL = 2   deal with that case outside   ', &
    'routine called for IREL = 3   and   even  NK                ', &
    'anti-unitary symmetry matrices created for IREL = 2         ']

    if (k > 0 .and. k <= kmax) then
      text = t(k)
    else
      text = 'unknown reason   key out of bounds'
    endif

    if (istop == 1) then
      write(6, fmt="(/,1x,79('*'),/,1x,79('*'),/,5x,'STOP in subroutine <',a,'>',/,5x,a,/,1x,79('*'),/,1x,79('*'),/)") routine, text
      die_here(" in subroutine"+routine-", message=''"-text-"''.")
    else
      write(6, fmt="(/,1x,79('<'),/,5x,'warning from subroutine <',a,'>',/,5x,a,/,1x,79('>'),/)") routine, text
      warn(6, " in subroutine"+routine-", message=''"-text-"''.")
    endif

  endsubroutine errortrap
  
  
  logical function checkrmat(rmat, co1, si1, co2, si2, co3, si3, i, j)
!   ********************************************************************
!   *                                                                  *
!   *  check whether the values of the cosinus and sinus found for the *
!   *  euler angles tet1, tet2, tet3 are consistent with the           *
!   *  rotation matrix   rmat                                          *
!   *                                                                  *
!   ********************************************************************
    double precision, intent(in) :: co1, co2, co3, si1, si2, si3
    double precision, intent(in) :: rmat(3,3)
    integer, intent(in) :: i, j

    double precision :: a, b
    logical equal

    equal(a,b) = (abs(a-b) < 1.d-7)

    checkrmat = .false.

    if (i == 1) then
      if (j == 1) then
        checkrmat = equal(rmat(1,1), co3*co2*co1-si3*si1)
      elseif (j == 2) then
        checkrmat = equal(rmat(1,2), co3*co2*si1+si3*co1)
      elseif (j == 3) then
        checkrmat = equal(rmat(1,3), -co3*si2)
      endif
    elseif (i == 2) then
      if (j == 1) then
        checkrmat = equal(rmat(2,1), -si3*co2*co1-co3*si1)
      elseif (j == 2) then
        checkrmat = equal(rmat(2,2), -si3*co2*si1+co3*co1)
      elseif (j == 3) then
        checkrmat = equal(rmat(2,3), si3*si2)
      endif
    elseif (j == 1) then
      checkrmat = equal(rmat(3,1), si2*co1)
    elseif (j == 2) then
      checkrmat = equal(rmat(3,2), si2*si1)
    elseif (j == 3) then
      checkrmat = equal(rmat(3,3), co2)
    endif
  endfunction checkrmat
  
  
  logical function latvec(n, qlat, vec) 
!- Checks if a set of vectors are lattice vectors                       
! ----------------------------------------------------------------------
!i Inputs:                                                              
!i   n     :number of vectors                                           
!i   qlat  :primitive translation vectors in reciprocal space           
!i   vec   :double-precision vector                                     
!o Outputs:                                                             
!o   latvec:.true. if all vectors are lattice vectors                   
!r Remarks:                                                             
! ----------------------------------------------------------------------
    integer, intent(in) :: n 
    double precision, intent(in) :: qlat(3,3), vec(3,*) 

    integer :: i, m 
    double precision :: vdiff 
    double precision, parameter :: tol=1.d-3 
    
    latvec = .false. 
    do i = 1, n 
      do m = 1, 3 
        vdiff = vec(1,i)*qlat(1,m) + vec(2,i)*qlat(2,m) + vec(3,i)*qlat(3,m) 
        vdiff = dabs(vdiff - nint(vdiff)) 
        if (vdiff > tol) return ! false
      enddo ! m
    enddo ! i
    latvec = .true.
  endfunction latvec
  
endmodule Symmetry_mod      