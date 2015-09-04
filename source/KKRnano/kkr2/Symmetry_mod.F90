module Symmetry_mod
  implicit none
  private
  public :: pointgrp, findgroup, symtaumat

  contains
  
  subroutine findgroup(bravais, recbv, rbasis, nbasis, rsymat, rotname, isymindex, nsym, naezd)
! **********************************************************
! This subroutine finds the rotation matrices that leave the
! real lattice unchanged. 
! input:  bravais(i,j)    true bravais lattice vectors
!                         i = x,y,z ; j = A, B, C (a.u.)
!         recbv(i,j)      reciprocal basis vectors 
!         rbasis          coordinates of basis atoms
!         nbasis          number of basis atoms
!         rsymat          all 64 rotation matrices.
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
    double precision, intent(in) :: rsymat(64,3,3)
    character(len=*), intent(in) :: rotname(64)

    double precision, external :: ddot,ddet33
    logical, external :: test
    
    double precision :: r(3,4), rotrbas(3,naezd), bravais1(3,3)
    integer :: i,j,isym,i0,ia
    double precision mdotmp, mvecq(3,naezd), mvecqp(3,naezd)
    double precision mrotr(3,3), symdet, summdotmp
    double precision stet
    character(len=10) :: name(64)
    logical :: llatbas, latvec, lbulk

    write (6,fmt="(5x,'< FINDGROUP > : Finding symmetry operations',/)")

    nsym = 0
    bravais1(1:3,1:3) = bravais(1:3,1:3)
    
!     check for surface mode. if so, set bravais1(3,3) very large, so
!     that only the in-plane symmetries are found. not checked, be careful of z--> -z!
    lbulk=.true.
!     now check the bravais vectors if they have a z component 
    if (all(bravais(1:3,3) == 0.d0)) lbulk = .false.
!     
    do isym = 1, 64
!
!--------------------------------- store rotation matrix
!
      mrotr(1:3,1:3) = rsymat(isym,1:3,1:3)
      summdotmp = 0.d0
      symdet = ddet33(mrotr)
!
!     rotate bravais lattice vectors
!     
!     in the case of slab/interface geometry look only for 
!     symmetry opperations that preserve the z axis..
!     
      if (lbulk .or. rsymat(isym,3,3) == 1.d0) then 
!     do rotation only in case bulk or if slab and z axis is restored..  
          
        do i = 1, 3            ! loop on bravais vectors
          do j = 1, 3         ! loop on coordinates
            r(j,i) = rsymat(isym,j,1)*bravais1(1,i) + rsymat(isym,j,2)*bravais1(2,i) + rsymat(isym,j,3)*bravais1(3,i)
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
            rotrbas(j,ia) = rsymat(isym,j,1)*rbasis(1,ia) + rsymat(isym,j,2)*rbasis(2,ia) + rsymat(isym,j,3)*rbasis(3,ia)
            rotrbas(j,ia) = rotrbas(j,ia) - rbasis(j,ia)
            r(j,4) = rotrbas(j,ia) 
          enddo ! j
!     do j=1,4
!     do i=1,3
!     write(6,*) 'rrr',i,j,r(i,j)
!     enddo
!     enddo
!     write(6,*) 'latvec',latvec(4,recbv,r)
          if (.not. latvec(4,recbv,r)) llatbas = .false.
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
      endif ! (lbulk .or. (rsymat(isym,3,3) == 1))
    enddo ! isym = 1, nmatd
!     nsym symmetries were found
!     the isymindex array has the numbers of the symmetries found
 
    write(6,'(8x,60(1h-))')
    if (lbulk) then
      write(6,fmt="(8X,'3D symmetries',$)") 
    else  ! bulk
      write(6,fmt="(8X,'surface symmetries',$)") 
    endif ! bulk
    write(6,fmt="(' found for this lattice: ',i2,/,8x,60(1h-))") nsym
    do i = 1, nsym
      i0 = isymindex(i)
      name(i) = rotname(i0) 
    enddo ! i
    do i = 1, nsym/5 + 1
      isym = min(5, nsym - (i-1)*5)
      write(6,fmt="(8x,5(a10,2x))") (name(j),j=(i-1)*5+1,(i-1)*5+isym)
    enddo ! i
    write(6,fmt="(8x,60(1h-),/)")
    
  endsubroutine findgroup
  
  
  
  
  
  
  
!*==symtaumat.f    processed by SPAG 6.05Rc at 15:50 on 10 Dec 2002
      subroutine symtaumat(rotname, rotmat, drot, nsym, isymindex, nqmax, nkmmax, nq, nl, krel, iprint, nsymaxd)
!   ********************************************************************
!   *                                                                  *
!   *  Find the symmetry matrices DROT that act on t, tau, ....        *
!   *  KREL=0: for real spherical harmonics                            *
!   *  KREL=1: for relativistic represntation                          *
!   *                                                                  *
!   *  The NSYM allowed symmetry operations are indicated by ISYMINDEX *
!   *  in the table  ROTMAT. 
!   *                                                                  *
!   *  The routine determines first the Euler angles correponding      *
!   *  to a symmetry operation. Reflections are decomposed into        *
!   *  inversion + rotation for this reason.                           *
!   *                                                                  *
!   ********************************************************************
      integer, intent(in) :: iprint,krel,nkmmax,nl,nq,nqmax,nsym,nsymaxd
      double complex, intent(out) :: drot(nkmmax,nkmmax,48)
      integer, intent(in) :: isymindex(nsymaxd)
      double precision, intent(in) :: rotmat(64,3,3)
      character(len=*), intent(in) :: rotname(64)
      
      double complex, parameter :: ci=(0.d0,1.d0), c1=(1.d0,0.d0), c0=(0.d0,0.d0)
      double precision a,b,co1,co2,co3,det,fact(0:100),pi,rj,rmj,si1,si2,si3,sk, symeulang(3,48),tet(1:3)
      logical checkrmat
      double precision, external :: ddet33
      double complex dinv(nkmmax,nkmmax),dtim(nkmmax,nkmmax), rc(nkmmax,nkmmax),w1(nkmmax,nkmmax),w2(nkmmax,nkmmax)
      logical equal
      integer i,i1,i2,ind0q(nqmax),invflag(48),iq,irel,ireleff,isym, itop,j,k,l,loop,m,n,nk,nkeff,nkm,nlm,nok,ns
      double precision rmat(3,3)
      double precision w

      equal(a,b) = (abs(a-b) < 1d-7) ! inline function?

      write (6,99001)

      pi = 4.d0*atan(1.d0)

      irel = krel*3
      nk = (1-krel)*nl + krel*(2*nl-1)
      nkm = (1+krel)*nl**2
!
!-----------------------------------------------------------------------
      fact(0) = 1.d0
      do i = 1,100
         fact(i) = fact(i-1)*dble(i)
      enddo
!-----------------------------------------------------------------------
!
      ind0q(1) = 0
      do iq = 2,nq
         ind0q(iq) = ind0q(iq-1) + nkm
      enddo
!
! ----------------------------------------------------------------------
!    rc  transforms from  real to  complex (l,m,s) - representation
!                 |lc> = sum[lr] |lr> * rc(lr,lc)
! ----------------------------------------------------------------------
      if (krel == 0) then
         nlm = nkm
!
         call zcopy(nkmmax*nkmmax,c0,0,rc,1)
!
         w = 1.d0/sqrt(2.d0)
!
         do l = 0,(nl-1)
            do m = -l,l
               i = l*(l+1) + m + 1
               j = l*(l+1) - m + 1
!
               if (m < 0) then
                  rc(i,i) = -ci*w
                  rc(j,i) = w
               endif
               if (m == 0) then
                  rc(i,i) = c1
               endif
               if (m > 0) then
                  rc(i,i) = w*(-1.d0)**m
                  rc(j,i) = ci*w*(-1.d0)**m
               endif
            enddo
         enddo
      endif
!
!=======================================================================
!     the routine determines first the euler angles correponding
!     to a symmetry operation. reflections are decomposed into
!     inversion + rotation for this reason.
!=======================================================================
!
      do isym = 1,nsym
!
         do i1 = 1,3
            do i2 = 1,3
               rmat(i1,i2) = rotmat(isymindex(isym),i1,i2)
            enddo
         enddo
!
         det = ddet33(rmat)
!
         invflag(isym) = 0
         if (det < 0.d0) then
            call dscal(9,-1.d0,rmat,1)
            invflag(isym) = 1
         endif
!
!----------------------------------------------------------------------
         co2 = rmat(3,3)
         tet(2) = acos(co2)
         loop = 0
 50      continue
         if (loop == 1) tet(2) = -tet(2)
         si2 = sin(tet(2))
!
         if (equal(co2,1d0)) then
            tet(1) = acos(rmat(1,1))
            if (.not.equal(rmat(1,2),sin(tet(1)))) then
               tet(1) = -tet(1)
               if (.not.equal(rmat(1,2),sin(tet(1)))) write (*,*) '>>>>>>>>>>>>>>> STRAGE 1'
            endif
            tet(2:3) = 0.d0
         else if (equal(co2,-1d0)) then
            tet(1) = acos(-rmat(1,1))
            if (.not.equal(rmat(1,2),-sin(tet(1)))) then
               tet(1) = -tet(1)
               if (.not.equal(rmat(1,2),-sin(tet(1)))) write (*,*) '>>>>>>>>>>>>>>> STRAGE 2'
            endif
            tet(2:3) = [pi, 0.d0]
         else
            tet(1) = acos(rmat(3,1)/si2)
            if (.not.equal(rmat(3,2),si2*sin(tet(1)))) then
               tet(1) = -tet(1)
               if (.not.equal(rmat(3,2),si2*sin(tet(1)))) write (*,*) '>>>>>>>>>>>>>>> STRAGE 3'
            endif
!
            tet(3) = acos(-rmat(1,3)/si2)
            if (.not.equal(rmat(2,3),si2*sin(tet(3)))) then
               tet(3) = -tet(3)
               if (.not.equal(rmat(2,3),si2*sin(tet(3)))) write (*,*) '>>>>>>>>>>>>>>> STRAGE 4'
            endif

         endif

         co1 = cos(tet(1))
         si1 = sin(tet(1))
         co2 = cos(tet(2))
         si2 = sin(tet(2))
         co3 = cos(tet(3))
         si3 = sin(tet(3))
!
         nok = 0
         do i1 = 1,3
            do i2 = 1,3
               if (checkrmat(rmat,co1,si1,co2,si2,co3,si3,i1,i2)) then
                  nok = nok + 1
               else if (loop < 1) then
                  loop = loop + 1
                  goto 50
               endif
            enddo
         enddo
!
         symeulang(1,isym) = tet(1)*(180.d0/pi)
         symeulang(2,isym) = tet(2)*(180.d0/pi)
         symeulang(3,isym) = tet(3)*(180.d0/pi)
!
         if (nok /= 9) write (*,99009) nok
         write (*,99008) isym,rotname(isymindex(isym)),invflag(isym),(symeulang(i,isym),i=1,3)
!
      enddo
      write(6,'(8x,57(1h-),/)')
!
!-----------------------------------------------------------------------
!                    initialize all rotation matrices
!-----------------------------------------------------------------------
!
      call zcopy(nkmmax*nkmmax*nsym,c0,0,drot,1)
!
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
!
      i = 0
      if (irel > 2) then
         ns = 2
      else
         ns = 1
      endif
      do l = 0,(nl-1)
         do m = 1,ns*(2*l+1)
            i = i + 1
            dinv(i,i) = (-1.d0)**l
         enddo
      enddo
      itop = i
!
!-----------------------------------------------------------------------
!                         include inversion
!-----------------------------------------------------------------------
      do isym = 1,nsym
         if (invflag(isym) /= 0) then
!
            call zgemm('n','n',nkm,nkm,nkm,c1,drot(1,1,isym),nkmmax,dinv,nkmmax,c0,w2,nkmmax)
!
            do j = 1,nkm
               call zcopy(nkm,w2(1,j),1,drot(1,j,isym),1)
            enddo
         endif
      enddo
!
!-----------------------------------------------------------------------
!            add second spin-diagonal block for  irel=2
!            spin off-diagonal blocks have been initialized before
!-----------------------------------------------------------------------
      if (irel == 2) then
         nlm = nkm/2
         if (itop /= nlm) call errortrap('SYMTAUMAT',11,1)
         do isym = 1,nsym
!
            do j = 1,nlm
               call zcopy(nlm,drot(1,j,isym),1,drot(nlm+1,nlm+j,isym),1)
            enddo
         enddo
      endif
!-----------------------------------------------------------------------
!            transform to real spherical representation for  krel=0
!-----------------------------------------------------------------------
      n = nkm
      m = nkmmax
      if (krel == 0) then
         do isym = 1,nsym
            call zgemm('n','n',n,n,n,c1,rc,m,drot(1,1,isym),m,c0,w1,m)
            call zgemm('n','c',n,n,n,c1,w1,m,rc,m,c0,drot(1,1,isym),m)
         enddo
      endif
!-----------------------------------------------------------------------
!                     create matrix for time reversal
!-----------------------------------------------------------------------
      if (irel > 1) then
!
         call zcopy(nkmmax*nkmmax,c0,0,dtim,1)
!
         i = 0
         do k = 1,nk
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
               i1 = nint(2*l*(rj+0.5d0)+rj+rmj+1)
               i2 = nint(2*l*(rj+0.5d0)+rj-rmj+1)
               dtim(i1,i2) = sk*(-1)**nint(rmj+0.5d0)
            enddo
         enddo
!
      endif
!=======================================================================
!            set up of transformation matrices completed
!=======================================================================
!
!
!-----------------------------------------------------------------------
! for testing
!
!ccc      write (6,*) ' number of symmetries : ', nsym
!ccc      
!ccc      do isym = 1,nsym
!ccc         write(6,*) ' isym = ',isym
!ccc         call cmatstr('DROT',4,drot(1,1,isym),nkm,nkmmax,krel*3,krel*3,0,1d-12,6)
!ccc         write(6,*)
!ccc      enddo
!
!-----------------------------------------------------------------------
!
      if (iprint == 0) return
!
!=======================================================================
!       find the structure of the site-diagonal tau - matrices  tauq
!=======================================================================
!
      call taustruct(drot, nsym, nkm, nq, nqmax, nkmmax, iprint, irel)

      return
99001 format (5x,'<SYMTAUMAT> : rotation matrices acting on t/G/tau',//,8x,57(1h-),/,8x,'ISYM            INV          Euler angles      Unitarity',/,8x,57(1h-))
99008 format (8x,i2,3x,a,i3,3f10.5)
99009 format (50('>'),' trouble in <SYMTAUMAT>',i3,f10.5)
      endsubroutine symtaumat
  
  
  
  
  
!*==taustruct.f    processed by SPAG 6.05Rc at 15:50 on 10 Dec 2002
      subroutine taustruct(drot, nsym, nkm, nq, nqmax, nkmmax, iprint, irel)
!   ********************************************************************
!   *                                                                  *
!   *   find the structure of the site-diagonal tau - matrices  tauq   *
!   *                                                                  *
!   ********************************************************************
      integer, intent(in) :: iprint, irel, nkm, nkmmax, nq, nqmax, nsym
      double complex, intent(in) :: drot(nkmmax,nkmmax,48)
      
      double complex, parameter :: c0=(0.d0,0.d0)
      double precision :: abst,x
      integer :: i,i0,imweight,iq,isym,iw,iwr,j,k,l,lin,nelmt, nkmq(nqmax), nkmtop,nlin,non0(nqmax)
      double complex :: st(nkmmax,nkmmax), tauk(nkmmax,nkmmax,nqmax)

      do iq = 1,nq
        nkmq(iq) = nkm
      enddo

      imweight = 0
      nelmt = 0
      nlin = 0
      iw = 6

      do iq = 1,nq
         non0(iq) = 0
         nkmtop = nkmq(iq)

         if (iprint > 0) write (6,99004) iq
         do i = 1,nkmtop
            do j = 1,nkmtop
               st(i,j) = 0.d0

               call zcopy(nkmmax*nkmmax*nqmax,c0,0,tauk,1)

               do isym = 1,nsym
                  i0 = iq

                     do l = 1,nkmtop
                        do k = 1,nkmtop
                           tauk(k,l,i0) = tauk(k,l,i0) + drot(i,k,isym)*dconjg(drot(j,l,isym))
                        enddo
                     enddo
               enddo

               lin = 0
               iwr = 0
               do k = 1,nkmq(iq)
                  do l = 1,nkmq(iq)
                     abst = abs(tauk(k,l,iq))
                     st(i,j) = st(i,j) + abst
                     if (abst > 1d-8) then
                        if (dimag(tauk(k,l,iq)) > 1d-5) then
                           if (iprint > 0) write (*,*) ' Im(Weight) > 1D-5 ',i,j,k,l
                           imweight = 1
                        endif
                        x = dreal(tauk(k,l,iq))/dble(nsym)

                        if (iprint > 1) then
                           if (iwr == 0) then
                              iwr = 1
                              write (iw,99002) i,j,iq,x,k + (iq-1)*nkm, l + (iq-1)*nkm
                           else
                              write (iw,99003) x,k + (iq-1)*nkm, l + (iq-1)*nkm
                           endif
                        endif
                        lin = lin + 1
                     endif

                  enddo
               enddo

               if (lin > 0) then
                  nlin = nlin + lin
                  nelmt = nelmt + 1
                  non0(iq) = non0(iq) + 1
               endif

               if (abs(st(i,j)) > 1d-5) st(i,j) = 2

            enddo
         enddo
      enddo

      write (6,99005) nelmt,(non0(iq),iq=1,nq)
      write (6,99006) nlin

      if (imweight /= 0) write (*,99007)

      return
99002 format ('     TAUQ(',I2,',',I2,',',I2,') =  ','   ',F10.4,' * <',I3,'|T(K)|',I3,'>')
99003 FORMAT (23X,' + ',F10.4,' * <',I3,'|T(K)|',I3,'>')
99004 FORMAT (//,' ===========================================================',/,'   structure of  TAU-matrix   INT <i|t(k)|j>     IQ=',I3,/,' ===========================================================',/)
99005 FORMAT (/,5X,'non-0 TAU-elements          ',I5,'   Q:',80I4)
99006 FORMAT (5X,'terms to sum up             ',I5,/)
99007 FORMAT (/,5X,50('#'),/,5X,'WARNING: complex TAU weights found',/,5X,'this may occur for rotated magnetic moments',/,5X,'relevant only for tetrahedron BZ-integration',/,5X,50('#'),/)
      endsubroutine taustruct
  
  
  
  
  
  
  subroutine pointgrp(rotmat, rotname)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after 
! J.F. Cornwell (Group Theory??) second edition 
! Appendix D, p 324-325
! 
! *********************************************    
    double precision, intent(out) :: rotmat(64,3,3)
    character(len=*), intent(out) :: rotname(64)
    integer :: is
    double precision, parameter :: RTHREE = SQRT(3.d0)/2.d0, HALF=0.5d0, ONE=1.d0

    rotmat(:,:,:) =  0.d0 ! set all matrices to zero

    rotmat(1,1,1) =  ONE
    rotmat(1,2,2) =  ONE
    rotmat(1,3,3) =  ONE
    rotname(1) = 'E'

    rotmat(2,1,2) =  ONE
    rotmat(2,2,3) = -ONE
    rotmat(2,3,1) = -ONE
    rotname(2) = 'C3alfa'          

    rotmat(3,1,2) = -ONE
    rotmat(3,2,3) = -ONE
    rotmat(3,3,1) =  ONE
    rotname(3) = 'C3beta '

    rotmat(4,1,2) = -ONE
    rotmat(4,2,3) =  ONE
    rotmat(4,3,1) = -ONE
    rotname(4) = 'C3gamma'

    rotmat(5,1,2) = ONE
    rotmat(5,2,3) = ONE
    rotmat(5,3,1) = ONE
    rotname(5) = 'C3delta '

    rotmat(6,1,3) = -ONE
    rotmat(6,2,1) =  ONE
    rotmat(6,3,2) = -ONE
    rotname(6) = 'C3alfa-1'

    rotmat(7,1,3) =  ONE
    rotmat(7,2,1) = -ONE
    rotmat(7,3,2) = -ONE
    rotname(7) = 'C3beta-1 '

    rotmat(8,1,3) = -ONE
    rotmat(8,2,1) = -ONE
    rotmat(8,3,2) =  ONE
    rotname(8) = 'C3gamma-1'

    rotmat(9,1,3) =  ONE
    rotmat(9,2,1) =  ONE
    rotmat(9,3,2) =  ONE
    rotname(9) = 'C3delta-1'

    rotmat(10,1,1) =  ONE
    rotmat(10,2,2) = -ONE
    rotmat(10,3,3) = -ONE
    rotname(10) = 'C2x'
          
    rotmat(11,1,1) = -ONE
    rotmat(11,2,2) =  ONE
    rotmat(11,3,3) = -ONE
    rotname(11) = 'C2y'

    rotmat(12,1,1) = -ONE
    rotmat(12,2,2) = -ONE
    rotmat(12,3,3) =  ONE
    rotname(12) = 'C2z'

    rotmat(13,1,1) =  ONE
    rotmat(13,2,3) =  ONE
    rotmat(13,3,2) = -ONE
    rotname(13) = 'C4x'
          
    rotmat(14,1,3) = -ONE
    rotmat(14,2,2) =  ONE
    rotmat(14,3,1) =  ONE
    rotname(14) = 'C4y '

    rotmat(15,1,2) =  ONE
    rotmat(15,2,1) = -ONE
    rotmat(15,3,3) =  ONE
    rotname(15) = 'C4z'
          
    rotmat(16,1,1) =  ONE
    rotmat(16,2,3) = -ONE
    rotmat(16,3,2) =  ONE
    rotname(16) = 'C4x-1 '

    rotmat(17,1,3) =  ONE
    rotmat(17,2,2) =  ONE
    rotmat(17,3,1) = -ONE
    rotname(17) = 'C4y-1'

    rotmat(18,1,2) = -ONE
    rotmat(18,2,1) =  ONE
    rotmat(18,3,3) =  ONE
    rotname(18) = 'C4z-1'
          
    rotmat(19,1,2) =  ONE
    rotmat(19,2,1) =  ONE
    rotmat(19,3,3) = -ONE
    rotname(19) = 'C2a'

    rotmat(20,1,2) = -ONE
    rotmat(20,2,1) = -ONE
    rotmat(20,3,3) = -ONE
    rotname(20) = 'C2b'

    rotmat(21,1,3) =  ONE
    rotmat(21,2,2) = -ONE
    rotmat(21,3,1) =  ONE
    rotname(21) = 'C2c'

    rotmat(22,1,3) = -ONE
    rotmat(22,2,2) = -ONE
    rotmat(22,3,1) = -ONE
    rotname(22) = 'C2d'

    rotmat(23,1,1) = -ONE
    rotmat(23,2,3) =  ONE
    rotmat(23,3,2) =  ONE
    rotname(23) = 'C2e'

    rotmat(24,1,1) = -ONE
    rotmat(24,2,3) = -ONE
    rotmat(24,3,2) = -ONE
    rotname(24) = 'C2f'
    
    do is = 1, 24
      rotmat(is+24,1:3,1:3) = -rotmat(is,1:3,1:3)
      rotname(is+24) = 'I'//rotname(is)
    enddo ! i1      

!      
!*********************************************
! Trigonal and hexagonal groups
!*********************************************
!
    rotmat(49,1,1) = -HALF
    rotmat(49,1,2) =  RTHREE
    rotmat(49,2,1) = -RTHREE
    rotmat(49,2,2) = -HALF
    rotmat(49,3,3) =  ONE
    rotname(49) = 'C3z'  

    rotmat(50,1,1) = -HALF
    rotmat(50,1,2) = -RTHREE
    rotmat(50,2,1) =  RTHREE
    rotmat(50,2,2) = -HALF
    rotmat(50,3,3) =  ONE
    rotname(50) = 'C3z-1'

    rotmat(51,1,1) =  HALF
    rotmat(51,1,2) =  RTHREE
    rotmat(51,2,1) = -RTHREE
    rotmat(51,2,2) =  HALF
    rotmat(51,3,3) =  ONE
    rotname(51) = 'C6z'

    rotmat(52,1,1) =  HALF
    rotmat(52,1,2) = -RTHREE
    rotmat(52,2,1) =  RTHREE
    rotmat(52,2,2) =  HALF
    rotmat(52,3,3) =  ONE
    rotname(52) = 'C6z-1'

    rotmat(53,1,1) = -HALF
    rotmat(53,1,2) =  RTHREE
    rotmat(53,2,1) =  RTHREE
    rotmat(53,2,2) =  HALF
    rotmat(53,3,3) = -ONE
    rotname(53) = 'C2A'    

    rotmat(54,1,1) = -HALF
    rotmat(54,1,2) = -RTHREE
    rotmat(54,2,1) = -RTHREE
    rotmat(54,2,2) =  HALF
    rotmat(54,3,3) = -ONE
    rotname(54) = 'C2B'

    rotmat(55,1,1) =  HALF
    rotmat(55,1,2) = -RTHREE
    rotmat(55,2,1) = -RTHREE
    rotmat(55,2,2) = -HALF
    rotmat(55,3,3) = -ONE
    rotname(55) = 'C2C'

    rotmat(56,1,1) =  HALF
    rotmat(56,1,2) =  RTHREE
    rotmat(56,2,1) =  RTHREE
    rotmat(56,2,2) = -HALF
    rotmat(56,3,3) = -ONE
    rotname(56) = 'C2D'
    
    do is = 1, 8
      rotmat(56+is,1:3,1:3) = -rotmat(48+is,1:3,1:3)
      rotname(56+is) = 'I'//rotname(48+is) 
    enddo ! is
      
  endsubroutine pointgrp

endmodule Symmetry_mod      