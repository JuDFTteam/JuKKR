  ! **********************************************************************
  ! *                                                                    *
  ! *  generate lattice vectors of direct and reciprocal space from      *
  ! *  basic translation vectors br                                      *
  ! *                                                                    *
  ! *  alat            : lattice constant                                *
  ! *  br(i,j)         : i=x,y,z j= 1,2,3 bravais vectors                *
  ! *                    *** in a.u. ****                                *
  ! *  rmax            : maximum radius in real space        (input)     *
  ! *  gmax            : maximum radius in reciprocal space  (input)     *
  ! *  ngmax           : Number of reciprocal lattice vectors            *
  ! *  gn(3,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
  ! *  nrmax           : Number of real lattice vectors                  *
  ! *  rm(3,nmaxd)     : x,y,z  of real space vectors                    *
  ! *  nshlg           : shells in reciprocal space                      *
  ! *  nshlr           : shells in real space                            *
  ! *  nsg,nsr         : integer arrays, number of atoms in each shell   *
  ! *                                                                    *
  ! *  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
  ! *  one it is used only locally (GNR/RMR)       v.popescu May 2004    *
  ! *                                                                    *
  ! **********************************************************************

subroutine lattice3d(alat,bravais,recbv,ngmax,nrmax,nshlg,nshlr, &
     nsg,nsr,rmax,gmax,gn,rm,iprint,nmaxd,ishld, &
     print_info)

  implicit none
  !     ..
  !     .. Scalar arguments ..
  integer print_info
  integer iprint,ngmax,nrmax,nshlg,nshlr,nmaxd,ishld
  double precision alat
  !     ..
  !     .. Array arguments ..
  double precision bravais(3,3),recbv(3,3)
  double precision gn(3,nmaxd),rm(3,nmaxd)
  integer nsg(ishld),nsr(ishld)
  !     ..
  !     .. Local scalars ..
  double precision a,absgm,absrm,ag,ar,b,c,da,db,gmax,gx,gy,gz,pi, &
       rmax,rx,ry,rz,vmin
  integer i,k,l,m,n,n1,ng,nr,nsh,nshl,numg,numgh,numr, &
       numrh
  double precision dble
  integer idint
  !     ..
  !     .. Local arrays ..
  double precision absg(3),absr(3),bg(3,3),br(3,3),cj(4,nmaxd)
  double precision gnr(nmaxd),rmr(nmaxd)
  !     ..
  !     .. Intrinsic functions ..
  intrinsic abs,atan,dble,idint,max,mod,sqrt

  pi = 4.0d0*atan(1.0d0)
  !
  if (print_info.eq.0) then
     write (6,'(5X,2A,/)') '< LATTICE3D > : ', &
          'generating direct/reciprocal lattice vectors'
  endif
  !
  ! ======================================================================
  !
  ! --> basic trans. vectors and basis vectors
  !
  do i = 1,3
     br(1,i) = bravais(1,i)*alat
     br(2,i) = bravais(2,i)*alat
     br(3,i) = bravais(3,i)*alat
  end do
  ! ======================================================================
  !
  ! --> generate primitive vectors BG of reciprocal space
  !
  do i = 1,3
     bg(1,i) = recbv(1,i)*2d0*pi/alat
     bg(2,i) = recbv(2,i)*2d0*pi/alat
     bg(3,i) = recbv(3,i)*2d0*pi/alat
  end do
  ! ======================================================================
  !
  ! --> estimate no. of lattice vectors
  !
  do i = 1,3
     absr(i) = sqrt(br(1,i)**2+br(2,i)**2+br(3,i)**2)
     absg(i) = sqrt(bg(1,i)**2+bg(2,i)**2+bg(3,i)**2)
  end do
  !
  absrm = max(absr(1),absr(2),absr(3))
  absgm = max(absg(1),absg(2),absg(3))
  absrm = 2.0d0*pi/absrm
  absgm = 2.0d0*pi/absgm
  numr = 2*(idint(rmax/absgm)+1) + 1
  numg = 2*(idint(gmax/absrm)+1) + 1
  numrh = numr/2 + 1
  numgh = numg/2 + 1
  !
  ! **********************************************************************
  !                 generate lattice vectors of real space
  ! **********************************************************************
  !
  nr = 0
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do l = 1,numr
     a = dble(l-numrh)
     do m = 1,numr
        b = dble(m-numrh)
        do n = 1,numr
           c = dble(n-numrh)
           ! ----------------------------------------------------------------------
           rx = a*br(1,1) + b*br(1,2) + c*br(1,3)
           ry = a*br(2,1) + b*br(2,2) + c*br(2,3)
           rz = a*br(3,1) + b*br(3,2) + c*br(3,3)
           ar = sqrt(rx*rx+ry*ry+rz*rz)
           ! ----------------------------------------------------------------------
           if ( ar.le.rmax ) then
              nr = nr + 1
              if ( nr.gt.nmaxd ) then
                 write (6,*) &
                      ' ERROR: Dimension NMAXD in inc.p too small' &
                      ,nr,nmaxd
                 stop
              end if
              cj(1,nr) = rx
              cj(2,nr) = ry
              cj(3,nr) = rz
              cj(4,nr) = ar
           end if
        end do
     end do
  end do
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  nrmax = nr
  ! ======================================================================
  !
  ! --> sort vectors in order of increasing absolute value
  !
  da = 1.d-06
  nsh = 0
  nshl = -1
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do k = 1,nr
     vmin = rmax + 1.0d0
     do n = 1,nr
        if ( cj(4,n)-vmin.lt.0d0 ) then
           vmin = cj(4,n)
           n1 = n
        end if
     end do
     !
     nshl = nshl + 1
     rm(1,k) = cj(1,n1)
     rm(2,k) = cj(2,n1)
     rm(3,k) = cj(3,n1)
     rmr(k)  = cj(4,n1)
     db = vmin
     ! ----------------------------------------------------------------------
     if ( db.gt.da+1.d-06 ) then
        nsh = nsh + 1
        if ( nsh.gt.ishld ) then
           write (6,*) ' ERROR: Dimension ISHLD in inc.p too small', &
                nsh,ishld
           stop
        end if
        !
        nsr(nsh) = nshl
        nshl = 0
        da = db
     end if
     ! ----------------------------------------------------------------------
     cj(4,n1) = rmax + 1.0d0
  end do
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  nsh = nsh + 1
  nshl = nshl + 1
  if ( nsh.gt.ishld ) then
     write (6,*) ' ERROR: Dimension ISHLD in inc.p too small',nsh, &
          ishld
     stop
  end if
  !
  nsr(nsh) = nshl
  nshlr = nsh
  if ( nshlr.le.1 ) stop ' ERROR: cut-off radius RMAX too small '
  !
  ! **********************************************************************
  !
  ! **********************************************************************
  !                 generate lattice vectors of real space
  ! **********************************************************************
  !
  ng = 0
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do l = 1,numg
     a = dble(l-numgh)
     do m = 1,numg
        b = dble(m-numgh)
        do n = 1,numg
           c = dble(n-numgh)
           ! ----------------------------------------------------------------------
           gx = a*bg(1,1) + b*bg(1,2) + c*bg(1,3)
           gy = a*bg(2,1) + b*bg(2,2) + c*bg(2,3)
           gz = a*bg(3,1) + b*bg(3,2) + c*bg(3,3)
           ag = sqrt(gx*gx+gy*gy+gz*gz)
           ! ----------------------------------------------------------------------
           if ( ag.le.gmax ) then
              ng = ng + 1
              if ( ng.gt.nmaxd ) then
                 write (6,*) &
                      ' ERROR: Dimension NMAXD in inc.p too small' &
                      ,ng,nmaxd
                 stop
              end if
              cj(1,ng) = gx
              cj(2,ng) = gy
              cj(3,ng) = gz
              cj(4,ng) = ag
           end if
        end do
     end do
  end do
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ngmax = ng
  ! ======================================================================
  !
  ! --> sort vectors in order of increasing abs. value
  !
  da = 1.d-06
  nsh = 0
  nshl = -1
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do k = 1,ng
     vmin = gmax + 1.0d0
     do n = 1,ng
        if ( cj(4,n).lt.vmin ) then
           vmin = cj(4,n)
           n1 = n
        end if
     end do
     !
     nshl = nshl + 1
     gn(1,k) = cj(1,n1)
     gn(2,k) = cj(2,n1)
     gn(3,k) = cj(3,n1)
     gnr(k)  = cj(4,n1)
     db = vmin
     ! ----------------------------------------------------------------------
     if ( db.gt.da+1.d-07 ) then
        nsh = nsh + 1
        if ( nsh.gt.ishld ) then
           write (6,*) ' ERROR: Dimension ISHLD in inc.p too small', &
                nsh,ishld
           stop
        end if
        !
        nsg(nsh) = nshl
        nshl = 0
        da = db
     end if
     ! ----------------------------------------------------------------------
     cj(4,n1) = gmax + 1.0d0
  end do
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  nsh = nsh + 1
  nshl = nshl + 1
  if ( nsh.gt.ishld ) then
     write (6,*) ' ERROR: Dimension ISHLD in inc.p too small',nsh, &
          ishld
     stop
  end if
  !
  nsg(nsh) = nshl
  nshlg = nsh
  if ( nshlg.le.1 ) stop ' ERROR: cut-off radius GMAX too small '
  ! **********************************************************************
  !
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (print_info.eq.0) then
     write (6,fmt=99002)
     write (6,fmt=99003) 'Direct  lattice',nrmax,nshlr,rmr(nrmax)
     write (6,fmt=99003) 'Recipr. lattice',ngmax,nshlg,gnr(ngmax)
     write (6,fmt=99004)
  endif
  !
  if ( iprint.lt.3 ) return
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  !
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  k = 0
  write (6,fmt=99005) 'real-space'
  do l = 1,nshlr
     write (6,99006) l,nsr(l),rmr(k+1),(rm(m,k+1),m=1,3)
     do n = 2,nsr(l)
        write (6,fmt=99007) (rm(m,k+n),m=1,3)
     end do
     if ( l.ne.nshlr ) write (6,99008)
     k = k + nsr(l)
  end do
  write (6,99009)
  k = 0
  write (6,fmt=99005) 'reciprocal'
  do l = 1,nshlg
     write (6,99006) l,nsg(l),gnr(k+1),(gn(m,k+1),m=1,3)
     do n = 2,nsg(l)
        write (6,fmt=99007) (gn(m,k+n),m=1,3)
     end do
     if ( l.ne.nshlg ) write (6,99008)
     k = k + nsg(l)
  end do
  write (6,99009)
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  !
99002 format (10x,'               vectors  shells  max. R ',/,10x, &
       '               ------------------------------')
99003 format (10x,a,i7,2x,i6,2x,f9.5)
99004 format (10x,'               ------------------------------',/)
99005 format (10x,55('+'),/,18x,'generated ',a,' lattice vectors',/,10x, &
       55('+'),/,10x, &
       'shell Nvec    radius          x         y         z',/, &
       10x,55('-'))
99006 format (10x,i5,i5,f12.6,2x,3f10.5)
99007 format (34x,3f10.5)
99008 format (13x,52('-'))
99009 format (10x,55('+'),/)
end subroutine lattice3d
