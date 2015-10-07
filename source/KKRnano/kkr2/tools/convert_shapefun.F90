#ifndef nAtoms
#define nAtoms 1
#endif
program convert_shapefun
  implicit none
  
  integer :: ia
  
  do ia = 1, nAtoms
    call read_shapefun_file(ia)
  enddo ! ia
  
  contains
  
  subroutine read_shapefun_file(atom_id)
    integer, intent(in) :: atom_id
#ifndef LMAX_LIMIT
    integer, parameter :: LMAX_LIMIT = -1 ! -1:no_limit
#endif  
    
    integer :: ir, ifun, npan, meshn, nfun, lmmax, lmax, ilm, mfun, iq
    character(len=32) :: filename
    character(len=8192) :: frmt
    character(len=*), parameter :: F9000="(16i5)", F9010="(4d20.12)" ! todo: define these formats with meaningful names in ShapefunData_mod
    double precision, allocatable :: thetas(:,:), xrn(:), drn(:)
    integer, allocatable :: ilmsp(:), nm(:), ndigits(:)
    integer, parameter :: nq = 64, ells(16) = [0, 1, 1, 1, 2,2, 2, 2,2, 3,3,3, 3, 3,3,3]
    double precision :: factor, dr, dq, fq(16,0:nq-1), q, q2, r, r2, theta(16), fr(16)
    double precision, parameter :: pi = 4*atan(1.d0), SQRT2OVERPI=sqrt(2.d0/pi)
    character(len=*), parameter :: F01(0:1) = ['," 0"', '," ",f0.3']
    
    factor = 100./sqrt(4*pi) ! 1/sqrt(4pi) --> this lets the l=0 m=0-shapefunction become 100% inside the MT radius 
    
    write(filename, fmt="(a,i7.7)") "shape.",atom_id

    open(14, file=filename, form="formatted", action='read')
    read(14, fmt=F9000) npan, meshn
    allocate(nm(npan), xrn(meshn), drn(meshn))
    read(14, fmt=F9000) nm(1:npan)
    read(14, fmt=F9010) (xrn(ir), drn(ir), ir=1,meshn) ! interleaved pairs (r, dr)
    read(14, fmt=F9000) nfun
    allocate(ilmsp(nfun), thetas(1:meshn,nfun))
    do ifun = 1, nfun
      read(14, fmt=F9000) ilmsp(ifun)
      read(14, fmt=F9010) thetas(1:meshn,ifun)
    enddo ! ifun
    close(14)

    lmmax = maxval(ilmsp)
    lmax = floor(sqrt(0.999*lmmax))
    mfun = nfun
    
    if (LMAX_LIMIT >= 0 .and. lmax > LMAX_LIMIT) then
      lmax = LMAX_LIMIT ! limit to
      mfun = 0
      do ifun = 1, nfun
        if (ilmsp(ifun) <= (lmax+1)**2) mfun = max(mfun, ifun)
      enddo ! ifun
    endif ! limit lmax
    
    lmmax = (lmax+1)**2 ! correct
    write(*, fmt=*) lmax
    allocate(ndigits(lmmax))
    ndigits = 0        ! format string ', 0'
    ndigits(ilmsp(:mfun)) = 1 ! format string '," ",f0.6'
    
    write(filename, fmt="(a,i7.7)") "plot_shape.",atom_id
    
    write(unit=frmt, fmt='(999a)') '(F9.6', (trim(F01(ndigits(ilm))), ilm=1,lmmax), ',9(" ",f0.6))'
!   write(*, fmt=*) trim(frmt) ! show format
    
    open(15, file=filename, form="formatted", action='write')
    !! output to stdout
    write(15,'(/,9(a,i0))') '# atom_id=',atom_id,' mask values in percent'
      write(15, fmt='(F9.6,999a)') 0.0, ' 100.000', (' 0', ilm=2,lmmax) ! all shape function values should be zero here except the first
    do ir = 1, meshn
      write(15, fmt=frmt) xrn(ir), thetas(ir,:mfun)*factor
!       write(8, fmt=*) xrn(ir), 1./drn(ir)
    enddo ! ir
      write(15, fmt='(F9.6,999a)') xrn(meshn)+drn(meshn), (' 0', ilm=1,lmmax) ! all shape function values should be zero here
      write(15, fmt='(F9.6,999a)') 2*xrn(meshn),          (' 0', ilm=1,lmmax) ! all shape function values should be zero here
    close(15)

    
    if (lmax > 3) then
      lmax = 3 ! limit to
      mfun = 0
      do ifun = 1, nfun
        if (ilmsp(ifun) <= (lmax+1)**2) mfun = max(mfun, ifun)
      enddo ! ifun
    endif ! limit lmax
    
    dr = (xrn(meshn) - xrn(1))/nq
    dq = pi/dr
    
    fq(:,:) = 0.d0
!     ! inner region
!     do ir = 1, xrn(ir)/dr -1
!       r = ir*dr
!       r2 = r*r
!       do iq = 0, nq-1
!         q = iq*dq
!         fq(1,iq) = fq(1,iq) + thetas(ir,1)*r2*dr*Jell(0, x=q*r)
!       enddo ! iq
!     enddo ! ir
    
    ! interstitial region
    do ir = 1, meshn
      theta = 0.d0
      theta(ilmsp(:mfun)) = thetas(ir,:mfun)
      r = xrn(ir) - xrn(1)
      r2 = r*r
      do iq = 0, nq-1
        q = iq*dq
        fq(:,iq) = fq(:,iq) + theta(:)*r2*drn(ir)*Jell(ells, x=q*r)
      enddo ! iq
    enddo ! ir
    
    fq = fq*SQRT2OVERPI ! scale
    
    do iq = 0, nq-1
      write(*,'(99F16.6)') iq*dq, fq(:,iq) ! show coefficients in Bessel space
    enddo ! iq
    
!       fq(iq) = sum( rf(:nr)*mrg(:)*rg(:nr)**(2+2*ell-ellin)*Jellxmell( ell, x=q*rg(:nr) )*drg(:nr) )*q**ell*SQRT2OVERPI
!       fr(ir) = sum( fq(:)*qgd(:)**(2+ell)*Jellxmell( ell, x=qgd(:)*r ) )*dq*SQRT2OVERPI
    
    write(filename, fmt="(a,i7.7)") "filtered_shape.",atom_id
    
    open(17, file=filename, form="formatted", action='write')
    !! output to stdout
    write(17,'(/,9(a,i0))') '# atom_id=',atom_id,' mask values in percent'
    do ir = 1, meshn
      r = dr*(ir-1) ! synthetic equidistant grid
      fr = 0.d0
      do iq = 0, nq-1
        q = iq*dq
        q2 = q*q
        fr(:) = fr(:) + fq(:,iq)*q2*dq*Jell(ells, x=q*r)
      enddo ! iq
      write(17, fmt='(999(f0.6," "))') xrn(1) + r, fr*SQRT2OVERPI
    enddo ! ir
    close(17)
    
    
  endsubroutine ! write_shapefun_file


!   !! Bessel function J_ell(x) multipied with x^(-ell)
!   double precision elemental function Jellxmell( ell, x ) result( J )
!     integer, intent(in) :: ell
!     double precision, intent(in) :: x
! 
!     double precision, parameter :: THRESHOLD = 1.E-8, STARTVAL(0:3) = [1., 3., 15., 105.]
! 
!     if( abs( x ) < THRESHOLD ) then
!       j = 1./STARTVAL(ell)
!       return
!     endif
!     selectcase( ell )
! !     case(:-1) ; stop 'tBFUN Jellxmell: ell < 0 unphysical!'
!     case( 0 ) ; J = sin(x)/x                                           ! j0(x)
!     case( 1 ) ; J = (sin(x)-x*cos(x))/x**3                             ! j1(x) / x
!     case( 2 ) ; J = ((3.-x*x)*sin(x)-3.*x*cos(x))/x**5                 ! j2(x) / x^2
!     case( 3 ) ; J = ((15.-6.*x*x)*sin(x)-(15.-x*x)*x*cos(x))/x**7      ! j3(x) / x^3
! !     case( 4:) ; stop 'tBFUN Jellxmell: ell > 3 not implemented!'
!     endselect ! ell
!   endfunction ! Jellxmell


  !! Bessel function J_ell(x) multipied with x^(-ell)
  double precision elemental function Jell( ell, x ) result( J )
    integer, intent(in) :: ell
    double precision, intent(in) :: x

    double precision, parameter :: THRESHOLD = 1.E-8, STARTVAL(0:3) = [1., 3., 15., 105.]

    if( abs( x ) < THRESHOLD ) then
      j = x**ell/STARTVAL(ell)
      return
    endif
    selectcase( ell )
!     case(:-1) ; stop 'tBFUN Jell: ell < 0 unphysical!'
    case( 0 ) ; J = sin(x)/x                                           ! j0(x)
    case( 1 ) ; J = (sin(x)-x*cos(x))/x**2                             ! j1(x)
    case( 2 ) ; J = ((3.-x*x)*sin(x)-3.*x*cos(x))/x**3                 ! j2(x)
    case( 3 ) ; J = ((15.-6.*x*x)*sin(x)-(15.-x*x)*x*cos(x))/x**4      ! j3(x)
!     case( 4:) ; stop 'tBFUN Jell: ell > 3 not implemented!'
    endselect ! ell
  endfunction ! Jell
  
endprogram

