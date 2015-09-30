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

    integer, parameter :: LMAX_LIMIT = 3 ! -1:no_limit
    
    integer :: ir, ifun, npan, meshn, nfun, lmmax, lmax, ilm, mfun
    character(len=32) :: filename
    character(len=8192) :: frmt
    character(len=*), parameter :: F9000="(16i5)", F9010="(4d20.12)" ! todo: define these formats with meaningful names in ShapefunData_mod
    double precision, allocatable :: thetas(:,:), xrn(:), drn(:), theta(:)
    integer, allocatable :: ilmsp(:), nm(:), ndigits(:)
    double precision :: factor
    character(len=*), parameter :: F01(0:1) = ['," 0"', '," ",f0.3']
    
    factor = 100./sqrt(16*atan(1.d0)) ! 1/sqrt(4pi) --> this lets the l=0 m=0-shapefunction become 100% inside the MT radius 
    
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
    allocate(ndigits(lmmax), theta(lmmax))
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
    enddo ! ir
      write(15, fmt='(F9.6,999a)') xrn(meshn)+drn(meshn), (' 0', ilm=1,lmmax) ! all shape function values should be zero here
      write(15, fmt='(F9.6,999a)') 2*xrn(meshn),          (' 0', ilm=1,lmmax) ! all shape function values should be zero here
    close(15)
    
  endsubroutine ! write_shapefun_file

endprogram

