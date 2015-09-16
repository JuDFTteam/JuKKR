module PositionReader_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: getAtomData !, readXYZfile, PSE, atomic_number_by_symbol

  ! element symbols
  character(len=2), parameter :: PSE(-1:116) = [ & !! element symbol from the periodic tabel of elements
  'e ', '__', &                                                    ! vacuum
  'H ',                              'He', &                       ! 1s
  'Li','Be','B ','C ','N ','O ','F ','Ne', &                       ! 2s, 2p
  'Na','Mg','Al','Si','P ','S ','Cl','Ar',&                        ! 3s, 3p
  'K ','Ca', &                                                     ! 4s
        'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &       ! 3d
            'Ga','Ge','As','Se','Br','Kr', &                       ! 4p
  'Rb','Sr', &                                                     ! 5s
        'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', &       ! 4d
            'In','Sn','Sb','Te','I ','Xe', &                       ! 5p
  'Cs','Ba', &                                                     ! 6s
        'La', &                                                    ! 5d
                            'Ce','Pr','Nd','Pm','Sm','Eu','Gd', &  ! 4f Lanthanides
                            'Tb','Dy','Ho','Er','Tm','Yb','Lu', &  ! 4f Lanthanides
            'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &        ! 5d
            'Tl','Pb','Bi','Po','At','Rn', &                       ! 6p
  'Fr','Ra', &                                                     ! 7s
        'Ac', &                                                    ! 7p
                            'Th','Pa','U ','Np','Pu','Am','Cm', &  ! 5f Actinides
                            'Bk','Cf','Es','Fm','Md','No','Lr', &  ! 5f Actinides
            'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn', &        ! 6d
            '+-','Fl','++','Lv']                                   ! 7p and custom
  
  
  integer(kind=1), parameter :: & ! enum-replacement
    MODIFIED_INITIALIZED  = 0,  & ! set to zero
    MODIFIED_AUTO_DEFAULT = 1,  & ! use hard coded default values for elements
    MODIFIED_ELEM_DEFAULT = 2,  & ! use hard coded default values for elements
    MODIFIED_CONF_DEFAULT = 3,  & ! use default values for elements which have been adjusted in the input file
    MODIFIED_SAME_DEFAULT = 4,  & ! values given in the .xyz file have the same value as the MODIFIED_CONF_DEFAULT level
    MODIFIED_SPECIALIZED  = 5     ! use values given in the .xyz file
  character(len=*), parameter :: MODIFY_STRING(0:5) = ['initialized', 'automatic', 'element', 'configured', 'same', 'specialized']
  
  
  contains
  
  
#define MPI_master(COMM) .true.

  integer function getAtomData(filename, natoms, pos, comm) result(ist)
!   use ConfigReader_mod, only: ConfigReader, create, destroy, parseFile, getValue
    character(len=*), intent(in) :: filename
    integer, intent(out) :: natoms
    double precision, allocatable, intent(out) :: pos(:,:) ! pos(0:3,natoms)
    integer, intent(in) :: comm ! MPI communicator

    integer,          parameter :: nKeys=6, nDoF=12
    character(len=*), parameter :: keywords(nKeys) = ['Z=','Vw=','rMT=','Mag=','core=','start='] ! all strings need to tail with '='
    integer         , parameter :: ikeypos(nKeys)  = [ 0,    4,     5,     6,      7,      11  ] ! position in params(0:)
    double precision :: z_defaults(0:nDoF-1,-1:116)
    integer(kind=1)  :: z_modified(0:nDoF-1,-1:116)
    integer :: iZ
    
!   type(ConfigReader) :: cr
    character(len=16)  :: element_keyword
    character(len=256) :: configuration_line, elem_file

    z_defaults(:,:) = 0.d0; z_modified(:,:) = MODIFIED_INITIALIZED
    
    elem_file = 'elem' ! try to find the element configurations in the .xyz-file
    
!     call create(cr)
!     if (parseFile(cr, elem_file) /= 0) die_here('parsing file "'-elem_file-'" failed!')
    
    do iZ = -1, 116
      z_defaults(0,:) = dble(iZ); z_modified(0,:) = MODIFIED_AUTO_DEFAULT
      ! add more default functions of iZ here
      ! ...
      write(unit=element_keyword, fmt='(9a)') 'element-',PSE(iZ)
      configuration_line = '' ! init since getValue is intent(inout)
!     ist = getValue(cr, element_keyword, configuration_line, def='')
      if (ist == 0) then
        ist = parseParams(configuration_line, z_defaults(0:,iZ), z_modified(0:,iZ), keywords, ikeypos, &
                           linenumber=0, filename=elem_file, modify_states=[MODIFIED_ELEM_DEFAULT, MODIFIED_CONF_DEFAULT])
      endif
    enddo ! iZ

!     call destroy(cr)
    
    ist = readXYZfile(filename, natoms, pos, keywords, ikeypos, z_defaults, z_modified, comm)
    
  endfunction ! getAtomData


  integer function readXYZfile(filename, natoms, pos, keywords, ikeypos, z_defaults, z_modified, comm) result(ist)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: natoms
    double precision, allocatable, intent(out) :: pos(:,:) ! pos(0:3,natoms)
    double precision, intent(in) :: z_defaults(0:,-1:)
    integer(kind=1),  intent(in) :: z_modified(0:,-1:)
    integer, intent(in) :: comm ! MPI communicator

    character(len=*), intent(in) :: keywords(:) ! e.g. ['Z=','Vw=','rMT=','Mag=','core=','start='] ! all strings need to tail with '='
    integer         , intent(in) :: ikeypos(:)  ! e.g. [ 0,    4,     5,     6,      7,      11  ] ! position in params(0:)
    
    integer, parameter :: FU=23
    integer, parameter :: L=128
    character(len=L) :: line
    character(len=L+32) :: longline
    character(len=4)   :: symbol ! chemical symbol
    double precision   :: xyz(3) ! coordinates
    integer :: ia, iZ, n_warn_line_truncated, n_warn_ignored_param, ndof, ios(4), nadd, im
    logical, parameter :: checks = .true.
    integer(kind=1), allocatable :: modified(:,:) ! (0:ndof-1,natoms)
    double precision, parameter :: IMPOSSIBLE_VALUE = -(2.d0**17)
    
    n_warn_line_truncated = 0
    n_warn_ignored_param = 0
    
    deallocate(pos, stat=ist) ! ignore status
    
    if (MPI_master(comm)) then
      open(unit=FU, file=filename, status='old', action='read', iostat=ios(1))
      if (ios(1) /= 0) die_here('failed to open "'-filename-'"!')
      
      ! read number of atoms (first line in an .xyz file)
      read(unit=FU, fmt='(a)', iostat=ios(2)) line
      if (ios(2) /= 0) die_here('failed to read the first line in "'-filename-'"!')
      
      read(unit=line, fmt=*, iostat=ios(3)) natoms
      write(*,'(9(a,i0))') 'nAtoms = ',natoms
      if (ios(3) /= 0) die_here('failed to find the number of atoms in the first line of "'-filename-'", line reads "'-line-'"!')
      
      read(unit=FU, fmt='(a)', iostat=ios(4)) line ! second line in an .xyz file
      write(*,'(9a)') 'comment = "',trim(line),'"' ! echo comment line
      if (ios(4) /= 0) warn(6, 'unable to read a comment in line #2 of file "'-filename-'"!')
      
    endif ! master

    ndof = size(z_defaults, 1)
      
    ! broadcast ios(1:4)
    ! broadcast natoms
    ! broadcast ndof
    
    if (any(ios(1:3) /= 0) .or. natoms < 1) return
    
    allocate(pos(0:ndof-1,natoms), modified(0:ndof-1,natoms), stat=ist)
    if (ist /= 0) die_here('failed allocate position array with'+natoms+'x'+ndof+'entries (reading file "'-filename-'")!')
    
    if (MPI_master(comm)) then
      
      allocate(modified(0:ndof-1,natoms), stat=ist)
      pos(:,:) = 0.d0; modified(:,:) = MODIFIED_INITIALIZED

      do ia = 1, natoms
        longline = ''
        read(unit=FU, fmt='(a)', iostat=ist) longline
        if (ist /= 0) return ! error

        line = adjustl(longline) ! possible data loss due to string truncation here
        if (line /= adjustl(longline)) then
          write(0,'(a,i0,9a)') 'Warning: line #',ia+2,' in file "',trim(filename),'" was truncated!' ! DEBUG
          n_warn_line_truncated = n_warn_line_truncated + 1
        endif
        
        symbol = ''
        xyz = 0.d0
        read(unit=line, fmt=*, iostat=ist) symbol, xyz(1:3)
        if (ist /= 0) then
          warn(6, 'expected a chemical symbol and 3 positions in line #'-(ia+2)+'of file "'-filename-'"!')
          return ! error
        endif
        
!       write(*,'(9(a,i0,2x))') 'atom #',ia,  ' ist=',ist ! status will be -1 if no configuration string is shown
!       write(*,'(a,3F16.9,2x,a)') symbol, xyz(1:3)
        
        iZ = atomic_number_by_symbol(symbol)
        if (iZ < -1) die_here('cannot read chemical symbol "'-symbol-'" in line #'-(ia+2)+'of file "'-filename-'"!')
        
        if (iZ > ubound(z_defaults, 2)) then
          pos(0,ia) = dble(iZ)           ; modified(0,ia) = MODIFIED_AUTO_DEFAULT
          die_here('chemical symbol "'-symbol-'" in line #'-(ia+2)+'of file "'-filename-'" is higher than defaults are given!')
        else
          pos(0:,ia) = z_defaults(0:,iZ) ; modified(0:,ia) = z_modified(0:,iZ)
        endif
        
        ist = parseParams(line, pos(0:,ia), modified(0:,ia), keywords, ikeypos, linenumber=ia+2, filename=filename, n_warn_ignored=n_warn_ignored_param)
        if (ist /= 0) then
          warn(6, 'parsing of additional atom data in line #'-(ia+2)+'of file "'-filename-'"!')
          return ! error
        endif
        
        pos(1:3,ia) = xyz(1:3); modified(1:3,ia) = MODIFIED_SPECIALIZED ! set Cartesian positions

!       write(*, '(i0,99("  ",f0.2))') ia, pos(:,ia) ! show all params  
        
      enddo ! ia

      if (n_warn_line_truncated > 0) warn(6, 'in file "'-filename-'",'+n_warn_line_truncated+'lines have been truncated!')
      if (n_warn_ignored_param  > 0) warn(6, 'in file "'-filename-'",'+n_warn_ignored_param+'expressions have been ignored!')
      
      if (checks) then
        nadd = 0
        ia = natoms
          read(unit=FU, fmt=*, iostat=ist) symbol, xyz(1:3) ! read one more
        do while (ist == 0)
          nadd = nadd + 1
          read(unit=FU, fmt=*, iostat=ist) symbol, xyz(1:3) ! read even more lines
        enddo ! ia
        if (nadd > 0) warn(6, 'at least'+nadd+'additional lines in file "'-filename-'" seem to be valid atom entries!')
      endif ! checks
      
      close(FU, iostat=ist) ! ignore status
      write(*,'(9a)') 'file "',trim(filename),'" has been read in.' ! success message, todo suppress according to verbosity level
      
      ! do some statistics on modified before its deallocation
      do im = MODIFIED_INITIALIZED, MODIFIED_SPECIALIZED
        write(*,'(9(3a,i9))') 'in file "',trim(filename),'"',count(modified == im),' atom data items are at level "',trim(MODIFY_STRING(im)),'".'
      enddo ! im
      deallocate(modified, stat=ist)
      
    else
      pos = IMPOSSIBLE_VALUE ! DEBUG
    endif ! master
    
    ! todo: broadcast positions here
    
    ! DEBUG: assert( all(pos /= IMPOSSIBLE_VALUE) )
    
  endfunction ! read
  
  integer function parseParams(line, params, modi, keyw, ikey, linenumber, filename, n_warn_ignored, modify_states) result(ist)
    character(len=*),     intent(in) :: line
    double precision,  intent(inout) :: params(0:)
    integer(kind=1),     intent(out) :: modi(0:) ! has been modified? which level
    character(len=*),     intent(in) :: keyw(:) ! e.g. ['Z=','Vw=','rMT=','Mag=','core=','start='] ! all strings need to tail with '='
    integer,              intent(in) :: ikey(:) ! e.g. [ 0,    4,     5,     6,      7,      11  ] ! position in params(0:)
    integer,              intent(in) :: linenumber
    character(len=*),     intent(in) :: filename
    integer, intent(inout), optional :: n_warn_ignored
    integer(kind=1), intent(in), optional :: modify_states(1:2)

    character(len=16) :: bad
    double precision :: value
    integer :: ieq(size(keyw)), ii, ik, jc, nk
    integer(kind=1) :: ms(1:2)

    ms(1:2) = [4, 5]; if (present(modify_states)) ms = modify_states
    
    nk = size(keyw)
    assert( nk == size(ikey) ) ! the same number has to be passed
    if (any(ikey(1:nk) > ubound(params, 1))) die_here('invalid entries for iKey could lead to array-out-of-bounds! iKey='+ikey)
    
    ist = 0
    ieq(:) = 0
    do ik = 1, nk
      ii = index(line, trim(keyw(ik)))
      if (ii > 0) then ! the keyword has been found in line
        ii = ii + len_trim(keyw(ik)) ! start to read after the keyword
        read(unit=line(ii:), fmt=*, iostat=ist) value
        if (ist /= 0) then 
          ! reading error
          write(*, fmt='(/,9a)')  'string: "',trim(line),'"'
          write(*, fmt='(999a)')  '  invalid ',('-', jc=3,ii),'^ here for key "',trim(keyw(ik)),'"!'
          bad = adjustl((line(max(1, ii-7):min(ii+7, len(line)))))
          write(*, fmt='(3a,i0)') '  found in file "',trim(filename),'" at line #',linenumber
          die_here('value cannot be parsed at "...'-bad-'..." in line #'-linenumber+'in file "'-filename-'"!')
          return ! if error was treated as soft error
        endif
!!!     write(*, fmt='(a,f0.3)') trim(keyw(ik)), value ! show found values
        modi(ikey(ik)) = ms(2)
        if (params(ikey(ik)) == value) modi(ikey(ik)) = ms(1)
        params(ikey(ik)) = value ! store
!!!     write(*,*) mk, line(ieq(ik):ieq(ik)) ! confirm that ieq(ik) is the position of an equality char
        ieq(ik) = ii - 1 ! store the position of the equality char in line
      endif ! found
    enddo ! ik

    do ii = 1, len(line) ! loop through the string
      if (line(ii:ii) /= '=') cycle ! look only for equality chars
      if (any(ieq(:) == ii)) cycle ! this equality char has been found with a keyword
      ! launch a warning
      write(*, fmt='(/,9a)')     'string: "',trim(line),'"' ! WARNING
      write(*, fmt='(999a)')     '  ignored ',('-', jc=3,ii),'^ here! Warning!'
      bad = adjustl((line(max(1, ii-7):min(ii+7, len(line)))))
      write(*, fmt='(5a,i0,a,/)') '  expression "...',trim(bad),'..." found in file "',trim(filename),'" at line #',linenumber,' was ignored!'
      if (present(n_warn_ignored)) then
        n_warn_ignored = n_warn_ignored + 1 ! warnings are collected
      else ! warnings are launched
        warn(6, 'expression "...'-bad-'..." found in file "'-filename-'" at line #'-linenumber+'was ignored!')
      endif
    enddo ! ii

!   write(*, fmt='(99(2a,f0.3))') ('  ',trim(keyw(ik)),params(ikey(ik)), ik=1,nk) ! echo
    
  endfunction ! parse
  
  integer(kind=1) function atomic_number_by_symbol(Sym) result(Z)
  !! retrieves the atomic number from an element symbol of length 3
  !! valid input are all chemical symbols of the periodic table
  !! in correct case, i.e. first letter capitalized
  !!  second letter (if any) lower case
  !! also accepts (integer) numbers in the range [0:121]
  !! errors: Z=-6
    character(len=3), intent(in)    :: Sym

    integer(kind=1), parameter      :: Z_ERROR = -6
    integer                         :: ios
    character                       :: y ! 2nd character of Sym
    integer(kind=1)                 :: t(0:15)
   
    read(unit=Sym, fmt=*, iostat=ios) Z ! try integer reading
    if (ios == 0) then ! an integer number could be read
      if (Z < -1 .or. Z > 116) Z = Z_ERROR
      return ! error
    endif ! ios == 0

    Z = Z_ERROR
    t(0:) = Z_ERROR

    y = Sym(2:2) ! abbrev.

    selectcase(Sym(1:1))
    case('C'); t(1:12)= [  6, 20, 48, 58, 98, 17, 96,112, 27, 24, 55, 29]; Z = t(scan(' adeflmnorsu', y))
    case('P'); t(1:9) = [ 15, 91, 82, 46, 61, 84, 59, 78, 94]; Z = t(scan(' abdmortu', y))
    case('S'); t(1:9) = [ 16, 51, 21, 34,106, 14, 62, 50, 38]; Z = t(scan(' bcegimnr', y))
    case('N'); t(1:8) = [  7, 11, 41, 60, 10, 28,102, 93]; Z = t(scan(' abdeiop', y))
    case('A'); t(1:8) = [ 89, 47, 13, 95, 18, 33, 85, 79]; Z = t(scan('cglmrstu', y))
    case('R'); t(1:8) = [ 88, 37, 75,104,111, 45, 86, 44]; Z = t(scan('abefghnu', y))
    case('T'); t(1:8) = [ 73, 65, 43, 52, 90, 22, 81, 69]; Z = t(scan('abcehilm', y))
    case('B'); t(1:7) = [  5, 56,  4,107, 83, 97, 35]; Z = t(scan(' aehikr', y))
    case('H'); t(1:6) = [  1,  2, 72, 80, 67,108]; Z = t(scan(' efgos', y))
    case('M'); t(1:5) = [101, 12, 25, 42,109]; Z = t(scan('dgnot', y))
    case('F'); t(1:5) = [  9, 26,114,100, 87]; Z = t(scan(' elmr', y))
    case('L'); t(1:5) = [ 57,  3,103, 71,116]; Z = t(scan('airuv', y))
    case('D'); t(1:3) = [105,110, 66]; Z = t(scan('bsy', y))
    case('E'); t(1:3) = [ 68, 99, 63]; Z = t(scan('rsu', y))
    case('G'); t(1:3) = [ 31, 64, 32]; Z = t(scan('ade', y))
    case('I'); t(1:3) = [ 53, 49, 77]; Z = t(scan(' nr', y))
    case('K'); t(1:2) = [ 19, 36]; Z = t(scan(' r', y))
    case('O'); t(1:2) = [  8, 76]; Z = t(scan(' s', y))
    case('Y'); t(1:2) = [ 39, 70]; Z = t(scan(' b', y))
    case('Z'); t(1:2) = [ 30, 40]; Z = t(scan('nr', y))
    case('U'); if (y == ' ') Z =  92 ! Uranium
    case('V'); if (y == ' ') Z =  23 ! Vanadium
    case('W'); if (y == ' ') Z =  74 ! Tungsten (Wolfram)
    case('X'); if (y == 'e') Z =  54 ! Xenon
    ! ==== specialties ====
               if (y == ' ') Z = 115 ! custom element X (jmol style)
    case('_'); if (y == '_') Z =   0 ! "__" vacuum
    case('e'); if (y == ' ') Z =  -1 ! electron
    case('+'); t(1:2) = [113,115]; Z = t(scan('-+', y)) ! custom
    endselect ! S

  endfunction ! atomic_number_by_symbol
  
endmodule PositionReader_mod

#ifdef __MAIN__
program test_PositionsReader_mod
  use PositionReader_mod, only: readXYZfile, PSE, atomic_number_by_symbol
  implicit none

  integer :: ios, na, iZ, jZ
  double precision, allocatable :: apos(:,:)
  character(len=3) :: sym
  
  do iZ = lbound(PSE, 1), ubound(PSE, 1)
    sym = PSE(iZ)
    jZ = atomic_number_by_symbol(sym)
    if (iZ /= jZ) write(0,*) iZ,jZ,sym, '  differ!'
  enddo ! iZ
  
  !! ios = readXYZfile('pos.xyz', na, apos, comm=0) ! interface has changed
  
endprogram
#endif

! #if 0
!   integer function parseOptionalParams(line, params, startatword) result(ist)
!     character(len=*), intent(in) :: line
!     double precision, intent(out) :: params(:)
!     integer, intent(in), optional :: startatword
!     
!     integer, parameter :: ErrorKey=-2
!     integer, parameter :: NumberKey=-1
!     integer, parameter :: IgnoreKey=0
!     integer, parameter :: NK = 4
!      
!     character(len=*), parameter :: keyword(-2:NK) = ['???','<f>','___','vw','rmt','mag','core']
!     integer(kind=1),  parameter :: expects( 1:NK) =                   [  1,   1,    3,    4   ]
!     
!     integer :: ik, ic, iw, is
!     character(len=len(line)) :: copy
!     character(len=16) :: wrd(2*NK)
!     integer           :: key(2*NK)
!     double precision  :: val(2*NK)
!     
!     is = 1; if (present(startatword)) is = max(1, startatword)
!     
!     do ic = 1, len(line)
!       copy(ic:ic) = line(ic:ic)
!       if (line(ic:ic) == '=') copy(ic:ic) = ' ' ! replace '=' by blank
!     enddo ! ic
!     
!     wrd = ''
!     read(unit=copy, fmt=*, iostat=ist) wrd
!     
!     key(:) = IgnoreKey ! init
!     val(:) = 0.d0 ! init
!     
!     do iw = is, 2*NK
!       do ik = 1, NK
!         if (wrd(iw) == '') then
!           key(iw) = IgnoreKey
!         elseif (wrd(iw) == keyword(ik)) then
!           key(iw) = ik
!           write(unit=wrd(iw), fmt='(9a)') trim(keyword(ik)),'=' ! reformat
!         else
!           read(unit=wrd(iw), fmt=*, iostat=ist) val(iw)
!           if (ist /= 0) then
!             key = ErrorKey ! cannot be interpreted
!             wrd(iw) = keyword(ErrorKey)
!           else
!             key = NumberKey
!             write(unit=wrd(iw), fmt='(f16.6)') val(iw) ! reformat
!           endif
!         endif
!       enddo ! ik
!     enddo ! iw
! 
!     if (any(key == ErrorKey)) then
!       write(*, fmt='(99a)') 'string: "',trim(line),'"'
!       write(*, fmt='(99a)') 'detected sequence: ',(trim(keyword(key(iw))),' ',iw=1,2*NK), &
!                             ' where ',trim(keyword(NumberKey)),' is of numeric type'
!       stop 'Error parsing!'
!     endif
!     write(*,'(99(a,1x))') 'string:', (trim(wrd(iw)), iw=1,2*NK)
!     
!     ist = 0
! 
!   endfunction ! parse
! #endif
