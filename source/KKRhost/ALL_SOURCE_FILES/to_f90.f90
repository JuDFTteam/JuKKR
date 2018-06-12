! Takes Fortran 77 code in standard format and makes some changes to produce
! free-format Fortran 90 code.
! N.B. It expects STANDARD F77 code.   Non-standard extensions such as
! DO .. END DO (i.e. no label) or in-line comments may cause havoc!

! Changes included are:
! C or c in column 1 replaced with !
! Continuation denoted by a character in column 6 replaced with & at the
! end of the previous line.
! Indenting of code for DO-loops and IF blocks.
! END of program unit replaced by END SUBROUTINE (/PROGRAM/FUNCTION) name
! Fortran `keywords' are in upper case, all other words other than those
! in character strings are converted to lower case.
! .LT., .EQ., etc. replaced with <, ==, etc.
! Labels removed from DO loops; all of which will end with END DO.
! If labels are not referenced, they are removed.
! Short continued lines are adjoined to the previous line.
! ENDIF, ELSEIF & GOTO split into separate words.
! 3-way arithmetic IF constructs are converted to IF .. ELSE IF form.
! Embedded blanks are removed from numbers in DATA statements.
! INTENT declarations are added for dummy arguments.
! Some GO TOs are converted to CYCLE or EXIT.
! Converts CHARACTER * to CHARACTER (LEN=xx) ::.
! Converts computed GO TOs to SELECT CASE.

! To be done:
! DATA statements to be replaced by assignments on the declaration line.
! IMPLICIT NONE statements to be included.
! Declaration of types of unlisted variables.
! Functions to be converted to ELF90 form, i.e. REAL FUNCTION XYZ(arg)
! converted to FUNCTION xyz(arg) RESULT(fn_val).

! Known problems
! Cannot handle character strings or names broken at the end of lines.
! No attempt to convert BLOCKDATA, COMMON or EQUIVALENCE.
! Does not convert Hollerith strings, e.g. 31HTHIS IS A COMMENT ...
! May do the wrong thing if variable names start with IF or end with DO.
! INTENTs are sometimes wrong.  In particular, INTENT(IN) arguments are
! often shown as INTENT(IN OUT).
! Cannot handle comment lines in the middle of continued instructions.
! Can handle 'character*(*) str' but not 'character str*(*)'.

! The default extension for the name of the input file is `for'; this can be
! over-ruled by giving the full name (e.g. myprog.f77).   The output file name
! will be the input name (and directory) with extension `.f90'.

! Added conversion of `enddo' to END DO - 13 March 1997
! Corrected bug which occurred when an arithmetic IF within a DO-loop involved
! a jump to the end of the DO-loop - 17 August 1997.

! ELSEIF, ENDIF & ELSEIF were being split into 2 separate words, and then the
! last letter converted back to lower case - corrected 17 August 1997.
! Corrected bug which occurred when .LT. (or other comparison) had a blank
! before and/or after, followed on the same line by a text string, followed
! by a Fortran word such as THEN or GO TO - 8 December 1997.
! Added (LEN=1) after CHARACTER if length not specified - 9 December 1997.
! Embedded blanks are removed from numerical constants in DATA statements.
! Added 9 December 1997.
! Added INTENTs and TYPE declarations for dummy arguments - 23 December 1997.
! Corrected problem when DO statement contains a comma immediately after DO,
! and improved the detection of INTENTs when a dummy argument appears in an
! IF-expression.  Added extra indentation on continuation lines to improve
! readability - 13 January 1998
! Corrected a bug which could occur when the last type declaration was matched
! to a dummy variable and the line deleted - 5 June 1998
! Corrected jumps out of inner nested DO loops, and replaced GO TOs, out of
! DO loops to the next executable line, with EXIT - 8 June 1998
! Added conversion of CHARACTER * to CHARACTER (LEN=xx) ::
! including CHARACTER*10 a, d, c*50, d   - 21 June 1998.
! Corrected for case of final command of a DO loop which is not CONTINUE and
! which flows onto the next line - 29 June 1998.
! Added conversion of computed GO TOs to SELECT CASE form, and
! fixed problem when a CHARACTER dummy variable had '*' as its last
! dimension - 26 November 1998.
! Fixed problems when the dimensions of a dummy variable involved another
! dummy variable, e.g. wk(nrow*(ncols+1)) - 25 December 1998
! Added date & time stamp - 27 January 1999
! Finally fixed the problems with CYCLE & EXIT, I hope! - 2 February 1999
! Fixed a problem when a type declaration was continued and the next line
! declared the type(s) of dummy arguments - 3 February 1999
! Added conversion of PARAMETER statements from PARAMETER (name1=v1, .. )
! to TYPE1, PARAMETER :: name1=v1  - 8 February 1999
! Added EQV to the list of FORTRAN `words' - 11 February 1999
! Partially corrected problems with the construct:
! IF (condition) GO TO (10, 20, ..
! ..., 99), next
! i.e. with IF & computed GOTO in the same statement (without a THEN), and
! continued onto the next line.
! Also changed a DATA statement to a PARAMETER statement to make the code
! compatible with ELF90 (Thanks to David Ormerod) - 20 May 1999
! Added test for existence of source file.  Program crashed previously if
! the file was not found - 3 August 1999
! Corrected SUBROUTINE fix_3way_IF so that it does not interpret IFIX (or
! similar) as an IF - 23 January 2000.
! At last fixed up strings in quotes which flowed from one line to the next
! - 24 January 2000
! Fixed an error which sometimes caused GOTOs to be wrongly converted to CYCLE
! - 21 March 2000
! Made minor change which helps finding the INTENT of variables which occur
! at the beginning or end of conditions in IF(condition expression)
! - 28 April 2000

! Latest revision - 28 April 2000
! Author - Alan.Miller @ vic.cmis.csiro.au
! WWW-page:  www.ozemail.com.au/~milleraj
! Overflow:  users.bigpond.net.au/amiller/


module implicit
  ! Module to set and reset implicit variable types for use by to_f90.

  implicit none
  integer, save :: var_type(26) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  ! a b c d e f g h i j k l m n o p q r s t u v w x y z
  character (len=24), save :: vt(0:7) = (/ 'NO TYPE                 ', &
    'REAL                    ', 'INTEGER                 ', &
    'real (kind=dp)        ', 'LOGICAL                 ', &
    'COMPLEX                 ', 'CHARACTER               ', &
    'OTHER TYPE              ' /)

contains


  subroutine reset_defaults

    var_type(1:8) = 1              ! REAL (A-H)
    var_type(9:14) = 2             ! INTEGER (I-N)
    var_type(15:26) = 1            ! REAL (O-Z)

    return
  end subroutine reset_defaults



  subroutine set_implicit_types(text)
    ! Read in implicit statement and interpret.

    character (len=*), intent (inout) :: text

    ! Local variables
    integer :: ivt, length, start, i, j, pos, left, right
    logical :: first

    i = index(text, 'IMPLICIT')
    if (i>0) text = text(i+8:)
    text = adjustl(text)

    do
      if (text(1:4)=='NONE') then
        var_type = 0
        return
      else if (text(1:4)=='REAL') then
        ivt = 1
      else if (text(1:7)=='INTEGER') then
        ivt = 2
      else if (text(1:24)=='real (kind=dp) COMPLEX') then
        ivt = 7
        vt(7) = 'real (kind=dp) COMPLEX'
      else if (text(1:16)=='real (kind=dp)') then
        ivt = 3
      else if (text(1:7)=='LOGICAL') then
        ivt = 4
      else if (text(1:7)=='COMPLEX') then
        ivt = 5
      else if (text(1:9)=='CHARACTER') then
        ivt = 6
      else
        ivt = 7
        i = index(text, ' ')
        vt(7) = text(1:i-1)
      end if

      ! Interpret the part in brackets, e.g. (a - h, o - z)

      length = len_trim(text)
      start = 5
      left = index(text(start:length), '(') + start - 1
      if (left<start) return
      right = index(text(start:length), ')') + start - 1
      if (right<left) return
      ! Interpret text(left+1:right-1)
      first = .true.
      do pos = left + 1, right
        select case (text(pos:pos))
        case (' ')
          cycle
        case ('-')
          first = .false.
        case (',', ')')
          if (first) then
            var_type(i) = ivt
          else
            var_type(i:j) = ivt
            first = .true.
          end if
        case default
          if (first) then
            i = ichar(text(pos:pos)) - ichar('a') + 1
            if (i<1) then
              i = ichar(text(pos:pos)) - ichar('A') + 1
            end if
          else
            j = ichar(text(pos:pos)) - ichar('a') + 1
            if (j<1) then
              j = ichar(text(pos:pos)) - ichar('A') + 1
            end if
          end if
        end select
      end do

      start = right + 1
      if (start>=length) return
      text = text(start:length)
      do
        if (text(1:1)==',' .or. text(1:1)==' ') then
          text = text(2:)
        else
          exit
        end if
      end do
    end do

    return
  end subroutine set_implicit_types



  function implicit_type(ch) result (vtype)
    ! Return the variable type given the first character of its name.
    ! The first character is expected to be lower case, but just in case ..

    character (len=1), intent (in) :: ch
    character (len=24) :: vtype

    ! Local variable
    integer :: i, j

    i = ichar(ch) - ichar('a') + 1
    if (i>=1 .and. i<=26) then
      j = var_type(i)
      vtype = vt(j)
    else
      i = ichar(ch) - ichar('A') + 1
      if (i>=1 .and. i<=26) then
        j = var_type(i)
        vtype = vt(j)
      else
        vtype = ' '
      end if
    end if

    return
  end function implicit_type

end module implicit



program to_f90
  use :: implicit
  implicit none

  type :: code
    character (len=140) :: text
    character (len=5) :: label
    type (code), pointer :: next
  end type code

  type :: argument
    character (len=10) :: name
    integer :: intention           ! IN = 1, OUT = 2, IN OUT = 3
    character (len=24) :: var_type                             ! Room for real (kind=dp) COMPLEX


    integer :: dim                 ! DIM = 0 for scalars
    character (len=24) :: dimensions                             ! Not used if DIM = 0


    type (argument), pointer :: next
  end type argument

  character (len=60) :: f77_name, f90_name
  character (len=1) :: tab = char(9), ch
  character (len=50) :: prog_unit_name = ' ', blank = ' ', case_expr
  character (len=9) :: delimiters = ' =+-*/,()'
  character (len=10) :: numbers = '1234567890'
  character (len=5) :: lab
  character (len=30) :: text, vtype
  character (len=140) :: statement
  character (len=8) :: date
  character (len=10) :: time
  integer :: iostatus, pos, count, last, n_marks, pos1(20), pos2(20), &
    lab_length, indent, i, i1, i2, length, numb_arg, i3, i4
  type (code), pointer :: head, current, tail, last_line, next_line, &
    first_decl, last_decl, start_prog_unit, end_prog_unit
  logical :: asterisk, ok, data_stmnt, first_arg, continuation
  type (argument), pointer :: arg_start, arg, last_arg

  interface
    subroutine mark_text(text, n_marks, pos1, pos2, continuation)
      implicit none
      character (len=*), intent (in) :: text
      integer, intent (out) :: n_marks, pos1(:), pos2(:)
      logical, intent (in) :: continuation
    end subroutine mark_text

    subroutine convert_text(text, n_marks, pos1, pos2)
      implicit none
      character (len=*), intent (inout) :: text
      integer, intent (in) :: n_marks
      integer, intent (inout) :: pos1(:), pos2(:)
    end subroutine convert_text

    subroutine remove_data_blanks(text)
      implicit none
      character (len=*), intent (inout) :: text
    end subroutine remove_data_blanks

    function last_char(text) result (ch)
      implicit none
      character (len=*), intent (in) :: text
      character (len=1) :: ch
    end function last_char

    function find_delimited_name(text, name) result (pos)
      implicit none
      character (len=*), intent (in) :: text, name
      integer :: pos
    end function find_delimited_name
  end interface

  open (7777, file='tof90_filename.txt', form='formatted')
  do
    write (*, '(a)', advance='NO') ' Enter name of Fortran source file: '
    write (*, *) 'read filename from tof90_filename.txt'
    read (7777, '(a)', iostat=iostatus) f77_name
    ! READ(*, '(a)', IOSTAT=iostatus) f77_name
    if (iostatus<0) stop           ! Halts gracefully when the names are
    ! read from a file and the end is reached

    if (len_trim(f77_name)==0) cycle
    if (index(f77_name,'.')==0) then
      last = len_trim(f77_name)
      f77_name(last+1:last+4) = '.for'
    end if
    open (8, file=f77_name, status='old', iostat=iostatus)
    if (iostatus/=0) then
      write (*, *) '** Unable to open file: ', f77_name
      cycle
    end if

    pos = index(f77_name, '.', back=.true.) ! Added BACK=.TRUE. for Unix
    ! names e.g. prog.test.f
    f90_name = f77_name(1:pos) // 'f90'
    open (9, file=f90_name)

    ! Set up a linked list containing the lines of code

    nullify (head, tail)
    allocate (head)
    tail => head
    read (8, '(a)') head%text
    if (head%text(1:1)=='C' .or. head%text(1:1)=='c' .or. head%text(1:1)=='*') &
      then
      head%text(1:1) = '!'
    else if (head%text(1:1)==tab) then
      head%text = '      ' // head%text(2:)
    end if
    head%label = ' '
    count = 1

    do
      nullify (current)
      allocate (current)
      read (8, '(a)', iostat=iostatus) current%text
      if (iostatus/=0) exit

      ! Change C, c or * in column 1 to !
      if (current%text(1:1)=='C' .or. current%text(1:1)=='c' .or. &
        current%text(1:1)=='*') then
        if (len_trim(current%text)>1) then
          current%text(1:1) = '!'
        else
          current%text = ' '       ! Leave blank if nothing else on line
        end if
        current%label = ' '
      else
        current%label = adjustl(current%text(1:5))
      end if

      count = count + 1
      if (current%label(1:1)==tab) then ! Expand tabs
        current%label = ' '
        current%text = '      ' // current%text(2:)
      else if (current%label(1:1)=='!') then
        current%label = ' '
      else
        current%label = adjustl(current%label)
      end if

      nullify (current%next)
      tail%next => current
      tail => current
    end do

    write (*, *) 'No. of lines read =', count

    ! ---------------------------------------------------------------------------

    current => head
    nullify (last_line)
    data_stmnt = .false.

    do
      ! Look for blanks in columns 1-5 followed by non-blank in column 6.
      ! If found, add an ampersand at the end of the previous line.

      if (current%label=='     ' .and. current%text(6:6)/=' ' .and. &
        current%text(1:1)/='!') then
        last = len_trim(last_line%text)
        last_line%text(last+3:last+3) = '&'
        current%text(6:6) = ' '
        continuation = .true.
      else
        data_stmnt = .false.
        continuation = .false.
      end if

      ! Replace tabs with single spaces
      do
        pos = index(current%text, tab)
        if (pos==0) exit
        current%text(pos:pos) = ' '
      end do

      ! Remove leading blanks
      current%text = adjustl(current%text)

      ! Mark regions of text which must not have their case changed.
      call mark_text(current%text, n_marks, pos1, pos2, continuation)

      ! Convert cases of regions which are not protected.
      call convert_text(current%text, n_marks, pos1, pos2)

      ! If line is start of a program unit, record its name
      if (current%text(1:7)=='PROGRAM') then
        prog_unit_name = current%text(1:50)
      else if (current%text(1:10)=='SUBROUTINE') then
        pos = index(current%text, '(') - 1
        if (pos<0) pos = len_trim(current%text)
        prog_unit_name = current%text(1:pos)
      else if (current%text(1:9)=='BLOCKDATA') then
        prog_unit_name = current%text(1:50)
      else
        ! N.B. 'FUNCTION' could be part of a comment
        pos = index(current%text, 'FUNCTION')
        if (pos>0 .and. index(current%text,'!')==0 .and. &
          index(current%text,'''')==0) then
          last = index(current%text, '(') - 1
          if (last<0) last = len_trim(current%text)
          prog_unit_name = current%text(pos:last)
        end if
      end if

      ! If first word is one of INTEGER, REAL, real (kind=dp), CHARACTER ,
      ! LOGICAL or COMPLEX, add :: unless FUNCTION appears on the same line
      ! or next non-blank character is '*' as in REAL*8.
      if (index(current%text,'FUNCTION')==0) then
        pos = 0
        if (index(current%text,'INTEGER')==1) then
          pos = 9
        else if (index(current%text,'REAL')==1) then
          pos = 6
        else if (index(current%text,'real (kind=dp)')==1) then
          pos = 18
        else if (index(current%text,'CHARACTER')==1) then
          pos = 11
        else if (index(current%text,'COMPLEX')==1) then
          pos = 9
        else if (index(current%text,'LOGICAL')==1) then
          pos = 9
        end if

        if (pos>0) then
          asterisk = index(current%text(pos-1:pos), '*') > 0
          if (.not. asterisk) then
            if (pos/=11) then
              current%text = current%text(1:pos-1) // ':: ' // &
                adjustl(current%text(pos:))
            else                   ! CHARACTER type, default length = 1
              current%text = 'CHARACTER (LEN=1) :: ' // &
                adjustl(current%text(pos:))
            end if
          else
            if (pos==11) then      ! CHARACTER * found
              i1 = index(current%text, '*') + 1
              length = len_trim(current%text)
              ! Get length, could be (*)
              do
                if (current%text(i1:i1)/=' ') exit
                if (i1>=length) exit
                i1 = i1 + 1
              end do
              if (current%text(i1:i1)=='(') then
                i1 = i1 + 1
                i2 = index(current%text, ')') - 1
              else
                i2 = index(current%text(i1:), ' ') + i1 - 2
              end if
              current%text = 'CHARACTER (LEN=' // current%text(i1:i2) // &
                ') :: ' // adjustl(current%text(i2+2:))
            end if
          end if
          ! Check for 2 or more lengths in CHARACTER declaration.
          ! e.g. CHARACTER a, b, c*10, d
          ! Put 2nd (& later) declarations on separate lines:
          ! CHARACTER*10 c
          ! But check for CHARACTER*10 a(*) where last * is not a
          ! length but a dimension
          if (pos==11) then
            pos = index(current%text, '::') + 2
            do
              i = index(current%text(pos:), '*')
              if (i==0) exit
              i = i + pos - 1
              length = len_trim(current%text)
              i1 = index(current%text(:i-1), ',', back=.true.)
              i1 = max(pos, i1)
              i2 = index(current%text(i+1:), ',')
              if (i2==0) then
                i2 = length + 1
              else
                i2 = i2 + i
              end if
              ! i1, i2 mark commas at beginning & end of `, name*xx,'
              ! but we could have `name(xx, *), '
              ! Test for * after ( , or ) before ,
              i3 = index(current%text(i1:i2), '(')
              i4 = index(current%text(i1:), ')')
              if (i3>0) then
                i4 = i4 + i1 - 1
                i2 = index(current%text(i4+1:), ',')
                if (i2==0) then
                  i2 = length + 1
                else
                  i2 = i2 + i4
                end if
                pos = i2 + 1
                cycle
              else if (i4>0) then
                i4 = i4 + i1 - 1
                i2 = index(current%text(i4+1:), ',')
                if (i2==0) then
                  i2 = length + 1
                else
                  i2 = i2 + i4
                end if
                pos = i2 + 1
                cycle
              end if

              if (i1==pos .and. i2==length+1) then
                ! Only one declaration left on line, e.g.
                ! CHARACTER :: name*50
                current%text = 'CHARACTER (LEN=' // current%text(i+1:length) &
                  // ') :: ' // adjustl(current%text(i1:i-1))
                exit
              end if

              allocate (next_line)
              next_line%next => current%next
              current%next => next_line
              next_line%text = 'CHARACTER' // current%text(i:i2-1) // ' ' // &
                current%text(i1+1:i-1)
              if (i2<length) then
                current%text = current%text(:i1) // current%text(i2+1:length)
              else
                current%text = current%text(:i1-1)
              end if
            end do
          end if
        end if
      end if

      ! If this is in a DATA statement, eliminate any blanks within numbers
      if (data_stmnt .or. current%text(1:4)=='DATA') then
        call remove_data_blanks(current%text)
        last = len_trim(current%text)
        data_stmnt = .true.
      end if

      ! If line only contains 'END', add the program unit name
      if (len_trim(current%text)==3 .and. current%text(1:3)=='END') then
        current%text = current%text(1:3) // ' ' // prog_unit_name
        prog_unit_name = ' '

        ! Convert `enddo' to 'END DO'
      else if (current%text(1:5)=='enddo') then
        current%text = 'END DO' // current%text(6:)
      end if

      last_line => current
      if (associated(current,tail)) exit
      if (.not. associated(current)) exit
      current => current%next
    end do

    ! -------------------------------------------------------------------------

    ! Now convert Do-loops

    current => head
    write (*, *) '      Converting DO-loops, 3-way IFs, & computed GO TOs'
    do
      if (current%text(1:1)/='!' .and. current%text(1:1)/=' ') then
        pos = index(current%text, 'DO')
        if (pos>0 .and. (current%text(pos+2:pos+2)==' ' .or. current%text(pos+ &
          2:pos+2)==',')) then
          if (current%text(pos+2:pos+2)==',') current%text(pos+2:pos+2) = ' '
          if (pos==1) then
            ok = .true.
          else if (scan(current%text(pos-1:pos-1),delimiters)>0) then
            ok = index(current%text(:pos-1), 'END ') == 0
          else
            ok = .false.
          end if
          if (ok) then
            text = adjustl(current%text(pos+3:))
            last = index(text, ' ')
            lab = text(:last-1)
            if (scan(lab(1:1),numbers)==0) lab = ' '
            lab_length = len_trim(lab)
            if (lab_length>0) then
              pos = index(lab, ',') ! Check for a comma after label
              if (pos>0) then
                lab(pos:) = ' '
                i = index(current%text, ',')
                current%text(i:i) = ' '
                lab_length = pos - 1
              end if
              call do_loop_fixup(current, lab)
            end if
          end if

          ! Test for computed GO TO
        else if (index(current%text,'GO TO')>0) then
          i1 = index(current%text, 'GO TO')
          statement = adjustl(current%text(i1+5:))
          ! Test for a `('
          if (statement(1:1)=='(') then
            ok = .true.
            ! If current line is continued, try appending
            ! the next line
            if (last_char(statement)=='&') then
              next_line => current%next
              length = len_trim(statement) + len_trim(next_line%text)
              ok = (length<=141) .and. (last_char(next_line%text)/='&')
              if (ok) then
                pos = len_trim(statement)
                statement = trim(statement(:pos-1)) // trim(next_line%text)
                current%next => next_line%next
                deallocate (next_line)
              end if
            end if

            if (ok) then
              ! Check for comma between ( and )
              pos = index(statement, ')')
              if (index(statement(2:pos-1),',')>0) then
                ! We could have something like:
                ! IF (condition) GO TO (100, 200, 300) ivar
                ! Before doing any more, split into 3 lines:
                ! IF (condition) THEN
                ! GO TO (100, 200, 300) ivar
                ! END IF
                if (current%text(1:2)=='IF') then
                  if (current%text(3:3)==' ' .or. current%text(3:3)=='(') then
                    current%text = current%text(:i1-1) // 'THEN'
                    i1 = 2
                    call insert_and_moveto_newline(current)
                    current%text = ' '
                    next_line => current
                    call insert_and_moveto_newline(next_line)
                    next_line%text = 'END IF'
                  end if
                end if
                ! Get the CASE variable or expression
                case_expr = adjustl(statement(pos+1:))
                if (case_expr(1:1)==',') case_expr = adjustl(case_expr(2:))
                current%text = current%text(:i1-1) // 'SELECT CASE ( ' // &
                  trim(case_expr) // ' )'
                ! Put in pairs of lines:  CASE ( i )
                ! GO TO i-th label
                call goto_cases(statement(2:pos-1))
              end if
            end if
          end if

          ! Look for IF, then a number as last non-blank character
        else
          pos = index(current%text, 'IF')
          if (pos>0) then
            last = len_trim(current%text)
            if (scan(current%text(last:last),numbers)>0) then
              call fix_3way_if(current)
            end if
          end if
        end if
      end if

      if (associated(current,tail)) exit
      if (.not. associated(current)) exit
      current => current%next
    end do

    ! -------------------------------------------------------------------------

    ! Determine INTENTs for dummy arguments

    write (*, *) '      Determining INTENTs of dummy arguments'

    ! Search for either FUNCTION or SUBROUTINE.
    ! Extract name of program unit.

    current => head
    nullify (last_line)
outer_loop: do
      do
        if (current%text(1:1)/='!' .and. current%text(1:1)/=' ') then
          if (current%text(1:10)=='SUBROUTINE') then
            pos = index(current%text, '(') - 1
            if (pos<0) pos = len_trim(current%text)
            prog_unit_name = current%text(1:pos)
            exit
          else
            pos = index(current%text, 'FUNCTION')
            if (pos>0) then
              last = index(current%text, '(') - 1
              if (last<0) last = len_trim(current%text)
              prog_unit_name = current%text(pos:last)
              exit
            end if
          end if
        end if

        last_line => current
        current => current%next
        if (associated(current,tail)) exit outer_loop
      end do

      ! If there is no blank line between this program unit and the previous
      ! one, then insert one.

      if (associated(last_line)) then
        if (len_trim(last_line%text)>0) then
          call insert_and_moveto_newline(last_line)
          last_line%text = ' '
        end if
      end if

      allocate (start_prog_unit)
      start_prog_unit => current

      ! Find end of program unit

      do
        current => current%next
        if (current%text(1:1)/='!' .and. current%text(1:1)/=' ') then
          if (current%text(1:3)=='END') then
            if (index(current%text(5:),prog_unit_name)>0) then
              allocate (end_prog_unit)
              end_prog_unit => current
              exit
            end if
          end if
        end if
        if (associated(current,tail)) exit outer_loop
      end do

      ! Find first & last declarations

      allocate (first_decl, last_decl)
      call find_declarations(start_prog_unit, end_prog_unit, first_decl, &
        last_decl)
      if (.not. associated(last_decl)) go to 100

      ! Extract list of dummy arguments

      call get_arg_list
      if (numb_arg==0) go to 100

      ! See if the declarations contain any IMPLICIT statements

      call reset_defaults
      current => first_decl
      do
        if (current%text(1:8)=='IMPLICIT') then
          statement = current%text(10:)
          call set_implicit_types(statement)
        end if
        if (associated(current,last_decl)) exit
        current => current%next
      end do

      ! Search through the declarations for variable types & dimensions

      call get_var_types

      ! Search through rest of code to try to determine the INTENTs

      call get_intents

      ! Insert INTENT statements

      statement = first_decl%text
      first_decl%text = ' '
      current => first_decl
      arg => arg_start
      do
        call insert_and_moveto_newline(current)
        current%text = arg%var_type
        select case (arg%intention)
        case (0, 3)
          current%text = trim(current%text) // ', INTENT(IN OUT)'
        case (1)
          current%text = trim(current%text) // ', INTENT(IN)'
        case (2)
          current%text = trim(current%text) // ', INTENT(OUT)'
        end select
        current%text = current%text(:41) // ':: ' // arg%name
        if (arg%dim>0) current%text = trim(current%text) // arg%dimensions

        if (associated(arg,last_arg)) exit
        arg => arg%next
      end do
      call insert_and_moveto_newline(current)
      current%text = statement

      ! Search for, and convert, any PARAMETER statements

      current => first_decl
      do
        if (current%text(1:9)=='PARAMETER') then
          call convert_parameter(current)
        end if
        if (associated(current,last_decl)) exit
        current => current%next
      end do

      ! Insert a blank line after the last declaration if there is not one
      ! there already, or a comment.

      next_line => last_decl%next
      if (next_line%text(1:1)/=' ' .and. next_line%text(1:1)/='!') then
        call insert_and_moveto_newline(last_decl)
        last_decl%text = ' '
      end if

      ! Move onto the next SUBROUTINE or FUNCTION

100   current => end_prog_unit
      if (associated(current,tail)) exit
      last_line => current
      current => current%next
      if (associated(current,tail)) exit
    end do outer_loop

    ! -------------------------------------------------------------------------

    ! Indenting and writing output file

    ! Output header line & any continuation lines

    current => head
    continuation = .false.
    do
      if (continuation) then
        write (9, '(t9, a)') trim(current%text)
      else
        write (9, '(a)') trim(current%text)
      end if
      ch = last_char(current%text)
      current => current%next
      if (ch/='&') exit
      continuation = .true.
    end do
    ! Date & time stamp
    call date_and_time(date, time)
    if (ch/=' ') write (9, *)
    write (9, '("! Code converted using TO_F90 by Alan Miller")')
    write (9, '("! Date: ", a4, "-", a2, "-", a2, "  Time: ", a2, &
      &":", a2,      ":", a2)') date(1:4), date(5:6), date(7:8), time(1:2), &
      time(3:4), time(5:6)
    if (len_trim(current%text)>0) write (9, *)

    indent = 0
    continuation = .false.
    write (*, *) '      Writing file: ', f90_name

    do
      if (current%text(1:1)/='!') then
        if (index(current%text,'END ')>0) then
          if (index(current%text,'END SELECT')==0) indent = max(indent-2, 0)
          write (9, '(a)') blank(:indent) // trim(current%text)
          continuation = (last_char(current%text)=='&')
        else if (index(current%text,'DO ')>0) then
          write (9, '(a)') blank(:indent) // trim(current%text)
          continuation = (last_char(current%text)=='&')
          indent = indent + 2
          ! Temporary reduction in
          ! indentation for `ELSE'
        else if (index(current%text,'ELSE')>0) then
          last = max(0, indent-2)
          write (9, '(a)') blank(:last) // trim(current%text)
          continuation = (last_char(current%text)=='&')
          ! Indent increased if `IF'
          ! is followed by `THEN'
        else if (index(current%text,'IF ')>0 .or. index(current%text,'IF(')>0) &
            then
          current%text = blank(:indent) // trim(current%text)
          ! If IF statement runs onto
          ! next line, try joining
          last = len_trim(current%text)
          if (current%text(last:last)=='&') then
            next_line => current%next
            if (last+len_trim(next_line%text)<80) then
              current%text(last:last) = ' '
              current%text = trim(current%text) // ' ' // trim(next_line%text)
              current%next => next_line%next
            end if
          end if

          write (9, '(a)') trim(current%text)
          continuation = (last_char(current%text)=='&')
          next_line => current
          do
            if (index(next_line%text,' THEN')>0 .or. &
              index(next_line%text,')THEN')>0) then
              indent = indent + 2
              exit
            else
              if (last_char(next_line%text)/='&') exit
            end if
            next_line => next_line%next
          end do
        else

          ! If line ends with '&', attempt to join on the next line if it is
          ! short.

          last = len_trim(current%text)
          if (last>0) then
            if (current%text(last:last)=='&') then
              last = len_trim(current%text(:last-1))
              next_line => current%next
              if (last+indent+len_trim(next_line%text)<78) then
                current%text = current%text(:last) // ' ' // &
                  trim(next_line%text)
                current%next => next_line%next
                deallocate (next_line)
              end if
            end if
          end if

          if (continuation) then
            write (9, '(a)') blank(:indent+4) // trim(current%text)
          else
            write (9, '(a)') blank(:indent) // trim(current%text)
          end if
          continuation = (last_char(current%text)=='&')
        end if
        ! Comment line (unchanged)
      else
        write (9, '(a)') trim(current%text)
        continuation = .false.
      end if
      if (associated(current,tail)) exit
      if (.not. associated(current)) exit
      current => current%next
    end do

    close (8)
    close (9)
  end do

  stop


contains


  subroutine do_loop_fixup(start, lab)

    ! Convert DO-loops from:    DO xxx i=1,n    To:   DO i=1,n
    ! xxx CONTINUE              END DO

    ! `start' points to the first line of the DO loop
    ! `lab' is the label

    type (code), pointer :: start
    character (len=*), intent (in) :: lab

    ! Local variables

    type (code), pointer :: current, end_loop
    integer :: i, j, level, nmult, nl_length
    logical :: continued, jump_from_inner, referenced
    character (len=5) :: label(20), next_label, text
    character (len=10) :: loop_name

    ! -------------------------------------------------------------------
    ! PASS 1. Analysis
    ! Find end of loop (end_loop)
    ! Test for multiple loops using same label
    ! Test for jumps to end of this loop from this DO loop (referenced)
    ! or from inner loops (jump_from_inner)
    ! Find if label is on a statement other than CONTINUE
    ! Find if next executable line beyond loop is labelled (for EXIT)

    current => start%next
    nmult = 1
    level = 0
    jump_from_inner = .false.
    referenced = .false.
    do
      if (current%label==lab) then
        continued = (index(current%text,'CONTINUE')>0)
        exit
      end if

      ! Check for nested DO loop or multiple use of current loop

      if (current%text(1:1)=='!' .or. current%text(1:1)==' ') go to 100
      i = index(current%text, 'DO ')
      if (i>0 .and. index(current%text,'END DO')==0) then
        text = adjustl(current%text(i+3:))
        if (scan(text(1:1),numbers)>0) then
          if (text(:lab_length)==lab) then
            nmult = nmult + 1
          else
            level = level + 1
            i = scan(text, ' ,')
            if (i>0) text = text(:i-1)
            label(level) = text
          end if
        end if
      end if

      ! Check for end of nested loop

      if (current%label/='     ' .and. level>0) then
        do
          if (current%label==label(level)) then
            level = level - 1
            if (level<=0) exit
          else
            exit
          end if
        end do
      end if

      ! Test for GO TO current loop label

      i = index(current%text, 'GO TO')
      if (i>0) then
        text = adjustl(current%text(i+5:))
        if (text(:lab_length)==lab) then
          if (level>0) then
            jump_from_inner = .true.
          else
            referenced = .true.
          end if
        end if
      end if

      ! Get next line

100   if (.not. associated(current)) return
      current => current%next
    end do

    end_loop => current

    ! Find label of next executable line.
    ! First advance past any continuation lines after the end of the DO loop.

    next_label = ' '
    do
      if (last_char(current%text)/='&') exit
      if (.not. associated(current)) go to 110
      current => current%next
    end do

    do
      current => current%next
      if (current%text(1:1)/='!') exit
      if (.not. associated(current)) go to 110
    end do
    next_label = current%label
    nl_length = len_trim(next_label)

    ! -------------------------------------------------------------------
    ! PASS 2. Transform beginning & end of loop

110 current => start

    ! Remove label from DO line
    ! There may be a comma after the label, if so, remove it.

    i = index(current%text, lab(:lab_length))
    current%text = current%text(:i-1) // current%text(i+lab_length:)
    length = len_trim(current%text)
    do j = i, length
      if (current%text(j:j)==' ') cycle
      if (current%text(j:j)==',') current%text(j:j) = ' '
      exit
    end do

    ! Jump out of inner loop detected, set up DO construct.

    if (jump_from_inner) then
      loop_name = 'loop' // lab
      current%text = trim(loop_name) // ':  ' // current%text
      current%label = ' '
    end if

    ! Insert `END DO' at end of loop

    current => end_loop
    if (continued) then
      current%text = 'END DO'
      current%label = ' '
    else
      if (.not. referenced) then
        current%label = ' '
        i = index(current%text, lab(:lab_length))
        if (i>0) current%text = adjustl(current%text(i+lab_length:))
      end if
      ! If there are continuation lines, advance to last one
      do
        if (last_char(current%text)=='&') then
          current => current%next
        else
          exit
        end if
      end do
      call insert_and_moveto_newline(current)
      end_loop => current
      current%text = 'END DO'
    end if
    if (jump_from_inner) current%text = trim(current%text) // ' ' // loop_name

    ! Insert multiple CONTINUE's if necessary

    if (nmult>1) then
      call insert_and_moveto_newline(current)
      end_loop => current
      current%text = lab // ' CONTINUE'
      current%label = lab
    end if

    ! -------------------------------------------------------------------
    ! PASS 3. Change GO TOs to CYCLE or EXIT where appropriate

    current => start%next
    if (continued) then
      do
        if (current%text(1:1)=='!' .or. current%text(1:1)==' ') go to 120
        i = index(current%text, 'GO TO')
        if (i>0) then
          text = adjustl(current%text(i+5:))
          if (text(:5)==lab) then
            current%text(i:) = 'CYCLE'
            if (jump_from_inner) current%text = trim(current%text) // ' ' // &
              loop_name
          else if (nl_length>0 .and. text(:nl_length)==next_label) then
            current%text(i:) = 'EXIT'
            if (jump_from_inner) current%text = trim(current%text) // ' ' // &
              loop_name
          end if
        end if

        ! Get next line

120     current => current%next
        if (associated(current,end_loop)) exit
        if (.not. associated(current)) exit
      end do
    end if

    return
  end subroutine do_loop_fixup



  subroutine fix_3way_if(start)
    ! Convert 3-way IFs to IF () THEN .. ELSE IF () THEN .. ELSE

    type (code), pointer :: start

    ! Local variables

    type (code), pointer :: current
    integer :: pos1, count, length, pos2, i, lab1, lab2, lab3, lenq, &
      next_label, lenz
    character (len=1) :: ch
    character (len=128) :: quantity
    character (len=3) :: zero_txt

    current => start
    length = len_trim(current%text)

    ! Find closing bracket to match the opening bracket.
    ! Only cases with the closing bracket on the same line are converted.

    pos1 = index(current%text, 'IF')

    ! Check that next non-blank character after 'IF' is '('.
    i = pos1 + 2
    do
      ch = current%text(i:i)
      if (ch/=' ') exit
      i = i + 1
      if (i>length) return
    end do
    if (ch/='(') return

    pos1 = i
    count = 1
    pos2 = pos1 + 1
    do
      i = scan(current%text(pos2:length), '()')
      if (i==0) return
      pos2 = i + pos2 - 1
      if (current%text(pos2:pos2)=='(') then
        count = count + 1
      else
        count = count - 1
      end if
      if (count==0) exit
      pos2 = pos2 + 1
    end do

    ! See if there are 3 labels after the closing bracket.

    read (current%text(pos2+1:), *, err=100) lab1, lab2, lab3

    ! As it is probably very old code, the first alphabetic character in the
    ! expression should tell us whether the quantity is REAL or INTEGER.

    do i = pos1 + 1, pos2 - 1
      ch = current%text(i:i)
      if (ch>='i' .and. ch<='n') then
        zero_txt = '0'
        lenz = 1
        exit
      else if (ch>='a' .and. ch<='z') then
        zero_txt = '0.0'
        lenz = 3
        exit
      else if (i==pos2-1) then
        return
      end if
    end do

    quantity = current%text(pos1:pos2)
    lenq = len_trim(quantity)

    ! Find the next executable line to see if it is labelled.
    next_label = 0
    do
      if (.not. associated(current)) exit
      current => current%next
      if (current%text(1:1)=='!' .or. len_trim(current%text)==0) cycle
      if (len_trim(current%label)>0) read (current%label, *) next_label
      exit
    end do
    current => start

    if (lab1==lab2) then
      current%text = current%text(:pos2-1) // ' > ' // zero_txt(:lenz) // &
        ') THEN'
      call insert_and_moveto_newline(current)
      current%text = ' '
      write (current%text, '(a, i5)') 'GO TO ', lab3
      if (lab1/=next_label) then
        call insert_and_moveto_newline(current)
        current%text = 'ELSE'
        call insert_and_moveto_newline(current)
        current%text = ' '
        write (current%text, '(a, i5)') 'GO TO ', lab1
      end if
      call insert_and_moveto_newline(current)
      current%text = 'END IF'

    else if (lab2==lab3) then
      current%text = current%text(:pos2-1) // ' < ' // zero_txt(:lenz) // &
        ') THEN'
      call insert_and_moveto_newline(current)
      current%text = ' '
      write (current%text, '(a, i5)') 'GO TO ', lab1
      if (lab2/=next_label) then
        call insert_and_moveto_newline(current)
        current%text = 'ELSE'
        call insert_and_moveto_newline(current)
        current%text = ' '
        write (current%text, '(a, i5)') 'GO TO ', lab2
      end if
      call insert_and_moveto_newline(current)
      current%text = 'END IF'

    else if (lab1==lab3) then
      current%text = current%text(:pos2-1) // ' == ' // zero_txt(:lenz) // &
        ') THEN'
      call insert_and_moveto_newline(current)
      current%text = ' '
      write (current%text, '(a, i5)') 'GO TO ', lab2
      if (lab1/=next_label) then
        call insert_and_moveto_newline(current)
        current%text = 'ELSE'
        call insert_and_moveto_newline(current)
        current%text = ' '
        write (current%text, '(a, i5)') 'GO TO ', lab1
      end if
      call insert_and_moveto_newline(current)
      current%text = 'END IF'

    else
      current%text = current%text(:pos2-1) // ' < ' // zero_txt(:lenz) // &
        ') THEN'
      call insert_and_moveto_newline(current)
      current%text = ' '
      write (current%text, '(a, i5)') 'GO TO ', lab1
      call insert_and_moveto_newline(current)
      current%text = 'ELSE IF ' // quantity(1:lenq-1) // ' == ' // &
        zero_txt(:lenz) // ') THEN'
      call insert_and_moveto_newline(current)
      current%text = ' '
      write (current%text, '(a, i5)') 'GO TO ', lab2
      if (lab3/=next_label) then
        call insert_and_moveto_newline(current)
        current%text = 'ELSE'
        call insert_and_moveto_newline(current)
        current%text = ' '
        write (current%text, '(a, i5)') 'GO TO ', lab3
      end if
      call insert_and_moveto_newline(current)
      current%text = 'END IF'

    end if

100 return
  end subroutine fix_3way_if



  subroutine insert_and_moveto_newline(current)
    ! Insert a new line AFTER the current line, and move `current' to point to
    ! it.

    type (code), pointer :: current

    ! Local variable
    type (code), pointer :: new_line

    allocate (new_line)
    new_line%next => current%next
    current%next => new_line
    current => new_line

    return
  end subroutine insert_and_moveto_newline



  subroutine find_declarations(start, tail, first_decl, last_decl)
    ! Find the first & last declaration lines in a program unit.

    type (code), pointer :: start, tail
    type (code), pointer :: first_decl, last_decl

    ! Local variables
    character (len=9), parameter :: declaration(13) = (/ 'IMPLICIT ', &
      'INTEGER  ', 'REAL     ', 'DOUBLE   ', 'LOGICAL  ', 'COMPLEX  ', &
      'DIMENSION', 'EXTERNAL ', 'DATA     ', 'COMMON   ', 'PARAMETER', &
      'SAVE     ', 'CHARACTER' /)
    type (code), pointer :: current
    integer :: pos, length, i

    nullify (first_decl, last_decl)

    ! Search for first declaration
    current => start%next
search1: do
      if (current%text(1:1)/='!' .and. current%text(1:1)/=' ') then
        pos = scan(current%text(1:13), delimiters)
        if (pos>0) then
          length = min(9, pos-1)
          if (length>=4) then
            do i = 1, 13
              if (current%text(:length)==declaration(i)(:length)) then
                first_decl => current
                exit search1
              end if
            end do
          end if
        end if
      end if

      current => current%next
      if (associated(current,tail)) return
    end do search1

    ! Search for last declaration

    last_decl => first_decl
    do
      if (current%text(1:1)/='!' .and. current%text(1:1)/=' ') then
        pos = index(current%text, '=')
        if (pos>0) then
          if (pos<12) return
          if (current%text(1:9)/='PARAMETER' .and. current%text(1:9)/= &
            'CHARACTER') return
        end if

        if (current%text(1:4)=='CALL') return

        if (current%text(1:2)=='IF') then
          if (current%text(3:3)==' ') return
          if (current%text(3:3)=='(') return
        end if

        if (current%text(1:3)=='DO ') return

        ! Skip continuation lines

        do
          if (last_char(current%text)/='&') exit
          current => current%next
        end do

        last_decl => current
      end if

      current => current%next
      if (associated(current,tail)) return
    end do

    return
  end subroutine find_declarations


  subroutine get_arg_list
    ! Extract list of dummy arguments

    ! Local variables
    integer :: pos, last

    current => start_prog_unit
    numb_arg = 0
    do                             ! Find '(' if there are any arguments
      pos = index(current%text, '(')
      if (pos==0) then
        if (last_char(current%text)/='&') return
        current => current%next
      else
        exit
      end if
    end do
    pos = pos + 1

    nullify (arg_start)
    allocate (arg_start)
    first_arg = .true.
    do                             ! Loop through lines of arguments
      last = scan(current%text(pos:), ',)')
      if (last==0) then
        if (last_char(current%text)/='&') exit
        current => current%next
        pos = 1
      else
        last = last + pos - 1
        nullify (arg)
        allocate (arg)
        if (first_arg) then
          if (len_trim(current%text(pos:last-1))==0) exit
          arg_start => arg
          first_arg = .false.
          nullify (last_arg)
          allocate (last_arg)
        else
          last_arg%next => arg
        end if
        numb_arg = numb_arg + 1
        last_arg => arg

        arg%name = adjustl(current%text(pos:last-1))
        arg%intention = 0
        arg%var_type = ' '
        arg%dim = 0
        pos = last + 1
      end if
    end do

    return
  end subroutine get_arg_list



  subroutine get_var_types
    ! Search thru the declarations for the types of dummy arguments

    current => first_decl
    do
      text = current%text(:30)
      if (text(:4)=='REAL' .or. text(:7)=='INTEGER' .or. &
        text(:6)=='DOUBLE' .or. text(:9)=='CHARACTER' .or. &
        text(:7)=='LOGICAL' .or. text(:7)=='COMPLEX') then
        ! Copy the variable type to vtype
        last = index(text, ' ::') - 1
        if (last<0) then
          last = index(text, '*')
          if (last==0) then
            last = 24
          else
            last = index(text(last+2:), ' ') + last
          end if
          i1 = last + 2
        else
          i1 = last + 4
        end if
        vtype = text(:last)
        call extract_declarations(i1)

      else if (text(:9)=='DIMENSION') then
        i1 = 11
        vtype = ' '
        call extract_declarations(i1)
      end if

      if (associated(current,last_decl)) exit
      current => current%next
    end do

    ! If there are any arguments for which the type has not been determined,
    ! use the implicit types

    arg => arg_start
    do
      if (arg%var_type==' ') arg%var_type = implicit_type(arg%name(1:1))
      if (associated(arg,last_arg)) exit
      arg => arg%next
    end do

    return
  end subroutine get_var_types


  subroutine get_intents
    ! Search thru the body of the current program unit to try to determine
    ! the intents of dummy arguments.
    ! arg % intention = 0 unknown
    ! = 1 IN
    ! = 2 OUT
    ! = 3 IN OUT

    character (len=80) :: last_part
    integer :: j, nbrac

    do
      if (current%text(1:1)/='!' .and. current%text(1:1)/=' ') then
        statement = current%text
        if (statement(1:3)=='IF ' .or. statement(1:3)=='IF(') then
          ! Split line into two parts
          ! IF (condition) | last_part
          i = index(statement, '(')
          length = len_trim(statement)
          nbrac = 1
          do j = i + 1, length - 1
            if (statement(j:j)==')') then
              nbrac = nbrac - 1
              if (nbrac==0) exit
            else if (statement(j:j)=='(') then
              nbrac = nbrac + 1
            end if
          end do
          if (j<length) then
            last_part = statement(j+1:)
            if (adjustl(last_part)=='THEN') last_part = ' '
          else
            last_part = ' '
          end if
          statement = statement(i+1:j-1)
          ! It is assumed that a variable cannot
          ! be altered inside an IF-expression
          arg => arg_start
          do
            i = find_delimited_name(statement, arg%name)
            if (i>0) then
              if (arg%intention==0) arg%intention = 1
            end if
            if (associated(arg,last_arg)) exit
            arg => arg%next
          end do
          statement = last_part
        end if

        pos = index(statement, '=', back=.true.)
        if (pos>0) then
          if (statement(pos-1:pos-1)/='=' .and. statement(pos-1:pos-1)/='/' &
            .and. statement(pos-1:pos-1)/='<' .and. statement(pos-1:pos-1)/= &
            '>') then

            ! Look for each argument name;
            ! is it before or after '='?
            arg => arg_start
            do
              i = find_delimited_name(statement, arg%name)
              if (i>0) then
                if (i<pos) then
                  arg%intention = ior(arg%intention, 2)
                else
                  if (arg%intention==0) arg%intention = 1
                end if
              end if
              if (associated(arg,last_arg)) exit
              arg => arg%next
            end do
          end if
        end if
      end if

      if (associated(current,end_prog_unit)) exit
      current => current%next
    end do

    return
  end subroutine get_intents



  subroutine goto_cases(text)
    ! Inserts pairs:
    ! CASE (i)
    ! GO TO i-th label
    ! Terminated with:
    ! END SELECT

    character (len=*), intent (inout) :: text

    integer :: case_number, pos, i2

    case_number = 1

    do
      pos = index(text, ',')
      if (pos>0) then
        i2 = pos - 1
      else
        i2 = len_trim(text)
      end if
      call insert_and_moveto_newline(current)
      write (current%text, '("  CASE (", i5, ")")') case_number
      call insert_and_moveto_newline(current)
      current%text = '    GO TO ' // trim(text(:i2))
      if (pos==0) exit
      text = text(pos+1:)
      case_number = case_number + 1
    end do

    call insert_and_moveto_newline(current)
    current%text = 'END SELECT'

    return
  end subroutine goto_cases


  subroutine extract_declarations(start_pos)
    ! Take the current line, and any continuations, look for dummy variables,
    ! and remove them, after storing any relevant type & dimension info.

    integer, intent (in) :: start_pos

    ! Local variables

    integer :: i, i1, j, ndim
    character (len=70) :: text

    i1 = start_pos
    do
      i = scan(current%text(i1:), '(,') ! Find next ( or ,
      ndim = 0
      if (i==0) then               ! No comma or ( on this line
        if (last_char(current%text)=='&') then
          current => current%next
          i1 = 1
          cycle
        else
          text = adjustl(current%text(i1:))
          ! Just in case there is an in-line
          pos = index(text, '!')   ! comment (though illegal in F77)
          if (pos>0) text = text(:pos-1)

          if (len_trim(text)==0) return
          pos = len_trim(current%text)
        end if
      else
        pos = i + i1 - 1
        if (current%text(pos:pos)==',') then ! Comma found
          text = current%text(i1:pos-1)
        else                       ! ( found; find matching )
          count = 1
          ndim = 1
          pos = pos + 1
          do
            j = scan(current%text(pos:), '(,)')
            if (j==0) then         ! No bracket or comma
              if (last_char(current%text)=='&') then
                length = len_trim(current%text)
                next_line => current%next
                current%text = trim(current%text(:length-1)) // ' ' // &
                  adjustl(next_line%text)
                if (associated(next_line,last_decl)) last_decl => current
                current%next => next_line%next
                cycle
              else
                return
              end if
            end if

            pos = pos + j - 1
            select case (current%text(pos:pos))
            case ('(')
              count = count + 1
            case (')')
              count = count - 1
              if (count<=0) then
                text = current%text(i1:pos)
                exit
              end if
            case (',')
              ndim = ndim + 1
            end select
            pos = pos + 1
          end do                   ! End matching ) search
        end if
      end if

      ! Variable name isolated, with ndim dimensions
      ! Now see if it matches a dummy argument

      arg => arg_start
      text = adjustl(text)
      if (ndim<=0) then
        length = len_trim(text)
      else
        length = index(text, '(') - 1
      end if
      do
        if (text(:length)==arg%name) then ! Argument matched
          ! Insert variable type
          if (arg%var_type==' ') arg%var_type = vtype
          if (ndim>arg%dim) then
            arg%dim = ndim
            i = index(text, '(')
            arg%dimensions = text(i:)
          end if
          ! Remove variable ( & comma)
          text = adjustl(current%text(pos+1:))
          if (len_trim(text)==0) then
            if (i1>1) then
              current%text(i1-1:) = ' '
            else
              current%text = ' '
            end if
            if (i1==start_pos) current%text = ' '
            return
          else
            if (text(1:1)==',') text = adjustl(text(2:))
            if (text(1:1)=='&') then
              next_line => current%next
              if (i1==start_pos) then
                current%text = current%text(:i1-1) // ' ' // &
                  adjustl(next_line%text)
                if (associated(next_line,last_decl)) last_decl => current
                current%next => next_line%next
              else
                current%text = current%text(:i1-1) // '  &'
                current => next_line
                i1 = 1
              end if
            else
              current%text = current%text(:i1-1) // ' ' // text
            end if
          end if
          exit
        end if

        if (associated(arg,last_arg)) then
          i1 = pos + 1             ! Skip over comma, if present
          exit
        end if
        arg => arg%next
      end do
    end do

    return
  end subroutine extract_declarations



  subroutine convert_parameter(start)

    ! Convert PARAMETER statements from:
    ! PARAMETER (name1 = value1, name2 = value2, ... )
    ! to:
    ! TYPE1, PARAMETER :: name1 = value1
    ! TYPE2, PARAMETER :: name2 = value2

    type (code), pointer :: start

    ! Local variables

    type (code), pointer :: current, next_line
    integer :: count, i, j, length, pos
    character (len=10) :: text
    character (len=30) :: vtype

    current => start

    ! Replace opening ( with ::

    i = index(current%text, '(')
    if (i==0) return
    current%text = trim(current%text(:i-1)) // ' :: ' // &
      adjustl(current%text(i+1:))
    i = index(current%text, '::') + 3
    do
      j = index(current%text(i:), '=')
      if (j==0) then
        if (last_char(current%text)/='&') return
        next_line => current%next
        j = len_trim(current%text)
        current%text = trim(current%text(:j-1)) // next_line%text
        current%next => next_line%next
        j = index(current%text(i:), '=')
        if (j==0) return
      end if
      j = i + j - 1
      text = adjustl(current%text(i:j-1))
      call find_type(text, vtype, first_decl, start)

      current%text = trim(vtype) // ', ' // current%text
      j = j + 2 + len_trim(vtype)

      ! Is there another value set in this statement?
      ! Find end of the expression for the value, which may involve brackets
      ! and commas.

100   length = len_trim(current%text)
      count = 0
      do i = j + 1, length
        select case (current%text(i:i))
        case ('(')
          count = count + 1
        case (')')
          count = count - 1
          if (count<0) then
            ! Remove final ) and return
            current%text = current%text(:i-1)
            return
          end if
        case (',')
          ! If count = 0, there is another declaration
          if (count==0) then
            ! Break line, check for '&' as first character
            text = adjustl(current%text(i+1:))
            if (text(1:1)=='&') then
              current%text = current%text(:i-1)
              current => current%next
              current%text = 'PARAMETER :: ' // adjustl(current%text)
            else
              allocate (next_line)
              next_line%next => current%next
              current%next => next_line
              next_line%text = 'PARAMETER :: ' // adjustl(current%text(i+1:))
              current%text = current%text(:i-1)
              if (associated(current,last_decl)) last_decl => next_line
              current => next_line
              start => start%next
            end if
            exit
          end if
        case ('&')
          ! Expression continued on next line, merge lines
          next_line => current%next
          pos = len_trim(current%text(:i-1))
          current%text = current%text(:pos) // next_line%text
          current%next => next_line%next
          go to 100
        end select
      end do

      if (i>length) exit
      i = 14
    end do

    return
  end subroutine convert_parameter



  subroutine find_type(vname, vtype, first_decl, last_decl)

    ! Find the type of variable 'vname'

    character (len=*), intent (in) :: vname
    character (len=*), intent (out) :: vtype
    type (code), pointer :: first_decl, last_decl

    ! Local variables

    type (code), pointer :: current
    character (len=30) :: text
    integer :: i1, last, length, pos

    current => first_decl
    length = len_trim(vname)
    if (length==0) return
    do
      text = current%text(:30)
      if (text(:4)=='REAL' .or. text(:7)=='INTEGER' .or. &
        text(:6)=='DOUBLE' .or. text(:9)=='CHARACTER' .or. &
        text(:7)=='LOGICAL' .or. text(:7)=='COMPLEX') then
        ! Copy the variable type to vtype
        last = index(text, ' ::') - 1
        if (last<0) then
          last = index(text, '*')
          if (last==0) then
            last = 24
          else
            last = index(text(last+2:), ' ') + last
          end if
          i1 = last + 2
        else
          i1 = last + 4
        end if
        vtype = text(:last)

        ! See if variable is declared on this line (& any continuation)

        do
          pos = find_delimited_name(current%text(i1:), vname(:length))
          if (pos==0) then
            if (last_char(current%text)=='&') then
              current => current%next
              i1 = 1
              cycle
            end if
          end if
          exit
        end do

        ! Variable name found if pos > 0.

        if (pos>0) then            ! Remove variable name
          pos = pos + i1 - 1
          current%text = current%text(:pos-1) // current%text(pos+length:)
          ! Delete line if only TYPE :: remains
          if (last_char(current%text)==':') then
            current%text = ' '
            return
          end if
          ! Remove any following comma
          i = pos
          length = len_trim(current%text)
          do
            if (i>length) then
              return
            else if (current%text(i:i)==',') then
              current%text = current%text(:i-1) // current%text(i+1:)
              return
            else if (current%text(i:i)/=' ') then
              return
            end if
            i = i + 1
          end do
        end if

      end if

      ! If last declaration has been reached, return default type.
      ! Otherwise proceed to next line.

      if (associated(current,last_decl)) then
        vtype = implicit_type(vname(1:1))
        exit
      else
        current => current%next
      end if
    end do

    return
  end subroutine find_type

end program to_f90



subroutine mark_text(text, n_marks, pos1, pos2, continuation)

  ! Look for exclamation marks or quotes to find any text which must be
  ! protected from case changes.
  ! It is assumed that strings are NOT continued from one line to the next.
  implicit none

  character (len=*), intent (in) :: text
  logical, intent (in) :: continuation
  integer, intent (out) :: n_marks, pos1(:), pos2(:)

  ! Local variables
  integer :: mark, start, pos_exclaim, pos_sngl_quote, pos_dbl_quote, pos, &
    endpos
  character (len=1), save :: quote
  logical, save :: protect = .false.

  mark = 1
  start = 1
  if (continuation .and. protect) then
    pos1(mark) = 1
    pos = 0
    go to 110
  end if

  ! Find next opening quote or exclamation mark

100 protect = .false.
  pos_exclaim = index(text(start:80), '!')
  pos_sngl_quote = index(text(start:80), '''')
  pos_dbl_quote = index(text(start:80), '"')
  if (pos_exclaim==0) pos_exclaim = 81
  if (pos_sngl_quote==0) pos_sngl_quote = 81
  if (pos_dbl_quote==0) pos_dbl_quote = 81
  pos1(mark) = min(pos_exclaim, pos_sngl_quote, pos_dbl_quote)

  if (pos1(mark)==81) then         ! No more protected regions
    n_marks = mark - 1
    return
  else if (pos_exclaim==pos1(mark)) then ! Rest of line is a comment
    pos1(mark) = pos1(mark) + start - 1
    pos2(mark) = 80
    n_marks = mark
    return
  end if

  pos = start - 1 + pos1(mark)
  pos1(mark) = pos
  quote = text(pos:pos)

  ! Search for matching quote

110 endpos = index(text(pos+1:), quote)
  if (endpos>0) then
    pos2(mark) = pos + endpos
    start = pos2(mark) + 1
    mark = mark + 1
    go to 100
  end if

  ! No matching end quote - it should be on the next line

  pos2(mark) = 80
  n_marks = mark
  protect = .true.

  return
end subroutine mark_text


subroutine convert_text(text, n_marks, pos1, pos2)

  ! Convert unprotected text to upper case if it is a FORTRAN word,
  ! otherwise convert to lower case.
  implicit none

  character (len=*), intent (inout) :: text
  integer, intent (in) :: n_marks
  integer, intent (inout) :: pos1(:), pos2(:)

  ! Local variables

  integer :: length, inc = ichar('A') - ichar('a'), pos, mark, i, i1, j, j1, &
    j2, ptr
  logical :: matched
  character (len=11) :: fortran_word(186) = (/ 'ABS        ', 'ACCESS     ', &
    'ACOS       ', 'AIMAG      ', 'AINT       ', 'ALOG       ', 'ALOG10     ', &
    'AMAX0      ', 'AMAX1      ', 'AMIN0      ', 'AMIN1      ', 'AMOD       ', &
    'AND        ', 'ANINT      ', 'APPEND     ', 'ASIN       ', 'ASSIGN     ', &
    'ATAN       ', 'ATAN2      ', 'BACKSPACE  ', 'BLANK      ', 'BLOCK      ', &
    'BLOCKDATA  ', 'BLOCKSIZE  ', 'CALL       ', 'CCOS       ', 'CDABS      ', &
    'CDCOS      ', 'CDEXP      ', 'CDLOG      ', 'CDSIN      ', 'CDSQRT     ', &
    'CEXP       ', 'CHAR       ', 'CHARACTER  ', 'CLOG       ', 'CLOSE      ', &
    'CMPLX      ', 'COMMON     ', 'COMPLEX    ', 'CONJG      ', 'CONTINUE   ', &
    'COS        ', 'COSH       ', 'CSIN       ', 'CSQRT      ', 'DABS       ', &
    'DACOS      ', 'DASIN      ', 'DATA       ', 'DATAN      ', 'DATAN2     ', &
    'DBLE       ', 'DCMPLX     ', 'DCONJG     ', 'DCOS       ', 'DCOSH      ', &
    'DELETE     ', 'DEXP       ', 'DIMAG      ', 'DINT       ', 'DIRECT     ', &
    'DLOG       ', 'DLOG10     ', 'DMAX1      ', 'DIMENSION  ', 'DMIN1      ', &
    'DMOD       ', 'DNINT      ', 'DO         ', 'DOUBLE     ', 'DSIGN      ', &
    'DSIN       ', 'DSINH      ', 'DSQRT      ', 'DTAN       ', 'DTANH      ', &
    'ELSE       ', 'ELSEIF     ', 'END        ', 'ENDFILE    ', 'ENDIF      ', &
    'ENTRY      ', 'EQ         ', 'EQUIVALENCE', 'EQV        ', 'ERR        ', &
    'EXIST      ', 'EXIT       ', 'EXP        ', 'EXTERNAL   ', 'FILE       ', &
    'FLOAT      ', 'FMT        ', 'FORM       ', 'FORMAT     ', 'FORMATTED  ', &
    'FUNCTION   ', 'GE         ', 'GOTO       ', 'GO         ', 'GT         ', &
    'IABS       ', 'IAND       ', 'ICHAR      ', 'IDINT      ', 'IDNINT     ', &
    'IEOR       ', 'IF         ', 'IFIX       ', 'IMPLICIT   ', 'INCLUDE    ', &
    'INDEX      ', 'INPUT      ', 'INQUIRE    ', 'INT        ', 'INTEGER    ', &
    'INTRINSIC  ', 'IOSTAT     ', 'ISIGN      ', 'KEEP       ', 'LE         ', &
    'LEN        ', 'LGE        ', 'LGT        ', 'LLE        ', 'LLT        ', &
    'LOG        ', 'LOG10      ', 'LOGICAL    ', 'LT         ', 'MAX        ', &
    'MAX0       ', 'MAX1       ', 'MIN        ', 'MIN0       ', 'MIN1       ', &
    'MOD        ', 'NAME       ', 'NAMELIST   ', 'NAMED      ', 'NE         ', &
    'NEQV       ', 'NEW        ', 'NEXTREC    ', 'NONE       ', 'NOT        ', &
    'NUMBER     ', 'OLD        ', 'OPEN       ', 'OPENED     ', 'OR         ', &
    'PARAMETER  ', 'PAUSE      ', 'POSITION   ', 'PRECISION  ', 'PRINT      ', &
    'PROGRAM    ', 'READ       ', 'REAL       ', 'REC        ', 'RECL       ', &
    'RETURN     ', 'REWIND     ', 'SAVE       ', 'SCRATCH    ', 'SEQUENTIAL ', &
    'SIGN       ', 'SIN        ', 'SINH       ', 'SNGL       ', 'SPACE      ', &
    'SQRT       ', 'STATUS     ', 'STOP       ', 'SUBROUTINE ', 'TAN        ', &
    'TANH       ', 'THEN       ', 'TO         ', 'TYPE       ', 'UNFORMATTED', &
    'UNIT       ', 'UNKNOWN    ', 'WHILE      ', 'WRITE      ' /)
  character (len=4) :: compare(6) = (/ '.LT.', '.LE.', '.EQ.', '.GE.', '.GT.', &
    '.NE.' /)
  character (len=2) :: replacement(6) = (/ '< ', '<=', '==', '>=', '> ', &
    '/=' /)

  ! A   B   C   D   E   F   G    H    I    J    K    L    M    N    O
  ! P    Q    R    S    T    U    V    W    X    Y    Z
  integer, parameter :: indx(27) = (/ 1, 20, 25, 47, 78, 92, 99, 103, 103, 121 &
    , 121, 122, 132, 139, 149, 153, 159, 159, 165, 177, 182, 185, 185, 187, &
    187, 187, 187 /)

  if (pos1(1)==1 .and. pos2(1)==80) return ! Entire line protected

  pos = 1
  mark = 1
  length = len_trim(text)
  do                               ! Convert to upper case
    if (n_marks>=mark .and. pos==pos1(mark)) then
      pos = pos2(mark) + 1
      mark = mark + 1
      if (pos>=length) exit
    end if
    if (text(pos:pos)>='a' .and. text(pos:pos)<='z') text(pos:pos) &
      = char(ichar(text(pos:pos))+inc)
    pos = pos + 1
    if (pos>length) exit
  end do

  ! Search for `words' in text.
  ! Convert to lower case if they are not FORTRAN words.
  i1 = 1
  pos = 1
  mark = 1
  do
    if (pos>length) exit
    if (n_marks>=mark .and. pos>=pos1(mark)) then
      pos = pos2(mark) + 1
      i1 = pos
      mark = mark + 1
      if (pos>=length) exit
    end if

    do
      if ((text(pos:pos)>='A' .and. text(pos:pos)<='Z') .or. (text( &
        pos:pos)>='0' .and. text(pos:pos)<='9') .or. text(pos:pos)=='_') then
        pos = pos + 1
        cycle
      else
        exit
      end if
    end do

    pos = pos - 1
    ! Now i1 & pos = positions of 1st & last characters of current string

    if (pos<i1) then               ! Single non-alphanumeric character
      pos = i1 + 1
      i1 = pos
      cycle
    end if

    ptr = ichar(text(i1:i1)) - ichar('A') + 1
    if (ptr<1 .or. ptr>26) then
      pos = pos + 1
      if (pos>length) exit
      i1 = pos
      cycle
    end if

    matched = .false.
    if (pos>i1) then
      j1 = indx(ptr)
      j2 = indx(ptr+1) - 1
      do j = j1, j2
        if (text(i1:pos)==fortran_word(j)) then
          matched = .true.
          exit
        end if
      end do
    end if

    ! Replace .LT. with <, etc.
    if (matched .and. i1>1) then
      if (text(i1-1:i1-1)=='.') then
        do j = 1, 6
          if (text(i1-1:pos+1)==compare(j)) then
            text(i1-1:pos+1) = ' ' // replacement(j) // ' '
            exit
          end if
        end do
        do                         ! Remove excess blanks
          i1 = max(i1, 3)
          j1 = index(text(i1-2:pos+2), '  ')
          if (j1==0) exit
          j1 = j1 + i1 - 3
          text(j1:) = text(j1+1:)
          pos2(mark) = pos2(mark) - 1 ! Adjust mark positions
          do i = mark + 1, n_marks
            pos1(i) = pos1(i) - 1
            pos2(i) = pos2(i) - 1
          end do
          pos = pos - 1
        end do
      end if
    end if

    ! Output line of text to screen if it contains SUBROUTINE or FUNCTION.
    ! Convert ENDIF to END IF, ELSEIF to ELSE IF, and GOTO to GO TO.
    if (matched) then
      if (text(i1:pos)=='SUBROUTINE' .or. text(i1:pos)=='FUNCTION') then
        write (*, '(1x, a)') text(1:length)
      else if (text(i1:pos)=='ENDIF') then
        text(i1:) = 'END IF' // text(pos+1:)
        pos = pos + 1
      else if (text(i1:pos)=='ELSEIF') then
        text(i1:) = 'ELSE IF' // text(pos+1:)
        pos = pos + 1
      else if (text(i1:pos)=='GOTO') then
        text(i1:) = 'GO TO' // text(pos+1:)
        pos = pos + 1
      end if
    end if

    ! If text is not matched, convert to lower case, if necessary.
    if (.not. matched) then
      do j = i1, pos
        if (text(j:j)>='A' .and. text(j:j)<='Z') text(j:j) = char(ichar(text( &
          j:j))-inc)
      end do
    end if

    pos = pos + 1
    if (pos>length) exit
    i1 = pos
  end do

  return
end subroutine convert_text



subroutine remove_data_blanks(text)
  ! Remove any blanks embedded between numerical digits in DATA statements

  implicit none
  character (len=*), intent (inout) :: text

  ! Local variables
  integer :: length, pos, i1
  character (len=10) :: numbers = '1234567890'

  length = len_trim(text)
  i1 = 2
  do
    pos = index(text(i1:length), ' ')
    if (pos==0) exit
    i1 = i1 + pos - 1
    if (scan(text(i1-1:i1-1),numbers)>0 .and. scan(text(i1+1:i1+ &
      1),numbers)>0) then
      text = text(:i1-1) // text(i1+1:length)
      length = length - 1
    end if
    i1 = i1 + 2
    if (i1>length) exit
  end do

  return
end subroutine remove_data_blanks


function last_char(text) result (ch)
  ! Return the last character on a line
  implicit none

  character (len=*), intent (in) :: text
  character (len=1) :: ch

  ! Local variable
  integer :: last

  last = len_trim(text)
  if (last==0) then
    ch = ' '
  else
    ch = text(last:last)
  end if

  return
end function last_char


function find_delimited_name(text, name) result (pos)
  ! Find a name in a character string with delimiters either side of it,
  ! or after it if it starts at position 1.
  ! An extended version of the intrinsic INDEX.
  ! pos = the position of the first character of name in text (= 0 if not
  ! found).
  ! N.B. When the name is short (e.g. i or n) it could occur as part of some
  ! other name.

  implicit none
  character (len=*), intent (in) :: text, name
  integer :: pos

  ! Local variables
  integer :: i1, ltext, lname

  i1 = 1
  ltext = len_trim(text)
  lname = len_trim(name)
  do
    pos = index(text(i1:ltext), trim(name))
    if (pos==0) return
    pos = pos + i1 - 1
    if (pos>1) then
      if (scan(text(pos-1:pos-1),' <=+-/*,')>0) then
        if (scan(text(pos+lname:pos+lname),' >=(+-/*,')>0) return
      end if
    else
      if (scan(text(pos+lname:pos+lname),' >=(+-/*,')>0) return
    end if
    i1 = pos + lname
    if (i1+lname>ltext) exit
  end do

  pos = 0

  return
end function find_delimited_name
