module ConfigReader_mod
!-------------------------------------------------------------------------------
!> Summary: Reader for the configuration file
!> Author: Elias Rabel, Paul F Baumeister
!> Category: KKRnano, input-output, initialization
!-------------------------------------------------------------------------------
  use ConfigReaderDictionary_mod, only: Dictionary
  implicit none
  private
  public :: ConfigReader, create, destroy
  public :: parseFile
  public :: getUnreadVariable
  public :: getValue ! (interfaced)
  
  ! Public constants, error codes
  !  Parse errors
  integer, parameter, public :: CONFIG_READER_ERR_NO_ERROR = 0
  integer, parameter, public :: CONFIG_READER_ERR_BAD_CHAR = 1
  integer, parameter, public :: CONFIG_READER_ERR_NO_EQUAL = 2
  integer, parameter, public :: CONFIG_READER_ERR_NO_VALUE = 3
  integer, parameter, public :: CONFIG_READER_ERR_EMPTY_STR = 4
  integer, parameter, public :: CONFIG_READER_ERR_END_LINE = 5

  !  Logical errors
  integer, parameter, public :: CONFIG_READER_ERR_VAR_NOT_UNIQUE = 10
  integer, parameter, public :: CONFIG_READER_ERR_VAR_NOT_FOUND = 11

  ! Datatype errors
  integer, parameter, public :: CONFIG_READER_ERR_NOT_INTEGER = 100
  integer, parameter, public :: CONFIG_READER_ERR_NOT_DOUBLE = 101
  integer, parameter, public :: CONFIG_READER_ERR_NOT_LOGICAL = 102

  ! I/O Errors
  integer, parameter, public :: CONFIG_READER_ERR_NO_FILE = 1000
  integer, parameter, public :: CONFIG_READER_ERR_IO_FAIL = 1001

  
  integer, parameter, public :: CONFIG_READER_USE_DEFAULT_VALUE = -1
  
  type ConfigReader
    private
    type(Dictionary) :: parse_dict
  endtype

  integer, parameter :: MAX_LINE_LENGTH = 150
  integer, parameter :: MAX_FILENAME_LENGTH = 255
  integer, parameter :: FILE_HANDLE = 242
  character(len=*), parameter :: COMMENT_CHARS = '#!'
  character(len=*), parameter :: EQUAL_CHARS = '='
  character(len=*), parameter :: STRING_DELIM = '"' // "'"

  character(len=*), parameter :: ALLOWED_CHARS = &
   'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890+-._'

   
  interface create
    module procedure createConfigReader
  endinterface
  
  interface destroy
    module procedure destroyConfigReader
  endinterface
  
  interface getValue
    module procedure getValueDouble, getValueReal, getValueDoubleVector, getValueInteger, getValueIntVector, getValueLogical, getValueString
  endinterface
  
  contains

!---------------------------------------------------------------------
  subroutine createConfigReader(this)
    use ConfigReaderDictionary_mod, only: create, CONFIG_READER_DICT_VALUE_LENGTH
    type(ConfigReader), intent(inout) :: this

    if (CONFIG_READER_DICT_VALUE_LENGTH < MAX_LINE_LENGTH) then
      write(*,*) "CONFIG_READER: << FATAL ERROR >>"
      write(*,*) "Constant CONFIG_READER_DICT_VALUE_LENGTH has to be >= MAX_LINE_LENGTH"
      stop
    endif

    call create(this%parse_dict) ! createDictionary

  endsubroutine ! create


!---------------------------------------------------------------------
  elemental subroutine destroyConfigReader(this)
    use ConfigReaderDictionary_mod, only: destroy
    type(ConfigReader), intent(inout) :: this

    call destroy(this%parse_dict)
  endsubroutine ! destroy

!---------------------------------------------------------------------
  integer function parseFile(this, filename) result(ierror) ! for parse errors
    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: filename

    integer :: line_number
    integer :: ios ! for I/O errors
    character(len=MAX_LINE_LENGTH) :: line_buf

    ierror = CONFIG_READER_ERR_NO_ERROR
    ios = 0

    line_number = 1

    open(unit=FILE_HANDLE, file=filename, status='old', action='read', iostat=ios)

    if (ios /= 0) then
      write(*,*) "CONFIG_READER: Could not open file ", filename
      ierror = CONFIG_READER_ERR_NO_FILE
      return
    endif

    do while (ios == 0 .and. ierror == CONFIG_READER_ERR_NO_ERROR)
      read(unit=FILE_HANDLE, fmt='(A)', iostat=ios) line_buf

      if (ios > 0) then
        write(*,*) "CONFIG_READER: Error reading file ", filename
        ierror = CONFIG_READER_ERR_IO_FAIL
      endif

      if (ios == 0) then
        !write(*,*) line_buf
        ierror = parseLine(this, line_buf, line_number)
        line_number = line_number + 1
      endif
    enddo ! while

    close(unit=FILE_HANDLE, iostat=ios)

  endfunction ! parseFile

!---------------------------------------------------------------------
  !> parse a line with some simple syntax rules
  !>
  !> e.g. to assign value 5 to the variable VAR, write
  !> VAR = 5
  !> assigning string values
  !> no special characters or spaces - no " " necessary
  !> VAR = yes
  !> for strings with special characters and spaces
  !> VAR = "Hello world!"
  !> the comment characters are allowed in strings
  !> VAR = "#!"
  integer function parseLine(this, line_buf, line_number) result(ierror) ! for parse errors

    use ConfigReaderDictionary_mod, only: pushBackDictionary, &
      CONFIG_READER_DICT_VAR_LENGTH, CONFIG_READER_DICT_VALUE_LENGTH, CONFIG_READER_DICT_NOT_UNIQUE

    type(ConfigReader), intent(inout) :: this
    character(len=MAX_LINE_LENGTH), intent(in) :: line_buf
    integer, intent(in) :: line_number

    character :: parse_char
    character :: string_delim_used
    character(len=CONFIG_READER_DICT_VAR_LENGTH) :: variable
    character(len=CONFIG_READER_DICT_VALUE_LENGTH) :: value

    integer :: dict_error ! for errors such as duplicate variable etc.
    integer :: column
    integer :: variable_counter, value_counter

    integer :: state
    integer :: action

    integer, parameter :: START_MODE = 0
    integer, parameter :: READ_VARIABLE_MODE = 1
    integer, parameter :: SEARCH_EQUAL_MODE = 2
    integer, parameter :: SEARCH_VALUE_MODE = 3
    integer, parameter :: READ_VALUE_MODE = 4
    integer, parameter :: SEARCH_STRING_MODE = 5
    integer, parameter :: READ_STRING_MODE = 6
    integer, parameter :: FINISHED = 7

    integer, parameter :: NO_ACTION = 0
    integer, parameter :: ADD_TO_VARIABLE = 1
    integer, parameter :: ADD_TO_VALUE = 2

    logical :: is_blank
    logical :: is_allowed_char
    logical :: is_equal_char
    logical :: is_comment_char
    logical :: is_string_delim

    ! set error codes to zero
    ierror = 0
    dict_error = 0

    ! initial values
    string_delim_used = ' '
    variable_counter = 0
    value_counter = 0
    state = START_MODE

    do column = 1, len(variable)
      variable(column:column) = ' '
    enddo ! column

    do column = 1, len(value)
      value(column:column) = ' '
    enddo ! column

    char_loop: do column = 1, len(line_buf)

      parse_char = line_buf(column:column)

      is_blank = (parse_char == ' ')
      is_allowed_char = (index(ALLOWED_CHARS, parse_char) /= 0)
      is_equal_char = (index(EQUAL_CHARS, parse_char) /= 0)
      is_comment_char = (index(COMMENT_CHARS, parse_char) /= 0)
      is_string_delim = (index(STRING_DELIM, parse_char) /= 0)

      action = NO_ACTION ! default: do nothing
      ierror = CONFIG_READER_ERR_NO_ERROR

      selectcase (state)

      case (START_MODE)
        if (is_blank) then
          action = NO_ACTION
        elseif (is_allowed_char) then
          state = READ_VARIABLE_MODE
          action = ADD_TO_VARIABLE
        elseif (is_comment_char) then
          state = FINISHED
        else
          ! Unexpected Character found
          state = FINISHED
          ierror = CONFIG_READER_ERR_BAD_CHAR
        endif 

      case (READ_VARIABLE_MODE)
        if (parse_char == ' ') then
          state = SEARCH_EQUAL_MODE
        elseif (is_equal_char) then
          state = SEARCH_VALUE_MODE
        elseif (is_allowed_char) then
          action = ADD_TO_VARIABLE
        elseif (is_comment_char) then
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_VALUE
        else
          state = FINISHED
          ierror = CONFIG_READER_ERR_BAD_CHAR
        endif

      case (SEARCH_EQUAL_MODE)
        if (is_blank) then
          action = NO_ACTION

        elseif (is_equal_char) then
          state = SEARCH_VALUE_MODE

        elseif (is_comment_char) then
          ! no equal sign found
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_EQUAL
        else
          ! expected =
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_EQUAL
        endif

      case (SEARCH_VALUE_MODE)
        if (is_blank) then
          action = NO_ACTION

        elseif (is_allowed_char) then
          state = READ_VALUE_MODE
          action = ADD_TO_VALUE

        elseif (is_string_delim) then
          state = SEARCH_STRING_MODE
          ! remember string delimiter used
          string_delim_used = parse_char

        elseif (is_comment_char) then
        ! expected Value - none found
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_VALUE

        else
        ! unexpected character
          state = FINISHED
          ierror = 1
        endif

       case (READ_VALUE_MODE)
         if (is_blank) then

           action = ADD_TO_VALUE
           ! Uncomment if reading should stop after one value
           !state = FINISHED
           !action = NO_ACTION

         elseif (is_allowed_char) then
           action = ADD_TO_VALUE

         elseif (is_comment_char) then
           state = FINISHED
         else
         ! unexpected character
           state = FINISHED
           ierror = CONFIG_READER_ERR_BAD_CHAR
         endif

       ! check for non-empty string
       case (SEARCH_STRING_MODE)
         if (parse_char == string_delim_used) then
           ! error: string was empty
           state = FINISHED
           ierror = CONFIG_READER_ERR_EMPTY_STR
         else
           state = READ_STRING_MODE
           action = ADD_TO_VALUE
         endif

       case (READ_STRING_MODE)
         if (parse_char == string_delim_used) then
         ! end reading string
           state = FINISHED
         else
           action = ADD_TO_VALUE
         endif

       case default
         write (*,*) "FATAL ERROR in simple_parser. Parser is in unknown state."
         STOP
      endselect

      !write(*,*) state, action

      selectcase (action)

        case (ADD_TO_VARIABLE)

          variable_counter = variable_counter + 1
          variable(variable_counter:variable_counter) = parse_char

        case (ADD_TO_VALUE)

          value_counter = value_counter + 1
          value(value_counter:value_counter) = parse_char

      endselect

      if (state == FINISHED) exit

    enddo char_loop

    ! check the state at end of line

    selectcase (state)
      case (READ_VARIABLE_MODE);  ierror = CONFIG_READER_ERR_NO_VALUE
      case (SEARCH_EQUAL_MODE);   ierror = CONFIG_READER_ERR_NO_VALUE
      case (SEARCH_VALUE_MODE);   ierror = CONFIG_READER_ERR_NO_VALUE
      case (SEARCH_STRING_MODE);  ierror = CONFIG_READER_ERR_END_LINE
      case (READ_STRING_MODE);    ierror = CONFIG_READER_ERR_END_LINE
      case (START_MODE) ! do nothing
      case (READ_VALUE_MODE) ! do nothing
      case default ! do nothing
    endselect

    !if (variable_counter /= 0) write(*,*) variable
    !if (value_counter /= 0) write(*,*) value

    if (variable_counter /= 0 .and. value_counter /= 0) then
      call pushBackDictionary(this%parse_dict, variable, value, .false., dict_error)
      if (dict_error == CONFIG_READER_DICT_NOT_UNIQUE) then
        ierror = CONFIG_READER_ERR_VAR_NOT_UNIQUE
      endif
    endif

    if (ierror /= 0) call displayParserError(line_number, column, ierror)

  endfunction ! parseLine

!---------------------------------------------------------------------
  subroutine displayParserError(line_number, column, ierror)
    integer, intent(in) :: line_number, column, ierror

    if (ierror == 0) return
    write(*,'(9(A,I0))') 'CONFIG_READER: Error in line ', line_number, ' in column ', column
    selectcase (ierror)
      case (CONFIG_READER_ERR_BAD_CHAR);       write (*,*) "Character not allowed."
      case (CONFIG_READER_ERR_NO_EQUAL);       write (*,*) "Expected equal sign not found."
      case (CONFIG_READER_ERR_NO_VALUE);       write (*,*) "No value assigned to variable."
      case (CONFIG_READER_ERR_EMPTY_STR);      write (*,*) "Value is empty string."
      case (CONFIG_READER_ERR_END_LINE);       write (*,*) "Unexpected end of line."
      case (CONFIG_READER_ERR_VAR_NOT_UNIQUE); write (*,*) "Parameter is already defined."
      case default;                            write (*,*) "Unknown error."
    endselect
    
  endsubroutine ! displayParserError

!---------------------------------------------------------------------
  integer function getValueString(this, variable, value, def) result(ierror)
    use ConfigReaderDictionary_mod, only: getDictionaryValue

    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: variable
    character(len=*), intent(inout) :: value
    character(len=*), intent(in), optional :: def

    logical :: tag

    tag = .true.
    ierror = 0
    call getDictionaryValue(this%parse_dict, variable, value, tag, ierror)
    
    if (ierror /= 0) then
      if (present(def)) then
        value = def
        ierror = CONFIG_READER_USE_DEFAULT_VALUE
      else
        ierror = CONFIG_READER_ERR_VAR_NOT_FOUND
      endif
    endif
    
  endfunction ! getValue

!---------------------------------------------------------------------
!> Get value of 'variable' and interpret it as integer value.
!>
!> a default value can be passed as int_value, which does not change on exit
!> if the variable is not found
  integer function getValueInteger(this, variable, value, def) result(ierror)
    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: variable
    integer, intent(inout) :: value
    integer, intent(in), optional :: def

    character(len=MAX_LINE_LENGTH) :: value_string
    integer :: value_read
    integer :: ios

    ierror = 0
    ios = 0

    ierror = getValueString(this, variable, value_string)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) value_read
      ! test if an integer was read
      if (ios == 0) then
        value = value_read
      else
        ierror = CONFIG_READER_ERR_NOT_INTEGER
      endif
    endif

    if (ierror /= 0 .and. present(def)) then
      value = def
      ierror = CONFIG_READER_USE_DEFAULT_VALUE
    endif
    
  endfunction ! getValue

!---------------------------------------------------------------------
!> Get value of 'variable' and interpret it as double prec. value.
!>
!> a default value can be passed as value, which does not change on exit
!> if the variable is not found
  integer function getValueDouble(this, variable, value, def) result(ierror)
    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: variable
    double precision, intent(inout) :: value
    double precision, intent(in), optional :: def

    character(len=MAX_LINE_LENGTH) :: value_string
    double precision :: value_read
    integer :: ios

    ierror = 0
    ios = 0

    ierror = getValueString(this, variable, value_string)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) value_read
      ! test if a double was read
      if (ios == 0) then
        value = value_read
      else
        ierror = CONFIG_READER_ERR_NOT_DOUBLE
      endif
    endif
    
    if (ierror /= 0 .and. present(def)) then
      value = def
      ierror = CONFIG_READER_USE_DEFAULT_VALUE
    endif

  endfunction ! getValue

  !! warpper for single precision
  integer function getValueReal(this, variable, value, def) result(ierror)
    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: variable
    double precision, intent(inout) :: value
    real, intent(in) :: def
    ierror = getValueDouble(this, variable, value, def=dble(def))
  endfunction ! getValue
  
  
!---------------------------------------------------------------------
!> Get value of 'variable' and interpret it as logical value.
!>
! a default value can be passed as value, which does not change on exit
! if the variable is not found
  integer function getValueLogical(this, variable, value, def) result(ierror)
    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: variable
    logical, intent(inout) :: value
    logical, intent(in), optional :: def

    character(len=MAX_LINE_LENGTH) :: value_string
    logical :: value_read
    integer :: ios

    ierror = 0
    ios = 0

    ierror = getValueString(this, variable, value_string)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) value_read
      ! test if a logical was read
      if (ios == 0) then
        value = value_read
      else
        ierror = CONFIG_READER_ERR_NOT_LOGICAL
      endif
    endif

    if (ierror /= 0 .and. present(def)) then
      value = def
      ierror = CONFIG_READER_USE_DEFAULT_VALUE
    endif
    
  endfunction ! getValue

!---------------------------------------------------------------------
!> Reads a double precision vector of fixed length from config-File.
!>
!> The dimension of the vector has to be passed as argument 'length'
!> The routine does not check if more then 'length' values are present
!> If the vector is too short, however, ierror contains an error-code
!
!> A default value can be passed as vector, which does not change on exit
!> if the variable is not found
  integer function getValueDoubleVector(this, variable, value, def) result(ierror)
    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: variable
    double precision, intent(inout) :: value(:)
    double precision, intent(in), optional :: def

    character(len=MAX_LINE_LENGTH) :: value_string
    double precision :: vector_read(size(value))
    integer :: ios

    ierror = 0
    ios = 0

    ierror = getValueString(this, variable, value_string)

    if (ierror == 0) then
      if (present(def)) vector_read(:) = def
      read(unit=value_string, fmt=*, iostat=ios) vector_read(:)
      ! test if read was successful
      if (ios == 0) then
        value = vector_read
      else
        write(*,*) value_string
        write(*,*) vector_read
        value = vector_read
        write(*,*) ios
        ierror = CONFIG_READER_ERR_NOT_DOUBLE
      endif
    endif

    if (ierror /= 0 .and. present(def)) then
      value = def
      ierror = CONFIG_READER_USE_DEFAULT_VALUE
    endif
    
  endfunction ! getValue

!---------------------------------------------------------------------
!> Reads an integer vector of fixed length from config-File.
!>
!> The dimension of the vector has to be passed as argument 'length'
!> The routine does not check if more then 'length' values are present
!> If the vector is too short, however, ierror contains an error-code
!
!> A default value can be passed as vector, which does not change on exit
!> if the variable is not found
  integer function getValueIntVector(this, variable, value, def) result(ierror)
    type(ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: variable
    integer, intent(inout) :: value(:)
    integer, intent(in), optional :: def

    character(len=MAX_LINE_LENGTH) :: value_string
    integer :: vector_read(size(value))
    integer :: ios

    ierror = 0
    ios = 0

    ierror = getValueString(this, variable, value_string)

    if (ierror == 0) then
      if (present(def)) vector_read(:) = def
      read(unit=value_string, fmt=*, iostat=ios) vector_read(:)
      ! test if read was successful
      if (ios == 0) then
        value = vector_read
      else
        write(*,*) value_string
        write(*,*) vector_read
        value = vector_read
        write(*,*) ios
        ierror = CONFIG_READER_ERR_NOT_INTEGER
      endif
    endif

    if (ierror /= 0 .and. present(def)) then
      value = def
      ierror = CONFIG_READER_USE_DEFAULT_VALUE
    endif
    
  endfunction ! getValue

!---------------------------------------------------------------------
!> This subroutine is used to get variables, which have not been read (yet).
!>
!> Initially an integer with value 1 has to be passed as next_ptr
!> next_ptr is an integer which has to be kept for the next call to this
!> subroutine, which gives the next unread variable
!> If no unread variable was found then ierror = CONFIG_READER_ERR_VAR_NOT_FOUND
  integer function getUnreadVariable(this, variable, next_ptr) result(ierror)
    use ConfigReaderDictionary_mod, only: getTaggedVariable, CONFIG_READER_DICT_NOT_FOUND

    type(ConfigReader), intent(in) :: this
    character(len=*),  intent(out) :: variable
    integer, intent(inout) :: next_ptr

    integer :: dict_error

    variable = ' '
    ierror = 0
    dict_error = 0

    call getTaggedVariable(this%parse_dict, variable, .false., next_ptr, dict_error)

    if (dict_error == CONFIG_READER_DICT_NOT_FOUND) ierror = CONFIG_READER_ERR_VAR_NOT_FOUND

  endfunction ! get

endmodule ! ConfigReader_mod

