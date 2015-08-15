!> @author Elias Rabel

module ConfigReader_mod
  use ConfigReaderDictionary_mod, only: Dictionary
  implicit none
  private
  public :: ConfigReader, create, destroy
  public :: createConfigReader, destroyConfigReader ! deprecated
  public :: parseFile
  public :: getUnreadVariable
  public :: getValueDouble, getValueDoubleVector, getValueInteger, getValueIntVector, getValueLogical, getValueString  
  
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

  type ConfigReader
    private
    type (Dictionary) :: parse_dict
  end type ConfigReader

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
   
  contains

!---------------------------------------------------------------------
  subroutine createConfigReader(this)
    use ConfigReaderDictionary_mod, only: createDictionary, CONFIG_READER_DICT_VALUE_LENGTH
    type (ConfigReader), intent(inout) :: this

    if(CONFIG_READER_DICT_VALUE_LENGTH < MAX_LINE_LENGTH) then
      write(*,*) "CONFIG_READER: << FATAL ERROR >>"
      write(*,*) "Constant CONFIG_READER_DICT_VALUE_LENGTH has to be >= MAX_LINE_LENGTH"
      stop
    end if

    call createDictionary(this%parse_dict)

  end subroutine createConfigReader


!---------------------------------------------------------------------
  subroutine destroyConfigReader(this)
    use ConfigReaderDictionary_mod, only: destroyDictionary
    type (ConfigReader), intent(inout) :: this

    call destroyDictionary(this%parse_dict)

  end subroutine destroyConfigReader

!---------------------------------------------------------------------
  subroutine parseFile(this, filename, ierror)

    type (ConfigReader), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierror ! for parse errors

    integer :: line_number
    integer :: ios ! for I/O errors
    character(len = MAX_LINE_LENGTH) :: line_buf

    ierror = CONFIG_READER_ERR_NO_ERROR
    ios = 0

    line_number = 1

    open(unit = FILE_HANDLE, file = filename, status='old', & 
         action='read', iostat = ios)

    if (ios /= 0) then
      write(*,*) "CONFIG_READER: Could not open file ", filename
      ierror = CONFIG_READER_ERR_NO_FILE
    end if

      do while (ios == 0 .and. ierror == CONFIG_READER_ERR_NO_ERROR)
        read(FILE_HANDLE, '(A)', iostat = ios) line_buf

        if (ios > 0) then
          write(*,*) "CONFIG_READER: Error reading file ", filename
          ierror = CONFIG_READER_ERR_IO_FAIL
        end if

        if (ios == 0) then
          !write(*,*) line_buf
          call parseLine(this, line_buf, line_number, ierror)
          line_number = line_number + 1
        end if
      end do


    close(FILE_HANDLE)

  end subroutine parseFile

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
  subroutine parseLine(this, line_buf, line_number, ierror)

    use ConfigReaderDictionary_mod, only: CONFIG_READER_DICT_VAR_LENGTH, &
      CONFIG_READER_DICT_VALUE_LENGTH, CONFIG_READER_DICT_NOT_UNIQUE, &
      pushBackDictionary

    type (ConfigReader), intent(inout) :: this
    character(len = MAX_LINE_LENGTH), intent(in) :: line_buf
    integer, intent(in) :: line_number
    integer, intent(out) :: ierror ! for parse errors

    character :: parse_char
    character :: string_delim_used
    character(len = CONFIG_READER_DICT_VAR_LENGTH) :: variable
    character(len = CONFIG_READER_DICT_VALUE_LENGTH) :: value

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
    end do

    do column = 1, len(value)
      value(column:column) = ' '
    end do

    char_loop: do column = 1, len(line_buf)

      parse_char = line_buf(column:column)

      is_blank = (parse_char == ' ')
      is_allowed_char = (index(ALLOWED_CHARS, parse_char) /= 0)
      is_equal_char = (index(EQUAL_CHARS, parse_char) /= 0)
      is_comment_char = (index(COMMENT_CHARS, parse_char) /= 0)
      is_string_delim = (index(STRING_DELIM, parse_char) /= 0)

      action = NO_ACTION ! default: do nothing
      ierror = CONFIG_READER_ERR_NO_ERROR

      select case (state)

      case(START_MODE)
        if (is_blank) then
          action = NO_ACTION
        else if (is_allowed_char) then
          state = READ_VARIABLE_MODE
          action = ADD_TO_VARIABLE
        else if (is_comment_char) then
          state = FINISHED
        else
          ! Unexpected Character found
          state = FINISHED
          ierror = CONFIG_READER_ERR_BAD_CHAR
        end if 

      case (READ_VARIABLE_MODE)
        if (parse_char == ' ') then
          state = SEARCH_EQUAL_MODE
        else if (is_equal_char) then
          state = SEARCH_VALUE_MODE
        else if (is_allowed_char) then
          action = ADD_TO_VARIABLE
        else if (is_comment_char) then
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_VALUE
        else
          state = FINISHED
          ierror = CONFIG_READER_ERR_BAD_CHAR
        end if

      case (SEARCH_EQUAL_MODE)
        if (is_blank) then
          action = NO_ACTION

        else if (is_equal_char) then
          state = SEARCH_VALUE_MODE

        else if (is_comment_char) then
          ! no equal sign found
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_EQUAL
        else
          ! expected =
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_EQUAL
        end if

      case (SEARCH_VALUE_MODE)
        if (is_blank) then
          action = NO_ACTION

        else if (is_allowed_char) then
          state = READ_VALUE_MODE
          action = ADD_TO_VALUE

        else if (is_string_delim) then
          state = SEARCH_STRING_MODE
          ! remember string delimiter used
          string_delim_used = parse_char

        else if (is_comment_char) then
        ! expected Value - none found
          state = FINISHED
          ierror = CONFIG_READER_ERR_NO_VALUE

        else
        ! unexpected character
          state = FINISHED
          ierror = 1
        end if

       case (READ_VALUE_MODE)
         if (is_blank) then

           action = ADD_TO_VALUE
           ! Uncomment if reading should stop after one value
           !state = FINISHED
           !action = NO_ACTION

         else if (is_allowed_char) then
           action = ADD_TO_VALUE

         else if (is_comment_char) then
           state = FINISHED
         else
         ! unexpected character
           state = FINISHED
           ierror = CONFIG_READER_ERR_BAD_CHAR
         end if

       ! check for non-empty string
       case (SEARCH_STRING_MODE)
         if (parse_char == string_delim_used) then
           ! error: string was empty
           state = FINISHED
           ierror = CONFIG_READER_ERR_EMPTY_STR
         else
           state = READ_STRING_MODE
           action = ADD_TO_VALUE
         end if

       case (READ_STRING_MODE)
         if (parse_char == string_delim_used) then
         ! end reading string
           state = FINISHED
         else
           action = ADD_TO_VALUE
         end if

       case default
         write (*,*) "FATAL ERROR in simple_parser. Parser is in unknown state."
         STOP
      end select

      !write(*,*) state, action

      select case (action)

        case (ADD_TO_VARIABLE)

          variable_counter = variable_counter + 1
          variable(variable_counter:variable_counter) = parse_char

        case (ADD_TO_VALUE)

          value_counter = value_counter + 1
          value(value_counter:value_counter) = parse_char

      end select

      if (state == FINISHED) exit

    end do char_loop

    ! check the state at end of line

    select case (state)

      case (START_MODE)
        ! do nothing

      case (READ_VALUE_MODE)
        ! do nothing

      case (READ_VARIABLE_MODE)
        ierror = CONFIG_READER_ERR_NO_VALUE

      case (SEARCH_EQUAL_MODE)
        ierror = CONFIG_READER_ERR_NO_VALUE

      case (SEARCH_VALUE_MODE)
        ierror = CONFIG_READER_ERR_NO_VALUE

      case (SEARCH_STRING_MODE)
        ierror = CONFIG_READER_ERR_END_LINE

      case (READ_STRING_MODE)
        ierror = CONFIG_READER_ERR_END_LINE

    end select

    !if (variable_counter /= 0) write(*,*) variable
    !if (value_counter /= 0) write(*,*) value

    if (variable_counter /= 0 .and. value_counter /= 0) then
      call pushBackDictionary(this%parse_dict, variable, value, .false., dict_error)
      if (dict_error == CONFIG_READER_DICT_NOT_UNIQUE) then
        ierror = CONFIG_READER_ERR_VAR_NOT_UNIQUE
      end if
    end if

    if (ierror /= 0) then
      call displayParserError(line_number, column, ierror)
    end if

  end subroutine parseLine

!---------------------------------------------------------------------
  subroutine displayParserError(line_number, column, ierror)
    integer, intent(in) :: line_number, column, ierror

    if (ierror /= 0) then
      write(*,'(9(A,I0))') 'CONFIG_READER: Error in line ', line_number, ' in column ', column
      select case (ierror)
        case (CONFIG_READER_ERR_BAD_CHAR)      
      write (*,*) "Character not allowed."
        case (CONFIG_READER_ERR_NO_EQUAL)      
          write (*,*) "Expected equal sign not found."
        case (CONFIG_READER_ERR_NO_VALUE)      
          write (*,*) "No value assigned to variable."
        case (CONFIG_READER_ERR_EMPTY_STR)     
          write (*,*) "Value is empty string."
        case (CONFIG_READER_ERR_END_LINE)      
          write (*,*) "Unexpected end of line."
        case (CONFIG_READER_ERR_VAR_NOT_UNIQUE);
          write (*,*) "Parameter is already defined."
        case default
          write (*,*) "Unknown error."
      end select
    end if
  end subroutine displayParserError

!---------------------------------------------------------------------
  subroutine getValueString(this, variable, value, ierror)
    use ConfigReaderDictionary_mod, only: getDictionaryValue

    type (ConfigReader), intent(inout) :: this
    character(len = *), intent(in) :: variable
    character(len = *), intent(inout) :: value
    integer, intent(out) :: ierror

    logical :: tag
    integer :: dict_error

    tag = .true.
    dict_error = 0
    ierror = 0
    call getDictionaryValue(this%parse_dict, variable, value, tag, dict_error)

    if (dict_error /= 0) then
      ierror = CONFIG_READER_ERR_VAR_NOT_FOUND
    end if
  end subroutine getValueString

!---------------------------------------------------------------------
!> Get value of 'variable' and interpret it as integer value.
!>
!> a default value can be passed as int_value, which does not change on exit
!> if the variable is not found
  subroutine getValueInteger(this, variable, int_value, ierror)
    type (ConfigReader), intent(inout) :: this
    character(len = *), intent(in) :: variable
    integer, intent(inout) :: int_value
    integer, intent(out) :: ierror

    character(len = MAX_LINE_LENGTH) :: value_string
    integer :: value_read
    integer :: ios

    ierror = 0
    ios = 0

    call getValueString(this, variable, value_string, ierror)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) value_read
      ! test if an integer was read
      if (ios == 0) then
        int_value = value_read
      else
        ierror = CONFIG_READER_ERR_NOT_INTEGER
      end if
    end if

  end subroutine getValueInteger

!---------------------------------------------------------------------
!> Get value of 'variable' and interpret it as double prec. value.
!>
!> a default value can be passed as double_value, which does not change on exit
!> if the variable is not found
  subroutine getValueDouble(this, variable, double_value, ierror)
    type (ConfigReader), intent(inout) :: this
    character(len = *), intent(in) :: variable
    real(kind = kind(1.0d0)), intent(inout) :: double_value
    integer, intent(out) :: ierror

    character(len = MAX_LINE_LENGTH) :: value_string
    real(kind = kind(1.0d0)) :: value_read
    integer :: ios

    ierror = 0
    ios = 0

    call getValueString(this, variable, value_string, ierror)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) value_read
      ! test if a double was read
      if (ios == 0) then
        double_value = value_read
      else
        ierror = CONFIG_READER_ERR_NOT_DOUBLE
      end if
    end if

  end subroutine getValueDouble

!---------------------------------------------------------------------
!> Get value of 'variable' and interpret it as logical value.
!>
! a default value can be passed as logical_value, which does not change on exit
! if the variable is not found
  subroutine getValueLogical(this, variable, logical_value, ierror)
    type (ConfigReader), intent(inout) :: this
    character(len = *), intent(in) :: variable
    logical, intent(inout) :: logical_value
    integer, intent(out) :: ierror

    character(len = MAX_LINE_LENGTH) :: value_string
    logical :: value_read
    integer :: ios

    ierror = 0
    ios = 0

    call getValueString(this, variable, value_string, ierror)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) value_read
      ! test if a logical was read
      if (ios == 0) then
        logical_value = value_read
      else
        ierror = CONFIG_READER_ERR_NOT_LOGICAL
      end if
    end if

  end subroutine getValueLogical

!---------------------------------------------------------------------
!> Reads a double precision vector of fixed length from config-File.
!>
!> The dimension of the vector has to be passed as argument 'length'
!> The routine does not check if more then 'length' values are present
!> If the vector is too short, however, ierror contains an error-code
!
!> A default value can be passed as vector, which does not change on exit
!> if the variable is not found
  subroutine getValueDoubleVector(this, variable, vector, length, ierror)
    type (ConfigReader), intent(inout) :: this
    character(len = *), intent(in) :: variable
    double precision, dimension(length), intent(inout) :: vector
    integer, intent(in) :: length
    integer, intent(out) :: ierror

    character(len = MAX_LINE_LENGTH) :: value_string
    double precision, dimension(length) :: vector_read
    integer :: ios
    integer :: ii

    ierror = 0
    ios = 0

    call getValueString(this, variable, value_string, ierror)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) (vector_read(ii), ii = 1, length)
      ! test if read was successful
      if (ios == 0) then
        vector = vector_read
      else
        write(*,*) value_string
        write(*,*) vector_read
        vector = vector_read
        write(*,*) ios
        ierror = CONFIG_READER_ERR_NOT_DOUBLE
      end if
    end if

  end subroutine getValueDoubleVector

!---------------------------------------------------------------------
!> Reads an integer vector of fixed length from config-File.
!>
!> The dimension of the vector has to be passed as argument 'length'
!> The routine does not check if more then 'length' values are present
!> If the vector is too short, however, ierror contains an error-code
!
!> A default value can be passed as vector, which does not change on exit
!> if the variable is not found
  subroutine getValueIntVector(this, variable, vector, length, ierror)
    type (ConfigReader), intent(inout) :: this
    character(len = *), intent(in) :: variable
    integer, dimension(length), intent(inout) :: vector
    integer, intent(in) :: length
    integer, intent(out) :: ierror

    character(len = MAX_LINE_LENGTH) :: value_string
    integer, dimension(length) :: vector_read
    integer :: ios
    integer :: ii

    ierror = 0
    ios = 0

    call getValueString(this, variable, value_string, ierror)

    if (ierror == 0) then
      read(unit=value_string, fmt=*, iostat=ios) (vector_read(ii), ii = 1, length)
      ! test if read was successful
      if (ios == 0) then
        vector = vector_read
      else
        write(*,*) value_string
        write(*,*) vector_read
        vector = vector_read
        write(*,*) ios
        ierror = CONFIG_READER_ERR_NOT_INTEGER
      end if
    end if

  end subroutine getValueIntVector

!---------------------------------------------------------------------
!> This subroutine is used to get variables, which have not been read (yet).
!>
!> Initially an integer with value 1 has to be passed as next_ptr
!> next_ptr is an integer which has to be kept for the next call to this
!> subroutine, which gives the next unread variable
!> If no unread variable was found then ierror = CONFIG_READER_ERR_VAR_NOT_FOUND
  subroutine getUnreadVariable(this, variable, next_ptr, ierror)
    use ConfigReaderDictionary_mod, only: getTaggedVariable, CONFIG_READER_DICT_NOT_FOUND

    type (ConfigReader), intent(in) :: this
    character(len = *),  intent(out) :: variable
    integer, intent(inout) :: next_ptr
    integer, intent(out) :: ierror

    integer :: dict_error

    variable = ' '
    ierror = 0
    dict_error = 0

    call getTaggedVariable(this%parse_dict, variable, .false., next_ptr, dict_error)

    if (dict_error == CONFIG_READER_DICT_NOT_FOUND) then
      ierror = CONFIG_READER_ERR_VAR_NOT_FOUND
    end if

  end subroutine getUnreadVariable

end module ConfigReader_mod

