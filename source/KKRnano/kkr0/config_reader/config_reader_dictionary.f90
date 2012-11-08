  !----------------------------------------------------------------------------
!> Implementation of a Dictionary
!> given a name 'variable', look up the corresponding 'value'
!> This is a slow implementation, using a dynamic array
!> (a faster implementation should i.e. use a hash table)

!> @author Elias Rabel, September 2011

module Config_Reader_Dictionary
  implicit none

  ! Status flags and error codes
  integer, parameter, public :: CONFIG_READER_DICT_NOT_UNIQUE = 1
  integer, parameter, public :: CONFIG_READER_DICT_NOT_FOUND = 2

  ! Maximal length of variable names and value string.
  ! Modify to allow different lengths.
  integer, parameter, public :: CONFIG_READER_DICT_VAR_LENGTH = 64
  integer, parameter, public :: CONFIG_READER_DICT_VALUE_LENGTH = 192

  type Dictionary
    private
    type (DictionaryEntry), dimension(:), pointer :: dict
    integer :: counter
  end type Dictionary

! private declarations

  integer, parameter, private :: INITIAL_SIZE = 50

  type, private :: DictionaryEntry
    character(len=CONFIG_READER_DICT_VAR_LENGTH) :: variable
    character(len=CONFIG_READER_DICT_VALUE_LENGTH) :: value
    logical :: tag
  end type DictionaryEntry

contains

  subroutine fatalErrorDictionary(message)
    character(len=*), intent(in), optional :: message
    write(*,*) "Config_Reader_Dictionary: << Fatal Error >>"

    if (present(message)) then
      write(*,*) message
    end if

    stop
  end subroutine

  subroutine createDictionary(this)
    type (Dictionary), intent(inout) :: this

    integer :: ierror

    this%counter = 0

    if (associated(this%dict)) then
      call fatalErrorDictionary("It seems that the dictionary was already created.")
    end if

    allocate(this%dict(INITIAL_SIZE), stat=ierror)
    if (ierror /= 0) then
      call fatalErrorDictionary()
    end if
  end subroutine createDictionary


  !----------------------------------------------------------------------------
  !> if an entry with the same name as in 'variable'
  !> exists then ierror = CONFIG_READER_DICT_NOT_UNIQUE and the variable/value
  !> pair is not inserted in the dictionary
  subroutine pushBackDictionary(this, variable, value, tag, ierror)
    type (Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: variable
    character(len=*), intent(in) :: value
    logical, intent(in) :: tag

    integer, intent(out) :: ierror

    type (DictionaryEntry), dimension(:), pointer :: dict_new

    integer :: ind
    integer :: capacity
    integer :: ios
    integer :: loop_ind

    if (.not. associated(this%dict)) then
      call fatalErrorDictionary("Dictionary was not created.")
    end if

    ierror = 0

    ! check if variable already exists
    do ind = 1, this%counter
      if (variable .eq. this%dict(ind)%variable) then
        ierror = CONFIG_READER_DICT_NOT_UNIQUE
        return
      end if
    end do

    ! index of next entry
    ind = this%counter + 1

    capacity = size(this%dict)

    ! if array is full, reallocate memory
    if (ind > capacity) then
      allocate(dict_new(capacity * 2), stat = ios)
      if (ios /= 0) then
        call fatalErrorDictionary()
      end if

      do loop_ind = 1, capacity
        dict_new(loop_ind) = this%dict(loop_ind)
      end do

      deallocate(this%dict, stat = ios)

      this%dict => dict_new
      if (ios /= 0) then
        call fatalErrorDictionary()
      end if
    end if

    this%dict(ind)%variable = variable
    this%dict(ind)%value = value
    this%dict(ind)%tag = tag

    ! don't forget to increment counter
    this%counter = this%counter + 1

  end subroutine


  !----------------------------------------------------------------------------
  !> retrieves the value corresponding to 'variable' and sets the
  !> tag to the given logical value (true/false)
  subroutine getDictionaryValue(this, variable, value, tag, ierror)
    type (Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: variable
    character(len=*), intent(inout) :: value
    logical, intent(in) :: tag
    integer, intent(out) :: ierror

    integer :: ind

    ierror = CONFIG_READER_DICT_NOT_FOUND

    do ind = 1, this%counter
      if (variable .eq. this%dict(ind)%variable) then
        ierror = 0
        value = this%dict(ind)%value
        this%dict(ind)%tag = tag     ! set the tag
        return
      end if
    end do
  end subroutine

  !----------------------------------------------------------------------------
  !> This routine can be used to find variables that have never been looked up.
  !>
  !> Use this to get all variables name which are tagged (with logical
  !> value given in 'tag')
  !> next_ptr is an integer which points to the next value that should be
  !> checked, set it to 1 if you want to start searching at the beginning
  !> if next_ptr has an illegal value or no tagged variable was found
  !> then ierror = CONFIG_READER_DICT_NOT_FOUND
  !> variable then has the value it had on entry!
  subroutine getTaggedVariable(this, variable, tag, next_ptr, ierror)
    type (Dictionary), intent(in) :: this
    character(len=*), intent(inout) :: variable
    logical, intent(in) :: tag
    integer, intent(inout) :: next_ptr
    integer, intent(out) :: ierror

    integer :: ind

    ierror = CONFIG_READER_DICT_NOT_FOUND

    if (next_ptr < 1 .or. next_ptr > this%counter) then
      return
    end if

    do ind = next_ptr, this%counter
      if (this%dict(ind)%tag .eqv. tag) then
        ierror = 0
          variable = this%dict(ind)%variable
        exit
      end if
    end do

    next_ptr = ind + 1

  end subroutine

  !> output of dictionary content to stdout
  !> useful for testing purposes
  subroutine printDictionary(this)
    type (Dictionary), intent(in) :: this

    integer :: ind

    do ind = 1, this%counter
      write (*,*)  trim(this%dict(ind)%variable), " ", &
                 & trim(this%dict(ind)%value), &
                 & this%dict(ind)%tag
    end do
  end subroutine

  !----------------------------------------------------------------------------
  !> Frees resources used by dictionary.
  !>
  !> This routine has to be called after use of the dictionary
  subroutine destroyDictionary(this)
    type (Dictionary), intent(inout) :: this

    integer :: ierror

    if (.not. associated(this%dict)) then
      call fatalErrorDictionary("Dictionary was already destroyed.")
    end if

    deallocate(this%dict, stat=ierror)
    if (ierror /= 0) then
      call fatalErrorDictionary()
    end if
  end subroutine destroyDictionary

end module Config_Reader_Dictionary
