!----------------------------------------------------------------------------
!> Implementation of a Dictionary
!> given a name 'variable', look up the corresponding 'value'
!> This is a slow implementation, using a dynamic array
!> (a faster implementation should i.e. use a hash table)

!> @author Elias Rabel, September 2011

module ConfigReaderDictionary_mod
  implicit none
  private
  public :: Dictionary, create, destroy
  public :: createDictionary, destroyDictionary ! deprecated
  public :: pushBackDictionary, getDictionaryValue, getTaggedVariable
  
  ! Status flags and error codes
  integer, parameter, public :: CONFIG_READER_DICT_NOT_UNIQUE=1, CONFIG_READER_DICT_NOT_FOUND=2

  ! Maximal length of variable names and value string.
  ! Modify to allow different lengths.
  integer, parameter, public :: CONFIG_READER_DICT_VAR_LENGTH=64, CONFIG_READER_DICT_VALUE_LENGTH=191

  type, private :: DictionaryEntry
    character(len=CONFIG_READER_DICT_VAR_LENGTH) :: variable
    character(len=CONFIG_READER_DICT_VALUE_LENGTH) :: value
    logical(kind=1) :: tag
  endtype
  
  type Dictionary
    private
    integer :: counter
    type(DictionaryEntry), allocatable :: dict(:)
  endtype 
  

  integer, parameter, private :: INITIAL_SIZE = 256
  
  interface create
    module procedure createDictionary
  endinterface
  
  interface destroy
    module procedure destroyDictionary
  endinterface
  
  contains

  subroutine fatalErrorDictionary(message)
    character(len=*), intent(in), optional :: message
    write(*,*) "ConfigReaderDictionary_mod: << Fatal Error >>"

    if (present(message)) write(*,*) message

    stop
  endsubroutine

  subroutine createDictionary(this)
    type(Dictionary), intent(inout) :: this

    integer :: ierror

    this%counter = 0

    if (allocated(this%dict)) call fatalErrorDictionary("It seems that the dictionary was already created.")

    allocate(this%dict(INITIAL_SIZE), stat=ierror)
    if (ierror /= 0) call fatalErrorDictionary()
  endsubroutine createDictionary


  !----------------------------------------------------------------------------
  !> if an entry with the same name as in 'variable'
  !> exists then ierror = CONFIG_READER_DICT_NOT_UNIQUE and the variable/value
  !> pair is not inserted in the dictionary
  subroutine pushBackDictionary(this, variable, value, tag, ierror)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: variable
    character(len=*), intent(in) :: value
    logical, intent(in) :: tag

    integer, intent(out) :: ierror

    type(DictionaryEntry), allocatable :: dict_tmp(:)

    integer :: ind, capacity, new_capacity, ios

    if (.not. allocated(this%dict)) call fatalErrorDictionary("Dictionary was not created.")

    ierror = 0

    ! check if variable already exists
    do ind = 1, this%counter
      if (variable == this%dict(ind)%variable) then
        ierror = CONFIG_READER_DICT_NOT_UNIQUE
        return
      endif
    enddo ! ind

    ind = this%counter + 1 ! index of next entry

    capacity = size(this%dict)

    ! if array is full, reallocate memory
    if (ind > capacity) then
      new_capacity = 2*capacity ! size doubling
      allocate(dict_tmp(capacity), stat=ios)
      if (ios /= 0) call fatalErrorDictionary()

      dict_tmp(1:capacity) = this%dict(:)
      
      deallocate(this%dict, stat=ios) ! ignore status
      allocate(this%dict(new_capacity), stat=ios)
      if (ios /= 0) call fatalErrorDictionary()
      
      this%dict(1:capacity) = dict_tmp(:)
      ! this%dict(capacity+1:) = DictionaryEntry("", "", .true.) ! do we need to init this array part?

      deallocate(dict_tmp, stat=ios) ! ignore status
    endif

    this%dict(ind)%variable = variable
    this%dict(ind)%value = value
    this%dict(ind)%tag = tag
    
    this%counter = this%counter + 1 ! increment the counter

  endsubroutine


  !----------------------------------------------------------------------------
  !> retrieves the value corresponding to 'variable' and sets the
  !> tag to the given logical value (true/false)
  subroutine getDictionaryValue(this, variable, value, tag, ierror)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: variable
    character(len=*), intent(inout) :: value
    logical, intent(in) :: tag
    integer, intent(out) :: ierror

    integer :: ind

    ierror = CONFIG_READER_DICT_NOT_FOUND

    do ind = 1, this%counter
      if (variable == this%dict(ind)%variable) then
        ierror = 0
        value = this%dict(ind)%value
        this%dict(ind)%tag = tag ! set the tag
        return
      endif
    enddo ! ind
    
  endsubroutine

  !----------------------------------------------------------------------------
  !> This routine can be used to find variables that have never been looked up.
  !>
  !> Use this to get all variable names which are tagged (with logical value given in 'tag')
  !> next_ptr is an integer which points to the next value that should be
  !> checked, set it to 1 if you want to start searching at the beginning
  !> if next_ptr has an illegal value or no tagged variable was found
  !> then ierror = CONFIG_READER_DICT_NOT_FOUND
  !> variable then has the value it had on entry!
  subroutine getTaggedVariable(this, variable, tag, next_ptr, ierror)
    type(Dictionary), intent(in) :: this
    character(len=*), intent(inout) :: variable
    logical, intent(in) :: tag
    integer, intent(inout) :: next_ptr
    integer, intent(out) :: ierror

    integer :: ind

    ierror = CONFIG_READER_DICT_NOT_FOUND

    if (next_ptr < 1 .or. next_ptr > this%counter) return

    do ind = next_ptr, this%counter
      if (this%dict(ind)%tag .eqv. tag) then
        ierror = 0
        variable = this%dict(ind)%variable
        exit
      endif
    enddo ! ind

    next_ptr = ind + 1

  endsubroutine

  !> output of dictionary content to stdout
  !> useful for testing purposes
  subroutine printDictionary(this)
    type(Dictionary), intent(in) :: this

    integer :: ind
    
    do ind = 1, this%counter
      write(*,*) trim(this%dict(ind)%variable), " ", trim(this%dict(ind)%value), this%dict(ind)%tag
    enddo ! ind
    
  endsubroutine ! print

  !----------------------------------------------------------------------------
  !> Frees resources used by dictionary.
  !>
  !> This routine has to be called after use of the dictionary
  subroutine destroyDictionary(this)
    type(Dictionary), intent(inout) :: this

    integer :: ierror

    if (.not. allocated(this%dict)) &
      call fatalErrorDictionary("Dictionary was already destroyed.")

    deallocate(this%dict, stat=ierror)
    if (ierror /= 0) call fatalErrorDictionary()
  endsubroutine destroyDictionary

endmodule ConfigReaderDictionary_mod
