#define status_t integer

module StringHelpers_mod
implicit none
  private ! default visibility

  public :: operator(+), operator(-)
  public :: test, to_lowercase, to_int, replace_underscore
  
  integer, parameter, public :: DEFAULT_STRING_LENGTH=128 !< change length of default strings here
  integer, parameter :: SLEN =  DEFAULT_STRING_LENGTH
#define string_t character(len=SLEN)
  
  interface operator(+) ! left," ",right
    module procedure concat_string_whitespace_string, &
      concat_string_whitespace_integer, concat_string_whitespace_real, concat_string_whitespace_double, &
      concat_string_whitespace_integers, concat_string_whitespace_doubles
  endinterface

  interface operator(-) ! left,right ! without space in between
    module procedure concat_string_string, concat_string_integer, concat_integer_string
  endinterface

  contains

  string_t function concat_string_whitespace_string(left, right) result(str)
    character(len=*), intent(in) :: left, right
    status_t :: ios
    write(unit=str,fmt='(9A)',iostat=ios) trim(left)," ",trim(adjustl(right))
  endfunction ! string," ",string

  string_t function concat_string_string(left, right) result(str)
    character(len=*), intent(in) :: left, right
    status_t :: ios
    write(unit=str,fmt='(9A)',iostat=ios) trim(left),trim(adjustl(right))
  endfunction ! string,string

  string_t function concat_string_whitespace_integer(left, right) result(str)
    character(len=*), intent(in) :: left
    integer, intent(in)          :: right
    status_t :: ios
    write(unit=str,fmt='(2A,I0)',iostat=ios) trim(left),' ',right
  endfunction ! string," ",integer
  
  string_t function concat_string_whitespace_integers(left, right) result(str)
    character(len=*), intent(in) :: left
    integer, intent(in)          :: right(:)
    status_t :: ios
    write(unit=str,fmt='(A,9999(" ",I0))',iostat=ios) trim(left),right
  endfunction ! string," ",integer(1)," ",integer(2), ...
  
  string_t function concat_string_integer(left, right) result(str)
    character(len=*), intent(in) :: left
    integer, intent(in)          :: right
    status_t :: ios
    write(unit=str,fmt='(A,I0)',iostat=ios) trim(left),right
  endfunction ! string,integer

  string_t function concat_integer_string(left, right) result(str)
    integer, intent(in)          :: left
    character(len=*), intent(in) :: right
    status_t :: ios
    write(unit=str,fmt='(A,I0)',iostat=ios) left,adjustl(right)
  endfunction ! integer,string

  string_t function concat_string_whitespace_real(left, right) result(str)
    character(len=*), intent(in) :: left
    real(kind=4), intent(in) :: right
    status_t :: ios
    write(unit=str,fmt='(2A,F0.6)',iostat=ios) trim(left),' ',right
  endfunction ! string," ",real

  string_t function concat_string_whitespace_double(left, right) result(str)
    character(len=*), intent(in) :: left
    double precision, intent(in) :: right
    status_t :: ios
    write(unit=str,fmt='(2A,F0.6)',iostat=ios) trim(left),' ',right
  endfunction ! string," ",double
  
  string_t function concat_string_whitespace_doubles(left, right) result(str)
    character(len=*), intent(in) :: left
    double precision, intent(in) :: right(:)
    status_t :: ios
    write(unit=str,fmt='(A,9999(" ",F0.6))',iostat=ios) trim(left),right
  endfunction ! string," ",double(1)," ",double(2), ...
  
  integer function to_lowercase(mixedcase, lowercase) result(n_chars_converted)
    character(len=*), intent(in)  :: mixedcase !! string containing possible mixed case characters
    character(len=*), intent(out) :: lowercase !! string with only lower case characters

    integer, parameter :: A=ichar('A'), Z=ichar('Z'), aA=ichar('a')-A 
    integer            :: ip, ic, lm
    n_chars_converted = 0 ! init
    lowercase = mixedcase ! string copy
    lm = len_trim(mixedcase)
    if (len(lowercase) < lm) stop 'to_lowercase: string is cropped!'
    do ip = 1, lm
      ic = ichar(mixedcase(ip:ip))
      if (ic < A) cycle
      if (ic > Z) cycle
      lowercase(ip:ip) = achar(ic + aA)
      n_chars_converted = n_chars_converted+1 ! count up
    enddo ! ip
  endfunction ! to_lowercase

  integer function to_int(str, iostat)
    character(len=*), intent(in) :: str
    integer, intent(out), optional :: iostat 
    status_t :: ios
    read(unit=str,fmt=*,iostat=ios) to_int
    if (present(iostat)) iostat = ios
  endfunction ! to_int
  
  string_t function replace_underscore(s) result(r)
    character(len=*), intent(in) :: s
    integer :: i, n
    r = s ! copy
    n = 0
    do i = 1, len_trim(s)
      if (s(i:i) == '_') then
        if(i+n+2 > len(r)) return ! avoid out-of-bounds
        r(i+n+2:) = s(i+1:) ! copy
        r(i+n:i+n+1) = '\_'
        n = n+1
      endif
    enddo ! i
  endfunction ! replace_underscore

  status_t function test()
    write(*,'(A)',iostat=test) "s"+"s" , "i"+1, "i"+0.123d0
    write(*,'(A)',iostat=test) "s"-"s"!, "i"-1, "i"-0.123d0
    write(*,'(9A)',iostat=test) 'replace_underscore("_a_Bc_DeF_gHiJ_kLmNo_PqR...") = "',trim(replace_underscore("_a_Bc_DeF_gHiJ_kLmNo_PqR...")),'"'
    write(*,'(9A)',iostat=test) 'replace_underscore("a_Bc_DeF_gHiJ_kLmNo_PqR...") = "',trim(replace_underscore("a_Bc_DeF_gHiJ_kLmNo_PqR...")),'"'
  endfunction ! test

endmodule ! StringHelpers_mod