module m_npy
   implicit none

   integer(4), parameter               :: p_un = 23
   character, parameter                :: magic_num = achar(147) ! x93
   character, parameter                :: major = achar(2)   !major *.npy version
   character, parameter                :: minor = achar(0)   !minor *.npy version
   character(len=*), parameter         :: zip_flag = "-q0"
   character(len=*), parameter         :: magic_str = "NUMPY"

   interface save_npy
      module procedure write_int64_vec, write_int64_mtx, &
         write_int32_vec, write_int32_mtx, write_int32_3d, &
         write_int16_vec, write_int16_mtx, &
         write_int8_vec, write_int8_mtx, write_int8_3d, &
         write_dbl_vec, write_dbl_mtx, &
         write_sng_vec, write_sng_mtx, &
         write_cmplx_sgn_vec, write_cmplx_sgn_mtx, &
         write_cmplx_dbl_vec, write_cmplx_dbl_mtx, &
         write_sng_3dT, write_dbl_3dT, &
         write_sng_4dT, write_dbl_4dT, &
         write_dbl_5dT, &
         write_cmplx_dbl_3dT, &
         write_cmplx_dbl_4dT, &
         write_cmplx_dbl_5dT, &
         write_cmplx_dbl_6dT

   end interface save_npy
   interface add_npz
      module procedure addrpl_int8_vec, addrpl_int8_mtx, &
         addrpl_int16_vec, addrpl_int16_mtx, &
         addrpl_int32_vec, addrpl_int32_mtx, &
         addrpl_int64_vec, addrpl_int64_mtx, &
         addrpl_sng_vec, addrpl_sng_mtx, &
         addrpl_dbl_vec, addrpl_dbl_mtx, &
         addrpl_cmplx_dbl_vec, addrpl_cmplx_dbl_mtx, &
         addrpl_cmplx_sng_vec, addrpl_cmplx_sng_mtx
   end interface add_npz

contains
   subroutine run_sys(cmd, stat)
      implicit none
      character(len=*), intent(in)     :: cmd
      integer(4), intent(out)          :: stat

      call execute_command_line(cmd, wait=.True., exitstat=stat)
   end subroutine run_sys

   subroutine addrpl_cmplx_sng_vec(zipfile, var_name, vec)
      implicit none
      complex(4), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_sng_vec

   subroutine addrpl_cmplx_sng_mtx(zipfile, var_name, mtx)
      implicit none
      complex(4), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_sng_mtx

   subroutine addrpl_cmplx_dbl_vec(zipfile, var_name, vec)
      implicit none
      complex(8), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_dbl_vec

   subroutine addrpl_cmplx_dbl_mtx(zipfile, var_name, mtx)
      implicit none
      complex(8), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_dbl_mtx

   subroutine addrpl_dbl_vec(zipfile, var_name, vec)
      implicit none
      real(8), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_dbl_vec

   subroutine addrpl_dbl_mtx(zipfile, var_name, mtx)
      implicit none
      real(8), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_dbl_mtx

   subroutine addrpl_sng_vec(zipfile, var_name, vec)
      implicit none
      real(4), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_sng_vec

   subroutine addrpl_sng_mtx(zipfile, var_name, mtx)
      implicit none
      real(4), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_sng_mtx

   subroutine addrpl_int8_vec(zipfile, var_name, vec)
      implicit none
      integer(1), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int8_vec

   subroutine addrpl_int8_mtx(zipfile, var_name, mtx)
      implicit none
      integer(1), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int8_mtx

   subroutine addrpl_int16_vec(zipfile, var_name, vec)
      implicit none
      integer(2), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int16_vec

   subroutine addrpl_int16_mtx(zipfile, var_name, mtx)
      implicit none
      integer(2), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int16_mtx

   subroutine addrpl_int32_vec(zipfile, var_name, vec)
      implicit none
      integer(4), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int32_vec

   subroutine addrpl_int32_mtx(zipfile, var_name, mtx)
      implicit none
      integer(4), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int32_mtx

   subroutine addrpl_int64_vec(zipfile, var_name, vec)
      implicit none
      integer(8), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int64_vec

   subroutine addrpl_int64_mtx(zipfile, var_name, mtx)
      implicit none
      integer(8), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int64_mtx

   Subroutine write_cmplx_sgn_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(4), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<c8"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor
      write (p_un) header_len
      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_cmplx_sgn_mtx

   Subroutine write_cmplx_sgn_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(4), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<c8"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor
      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_cmplx_sgn_vec

   Subroutine write_cmplx_dbl_6dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :, :, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_6dT

   Subroutine write_cmplx_dbl_5dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_5dT

   Subroutine write_cmplx_dbl_4dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_4dT

   Subroutine write_cmplx_dbl_3dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_3dT

   Subroutine write_cmplx_dbl_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_cmplx_dbl_mtx

   Subroutine write_cmplx_dbl_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_cmplx_dbl_vec

   Subroutine write_sng_3dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: tensor(:, :, :)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_sng_3dT

   Subroutine write_sng_4dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: tensor(:, :, :, :)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_sng_4dT

   Subroutine write_sng_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_sng_mtx

   Subroutine write_sng_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: vec(:)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_sng_vec

   Subroutine write_dbl_3dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: tensor(:, :, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_dbl_3dT

   Subroutine write_dbl_4dT(filename, tensor4)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: tensor4(:, :, :, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor4)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor4))
      write (p_un) tensor4
      close (unit=p_un)
   End Subroutine write_dbl_4dT

   Subroutine write_dbl_5dT(filename, tensor5)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: tensor5(:, :, :, :, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor5)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor5))
      write (p_un) tensor5
      close (unit=p_un)
   End Subroutine write_dbl_5dT

   Subroutine write_dbl_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_dbl_mtx

   Subroutine write_dbl_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: vec(:)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_dbl_vec

   Subroutine write_int64_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(8), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i8"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int64_mtx

   Subroutine write_int64_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(8), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<i8"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_int64_vec

   Subroutine write_int32_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(4), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i4"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int32_mtx

   Subroutine write_int32_3d(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(4), intent(in)           :: mtx(:,:,:)
      character(len=*), parameter      :: var_type = "<i4"
      integer(4)                       :: header_len, s_mtx(3), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int32_3d

   Subroutine write_int32_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(4), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<i4"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_int32_vec

   Subroutine write_int16_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(2), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i2"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int16_mtx

   Subroutine write_int16_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(2), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<i2"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_int16_vec

   Subroutine write_int8_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(1), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i1"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx


module mod_write_gflle 

  implicit none

contains

  !-------------------------------------------------------------------------------
  !> Summary: Write gflle file out in npy format
  !> Author: Philipp Rüßmann
  !> Category: writeout
  !> Deprecated: False 
  !> Creates one file per atom and energy, otherwise file can be very large which
  !> might be problematic in post-processing
  !-------------------------------------------------------------------------------
  subroutine write_gflle_to_npy(lmmaxd, ielast, nqdos, i1, gflle)

    use mod_datatypes, only: dp
    use m_npy, only: save_npy
    implicit none
    integer, intent(in) :: lmmaxd, ielast, nqdos, i1
    complex (kind=dp) :: gflle(lmmaxd,lmmaxd,ielast,nqdos)
    character (len=100) :: filename
    integer :: ie
    do ie = 1, ielast
      write(filename, "(A,1I0.3,A,1I0.3,A)") "gllke.", I1, ".", IE, ".npy"
      call save_npy(trim(filename), gflle(:,:,ie, :))
    end do

  end subroutine write_gflle_to_npy

end module mod_write_gflle 
