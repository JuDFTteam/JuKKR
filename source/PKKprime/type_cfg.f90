module type_cfg

use mpi

implicit none

  type :: cfg_TYPE

!    sequence

     !init
     integer          :: N = 30

     !general output parameter
     logical          :: debug    = .false.
     logical          :: verbose  = .false.

     !parameter for Fermi-surface calculation
     integer          :: nOct(3) = 1
     double precision :: filter_eigval  = 1d0
     double precision :: eps_degenerate = 1d-04
     double precision :: max_realW      = 1d-01

     !parameter for calculation on the Fermi surface
     logical          :: log_saveeigv = .false.
     logical          :: log_pdos = .false.
     logical          :: log_atomresolved = .false.
     integer          :: nAniso  = 1
     integer          :: spincomb   = -1
     integer          :: nspincomb  = -1

     !parameter for integration on Fermi surface
     integer          :: intpmeth   = -1
     integer          :: iregion    = -1
     integer          :: hist_nbins = -1
     integer          :: hist_logscale = -1
     double precision :: hist_xmin  = -1e16
     double precision :: hist_xmax  =  1e16

     !switches for scattering-mode
     logical          :: log_visdata      = .false.
     logical          :: log_Pkkfixed     = .false.
     logical          :: log_Pkksave      = .false.
     logical          :: log_lifetime     = .false.
     logical          :: log_conductivity = .false.

     !mode selection
     logical          :: log_band = .false.
     logical          :: log_fs   = .true.
     logical          :: log_fv   = .false.
     logical          :: log_spin = .false.
     logical          :: log_intg = .false.
     logical          :: log_scattering   = .false.

  end type cfg_TYPE


  type :: cfga_TYPE

     integer :: N = 2

     integer, allocatable :: ispincomb(:)

  end type cfga_TYPE

contains



  subroutine myMPI_create_struct_cfg(cfg,myMPItype,iout)

    type(cfg_type), intent(in)  :: cfg
    integer,        intent(out) :: myMPItype, iout

    integer :: N, blocklen(cfg%N), etype(cfg%N), ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp(cfg%N), base

    N = cfg%N

    !init
    call MPI_Get_address(cfg%N, disp(1), ierr)

    !general output parameter
    call MPI_Get_address(cfg%debug,   disp(2), ierr)
    call MPI_Get_address(cfg%verbose, disp(3), ierr)

    !parameter for Fermi-surface calculation
    call MPI_Get_address(cfg%nOct,           disp(4), ierr)
    call MPI_Get_address(cfg%filter_eigval,  disp(5), ierr)
    call MPI_Get_address(cfg%eps_degenerate, disp(6), ierr)
    call MPI_Get_address(cfg%max_realW,      disp(7), ierr)

    !parameter for calculation on the Fermi surface
    call MPI_Get_address(cfg%log_saveeigv,     disp(8),  ierr)
    call MPI_Get_address(cfg%log_pdos,         disp(9),  ierr)
    call MPI_Get_address(cfg%log_atomresolved, disp(10), ierr)
    call MPI_Get_address(cfg%nAniso,           disp(11), ierr)
    call MPI_Get_address(cfg%spincomb,         disp(12), ierr)
    call MPI_Get_address(cfg%nspincomb,        disp(13), ierr)

    !parameter for integration on Fermi surface
    call MPI_Get_address(cfg%intpmeth,      disp(14), ierr)
    call MPI_Get_address(cfg%iregion,       disp(15), ierr)
    call MPI_Get_address(cfg%hist_nbins,    disp(16), ierr)
    call MPI_Get_address(cfg%hist_logscale, disp(17), ierr)
    call MPI_Get_address(cfg%hist_xmin,     disp(18), ierr)
    call MPI_Get_address(cfg%hist_xmax,     disp(19), ierr)

    !switches for scattering-mode
    call MPI_Get_address(cfg%log_visdata,      disp(20), ierr)
    call MPI_Get_address(cfg%log_Pkkfixed,     disp(21), ierr)
    call MPI_Get_address(cfg%log_Pkksave,      disp(22), ierr)
    call MPI_Get_address(cfg%log_lifetime,     disp(23), ierr)
    call MPI_Get_address(cfg%log_conductivity, disp(24), ierr)

    !mode selection
    call MPI_Get_address(cfg%log_band,       disp(25), ierr)
    call MPI_Get_address(cfg%log_fs,         disp(26), ierr)
    call MPI_Get_address(cfg%log_fv,         disp(27), ierr)
    call MPI_Get_address(cfg%log_spin,       disp(28), ierr)
    call MPI_Get_address(cfg%log_intg,       disp(29), ierr)
    call MPI_Get_address(cfg%log_scattering, disp(30), ierr)


    base = disp(1)
    disp = disp - base

    blocklen(1:N) = 1
    blocklen(4)   = 3

    !init
    etype(1) = MPI_INTEGER

    !general output parameter
    etype(2:3) = MPI_LOGICAL

    !parameter for Fermi-surface calculation
    etype(4) = MPI_INTEGER
    etype(5:7)   = MPI_DOUBLE_PRECISION

    !parameter for calculation on the Fermi surface
    etype(8:10)  = MPI_LOGICAL
    etype(11:13) = MPI_INTEGER

    !parameter for integration on Fermi surface
    etype(14:17) = MPI_INTEGER
    etype(18:19) = MPI_DOUBLE_PRECISION

    !switches for scattering-mode
    etype(20:24)  = MPI_LOGICAL

    !mode selection
    etype(25:30)  = MPI_LOGICAL

    call MPI_Type_create_struct(N, blocklen, disp, etype, myMPItype, iout)

  end subroutine myMPI_create_struct_cfg





  subroutine myMPI_create_struct_cfga(cfg,cfga,myMPItype,iout)
    type(cfg_type),  intent(in)    :: cfg
    type(cfga_type), intent(inout) :: cfga
    integer,         intent(out)   :: myMPItype, iout

    integer :: N, blocklen(cfga%N), etype(cfga%N), ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp(cfga%N), base

    N = cfga%N

    if(.not.allocated(cfga%ispincomb)) then
      allocate( cfga%ispincomb(cfg%nspincomb), STAT=ierr )
      if(ierr/=0) stop 'problem allocating cfga%ispincomb in type_cfg'
    end if

    call MPI_Get_address(cfga%N,              disp(1),  ierr)
    call MPI_Get_address(cfga%ispincomb,      disp(2),  ierr)

    base = disp(1)
    disp = disp - base

    blocklen(1)   = 1
    blocklen(2)   = size(cfga%ispincomb)

    etype(1:N)   = MPI_INTEGER

    call MPI_Type_create_struct(N, blocklen, disp, etype, myMPItype, iout)

  end subroutine myMPI_create_struct_cfga

end module type_cfg
