module mod_ioformat

  implicit none

    integer, parameter :: MODE_INT=1, MODE_VIS=2

    character(len=*), parameter :: ext_formatted   = '.txt'
    character(len=*), parameter :: ext_unformatted = '.dat'
    character(len=*), parameter :: ext_mpiio       = '.mpidat'
    character(len=*), parameter :: ext_vtkxml      = '.vtp'
    character(len=*), parameter :: ext_new         = '.new'
    character(len=*), parameter :: ext_refined     = '.refined'
    character(len=*), parameter :: ext_orig        = '.orig'

    character(len=*), parameter :: fmt_fn_ext='(A,A)'
    character(len=*), parameter :: fmt_fn_sub_ext='(A,A,A)'
    character(len=*), parameter :: fmt_fn_rank_ext='(A,"_",I5.5,A)'

    character(len=*), parameter :: filemode_vis = '.vis'
    character(len=*), parameter :: filemode_int = '.int'
    character(len=*), parameter :: filemode_rot = '.rot'
    character(len=*), parameter :: filemode_ref = '.ref'

    character(len=*), parameter :: filename_outinfo        = 'outfile'
    character(len=*), parameter :: filename_tbkkrparams    = 'TBkkr_params'
    character(len=*), parameter :: filename_tbkkrcontainer = 'TBkkr_container'
    character(len=*), parameter :: filename_tbkkrrhod      = 'TBkkr_rhod'
    character(len=*), parameter :: filename_cubesinfo      = 'cubesinfo'
    character(len=*), parameter :: filename_cubesrefine    = 'cubesrefine'
    character(len=*), parameter :: filename_fsdata         = 'FSpoints'
    character(len=*), parameter :: filename_fvel           = 'fermivelocity'
    character(len=*), parameter :: filename_spin           = 'spinvalue'
    character(len=*), parameter :: filename_eigvect        = 'eigenvectors'
    character(len=*), parameter :: filename_weights        = 'weights'
    character(len=*), parameter :: filename_scattAmat      = 'AMAT'

    character(len=*), parameter :: filename_vtktest  = 'fstest'
    character(len=*), parameter :: filename_vtkonfs  = 'visdata'
    character(len=*), parameter :: filename_scattfixk= 'scattfix'
    character(len=*), parameter :: filename_lifetime = 'lifetime'
    character(len=*), parameter :: filename_scattmat = 'scatteringmatrix'

! old files - delete when coding done
    character(len=*), parameter :: filename_cubes_IBZ     = 'cubes_irrbz'
    character(len=*), parameter :: filename_kpoints       = 'FSpts_points'
    character(len=*), parameter :: filename_connectivity  = 'FSpts_connect'
    character(len=*), parameter :: filename_realpart      = 'realpart'
    character(len=*), parameter :: filename_spinvalue     = 'spinvalue'
    character(len=*), parameter :: filename_spinvalue_alm = 'spinvalue_alm'
    character(len=*), parameter :: filename_fermivelocity = 'fermivelocity'
    character(len=*), parameter :: filename_eigv_rot      = 'eigenvectors_rotated'
    character(len=*), parameter :: filename_kpoints_triangles = 'kpoints_triangles'
    character(len=*), parameter :: filename_scattering_fix = 'scattering_fix'
    character(len=*), parameter :: filename_invlifetime    = 'invlifetime'

    character(len=*), parameter :: filename_vtkinpuc     = 'VTK_puc'
    character(len=*), parameter :: filename_vtkinpuc_rot = 'VTK_puc_rot'
    character(len=*), parameter :: filename_vtkinBZ      = 'VTK_BZ'
! old files - delete when coding done


end module mod_ioformat
