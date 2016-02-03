module mod_rhoq

#ifdef CPP_MPI
use mpi
#endif

implicit none

type type_rhoq

  ! number of scalars in this type, needed for MPI communication
  integer :: Nscalars = 7
  ! number of arrays in this type, needed for MPI communication
  integer :: Narrays = 11
  
  ! impurity cluster
  integer :: Nscoef ! number of impurities in impurity cluster
  integer :: Nlayer ! number of different layers in impurity cluster
  integer, allocatable :: ilay_scoef(:) ! (Nscoef), layer index for all positions in impurity cluster
  double precision, allocatable :: r_scoef(:,:) ! (3,Nscoef), position vector to all positions in impurity cluster
  
  ! kmesh
  integer :: nkpt ! total number of kpoints
  double precision :: volbz ! Brillouin zone volume -> integration weight
  double precision, allocatable :: kpt(:,:) ! (3,nkpt), coordinates of kpoints
  
  ! geometry etc.
  integer :: mu_0 ! layer index of probe position
  integer :: natyp ! number of atoms in host system
  double precision, allocatable :: r_basis(:,:) ! (3,natyp), real space positions of all atoms in host system
  double precision, allocatable :: L_i(:,:) ! (3,Nscoef+1), lattice vectors for all atoms in impurity cluster and the probe layer
  integer :: lmmaxso ! size in l,m(,s) (if SOC calculation) subblocks
  
  ! Green functions etc.
  logical :: Ghost_k_memsave ! logical switch which determines if Ghost_k is stored in a file or kept in memory
  double complex, allocatable :: Ghost(:,:,:) ! (Nlayer,lmmaxso,lmmaxso), 
  double complex, allocatable :: Ghost_k(:,:,:,:) ! (Nlayer,nkpt,lmmaxso,lmmaxso), 
  double complex, allocatable :: Dt(:,:) ! (Nscoef,Nscoef,lmmaxso,lmmaxso), 
  double complex, allocatable :: Gimp(:,:) ! (Nscoef,Nscoef,lmmaxso,lmmaxso), 
  double complex, allocatable :: tau(:,:) ! (Nscoef,Nscoef,lmmaxso,lmmaxso), 
  
end type type_rhoq


type (type_rhoq), save :: t_rhoq


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef CPP_MPI
subroutine bcast_scalars_rhoq(t_rhoq)
  ! broadcast scalars in t_rhoq over all ranks
  use mpi
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  ! local variables
  integer :: ierr
  integer :: blocklen1(t_rhoq%Nscalars),etype1(t_rhoq%Nscalars),myMPItype1
  integer(kind=MPI_ADDRESS_KIND) :: disp1(t_rhoq%Nscalars), base
  
  call MPI_Get_address(t_rhoq%Nscoef,          disp1(1), ierr)
  call MPI_Get_address(t_rhoq%Nlayer,          disp1(2), ierr)
  call MPI_Get_address(t_rhoq%nkpt,            disp1(3), ierr)
  call MPI_Get_address(t_rhoq%mu_0,            disp1(4), ierr)
  call MPI_Get_address(t_rhoq%natyp,           disp1(5), ierr)
  call MPI_Get_address(t_rhoq%lmmaxso,          disp1(6), ierr)
  call MPI_Get_address(t_rhoq%Ghost_k_memsave, disp1(7), ierr)

  base  = disp1(1)
  disp1 = disp1 - base

  blocklen1(1:7) = 1

  etype1(1:6) = MPI_INTEGER
  etype1(7)   = MPI_LOGICAL

  call MPI_Type_create_struct(t_rhoq%Nscalars, blocklen1, disp1, etype1, myMPItype1, ierr)
  if(ierr/=MPI_SUCCESS) stop '[bcast_scalars_rhoq] Problem in create_mpimask_t_rhoq'
  
  call MPI_Type_commit(myMPItype1, ierr)
  if(ierr/=MPI_SUCCESS) stop '[bcast_scalars_rhoq] error comiting create_mpimask_t_rhoq'
  
  call MPI_Bcast(t_rhoq%Nscalars, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[bcast_scalars_rhoq] error brodcasting scalars in t_rhoq'
  
  call MPI_Type_free(myMPItype1, ierr)

end subroutine bcast_scalars_rhoq
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef CPP_MPI
subroutine bcast_arrays_rhoq(t_rhoq)
  ! broadcast arrays in t_rhoq over all ranks
  ! has to be done after bacst scalars and then init_t_rhoq on all other processes that are not the master
  use mpi
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  ! local variables
  integer :: ierr
  integer :: blocklen1(t_rhoq%Narrays),etype1(t_rhoq%Narrays),myMPItype1
  integer(kind=MPI_ADDRESS_KIND) :: disp1(t_rhoq%Narrays), base

  call MPI_Get_address(t_rhoq%ilay_scoef, disp1(1), ierr)
  call MPI_Get_address(t_rhoq%r_scoef,    disp1(2), ierr)
  call MPI_Get_address(t_rhoq%volbz,      disp1(3), ierr)
  call MPI_Get_address(t_rhoq%kpt,        disp1(4), ierr)
  call MPI_Get_address(t_rhoq%r_basis,    disp1(5), ierr)
  call MPI_Get_address(t_rhoq%L_i,        disp1(6), ierr)
  call MPI_Get_address(t_rhoq%Ghost,      disp1(7), ierr)
  call MPI_Get_address(t_rhoq%Ghost_k,    disp1(8), ierr)
  call MPI_Get_address(t_rhoq%Dt,         disp1(9), ierr)
  call MPI_Get_address(t_rhoq%Gimp,      disp1(10), ierr)
  call MPI_Get_address(t_rhoq%tau,       disp1(11), ierr)

  base  = disp1(1)
  disp1 = disp1 - base

  blocklen1(1)  = t_rhoq%Nscoef
  blocklen1(2)  = 3*t_rhoq%Nscoef
  blocklen1(3)  = 1
  blocklen1(4)  = 3*t_rhoq%nkpt
  blocklen1(5)  = 3*t_rhoq%natyp
  blocklen1(6)  = 3*(t_rhoq%Nscoef+1)
  blocklen1(7)  = t_rhoq%Nlayer*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(8)  = t_rhoq%Nlayer*t_rhoq%nkpt*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(9)  = t_rhoq%Nscoef*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(10) = t_rhoq%Nscoef*t_rhoq%Nscoef*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(11) = t_rhoq%Nscoef*t_rhoq%Nscoef*t_rhoq%lmmaxso*t_rhoq%lmmaxso

  etype1(1)    = MPI_INTEGER
  etype1(2:6)  = MPI_DOUBLE_PRECISION
  etype1(7:11) = MPI_DOUBLE_COMPLEX

  call MPI_Type_create_struct(t_rhoq%Narrays, blocklen1, disp1, etype1, myMPItype1, ierr)
  if(ierr/=MPI_SUCCESS) stop '[bcast_arrays_rhoq] Problem in create_mpimask_t_rhoq'
  
  call MPI_Type_commit(myMPItype1, ierr)
  if(ierr/=MPI_SUCCESS) stop '[bcast_arrays_rhoq] error comiting create_mpimask_t_rhoq'
  
  call MPI_Bcast(t_rhoq%Narrays, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[bcast_arrays_rhoq] error brodcasting scalars in t_rhoq'
  
  call MPI_Type_free(myMPItype1, ierr)


end subroutine bcast_arrays_rhoq
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine init_t_rhoq(t_rhoq)
  ! initialize allocatable arrays in t_rhoq
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer :: ierr, N
  
  N = t_rhoq%Nscoef*t_rhoq%lmmaxso

  ierr = 0
  if(.not.allocated(t_rhoq%ilay_scoef)) allocate(t_rhoq%ilay_scoef(t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating ilay_scoef in rhoq'
  if(.not.allocated(t_rhoq%r_scoef)) allocate(t_rhoq%r_scoef(3,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating r_scoef in rhoq'
  if(.not.allocated(t_rhoq%kpt)) allocate(t_rhoq%kpt(3,t_rhoq%nkpt), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating kpt in rhoq'
  if(.not.allocated(t_rhoq%r_basis)) allocate(t_rhoq%r_basis(3,t_rhoq%natyp), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating r_basis in rhoq'
  if(.not.allocated(t_rhoq%L_i)) allocate(t_rhoq%L_i(3,t_rhoq%Nscoef+1), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Nscoef in rhoq'
  if(.not.allocated(t_rhoq%Ghost)) allocate(t_rhoq%Ghost(t_rhoq%lmmaxso,t_rhoq%lmmaxso,t_rhoq%Nlayer), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Nscoef in rhoq'
  if(.not.allocated(t_rhoq%Ghost_k)) allocate(t_rhoq%Ghost_k(t_rhoq%lmmaxso,t_rhoq%lmmaxso,t_rhoq%Nlayer,t_rhoq%nkpt), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Nscoef in rhoq'
  if(.not.allocated(t_rhoq%Dt)) allocate(t_rhoq%Dt(N, N), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Nscoef in rhoq'
  if(.not.allocated(t_rhoq%Gimp)) allocate(t_rhoq%Gimp(N, N), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Nscoef in rhoq'
  if(.not.allocated(t_rhoq%tau)) allocate(t_rhoq%tau(N, N) , stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Nscoef in rhoq'

end subroutine init_t_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine read_scoef_rhoq(t_rhoq)
  ! read in impurity cluster information and save this in t_rhoq
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  !local variables
  integer natomimp, iatom, ierr, Nlayer, tmp
  double precision, allocatable :: ratomimp(:,:)
  integer, allocatable :: atomimp(:)

  open(unit=32452345,file='scoef',iostat=ierr)
  if (ierr/=0) stop '[read_scoef_rhoq] file not found'
  read(32452345,*) natomimp
  write(*,*) 'natomimp',natomimp
  allocate(ratomimp(3,natomimp))
  allocate(atomimp(natomimp))
  do iatom=1,natomimp
    read(32452345,*) ratomimp(:,iatom),atomimp(iatom)
    write(*,'(A,I,A,3F)') 'IMPATOM ',iatom,' :',ratomimp(:,iatom)
  end do
  close(32452345)
  
  ! save natomimp in t_rhoq
  t_rhoq%Nscoef     = natomimp
  
  ! allocate ilay_scoef and r_scoef arrays from t_rhoq and save information in there
  ierr = 0
  if(.not.allocated(t_rhoq%ilay_scoef)) allocate(t_rhoq%ilay_scoef(t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[read_scoef_rhoq] error allocating ilay_scoef in read_scoef_rhoq'
  if(.not.allocated(t_rhoq%r_scoef)) allocate(t_rhoq%r_scoef(3,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[read_scoef_rhoq] error allocating r_scoef in read_scoef_rhoq'
  
  t_rhoq%ilay_scoef = atomimp
  t_rhoq%r_scoef    = ratomimp
  
  ! find parameter Nlayer
  Nlayer = 0
  tmp = 0
  do iatom=1,natomimp
    if(tmp/=atomimp(iatom)) then 
      tmp = atomimp(iatom)
      Nlayer = Nlayer + 1
    end if
  end do
  
  t_rhoq%Nlayer = Nlayer
  
  ! deallocate unused arrays
  deallocate(ratomimp, atomimp)
  

end subroutine read_scoef_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine read_input_rhoq(t_rhoq, uio)
  ! read input such as proble layer mu_0, and parameters such as lmmaxso etc.
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  character*256, intent(in) :: uio ! unit in which to find inputcard
  ! local
  integer :: ier
  
!   open(32452345, file='inputcard', status='old', iostat=ier)
!   if(ier/=0) stop '[read_input_rhoq] error inputcard not found'
  
  call ioinput('mu0_rhoq        ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) t_rhoq%mu_0
    write(*,*) 'Proble layer mu_0= ',t_rhoq%mu_0
  else
    stop '[read_input_rhoq] error "mu0_rhoq" not found in inputcard'
  endif

  call ioinput('NATYP           ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) t_rhoq%natyp
    write(*,*) 'number of layers in system (natyp)= ',t_rhoq%natyp
  else
    stop '[read_input_rhoq] error "NATYP" not found in inputcard'
  endif

  call ioinput('Ghost_k_memsave ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) t_rhoq%Ghost_k_memsave
    write(*,*) 'save Ghost in memory? ',t_rhoq%Ghost_k_memsave
  else
    stop '[read_input_rhoq] error "Ghost_k_memsave" not found in inputcard'
  endif
  
!   close(32452345)

end subroutine read_input_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine save_kmesh_rhoq(t_rhoq,nkpt,kpt,volbz)
  ! save kmesh information in t_rhoq
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: nkpt ! total number of kpoints
  double precision, intent(in) :: volbz ! Brillouin zone volume -> integration weight
  double precision, intent(in) :: kpt(3,nkpt) ! coordinates in reciprocal space of the kpoints
  !local
  integer :: ierr
  
  t_rhoq%nkpt = nkpt
  ierr = 0
  t_rhoq%volbz = volbz
  if(.not.allocated(t_rhoq%kpt)) allocate(t_rhoq%kpt(3,nkpt), stat=ierr)
  if(ierr/=0) stop '[save_kmesh_rhoq] error allocating kpt in rhoq'
  t_rhoq%kpt = kpt
  
end subroutine save_kmesh_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine save_geometry_rhoq(t_rhoq,r_basis,lmmaxso,natyp)
  ! save geometry information
  ! Attention: t_rhoq has to contain mu_0 etc from inputcard (done with read_input_rhoq)
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  double precision, intent(in)   :: r_basis(3,natyp) ! real space positions of all atoms in host system
  integer, intent(in)            :: natyp, lmmaxso ! number of atoms in impcluster, size in l,m(,s) (if SOC calculation) subblocks
  !local
  integer :: ierr
  
  ! save lmmaxso (needed later)
  t_rhoq%lmmaxso = lmmaxso

  ! allocate and save r_basis
  ierr = 0
  if(.not.allocated(t_rhoq%r_basis)) allocate(t_rhoq%r_basis(3,t_rhoq%natyp), stat=ierr)
  if(ierr/=0) stop '[save_geometry_rhoq] error allocating r_basis in rhoq'
  t_rhoq%r_basis = r_basis
  
  ! calculate lattice vectors needed later in calc_rhoq
  ! x = R^mu_n + r with R^mu_n = X_mu + L_i where L_i is the lattice vector that appears in the Fourier transform
  call get_L_vecs_rhoq(t_rhoq)

end subroutine save_geometry_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_L_vecs_rhoq(t_rhoq)
  ! calculate L_i vectors needed in calculate_rhoq from information of impurity cluster etc.
  ! L_i = R_mu - R_cls - r^cls_i, here R_mu and R_cls are the positions in the host systems basis and r^cls_i is the position inside the impurity cluster
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  ! local
  integer :: mu_0, mu_cls, mu_orig ! layer indices for probe position, cluster center and origin position, respectively
  double precision :: R_mu(3), R_cls(3) ! position in host system according to mu_orig and mu_cls
  double precision, allocatable :: rcls_i(:,:) ! size=(3,Nscoef), positions inside the impurity cluster
  double precision, allocatable :: r_basis(:,:) ! size=(3,natyp), positions in host system
  double precision, allocatable :: L_i(:,:) ! size=(3,Nscoef), lattice vectors
  integer :: ilayer, irun ! loop parameter
  double precision, parameter :: eps = 1.0D-6 ! small epsilon -> threashold
  integer :: ierr
  
  ! allocations etc.
  allocate(rcls_i(3,t_rhoq%Nscoef), r_basis(3,t_rhoq%natyp), L_i(3,t_rhoq%Nscoef))
  r_basis = t_rhoq%r_basis
  rcls_i = t_rhoq%r_scoef
  
  write(*,*) 'rbasis',r_basis
  write(*,*) 'rcls_i', rcls_i
    
  ! find mu_orig
  ilayer = 0
  irun = 1
  do while(irun==1)
    ilayer = ilayer+1
    if(dsqrt(sum(r_basis(:,ilayer)**2))<eps) irun = 0
  end do
  mu_orig = ilayer
  write(*,*) 'mu_orig', mu_orig
  
  ! find mu_cls
  ilayer = 0
  irun = 1
  do while(irun==1)
    ilayer = ilayer+1
    if(dsqrt(sum(rcls_i(:,ilayer)**2))<eps) irun = 0
  end do
  mu_cls = t_rhoq%ilay_scoef(ilayer)
  write(*,*) 'mu_cls', mu_cls
  
  ! set R_mu and R_cls from mu_orig and mu_cls
  R_mu(:) = r_basis(:,mu_orig)
  R_cls(:) = r_basis(:,mu_cls)
  
  ! L_i = R_mu - R_cls - r^cls_i
  do ilayer=1,t_rhoq%Nscoef
    L_i(:,ilayer) = R_mu(:) - R_cls(:) - rcls_i(:,ilayer)
    write(*,*) 'L_i', L_I(:,ilayer)
  end do
  
  ! save result and exit
  ierr = 0
  if(.not.allocated(t_rhoq%L_i)) allocate(t_rhoq%L_i(3,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[get_L_vecs_rhoq] error allocating L_i in rhoq'
  t_rhoq%L_i(:,:) = L_i(:,:)
  deallocate(rcls_i, r_basis, L_i)

end subroutine get_L_vecs_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine save_Ghost_rhoq(t_rhoq, Ghost, lmmaxso, Nlayer)
  ! save structural GF  in t_rhoq
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: Nlayer, lmmaxso
  double complex, intent(in) :: Ghost(lmmaxso,lmmaxso,Nlayer)
  
  ! some checks
  if( .not.allocated(t_rhoq%Ghost) ) stop '[save_Ghost_rhoq] error: Ghost not allocated in rhoq'
  if( any( shape(t_rhoq%Ghost)/=shape(Ghost) ) ) stop '[save_Ghost_rhoq] error: input Ghost and allocated Ghost in rhoq do not match in shape'
  
  t_rhoq%Ghost(:,:,:) = Ghost(:,:,:)
  
end subroutine save_Ghost_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_Ghost_k_rhoq(t_rhoq,alat,cls,nacls,rr,ezoa,atom,rcls,lmax,naez,ncls,nr,nemb,naclsmax)
  ! calculate Ghost(k;E) from Ghost(E) and kmesh information by fourier transform
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  ! parameters
  integer, intent(in) :: lmax,naez,ncls,nr,nemb
  ! scalars
  double precision, intent(in) :: alat, naclsmax
  ! arrays
  double precision, intent(in) :: rcls(3,ncls,ncls),rr(3,0:nr)
  integer, intent(in)          :: atom(ncls,naez+nemb),cls(naez+nemb),ezoa(ncls,naez+nemb),nacls(ncls)
  ! local
  double complex, allocatable :: Ghost_k_tmp(:,:,:)
  integer :: ierr, ik


  ! allocations etc.
  if(t_rhoq%Ghost_k_memsave) then
    ierr = 0
    if(.not.allocated(t_rhoq%Ghost_k)) allocate(t_rhoq%Ghost_k(t_rhoq%lmmaxso,t_rhoq%lmmaxso,t_rhoq%Nlayer,t_rhoq%nkpt), stat=ierr)
    if(ierr/=0) stop '[calc_Ghost_k_rhoq] error in allocating Ghost_k in rhoq'
  end if ! t_rhoq%Ghost_k_memsave
  allocate(Ghost_k_tmp(t_rhoq%lmmaxso,t_rhoq%lmmaxso,t_rhoq%Nlayer), stat=ierr)
  if(ierr/=0) stop '[calc_Ghost_k_rhoq] error in allocating Ghost_k_tmp'
  
  ! do fourier transform of Ghost
  do ik=1,t_rhoq%nkpt
    call dlke0(Ghost_k_tmp,alat,t_rhoq%natyp,cls,nacls,naclsmax,rr,ezoa,atom,t_rhoq%kpt,rcls,t_rhoq%Ghost)
    if(t_rhoq%Ghost_k_memsave) t_rhoq%Ghost_k(:,:,:,ik) = Ghost_k_tmp(:,:,:)
  end do
  
  
  deallocate(Ghost_k_tmp)

end subroutine calc_Ghost_k_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine read_Dt_Gimp_rhoq(t_rhoq, lmmaxso, ncls)
  ! read in DTMTRX and GMATLL_GES, provided by zulapi code
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: lmmaxso, ncls

  ! read in 
  call read_DTMTRX( t_rhoq%Dt, lmmaxso, ncls)
  call read_green_ll(ncls, lmmaxso, t_rhoq%Gimp)

end subroutine read_Dt_Gimp_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! helping functions from Bernds FS code, in mod_scattering
! needed to read in DTMTRX and GMATLL_GES files

subroutine read_DTMTRX( tmat, lmmaxso, ncls)
  implicit none
  ! .... Arguments ....
  integer,           intent(in)  :: lmmaxso, ncls
  double complex,    intent(out) :: tmat(ncls*lmmaxso,ncls*lmmaxso) ! tmat file actually too large? not vector but diagonal matrix is stored???
  ! local
  integer :: ierr, icls, icls2, lm1, lm2, idummy1, idummy2, nClustertest
  double precision :: Rclstest(3), Zatom, Rdist
  double complex :: deltamat_dummy
  !parameter
  double precision, parameter :: eps = 1.0D-6 ! small epsilon -> threshold


  !============ BEGIN ===========!
  !=== read the file 'DTMTRX' ===!
  open(unit=562, file='DTMTRX', form='formatted', action='read')
  write(*,*) 'reading file DTMTRX'

  !read the number of atoms in the impurity-cluster
  read(562, fmt=*) nClustertest
  if( nClustertest/=t_rhoq%Nscoef ) stop 'Error: cluster information from DTMTRX and scoef file do not match, different number of atoms!'
  
  !read in the position of the atoms in the impurity cluster and compare to scoef file info
  do icls=1,nClustertest
    read(562, fmt=*) RClstest(:)
    if( any( abs(Rclstest(:)-t_rhoq%r_scoef(:,icls))>eps ) ) stop 'Error: cluster information from DTMTRX and scoef file do not match!'
  end do!icls

  !read in the t-matrix
  do lm1=1,ncls*lmmaxso
    do lm2=1,ncls*lmmaxso
!       write(*,*) lm1,lm2, ncls, lmmaxso, ncls*lmmaxso
      read(562,"(2I5,4e17.9)") idummy1, idummy2, tmat(lm2,lm1), deltamat_dummy
!       write(*,*) idummy1, idummy2
    end do!lm2
  end do!lm1
  
  close(562)
  write(*,*) 'done reading DTMTRX'
  !=== read the file 'DTMTRX' ===!
  !============  END  ===========!

end subroutine read_DTMTRX


subroutine read_green_ll(ncls,lmmaxso, Gll0)

  implicit none
  integer,          intent(in)  :: lmmaxso, ncls
  double complex, intent(out) :: Gll0(ncls*lmmaxso, ncls*lmmaxso)
  ! local
  double precision :: dEimag
  double complex :: energy(3), dE1, dE2, ctmp1, ctmp2
  double complex :: Gll_3(ncls*lmmaxso, ncls*lmmaxso,3), &
                  &  dGll(ncls*lmmaxso, ncls*lmmaxso),    &
                  & d2Gll(ncls*lmmaxso, ncls*lmmaxso)
  integer :: ienergy, lm1, lm2, id1, id2, ierr
  double complex, parameter :: CI=(0d0,1d0)

  ! Read the gll_3 (for the three energies)
  open(unit=1283, file='GMATLL_GES', form='formatted', action='read')
  write(*,*) 'reading file GMATLL_GES'
  
  do ienergy=1,3
  
    write(*,*) ' reading energy',ienergy
    read(1283,"(2(e17.9,X))") energy(ienergy)

    do lm1=1,ncls*lmmaxso
      do lm2=1,ncls*lmmaxso
        read(1283,"((2I5),(2e17.9))") id1, id2, Gll_3(lm2, lm1, ienergy)
      end do!icls2
    end do!icls1

  end do!ienergy

  ! Checks
  dE1 = energy(2)-energy(1)
  dE2 = energy(3)-energy(2)
  if(abs(dE1-dE2)>1d-8) stop '3 Energy points not equidistant'

  ! Construct first and second derivative of G(E)
  ctmp1 = 0.5d0/dE1
  ctmp2 = 1d0/dE1**2
  dGll  = ctmp1*( Gll_3(:,:,3)                    - Gll_3(:,:,1) )
  d2Gll = ctmp2*( GLL_3(:,:,3) - 2d0*Gll_3(:,:,2) + Gll_3(:,:,1) )

  ! extrapolate to energy with imag(E)=0
  dEimag = aimag(energy(2))
  Gll0 = Gll_3(:,:,2) - CI*dEimag*dGll(:,:) - 0.5d0*dEimag**2 * d2Gll(:,:)

end subroutine read_green_ll


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_tau_rhoq(t_rhoq)
  ! calculate tau from Dt and Gimp according to
  ! tau_i,i' = Dt_i delta_i,i' + Dt_i Gimp_i,i' Dt_i'
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  ! locals
  double complex, allocatable :: temp(:,:)
  integer :: ierr, N
  double complex, parameter :: C0=(0.0D0,0.0D0), C1=(1.0D0,0.0D0)

  ! matrix sizes
  N = t_rhoq%Nscoef*t_rhoq%lmmaxso
  
  ! now do matrix matrix operations:
  ! calculate first Gimp.Dt, here Dt is already diagonal matrix: temp = Gimp.Dt
  allocate(temp(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_tau_rhoq] error allocating temp'
  ! temp = Gimp.Dt
  call ZGEMM('n','n',N,N,N,C1,t_rhoq%Gimp,N,t_rhoq%Dt,N,C0,temp,N)
  
  ierr = 0
  if(.not.allocated(t_rhoq%tau)) allocate(t_rhoq%tau(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_tau_rhoq] error allocating tau in rhoq'

  ! tau = Dt + Dt.Gimp.Dt
  !     = tau + Dt.temp
  t_rhoq%tau = t_rhoq%Dt
  call ZGEMM('n','n',N,N,N,C1,t_rhoq%Dt,N,temp,N,C1,t_rhoq%tau,N)
  
  
  deallocate(temp)

end subroutine calc_tau_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_Q_mu_rhoq(lmax, ntotd, npan_tot, ncheb, npan_log, npan_eq, &
  & nsra, lmmaxso, Rll, Rllleft, rpan_intervall, ipan_intervall, ncleb,            &
  & lm2d, lmaxd, iend, cleb, icleb, loflm, nfund, lmxsp, ifunm,    &
  & nrmaxd, thetasnew, irmdnew, q_mu )
  ! calculate prefactor needed in rhoq calculation:
  ! Q^mu_Lambda,Lamabda' = Int d^3r Tr[ R^mu_Lmabda(\vec{r}) R^left,mu_Lmabda'(\vec{r}) ] where R^left denotes the left solution
  implicit none
  integer, intent(in) :: lmax, ntotd, npan_tot, ncheb, npan_log, npan_eq
  integer, intent(in) :: nsra, lmmaxso, irmdnew ! left/right solutions (wave functions)
  double complex, intent(in) :: Rll(nsra*lmmaxso,lmmaxso,irmdnew), Rllleft(nsra*lmmaxso,lmmaxso,irmdnew) ! wave functions in new mesh
  double precision, intent(in) :: rpan_intervall(0:ntotd)
  integer, intent(in) :: ipan_intervall(0:ntotd)
  integer, intent(in) :: ncleb, lm2d, lmaxd, iend ! Gaunt coefficients
  double precision, intent(in) :: cleb(ncleb) ! Gaunt coefficients
  integer, intent(in) :: icleb(ncleb,4), loflm(lm2d) ! Gaunt coefficients
  integer, intent(in) :: nfund, lmxsp, ifunm(lmxsp), nrmaxd ! shape functions
  double precision, intent(in) :: thetasnew(nrmaxd,nfund) ! shape functions in Chebychev mesh (new mesh)
  double complex, intent(out) :: q_mu(lmmaxso,lmmaxso)
  
  double precision, parameter :: c0ll = 1.0d0/sqrt(16.0d0*datan(1.0d0)) ! Y_00 = 1/sqrt(4*pi) needed for theta_00 * Y_00 = 1, other spherical harmonics are in Gaunt coefficients cleb(:,:)
  
  ! local
  integer :: imt1, lm1, lm2, lm3, ir, ifun, j, lm01, lm02
  double complex :: q_of_r(nrmaxd,lmmaxso,lmmaxso)
  

  ! get indices of imt and irmd in new mesh -> boundaries for ir loops
  imt1=ipan_intervall(npan_log+npan_eq)+1

  ! first treat spherical block (diagonal in lm, lm' -> only lm'==lm)
  do lm01 = 1,lmmaxso
    ! fill diagonal
    do ir = 1,irmdnew
      q_of_r(ir, lm01, lm01) = q_of_r(ir, lm01, lm01) + rll(lm01, lm01, ir) * rllleft(lm01, lm01, ir) ! note: diagonal in lm01!!!
    enddo ! ir
    ! add shapefunction to points with r> R_MT
    do ir=imt1+1,irmdnew
      q_of_r(ir, lm01, lm01) = q_of_r(ir, lm01, lm01)*thetasnew(ir,1)*c0ll ! add shapefunction on diagonal for last points outside of MT radius
    enddo ! ir
  enddo ! lm01

!   write(*,*) 'iend',iend
!   write(*,*) 'cleb',cleb
!   write(*,*) 'icleb',icleb
  
  ! treat non-spherical components -> off-diagonal blocks in lm, lm' matrices
  do lm01 = 1, lmmaxso ! loop ofer outer lm-component
    do lm02 = 1, lmmaxso ! loop ofer outer lm-component
      do j = 1,iend ! loop over number of non-vanishing Gaunt coefficients
        ! set lm indices for Gaunt coefficients
        lm1 = icleb(j,1)
        lm2 = icleb(j,2)
        lm3 = icleb(j,3)
!         write(*,*) lm1, lm2, lm3, j, lm01, lm02
        ! get shape function index for lm3
        ifun = ifunm(lm3)
        write(*,*) j, ifun, lm3
        ! do loop over non-spherical points here
        do ir=imt1+1,irmdnew
          q_of_r(ir, lm02, lm01) = q_of_r(ir, lm02, lm01) + cleb(j) * rll(lm01, lm1, ir) * rllleft(lm02, lm2, ir) * thetasnew(ir,ifun) ! ifun gives lm3 component
        enddo ! ir
      enddo ! j -> sum of Gaunt coefficients
      
      ! do radial integration
      call intcheb_cell(q_of_r(:,lm02,lm01),q_mu(lm02,lm01),rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
    end do ! lm02
  end do ! lm01

end subroutine calc_Q_mu_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_rhoq(t_rhoq)
  ! calculate Delta rho^mu(q) = Int{ Delta rho(q, X_mu+r), dr }
  !                           = -1/(2*pi*i) Tr[ Q^mu Int{ Ghost(k) tau Ghost(k+q),dk } - Q^mu,* Int{ Ghost^*(k) tau^* Ghost^*(k-q),dk }]
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  
  
  
  
end subroutine calc_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module mod_rhoq

! for testing only, comment out when ready
!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!
program test

  use mod_rhoq
  
  character*256 :: uio ! unit in which to find inputcard
  integer :: nkpt ! total number of kpoints
  integer :: lmmaxso, natyp, nfund, lmxsp, nrmaxd
  integer :: naez,ncls,nr,nemb
  integer :: lmax, ntotd, npan_tot, ncheb, npan_log, npan_eq
  integer :: nsra, irmdnew ! left/right solutions (wave functions)
  integer :: ncleb, lm2d, lmaxd, iend ! Gaunt coefficients
  double precision :: alat, naclsmax
  double precision :: volbz ! Brillouin zone volume -> integration weight
  
  double precision, allocatable :: kpt(:,:) ! coordinates in reciprocal space of the kpoints
  double precision, allocatable :: r_basis(:,:) ! real space positions of all atoms in host system
  double complex, allocatable   :: Ghost(:,:,:)
  double precision, allocatable :: rcls(:,:,:),rr(:,:)
  integer, allocatable          :: atom(:,:),cls(:),ezoa(:,:),nacls(:)
  double complex, allocatable   :: Rll(:,:,:), Rllleft(:,:,:) ! wave functions in new mesh
  double precision, allocatable :: rpan_intervall(:)
  integer, allocatable          :: ipan_intervall(:)
  double precision, allocatable :: cleb(:) ! Gaunt coefficients
  integer, allocatable          :: icleb(:,:), loflm(:) ! Gaunt coefficients
  integer, allocatable          :: ifunm(:) ! shape functions
  double precision, allocatable :: thetasnew(:,:) ! shape functions in Chebychev mesh (new mesh)
  double complex, allocatable   :: q_mu(:,:)
  
  integer :: lm1, lm2, N, i, j
  
  ! read in scalars
  uio = 'inputcard'
  open(9999, file='params.txt')
  read(9999,*) lmmaxso, natyp
  lmmaxso = lmmaxso/2
  read(9999,*) naez, ncls, nr, nemb, lmax
  read(9999,*) alat, naclsmax
  close(9999)

  write(*,*) 'params', lmmaxso, natyp, naez, ncls, nr, nemb, lmax, alat, naclsmax
  
  open(9999, file='kpoints.txt')
  read(9999,*) nkpt
  allocate( kpt(3,nkpt) )
  read(9999,*) volbz, kpt(1:3,1:nkpt)
  close(9999)

  write(*,*) 'kpoints', nkpt, volbz
  
  open(9999, file='host.txt')
  allocate( r_basis(3,natyp), Ghost(lmmaxso,lmmaxso,t_rhoq%Nlayer) )
  allocate( rcls(3,ncls,ncls), rr(3,0:nr) )
  allocate( atom(ncls,naez+nemb), cls(naez+nemb), ezoa(ncls,naez+nemb), nacls(ncls) )
  read(9999,*) r_basis(1:3,1:natyp)!, Ghost(1:lmmaxso,1:lmmaxso,)
  read(9999,*) rcls(1:3,1:ncls,1:ncls), rr(1:3,0:nr), atom(1:ncls,1:naez+nemb)
  read(9999,*) cls(1:naez+nemb), ezoa(1:ncls,1:naez+nemb), nacls(1:ncls)
  close(9999)
  
           open(9999, file='wavefunctions.txt')
           read(9999,'(100I9)') ntotd, npan_tot, ncheb, npan_log, npan_eq, nsra, irmdnew

           allocate( Rll(1:nsra*lmmaxso,1:lmmaxso,1:irmdnew),          &
     &               Rllleft(1:nsra*lmmaxso,1:lmmaxso,1:irmdnew) ,     &
     &               rpan_intervall(0:ntotd),ipan_intervall(0:ntotd) )
           do lm1=1,irmdnew
             read(9999,'(20000E16.7)') Rll(1:nsra*lmmaxso,1:lmmaxso,lm1)      
             read(9999,'(20000E16.7)') Rllleft(1:nsra*lmmaxso,1:lmmaxso,lm1)
           enddo
           do lm1=0,ntotd
             read(9999,'(E16.7,I9)') rpan_intervall(lm1), ipan_intervall(lm1)
           enddo
           close(9999)
  
  
  open(9999, file='cleb_shapefun.txt')
  read(9999,*) ncleb, lm2d, lmaxd, iend, nfund, lmxsp, nrmaxd
  allocate( icleb(ncleb,4), loflm(lm2d), cleb(ncleb) )
  allocate( ifunm(lmxsp) , thetasnew(nrmaxd,nfund), q_mu(lmmaxso,lmmaxso) )
  do lm1=1,ncleb
    read(9999,'(E16.7,4I9)') cleb(lm1), icleb(lm1,1:4)
  end do
  do lm1=1,lm2d
    read(9999,'(I9)') loflm(lm1)
    write(*,*) 'loflm',loflm(lm1)
  end do
  do lm1=1,lmxsp
    read(9999,'(I9)') ifunm(lm1)
    write(*,*) 'ifunm',ifunm(lm1)
  end do
  do lm1=1,nrmaxd
    do lm2=1,nfund
      read(9999,'(E16.7)') thetasnew(lm1,lm2)
    end do
  end do
  close(9999)
  

  call read_scoef_rhoq(t_rhoq)

  call read_input_rhoq(t_rhoq, uio)

  call save_kmesh_rhoq(t_rhoq,nkpt,kpt,volbz)

  call save_geometry_rhoq(t_rhoq,r_basis,lmmaxso,natyp)
  
!   call save_Ghost_rhoq(t_rhoq, Ghost, lmmaxso, t_rhoq%Nlayer)
! 
! #ifdef CPP_MPI
!   call bcast_scalars_rhoq(t_rhoq)
! #endif
  call init_t_rhoq(t_rhoq)
! #ifdef CPP_MPI
!   call bcast_arrays_rhoq(t_rhoq)
! #endif
! 
!   call calc_Ghost_k_rhoq(t_rhoq,alat,cls,nacls,rr,ezoa,atom,rcls,lmax,naez,ncls,nr,nemb,naclsmax)
! 

!   write(*,*) 'before Gimp,Dt',t_rhoq%Nscoef,allocated(t_rhoq%Dt), allocated(t_rhoq%Gimp)
!   write(*,*) 'shape Gimp,Dt',shape(t_rhoq%Dt), shape(t_rhoq%Gimp)
  call read_Dt_Gimp_rhoq(t_rhoq, lmmaxso, t_rhoq%Nscoef)
  
  call calc_tau_rhoq(t_rhoq)
  
!   N = t_rhoq%Nscoef*t_rhoq%lmmaxso
!   do i=1,N
!     do j=1,N
!       write(425,'(2E16.7)') t_rhoq%Dt(i,j)
!     end do
!   end do
!   do i=1,N
!     do j=1,N
!       write(426,'(2E16.7)') t_rhoq%Gimp(i,j)
!     end do
!   end do
!   do i=1,N
!     do j=1,N
!       write(427,'(2E16.7)') t_rhoq%tau(i,j)
!     end do
!   end do
      
      
!   write(*,*) 'calc Q'
!     write(*,*) icleb
!   call calc_Q_mu_rhoq(lmax, ntotd, npan_tot, ncheb, npan_log, npan_eq, &
!        & nsra, lmmaxso, Rll, Rllleft, rpan_intervall, ipan_intervall, ncleb,   &
!        & lm2d, lmaxd, iend, cleb, icleb, loflm, nfund, lmxsp, ifunm,           &
!        & nrmaxd, thetasnew, irmdnew, q_mu )
!     
!   call calc_rhoq(t_rhoq) 

end program test
!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!
