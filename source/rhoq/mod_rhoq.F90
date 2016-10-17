module mod_rhoq

#ifdef CPP_MPI
use mpi
#endif

implicit none

type type_rhoq

  ! number of scalars in this type, needed for MPI communication
  integer :: Nscalars = 7
  ! number of arrays in this type, needed for MPI communication
  integer :: Narrays = 12
  
  ! impurity cluster
  integer :: Nscoef ! number of impurities in impurity cluster
  integer :: Nlayer ! number of different layers in impurity cluster
  integer, allocatable :: ilay_scoef(:) ! (Nscoef), layer index for all positions in impurity cluster
  double precision, allocatable :: r_scoef(:,:) ! (3,Nscoef), position vector to all positions in impurity cluster
  
  ! kmesh
  integer :: nkpt ! total number of kpoints
  double precision :: volbz ! Brillouin zone volume -> integration weight
  double precision, allocatable :: kpt(:,:) ! (3,nkpt), coordinates of kpoints
  double precision, allocatable :: volcub(:) ! (nkpt), volume of cube associated with kpoint
  
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
  call MPI_Get_address(t_rhoq%volcub,     disp1(5), ierr)
  call MPI_Get_address(t_rhoq%r_basis,    disp1(6), ierr)
  call MPI_Get_address(t_rhoq%L_i,        disp1(7), ierr)
  call MPI_Get_address(t_rhoq%Ghost,      disp1(8), ierr)
  call MPI_Get_address(t_rhoq%Ghost_k,    disp1(9), ierr)
  call MPI_Get_address(t_rhoq%Dt,        disp1(10), ierr)
  call MPI_Get_address(t_rhoq%Gimp,      disp1(11), ierr)
  call MPI_Get_address(t_rhoq%tau,       disp1(12), ierr)

  base  = disp1(1)
  disp1 = disp1 - base

  blocklen1(1)  = t_rhoq%Nscoef
  blocklen1(2)  = 3*t_rhoq%Nscoef
  blocklen1(3)  = 1
  blocklen1(4)  = 3*t_rhoq%nkpt
  blocklen1(4)  = t_rhoq%nkpt
  blocklen1(6)  = 3*t_rhoq%natyp
  blocklen1(7)  = 3*(t_rhoq%Nscoef+1)
  blocklen1(8)  = t_rhoq%Nlayer*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(9)  = t_rhoq%Nlayer*t_rhoq%nkpt*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(10) = t_rhoq%Nscoef*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(11) = t_rhoq%Nscoef*t_rhoq%Nscoef*t_rhoq%lmmaxso*t_rhoq%lmmaxso
  blocklen1(12) = t_rhoq%Nscoef*t_rhoq%Nscoef*t_rhoq%lmmaxso*t_rhoq%lmmaxso

  etype1(1)    = MPI_INTEGER
  etype1(2:7)  = MPI_DOUBLE_PRECISION
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
  if(.not.allocated(t_rhoq%volcub)) allocate(t_rhoq%volcub(t_rhoq%nkpt), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating volcub in rhoq'
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
  if (ierr/=0) stop '[read_scoef_rhoq] file "scoef" not found'
  read(32452345,*) natomimp
  ! exclude probe layer:
  natomimp = natomimp-1
  if(natomimp==0) stop '[read_scoef_rhoq] Error only a single line found in scoef file! Did you include the probe layer?'
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
  
  open(9999, file='mu0', form='formatted')
  read(9999,*) tmp, Nlayer
  close(9999)
  
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
  
end subroutine read_input_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine save_kmesh_rhoq(t_rhoq,nkpt,kpt,volcub,volbz)
  ! save kmesh information in t_rhoq
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: nkpt ! total number of kpoints
  double precision, intent(in) :: volbz ! Brillouin zone volume -> integration weight
  double precision, intent(in) :: kpt(3,nkpt) ! coordinates in reciprocal space of the kpoints
  double precision, intent(in) :: volcub(nkpt) ! volume of kpoint cube
  !local
  integer :: ierr
  
  t_rhoq%nkpt = nkpt
  ierr = 0
  t_rhoq%volbz = volbz
  if(.not.allocated(t_rhoq%kpt)) allocate(t_rhoq%kpt(3,nkpt), stat=ierr)
  if(ierr/=0) stop '[save_kmesh_rhoq] error allocating kpt in rhoq'
  t_rhoq%kpt = kpt
  if(.not.allocated(t_rhoq%volcub)) allocate(t_rhoq%volcub(nkpt), stat=ierr)
  if(ierr/=0) stop '[save_kmesh_rhoq] error allocating volcub in rhoq'
  t_rhoq%volcub = volcub
  
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
  integer :: mu_cls, mu_orig ! layer indices for probe position, cluster center and origin position, respectively
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
  
!   write(*,*) 'rbasis',r_basis
!   write(*,*) 'rcls_i', rcls_i
    
  ! find mu_orig
  ilayer = 0
  irun = 1
  do while(irun==1 .and. ilayer<t_rhoq%natyp)
    ilayer = ilayer+1
    if(dsqrt(sum(r_basis(:,ilayer)**2))<eps) irun = 0
  end do
  mu_orig = ilayer
!   write(*,*) 'mu_orig', mu_orig
  
  ! find mu_cls
  ilayer = 0
  irun = 1
  do while(irun==1)
    ilayer = ilayer+1
    if(dsqrt(sum(rcls_i(:,ilayer)**2))<eps) irun = 0
  end do
  mu_cls = t_rhoq%ilay_scoef(ilayer)
!   write(*,*) 'mu_cls', mu_cls
  
  ! set R_mu and R_cls from mu_orig and mu_cls
  R_mu(:) = r_basis(:,t_rhoq%mu_0)
  R_cls(:) = r_basis(:,mu_cls)
  
  ! L_i = r^cls_i - \Chi^{\mu_i} + \Chi^{\mu_{cls}} - 
  do ilayer=1,t_rhoq%Nscoef
    L_i(:,ilayer) = rcls_i(:,ilayer)+R_cls(:)-R_mu(:)+r_basis(:,t_rhoq%mu_0)-r_basis(:,t_rhoq%ilay_scoef(ilayer))
!     write(*,'(A, 100ES15.7)') 'L_i', L_I(:,ilayer),r_basis(:,t_rhoq%ilay_scoef(ilayer)),rcls_i(:,ilayer),R_cls(:),R_mu(:)
  end do
!   write(*,*) 'R_mu',R_mu
!   write(*,*) 'R_cls',R_cls
!   write(*,*) 'rcls_i',rcls_i(:,:)
  
  ! save result and exit
  ierr = 0
  if(.not.allocated(t_rhoq%L_i)) allocate(t_rhoq%L_i(3,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[get_L_vecs_rhoq] error allocating L_i in rhoq'
  t_rhoq%L_i(:,:) = L_i(:,:)
  deallocate(rcls_i, r_basis, L_i)

end subroutine get_L_vecs_rhoq


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
  integer :: icls, lm1, lm2, idummy1, idummy2, nClustertest
  double precision :: Rclstest(3)
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
      read(562,"(2I5,4e17.9)") idummy1, idummy2, tmat(lm2,lm1), deltamat_dummy
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
  integer :: ienergy, lm1, lm2, id1, id2
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


subroutine calc_Q_mu_rhoq(lmax, ntotd, npan_tot, &
            & nsra, lmmaxso, Rll, Rllleft, ipan_intervall, &
            & irmdnew, trq_of_r )
  ! calculate prefactor needed in rhoq calculation:
  ! Q^mu_Lambda,Lamabda' = Int d^3r Tr[ R^mu_Lmabda(\vec{r}) R^left,mu_Lmabda'(\vec{r}) ] where R^left denotes the left solution
  implicit none
  integer, intent(in) :: lmax, ntotd, npan_tot
  integer, intent(inout) :: nsra, lmmaxso, irmdnew ! left/right solutions (wave functions)
!   integer, intent(in) :: nsra, lmmaxso, irmdnew ! left/right solutions (wave functions)
  double complex, intent(in) :: Rll(nsra*lmmaxso,lmmaxso,irmdnew), Rllleft(nsra*lmmaxso,lmmaxso,irmdnew) ! wave functions in new mesh
  integer, intent(in) :: ipan_intervall(0:ntotd)
  double complex, intent(out) :: trq_of_r(lmmaxso,lmmaxso, 725, irmdnew)
  
  double precision, parameter :: c0ll = 1.0d0/sqrt(16.0d0*datan(1.0d0)) ! Y_00 = 1/sqrt(4*pi) needed for theta_00 * Y_00 = 1, other spherical harmonics are in Gaunt coefficients cleb(:,:)
  
  ! local
  integer :: imt1, lm1, lm2, lm3, ir, j, lm01, lm02, lmshift(2), ispin
  integer :: ncleb, lpot, iend, ierr ! Gaunt coefficients
  double precision, allocatable :: cleb(:), WG(:), YRG(:,:,:) ! Gaunt coefficients
  integer, allocatable :: icleb(:,:), loflm(:), JEND(:,:,:) ! Gaunt coefficients
  
  ! parameter
  double complex, parameter :: C0=( 0.0d0, 0.0d0 ), C1=( 1.0d0, 0.0d0 )
  
  
  ! compute gaunt coefficients
  NCLEB = (lmax*2+1)**2 * (lmax+1)**2
  LPOT = 4*lmax
  
  allocate( WG(4*LMAX), YRG(4*LMAX,0:4*LMAX,0:4*LMAX), CLEB(NCLEB), stat=ierr)
  if(ierr/=0) stop '[calc_Q_mu_rhoq] Error allocating wg etc for gaunt2, gaunt'
  allocate( ICLEB(NCLEB,3), JEND((LPOT+1)**2,0:LMAX,0:LMAX), LOFLM((2*LPOT+1)**2), stat=ierr )
  if(ierr/=0) stop '[calc_Q_mu_rhoq] Error allocating icleb etc for gaunt'
  
  CALL GAUNT2( WG, YRG, 4*LMAX )
  
  CALL GAUNT_new( LMAX, LMAX, LPOT, WG, YRG, CLEB, LOFLM, ICLEB, IEND, JEND, NCLEB, LMAX, (LPOT+1)**2 )
  ! done computing gaunts
  
  ! initialization
  trq_of_r = C0
  
  ! get indices of imt and irmd in new mesh -> boundaries for ir loops
  imt1=ipan_intervall(npan_tot)+1
  

  ! define lmshifts for loop over spins
  lmshift(1) = 0
  lmshift(2) = lmmaxso/2

  ! compute Tr[ R^mu_Lmabda(\vec{r}) R^left,mu_Lmabda'(\vec{r}) ]_L
  do ispin=1,2 ! sum over spins
   do lm01=1,lmmaxso/2
   do lm02=1,lmmaxso/2
   
    ! diagonal in lm1, lm2, lm3 -> spherical part
    do lm1=1,lmmaxso/2
      do ir=1,irmdnew
        ! sum of Rll*Rllleft with Gaunt coefficients, with trace over big and small component (lmmaxso shift in first index of Rll, Rllleft)
        trq_of_r(lm01, lm02, lm1, ir) = trq_of_r(lm01, lm02, lm1, ir) + ( Rll(lm01, lm1, ir)*Rllleft(lm02, lm1, ir) + Rll(lm01+lmmaxso, lm1, ir)*Rllleft(lm02+lmmaxso, lm1, ir) )
      enddo ! ir
    end do ! lm1
   
    ! non-spherical part
    do j = 1,iend ! loop over number of non-vanishing Gaunt coefficients
      ! set lm indices for Gaunt coefficients
      lm1 = icleb(j,1) + lmshift(ispin) ! lmshift shifts between up and down spins
      lm2 = icleb(j,2) + lmshift(ispin) ! same entry since under trace the spinors can be written as: X_s1*X_s2^\dagger = X_s2^\dagger*X_s1 = \delta_{s1,s2}
      lm3 = icleb(j,3)
      ! do loop over non-spherical points here
!       do ir=imt1+1,irmdnew
      do ir=1,irmdnew
        ! sum of Rll*Rllleft with Gaunt coefficients, with trace over big and small component (lmmaxso shift in first index of Rll, Rllleft)
        trq_of_r(lm01, lm02, lm3, ir) = trq_of_r(lm01, lm02, lm3, ir) + cleb(j) * ( Rll(lm01, lm1, ir)*Rllleft(lm02, lm2, ir) + Rll(lm01+lmmaxso, lm1, ir)*Rllleft(lm02+lmmaxso, lm2, ir) )
      enddo ! ir
    enddo ! j -> sum of Gaunt coefficients
    
   enddo ! lm02
   enddo ! lm01
  enddo ! ispin -> sum over spins
  
  
  deallocate( WG, YRG, CLEB )
  deallocate( ICLEB, JEND, LOFLM )
  
  
end subroutine calc_Q_mu_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_G0_k(tau0_k, tinv_i, tinv_j, G0_k, mu_i, mu_j, N)
  ! Calculate G0 from tau_0 and t
  ! G0^mu,mu_i (k) = -t_mu^-1 delta(mu-mu_i) + t_mu^-1 tau0^mu,mu_i(k) t_mu_i^-1
  ! 
  ! This comes from 
  ! $G^0_{i,i'} = -t_i^{-1} \delta_{i,i'} + t_i^{-1} \tau^0_{i,i'} t_{i'}^{-1}$
  ! and the definition of the fourier transform
  ! $\tau^0_{i,i'} = \frac{1}{\Omega_{BZ}}\int d^2k \tau^0_{mu_i,mu_{i'}}(\vec{k}) e^{-i\vec{k}\cdot(\vec{L_i}-\vec{L_{i'}})}$.
  ! Then this gives
  ! $G^0_{mu_i,mu_{i'}}(\vec{k}) = -t_{mu_i}^{-1} \delta_{mu_i,mu_{i'}} + t_{mu_i}^{-1} \tau^0_{mu_i,mu_{i'}}(\vec{k}) t_{mu_{i'}}^{-1}$,
  ! where $t_{mu_i}\equiv t_{i}$.
  !                                       -> expfac
  ! finally construct G0^mu,mu_i (k) * exp(-i k*L_i)
  implicit none
  integer, intent(in) :: N, mu_i, mu_j
  double complex, intent(in)  :: tinv_i(N,N), tinv_j(N,N), tau0_k(N,N)
  double complex, intent(out) :: G0_k(N,N)
  ! local
  double complex :: temp(N,N)
  double complex, parameter :: C0=(0.0d0, 0.0d0), Ci=(0.0d0, 1.0d0), C1=(1.0d0, 0.0d0)
  
  ! initialize G0_k
  G0_k = C0
  temp = C0
  
  ! treat diagonal element from first term:
  if(mu_i==mu_j) G0_k = -tinv_i
  
  !treat second term
  !temp = tau0_k*tinv_j
  call ZGEMM('n','n',N,N,N,C1,tau0_k,N,tinv_j,N,C0,temp,N)
  !G0_k = G0_k + tinv_i*tau0_k*tinv_j = G0_k + tinv_i*temp
  call ZGEMM('n','n',N,N,N,C1,tinv_i,N,temp,N,C1,G0_k,N)
 
end subroutine calc_G0_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine red_Q(t_rhoq, rb, qvec, iq)
  ! folding Q back into Brillouin zone, i.e. find q such that Q=G+q with a reciprocal lattice vector G
  ! this is done with 
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  double precision, intent(in) :: qvec(2), rb(2,2)
  integer, intent(out) :: iq
  ! local
  double precision :: rbt(2,2), irbt(2,2), rq(2), kpt(2), det, qv(2), minimum(3)
  integer :: i, j
  double precision, parameter :: eps=1.0d-04
  
  ! compute inverse of transpose of reciprocal Bravais matrix
  ! note that transpose is needed because vectors are stored in rows and not in columns in the code
  ! fill transpose
  do i=1,2
    rbt(i,i)=rb(i,i)
    do j=1,2
      if(j/=i) rbt(i,j) = rb(j,i)
    end do
  end do
  ! compute 2x2 inverse
  det = 1.0d0/(rbt(1,1)*rbt(2,2) - rbt(1,2)*rbt(2,1))
  irbt(1,1) =  det*rbt(2,2)
  irbt(1,2) = -det*rbt(1,2)
  irbt(2,1) = -det*rbt(2,1)
  irbt(2,2) =  det*rbt(1,1)
  
  ! compute rq = (RB^T)^-1.Q
  ! this gives coordinates in multiples of the reciprocal lattice vectors
  do i=1,2
    rq(i) = irbt(1,i)*qvec(1) + irbt(2,i)*qvec(2)
  end do
  
  ! take modulus rq = mod(rq,1)
  do i=1,2
    rq(i) = rq(i)+100-int(rq(i)+100.0d0)
    if(1-rq(i)<eps) rq(i) = 0.0d0
  end do
  
  ! convert coordinates back
  ! qv = RB^T.rq
  do i=1,2
    qv(i) = rb(i,1)*rq(1) + rb(i,2)*rq(2)
  end do
  
  ! find iq from projection
  iq = -1
!   minimum = 100.0d0
  do i=1,t_rhoq%Nkpt
    kpt(:) = t_rhoq%kpt(:,i)
    if(dabs(dsqrt((kpt(1)-qv(1))**2+(kpt(2)-qv(2))**2))<eps) then
      iq=i
    end if
!     if(minimum(1)>dabs(dsqrt((kpt(1)-qv(1))**2+(kpt(2)-qv(2))**2))) then
!       minimum(1) = dabs(dsqrt((kpt(1)-qv(1))**2+(kpt(2)-qv(2))**2))
!       minimum(2:3) = kpt(1:2)
!     end if
  end do
  
  !check if iq was found
  if(iq<1) write(*,*) qv!, minimum
  if(iq<1) stop '[red_Q] Error: no iq found'
  
  
end subroutine red_Q


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_rhoq(t_rhoq, lmmaxso, Nkp, trq_of_r, rhoq, recbv, lmax,   &
        &        ntotd, npan_tot, ncheb, rpan_intervall, ipan_intervall, irmdnew, alat )
  ! calculate Delta rho^mu(q) = Int{ Delta rho(q, X_mu+r), dr }
  !                           = -1/(2*pi*i) Tr[ Q^mu Int{ Ghost(k) tau Ghost(k+q),dk } - Q^mu,* Int{ Ghost^*(k) tau^* Ghost^*(k-q),dk }]
  use omp_lib
  use mod_timing
  use IFPORT ! random numbers
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: lmmaxso, Nkp, irmdnew, lmax
  double complex, intent(inout) :: trq_of_r(lmmaxso,lmmaxso, 725, irmdnew)
  double precision, intent(in) :: recbv(3,3), alat
  double complex, allocatable, intent(out) :: rhoq(:)
  
  integer :: ntotd, npan_tot, ncheb, npan_log, npan_eq
  integer npan_lognew,npan_eqnew
  double precision, intent(in) :: rpan_intervall(0:ntotd)
  integer, intent(in) :: ipan_intervall(0:ntotd)
  
  ! local
  integer :: lrec, i, irec, k, mu, mu_i, ie, N, q, j, kpq, ierr, lm1, ix, mythread, nthreads, Ni, Nj, Nqpt, ixyz, lm01, lm02, ir, l1, m1, imt1, lm2, lm3, ifun, imin, ielast, nref, ikx, iky, Nx, Ny
  integer, allocatable :: ipvt(:), refpot(:)
  double complex, allocatable :: tinv(:,:,:), tau0_k(:,:), G0ij_k(:,:,:,:), G0ji_k(:,:,:,:), tau(:,:,:,:), tmpG0(:,:), tmpG02(:,:), tll(:,:), trefll(:,:,:), tmp(:,:), tmpsum1(:,:), tmpsum2(:,:), exG0_tmp(:,:), qint(:,:,:), eiqr_lm(:,:), q_mu(:,:), cYlm(:), jl(:)
  double complex :: tr_tmpsum1, tr_tmpsum2, kweight, Z
  double precision, allocatable :: kpt(:,:), L_i(:,:), Qvec(:,:), Ylm(:), rnew(:), qvec_tmp(:,:)
  integer, allocatable :: qvec_index_tmp(:,:), qvec_index(:,:)
  double precision :: QdotL, tmpk(3), tmpr(3), Rq, costheta, phi, box(3)
  integer, allocatable :: ifunm(:), lmsp(:), ntcell(:) ! shape functions
  double precision, allocatable :: thetasnew(:,:) ! shape functions in Chebychev mesh (new mesh)
  logical, allocatable :: kmask(:)

  ! parameters
  double complex, parameter :: C0=(0.0d0, 0.0d0), Ci=(0.0d0, 1.0d0), C1=(1.0d0, 0.0d0)
  double precision, parameter :: pi = 4.0d0*datan(1.0d0) ! pi/4 = arctan(1)
  double precision, parameter :: eps = 1.0d-10 ! small epsilon
  double precision, parameter :: c0ll = 1.0d0/sqrt(16.0d0*datan(1.0d0)) ! Y_00 = 1/sqrt(4*pi) needed for theta_00 * Y_00 = 1, other spherical harmonics are in Gaunt coefficients cleb(:,:)
  
  INTEGER IEND, NCLEB, LMAX_1, LMAX_2, LPOT_2, lm0    
  INTEGER, allocatable :: ICLEB(:,:),JEND(:,:,:),LOFLM(:), IRCUT(:)
  DOUBLE PRECISION, allocatable :: WG(:),YRG(:,:,:),CLEB(:)
  integer :: irid, irmd, ipand, irmin, irws, ipan, mu_0, nfund
  double precision, allocatable :: THETAS(:,:), R(:)
  double precision :: r_log
  
  ! allocate kpoint arrays
  ierr = 0
  N = t_rhoq%lmmaxso
  if(N/=lmmaxso) stop '[calc_rhoq] lmmaxso input does not match value in t_rhoq'
  if(Nkp/=t_rhoq%nkpt) stop '[calc_rhoq] Nkp input does not match value in t_rhoq'
  
  allocate(kpt(3,Nkp), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating kpt'
  kpt = t_rhoq%kpt
  
  allocate(L_i(3,T_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating L_i'
  L_i = t_rhoq%L_i
  
  allocate(kmask(Nkp), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating kmask'
  kmask(:) = .true.
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate scattering path operator tau from host, tmat and impurty GF
  
  !allocate temporary arrays
  allocate(tll(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tll'
  allocate(ipvt(N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating ipvt'
  allocate(tinv(N,N,t_rhoq%Nscoef+1), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tinv'
  allocate(tau0_k(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tau0_k'
  allocate(G0ij_k(N,N,t_rhoq%Nscoef,Nkp), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating G0ij_k'
  allocate(G0ji_k(N,N,t_rhoq%Nscoef,Nkp), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating G0ji_k'
  
  tll = C0
  tinv = C0
  tau0_k = C0
  G0ij_k = C0
  G0ji_k = C0
       
  ie = 1
  ielast = 1
  
  open(9999, file='refinfo', form='formatted')
  read(9999,'(1I9)') NREF
  allocate(refpot(1:t_rhoq%NATYP), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] ERROR in allocating refpot'
  read(9999, '(1000I9)') refpot(1:t_rhoq%NATYP)
  close(9999)
  
  
  allocate(trefll(N,N, NREF), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating trefll'
  open(9999, file='tref', access='direct', recl=4*lmmaxso**2 )
  ie = 2
  ielast = 3
  DO I = 1,NREF
     irec = ie + ielast*(i-1)
     read(9999, rec=irec) trefll(1:N,1:N,i)
     write(9992,'(2ES15.7)') trefll(1:N,1:N,i)
  ENDDO ! I1
  close(9999)
  
  
  ! read tmat and compute tinv = tmat^-1
  lrec = 4*N**2
  open(9999, file='tmat', access='direct', recl=lrec, form='unformatted')
  
  do i=1,t_rhoq%Nscoef
  
    irec = ie+ielast*(t_rhoq%ilay_scoef(i)-1) ! layer index 
    read(9999,rec=irec) tll(1:N,1:N)
    ! t-t_ref
    tll = tll - trefll(:,:,refpot(t_rhoq%ilay_scoef(i)))

    ! initialize tinv
    tinv(1:N,1:N,i)=C0
    do lm1=1,lmmaxso
      tinv(lm1,lm1,i)=C1
    enddo
    
    ! LU decomposition to invert tll
    call zgetrf(lmmaxso,lmmaxso,tll,lmmaxso,ipvt,ierr)
    call zgetrs('n',lmmaxso,lmmaxso,tll,lmmaxso,ipvt,tinv(:,:,i),lmmaxso,ierr)
    
  end do ! i
  
  irec =  ie+ielast*(t_rhoq%mu_0-1) ! layer index for probe layer mu_0
  !write(*,*) 'read tmat_mu0', irec
  read(9999,rec=irec) tll(1:N,1:N)
  ! t-t_ref
  tll = tll - trefll(:,:,refpot(t_rhoq%mu_0))

  ! initialize tinv
  tinv(1:N,1:N,t_rhoq%Nscoef+1)=C0
  do lm1=1,lmmaxso
    tinv(lm1,lm1,t_rhoq%Nscoef+1)=C1
  enddo
  ! LU decomposition to invert tll
  call zgetrf(lmmaxso,lmmaxso,tll,lmmaxso,ipvt,ierr)
  call zgetrs('n',lmmaxso,lmmaxso,tll,lmmaxso,ipvt,tinv(:,:,t_rhoq%Nscoef+1),lmmaxso,ierr)

  ! close file 'tmat'
  close(9999)

!   close(6699)
  
  
  deallocate(tll, ipvt, stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error deallocating tll etc.'
  
  
  ! open tau0_k file, written in writegreen of zulapi code
!   lrec = 4*(N**2)
  lrec = 4*(N**2+1)
!   lrec = 4*(N**2+2)
!   open(9999, file='tau0_k', access='direct', form='unformatted', recl=lrec)
  
  ! allocate temporary arrays
  allocate(tmp(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmp'
  allocate(tmpG0(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmpG0'
  allocate(tmpG02(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmpG02'
  allocate(tmpsum1(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmpsum1'
  allocate(tmpsum2(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmpsum2'
  tmpsum1 = C0
  tmpsum2 = C0
  
  ! find imin
  imin = 1000
  do i=1,t_rhoq%Nscoef
    if(t_rhoq%ilay_scoef(i)<imin) imin = t_rhoq%ilay_scoef(i)
    write(*,*) 'find imin' ,imin,t_rhoq%ilay_scoef(i), t_rhoq%Nlayer
  end do
  
  lm1 = 0
  lm2 = 0
  
  call timing_start('calc rhoq - comp tau')
  ! calculate tau
  write(*,*) 'calculate G0_k from tau0_k'
  write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
  write(*,FMT=190) !beginning of statusbar
  do k=1,Nkp
    ! find exG0_i(k) vector of lm-blocks (i is one component)

    tmpsum1 = C0
    tmpsum2 = C0

    do i=1,t_rhoq%Nscoef
    
      ! set mu, mu_i to calculate G0ij_k
      mu = imin !t_rhoq%ilay_scoef(t_rhoq%Nlayer)
      mu_i = t_rhoq%ilay_scoef(i)
      
      ! read in tau0(k0) from file
      ie = 2 ! maybe interpolation to real axis?
      irec = ((t_rhoq%Nlayer-1)*2)*(k-1) + ((t_rhoq%Nlayer-1)*2)*Nkp*(ie-1-1)
      ix = mu_i - mu + 1
      ! first mu_0,mu_i element
      irec = irec + ix + (t_rhoq%Nlayer-1)
      read(998899,'(10000ES15.7)') tmpk(1:2), tau0_k(1:N,1:N)
      
!       call red_Q(t_rhoq, recbv(1:2,1:2), tmpk(1:2), q)
      q = k


      ! set kmask
      if(dsqrt(dreal(tau0_k(1,1))**2+dimag(tau0_k(1,1))**2)<eps) then
        kmask(q) = .false.
        G0ij_k(1:N,1:N,i,q) = C0
        G0ji_k(1:N,1:N,i,q) = C0
        if(i==1) lm2 = lm2+1
      else
        kmask(q) = .true.
        if(i==1) lm1 = lm1+1
      end if      

      
      if(kmask(q)) then
        ! find G0ij_k from tau0_k and tinv
        call calc_G0_k(tau0_k, tinv(1:N,1:N,t_rhoq%Nscoef+1), tinv(1:N,1:N,i), tmpG0(1:N,1:N), t_rhoq%mu_0, mu_i, N)
        G0ij_k(1:N,1:N,i,q) = tmpG0(1:N,1:N)
!         G0ij_k(1:N,1:N,i,q) = tau0_k(1:N,1:N)
      end if
      
      ! now mu_i,mu_0 element
      irec = ((t_rhoq%Nlayer-1)*2)*(k-1) + ((t_rhoq%Nlayer-1)*2)*Nkp*(ie-1-1)
      ix = mu_i - mu + 1
      irec = irec + ix
      read(998888,'(10000ES15.7)') tmpk(1:2), tau0_k(1:N,1:N)
      
      if(kmask(q)) then
        !find G0ji_k from tau0_k and tinv
        call calc_G0_k(tau0_k, tinv(1:N,1:N,i), tinv(1:N,1:N,t_rhoq%Nscoef+1), tmpG02(1:N,1:N), mu_i, t_rhoq%mu_0, N)
        G0ji_k(1:N,1:N,i,q) = tmpG02(1:N,1:N)
!         G0ji_k(1:N,1:N,i,q) = tau0_k(1:N,1:N)
      end if

      tmpsum1 = tmpsum1 + G0ij_k(1:N,1:N,i,q)
      tmpsum2 = tmpsum2 + G0ji_k(1:N,1:N,i,q)
      
    end do ! i

    ! test writeout of qdos (w/o single site contribution which is the same for all kpts)
    tmp(1,1) = C0
    do i=1,lmmaxso
       tmp(1,1) = tmp(1,1) + tmpsum1(i,i)
    end do
    write(22345, '(1000ES15.7)') tmpk(1:2), tmp(1,1)
    !write(12345, '(1000ES15.7)') tmpk(1:2), tmp(1,1)
    tmp(1,1) = C0
    do i=1,lmmaxso
       tmp(1,1) = tmp(1,1) + tmpsum2(i,i)
    end do
    write(22346, '(1000ES15.7)') tmpk(1:2), tmp(1,1)
    !write(12346, '(1000ES15.7)') tmpk(1:2), tmp(1,1)
    
    
    !update statusbar
    if(Nkp>=50.and.mod(k,Nkp/50)==0) write(*,FMT=200)
    
  end do ! k
  
  write(*,*) 'kmask info:', lm1, lm2
  
!   close(998899)
!   close(998888)
  deallocate(tmp, tmpG0, tmpG02)
  deallocate(tmpsum1, tmpsum2)
  
  
  ! take tau from t_rhoq and reshape it 
  allocate( tau(lmmaxso,lmmaxso,t_rhoq%Nscoef,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tau imp'
  
  ! reshape tau
  do j=1,t_rhoq%Nscoef
    do i=1,t_rhoq%Nscoef
      do lm2=1,lmmaxso
        do lm1=1,lmmaxso
          tau(lm1, lm2 ,i, j) = t_rhoq%tau(lm1+lmmaxso*(i-1),lm2+lmmaxso*(j-1))
        end do
      end do
    end do
  end do
  
  write(*,*) 'done with tau'

  ! done calculating tau
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call timing_stop('calc rhoq - comp tau')
  
    
  
  call timing_start('calc rhoq - Gaunt')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! find second set of gaunt coefficients for sum_LL'L" [exp(-iqr)_L * Tr(R*Rleft)_L' * Theta)L" * C_LL'L"] which is then integrated radially
  
  NCLEB =   ((lmax*2+1)**2 * (lmax+1)**2 )*100
  LMAX_1 = 6
  LMAX_2 = 18
  LPOT_2 = 2*(LMAX_1+LMAX_2)
  
  allocate( WG(4*LMAX_2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating wg '
  allocate( YRG(4*LMAX_2,0:4*LMAX_2,0:4*LMAX_2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating yrg'
  allocate( CLEB(NCLEB), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating cleb'
  allocate( ICLEB(NCLEB,3), JEND((LPOT_2+1)**2,0:LMAX_1,0:LMAX_2), LOFLM((2*LPOT_2+1)**2), stat=ierr )
  if(ierr/=0) stop '[calc_rhoq] Error allocating icleb_2 etc for gaunt'
  
  !          ( > |  > |     <    )
  CALL GAUNT2( WG, YRG, 4*LMAX_2)
  
  !             (    <  |    <  |   <   | < |  < |  >  |  >   |   >  |  >  |  >  |  >   |   <   |     <         )
  CALL GAUNT_new( LMAX_1, LMAX_2, LPOT_2, WG, YRG, CLEB, LOFLM, ICLEB, IEND, JEND, NCLEB, LMAX_2  , (LPOT_2+1)**2 )
  
  
  write(*,*) 'done with gaunts'
  
  ! done finding Gaunts for larger lmax values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call timing_stop('calc rhoq - Gaunt')
  
  
  call timing_start('calc rhoq - shape')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! find shapefunctions
  
  ! some dimensions, needed in read in routine
  irid = 900
  irmd = 900
  ipand = 50
  nfund = 500
  
  allocate( ifunm( (2*lpot_2+1)**2 ), lmsp( (2*lpot_2+1)**2 ) ) ! shape functions
  allocate( thetasnew(NTOTD*(NCHEB+1),nfund) ) ! shape functions in Chebychev mesh (new mesh)
  allocate( NTCELL(1:t_rhoq%NATYP) )
  allocate( thetas(irid,nfund) )
  allocate( rnew(irmdnew) )
  
  rnew = 0.0d0
  
  thetas = 0.0d0
  thetasnew = 0.0d0
  
  ! read in corresponding mesh information
  open(9999, file='rmesh_mu0.txt', form='formatted')
  read(9999,*) !'# mu_0, IRMIN, IRWS, IRNS'
  read(9999,'(4I9)') mu_0, irmin, irws, ipan
  if(mu_0/=t_rhoq%mu_0) stop 'Error: mu0 value does not match!!!'
  allocate( R(1:irws) )
  allocate( ircut(1:ipan) )
  read(9999,*) !'# R(1:IRWS)'
  read(9999,'(1000E22.15)') R(1:irws)
  read(9999,*) !'# NTCELL(1:NATYP)'
  read(9999,'(1000I9)') ntcell(1:t_rhoq%natyp)
  read(9999,'(A)') !'# IRCUT(1:IPAN)'
  read(9999,'(1000I9)') ircut(1:ipan)
  read(9999,'(A)') !'# R_LOG, NPAN_LOG, NPAN_EQ'
  read(9999,'(E22.15,2I9)') R_LOG, NPAN_LOG, NPAN_EQ
  close(9999)
  
  ! read shapefunction for atom mu_0 from file shapefun (Note: this file needs to be )
  call read_shape(thetas, lmsp, ifunm, irid, irmd, nfund, lpot_2, t_rhoq%mu_0, ipan, ntcell, t_rhoq%natyp)
     
  ! interpolate shapefunction from old mesh (stored in r(1:irws)) to chebychev mesh (output as rnew)
  call interpol_shape(R,irmin,irws,ipan,ircut,r_log,npan_log,npan_eq,ncheb,           &
     &                npan_tot,rnew,rpan_intervall,ipan_intervall,thetas,thetasnew,     &
     &                nfund, npan_lognew, npan_eqnew)
     
  ! deallocate arrays in old mesh, from here on only the new (chebychev) mesh is used
  deallocate( R, ircut, thetas, ntcell )
  
  
  write(*,*) 'done with shapes'
  
  ! done finding shapefunction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call timing_stop('calc rhoq - shape')
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !find q-mesh from k-kpoints
  
  ! parameters for bigger q mesh
!   Ni = 1
!   Nj = 1
  Ni = 2
  Nj = 2
!   Ni = 3
!   Nj = 3
!   Ni = 4
!   Nj = 4
  Nqpt = Ni*Nj*Nkp
  
  Nx = int(dsqrt(dfloat(Nkp)))
  Ny = Nx
  if(Nx*Ny/=Nkp) stop '[calc_rhoq] Error: Nx*Ny/=Nkp'
  
  
!   write(*,*) 'alloc qvec', nkp, nqpt
  allocate(Qvec(3,Nqpt), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating Qvec'
  allocate(Qvec_index(3,Nqpt), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating Qvec_index'
  ! temporary arrays for reduces set of qpts in box (see below)
  allocate(Qvec_tmp(3,Nqpt), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating Qvec_tmp'
  allocate(Qvec_index_tmp(3,Nqpt), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating Qvec_index_tmp'

  ! find qvecs based on kpts from integration
  
!   open(888, file='qvec_index.dat', form='formatted')
  
  do k=1,Nkp
    do i=1,Ni
      do j=1,Nj
        irec = k + Nkp*(i-1) + Nkp*Ni*(j-1)
        Qvec(1:3,irec) = kpt(1:3,k) + (i-Ni)*recbv(1:3,1) + (j-Nj)*recbv(1:3,2)
        
        ! find kx, ky indices from k=iky+(ikx-1)*Ny = 1...Nx*Ny
        ikx = mod(k,Nx)
        if (ikx<=0) ikx = ikx+Nx
        iky = (k-ikx)/Nx+1
        if (iky<=0) iky = iky+Ny
        ! store these in qvec_index
        qvec_index(1,irec) = ikx
        qvec_index(2,irec) = iky
        
!         write(888,'(2F14.7,6i9)') qvec(1:2,irec), irec, k, i,j, ikx, iky 
        
      end do !j=1,Nj
    end do !i=1,Ni
    
  end do !k=1,Nkp
  
!   close(888)

  ! reduce qpts to box around Gamma

  box(1) = 0.4d0!1.0d0!0.7d0
  box(2) = 0.4d0!1.0d0!0.7d0
  box(3) = 0.4d0!1.0d0!0.1d0
  
  write(*,*) 'reduce kpts to:', 0.12, box(1)

  k = 0
  do q=1,Nqpt
    if((abs(qvec(1,q))-box(1)<=eps).and.(abs(qvec(2,q))-box(2)<=eps).and.(abs(qvec(3,q))-box(3)<=eps)) then
      if(dsqrt((qvec(1,q))**2+(qvec(2,q))**2+(qvec(3,q))**2)>0.12) then
        k = k+1
        qvec_tmp(1:3,k) = qvec(1:3,q)
        qvec_index_tmp(1:2,k) = qvec_index(1:2,q)
      end if
    end if
  end do

  Nqpt = k
  write(*,*) 'red qvecs box',k, box

  ! change allocation 
  deallocate(qvec)
  allocate(Qvec(3,Nqpt), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating Qvec'
  deallocate(qvec_index)
  allocate(Qvec_index(2,Nqpt), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] Error allocating Qvec'

  qvec(1:3, 1:Nqpt) = qvec_tmp(1:3, 1:Nqpt)
  qvec_index(1:2, 1:Nqpt) = qvec_index_tmp(1:2, 1:Nqpt)

  ! deallocate temporary array
  deallocate(qvec_tmp, qvec_index_tmp)
    
  ! done finding q-mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!   write(*,*) 'alloc rhoq'
  allocate(rhoq(Nqpt), stat=ierr)
  rhoq = C0
  if(ierr/=0) stop '[calc_rhoq] Error allocating rhoq'
  
  
  !print header of statusbar
190      FORMAT('                 |'$)   ! status bar
200      FORMAT('|'$)                    ! status bar
  write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
  write(*,FMT=190) !beginning of statusbar
  
  call timing_start('calc rhoq - q-loop')
  !$omp parallel default(shared) private(q, tmpsum1, tmpsum2, k, kpq, kweight, i, j, tmp, tr_tmpsum1, tr_tmpsum2, mythread, nthreads, irec, ixyz, exG0_tmp, QdotL, tmpk, tmpr, q_mu, ifun, lm1, lm2, lm3, ir, qint, l1, m1, lm01, lm02, lm0, Z, Ylm, cYlm, phi, costheta, Rq, eiqr_lm, jl)
!   write(*,*) 'alloc tmpsum'
  allocate(tmpsum1(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmpsum1'
  allocate(tmpsum2(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmpsum2'
  allocate(tmp(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmp'
  allocate(exG0_tmp(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating exG0_tmp'
  
  mythread = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate rho(q)
    
  ! get indices of imt and irmd in new mesh -> boundaries for ir loops
  imt1=ipan_intervall(npan_tot)+1
  allocate( qint(irmdnew, lmmaxso, lmmaxso), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating qint'
  allocate( eiqr_lm((LMAX_2+1)**2, irmdnew), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating eiqr_lm'
  allocate( jl(0:lmax_2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating jl'
  allocate( Ylm((LMAX_2+1)**2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating Ylm'
  allocate( cYlm((LMAX_2+1)**2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating cYlm'
  allocate( q_mu(lmmaxso,lmmaxso), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating q_mu'
  
  
  !$omp do 
  do q=1,Nqpt !q-loop     
   
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> k-loop >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! initialize kpt integration
     tmpsum1 = C0
     tmpsum2 = C0
     if(mythread==0) call timing_start('calc rhoq - q>k-loop')
     do k=1,Nkp !k-loop integration
       ! new: find kpq index from kvec+Qvec, compared with kpts in BZ
       ! find q from Q(q)
       if(mythread==0) call timing_start('calc rhoq - q>k-loop>red_Q')       
!        call red_Q(t_rhoq, recbv(1:2,1:2), Qvec(1:2,q)+kpt(1:2,k), kpq)
       ! find ikx, iky indices from combined index k=iky+(ikx-1)*Ny
       ikx = mod(k,Nx)
       if (ikx<=0) ikx = ikx+Nx
       iky = (k-ikx)/Nx+1
       if (iky<=0) iky = iky+Ny
       ! create i=(iqx+ikx)%Nx and j=(iqy+iky)%Ny
       i = mod(qvec_index(1,q)+ikx, Nx)
       if (i<=0) i = i+Nx
       j = mod(qvec_index(2,q)+iky, Ny)
       if (j<=0) j = j+Ny
       ! then k+q is given by j+(i-1)*Ny
       kpq = i+(j-1)*Nx
!        if(mythread==0) write(*,*) k,q,kpq
       if(mythread==0) call timing_pause('calc rhoq - q>k-loop>red_Q')
       
       if(kmask(kpq) .and. kmask(k)) then
     
         ! read kweight from memory
         kweight = dcmplx(t_rhoq%volcub(k), 0.0d0) ! complex number needed for zgemm later on
         
         ! kweight = kweight * exp(i q*L_i)
         
         do i=1,t_rhoq%Nscoef
           do j=1,t_rhoq%Nscoef
         
             ! Sum( Int( exG0_i(k+q) tau_i,j exG0_j(k); dk ); i,j)
             
             ! collect phase factors
             if(mythread==0) call timing_start('calc rhoq - q>k-loop>phase1')
             ! exG0 = G0*exp(-i(k+q)*L_i)
             tmpk(:) = Qvec(:,q)+kpt(:,k)
             tmpr(:) = L_i(:,i)
             QdotL = tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3)
             ! exG0 = exG0*exp(+ik*L_j) -> G0tauG0*exp(-[(k+q)*L_i - k*L_j])
             tmpk(:) = kpt(:,k)
             tmpr(:) = L_i(:,j)
             QdotL = QdotL - (tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3))
             ! multiply with phase
             exG0_tmp(1:N,1:N) = G0ji_k(1:N,1:N,j,k)*exp(-2.0d0*pi*Ci*QdotL)
             if(mythread==0) call timing_pause('calc rhoq - q>k-loop>phase1')
             
             ! tmp = tau*exG0_k(k)
             tmp(1:N,1:N) = C0
             if(mythread==0) call timing_start('calc rhoq - q>k-loop>1.ZGEMM')
             call ZGEMM('n','n',N,N,N,C1,tau(1:N,1:N,i,j),N,exG0_tmp(1:N,1:N),N,C0,tmp(1:N,1:N),N)
             if(mythread==0) call timing_pause('calc rhoq - q>k-loop>1.ZGEMM')
             ! tmpsum1 = tmpsum1 + G0_k(k+q)*tmp*kweight
             if(mythread==0) call timing_start('calc rhoq - q>k-loop>2.ZGEMM')
             call ZGEMM('n','n',N,N,N,kweight,G0ij_k(1:N,1:N,i,kpq),N,tmp(1:N,1:N),N,C1,tmpsum1(1:N,1:N),N)
             if(mythread==0) call timing_pause('calc rhoq - q>k-loop>2.ZGEMM')
           
           
             ! Int( exG0(k)^*.tau^*.exG0(k+q)^*, dk )
             
             if(mythread==0) call timing_start('calc rhoq - q>k-loop>phase2')
             ! collect phase factors
             ! exG0 = G0*exp(-i(k+q)*L_j)
             tmpk(:) = Qvec(:,q)+kpt(:,k)
             tmpr(:) = L_i(:,j)
             QdotL = tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3)
             ! exG0 = exG0*exp(+ik*L_i) -> G0tauG0*exp(-[(k+q)*L_j - k*L_i])
             tmpk(:) = kpt(:,k)
             tmpr(:) = L_i(:,i)
             QdotL = QdotL - (tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3))
             ! multiply with phase
             exG0_tmp(1:N,1:N) = dconjg(G0ij_k(1:N,1:N,i,k))*exp(-2.0d0*pi*Ci*QdotL)
             if(mythread==0) call timing_pause('calc rhoq - q>k-loop>phase2')
             
             
             ! tmp = dconjg(tau)*dconjg(exG0_k(k+q))
             tmp(1:N,1:N) = C0
             if(mythread==0) call timing_start('calc rhoq - q>k-loop>3.ZGEMM')
             call ZGEMM('n','n',N,N,N,C1,dconjg(tau(1:N,1:N,i,j)),N,dconjg(G0ji_k(1:N,1:N,j,kpq)),N,C0,tmp(1:N,1:N),N)
             if(mythread==0) call timing_pause('calc rhoq - q>k-loop>3.ZGEMM')
             ! tmpsum2 = tmpsum2 + dconjg(exG0_k(k))*tmp*kweight
             if(mythread==0) call timing_start('calc rhoq - q>k-loop>4.ZGEMM')
             call ZGEMM('n','n',N,N,N,kweight,exG0_tmp(1:N,1:N),N,tmp(1:N,1:N),N,C1,tmpsum2(1:N,1:N),N)
             if(mythread==0) call timing_pause('calc rhoq - q>k-loop>4.ZGEMM')
                        
             end do !j
         end do ! i
         
       end if ! (kmask(kpq) .and. kmask(k))
       
     end do ! k-loop
     if(mythread==0) call timing_pause('calc rhoq - q>k-loop')
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< k-loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eiqr  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     if(mythread==0) call timing_start('calc rhoq - q>eiqr')
     !eiqr_lm generation: exp(-iq*r) = j_l(qr)*Y_L(-q)
     eiqr_lm(:,:) = C0
     if(dsqrt(Qvec(1,q)**2+Qvec(2,q)**2+Qvec(3,q)**2)>=eps) then
        Ylm = 0.0d0
        cYlm = C0
        Rq = dsqrt(Qvec(1,q)**2+Qvec(2,q)**2+Qvec(3,q)**2)
        costheta = -Qvec(3,q)/Rq ! =cos(theta)
        phi = datan2(Qvec(2,q), Qvec(1,q)) + pi
        do l1=0,LMAX_2
          lm0 = l1**2+1 ! = l1*(l1+1)+m1+1 for m1=-l1
          call cspher(cYlm(l1**2+1:(l1+1)**2), l1, costheta) ! computes real spherical harmonic for l=l1 without factor exp(i*m*phi)
          do m1=-l1,l1
            lm01 = l1* (l1+1) + m1 + 1
            cYlm(lm01) = cYlm(lm01) * exp(Ci*m1*phi)
          end do ! m1
        end do ! l1
        
        ! convert complex spherical harmonics to real spherical harmonics
        do l1=0,lmax_2
          do m1=-l1,l1
            if(m1<0) then
              lm0 = l1* (l1+1) + m1 + 1
              lm01 = l1* (l1+1) + abs(m1) + 1
              Ylm(lm0) = dsqrt(2.0d0) * (-1.0d0)**m1 * dimag(cYlm(lm01))
            elseif(m1==0) then
              lm0 = l1* (l1+1) + m1 + 1
              Ylm(lm0) = dreal(cYlm(lm0))
            elseif(m1>0) then
              lm0 = l1* (l1+1) + m1 + 1
              Ylm(lm0) = dsqrt(2.0d0) * (-1.0d0)**m1 * dreal(cYlm(lm0))
            end if
          end do
        end do
        
        ! construct exp(-i q.r)_L
        do ir=1,irmdnew
          Z = dcmplx(Rq*rnew(ir)/alat*2.0d0*pi, 0.0d0)
          call CALC_JLK(jl, Z, LMAX_2)
          do l1=0,lmax_2
            do m1=-l1, l1
              lm01 = L1* (L1+1) + M1 + 1
              eiqr_lm(lm01, ir) = jl(l1) * Ylm(lm01) * (Ci**l1)/c0ll/c0ll
            end do ! m1
          end do ! l1
        end do ! ir
        
     else
     
        eiqr_lm(:,:) = C0
        
     end if
     if(mythread==0) call timing_pause('calc rhoq - q>eiqr')
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  eiqr  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  qint  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     if(mythread==0) call timing_start('calc rhoq - q>qint')
     ! initialize qint
     qint = C0

     ! first treat spherical block (diagonal in lm, lm' -> only lm'==lm)
     do lm01 = 1,lmmaxso
       do lm02 = 1,lmmaxso
           do l1=0,lmax_2
           ! only diagonal (spherical)
           m1 = 0
           lm1 = l1 * (l1+1) + m1 + 1
           ! fill diagonal
           do ir = 1,irmdnew
             qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + trq_of_r(lm01, lm02, lm1, ir) * eiqr_lm(lm1, ir) ! note: diagonal in lm01!!!
!              qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + trq_of_r(lm01, lm02, lm1, ir) ! note: diagonal in lm01!!!
!              qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + eiqr_lm(lm1, ir) ! note: diagonal in lm01!!!
           enddo ! ir
           ! add shapefunction to points with r> R_MT
           do ir=imt1+1,irmdnew
             qint(ir, lm01, lm02) = qint(ir, lm01, lm02)*thetasnew(ir,1)*c0ll ! add shapefunction on diagonal for last points outside of MT radius, c0ll comes from gaunt coefficients
           enddo ! ir
         enddo ! lm1
       enddo ! lm02
     enddo ! lm01
     
     
     
     ! treat non-spherical components -> off-diagonal blocks in lm, lm' matrices
     do j = 1,iend ! loop over number of non-vanishing Gaunt coefficients
       ! set lm indices for Gaunt coefficients
       lm1 = icleb(j,1)
       lm2 = icleb(j,2)
       lm3 = icleb(j,3)
       ! get shape function index for lm3
       ifun = ifunm(lm3)
       ! do loop over non-spherical points here
       if(lmsp(lm3)/=0) then
         do ir=imt1+1,irmdnew
           do lm01=1,lmmaxso
             do lm02=1,lmmaxso
               qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + cleb(j) * trq_of_r(lm01, lm02, lm1, ir) * eiqr_lm(lm2, ir) * thetasnew(ir,ifun)
!                qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + cleb(j) * trq_of_r(lm01, lm02, lm1, ir) * thetasnew(ir,ifun)
!                qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + cleb(j) * eiqr_lm(lm2, ir) * thetasnew(ir,ifun)
             enddo
           enddo
         enddo ! ir
       end if ! lmsp(ifun)/=0
     enddo ! j -> sum of Gaunt coefficients
         
       ! do radial integration
     q_mu = C0
     do lm01 = 1, lmmaxso ! loop over outer lm-component
       do lm02 = 1, lmmaxso ! loop over outer lm-component
         call intcheb_cell(qint(:,lm01,lm02),q_mu(lm01,lm02),rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
       end do ! lm02
     end do ! lm01
     if(mythread==0) call timing_pause('calc rhoq - q>qint')
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  qint  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> rhoq(iq) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! calculate trace of Q^mu times k-kpoint integral
     tr_tmpsum1 = C0
     tr_tmpsum2 = C0
     ! tmpsum1 = Tr{ q_mu Int(G0.tau.G0.exp(...), dk) } = Tr{ q_mu.tmpsum1 }
     tmp = C0
     tmp(1:N,1:N) = tmpsum1(1:N,1:N)
     tmpsum1(1:N,1:N) = C0
     call ZGEMM('n','n',N,N,N,C1,q_mu(1:N,1:N),N,tmp(1:N,1:N),N,C0,tmpsum1(1:N,1:N),N)
     ! tmpsum2 = Tr{ q_mu^* Int(G0^*.tau^*.G0^*.exp(...), dk) } = Tr{ q_mu^*.tmpsum2 }
     tmp = C0
     tmp(1:N,1:N) = tmpsum2(1:N,1:N)
     tmpsum2(1:N,1:N) = C0
     call ZGEMM('n','n',N,N,N,C1,dconjg(q_mu(1:N,1:N)),N,tmp(1:N,1:N),N,C0,tmpsum2(1:N,1:N),N)
           
     ! take trace
     do i=1,t_rhoq%lmmaxso
       tr_tmpsum1 = tr_tmpsum1 + tmpsum1(i,i)
       tr_tmpsum2 = tr_tmpsum2 + tmpsum2(i,i)
!        tr_tmpsum1 = tr_tmpsum1 + q_mu(i,i)
     end do
     
     !$omp critical
     rhoq(q) = 1.d0/(2.d0*Ci) * (tr_tmpsum1 - tr_tmpsum2)
     !$omp end critical
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< rhoq(iq) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !update statusbar
     !$omp critical
     if(Nqpt>=50.and.mod(q,Nqpt/50)==0) write(*,FMT=200)
     !$omp end critical
       
  end do !q
  !$omp end do
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop>red_Q')
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop>phase1')
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop>1.ZGEMM')
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop>2.ZGEMM')
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop>phase2')
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop>3.ZGEMM')
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop>4.ZGEMM')
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop')
  if(mythread==0) call timing_stop('calc rhoq - q>eiqr')
  if(mythread==0) call timing_stop('calc rhoq - q>qint')
  
  deallocate(tmpsum1, tmpsum2, tmp, qint , eiqr_lm, jl, Ylm)
  
  !$omp end parallel
  call timing_stop('calc rhoq - q-loop')
  
  ! done calculating rho(q)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call timing_start('calc rhoq - writeout')
  ! write out result
  write(*,*) '' ! finish statusbar
  write(*,*) 'saving rhoq(q) to out_rhoq.txt'
  open(9999, file='out_rhoq.txt', form='formatted')
  do q=1,Nqpt
     write(9999,'(5E15.7)') Qvec(:,q),rhoq(q)
  end do
  close(9999)
  call timing_stop('calc rhoq - writeout')
  
  
!   deallocate(tau)
!   deallocate(tau0_k)
    
end subroutine calc_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module mod_rhoq


! for testing only, comment out when ready
!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!
program test
  
  use mod_types, only: t_inc
  use mod_rhoq
  use mod_timing
  use mod_version_info
  
  character*256 :: uio ! unit in which to find inputcard
  integer :: nkpt ! total number of kpoints
  integer :: lmmaxso, natyp
  integer :: naez,ncls,nr,nemb
  integer :: lmax, ntotd, npan_tot, ncheb
  integer :: nsra, irmdnew ! left/right solutions (wave functions)
  double precision :: alat
  double precision :: volbz ! Brillouin zone volume -> integration weight
  
  double precision, allocatable :: kpt(:,:), volcub(:) ! coordinates in reciprocal space of the kpoints and volume of kpoint cube
  double precision, allocatable :: r_basis(:,:) ! real space positions of all atoms in host system
  double complex, allocatable   :: Ghost(:,:,:)
  double precision, allocatable :: rcls(:,:,:),rr(:,:)
  integer, allocatable          :: atom(:,:),cls(:),ezoa(:,:),nacls(:)
  double complex, allocatable   :: Rll(:,:,:), Rllleft(:,:,:) ! wave functions in new mesh
  double precision, allocatable :: rpan_intervall(:)
  integer, allocatable          :: ipan_intervall(:)
  double complex, allocatable   :: trq_of_r(:,:,:,:), rhoq(:)
  double precision, allocatable :: recbv(:,:), bravais(:,:) ! lattice information > Bravais matrix in real and reciprocal space
  double precision, allocatable :: rnew(:)  ! r mesh
  
  double precision, parameter :: pi = 4.0d0*datan(1.0d0) ! pi/4 = arctan(1)

  
  integer :: lm1, lm2, i, j, ir
  
  ! set timing output
  call construct_serialnr()
  t_inc%i_time = 1
  call timing_init(0)
  call timing_start('total time')
  
  
  call timing_start('read stuff')
  
  ! read in scalars
  uio = 'inputcard'
  open(9999, file='params.txt')
  read(9999,*) lmmaxso, natyp
!   lmmaxso = lmmaxso/2
  read(9999,*) naez, ncls, nr, nemb, lmax
  read(9999,*) alat
  close(9999)

!   write(*,*) 'params', lmmaxso, natyp, naez, ncls, nr, nemb, lmax, alat
  

  open(9999, file='kpts.txt', form='formatted')
  read(9999,'(I9)') nkpt
  read(9999,'(E16.7)') volbz
  allocate( kpt(3,nkpt), volcub(nkpt) )
  do i=1,nkpt
    read(9999,'(4E16.7)') (kpt(j,i), j=1,3), volcub(i)
  end do
  allocate(recbv(3,3), bravais(3,3))
  read(9999,'(100E16.7)') recbv(1:3,1:3),bravais(1:3,1:3)
  close(9999)
  
  ! convert from alat units to a.u. units (used in the remainder of the code)
  recbv = recbv !* (2.0d0*pi)/alat  !* alat !2*pi! * alat
  kpt = kpt !* (2.0d0*pi)/alat      !* alat !2*pi! * alat
!   recbv = recbv * (2.0d0*pi)/alat  !* alat !2*pi! * alat
!   kpt = kpt * (2.0d0*pi)/alat      !* alat !2*pi! * alat

!   write(*,*) 'kpoints', nkpt, volbz
  
  open(9999, file='host.txt')
  allocate( r_basis(3,natyp), Ghost(lmmaxso,lmmaxso,t_rhoq%Nlayer) )
  allocate( rcls(3,ncls,ncls), rr(3,0:nr) )
  allocate( atom(ncls,naez+nemb), cls(naez+nemb), ezoa(ncls,naez+nemb), nacls(ncls) )
  read(9999,*) r_basis(1:3,1:natyp)!, Ghost(1:lmmaxso,1:lmmaxso,)
  read(9999,*) rcls(1:3,1:ncls,1:ncls), rr(1:3,0:nr), atom(1:ncls,1:naez+nemb)
  read(9999,*) cls(1:naez+nemb), ezoa(1:ncls,1:naez+nemb), nacls(1:ncls)
  close(9999)
  
  open(9999, file='wavefunctions.txt')
  read(9999,'(100I9)') ntotd, npan_tot, ncheb, nsra, irmdnew
  write(*,*) ntotd, npan_tot, ncheb, nsra, irmdnew, lmmaxso
  allocate( rnew(irmdnew) )
!   read(9999,'(1000E26.17)') rnew
!   read(9999,*) rnew(1:irmdnew)
         read(9999,'(1000E26.17)') rnew(1:irmdnew)
  
  allocate( trq_of_r(lmmaxso,lmmaxso, 725, irmdnew) )


  allocate( Rll(1:nsra*lmmaxso,1:lmmaxso,1:irmdnew),          &
 &          Rllleft(1:nsra*lmmaxso,1:lmmaxso,1:irmdnew) ,     &
 &          rpan_intervall(0:ntotd),ipan_intervall(0:ntotd) )
  do ir=1,irmdnew
    do lm1=1,nsra*lmmaxso
      do lm2=1,lmmaxso
        read(9999,'(20000E16.7)') Rll(lm1, lm2, ir) 
        read(9999,'(20000E16.7)') Rllleft(lm1, lm2, ir)
      end do
    end do
  enddo
  do lm1=0,ntotd
    read(9999,'(E16.7,I9)') rpan_intervall(lm1), ipan_intervall(lm1)
  enddo
  close(9999)
  
  call timing_stop('read stuff')
  
  

  call timing_start('init rhoq')
  

  call read_scoef_rhoq(t_rhoq)

  call read_input_rhoq(t_rhoq, uio)

  call save_kmesh_rhoq(t_rhoq,nkpt,kpt,volcub,volbz)

  call save_geometry_rhoq(t_rhoq,r_basis,lmmaxso,natyp)
  
! #ifdef CPP_MPI
!   call bcast_scalars_rhoq(t_rhoq)
! #endif
  call init_t_rhoq(t_rhoq)
! #ifdef CPP_MPI
!   call bcast_arrays_rhoq(t_rhoq)
! #endif

  call timing_stop('init rhoq')
  

!   write(*,*) 'before Gimp,Dt',t_rhoq%Nscoef,allocated(t_rhoq%Dt), allocated(t_rhoq%Gimp)
!   write(*,*) 'shape Gimp,Dt',shape(t_rhoq%Dt), shape(t_rhoq%Gimp)
  call timing_start('read Dt Gimp')
  call read_Dt_Gimp_rhoq(t_rhoq, lmmaxso, t_rhoq%Nscoef)
  call timing_stop('read Dt Gimp')
  
  
  ! calculate impurity scattering path operator
  call timing_start('calc tau')
  call calc_tau_rhoq(t_rhoq)
  call timing_start('calc tau')
  
  
  ! calculate prefactor Q^{\mu}_{LL'} = Tr{ \int{ R_{L}(\vec{r})*Rleft_{L'}(\vec{r}) d\vec{r} }
  call timing_start('calc Q_mu')
  call calc_Q_mu_rhoq(lmax, ntotd, npan_tot, &
        & nsra, lmmaxso, Rll, Rllleft, ipan_intervall, &
        & irmdnew, trq_of_r)
  call timing_stop('calc Q_mu')

  ! calculate Fourier transform: \rho(\vec{q}) = \int \Delta\rho(\vec{q};\Chi_\mu+\vec{r}) d\vec{r}
!   allocate(rhoq(t_rhoq%Nkpt), stat=ierr)
!   if(ierr/=0) stop '[test] Error allocating rhoq'
  call timing_start('calc rhoq')
  call calc_rhoq(t_rhoq, t_rhoq%lmmaxso, t_rhoq%Nkpt, trq_of_r, rhoq, recbv, lmax,   &
        &        ntotd, npan_tot, ncheb, rpan_intervall, ipan_intervall, irmdnew, alat )
  call timing_stop('calc rhoq')
  
  
  
  call timing_stop('total time')


end program test
!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!
