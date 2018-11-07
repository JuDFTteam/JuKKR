module mod_rhoq

#ifdef CPP_MPI
use mpi
#endif

implicit none

type type_rhoq
  
  ! impurity cluster
  integer :: Nscoef ! number of impurities in impurity cluster
  integer :: Nlayer ! number of different layers in impurity cluster
  integer, allocatable :: ilay_scoef(:) ! (Nscoef), layer index for all positions in impurity cluster
  double precision, allocatable :: r_scoef(:,:) ! (3,Nscoef), position vector to all positions in impurity cluster

  ! exclude cluster
  integer :: Nexcl ! number of position which are excluded (positions directly above impurity)
  double precision, allocatable :: r_excl(:,:) ! (3,Nexcl), position vector to all positions in exclude cluster
  
  ! kmesh
  integer :: nkpt, Nkx, Nky ! total number of kpoints, nkpt=Nkx*Nky
  double precision :: volbz ! Brillouin zone volume -> integration weight
  double precision, allocatable :: kpt(:,:) ! (3,nkpt), coordinates of kpoints
  double precision, allocatable :: volcub(:) ! (nkpt), volume of cube associated with kpoint

  ! reduce q-points fo this box
  double precision :: qbox(3)
  
  ! geometry etc.
  integer :: mu_0 ! layer index of probe position
  integer :: natyp ! number of atoms in host system
  double precision, allocatable :: r_basis(:,:) ! (3,natyp), real space positions of all atoms in host system
  double precision, allocatable :: L_i(:,:) ! (3,Nscoef+1), lattice vectors for all atoms in impurity cluster and the probe layer
  integer :: lmmaxso ! size in l,m(,s) (if SOC calculation) subblocks
  
  ! Green functions etc.
  logical :: Ghost_k_memsave ! logical switch which determines if Ghost_k is stored in a file or kept in memory
!   double complex, allocatable :: Ghost(:,:,:) ! (Nlayer,lmmaxso,lmmaxso) 
!   double complex, allocatable :: Ghost_k(:,:,:,:) ! (Nlayer,nkpt,lmmaxso,lmmaxso) 
  double complex, allocatable :: G0tauG0_excl(:,:,:) ! (2*Nexcl, lmmaxso, lmmaxso), precalculated maxtrix sum_jk G0ij.tau_jk.G0_ki where j,k in imp cluster and i in exclude cluster

  ! these are used later on (not included in Bcast with Narrays parameter)
  double complex, allocatable :: Dt(:,:)   ! (Nscoef*lmmaxso,Nscoef*lmmaxso) 
  double complex, allocatable :: Gimp(:,:) ! (Nscoef*lmmaxso,Nscoef*lmmaxso) 
  double complex, allocatable :: tau(:,:)  ! (Nscoef*lmmaxso,Nscoef*lmmaxso) 

  ! logical switches
  logical :: exclude_only ! flag to determin if dry run (i.e. without q-loop) should be performed to get only C_M contribution
  
end type type_rhoq


type (type_rhoq), save :: t_rhoq


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef CPP_MPI
subroutine bcast_scalars_rhoq(t_rhoq)
  ! broadcast scalars in t_rhoq over all ranks
  use mpi
  use mod_mympi, only: myrank, master, nranks
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  ! local variables
  integer :: ierr

  call MPI_Bcast(t_rhoq%Nscoef, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%Nlayer, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%nkpt, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%Nkx, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%Nky, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%mu_0, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%natyp, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%lmmaxso, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%volbz, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%Ghost_k_memsave, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%exclude_only, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, ierr)

end subroutine bcast_scalars_rhoq
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef CPP_MPI
subroutine bcast_arrays_rhoq(t_rhoq)
  ! broadcast arrays in t_rhoq over all ranks
  ! has to be done after bacst scalars and then init_t_rhoq on all other processes that are not the master
  use mpi
  use mod_mympi, only: myrank, master, nranks
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  ! local variables
  integer :: ierr,N

  call MPI_Bcast(t_rhoq%ilay_scoef, t_rhoq%Nscoef,       MPI_INTEGER,          master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%r_scoef,    3*t_rhoq%Nscoef,     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%kpt,        3*t_rhoq%nkpt,       MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%volcub,     t_rhoq%nkpt,         MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%r_basis,    3*t_rhoq%natyp,      MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%L_i,        3*(t_rhoq%Nscoef+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%qbox,       3,                   MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%r_excl,     3*t_rhoq%Nexcl,      MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

end subroutine bcast_arrays_rhoq
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine init_t_rhoq(t_rhoq)
  ! initialize allocatable arrays in t_rhoq
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer :: ierr, N
  

  ierr = 0
  !integer array
  if(.not.allocated(t_rhoq%ilay_scoef)) allocate(t_rhoq%ilay_scoef(t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating ilay_scoef in rhoq'
  
  !double precision
  if(.not.allocated(t_rhoq%r_scoef)) allocate(t_rhoq%r_scoef(3,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating r_scoef in rhoq'
  if(.not.allocated(t_rhoq%kpt)) allocate(t_rhoq%kpt(3,t_rhoq%nkpt), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating kpt in rhoq'
  if(.not.allocated(t_rhoq%volcub)) allocate(t_rhoq%volcub(t_rhoq%nkpt), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating volcub in rhoq'
  if(.not.allocated(t_rhoq%r_basis)) allocate(t_rhoq%r_basis(3,t_rhoq%natyp), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating r_basis in rhoq'
  if(.not.allocated(t_rhoq%L_i)) allocate(t_rhoq%L_i(3,t_rhoq%Nscoef+1), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating L_i in rhoq'

  N = t_rhoq%Nscoef*t_rhoq%lmmaxso
  
  !double complex
  if(.not.allocated(t_rhoq%Dt)) allocate(t_rhoq%Dt(N, N), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Dt in rhoq'
  if(.not.allocated(t_rhoq%Gimp)) allocate(t_rhoq%Gimp(N, N), stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating Gimp in rhoq'
  if(.not.allocated(t_rhoq%tau)) allocate(t_rhoq%tau(N, N) , stat=ierr)
  if(ierr/=0) stop '[init_t_rhoq] error allocating tau in rhoq'
  
!   if(.not.allocated(t_rhoq%Ghost)) allocate(t_rhoq%Ghost(t_rhoq%lmmaxso,t_rhoq%lmmaxso,t_rhoq%Nlayer), stat=ierr)
!   if(ierr/=0) stop '[init_t_rhoq] error allocating Nscoef in rhoq'
!   if(.not.allocated(t_rhoq%Ghost_k)) allocate(t_rhoq%Ghost_k(t_rhoq%lmmaxso,t_rhoq%lmmaxso,t_rhoq%Nlayer,t_rhoq%nkpt), stat=ierr)  

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
  use mod_ioinput, only: ioinput
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
    call ioinput('NAEZ            ',uio,1,7,ier)
    if (ier.eq.0) then
      read (unit=uio,fmt=*) t_rhoq%natyp
      write(*,*) 'number of layers in system (natyp)= ',t_rhoq%natyp
    else
      stop '[read_input_rhoq] error neither "NATYP" nor "NAEZ" found in inputcard'
    endif
  endif

  call ioinput('Ghost_k_memsave ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) t_rhoq%Ghost_k_memsave
    write(*,*) 'save Ghost in memory? ',t_rhoq%Ghost_k_memsave
  else
    write(*,*) '[read_input_rhoq] "Ghost_k_memsave" not found in inputcard, use default value: F'
    t_rhoq%Ghost_k_memsave = .false.
  endif

  call ioinput('rhoq_qbox       ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) t_rhoq%qbox
    write(*,*) 'qbox? ',t_rhoq%qbox
  else
    t_rhoq%qbox(1) = 0.8d0
    t_rhoq%qbox(2) = 0.8d0
    t_rhoq%qbox(3) = 0.0d0
    write(*,*) 'rhoq_qbox keyword not found, taking default value for qbox:', t_rhoq%qbox
  endif

  call ioinput('exclude_only    ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) t_rhoq%exclude_only
    write(*,*) 'Compute exclude region only? ', t_rhoq%exclude_only
  else
    t_rhoq%exclude_only = .false.
    write(*,*) 'exclude_only keyword not found, taking default value:', t_rhoq%exclude_only
  endif
  
end subroutine read_input_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine save_kmesh_rhoq(t_rhoq,nkpt,kpt,volcub,volbz, Nkx, Nky)
  ! save kmesh information in t_rhoq
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: nkpt, Nkx, Nky ! total number of kpoints
  double precision, intent(in) :: volbz ! Brillouin zone volume -> integration weight
  double precision, intent(in) :: kpt(3,nkpt) ! coordinates in reciprocal space of the kpoints
  double precision, intent(in) :: volcub(nkpt) ! volume of kpoint cube
  !local
  integer :: ierr
  
  t_rhoq%nkpt = nkpt
  t_rhoq%Nkx = Nkx
  t_rhoq%Nky = Nky
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
  double precision :: Chi_mu0(3), Chi_cls(3) ! position in host system according to mu_orig and mu_cls (Chi_mu vector)
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
    
  ! find mu_orig
  ilayer = 0
  irun = 1
  do while(irun==1 .and. ilayer<t_rhoq%natyp)
    ilayer = ilayer+1
    if(dsqrt(sum(r_basis(:,ilayer)**2))<eps) irun = 0
  end do
  mu_orig = ilayer
  
  ! find mu_cls
  ilayer = 0
  irun = 1
  do while(irun==1)
    ilayer = ilayer+1
    if(dsqrt(sum(rcls_i(:,ilayer)**2))<eps) irun = 0
  end do
  mu_cls = t_rhoq%ilay_scoef(ilayer)
  
  ! set Chi_mu0 and Chi_cls from mu_orig and mu_cls (i.e. tip layer and center of imp cluster)
  Chi_mu0(:) = r_basis(:,t_rhoq%mu_0)
  Chi_cls(:) = r_basis(:,mu_cls)
  
  ! L_i = r^cls_i 
  do ilayer=1,t_rhoq%Nscoef
    L_i(:,ilayer) = rcls_i(:,ilayer)+Chi_cls(:)-r_basis(:,t_rhoq%ilay_scoef(ilayer))
    write(*,'(A,I,3F)') 'L_i',ilayer, L_i(:,ilayer)
  end do
  
  ! save result and exit
  ierr = 0
  if(.not.allocated(t_rhoq%L_i)) allocate(t_rhoq%L_i(3,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[get_L_vecs_rhoq] error allocating L_i in rhoq'
  t_rhoq%L_i(:,:) = L_i(:,:)
  deallocate(rcls_i, r_basis, L_i)

end subroutine get_L_vecs_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine start_excl(t_rhoq)
  ! prepare stuff for exclude cluster:
  ! 1. read in cluster info and host Green function for exclude cluster
  ! 2. precalculate sum_jk G0ij(E).tau_jk.G0ki(E) (used in C_M construction)
  ! output: fills t_rhoq%Nexcl, t_rhoq%r_excl, t_rhoq%G0tauG0_excl
  !         for later use in q-loop
  use mod_mympi, only: myrank, master
#ifdef CPP_MPI
  use mpi
#endif
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  !local variables
  integer natomimp, iatom, itmp1, itmp2, icls1, icls2, lm1, lm2, N, ienergy
  integer :: ierr, N1, N2, irun, ilayer, mu_cls
  integer, allocatable :: ilay_excl(:)
  double complex, allocatable :: temp(:,:),temp2(:,:)
  double precision, allocatable :: ratomimp(:,:)
  double complex, parameter :: C0=(0.0D0,0.0D0), C1=(1.0D0,0.0D0), CI=(0d0,1d0)
  double complex, allocatable :: Gll0(:,:)     !(ncls*lmmaxso, ncls*lmmaxso)
  double complex, allocatable :: Gll_3(:,:,:)  !(ncls*lmmaxso, ncls*lmmaxso,3)
  double complex, allocatable :: dGll(:,:)     !(ncls*lmmaxso, ncls*lmmaxso)
  double complex, allocatable :: d2Gll(:,:)    !(ncls*lmmaxso, ncls*lmmaxso)
  double complex, allocatable :: G0ijexcl(:,:) ! (Nexcl*lmmaxso,Nscoef*lmmaxso), real-space GF between position in exclude region and position in impurity cluster
  double complex, allocatable :: G0jiexcl(:,:) ! (Nexcl*lmmaxso,Nscoef*lmmaxso), like G0ijexcl but propagator from imp cluster to exclude region
  double precision :: dEimag
  double complex :: energy(3), dE1, dE2, ctmp1, ctmp2
  double precision, parameter :: eps = 1.0D-6 ! small epsilon -> threshold

  logical :: l_exist
  double precision :: r_offset(3)



  ! read in etc on master rank only and distribute result
  if(myrank==master) then

     open(unit=32452345,file='scoef_excl',iostat=ierr)
     if (ierr/=0) stop '[read_excl] file "scoef_excl" not found'
     read(32452345,*) natomimp
     allocate(ratomimp(3,natomimp), ilay_excl(natomimp))
     do iatom=1,natomimp
       read(32452345,*) ratomimp(:,iatom), ilay_excl(iatom)
       if(iatom>t_rhoq%Nscoef) write(*,'(A,I,A,3F,3i)') 'exclude cluster ',iatom,' :',ratomimp(:,iatom), natomimp, t_rhoq%Nscoef, ilay_excl(iatom)
     end do
     close(32452345)
     
     ! save Nexcl in t_rhoq
     t_rhoq%Nexcl = natomimp - t_rhoq%Nscoef
     
     ! allocate r_excl array from t_rhoq and save information in there
     ierr = 0
     if(.not.allocated(t_rhoq%r_excl)) allocate(t_rhoq%r_excl(3,t_rhoq%Nexcl), stat=ierr)
     if(ierr/=0) stop '[read_excl] error allocating r_excl in read_excl'
     
     t_rhoq%r_excl = ratomimp(:,t_rhoq%Nscoef+1:t_rhoq%Nscoef+t_rhoq%Nexcl)

     ! do the same change in reference point as for r_scoef in L_i determination (see get_L_vecs_rhoq subroutine for more details)
     ilayer = 0
     irun = 1
     do while(irun==1)
       ilayer = ilayer+1
       if(dsqrt(sum(t_rhoq%r_scoef(:,ilayer)**2))<eps) irun = 0
     end do
     mu_cls = t_rhoq%ilay_scoef(ilayer)

     inquire(file='r_offset.dat', exist=l_exist)
     if(l_exist) then
         write(*,*) 'found r_offset.dat file'
         open(unit=1283, file='r_offset.dat', form='formatted')
         read(1283, *) r_offset(1), r_offset(2), r_offset(3)
         close(1283)
         write(*,*) r_offset
     else
         r_offset(:) = 0.
     end if
     
     do ilayer=1,t_rhoq%Nexcl
       !write(*,'(I,9F)') ilayer, t_rhoq%r_excl(:,ilayer), t_rhoq%r_basis(:,mu_cls), t_rhoq%r_basis(:,ilay_excl(ilayer+t_rhoq%Nscoef))
       t_rhoq%r_excl(:,ilayer) = t_rhoq%r_excl(:,ilayer) + r_offset -t_rhoq%r_basis(:,ilay_excl(ilayer+t_rhoq%Nscoef))  +t_rhoq%r_basis(:,mu_cls)
       !t_rhoq%r_excl(:,ilayer) = t_rhoq%r_excl(:,ilayer) + r_offset! -t_rhoq%r_basis(:,ilay_excl(ilayer+t_rhoq%Nscoef))  !+t_rhoq%r_basis(:,mu_cls)-
       write(*,'(I,3F)') ilayer, t_rhoq%r_excl(:,ilayer)
     end do
     write(*,*) t_rhoq%r_basis(:,mu_cls), r_offset
     
     ! deallocate unused arrays
     deallocate(ratomimp, ilay_excl)
     
     
     ! now read in green_host_excl file and fill G0ijexcl and G0jiexcl
     N = t_rhoq%lmmaxso*(t_rhoq%Nscoef+t_rhoq%Nexcl)
     allocate(Gll0(N,N), Gll_3(N,N,3), dGll(N,N), d2Gll(N,N), stat=ierr)
     if(ierr/=0) stop '[rhoq] Error allocating Gll0 etc in read_excl'
     
     ! Read the gll_3 (for the three energies)
     open(unit=1283, file='green_host_excl', form='formatted', action='read')
     write(*,*) 'reading file green_host_excl'
     
     do ienergy=1,3
     
       write(*,*) ' reading energy',ienergy
       read(1283,"(2(e17.9,X))") energy(ienergy)
     
       do icls1=1,t_rhoq%lmmaxso*(t_rhoq%Nscoef+t_rhoq%Nexcl)
         do icls2=1,t_rhoq%lmmaxso*(t_rhoq%Nscoef+t_rhoq%Nexcl)
           read(1283,"((2I5),(2e17.9))") lm1, lm2, Gll_3(icls2, icls1, ienergy)
         end do
       end do
     
     end do!ienergy
     
     close(1283)
     
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
     
     ! store in G0ijexcl and G0jiexcl arrays
     if(.not.allocated(G0ijexcl)) allocate(G0ijexcl(t_rhoq%Nexcl*t_rhoq%lmmaxso,t_rhoq%Nscoef*t_rhoq%lmmaxso), stat=ierr)
     if(ierr/=0) stop '[read_excl] error allocating G0ijexcl'
     if(.not.allocated(G0jiexcl)) allocate(G0jiexcl(t_rhoq%Nscoef*t_rhoq%lmmaxso,t_rhoq%Nexcl*t_rhoq%lmmaxso), stat=ierr)
     if(ierr/=0) stop '[read_excl] error allocating G0jiexcl'
     G0ijexcl = C0
     G0jiexcl = C0
     write(*,*) 'fill G0ijexcl'
     do icls1=1,t_rhoq%Nexcl
       do icls2=1,t_rhoq%Nscoef
         itmp1 = (icls1-1)*t_rhoq%lmmaxso+1
         itmp2 = (icls2-1)*t_rhoq%lmmaxso+1
         N1 = (icls1+t_rhoq%Nscoef-1)*t_rhoq%lmmaxso+1
         N2 = (icls2-1)*t_rhoq%lmmaxso+1
         !write(*,*) 'G0ijexcl',icls1,icls2,itmp1,itmp2
         G0ijexcl(itmp1:itmp1-1+t_rhoq%lmmaxso,itmp2:itmp2-1+t_rhoq%lmmaxso) = Gll0(N1:N1+t_rhoq%lmmaxso-1,N2:N2+t_rhoq%lmmaxso-1)
       end do !icls2
     end do !icls1
     write(*,*) 'fill G0jiexcl'
     do icls1=1,t_rhoq%Nscoef
       do icls2=1,t_rhoq%Nexcl
         itmp1 = (icls1-1)*t_rhoq%lmmaxso+1
         itmp2 = (icls2-1)*t_rhoq%lmmaxso+1
         N1 = (icls1-1)*t_rhoq%lmmaxso+1
         N2 = (icls2+t_rhoq%Nscoef-1)*t_rhoq%lmmaxso+1
         G0jiexcl(itmp1:itmp1-1+t_rhoq%lmmaxso,itmp2:itmp2-1+t_rhoq%lmmaxso) = Gll0(N1:N1+t_rhoq%lmmaxso-1,N2:N2+t_rhoq%lmmaxso-1)
       end do !icls2
     end do !icls1
     
     ! deallocate temporary arrays
     deallocate(dGll, d2Gll, Gll_3, Gll0)
     
     
     ! Now calculate sum_jk [ G0_ij(E) tau_jk G0_ki(E) ]
     ! this prefactor is used in C_M calculation
     
     if(.not.allocated(t_rhoq%G0tauG0_excl)) allocate(t_rhoq%G0tauG0_excl(2*t_rhoq%Nexcl,t_rhoq%lmmaxso,t_rhoq%lmmaxso)) 
     if(ierr/=0) stop '[precalc_G0tauG0] Error allocating t_rhoq%G0tauG0_excl'
     
     ! matrix sizes
     N1 = t_rhoq%Nscoef*t_rhoq%lmmaxso
     N2 = t_rhoq%Nexcl*t_rhoq%lmmaxso
     
     ! now do matrix matrix operations:
     ! calculate first temp = Gij.tau
     allocate(temp(N2,N1),temp2(N2,N2), stat=ierr)
     if(ierr/=0) stop '[precalc_G0tauG0_excl] error allocating temp,temp2'
     write(*,*) 'calc G0tauG01',shape(G0ijexcl),shape(t_rhoq%tau),shape(temp)
     write(*,*) 'dims',N2,N1, N1,N1, N2,N1
     ! temp = G0ij.tau
     call ZGEMM('n','n',N2,N1,N1,C1,G0ijexcl,N2,t_rhoq%tau,N1,C0,temp,N2)
     
     ! sum_jk [...] = Gij.tau.Gji
     !              = temp.Gji
     write(*,*) 'calc G0tauG02',shape(temp),shape(G0jiexcl),shape(temp2),shape(t_rhoq%G0tauG0_excl)
     write(*,*) 'dims',N2,N1, N1,N2, N2,N2
     call ZGEMM('n','n',N2,N2,N1,C1,temp,N2,G0jiexcl,N1,C0,temp2,N2)

     ! now restructure (take only diagonal lm-blocks)
     do icls1=1,t_rhoq%Nexcl
        N1 = (icls1-1)*t_rhoq%lmmaxso+1
        N2 = icls1*t_rhoq%lmmaxso
        write(*,*) 'restruc1', icls1, n1,n2
        t_rhoq%G0tauG0_excl(icls1,:,:) = temp2(N1:N2,N1:N2)
        !do lm1 = 1, t_rhoq%lmmaxso
        ! do lm2 = 1, t_rhoq%lmmaxso
        !  write(7894521,'(3i5,2f14.7)') icls1,lm1, lm2, t_rhoq%G0tauG0_excl(icls1,lm1,lm2)
        ! end do
        !end do
     end do


     ! matrix sizes
     N1 = t_rhoq%Nscoef*t_rhoq%lmmaxso
     N2 = t_rhoq%Nexcl*t_rhoq%lmmaxso

     ! temp = G0ij^*.tau^*
     call ZGEMM('n','n',N2,N1,N1,C1,dconjg(G0ijexcl),N2,dconjg(t_rhoq%tau),N1,C0,temp,N2)
     
     ! sum_jk [...] = Gij^*.tau^*.Gji^*
     !              = temp.Gji^*
     call ZGEMM('n','n',N2,N2,N1,C1,temp,N2,dconjg(G0jiexcl),N1,C0,temp2,N2)

     ! now restructure (take only diagonal lm-blocks)
     do icls1=1,t_rhoq%Nexcl
        N1 = (icls1-1)*t_rhoq%lmmaxso+1
        N2 = icls1*t_rhoq%lmmaxso
        write(*,*) 'restruc2', icls1, n1,n2
        t_rhoq%G0tauG0_excl(icls1+t_rhoq%Nexcl,:,:) = temp2(N1:N2,N1:N2)
        !do lm1 = 1, t_rhoq%lmmaxso
        ! do lm2 = 1, t_rhoq%lmmaxso
        !  write(7894521,'(3i5,2f14.7)') icls1+t_rhoq%Nexcl,lm1, lm2, t_rhoq%G0tauG0_excl(icls1+t_rhoq%Nexcl,lm1,lm2)
        ! end do
        !end do
     end do
     
     deallocate(temp,temp2)

  end if ! myrank==master

#ifdef CPP_MPI
  ! allocate for myrank/=master and broadcast data from master
  call MPI_Bcast(t_rhoq%Nexcl, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  if(.not.allocated(t_rhoq%G0tauG0_excl)) allocate(t_rhoq%G0tauG0_excl(t_rhoq%Nexcl*2,t_rhoq%lmmaxso,t_rhoq%lmmaxso)) 
  if(ierr/=0) stop '[precalc_G0tauG0] Error allocating t_rhoq%G0tauG0_excl'
  if(.not.allocated(t_rhoq%r_excl)) allocate(t_rhoq%r_excl(3,t_rhoq%Nexcl), stat=ierr)
  if(ierr/=0) stop '[read_excl] error allocating r_excl in read_excl'
  call MPI_Bcast(t_rhoq%r_excl, 3*t_rhoq%Nexcl, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_rhoq%G0tauG0_excl, 2*t_rhoq%Nexcl*t_rhoq%lmmaxso*t_rhoq%lmmaxso, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
#endif

end subroutine start_excl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine read_Dt_Gimp_rhoq(t_rhoq, lmmaxso, ncls)
  ! read in DTMTRX and GMATLL_GES, provided by zulapi code
  use mod_mympi, only: myrank, master
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: lmmaxso, ncls
#ifdef CPP_MPI
  integer :: ierr
#endif

  if(myrank==master) then
    ! read in DTMTRX and GMATLL_GES on master
    call read_DTMTRX( t_rhoq%Dt, lmmaxso, ncls)
    call read_green_ll(ncls, lmmaxso, t_rhoq%Gimp)
  endif

#ifdef CPP_MPI
  ! broadcast Dt and Gimp

  call MPI_Bcast(t_rhoq%Dt, (ncls*lmmaxso)**2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  if(ierr/=0) stop '[read_Dt_Gimp_rhoq] Error Bcast Dt'
  
  call MPI_Bcast(t_rhoq%Gimp, (ncls*lmmaxso)**2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  if(ierr/=0) stop '[read_Dt_Gimp_rhoq] Error Bcast Gimp'
#endif

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
  !     = tau + tau.temp using tau=Dt and temp=Gimp.Dt
  t_rhoq%tau = t_rhoq%Dt
  call ZGEMM('n','n',N,N,N,C1,t_rhoq%Dt,N,temp,N,C1,t_rhoq%tau,N)
  
  
  deallocate(temp)

end subroutine calc_tau_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_Q_mu_rhoq(lmax, ntotd, npan_tot, &
            & nsra, lmmaxso, Rll, Rllleft, ipan_intervall, &
            & irmdnew, trq_of_r, ncleb_max )
  ! calculate prefactor needed in rhoq calculation:
  ! Q^mu_Lambda,Lamabda' = Int d^3r Tr[ R^mu_Lmabda(\vec{r}) R^left,mu_Lmabda'(\vec{r}) ] where R^left denotes the left solution
  use mod_gaunt2, only: gaunt2
  implicit none
  integer, intent(in) :: lmax, ntotd, npan_tot
  integer, intent(inout) :: nsra, lmmaxso, irmdnew ! left/right solutions (wave functions)
!   integer, intent(in) :: nsra, lmmaxso, irmdnew ! left/right solutions (wave functions)
  double complex, intent(in) :: Rll(nsra*lmmaxso,lmmaxso,irmdnew), Rllleft(nsra*lmmaxso,lmmaxso,irmdnew) ! wave functions in new mesh
  integer, intent(in) :: ipan_intervall(0:ntotd)
  double complex, allocatable, intent(out) :: trq_of_r(:,:,:,:) !(lmmaxso,lmmaxso, 100, irmdnew)
  integer, intent(out) :: ncleb_max
  
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
  
  allocate( WG(4*LMAX), YRG(4*LMAX,0:4*LMAX,0:4*LMAX), CLEB(NCLEB), ICLEB(NCLEB,3), JEND((LPOT+1)**2,0:LMAX,0:LMAX), LOFLM((2*LPOT+1)**2), stat=ierr )
  if(ierr/=0) stop '[calc_Q_mu_rhoq] Error allocating icleb etc for gaunt'
  
  ! find wg, yrg
  CALL GAUNT2( WG, YRG, 4*LMAX )
  
  ! compute gaunt coefficients (normal ones with lmax_1 = lmax_2 = lmax)
  CALL GAUNT_new( LMAX, LMAX, LPOT, WG, YRG, CLEB, LOFLM, ICLEB, IEND, JEND, NCLEB, LMAX, (LPOT+1)**2 )
  open(777888, file='icleb.3')
  ncleb_max = (2*lmax+1)**2 !LM cutoff for trq_of_r, (49 for lmax=3); cutoff of trq_of_r due to trq_of_r(L) = R_L' * R_L'' * C_L'L''L (R_L has cutoff of lmax)
  allocate(trq_of_r(lmmaxso,lmmaxso, ncleb_max, irmdnew))
  write(777888, *) ncleb_max
  write(777888, *) icleb(:,3)
  close(777888)
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
 !          do ir=imt1+1,irmdnew
          do ir=1,irmdnew
            ! sum of Rll*Rllleft with Gaunt coefficients, with trace over big and small component (lmmaxso shift in first index of Rll, Rllleft)
            trq_of_r(lm01, lm02, lm3, ir) = trq_of_r(lm01, lm02, lm3, ir) + cleb(j) * ( Rll(lm01, lm1, ir)*Rllleft(lm02, lm2, ir) + Rll(lm01+lmmaxso, lm1, ir)*Rllleft(lm02+lmmaxso, lm2, ir) )
          enddo ! ir
        enddo ! j -> sum of Gaunt coefficients
         
        enddo ! lm02
    enddo ! lm01
  enddo ! ispin -> sum over spins
        
        
  deallocate( WG, YRG, CLEB, ICLEB, JEND, LOFLM, stat=ierr )
  if(ierr/=0) stop '[calc_Q_mu_rhoq] Error deallocating gaunt arrays'
  
  
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
  double precision :: rbt(2,2), irbt(2,2), rq(2), kpt(2), det, qv(2)!, minimum(3)
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
  do i=1,t_rhoq%Nkpt
    kpt(:) = t_rhoq%kpt(:,i)
    if(dabs(dsqrt((kpt(1)-qv(1))**2+(kpt(2)-qv(2))**2))<eps) then
      iq=i
    end if
  end do
  
  !check if iq was found
  if(iq<1) write(*,*) qv!, minimum
  if(iq<1) stop '[red_Q] Error: no iq found'
  
  
end subroutine red_Q


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_tinv(t_rhoq, tinv)
  use mod_mympi, only: myrank, master
  implicit none
  ! parameter
  double complex, parameter :: C0=(0.0d0, 0.0d0), C1=(1.0d0, 0.0d0)

  ! input
  type(type_rhoq), intent(in) :: t_rhoq
  ! output
  double complex, intent(out) :: tinv(:,:,:) !allocated outside with (N,N,Nscoef+1)
  ! local
  integer :: nref, ie, ielast, ierr, i, irec, lrec, N, lm1
  integer, allocatable :: ipvt(:), refpot(:)
  double complex, allocatable :: tll(:,:), trefll(:,:,:)

  !matrix dimensions (depending on lmax)
  N = t_rhoq%lmmaxso
  
  !read stuff in master only and pass to other ranks (is fastest)
  if(myrank==master) then
    open(9999, file='refinfo', form='formatted')
    read(9999,'(1I9)') NREF
    allocate(refpot(1:t_rhoq%NATYP), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] ERROR in allocating refpot'
    read(9999, '(1000I9)') refpot(1:t_rhoq%NATYP)
    close(9999)
  else
    allocate(refpot(1:t_rhoq%NATYP), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] ERROR in allocating refpot'
  endif
  
#ifdef CPP_MPI
  call MPI_Bcast(refpot, t_rhoq%NATYP, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[calc_rhoq] error brodcasting scalars in t_rhoq'
#endif
  
  if(myrank==master) then
    allocate(trefll(N,N, NREF), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] error allocating trefll'
    open(9999, file='tref', access='direct', recl=4*N**2 )
    ie = 2
    ielast = 3
    do i = 1,nref
       irec = ie + ielast*(i-1)
       read(9999, rec=irec) trefll(1:N,1:N,i)
    enddo ! i
    close(9999)
  endif
  
#ifdef CPP_MPI
  call MPI_Bcast(NREF, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  if(myrank/=master) then
    allocate(trefll(N,N, NREF), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] ERROR in allocating trefll'
  end if
  call MPI_Bcast(trefll, N*N*NREF, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[calc_rhoq] error brodcasting scalars in t_rhoq'
#endif
  
  if(myrank==master) then
    !temporary array
    allocate(ipvt(N), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] error allocating ipvt'
    ipvt = 0
    allocate(tll(N,N), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] error allocating tll'
    tll = C0

    ! read tmat and compute tinv = tmat^-1
    lrec = 4*N**2
    open(9999, file='tmat', access='direct', recl=lrec, form='unformatted')
  
    do i=1,t_rhoq%Nscoef
  
      irec = ie+ielast*(t_rhoq%ilay_scoef(i)-1) ! layer index 
      read(9999,rec=irec) tll(1:N,1:N)
      ! tau0_k is defined with respect to t_ref, thus: t-t_ref
      tll = tll - trefll(:,:,refpot(t_rhoq%ilay_scoef(i)))

      ! initialize tinv
      tinv(1:N,1:N,i)=C0
      do lm1=1,N
        tinv(lm1,lm1,i)=C1
      enddo
    
      ! LU decomposition to invert tll: tinv = tll^-1
      call zgetrf(N,N,tll,N,ipvt,ierr)
      call zgetrs('n',N,N,tll,N,ipvt,tinv(:,:,i),N,ierr)
    
    end do ! i
  
    irec =  ie+ielast*(t_rhoq%mu_0-1) ! layer index for probe layer mu_0
    read(9999,rec=irec) tll(1:N,1:N)
    ! t-t_ref
    tll = tll - trefll(:,:,refpot(t_rhoq%mu_0))

    ! initialize tinv
    tinv(1:N,1:N,t_rhoq%Nscoef+1)=C0
    do lm1=1,N
      tinv(lm1,lm1,t_rhoq%Nscoef+1)=C1
    enddo
    ! LU decomposition to invert tll
    call zgetrf(N,N,tll,N,ipvt,ierr)
    call zgetrs('n',N,N,tll,N,ipvt,tinv(:,:,t_rhoq%Nscoef+1),N,ierr)

    ! close file 'tmat'
    close(9999)
  
    deallocate(tll, ipvt, trefll, stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] error deallocating tll etc.'
  end if !myrank==master
  
#ifdef CPP_MPI
  call MPI_Bcast(tinv, N*N*(t_rhoq%Nscoef+1), MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[calc_rhoq] error brodcasting scalars in t_rhoq'
#endif
  
end subroutine get_tinv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_rhoq(t_rhoq, lmmaxso, Nkp, trq_of_r, recbv, lmax,   &
        &        ntotd, npan_tot, ncheb, rpan_intervall, ipan_intervall, irmdnew, alat, ncleb_max )
  ! calculate Delta rho^mu(q) = Int{ Delta rho(q, X_mu+r), dr }
  !                           = -1/(2*pi*i) Tr[ Q^mu Int{ Ghost(k) tau Ghost(k+q),dk } - Q^mu,* Int{ Ghost^*(k) tau^* Ghost^*(k-q),dk }]
  use omp_lib
  use mod_timing
  use IFPORT ! random numbers
#ifdef CPP_MPI
  use mod_mympi, only: master, myrank, nranks, distribute_linear_on_tasks
#else
  use mod_mympi, only: master, myrank, nranks
#endif
  use mod_intcheb_cell, only: intcheb_cell
  use mod_gaunt2, only: gaunt2
  implicit none
  type(type_rhoq), intent(inout) :: t_rhoq
  integer, intent(in) :: lmmaxso, Nkp, irmdnew, lmax, ncleb_max
  double complex, intent(inout) :: trq_of_r(lmmaxso,lmmaxso, ncleb_max, irmdnew)
  double precision, intent(in) :: recbv(3,3), alat
  
  integer :: ntotd, npan_tot, ncheb, npan_log, npan_eq
  integer npan_lognew,npan_eqnew
  double precision, intent(in) :: rpan_intervall(0:ntotd)
  integer, intent(in) :: ipan_intervall(0:ntotd)
  
  ! local
  integer :: lrec, i, irec, k, mu, mu_i, ie, N, q, j, kpq, ierr, lm1, ix, nthreads, Ni, Nj, Nqpt, lm01, lm02, ir, l1, m1, imt1, lm2, lm3, ifun, imin, ikx, iky, Nx, Ny
  
  integer, save :: mythread
  !$omp threadprivate(mythread)
  
  double complex, allocatable :: tinv(:,:,:), tau0_k(:,:), G0ij_k(:,:,:,:), G0ji_k(:,:,:,:), tau(:,:,:,:), tmpG0(:,:), rhoq(:), rhoq_excl(:)
  
  double complex, allocatable, save :: tmp(:,:), exG0_tmp(:,:), jl(:), tmpsum1(:,:), tmpsum2(:,:), eiqr_lm(:,:), cYlm(:), qint(:,:,:), q_mu(:,:)
  !$omp threadprivate(tmp, exG0_tmp, jl, tmpsum1, tmpsum2, eiqr_lm, cYlm, qint, q_mu)
  
  double complex :: tr_tmpsum1, tr_tmpsum2, kweight, Z
  double precision, allocatable :: kpt(:,:), L_i(:,:), Qvec(:,:), rnew(:), qvec_tmp(:,:)
  
  double precision, allocatable, save :: Ylm(:)
  !$omp threadprivate(Ylm)
  
  integer, allocatable :: qvec_index_tmp(:,:), qvec_index(:,:)
  double precision :: QdotL, tmpk(3), tmpr(3), Rq, costheta, phi, box(3)
  integer, allocatable :: ifunm(:), lmsp(:), ntcell(:) ! shape functions
  double precision, allocatable :: thetasnew(:,:) ! shape functions in Chebychev mesh (new mesh)
  logical, allocatable :: kmask(:)
  integer :: kmask_tmp
  logical :: no_kmask_given

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
  
#ifdef CPP_MPI
  integer, allocatable :: ntot_pT(:), ioff_pT(:)
#endif
  integer :: q_start, q_end
  integer, allocatable :: q_rand(:)
  integer, parameter :: Nrandomize=1000
  integer :: tmp_int, iq
  double precision :: rand_num
  
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
  allocate(tinv(N,N,t_rhoq%Nscoef+1), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tinv'
  allocate(G0ij_k(N,N,t_rhoq%Nscoef,Nkp), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating G0ij_k'
  allocate(G0ji_k(N,N,t_rhoq%Nscoef,Nkp), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating G0ji_k'
  
  G0ij_k = C0
  G0ji_k = C0
  
  !tinv = (t-tref)^-1
  call get_tinv(t_rhoq, tinv)
  
  
  
  call timing_start('calc rhoq - comp tau')
  if(myrank==master) then
  
    ! find imin
    imin = 1000
    do i=1,t_rhoq%Nscoef
      if(t_rhoq%ilay_scoef(i)<imin) imin = t_rhoq%ilay_scoef(i)
      if(myrank==master) write(*,*) 'find imin' ,imin,t_rhoq%ilay_scoef(i), t_rhoq%Nlayer
    end do
  
    ! open tau0_k file, written in writegreen of zulapi code
    lrec = 4*(N**2+1)
!     open(9999, file='tau0_k', access='direct', form='unformatted', recl=lrec)
  
    ! allocate temporary arrays
    allocate(tau0_k(N,N), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] error allocating tau0_k'
    allocate(tmpG0(N,N), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] error allocating tmpG0'
  
    !counter: how many kpts are inside/outside
    lm1 = 0
    lm2 = 0

    !FSqdos_rhoq.txt writeout
    open(776655, file='FSqdos_tau_rhoq.txt', form='formatted')
    write(776655, '(A)') '# first half is from tau_0ij_k, second half tau_0ji_k'
    write(776655, '(A)') '# ik, iscoef, kx, ky, tr(tau_0ij_k(:,:,iscoef, ik)'
    open(556677, file='FSqdos_rhoq.txt', form='formatted')
    write(556677, '(A)') '# first half is from G0ij_k, second half G0ji_k'
    write(556677, '(A)') '# ik, iscoef, kx, ky, tr(G0ij_k(:,:,iscoef, ik)'

    if(.not.t_rhoq%exclude_only) then

      inquire(file='kpts_mask.txt', exist=no_kmask_given)
      no_kmask_given = .not. no_kmask_given
      if(.not. no_kmask_given) open(8888, file='kpts_mask.txt', form='formatted')

      write(*,*) 'found kpts_mask.txt?', no_kmask_given

      ! calculate tau
      write(*,*) 'calculate G0_k from tau0_k'
      write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
      write(*,FMT=190) !beginning of statusbar
      do k=1,Nkp
        ! find exG0_i(k) vector of lm-blocks (i is one component)
      
        if(no_kmask_given) then
          kmask_tmp = 1
        else
          read(8888, *) kmask_tmp
        end if

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

          if (kmask_tmp>0 .or. no_kmask_given) then
             read(998899,'(10000ES15.7)') tmpk(1:2), tau0_k(1:N,1:N)
          else
             tau0_k(1:N,1:N) = C0
             read(998899, *) tmpk(1:2)
          end if

          q = k
      
      
          ! set kmask
          if( kmask_tmp==0 .or. (no_kmask_given .and. dsqrt(dreal(tau0_k(1,1))**2+dimag(tau0_k(1,1))**2)<eps) ) then
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
      
            !FStauFStauFStauFStauFStauFStauFStauFStauFStauFStauFS
            ! writeout trace of tau_0ij_k and tau_0ji_k
            tmpG0(1,1) = 0.0d0
            do ix=1,N
               tmpG0(1,1) = tmpG0(1,1) + tau0_k(ix, ix)
            end do
            write(776655, '(2i9,100ES15.7)') q, i, tmpk(1:2), tmpG0(1,1)/kmask_tmp
            !FStauFStauFStauFStauFStauFStauFStauFStauFStauFStauFS
      
          end if !(kmask(q))
          
          ! now mu_i,mu_0 element
          irec = ((t_rhoq%Nlayer-1)*2)*(k-1) + ((t_rhoq%Nlayer-1)*2)*Nkp*(ie-1-1)
          ix = mu_i - mu + 1
          irec = irec + ix
          if (kmask(q)) then
             read(998888,'(10000ES15.7)') tmpk(1:2), tau0_k(1:N,1:N)
          else
             tau0_k(1:N,1:N) = C0
             read(998888, *) tmpk(1:2)
          end if
          
          if(kmask(q)) then
      
            !find G0ji_k from tau0_k and tinv
            call calc_G0_k(tau0_k, tinv(1:N,1:N,i), tinv(1:N,1:N,t_rhoq%Nscoef+1), tmpG0(1:N,1:N), mu_i, t_rhoq%mu_0, N)
            G0ji_k(1:N,1:N,i,q) = tmpG0(1:N,1:N)
      
            !FStauFStauFStauFStauFStauFStauFStauFStauFStauFStauFS
            ! writeout trace of tau_0ij_k and tau_0ji_k
            tmpG0(1,1) = 0.0d0
            do ix=1,N
               tmpG0(1,1) = tmpG0(1,1) + tau0_k(ix, ix)
            end do
            write(776655, '(2i9,100ES15.7)') q, i, tmpk(1:2), tmpG0(1,1)/kmask_tmp
            !FStauFStauFStauFStauFStauFStauFStauFStauFStauFStauFS
      
            !FSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFS
            ! writeout trace of G0ij_k and G0ji_k
            tmpG0(1,1) = 0.0d0
            do ix=1,N
               tmpG0(1,1) = tmpG0(1,1) + G0ij_k(ix, ix, i, q)
            end do
            write(556677, '(2i9,100ES15.7)') q, i, tmpk(1:2), tmpG0(1,1)
            tmpG0(1,1) = 0.0d0
            do ix=1,N
               tmpG0(1,1) = tmpG0(1,1) + G0ji_k(ix, ix, i, q)
            end do
            write(556677, '(2i9,100ES15.7)') q, i, tmpk(1:2), tmpG0(1,1)
            !FSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFSFS
      
          end if !kmask
          
        end do ! i
        
        !update statusbar
        if(Nkp>=50.and.mod(k,Nkp/50)==0) write(*,FMT=200)
        
      end do ! k
      
      !FSqdos_rhoq.txt writeout
      close(556677)
      if(.not. no_kmask_given) close(8888)
      
!       close(998899)
!       close(998888)
      deallocate(tmpG0, tau0_k, stat=ierr)
      if(ierr/=0) stop '[calc_rhoq] error deallocating tmpG0 etc'
      
      write(*,*) !status bar
      write(*,*) 'kmask info (inside/outside):', lm1, lm2

    else ! t_rhoq%exclude_only
      ! set kmask to False to prevent rhoq calculation
      kmask(:) = .false.
    end if ! .not. t_rhoq%exclude_only
  
  endif !myrank==master

#ifdef CPP_MPI
  !if(myrank==master) write(*,*) 'bcast kmask', Nkp
  write(*,'(A,10I)') 'bcast kmask', myrank, Nkp
  call MPI_Bcast(kmask, Nkp, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[calc_rhoq] error brodcasting kmask in t_rhoq'
  !if(myrank==master) write(*,*) 'bcast G0ji_k', N, t_rhoq%Nscoef, shape(G0ji_k)
  write(*,'(A,10I)') 'bcast G0ji_k', myrank, N, t_rhoq%Nscoef, shape(G0ji_k)
  call MPI_Bcast(G0ji_k, N*N*t_rhoq%Nscoef*Nkp, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[calc_rhoq] error brodcasting G0ji_k in t_rhoq'
  !if(myrank==master) write(*,*) 'bcast G0ij_k', N, t_rhoq%Nscoef, shape(G0ij_k)
  write(*,'(A,10I)') 'bcast G0ij_k', myrank, N, t_rhoq%Nscoef, shape(G0ij_k)
  call MPI_Bcast(G0ij_k, N*N*t_rhoq%Nscoef*Nkp, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop '[calc_rhoq] error brodcasting G0ij_k in t_rhoq'
#endif
  
  
  ! take tau from t_rhoq and reshape it 
  allocate( tau(lmmaxso,lmmaxso,t_rhoq%Nscoef,t_rhoq%Nscoef), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tau imp'
  
  ! reshape tau
  !$omp parallel do default(shared) private(j,i,lm2,lm1)
  do j=1,t_rhoq%Nscoef
    do i=1,t_rhoq%Nscoef
      do lm2=1,lmmaxso
        do lm1=1,lmmaxso
          tau(lm1, lm2 ,i, j) = t_rhoq%tau(lm1+lmmaxso*(i-1),lm2+lmmaxso*(j-1))
        end do
      end do
    end do
  end do
  !$omp end parallel do
  
  if(myrank==master) write(*,*) 'done with tau'

  ! done calculating tau
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call timing_stop('calc rhoq - comp tau')
  
    
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! find second set of gaunt coefficients for sum_LL'L" [exp(-iqr)_L * Tr(R*Rleft)_L' * Theta)L" * C_LL'L"] which is then integrated radially
  call timing_start('calc rhoq - Gaunt')
  
  NCLEB =   ((lmax*2+1)**2 * (lmax+1)**2 )*100
  LMAX_1 = 6
  LMAX_2 = 18
  LPOT_2 = 2*(LMAX_1+LMAX_2)
  
  allocate( CLEB(NCLEB), ICLEB(NCLEB,3), JEND((LPOT_2+1)**2,0:LMAX_1,0:LMAX_2), LOFLM((2*LPOT_2+1)**2), stat=ierr )
  if(ierr/=0) stop '[calc_rhoq] Error allocating icleb etc for gaunt'
  
  if(myrank==master) then
  
    !temporary arrays
    allocate( WG(4*LMAX_2), YRG(4*LMAX_2,0:4*LMAX_2,0:4*LMAX_2), stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] Error allocating wg, yrg'
    
    ! first set wg and yrg arrays (input for gaunt_new)
    !          ( > |  > |     <    )
    CALL GAUNT2( WG, YRG, 4*LMAX_2)
    
    ! compute gaunt coefficients depending on two different cutoffs (lmax_1 and lamx_2)
    !             (    <  |    <  |   <   | < |  < |  >  |  >   |   >  |  >  |  >  |  >   |   <     |     <         )
    CALL GAUNT_new( LMAX_1, LMAX_2, LPOT_2, WG, YRG, CLEB, LOFLM, ICLEB, IEND, JEND, NCLEB, LMAX_2  , (LPOT_2+1)**2 )

  open(777888, file='icleb.32')
  ierr = 0
  do irid=1,iend
      if (icleb(irid,3)>ierr) ierr=icleb(irid,3)
  end do
  write(777888, *) ierr
  write(777888, *) icleb(:,3)
  close(777888)
    
    deallocate(WG, YRG, stat=ierr)
    if(ierr/=0) stop '[calc_rhoq] Error deallocating wg, yrg'
    
    if(myrank==master) write(*,*) 'done with gaunts'
    
  end if !myrank==master
#ifdef CPP_MPI
  call MPI_Bcast(IEND, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(CLEB, NCLEB, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ICLEB, NCLEB*3, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(JEND, (LPOT_2+1)**2*(LMAX_1+1)*(LMAX_2+1), MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(LOFLM, (2*LPOT_2+1)**2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
#endif
  call timing_stop('calc rhoq - Gaunt')
  ! done finding Gaunts for larger lmax values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
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
  allocate( rnew(irmdnew) ) ! new mesh
  
  if(myrank==master) then
  
    allocate( NTCELL(1:t_rhoq%NATYP) )
    allocate( thetas(irid,nfund) )
    
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
    
    
    if(myrank==master) write(*,*) 'done with shapes'
    
  end if ! myrank=master
#ifdef CPP_MPI
  !
  call MPI_Bcast(ifunm, (2*lpot_2+1)**2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  !arrays
  call MPI_Bcast(ifunm, (2*lpot_2+1)**2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(lmsp, (2*lpot_2+1)**2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(thetasnew, NTOTD*(NCHEB+1)*nfund, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(rnew, irmdnew, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
#endif
  
  ! done finding shapefunction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call timing_stop('calc rhoq - shape')
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !find q-mesh from k-kpoints
  
  ! parameters for bigger q mesh, define here up to second brillouin zone
  Ni = 2
  Nj = 2
  Nqpt = Ni*Nj*Nkp
  
  Nx = t_rhoq%Nkx
  Ny = t_rhoq%Nky
  if(Nx*Ny/=Nkp) stop '[calc_rhoq] Error: Nx*Ny/=Nkp'
  
  
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
  !$omp parallel do default(shared) private(k,i,j,irec,ikx, iky)
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
        qvec_index(1,irec) = ikx-1
        qvec_index(2,irec) = iky-1

      end do !j=1,Nj
    end do !i=1,Ni
    
  end do !k=1,Nkp
  !$omp end parallel do

  ! reduce qpts to box around Gamma
  box(1) = t_rhoq%qbox(1)   !0.6 !0.8   !0.50d0
  box(2) = t_rhoq%qbox(2)   !0.3 !0.8   !0.20d0
  box(3) = t_rhoq%qbox(3)   !0.0 !0.00d0 ! no kz-limit, use for inner radius cutoff
  
  if(myrank==master) write(*,*) 'reduce kpts to:', box(3), box(1)

  k = 0
  do q=1,Nqpt
    if((abs(qvec(1,q))-box(1)<=eps).and.(abs(qvec(2,q))-box(2)<=eps)) then
      if(.true.) then !(qvec(1,q)>=0).and.(qvec(2,q)>=0)) then
    !if((abs(qvec(1,q))-box(1)<=eps).and.(abs(qvec(2,q))-box(2)<=eps).and.(abs(qvec(3,q))-box(3)<=eps)) then
    !  if(dsqrt((qvec(1,q))**2+(qvec(2,q))**2+(qvec(3,q))**2)>box(3)) then
        k = k+1
        qvec_tmp(1:3,k) = qvec(1:3,q)
        qvec_index_tmp(1:2,k) = qvec_index(1:2,q)
      end if
    end if
  end do

  Nqpt = k
  if(myrank==master) write(*,*) 'red qvecs box',k, box

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


  allocate(rhoq(Nqpt),rhoq_excl(Nqpt), stat=ierr)
  rhoq = C0
  rhoq_excl = C0
  if(ierr/=0) stop '[calc_rhoq] Error allocating rhoq,rhoq_excl'
  
  call timing_start('calc rhoq - q-loop')
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate rho(q)
    
  ! get indices of imt and irmd in new mesh -> boundaries for ir loops
  imt1=ipan_intervall(npan_tot)+1
  
  
  
  !threadprivate ararys:
  !$omp parallel private(ierr)
  
  mythread = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  
  
  allocate( qint(irmdnew, lmmaxso, lmmaxso), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating qint'
  
  allocate( eiqr_lm((LMAX_2+1)**2, irmdnew), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating eiqr_lm'
  
  allocate(tmpsum1(N,N), tmpsum2(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmpsum1/2'
  
  allocate(tmp(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating tmp'
  allocate(exG0_tmp(N,N), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating exG0_tmp'
  allocate( jl(0:lmax_2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating jl'
  
  allocate( Ylm((LMAX_2+1)**2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating Ylm'
  allocate( cYlm((LMAX_2+1)**2), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating cYlm'
  
  allocate( q_mu(lmmaxso,lmmaxso), stat=ierr)
  if(ierr/=0) stop '[calc_rhoq] error allocating q_mu'
  
  !$omp end parallel
  
  allocate(q_rand(Nqpt), stat=ierr)
  if (ierr/=0) stop 'Error allocating q_rand'
  if(myrank==master) then
    !initialize q_rand with initial sequence of inttgers from 1 to Nqpt
    do k=1, Nqpt
      q_rand(k) = k
    end do
    ! randomize q_rand by Nrandomize swapping two elements at random
    do k=1, Nrandomize*Nqpt
      call random_number(rand_num)
      ikx = int(rand_num*Nqpt)+1
      call random_number(rand_num)
      iky = int(rand_num*Nqpt)+1
      tmp_int = q_rand(ikx)
      q_rand(ikx) = q_rand(iky)
      q_rand(iky) = tmp_int
    end do
  end if ! myrank==master
  
#ifdef CPP_MPI
  !communicate q_rand to all other ranks
  call MPI_Bcast(q_rand, Nqpt, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  if(ierr/=0) stop 'Error communicating q_rand'


  ! distribute work over ranks
  allocate(ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), stat=ierr)
  if(ierr/=0) stop 'Error allocating ntot_pT etc.'
  
  call distribute_linear_on_tasks(nranks, myrank, master, Nqpt, ntot_pT, ioff_pT, .true.)
  
  q_start = ioff_pT(myrank) + 1
  q_end   = ioff_pT(myrank) + ntot_pT(myrank)
  
  write(*,'(A,I9,A,I9,A,I9)') 'Rank ', myrank, ' does q-points ', q_start, ' to ', q_end
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  
#else

  q_start = 1
  q_end = Nqpt
  
#endif


  ! open big parallel section for k-loop etc. in q-loop, closed after q-loop
  !$omp parallel default(shared) private(q, k, kpq, kweight, i, j, tr_tmpsum1, tr_tmpsum2, nthreads, irec, QdotL, tmpk, tmpr, ifun, lm1, lm2, lm3, ir, l1, m1, lm01, lm02, lm0, Z, phi, costheta, Rq)
  !print header of statusbar
190      FORMAT('                 |'$)   ! status bar
200      FORMAT('|'$)                    ! status bar
  if(myrank==master .and. mythread==master) write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
  if(myrank==master .and. mythread==master) write(*,FMT=190) !beginning of statusbar
  !$omp do schedule(dynamic,1)
  do iq=q_start,q_end !q-loop

     ! take q from randomized array to minimize load imbalance betwee ranks
     q = q_rand(iq)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> k-loop >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! initialize kpt integration
     tmpsum1 = C0
     tmpsum2 = C0

     if(mythread==0) call timing_start('calc rhoq - q>k-loop')

     do k=1,Nkp !k-loop integration
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
       
       if(kmask(kpq) .and. kmask(k)) then
     
         ! read kweight from memory
         kweight = dcmplx(t_rhoq%volcub(k), 0.0d0) ! complex number needed for zgemm later on
         
         ! kweight = kweight * exp(i q*L_i) (exp factor is explicityl included into exG0 in the following)
         
         do i=1,t_rhoq%Nscoef
           do j=1,t_rhoq%Nscoef

             ! Sum( Int( exG0_i(k+q) tau_i,j exG0_j(k); dk ); i,j)
             
             ! collect phase factors
             ! exG0 = G0*exp(-i(k+q)*L_i)
             tmpk(:) = kpt(:,kpq) !Qvec(:,q)+kpt(:,k)
             tmpr(:) = L_i(:,i)
             QdotL = tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3)
             ! exG0 = exG0*exp(+ik*L_j) -> G0tauG0*exp(-[(k+q)*L_i - k*L_j])
             tmpk(:) = kpt(:,k)
             tmpr(:) = L_i(:,j)
             QdotL = QdotL - (tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3))
             ! multiply with phase
             exG0_tmp(1:N,1:N) = G0ji_k(1:N,1:N,j,k)*exp(-2.0d0*pi*Ci*QdotL) !/alat)
             

             ! tmp = tau*exG0_k(k)
             tmp(1:N,1:N) = C0
             call ZGEMM('n','n',N,N,N,C1,tau(1:N,1:N,i,j),N,exG0_tmp(1:N,1:N),N,C0,tmp(1:N,1:N),N)
             ! tmpsum1 = tmpsum1 + G0_k(k+q)*tmp*kweight
             call ZGEMM('n','n',N,N,N,kweight,G0ij_k(1:N,1:N,i,kpq),N,tmp(1:N,1:N),N,C1,tmpsum1(1:N,1:N),N)


             ! Second part:  Int( exG0(k)^*.tau^*.exG0(k+q)^*, dk )
             
             ! collect phase factors
             ! exG0 = G0*exp(-i(k+q)*L_j)
             tmpk(:) = kpt(:,kpq) !Qvec(:,q)+kpt(:,k)
             tmpr(:) = L_i(:,j)
             QdotL = tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3)
             ! exG0 = exG0*exp(+ik*L_i) -> G0tauG0*exp(-[(k+q)*L_j - k*L_i])
             tmpk(:) = kpt(:,k)
             tmpr(:) = L_i(:,i)
             QdotL = QdotL - (tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3))
             ! multiply with phase
             exG0_tmp(1:N,1:N) = dconjg(G0ij_k(1:N,1:N,i,k))*exp(-2.0d0*pi*Ci*QdotL) !/alat)
             
             
             ! tmp = dconjg(tau)*dconjg(exG0_k(k+q))
             tmp(1:N,1:N) = C0
             call ZGEMM('n','n',N,N,N,C1,dconjg(tau(1:N,1:N,i,j)),N,dconjg(G0ji_k(1:N,1:N,j,kpq)),N,C0,tmp(1:N,1:N),N)
             ! tmpsum2 = tmpsum2 + dconjg(exG0_k(k))*tmp*kweight
             call ZGEMM('n','n',N,N,N,kweight,exG0_tmp(1:N,1:N),N,tmp(1:N,1:N),N,C1,tmpsum2(1:N,1:N),N)

                        
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
     Rq = dsqrt(Qvec(1,q)**2+Qvec(2,q)**2+Qvec(3,q)**2)
     if(Rq>=eps) then

        Ylm = 0.0d0
        cYlm = C0
        ! maybe here a factor 2*pi/alat is missing in determination of cos(theta) and phi?
        ! probably this is not the case since Qvec and Rq are in the same units and everything is determined by deviding numbers
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
          Z = dcmplx(Rq*rnew(ir), 0.0d0)!/alat*2.0d0*pi, 0.0d0)
          call CALC_JLK(jl, Z, LMAX_2)
          do l1=0,lmax_2
            do m1=-l1, l1
              lm01 = L1* (L1+1) + M1 + 1
              eiqr_lm(lm01, ir) = jl(l1) * Ylm(lm01) * (Ci**l1)/c0ll/c0ll
            end do ! m1
          end do ! l1
        end do ! ir
        
     else
     
        ! q==0 case: set to zero instead of one since this is anyways cut out
        eiqr_lm(:,:) = C0
        
     end if !(Rq>=eps)
     
     
     if(mythread==0) call timing_pause('calc rhoq - q>eiqr')
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  eiqr  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  qint  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     if(mythread==0) call timing_start('calc rhoq - q>qint')
     ! initialize qint
     qint = C0
     
     ! first treat spherical block (diagonal in lm, lm' -> only lm'==lm)
     do lm02 = 1,lmmaxso
       do lm01 = 1,lmmaxso
           do l1=0,lmax !lmax_2 lmax instead of lmax_2 since trq_of_r has cutoff of 2*lmax and not lmax_2
           ! only diagonal (spherical), thus m=0
           m1 = 0
           lm1 = l1 * (l1+1) + m1 + 1
           ! fill diagonal
           do ir = 1,irmdnew
             qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + trq_of_r(lm01, lm02, lm1, ir) * eiqr_lm(lm1, ir) ! note: diagonal in lm01!!!
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
       if(lmsp(lm3)/=0 .and. lm1<=(2*lmax+1)**2) then ! check if theta/=0 and if trq_of_r element is smaller than curoff of 2*lmax (not lmax_2!)
         do ir=imt1+1,irmdnew
           do lm01=1,lmmaxso
             do lm02=1,lmmaxso
               qint(ir, lm01, lm02) = qint(ir, lm01, lm02) + cleb(j) * trq_of_r(lm01, lm02, lm1, ir) * eiqr_lm(lm2, ir) * thetasnew(ir,ifun)
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
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> C_M(iq) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! calculate trace of Q^mu times k-kpoint integral

     ! C_M(q) = sum_j exp(-i*q*(R_j+Chi_nu)) * q_mu.(G0_ji.tau_ii'.G0_i'j) (all are matrices in lms-space; i,i' are in imp cluster and j in exclude cluster)
     tmp = C0
     do i=1,t_rhoq%Nexcl
        tmpk(:) = Qvec(:,q)
        tmpr(:) = t_rhoq%r_excl(:,i)
        QdotL = tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3)
        kweight = exp(-2.0d0*pi*Ci*QdotL) !/alat)
        call ZGEMM('n','n',N,N,N,kweight,q_mu(1:N,1:N),N,t_rhoq%G0tauG0_excl(i,1:N,1:N),N,C1,tmp(1:N,1:N),N)
        !tmp(1:N,1:N) = tmp(1:N,1:N)+q_mu(1:N,1:N)!t_rhoq%G0tauG0_excl(i,1:N,1:N)
        !tmp(1,1) = tmp(1,1)+kweight
     end do
           
     ! take trace
     tr_tmpsum1 = C0
     do i=1,t_rhoq%lmmaxso
       tr_tmpsum1 = tr_tmpsum1 + tmp(i,i)
     end do

     ! C_M^*(-q) = sum_j exp(i*(-q)*(R_j+Chi_nu)) *  q_mu^*.(G0_ji.tau_ii'.G0_i'j)^* (all are matrices in lms-space; i,i' are in imp cluster and j in exclude cluster)
     tmp = C0
     do i=1,t_rhoq%Nexcl
        tmpk(:) = -Qvec(:,q)
        tmpr(:) = t_rhoq%r_excl(:,i)
        QdotL = tmpr(1)*tmpk(1)+tmpr(2)*tmpk(2)+tmpr(3)*tmpk(3)
        kweight = exp(2.0d0*pi*Ci*QdotL) !/alat)
        call ZGEMM('n','n',N,N,N,kweight,dconjg(q_mu(1:N,1:N)),N,t_rhoq%G0tauG0_excl(i+t_rhoq%Nexcl,1:N,1:N),N,C1,tmp(1:N,1:N),N)
        !tmp(1:N,1:N) = tmp(1:N,1:N)+dconjg(q_mu(1:N,1:N))!t_rhoq%G0tauG0_excl(i+t_rhoq%Nexcl,1:N,1:N)
        !tmp(1,1) = tmp(1,1)+kweight
     end do
           
     ! take trace
     tr_tmpsum2 = C0
     do i=1,t_rhoq%lmmaxso
       tr_tmpsum2 = tr_tmpsum2 + tmp(i,i)
     end do
     
     rhoq_excl(q) = -1.d0/(2.d0*Ci*pi) * (tr_tmpsum1 - tr_tmpsum2)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< C_M(iq) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> rhoq(iq) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! calculate trace of Q^mu times k-kpoint integral

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
     tr_tmpsum1 = C0
     tr_tmpsum2 = C0
     do i=1,t_rhoq%lmmaxso
       tr_tmpsum1 = tr_tmpsum1 + tmpsum1(i,i)
       tr_tmpsum2 = tr_tmpsum2 + tmpsum2(i,i)
     end do
     
     rhoq(q) = -1.d0/(2.d0*Ci*pi) * (tr_tmpsum1 - tr_tmpsum2)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< rhoq(iq) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !update statusbar
     if(myrank==master) then
       !$omp critical
       if( ((q_end-q_start)>=50.and.mod(q,(q_end-q_start)/50)==0) .or. ((q_end-q_start)<50) ) write(*,FMT=200)
       !$omp end critical
     end if
       
  end do !q
  !$omp end do
  
  !finish big parallel section
  !$omp end parallel
  
  if(mythread==0) call timing_stop('calc rhoq - q>k-loop')
  if(mythread==0) call timing_stop('calc rhoq - q>eiqr')
  if(mythread==0) call timing_stop('calc rhoq - q>qint')
  
  
  ! threadprivate arrays
  !$omp parallel private(ierr)
  deallocate(qint, eiqr_lm, Ylm, cYlm, tmpsum1, tmpsum2, tmp, jl, exG0_tmp, q_mu, stat=ierr)
  if(ierr/=0) stop 'Error deallocating threadprivate arrays'
  !$omp end parallel
  
  
  call timing_stop('calc rhoq - q-loop')
  
  ! done calculating rho(q)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(myrank==master) then
    write(*,*) '' ! finish statusbar
    write(*,*) 'saving rhoq(q) to out_rhoq.txt'
  end if
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write out result
  call timing_start('calc rhoq - writeout')
  
#ifdef CPP_MPI
  if(myrank==master) write(*,*) 'MPI communication ...'
  !Gather results on master
  allocate(tmp(Nqpt,1), stat=ierr)
  if(ierr/=0) stop 'Error allocating tmp for REDUCE of rhoq'

  call MPI_REDUCE(rhoq, tmp, Nqpt, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop 'ERROR in REDUCE for rhoq(q) before writeout'
  if(myrank==master) then
    rhoq(:) = tmp(:,1)
  endif
    
  call MPI_REDUCE(rhoq_excl, tmp, Nqpt, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
  if(ierr/=MPI_SUCCESS) stop 'ERROR in REDUCE for rhoq(q) before writeout'
  if(myrank==master) then
    rhoq_excl(:) = tmp(:,1)
  endif
    
  deallocate(tmp, stat=ierr)
  if(ierr/=0) stop 'Error deallocating tmp for REDUCE of rhoq'
#endif
  
  if(myrank==master) then
    write(*,'("writeout loop:   |",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
    write(*,FMT=190) !beginning of statusbar
    open(9999, file='out_rhoq.txt', form='formatted')
    do q=1,Nqpt
       write(9999,'(7E15.7)') Qvec(:,q),rhoq(q),rhoq_excl(q)
       if( (Nqpt>=50.and.mod(q,Nqpt/50)==0) .or. (Nqpt<50) ) write(*,FMT=200)
    end do
    close(9999)
    write(*,*) '' ! finish statusbar
  end if
  
  call timing_stop('calc rhoq - writeout')
  ! done with write out result
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
#ifdef CPP_MPI
  deallocate(ntot_pT, ioff_pT, stat=ierr)
  if(ierr/=0) stop 'Error deallocating ntot_pT etc.'
#endif
  
  
  deallocate(rhoq, rhoq_excl)
    
end subroutine calc_rhoq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module mod_rhoq


! standalone version for testing only, comment out when ready
!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!
program test
  
  
#ifdef CPP_MPI
  use mpi
#endif
  use mod_mympi, only: mympi_init, myrank, nranks, master
  use mod_types, only: t_inc
  use mod_rhoq
  use mod_timing
  use mod_version_info
  use mod_mympi, only: myrank, master, nranks

  implicit none
  
  character*256 :: uio ! unit in which to find inputcard
  integer :: nkpt, Nkx, Nky ! total number of kpoints, k-mesh in x/y
  integer :: lmmaxso, natyp
  integer :: naez,ncls,nr,nemb
  integer :: lmax, ntotd, npan_tot, ncheb
  integer :: nsra, irmdnew ! left/right solutions (wave functions)
  integer :: ncleb_max ! numer of lm-indices in cleb array (dimension of trq_of_r)
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
  double complex, allocatable   :: trq_of_r(:,:,:,:)
  double precision, allocatable :: recbv(:,:), bravais(:,:) ! lattice information > Bravais matrix in real and reciprocal space
  double precision, allocatable :: rnew(:)  ! r mesh
  
  double precision, parameter :: pi = 4.0d0*datan(1.0d0) ! pi/4 = arctan(1)

  
  integer :: lm1, lm2, i, j, ir, ierr
  
#ifdef CPP_MPI
  ! initialize MPI
  call MPI_Init ( ierr )
#endif
  ! set variables master, myrank and nranks for serial (myrank=master=0, nranks=1) as well as parallel execution
  call mympi_init()
  
  
  ! set timing output
  call construct_serialnr()
  t_inc%i_time = 1
  call timing_init(myrank)
  call timing_start('total time')
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>>>>>>>>>>>>>>>> read stuff in  >>>>>>>>>>>>>
    call timing_start('read stuff')
  if(myrank==master) then
    
    ! read in scalars
    uio = 'inputcard'
    open(9999, file='params.txt')
    read(9999,*) lmmaxso, natyp
    read(9999,*) naez, ncls, nr, nemb, lmax
    read(9999,*) alat
    close(9999)

    

    open(9999, file='kpts.txt', form='formatted')
    read(9999,'(3I9)') nkpt, Nkx, Nky
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
!     recbv = recbv * (2.0d0*pi)/alat  !* alat !2*pi! * alat
!     kpt = kpt * (2.0d0*pi)/alat      !* alat !2*pi! * alat
    
    open(9999, file='host.txt')
    allocate( r_basis(3,natyp), Ghost(lmmaxso,lmmaxso,t_rhoq%Nlayer) )
    allocate( rcls(3,ncls,ncls), rr(3,0:nr) )
    allocate( atom(ncls,naez+nemb), cls(naez+nemb), ezoa(ncls,naez+nemb), nacls(ncls) )
    read(9999,*) r_basis(1:3,1:natyp)
    read(9999,*) rcls(1:3,1:ncls,1:ncls), rr(1:3,0:nr), atom(1:ncls,1:naez+nemb)
    read(9999,*) cls(1:naez+nemb), ezoa(1:ncls,1:naez+nemb), nacls(1:ncls)
    close(9999)
    
    open(9999, file='wavefunctions.txt')
    read(9999,'(100I9)') ntotd, npan_tot, ncheb, nsra, irmdnew
    write(*,*) ntotd, npan_tot, ncheb, nsra, irmdnew, lmmaxso
    allocate( rnew(irmdnew) )
           read(9999,'(1000E26.17)') rnew(1:irmdnew)
    
    allocate( Rll(1:nsra*lmmaxso,1:lmmaxso,1:irmdnew),          &
   &          Rllleft(1:nsra*lmmaxso,1:lmmaxso,1:irmdnew) ,     &
   &          rpan_intervall(0:ntotd),ipan_intervall(0:ntotd) )
    do ir=1,irmdnew
      do lm1=1,nsra*lmmaxso
        do lm2=1,lmmaxso
          read(9999,'(20000E16.7)') Rll(lm1, lm2, ir), Rllleft(lm1, lm2, ir)
        end do
      end do
    enddo
    do lm1=0,npan_tot
      read(9999,'(E16.7,I9)') rpan_intervall(lm1), ipan_intervall(lm1)
    enddo
    close(9999)
  
  endif !myrank==master
  call timing_stop('read stuff')
  !<<<<<<<<<<<<<<<< read stuff in  <<<<<<<<<<<<<
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>>>>>>>>>>>>>>>>> init rhoq >>>>>>>>>>>>>>>>>
  call timing_start('init rhoq')
  
  if(myrank==master) then

    ! find ilay_scoef, r_scoef, Nlayer, Nscoef in t_rhoq from scoef
    call read_scoef_rhoq(t_rhoq)

    ! read from inputcard
    call read_input_rhoq(t_rhoq, uio)

    ! save kpt, volcub, volbz nkpt to t_rhoq
    call save_kmesh_rhoq(t_rhoq,nkpt,kpt,volcub,volbz, Nkx, Nky)

    ! save r_basis, lmmaxso, natyp to t_rhoq
    call save_geometry_rhoq(t_rhoq,r_basis,lmmaxso,natyp)

  endif !myrank==master
  
#ifdef CPP_MPI
  ! broadcast scalars in t_rhoq from master to other ranks
  call bcast_scalars_rhoq(t_rhoq)
#endif

  !allocate arrays in t_rhoq from scalars (used for myrank/=master)
  call init_t_rhoq(t_rhoq)
  
#ifdef CPP_MPI
  ! bradcast arrays of t_rhoq from master to all other ranks 
  call bcast_arrays_rhoq(t_rhoq)
#endif

  call timing_stop('init rhoq')
  !<<<<<<<<<<<<<<<<< init rhoq <<<<<<<<<<<<<<<<<
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  ! read DTMTRX and GMATLL_GES (scattering files), stored in t_rhoq%
  call timing_start('read Dt Gimp')
  call read_Dt_Gimp_rhoq(t_rhoq, t_rhoq%lmmaxso, t_rhoq%Nscoef)
  call timing_stop('read Dt Gimp')
  
  ! calculate impurity scattering path operator from Dt, Gimp
  ! tau_i,i' = Dt_i delta_i,i' + Dt_i Gimp_i,i' Dt_i'
  call timing_start('calc tau')
  call calc_tau_rhoq(t_rhoq)
  call timing_stop('calc tau')

  ! read exclude cluster and green_host_excl, then precaclulate sum_jk G0_ij(E).tau_jk.G0_ki(E)
  ! this is needed for exclude cluster of C_M calculation
  call timing_start('read in of exclude cluster stuff')
  call start_excl(t_rhoq)
  call timing_stop('read in of exclude cluster stuff')
  
  
  ! calculate prefactor Q^{\mu}_{LL'}(r) = R_{L}(\vec{r})*Rleft_{L'}(\vec{r})
  ! used in calc_rhoq to get Q^{\mu}_{LL'}(\vec{q}) = \Int R_{L}(\vec{r})*Rleft_{L'}(\vec{r}) \exp{-i\vec{q}\cdot\vec{r}} d\vec{r}
  call timing_start('calc Q_mu')
  if(myrank==master) then
    write(*,*) 'before calc_Q_mu', irmdnew
    call calc_Q_mu_rhoq(lmax, ntotd, npan_tot, &
          & nsra, lmmaxso, Rll, Rllleft, ipan_intervall, &
          & irmdnew, trq_of_r, ncleb_max)
    write(*,*) 'after calc_Q_mu', irmdnew, ncleb_max
  endif
#ifdef CPP_MPI
  call MPI_Bcast(ncleb_max, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(irmdnew, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  
  if(myrank/=master) then
    allocate(recbv(3,3), stat=ierr)
    if(ierr/=0) stop 'Error allocating recbv for myrank/=master'
  end if
  call MPI_Bcast(recbv, 3*3, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  
  call MPI_Bcast(lmax, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ntotd, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(npan_tot, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ncheb, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(alat, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  
  if(myrank/=master) then
    allocate(rpan_intervall(0:ntotd),ipan_intervall(0:ntotd), stat=ierr)
    if(ierr/=0) stop 'Error allocating rpan_ntervall for myrank/=master'
  end if
  call MPI_Bcast(rpan_intervall, ntotd+1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ipan_intervall, ntotd+1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
  
  
  if(myrank/=master) then
    allocate(trq_of_r(t_rhoq%lmmaxso,t_rhoq%lmmaxso,ncleb_max,irmdnew), stat=ierr)
    !allocate(trq_of_r(t_rhoq%lmmaxso,t_rhoq%lmmaxso,100,irmdnew), stat=ierr)
    if(ierr/=0) stop 'Error allocating trq_of_r for myrank/=master'
    trq_of_r(:,:,:,:) = (0.0d0, 0.0d0)
  end if
  write(*,*) myrank, shape(trq_of_r), t_rhoq%lmmaxso, irmdnew, ierr, master
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  write(*,*) 'going in bcast', myrank
  call MPI_Bcast(trq_of_r(:,:,:,:), t_rhoq%lmmaxso*t_rhoq%lmmaxso*ncleb_max*irmdnew, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  !all MPI_Allreduce(trq_of_r, t_rhoq%lmmaxso*t_rhoq%lmmaxso*100*irmdnew, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
  write(*,*) 'out of bcast', myrank
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  call timing_stop('calc Q_mu')


  ! calculate Fourier transform: \rho(\vec{q}) = \int \Delta\rho(\vec{q};\Chi_\mu+\vec{r}) d\vec{r}
  ! = Tr{ Q^{\mu}_{LL'}(q) \int [ \sum_{i,j}{\exp{-i\vec{q}\cdot (\vec{L}_i}-\vec{L}_j)} G0^{\mu_0,i}(k) \tau_i,j G0^{j,\mu_0}(k+q)} - \sum_{i,j}{\exp{-i\vec{q}\cdot (\vec{L}_i}-\vec{L}_j)} G0^{\mu_0,i}(k) \tau_i,j G0^{j,\mu_0}(k-q)}^* ] d^2k
  ! where G0^{\mu,\mu_i} (k) = -t_{\mu}^-1 \delta_{\mu-\mu_i} + t_{mu}^-1 \tau_0^{\mu,\mu_i}(k) t_{mu_i}^-1
  call timing_start('calc rhoq')
  call calc_rhoq(t_rhoq, t_rhoq%lmmaxso, t_rhoq%Nkpt, trq_of_r, recbv, lmax,   &
        &        ntotd, npan_tot, ncheb, rpan_intervall, ipan_intervall, irmdnew, alat, ncleb_max)
  call timing_stop('calc rhoq')
  
  
  
  call timing_stop('total time')
  
  
#ifdef CPP_MPI
  ! finalize MPI
  call MPI_Finalize ( ierr )
#endif


end program test
!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!TEST!
