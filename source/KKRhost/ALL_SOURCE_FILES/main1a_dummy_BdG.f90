!> this is dummy version of main1a
subroutine main1a_dummy
  use Profiling
  use Constants
  use global_variables
  Use mod_datatypes, Only: dp

  use mod_types, only: t_tgmat, t_inc, t_lloyd, t_dtmatJij,init_t_dtmatJij, init_t_dtmatJij_at
  use mod_mympi, only: nranks, master, myrank
  use mod_timing
  use mod_wunfiles
  use mod_jijhelp, only: set_Jijcalc_flags

  use mod_main0
  use mod_tmatnewsolver, only: tmat_newsolver


  implicit none

  ! .. Local variables
  integer :: I1
  integer :: IPOT
  integer :: ILTMP
  integer :: ISPIN
  integer :: ITMPDIR
  character(len=80) :: TMPDIR
  logical :: OPT
  logical :: TEST
  logical :: LREFSYS
  integer :: LRECTMT
  integer :: LRECTRA
  ! .. Local arrays
  real (kind=dp), dimension(NATYPD) :: PHI
  real (kind=dp), dimension(NATYPD) :: THETA
  real (kind=dp), dimension(:,:,:), allocatable :: VINSNEW

  integer :: i1_run, i1_start, i1_end, ierr,i_stat,i_all

  ! for data import:
  character(len=25) dummy
  integer :: ier
  character (len=256) :: uio                             ! NCOLIO=256

  ! BdG specific:
  logical :: use_BdG
  complex (kind=dp) :: delta_BdG


  !================
  ! start read-in

  write(*,*) 'start reading BdG-inputs from inputcard ...'

  call ioinput('use_BdG         ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) use_BdG
    write(*,*) 'use_BdG= ', use_BdG
  else
    stop '[main1a_dummy] error "use_BdG" not found in inputcard'
  endif

  call ioinput('delta_BdG       ',uio,1,7,ier)
  if (ier.eq.0) then
    read (unit=uio,fmt=*) delta_BdG
    write(*,*) 'delta_BdG= ', delta_BdG
  else
    stop '[main1a_dummy] error "delta_BdG" not found in inputcard'
  endif

  ! end read-in
  !================


  write(*,*) 'now read atom-specific input for tmat_newsolver'

  allocate(VINSNEW(NRMAXD,LMPOTD,NSPOTD),stat=i_stat)
  call memocc(i_stat,product(shape(VINSNEW))*kind(VINSNEW),'VINSNEW','main1a')
  VINSNEW=0.0D0

  i1_start = 1
  i1_end = NATYP

  do I1_run=i1_start,i1_end

    if (TEST('BdG_dev ')) then
      ! read out inputs for tmat_newsolver to extract first BdG 
      if (nranks>1) stop 'test option BdG_dev can only be used in serial!'
      if (i1_run==1) open(887766, file='BdG_tmat_inputs.txt', form='formatted')
      if (i1_run==1) then 
        read(887766, *) dummy
        read(887766, *)  
        read(887766, *) dummy, IELAST
        read(887766, *) dummy, NSPIN
        read(887766, *) dummy, LMAX
        read(887766, *) dummy, NSRA
        read(887766, *) dummy, IEND
        read(887766, *) dummy, LMPOTD
        read(887766, *) dummy, LLY
        read(887766, '(A,2ES21.9)') dummy, DELTAE
        read(887766, *) dummy, IDOLDAU
        read(887766, *) dummy, NCLEB
        read(887766, *) dummy, NCHEB
        read(887766, *) dummy, NTOTD
        read(887766, *) dummy, MMAXD
        read(887766, *) dummy, NSPIND
        read(887766, *) dummy, IEMXD
        read(887766, *) dummy, NRMAXD
        read(887766, *) dummy, NSPOTD
        read(887766, *) dummy
        read(887766, *) CLEB(:,1)
        read(887766, *) dummy
        read(887766, *) ICLEB(:,:)
        read(887766, *) dummy
        read(887766, '(2ES21.9)') EZ
        read(887766, *)
        read(887766, *) dummy
        read(887766, *)
      end if 
      read(887766, *) dummy, I1
      read(887766, *) dummy, IPOT
      read(887766, *) dummy, NPAN_TOT(I1)
      read(887766, *) dummy, LOPT(I1)
      read(887766, *) dummy, IPAN_INTERVALL(:,I1)
      read(887766, *) dummy, ZAT(I1)
      read(887766, *) dummy, PHI(I1)
      read(887766, *) dummy, THETA(I1)
      read(887766, *) dummy, SOCSCALE(I1)
      read(887766, *) dummy, RNEW(:,I1)
      read(887766, *) dummy, RPAN_INTERVALL(:,I1)
      read(887766, *) dummy, WLDAU(:,:,:,I1)
      read(887766, *) dummy, VINSNEW
      read(887766, *) 
      if (i1==i1_end) then
        close(887766)
        write(*,*) 'done reading tmat_newsolver input of test option BdG_dev'
      end if
    end if

    write(*,*) 'start tmat_newsolver ...'

    call init_t_dtmatJij(t_inc,t_dtmatJij)

    call TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT(I1),SOCSCALE(I1),EZ,  &
       NSRA,CLEB(:,1),ICLEB,IEND,NCHEB,NPAN_TOT(I1),                &
       RPAN_INTERVALL(:,I1),IPAN_INTERVALL(:,I1),RNEW(:,I1),        &
       VINSNEW,THETA(I1),PHI(I1),I1,IPOT,LMPOTD,LLY,DELTAE,IDOLDAU, &
       LOPT(I1),WLDAU(:,:,:,I1),t_dtmatJij(I1))

  enddo !I1, atom loop

end subroutine main1a_dummy
