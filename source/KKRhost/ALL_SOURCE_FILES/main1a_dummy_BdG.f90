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

  character(len=25) dummy

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
        read(887766, *) dummy, DELTAE
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
        write(*,*), 'EZ', EZ
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

      if (i1==i1_start) open(887766, file='check_BdG_tmat_inputs.txt', form='formatted')
      if (i1==1) then 
        write(887766, '(A25)') 'global parameters:'
        write(887766, *) 
        write(887766, '(A25,I9)') 'IELAST= ', IELAST
        write(887766, '(A25,I9)') 'NSPIN= ', NSPIN
        write(887766, '(A25,I9)') 'LMAX= ', LMAX
        write(887766, '(A25,I9)') 'NSRA= ', NSRA
        write(887766, '(A25,I9)') 'IEND= ', IEND
        write(887766, '(A25,I9)') 'LMPOTD= ', LMPOTD
        write(887766, '(A25,I9)') 'LLY= ', LLY
        write(887766, '(A25,2ES21.9)') 'DELTAE= ', DELTAE
        write(887766, '(A25,I9)') 'IDOLDAU= ', IDOLDAU
        write(887766, '(A25,I9)') 'NCLEB= ', NCLEB
        write(887766, '(A25,I9)') 'NCHEB= ', NCHEB
        write(887766, '(A25,I9)') 'NTOTD= ', NTOTD
        write(887766, '(A25,I9)') 'MMAXD= ', MMAXD
        write(887766, '(A25,I9)') 'NSPIND= ', NSPIND
        write(887766, '(A25,I9)') 'IEMXD= ', IEMXD
        write(887766, '(A25,I9)') 'NRMAXD= ', NRMAXD
        write(887766, '(A25,I9)') 'NSPOTD= ', NSPOTD
        write(887766, '(A25)') 'CLEB= '
        write(887766, '(999999999ES21.9)') CLEB(:,1)
        write(887766, '(A25)') 'ICLEB= '
        write(887766, '(999999999I9)') ICLEB(:,:)
        write(887766, '(A25)') 'EZ= '
        write(887766, '(2ES21.9)') EZ
        write(887766, *) 
        write(887766, '(A25)') 'atom-dependent input:'
        write(887766, *) 
      end if 
      write(887766, '(A25,I9)') 'I1= ', I1
      write(887766, '(A25,I9)') 'IPOT= ', IPOT
      write(887766, '(A25,I9)') 'NPAN_TOT= ', NPAN_TOT(I1)
      write(887766, '(A25,I9)') 'LPOT= ', LOPT(I1)
      write(887766, '(A25,999999999I9)') 'IPAN_INTERVALL= ', IPAN_INTERVALL(:,I1)
      write(887766, '(A25,ES21.9)') 'ZAT= ', ZAT(I1)
      write(887766, '(A25,ES21.9)') 'PHI= ', PHI(I1)
      write(887766, '(A25,ES21.9)') 'THETA= ', THETA(I1)
      write(887766, '(A25,ES21.9)') 'SOCSCALE= ', SOCSCALE(I1)
      write(887766, '(A25,999999999ES21.9)') 'RNEW= ', RNEW(:,I1)
      write(887766, '(A25,999999999ES21.9)') 'RPAN_INTERVALL= ', RPAN_INTERVALL(:,I1)
      write(887766, '(A25,999999999ES21.9)') 'WLDAU= ', WLDAU(:,:,:,I1)
      write(887766, '(A25,999999999ES21.9)') 'VINSNEW= ', VINSNEW
      write(887766, *)
      if (i1==i1_end) close(887766)
    end if

         call init_t_dtmatJij(t_inc,t_dtmatJij)
         if(OPT('XCPL    '))then
            call set_Jijcalc_flags(t_dtmatJij,NATYP,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
         end if !OPT('XCPL')

    call TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT(I1),SOCSCALE(I1),EZ,  &
       NSRA,CLEB(:,1),ICLEB,IEND,NCHEB,NPAN_TOT(I1),                &
       RPAN_INTERVALL(:,I1),IPAN_INTERVALL(:,I1),RNEW(:,I1),        &
       VINSNEW,THETA(I1),PHI(I1),I1,IPOT,LMPOTD,LLY,DELTAE,IDOLDAU, &
       LOPT(I1),WLDAU(:,:,:,I1),t_dtmatJij(I1))

  enddo !I1, atom loop

end subroutine main1a_dummy
