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

  integer :: i1_start, i1_end, ierr,i_stat,i_all

  do I1=i1_start,i1_end

    IPOT=NSPIN*(I1-1)+1

    if (TEST('BdG_dev ')) then
      ! write out inputs for tmat_newsolver to extract first BdG 
      if (nranks>1) stop 'test option BdG_dev can only be used in serial!'
      if (i1==i1_start) open(887766, file='BdG_tmat_inputs.txt', form='formatted')
      if (i1==1) then 
        write(887766, '(A)') 'global parameters:'
        write(887766, *) 
        write(887766, '(A,I9)') 'IELAST= ', IELAST
        write(887766, '(A,I9)') 'NSPIN= ', NSPIN
        write(887766, '(A,I9)') 'LMAX= ', LMAX
        write(887766, '(A,I9)') 'NSRA= ', NSRA
        write(887766, '(A,I9)') 'IEND= ', IEND
        write(887766, '(A,I9)') 'LMPOTD= ', LMPOTD
        write(887766, '(A,I9)') 'LLY= ', LLY
        write(887766, '(A,2ES21.9)') 'DELTAE= ', DELTAE
        write(887766, '(A,I9)') 'IDOLDAU= ', IDOLDAU
        write(887766, '(A,I9)') 'NCLEB= ', NCLEB
        write(887766, '(A,I9)') 'NCHEB= ', NCHEB
        write(887766, '(A,I9)') 'NTOTD= ', NTOTD
        write(887766, '(A,I9)') 'MMAXD= ', MMAXD
        write(887766, '(A,I9)') 'NSPIND= ', NSPIND
        write(887766, '(A,I9)') 'IEMXD= ', IEMXD
        write(887766, '(A,I9)') 'NRMAXD= ', NRMAXD
        write(887766, '(A,I9)') 'NSPOTD= ', NSPOTD
        write(887766, '(A,999999999ES21.9)') 'CLEB= ', CLEB(:,1)
        write(887766, '(A,999999999I9)') 'ICLEB= ', ICLEB
        write(887766, '(A,999999999ES21.9)') 'EZ= ', EZ
        write(887766, *) 
        write(887766, '(A)') 'atom-dependent input:'
        write(887766, *) 
      end if 
      write(887766, '(A,I9)') 'I1= ', I1
      write(887766, '(A,I9)') 'IPOT= ', IPOT
      write(887766, '(A,99999999I9)') 'IPOT= ', IPOT
      write(887766, '(A,I9)') 'NPAN_TOT= ', NPAN_TOT(I1)
      write(887766, '(A,I9)') 'LPOT= ', LOPT(I1)
      write(887766, '(A,999999999I9)') 'IPAN_INTERVALL= ', IPAN_INTERVALL(:,I1)
      write(887766, '(A,ES21.9)') 'ZAT= ', ZAT(I1)
      write(887766, '(A,ES21.9)') 'PHI= ', PHI(I1)
      write(887766, '(A,ES21.9)') 'THETA= ', THETA(I1)
      write(887766, '(A,ES21.9)') 'SOCSCALE= ', SOCSCALE(I1)
      write(887766, '(A,999999999ES21.9)') 'RNEW= ', RNEW(:,I1)
      write(887766, '(A,999999999ES21.9)') 'RPAN_INTERVALL= ', RPAN_INTERVALL(:,I1)
      write(887766, '(A,999999999ES21.9)') 'WLDAU= ', WLDAU(:,:,:,I1)
      write(887766, '(A,999999999ES21.9)') 'VINSNEW= ', VINSNEW
      write(887766, '(A,I9)') 'end atom=', I1
      write(887766, *)
      if (i1==i1_end) then
        close(887766)
        stop 'done writing tmat_newsolver input of test option BdG_dev'
      end if
    end if

    call TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT(I1),SOCSCALE(I1),EZ,  &
       NSRA,CLEB(:,1),ICLEB,IEND,NCHEB,NPAN_TOT(I1),                &
       RPAN_INTERVALL(:,I1),IPAN_INTERVALL(:,I1),RNEW(:,I1),        &
       VINSNEW,THETA(I1),PHI(I1),I1,IPOT,LMPOTD,LLY,DELTAE,IDOLDAU, &
       LOPT(I1),WLDAU(:,:,:,I1),t_dtmatJij(I1))

  enddo !I1, atom loop

end subroutine main1a_dummy
