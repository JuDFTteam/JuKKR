module mod_fermisurf_basic

  implicit none

  private
  public :: ROOT_IMAG, ROOT_ANY, ROOT_REAL, roots_along_edge, compare_two_eigv_in_substeps, &
          & connect_eigw_in_substeps, find_roots_any_eigw, mark_cubes_FScross, &
          & generate_cubevertices, generate_squarevertices, unroll_ixyz, &
          & read_cubesfile, save_cubesfile, find_kpoints_irredset, save_kpointsfile_vis, &
          & save_kpointsfile_int, testpath, get_cubesinfo_filename

  integer, parameter :: ROOT_ANY=0, ROOT_REAL=1, ROOT_IMAG=2

contains





  subroutine testpath(inc, lattice, cluster, tgmatrx)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_eigvects,   only: normeigv_new, compare_eigv
!   use mod_symmetries, only: symmetries_type, set_symmetries, get_IBZwedge_faces
    use mod_kkrmat,     only: compute_kkr_eigenvectors
    use mod_mathtools,  only: bubblesort
    use mod_ioinput,    only: IoInput
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(inout) :: tgmatrx

    integer          :: ltest, nsteps
    double precision :: k1(3), k2(3), dkp(3)
    double precision, allocatable :: ksub(:,:)
    double complex,   allocatable :: LVeig(:,:,:),&
                                   & RVeig(:,:,:),&
                                   & eigw(:,:)

!   type(symmetries_type) :: symmetries
    integer               :: nfaces
    double precision      :: bounds(3,2)
    double precision, allocatable :: nvec(:,:), dscal(:)

    integer          :: sorted(inc%almso)
    double precision :: eigwabs(inc%almso)


    integer          :: nrootsteps, nrootiter, roottype, nroot, lmroot(inc%nrootmax)
    double precision :: kends(3,2), rooteps
    double precision :: kroot(3,inc%nrootmax)
!   double complex   :: LVroot(inc%almso,inc%almso,inc%nrootmax),&
!                     & RVroot(inc%almso,inc%almso,inc%nrootmax),&
!                     & eigwroot(inc%almso,inc%nrootmax)
    double complex, allocatable :: LVroot(:,:,:),&
                                 & RVroot(:,:,:),&
                                 & eigwroot(:,:)
    character(len=8) :: rootinp


    integer :: itest, nlmtest, ikp, lm0, lm1, lm2, ierr
    character(len=80) :: filename
    character(len=80) :: uio
    integer,          allocatable :: ilmtest(:), connection(:,:)

    integer :: indum(6)
    integer :: nCub3(3), nmarked
    integer, allocatable :: imarked(:)

    !--------------------------------------!
    !--- Read in options from inputcard ---!
    !--------------------------------------!
    call IoInput('LTEST     ',uio,1,7,ierr)
    read(unit=uio,fmt=*) ltest
    if(ltest/=1) return

    !initialize symmetries
!   call set_symmetries(inc, lattice, symmetries)
!   call get_IBZwedge_faces(lattice%recbv, symmetries%nsym_used, symmetries%rotmat, symmetries%isym_used, nfaces, nvec, dscal, bounds)
!   call init_cube2tetralines()

    ! TESTPART 1 - PRINT EIGENVALUES
    call IoInput('TESTKPOI1 ',uio,1,7,ierr)
    read(unit=uio,fmt=*) k1
    call IoInput('TESTKPOI2 ',uio,1,7,ierr)
    read(unit=uio,fmt=*) k2
    call IoInput('TESTNKPTS ',uio,1,7,ierr)
    read(unit=uio,fmt=*) nsteps
    call IoInput('NLMTEST   ',uio,1,7,ierr)
    read(unit=uio,fmt=*) nlmtest
    allocate( ilmtest(nlmtest) )
    call IoInput('ILMTEST   ',uio,1,7,ierr)
    read(unit=uio,fmt=*) ilmtest

    allocate( ksub(3,nsteps+1),&
            & LVeig(inc%almso,inc%almso,nsteps+1),&
            & RVeig(inc%almso,inc%almso,nsteps+1),&
            & eigw(inc%almso,nsteps+1),&
            & connection(inc%almso,nsteps+1),&
            & STAT=ierr )
    if(ierr/=0) stop 'Problem allocating arrays in testpath'

    allocate( LVroot(inc%almso,inc%almso,inc%nrootmax),&
            & RVroot(inc%almso,inc%almso,inc%nrootmax),&
            & eigwroot(inc%almso,inc%nrootmax),        &
            & STAT=ierr                                )
    if(ierr/=0) stop 'Problem allocating LVroot etc. in testpath'

    !----------------------------------------------------!
    !--- create submesh to divide the input-intervall ---!
    !---  and find the KKR-matrix, -eigenvectors and  ---!
    !---      -eigenvalues on each sub-kpoint         ---!
    !----------------------------------------------------!
    dkp = (k2 - k1)/nsteps
    do ikp=1,nsteps+1
      ksub(:,ikp) = k1 + dble(ikp-1)*dkp
      call compute_kkr_eigenvectors( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                   & tgmatrx%tinvll(:,:,:,1), ksub(:,ikp),        &
                                   & eigw(:,ikp), LVeig(:,:,ikp), RVeig(:,:,ikp)  )
!     call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
!                                   & tgmatrx%tmat(:,:,1), ksub(:,ikp),            &
!                                   & eigw(:,ikp), LVeig(:,:,ikp), RVeig(:,:,ikp)  )
    end do!ikp


    !write the (sorted) eigenvalues to a file
    open(unit=135,file='test.eigw.txt',form='formatted',action='write')
    write(135,'("# real(eigw)             | imag(eigw)             | abs(eigw)              | ikp | lm")')
    do ikp=1,nsteps+1
      eigwabs = abs(eigw(:,ikp))
      call bubblesort(inc%almso, eigwabs, sorted)
      do lm0=1,inc%almso
        write(135,'(3ES25.16,2(2X,I0))') eigw(sorted(lm0),ikp), abs(eigw(sorted(lm0),ikp)), ikp, sorted(lm0)
      end do!lm0=1,inc%almso
    end do!ikp
    close(135)

    !for selected eigenvalues, write the eigenvalue-path to a file
    do itest=1,nlmtest
      lm1 = ilmtest(itest)

      write(filename,'("test.path_lm=",I0,".txt")') lm1
      open(unit=368, file=trim(filename), action='write', form='formatted')
      write(368,'("# real(eigw)             | imag(eigw)             | abs(eigw)              | ikp | lm")')
      write(368,'(3ES25.16,2X,I4,2X,I3)') eigw(lm1,1), abs(eigw(lm1,1)), 1, lm1

      write(filename,'("test.proj_lm=",I0,".txt")') lm1
      open(unit=468, file=trim(filename), action='write', form='formatted')
      write(468,'("#  max_{-n}(proj) .... max(proj) ")')

      do ikp=1,nsteps
        call compare_eigv(inc, LVeig(:,lm1,ikp), RVeig(:,:,ikp+1), lm2, 468)
        write(368,'(3ES25.16,2X,I4,2X,I3)') eigw(lm2,ikp+1), abs(eigw(lm2,ikp+1)), ikp+1, lm2
        lm1 = lm2
      end do!ikp
      close(368)
      close(468)
    end do!itest



!######################################################
!### test the subroutine 'connect_eigw_in_substeps' ###
!######################################################
    kends(:,1) = k1
    kends(:,2) = k2
    call connect_eigw_in_substeps( inc, lattice, cluster, tgmatrx, nsteps, kends, &
                                 & connection, ksub, eigw, LVeig, RVeig           )

    !write the (sorted) eigenvalues to a file
    open(unit=335,file='test.connect.eigw.txt',form='formatted',action='write')
    write(335,'("# real(eigw)             | imag(eigw)             | abs(eigw)              | ikp | lm")')
    do ikp=1,nsteps+1
      eigwabs = abs(eigw(:,ikp))
      call bubblesort(inc%almso, eigwabs, sorted)
      do lm0=1,inc%almso
        write(335,'(3ES25.16,2(2X,I0))') eigw(sorted(lm0),ikp), abs(eigw(sorted(lm0),ikp)), ikp, sorted(lm0)
      end do!lm0=1,inc%almso
    end do!ikp
    close(335)

    !for selected eigenvalues, write the eigenvalue-path to a file
    do itest=1,nlmtest
      lm1 = ilmtest(itest)
      write(filename,'("test.connect.lm=",I0,".txt")') lm1
      open(unit=368, file=trim(filename), action='write', form='formatted')

      write(368,'("# real(eigw)             | imag(eigw)             | abs(eigw)              | ikp | lm")')
      write(368,'(3ES25.16,2X,I4,2X,I3)') eigw(lm1,1), abs(eigw(lm1,1)), 1, lm1

      do ikp=1,nsteps
        lm2 = connection(lm1,ikp+1)
        write(368,'(3ES25.16,2X,I4,2X,I3)') eigw(lm2,ikp+1), abs(eigw(lm2,ikp+1)), ikp+1, lm2
      end do!ikp
      close(368)
    end do!itest



!##############################################
!### test the subroutine 'roots_along_edge' ###
!##############################################
    call IoInput('ROOTTSTPS ',uio,1,7,ierr)
    read(unit=uio,fmt=*) nrootsteps
    write(*,'("rootsteps= ",I0)') nrootsteps

    call IoInput('ROOTTITER ',uio,1,7,ierr)
    read(unit=uio,fmt=*) nrootiter
    write(*,'("rootiter= ",I0)') nrootiter

    call IoInput('ROOTTTYPE ',uio,1,7,ierr)
    read (unit=uio,fmt=*) rootinp
    select case (trim(rootinp))
      case('real'); roottype=ROOT_REAL
      case('imag'); roottype=ROOT_IMAG
      case default; stop 'case for roottype not known: real/imag'
    end select
    write(*,'("rootinp= ",(A),", FLAG=",I0)') rootinp, roottype

    call IoInput('ROOTTEPS  ',uio,1,7,ierr)
    read (unit=uio,fmt=*) rooteps
    write(*,'("rooteps= ",ES18.9)') rooteps

    kends(:,1) = k1
    kends(:,2) = k2
    write(*,'("kends(:,1)= (",3ES18.9,")")') kends(:,1)
    write(*,'("kends(:,2)= (",3ES18.9,")")') kends(:,2)
    open(unit=136,file='test.roots_along_edge.txt',form='formatted',action='write')
    call roots_along_edge( inc, lattice, cluster, tgmatrx, nrootsteps, kends, nrootiter, &
                         & roottype, rooteps, nroot, lmroot, kroot, LVroot, RVroot,      &
                         & eigwroot, 136                                                 )
    close(136)


  end subroutine testpath





  subroutine connect_eigw_in_substeps( inc, lattice, cluster, tgmatrx, nsteps, kends, &
                                     & connection, ksub, eigw, LVeig, RVeig           )
  !===========================================================================================!
  ! Connect the eigenvalues between two k-points, kends(:,1) and kends(:,2).                  !
  ! The connection is made by dividing the intervall into nsteps subintervals,                !
  !   and comparing the eigenvectors of each subinterval.                                     !
  !                                                                                           !
  ! The subroutine returns:                                                                   !
  !  - the connection-array, defined by:                                                      !
  !      eigw(lm1, ikp=1) == eigw(connection(lm1,1),1) <--> eigw(connection(lm1,ikp),ikp)     !
  !  - the kpoints of the submesh                                                             !
  !  - the eigenvalues and correcly normalized (left and right) eigenvectors                  !
  !===========================================================================================!

    use type_inc
    use type_data,    only: lattice_type, cluster_type, tgmatrx_type
    use mod_kkrmat,   only: compute_kkr_eigenvectors
    use mod_eigvects, only: compare_eigv
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx
    integer,            intent(in) :: nsteps
    double precision,   intent(in) :: kends(3,2)

    integer,          intent(out) :: connection(inc%almso,nsteps+1)
    double precision, intent(out) :: ksub(3,nsteps+1)
    double complex,   intent(out) :: LVeig(inc%almso,inc%almso,nsteps+1),&
                                   & RVeig(inc%almso,inc%almso,nsteps+1),&
                                   & eigw(inc%almso,nsteps+1)

    !--- Locals ---

    integer :: ikp, ierr, lm1, lmconn, istep

    ! helper variables
    double precision :: dkp(3)

    dkp = (kends(:,2) - kends(:,1))/nsteps
    do ikp=1,nsteps+1

      !--- create submesh to divide the input-intervall ---!
      ksub(:,ikp) = kends(:,1) + dble(ikp-1)*dkp

      !--- find the KKR-matrix, -eigenvectors and ---!
      !---  -eigenvalues on each sub-kpoint       ---!
      call compute_kkr_eigenvectors( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                   & tgmatrx%tinvll(:,:,:,1), ksub(:,ikp),        &
                                   & eigw(:,ikp), LVeig(:,:,ikp), RVeig(:,:,ikp)  )
!     call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
!                                   & tgmatrx%tmat(:,:,1), ksub(:,ikp),            &
!                                   & eigw(:,ikp), LVeig(:,:,ikp), RVeig(:,:,ikp)  )
    end do!ikp



    !-----------------------------------------------------!
    !--- for each step: make connection of eigenvalues ---!
    !--- between 2 k-points by comparing eigenvectors. ---!
    !-----------------------------------------------------!

    !initialize the connections
    connection = 0
    do lm1=1,inc%almso
        connection(lm1,1) = lm1
    end do!lm1

    !compare eigenvectors and make the connection
    do istep=1,nsteps
      do lm1=1,inc%almso
        lmconn = connection(lm1,istep)
        call compare_eigv(inc, LVeig(:,lmconn,istep), RVeig(:,:,istep+1), connection(lm1,istep+1))
      end do!lm1
    end do!istep


  end subroutine connect_eigw_in_substeps





  subroutine roots_along_edge( inc, lattice, cluster, tgmatrx, nsteps, kends, niterate, &
                             & roottype, rooteps, nroot, lmroot, kroot, LVroot,         &
                             & RVroot, eigwroot, uio, eigwends                          )
  !===========================================================================================!
  ! Find the roots of the KKR eigenvalues along a line from kends(:,1) to kends(:,2).         !
  !                                                                                           !
  ! In a first step, the line is divided into nedgesteps subintervals (see routine            !
  !   connect_eigw_in_substeps for details). Then, an iterative scheme (nested intervals) is  !
  !   used to refine the position of the root on the line.                                    !
  !   If uio>0, then information if written to the io-unit 'uio'.                             !
  !                                                                                           !
  ! The subroutine returns:                                                                   !
  !  - nroot:         the number of roots on the edge.                                        !
  !  - lmroot(j):     the number of eigenvalue belonging to the root at the j-th k-point.     !
  !  - kroot(:,:,j):  the j-th k-point where an eigenvalue becomes zero.                      !
  !  - LVroot(:,:,j): all left EVs at the j-th kpoint.                                        !
  !  - RVroot(:,:,j): all right EVs at the j-th kpoint.                                       !
  !  - eigwroot(:,j): all eigenvalues at the j-th kpoint.                                     !
  !  - eigwends(i,j): the eigenvalues of the band belonging to the j-th k-point at the start  !
  !===========================================================================================!

    use type_inc
    use type_data,    only: lattice_type, cluster_type, tgmatrx_type
    use mod_kkrmat,   only: compute_kkr_eigenvectors
    use mod_eigvects, only: compare_eigv
    use mod_mathtools,only: bubblesort
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx
    integer,            intent(in) :: nsteps, niterate, uio
    double precision,   intent(in) :: kends(3,2)
    integer,            intent(in) :: roottype
    double precision,   intent(in) :: rooteps
    integer,            intent(out) :: nroot, lmroot(inc%nrootmax)
    double precision,   intent(out) :: kroot(3,inc%nrootmax)
    double complex,     intent(out) :: LVroot(inc%almso,inc%almso,inc%nrootmax), &
                                     & RVroot(inc%almso,inc%almso,inc%nrootmax), &
                                     & eigwroot(inc%almso,inc%nrootmax)
    double complex, intent(out), optional :: eigwends(2,inc%nrootmax)

    integer          :: connection(inc%almso,nsteps+1), root_interval(inc%almso)
    double precision :: kpoints_steps(3,nsteps+1)
    double complex   :: LV_steps(inc%almso,inc%almso,nsteps+1),&
                      & RV_steps(inc%almso,inc%almso,nsteps+1),&
                      & eigw_steps(inc%almso,nsteps+1)

    ! for the iterative scheme
    integer          :: iter, step_intersect, counter
    double precision :: xtmp, dtmp, kleft(3), kright(3), kmiddle(3)
    double complex   :: LVleft(inc%almso), &
                      & RVleft(inc%almso), &
                      & LVright(inc%almso),&
                      & RVright(inc%almso),&
                      & eigwleft,          &
                      & eigwright,         &
                      & eigwmiddle(inc%almso),        &
                      & LVmiddle(inc%almso,inc%almso),&
                      & RVmiddle(inc%almso,inc%almso)

    integer :: lmsort(inc%almso), loopstep, lmi, lm1, ierr, newconnect
    double precision :: dtmpsort(inc%almso)

    call connect_eigw_in_substeps( inc, lattice, cluster, tgmatrx, nsteps, kends,            &
                                 & connection, kpoints_steps, eigw_steps, LV_steps, RV_steps )

    !determine which eigenvalues to scan
    if(inc%ndegen==2)then
      !in case of a degeneracy of each eigenvalue, take only one of each pair
      dtmpsort = abs(eigw_steps(:,1))
      call bubblesort(inc%almso, dtmpsort, lmsort) ! returns permutation in sortindx
      loopstep = 2 !skip every second eigenvalue (see loop below)
    else
      do lm1=1,inc%almso
        lmsort(lm1) = lm1
      end do!lm1
      loopstep = 1
    end if

    root_interval = 0
    !find which lm intersects with zero and store the interval-number
    do lmi=1,inc%almso,loopstep
      call find_roots_lm_eigw(inc%almso,nsteps,connection,eigw_steps,roottype,1d0,lmsort(lmi),root_interval(lmi))
    end do!lmi

    if(uio>0) write(uio,'("Number of possible roots=" ,I0)') count(root_interval>0)

    !find the precise k-point by iteration
    counter = 0
    lmloop: do lmi=1,inc%almso,loopstep

      !pick variables
      lm1 = lmsort(lmi)
      step_intersect = root_interval(lmi)
      if(uio>0) write(uio,'(2X,"lm= ",I0,", interval= ",I0,", abs(eigw)= ",ES18.9)') lm1, step_intersect, abs(eigw_steps(lm1,1))
      if(step_intersect==0) cycle !no intersection of this eigenvalue with the k-point

      ! pick kpoints, eigenvectors and eigenvalues of this subinterval
      kleft(:)  = kpoints_steps(:,step_intersect)
      LVleft(:) = LV_steps(:,connection(lm1,step_intersect),step_intersect)
      RVleft(:) = RV_steps(:,connection(lm1,step_intersect),step_intersect)
      eigwleft  = eigw_steps(connection(lm1,step_intersect),step_intersect)

      kright(:)  = kpoints_steps(:,step_intersect+1)
      LVright(:) = LV_steps(:,connection(lm1,step_intersect+1),step_intersect+1)
      RVright(:) = RV_steps(:,connection(lm1,step_intersect+1),step_intersect+1)
      eigwright  = eigw_steps(connection(lm1,step_intersect+1),step_intersect+1)


      !================ START ===============!
      !--- Refine the k-point on the edge ---!
      !======================================!
      do iter=1,niterate

        !--------------------------------------------------------!
        ! STEP b): determine the interpolated intersection-point !
        !--------------------------------------------------------!
        select case (roottype)
          case(ROOT_REAL)
            xtmp = dble(eigwleft)/(dble(eigwleft)-dble(eigwright))
          case(ROOT_IMAG)
            xtmp = aimag(eigwleft)/(aimag(eigwleft)-aimag(eigwright))
          case default; stop 'option not known: only real/imag'
        end select!roottype

        kmiddle = (1d0-xtmp)*kleft + xtmp*kright

        !--------------------------------------------------------------------!
        ! STEP c): determine the (true) eigenvalue at the interpolated point !
        !--------------------------------------------------------------------!

        call compute_kkr_eigenvectors( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                     & tgmatrx%tinvll(:,:,:,1), kmiddle,            &
                                     & eigwmiddle, LVmiddle, RVmiddle               )
!       call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
!                                     & tgmatrx%tmat(:,:,1), kmiddle,                &
!                                     & eigwmiddle, LVmiddle, RVmiddle               )

        ! make the connection to the left end of the interval
        newconnect=0
        call compare_eigv( inc, LVleft, RVmiddle, newconnect )



        !-------------------------------------!
        ! STEP c): determine the new interval !
        !-------------------------------------!
        select case (roottype)
          case(ROOT_REAL)
            dtmp = dble(eigwleft)*dble(eigwmiddle(newconnect))
          case(ROOT_IMAG)
            dtmp = aimag(eigwleft)*aimag(eigwmiddle(newconnect))
          case default; stop 'option not known: only real/imag'
        end select!roottype

        if(dtmp<0d0)then
          kright(:)  = kmiddle
          LVright(:) = LVmiddle(:,newconnect)
          RVright(:) = RVmiddle(:,newconnect)
          eigwright  = eigwmiddle(newconnect)
        else
          kleft(:)  = kmiddle
          LVleft(:) = LVmiddle(:,newconnect)
          RVleft(:) = RVmiddle(:,newconnect)
          eigwleft  = eigwmiddle(newconnect)
        end if

        if(uio>0) write(uio,'(4X,"iter= ",I0,", kmiddle= (",3ES25.16,"), eigw_min=(",2ES18.9,"), abs(eigw_min)= ",ES18.9,", lm_min=",I0)') iter, kmiddle, eigwmiddle(newconnect), abs(eigwmiddle(newconnect)), newconnect

      end do!iter
      !======================================!
      !--- Refine the k-point on the edge ---!
      !=============== END ==================!

      !save results
      if(abs(eigwmiddle(newconnect))<rooteps)then
        counter = counter+1
        if(counter>inc%nrootmax)then
          counter = counter-1
          write(*,*) 'Attention: counter>nrootmax is reached'
          exit lmloop
        end if!counter>inc%nrootmax

        kroot(:,counter)  = kmiddle
        LVroot(:,:,counter) = LVmiddle(:,:)
        RVroot(:,:,counter) = RVmiddle(:,:)
        eigwroot(:,counter) = eigwmiddle(:)
        lmroot(counter)     = newconnect
        if(present(eigwends))then
          eigwends(1,counter) = eigw_steps(connection(lm1,1),1)
          eigwends(2,counter) = eigw_steps(connection(lm1,nsteps+1),nsteps+1)
        end if
      end if!abs(dtmp)<rooteps

    end do lmloop!lmi

    nroot=counter
    if(uio>0) write(uio,'("Number of roots found matching the criteria: ", I0)') nroot
    if(uio>0 .and. nroot>0) write(uio,'(2X,"lm= ",8I8)') lmroot(1:nroot)

  end subroutine roots_along_edge





  subroutine find_roots_any_eigw(almso,nsteps,connection,eigw,roottype,rooteps,anyroot)

    implicit none

    integer,          intent(in) :: almso, nsteps, connection(almso,nsteps+1)
    double complex,   intent(in) :: eigw(almso,nsteps+1)
    integer,          intent(in) :: roottype
    double precision, intent(in) :: rooteps

    logical, intent(out) :: anyroot

    integer          :: istep, lm1
    double precision :: dtmp1, dtmp2, xscal
    double complex   :: ctmp1, ctmp2

    !search for an root on the edge
    anyroot=.false.
    outer: do istep=1,nsteps
      inner: do lm1=1,almso


        ctmp1 = eigw(connection(lm1,istep), istep)
        ctmp2 = eigw(connection(lm1,istep+1), istep+1)

        select case (roottype)
          case(ROOT_ANY)
            dtmp1 =  dble(ctmp1) *  dble(ctmp2)
            dtmp2 = aimag(ctmp1) * aimag(ctmp2)
            if(dtmp1<0d0 .or.  dtmp2<0d0 )then
              anyroot = .true.
              exit outer
            end if
          case(ROOT_REAL)
            dtmp1 =  dble(ctmp1) *  dble(ctmp2)
            xscal = -dble(ctmp1) / (dble(ctmp2) - dble(ctmp1))
            dtmp2 = (aimag(ctmp2)-aimag(ctmp1))*xscal + aimag(ctmp1)
            if(dtmp1<0d0 .and.  abs(dtmp2)<rooteps )then
              anyroot = .true.
              exit outer
            end if
          case(ROOT_IMAG)
            dtmp1 =  aimag(ctmp1) *  aimag(ctmp2)
            xscal = -aimag(ctmp1) / (aimag(ctmp2) - aimag(ctmp1))
            dtmp2 = (dble(ctmp2)-dble(ctmp1))*xscal + dble(ctmp1)
            if(dtmp1<0d0 .and.  abs(dtmp2)<rooteps )then
              anyroot = .true.
              exit outer
            end if
          case default; stop 'option not known'
        end select!roottype

      end do inner!lm1
    end do outer!istep

  end subroutine find_roots_any_eigw





  subroutine find_roots_lm_eigw(almso,nsteps,connection,eigw,roottype,rooteps,lm1,root_interval)

    implicit none

    integer,          intent(in) :: almso, nsteps, connection(almso,nsteps+1)
    double complex,   intent(in) :: eigw(almso,nsteps+1)
    integer,          intent(in) :: roottype, lm1
    double precision, intent(in) :: rooteps

    integer, intent(out) :: root_interval

    integer          :: istep, icount
    double precision :: dtmp1, dtmp2, xscal
    double complex   :: ctmp1, ctmp2

    !search for an root on the edge
    root_interval = 0
    icount = 0
    do istep=1,nsteps

        ctmp1 = eigw( connection(lm1,istep),   istep   )
        ctmp2 = eigw( connection(lm1,istep+1), istep+1 )

        select case (roottype)
          case(ROOT_ANY)
            dtmp1 =  dble(ctmp1) *  dble(ctmp2)
            dtmp2 = aimag(ctmp1) * aimag(ctmp2)
            if(dtmp1<0d0 .or.  dtmp2<0d0 ) then
              root_interval = istep
              icount = icount+1
            end if
          case(ROOT_REAL)
            dtmp1 =  dble(ctmp1) *  dble(ctmp2)
            xscal = -dble(ctmp1) / (dble(ctmp2) - dble(ctmp1))
            dtmp2 = (aimag(ctmp2)-aimag(ctmp1))*xscal + aimag(ctmp1)
            if(dtmp1<0d0 .and.  abs(dtmp2)<rooteps )then
              root_interval = istep
              icount = icount+1
            end if
          case(ROOT_IMAG)
            dtmp1 =  aimag(ctmp1) *  aimag(ctmp2)
            xscal = -aimag(ctmp1) / (aimag(ctmp2) - aimag(ctmp1))
            dtmp2 = (dble(ctmp2)-dble(ctmp1))*xscal + dble(ctmp1)
            if(dtmp1<0d0 .and.  abs(dtmp2)<rooteps )then
              root_interval = istep
              icount = icount+1
            end if
          case default; stop 'option not known'
        end select!roottype

    end do!istep

    if(icount>1) root_interval = 0

  end subroutine find_roots_lm_eigw





  subroutine compare_two_eigv_in_substeps( inc, lattice, cluster, tgmatrx, nsteps, &
                                         & kends, lm_in, LV1all, RV2ref, eigw2ref, &
                                         & match                                   )

    use type_inc
    use type_data,     only: lattice_type, cluster_type, tgmatrx_type
    use mod_kkrmat,    only: compute_kkr_eigenvectors
    use mod_eigvects,  only: compare_eigv
    use mod_mathtools, only: bubblesort
    implicit none

    type(inc_type),     intent(in)  :: inc
    type(lattice_type), intent(in)  :: lattice
    type(cluster_type), intent(in)  :: cluster
    type(tgmatrx_type), intent(in)  :: tgmatrx
    integer,            intent(in)  :: nsteps, lm_in
    double precision,   intent(in)  :: kends(3,2)
    double complex,     intent(in)  :: LV1all(inc%almso,inc%almso),&
                                     & RV2ref(inc%almso),&
                                     & eigw2ref
    logical,            intent(out) :: match

    ! kkr-matrix
    double precision :: ksub(3)
    double complex   :: LVleft(inc%almso),&
                      & LVeig(inc%almso,inc%almso),&
                      & RVeig(inc%almso,inc%almso),&
                      & eigw(inc%almso)

    !sorting
    integer          :: lm2, sorted(inc%almso)
    double complex   :: cproj
    double precision :: projs(inc%almso), diffs(inc%almso)

    integer          :: ierr, istep, lm0

!   write(999,'("kends(:,1)=",3ES18.9)') kends(:,1)
!   write(999,'("kends(:,2)=",3ES18.9)') kends(:,2)

    lm0 = lm_in
    LVleft = LV1all(:,lm0)
    LVeig = LV1all

    do istep=1,nsteps+1
      ksub = kends(:,1) + dble(istep)/(nsteps+1)*( kends(:,2) - kends(:,1) )

      call compute_kkr_eigenvectors( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                   & tgmatrx%tinvll(:,:,:,1), ksub, eigw, LVeig, RVeig)
!     call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1), &
!                                   & tgmatrx%tmat(:,:,1), ksub, eigw, LVeig, RVeig )

      ! find the corresponding eigenvctor to LVleft
      call compare_eigv(inc, LVleft, RVeig, lm0)

      LVleft = LVeig(:,lm0)
!     write(999,'(2ES18.9)') eigw(lm0)

    end do!istep


    !====================================!
    !=== find the FINAL corresponding ===!
    !===      eigenvctor to RV2       ===!
    !====================================!

!   !calculate projections on other eigenvectors
!   projs  = 0d0
!   do lm2=1,inc%almso
!       cproj = dot_product(LVeig(:,lm2),RV2ref(:))
!       projs(lm2)  = dble(cproj)**2 + dimag(cproj)**2
!   end do!lm2

!   !find if corresponding eigenvector is matching
!   call bubblesort(inc%almso, projs, sorted)
!   match = .false.
!   if(any(sorted(inc%almso+1-inc%ndegen:inc%almso)==lm0)) match=.true.

    !find if the corresponding eigenvalue is matching
    diffs = abs( eigw-eigw2ref )
    call bubblesort(inc%almso, diffs, sorted)

    match = .false.
    if(any(sorted(1:inc%ndegen)==lm0)) match=.true.

  end subroutine compare_two_eigv_in_substeps







  subroutine unroll_ixyz(ixyz, nCub3, ii3)
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !+ helper function to divide the UID to an integer triple +
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    implicit none

    integer, intent(in)  :: ixyz, nCub3(3)
    integer, intent(out) :: ii3(3)

    integer :: irest

    ii3(3) = int((ixyz-1)/(nCub3(1)*nCub3(2)))
    irest= ixyz-1-ii3(3)*(nCub3(1)*nCub3(2))
    ii3(2) = int(irest/nCub3(1))
    irest= ixyz-1-ii3(3)*(nCub3(1)*nCub3(2)) - ii3(2)*nCub3(1)
    ii3(1) = irest

  end subroutine unroll_ixyz


  function get_cubesinfo_filename(nCub3,lintermediate) result(filename)
    use mod_ioformat,   only: filename_cubesinfo, ext_formatted, fmt_fn_ext
    implicit none
    integer, intent(in) :: nCub3(3)
    logical, intent(in) :: lintermediate
    character(len=256)  :: filename
    character(len=256)  :: fileprefix

    !create filename
    if(lintermediate) then
      write(fileprefix,'(A,3("_",I0))') filename_cubesinfo, nCub3
    else!ltmp
      fileprefix= filename_cubesinfo
    endif!ltmp
     
    write(filename,fmt_fn_ext) trim(fileprefix), ext_formatted

  end function get_cubesinfo_filename


  subroutine save_cubesfile(nCub3, nmarked, imarked, lintermediate)

    implicit none
    integer, intent(in) :: nCub3(3), nmarked, imarked(nmarked)
    logical, intent(in), optional :: lintermediate
    character(len=256) :: filename
    logical :: ltmp
    integer, parameter :: iounit=1566

    ltmp=.false.
    if(present(lintermediate)) ltmp=lintermediate

    filename = get_cubesinfo_filename(nCub3,ltmp)

    open(unit=iounit, file=filename, form='formatted', action='write')
     write(iounit,'(3I8)')  nCub3
     write(iounit,'(I8)')   nmarked
     write(iounit,'(10I8)') imarked
    close(iounit)

  end subroutine save_cubesfile





  subroutine read_cubesfile(nCub3, nmarked, imarked, nCub3_intermediate)

    use mod_mympi, only: myrank, nranks, master
    use mod_ioformat,   only: filename_cubesinfo, ext_formatted, fmt_fn_ext
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer, intent(in), optional :: nCub3_intermediate(3)
    integer, intent(out) :: nCub3(3), nmarked
    integer, allocatable, intent(out) :: imarked(:)
    integer :: ierr, ntmp3(3)=0
    logical :: ltmp
    character(len=256)  :: filename
    integer, parameter  :: iounit=15668

    if(myrank==master)then

      if(present(nCub3_intermediate))then
        filename = get_cubesinfo_filename(nCub3_intermediate,.true.)
      else!present(nCub3_intermediate)
        filename = get_cubesinfo_filename(ntmp3,.false.)
      end if!present(nCub3_intermediate)

      open(unit=iounit, file=filename, form='formatted', action='read')
      read(iounit,'(3I8)')  nCub3
      read(iounit,'(I8)')   nmarked
      allocate(imarked(nmarked), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating imarked in readin on master'
      read(iounit,'(10I8)') imarked
      close(iounit)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(nCub3, 3, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in Bcast(nCub3)'
    call MPI_Bcast(nmarked, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in Bcast(nmarked)'
    if(myrank/=master)then
      allocate(imarked(nmarked), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating imarked in readin on slaves'
    end if!myrank/=master
    call MPI_Bcast(imarked, nmarked, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in Bcast(imarked)'
#endif

  end subroutine read_cubesfile







  subroutine mark_cubes_FScross(inc, lattice, cluster, tgmatrx, nCub3, nsteps, nverts, nedges, diag_ids, bounds, roottype, rooteps, nmarked, imarked)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_mympi,      only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx
    integer,            intent(in) :: nCub3(3), nsteps, nverts, nedges, diag_ids(2,nedges)
    double precision,   intent(in) :: bounds(3,2)
    integer,            intent(in) :: roottype
    double precision,   intent(in) :: rooteps
    integer,              intent(inout) :: nmarked
    integer, allocatable, intent(inout) :: imarked(:)

    integer :: ncubmarked_tmp, ioff, nkpt, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
    integer :: ncubmarked_pT(0:nranks-1), iioff(0:nranks-1)
    double precision :: kverts(3,nverts)
    integer, allocatable :: icubmarked_tmp(:)

    integer          :: connection(inc%almso,nsteps+1)
    double precision :: ksub(3,nsteps+1), kends(3,2)
    double complex   :: LVeig(inc%almso,inc%almso,nsteps+1),&
                      & RVeig(inc%almso,inc%almso,nsteps+1),&
                      & eigw(inc%almso,nsteps+1)

    integer :: iedge, irank, icub, ierr
    logical :: edgeroot(nedges)

    !Parallelize
    call distribute_linear_on_tasks(nranks,myrank,master, nmarked, ntot_pT, ioff_pT, .false.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !Create temp arrays
    allocate(icubmarked_tmp(nkpt), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating icubmarked_tmp in mark_cubes_FScross'
    ncubmarked_tmp = 0
    icubmarked_tmp = -1

    !Test whether cubes cross FS
    do icub=1,nkpt
      if(mod(icub,10)==0 .and. myrank==master) write(*,"(2X,F8.3,A)") dble(icub)/nkpt*100, ' percent done'
      select case (nverts)
        case (8); call generate_cubevertices(nCub3,imarked(icub+ioff), bounds, kverts)
        case (4); call generate_squarevertices(nCub3,imarked(icub+ioff), bounds, kverts)
        case default; stop 'Case nverts not allowed in mark_cubes_FScross'
      end select!
      edgeroot = .false.
      do iedge=1,nedges
        kends(:,1) = kverts(:,diag_ids(1,iedge))
        kends(:,2) = kverts(:,diag_ids(2,iedge))
        call connect_eigw_in_substeps( inc, lattice, cluster, tgmatrx, nsteps, kends, &
                                     & connection, ksub, eigw, LVeig, RVeig           )
        call find_roots_any_eigw( inc%almso, nsteps, connection, eigw, &
                                & roottype, rooteps, edgeroot(iedge)   )
      end do!iedge
      if(any(edgeroot))then
        ncubmarked_tmp = ncubmarked_tmp + 1
        icubmarked_tmp(ncubmarked_tmp) = imarked(icub+ioff)
      end if!l_cut
    end do!icub

    deallocate(imarked)

#ifdef CPP_MPI
    !Combine results
    call MPI_Allgather(ncubmarked_tmp, 1, MPI_INTEGER, ncubmarked_pT, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgather(ncubmarked_tmp) in mark_cubes_FScross'

    nmarked = sum(ncubmarked_pT)
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_cubes_FScross'

    iioff = 0
    do irank=1,nranks-1
      iioff(irank) = sum(ncubmarked_pT(0:irank-1))
    end do

    call MPI_Allgatherv( icubmarked_tmp, ncubmarked_tmp, MPI_INTEGER, &
                       & imarked, ncubmarked_pT, iioff, MPI_INTEGER,  &
                       & MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgatherv(icubmarked_tmp) in mark_cubes_FScross'

#else
    nmarked=ncubmarked_tmp
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_cubes_FScross'
    imarked = icubmarked_tmp(1:ncubmarked_tmp)
#endif

  end subroutine mark_cubes_FScross





  subroutine generate_cubevertices(nCub3,cubeid,bounds,kverts)
    !generates the vertices of a cube given by the cubeid
    implicit none

    integer, intent(in) :: nCub3(3), cubeid
    double precision, intent(in)  :: bounds(3,2)
    double precision, intent(out) :: kverts(3,8)

    integer :: ii3(3), ivert

    call unroll_ixyz(cubeid, nCub3, ii3)

    !generate vertices, limits are from 0 to 1
    kverts(:,1) = dble( (/ ii3(1)  , ii3(2)  , ii3(3)   /) )/nCub3
    kverts(:,2) = dble( (/ ii3(1)+1, ii3(2)  , ii3(3)   /) )/nCub3
    kverts(:,3) = dble( (/ ii3(1)  , ii3(2)+1, ii3(3)   /) )/nCub3
    kverts(:,4) = dble( (/ ii3(1)+1, ii3(2)+1, ii3(3)   /) )/nCub3
    kverts(:,5) = dble( (/ ii3(1)  , ii3(2)  , ii3(3)+1 /) )/nCub3
    kverts(:,6) = dble( (/ ii3(1)+1, ii3(2)  , ii3(3)+1 /) )/nCub3
    kverts(:,7) = dble( (/ ii3(1)  , ii3(2)+1, ii3(3)+1 /) )/nCub3
    kverts(:,8) = dble( (/ ii3(1)+1, ii3(2)+1, ii3(3)+1 /) )/nCub3

    !shift the vertices to the correct bounds
    do ivert=1,8
      kverts(:,ivert) = kverts(:,ivert)*(bounds(:,2)-bounds(:,1))+bounds(:,1)
    end do!ivert

  end subroutine generate_cubevertices




  subroutine generate_squarevertices(nCub3,cubeid,bounds,kverts)
    !generates the vertices of a square given by the id
    implicit none

    integer, intent(in) :: nCub3(3), cubeid
    double precision, intent(in)  :: bounds(3,2)
    double precision, intent(out) :: kverts(3,4)

    integer :: ii3(3), ivert

    call unroll_ixyz(cubeid, nCub3, ii3)

    !generate vertices, limits are from 0 to 1
    kverts(:,1) = dble( (/ ii3(1)  , ii3(2)  , 0 /) )/nCub3
    kverts(:,2) = dble( (/ ii3(1)+1, ii3(2)  , 0 /) )/nCub3
    kverts(:,3) = dble( (/ ii3(1)  , ii3(2)+1, 0 /) )/nCub3
    kverts(:,4) = dble( (/ ii3(1)+1, ii3(2)+1, 0 /) )/nCub3

    !shift the vertices to the correct bounds
    do ivert=1,4
      kverts(:,ivert) = kverts(:,ivert)*(bounds(:,2)-bounds(:,1))+bounds(:,1)
    end do!ivert

  end subroutine generate_squarevertices



  subroutine save_kpointsfile_vis(nkpts, nkpts_irr, kpoints, nsym, isym, kpt2irr, irr2kpt, vis2int, filenamein)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filemode_vis, filename_fsdata
    implicit none

    integer, intent(in) :: nkpts, nkpts_irr, nsym, isym(nsym), kpt2irr(nkpts), irr2kpt(nkpts_irr)
    double precision, intent(in) :: kpoints(3,nkpts)
    integer, intent(in), optional :: vis2int(nkpts)
    character(len=*), intent(in), optional :: filenamein

    double precision :: kpoints_irr(3,nkpts_irr)

    integer :: ii
    character(len=256) :: filename

    kpoints_irr = kpoints(:,irr2kpt)

    if(present(filenamein))then
      filename=filenamein
    else
      write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_vis, ext_formatted
    end if
    open(unit=326523, file=trim(filename), form='formatted', action='write')
    write(326523,'(3I8)') nkpts, nkpts_irr, nsym
    write(326523,*) '#==== isym ====#'
    write(326523,'(12I8)') isym
    write(326523,*) '#==== irr2kpt ====#'
    write(326523,'(10I8)') irr2kpt
    write(326523,*) '#==== kpt2irr ====#'
    write(326523,'(10I8)') kpt2irr
    write(326523,*) '#==== kpoints ====#'
    write(326523,'(3ES25.16)') kpoints_irr
    if(present(vis2int))then
      write(326523,*) '#==== vis2int ====#'
      write(326523,'(10I8)') vis2int
    end if
    close(326523)

  end subroutine save_kpointsfile_vis





  subroutine save_kpointsfile_int(nkpts, kpoints, areas, nsym, isym, filenamein)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filemode_int, filename_fsdata
    implicit none

    integer, intent(in) :: nkpts, nsym, isym(nsym)
    double precision, intent(in) :: kpoints(3,nkpts), areas(nkpts)
    character(len=*), intent(in), optional :: filenamein

    integer :: ii
    character(len=256) :: filename

    if(present(filenamein))then
      filename=filenamein
    else
      write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_int, ext_formatted
    end if
    open(unit=326523, file=trim(filename), form='formatted', action='write')
    write(326523,'(2I8)') nkpts, nsym
    write(326523,'(12I8)') isym
    do ii=1,nkpts
      write(326523,'(4ES25.16)') kpoints(:,ii), areas(ii)
    end do!ii
    close(326523)

  end subroutine save_kpointsfile_int





  subroutine find_kpoints_irredset( bounds, nkpts, kpoints, nkpts_irr, kpt2irr, irr2kpt )
    use mod_parutils, only: parallel_quicksort
    use mod_mympi, only: myrank, nranks, master
    implicit none

    integer, intent(in)  :: nkpts
    double precision, intent(in) :: bounds(3,2), kpoints(3,nkpts)
    integer, intent(out) :: nkpts_irr
    integer, allocatable, intent(out) :: kpt2irr(:), irr2kpt(:)

    integer :: isort(nkpts), itmplist(nkpts)
    double precision, allocatable :: values(:)

    integer :: ipoint, ierr, ikp
    double precision :: delk(3), kpoint0(3)
    double precision, parameter :: eps = 1.0d-12

    !sort kpoints
    allocate(values(nkpts), kpt2irr(nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating values'
    values = 1d0*(kpoints(1,:)-bounds(1,1)) + 3d0*(kpoints(2,:)-bounds(2,1)) + 5d0*(kpoints(3,:)-bounds(3,1))
    call parallel_quicksort(nranks, myrank, master, nkpts, values, isort)
    deallocate(values)

    !filter the irreducibe number
    nkpts_irr = 1
    ipoint  = isort(1)
    itmplist(1) = ipoint
    kpt2irr(ipoint) = 1
    kpoint0 = kpoints(:,ipoint)
    do ikp=2,nkpts
      ipoint = isort(ikp)
      delk   = kpoints(:,ipoint) - kpoint0
      if(any(abs(delk)>=eps)) then
        nkpts_irr = nkpts_irr+1
        kpoint0  = kpoints(:,ipoint)
        itmplist(nkpts_irr) = ipoint
      end if
      kpt2irr(ipoint) = nkpts_irr
    end do

    allocate( irr2kpt(nkpts_irr), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating irr2kpt etc.'

    irr2kpt = itmplist(1:nkpts_irr)

  end subroutine find_kpoints_irredset



end module mod_fermisurf_basic
