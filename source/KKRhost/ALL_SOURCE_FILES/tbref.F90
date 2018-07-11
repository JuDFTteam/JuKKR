!-------------------------------------------------------------------------------
! SUBROUTINE: TBREF
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine TBREF(EZ,IELAST,ALATC,VREF,IEND,LMAX,NCLS,NINEQ,NREF,CLEB,RCLS,ATOM,  &
   CLS,ICLEB,LOFLM,NACLS,REFPOT,RMTREF,TOLRDIF,TMPDIR,ITMPDIR,ILTMP,NAEZ,LLY)

   use mod_mympi, only: myrank, nranks, master
   use mod_types, only: t_tgmat, t_lloyd, t_inc
#ifdef CPP_MPI
   use mod_types, only: t_mpi_c_grid, init_tgmat, init_tlloyd
   use mpi
   use mod_mympi, only: myrank, nranks, master,find_dims_2d,distribute_linear_on_tasks
#endif
   use Constants
   use Profiling
   use global_variables
   use mod_datatypes, Only: dp

   implicit  none
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green
   !          function, set up in the spin-independent non-relativstic
   !          (l,m_l)-representation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     ..
   ! .. Input variables
   integer, intent(in) :: LLY    !< LLY <> 0 : alpply Lloyds formula
   integer, intent(in) :: IEND
   integer, intent(in) :: LMAX   !< Maximum l component in wave function expansion
   integer, intent(in) :: NCLS   !< Number of reference clusters
   integer, intent(in) :: NREF   !< Number of diff. ref. potentials
   integer, intent(in) :: NAEZ   !< Number of atoms in unit cell
   integer, intent(in) :: NINEQ  !< Number of ineq. positions in unit cell
   integer, intent(in) :: IELAST
   real (kind=dp), intent(in) :: ALATC
   real (kind=dp), intent(in) :: TOLRDIF !< For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.
   integer, dimension(NAEZD+NEMBD), intent(in) :: CLS      !< Cluster around atomic sites
   integer, dimension(LM2D), intent(in)      :: LOFLM    !< l of lm=(l,m) (GAUNT)
   integer, dimension(NCLSD), intent(in)     :: NACLS    !< Number of atoms in cluster
   integer, dimension(NAEZD+NEMBD), intent(in) :: REFPOT   !< Ref. pot. card  at position
   integer, dimension(NACLSD,NAEZD+NEMBD), intent(in)   :: ATOM  !< Atom at site in cluster
   integer, dimension(NCLEB,4), intent(in)            :: ICLEB !< Pointer array
   real (kind=dp), dimension(NREF), intent(in) :: VREF
   real (kind=dp), dimension(NREF), intent(in) :: RMTREF   !< Muffin-tin radius of reference system
   real (kind=dp), dimension(NCLEB,2), intent(in) :: CLEB  !< GAUNT coefficients (GAUNT)
   real (kind=dp), dimension(3,NACLSD,NCLSD), intent(in) :: RCLS   !< Real space position of atom in cluster
   ! .. In/Out variables
   integer, intent(inout) :: ILTMP
   integer, intent(inout) :: ITMPDIR
   character(len=80), intent(inout) :: TMPDIR
   complex (kind=dp), dimension(IEMXD), intent(inout) :: EZ
   ! .. Local variables
   integer :: I1,IC,ICLS,IE,LM1,NACLSMAX,LRECGRF1, i_stat,i_all
   complex (kind=dp) :: ERYD
   complex (kind=dp) :: LLY_G0TR_IE    ! LLY
   complex (kind=dp) :: lly_g0tr_dum   ! LLY dummy variable used if no LLY is chosen to save memory
   ! .. Parameters
   integer :: LRECGRF
   ! .. Local Arrays
   complex (kind=dp), dimension(0:LMAX,NREF) :: ALPHAREF  ! LLY Lloyd Alpha matrix
   complex (kind=dp), dimension(0:LMAX,NREF) :: DALPHAREF ! LLY Derivative of the Lloyd Alpha matrix
   complex (kind=dp), dimension(LMGF0D,LMGF0D,NREF) :: TREFLL    ! LLY
   complex (kind=dp), dimension(LMGF0D,LMGF0D,NREF) :: DTREFLL   ! LLY

   ! .. Local allocatable arrays
   complex (kind=dp), dimension(:,:), allocatable :: LLY_G0TR  ! LLY
   complex (kind=dp), dimension(:,:), allocatable :: dginp_dum ! LLY dummy variable used if no LLY is chosen to save memory
   complex (kind=dp), dimension(:,:,:), allocatable :: GINP
   complex (kind=dp), dimension(:,:,:), allocatable :: DGINP
#ifdef CPP_MPI
   ! ..
   ! .. MPI variables
   integer :: ntot1, idim
   integer, dimension(0:nranks-1) :: ntot_pT, ioff_pT
   complex (kind=dp), dimension(:,:,:), allocatable :: work
#endif
   integer :: IE_START, IE_END
   integer :: i1_start, i1_end
   integer :: ie_num, ierr
   !     ..
   !     .. External Functions ..
   logical :: TEST,OPT
   external :: TEST,OPT
   !     ..
   !     .. External Subroutines ..
   external :: CALCTREF13,GLL13
   !     ..

   NACLSMAX = 1
   do IC = 1,NCLS
      if (NACLS(IC).GT.NACLSMAX) NACLSMAX = NACLS(IC)
   enddo
   LRECGRF1 = WLENGTH*4*NACLSMAX*LMGF0D*LMGF0D*NCLS
   LRECGRF=WLENGTH*4*NACLSD*LMGF0D*LMGF0D*NCLSD

   ! allocate and initialize ginp
   allocate(GINP(NACLSMAX*LMGF0D,LMGF0D,NCLS), stat=i_stat)
   call memocc(i_stat,product(shape(GINP))*kind(GINP),'GINP','TBREF')
   GINP = CZERO
   allocate(dginp_dum(NACLSMAX*LMGF0D,LMGF0D), stat=i_stat)
   call memocc(i_stat,product(shape(dginp_dum))*kind(dginp_dum),'dginp_dum','TBREF')
   dginp_dum = CZERO

   if (LLY.NE.0) then
      ! allocate and initialize dginp and lly_g0tr
      allocate(DGINP(NACLSMAX*LMGF0D,LMGF0D,NCLS), stat=i_stat)
      call memocc(i_stat,product(shape(DGINP))*kind(DGINP),'DGINP','TBREF')
      DGINP = CZERO
      allocate(LLY_G0TR(IELAST,NCLSD), stat=i_stat)
      call memocc(i_stat,product(shape(LLY_G0TR))*kind(LLY_G0TR),'LLY_G0TR','TBREF')
      LLY_G0TR = CZERO
   endif

   if (t_tgmat%gref_to_file) then
      call OPENDAFILE(68,'gref',4,LRECGRF1,TMPDIR,ITMPDIR,ILTMP)
   end if
   if (LLY.NE.0) then
      if(t_lloyd%dgref_to_file) then
         call OPENDAFILE(681,'dgrefde',7,LRECGRF1,TMPDIR,ITMPDIR,ILTMP)
      end if
      if(t_lloyd%g0tr_to_file) then
         open(682,FILE='lly_g0tr_ie.ascii',FORM='FORMATTED')
      end if
   endif
   !
   !----------------------------------------------------------------------------
#ifdef CPP_MPI
   if(t_inc%i_write>0) write(1337,*) myrank, 'start tbref e-loop'

    call distribute_linear_on_tasks(t_mpi_c_grid%nranks_at,   &
               t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, &
               master,IELAST,ntot_pT,ioff_pT,.true.,.true.)
     
   ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
   ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

     t_mpi_c_grid%ntot2=ie_end   !t_mpi_c_grid%dims(1)
     if(.not. (allocated(t_mpi_c_grid%ntot_pT2) .or. allocated(t_mpi_c_grid%ioff_pT2))) then 
        allocate(t_mpi_c_grid%ntot_pT2(0:t_mpi_c_grid%nranks_at-1), t_mpi_c_grid%ioff_pT2(0:t_mpi_c_grid%nranks_at-1))
     end if
     t_mpi_c_grid%ntot_pT2 = ntot_pT
     t_mpi_c_grid%ioff_pT2 = ioff_pT

     ! now initialize arrays for tmat, gmat, and gref
     call init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
     if(lly.ne.0) call init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)
#else
   ie_start = 0
   ie_end = IELAST
#endif

   do ie_num = 1,ie_end
      IE = ie_start+ie_num
      if(t_inc%i_write>0) write(1337,FMT='(A23,I5,2F14.8)')'TBREF: GREF for energy:',IE,EZ(IE)

      ! reset arrays to zero (needed in case of MPIatom energy/atom distribution
      ginp(:,:,:) = CZERO
      if(LLY>0) then
         dginp(:,:,:) = CZERO
         LLY_G0TR(IE,:) = CZERO
         dginp_dum = CZERO
         lly_g0tr_dum = CZERO
      endif

      !
      !skip this part with GREENIMP option
      if(opt('GREENIMP')) then
         if(myrank==master) write(*,*)'Skipping atom loop in tbref for energy point',ie,'of',ielast
         i1_start = 1
         i1_end = 0
      else
         i1_start = 1
         i1_end = NREF
      end if
      ERYD = EZ(IE)
      DO I1 = i1_start, i1_end!1,NREF
         CALL CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1,&              ! LLY Lloyd
            TREFLL(1,1,I1),DTREFLL(1,1,I1),                &   ! LLY Lloyd
            ALPHAREF(0,I1),DALPHAREF(0,I1),LMAX+1,LMGF0D)    ! LLY Lloyd

         IF (TEST('flow    ').and.(t_inc%i_write>0))WRITE (1337,FMT=*) 'tll(ref),  i1 = ',I1
      END DO
      !
      IF (TEST('flow    ') .and. (t_inc%i_write>0))WRITE (1337,FMT=*) 't-mat(Ref) o.k.', IE
      ! ----------------------------------------------------------------------
#ifdef CPP_MPI
      ntot1 = NCLS

      if(myrank==master) write(1337,*) 'tbref NCLS loop:', NCLS,t_mpi_c_grid%nranks_ie
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie,  &
         t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,master, &
         ntot1,ntot_pT,ioff_pT,.true.)
      i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
      i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+ntot_pT(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)

      t_mpi_c_grid%ntot_pT1 = ntot_pT
      t_mpi_c_grid%ioff_pT1 = ioff_pT

#else
      i1_start = 1
      i1_end   = NCLS

#endif
      !skip this part with GREENIMP option
      if(opt('GREENIMP')) then
         i1_start = 1
         i1_end = 0
      end if

      do ICLS=i1_start,i1_end
         I1 = 1
         IC = 0
         do while (IC.EQ.0 .AND. I1.LE.NINEQ)
            if (CLS(I1).EQ.ICLS) IC = I1
            I1 = I1 + 1
         end do
         !
         if (IC.EQ.0) STOP 'Error in CLS(*) array in tbref'
         if (TEST('flow    ').and.(t_inc%i_write>0)) then
            write (1337,FMT=*) 'CLUSTER ',ICLS,' at ATOM ',IC
         endif
         !
         call GLL13(ERYD,CLEB(1,2),ICLEB,LOFLM,IEND,TREFLL,DTREFLL,ATOM(1,IC),&
            REFPOT,RCLS(1,1,ICLS),NACLS(ICLS),TOLRDIF,ALATC,0,GINP(1,1,ICLS), &
            dginp_dum,NACLSMAX,lly_g0tr_dum,LLY)

         ! copy from dummy variable to arrays
         if(LLY.ne.0) then
            DGINP(:,:,ICLS) = dginp_dum(:,:)
            LLY_G0TR(IE,ICLS) = lly_g0tr_dum
         end if

      END DO !icls


#ifdef CPP_MPI
      if(t_mpi_c_grid%nranks_ie>1) then
         !gather results of ginp, dginp and lly_g0tr from above parallel loop
         IDIM = NACLSMAX*LMGF0D*LMGF0D*NCLS
         allocate(work(NACLSMAX*LMGF0D,LMGF0D,NCLS), stat=i_stat)
         call memocc(i_stat,product(shape(work))*kind(work),'work','TBREF')

         call MPI_ALLREDUCE(GINP,work,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            t_mpi_c_grid%mympi_comm_ie,IERR)
         call ZCOPY(idim,work,1,ginp,1)

         i_all=-product(shape(work))*kind(work)
         deallocate(work,stat=i_stat)
         call memocc(i_stat,i_all,'work','TBREF')
         if(LLY.NE.0) THEN
            IDIM = NACLSMAX*LMGF0D*LMGF0D*NCLS
            allocate(work(NACLSMAX*LMGF0D,LMGF0D,NCLS), stat=i_stat)
            call memocc(i_stat,product(shape(work))*kind(work),'work','TBREF')

            call MPI_ALLREDUCE(DGINP,work,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM, &
               t_mpi_c_grid%mympi_comm_ie,IERR)
            call ZCOPY(idim,work,1,dginp,1)
            i_all=-product(shape(work))*kind(work)
            deallocate(work,stat=i_stat)
            call memocc(i_stat,i_all,'work','TBREF')
            IDIM = NCLS
            allocate(work(NCLSD,1,1), stat=i_stat)
            call memocc(i_stat,product(shape(work))*kind(work),'work','TBREF')
            call MPI_ALLREDUCE(LLY_G0TR(IE,:),work,IDIM,MPI_DOUBLE_COMPLEX,&
               MPI_SUM,t_mpi_c_grid%mympi_comm_ie,IERR)
            call ZCOPY(idim,work,1,lly_g0tr(ie,:),1)
            i_all=-product(shape(work))*kind(work)
            deallocate(work,stat=i_stat)
            call memocc(i_stat,i_all,'work','TBREF')
         endif
      endif!t_mpi_c_grid%nranks_ie>1
#endif
      if (LLY.NE.0) then                                    ! LLY Lloyd
         LLY_G0TR_IE = CZERO                                ! LLY Lloyd
         do I1 = 1,NAEZ                                     ! LLY Lloyd
            ICLS = CLS(I1)                                  ! LLY Lloyd
            LLY_G0TR_IE = LLY_G0TR_IE + LLY_G0TR(IE,ICLS)   ! LLY Lloyd
         enddo                                              ! LLY Lloyd
      endif                                                 ! LLY Lloyd
      ! ----------------------------------------------------------------------
      if (t_tgmat%gref_to_file) then
#ifdef CPP_MPI
         !make sure only one processor writes for one energy point
         if(t_mpi_c_grid%myrank_ie==0) then
            write (68,REC=IE) GINP                 ! Gref
         end if !if(t_mpi_c_grid%myrank_ie==0) then

         ! human readable writeout if test option is hit
         if(test('fileverb')) then
            ! test writeout
            do i1=0,nranks-1
               if(myrank==i1) then
                  write(686868+myrank,*) myrank,ie,ginp
                  !write(686868+myrank,'(2i9,200000F15.7)') myrank,ie,ginp
               endif
               call MPI_BARRIER(t_mpi_c_grid%mympi_comm_ie,ierr)
            end do
            ! end test writeout
         end if
#else
         WRITE (68,REC=IE) GINP                 ! Gref
#endif
      else   ! (t_tgmat%gref_to_file)
         t_tgmat%gref(:,:,:,ie_num) = GINP(:,:,:)
      end if ! (t_tgmat%gref_to_file)
      ! store result either to file or in memory
      if (LLY.NE.0) then                                ! LLY Lloyd
         if(t_lloyd%dgref_to_file) then                 ! LLY Lloyd
            write (681,REC=IE) DGINP ! dGref/dE         ! LLY Lloyd
            if (test('fileverb')) then
              write(681681,'(i9,200000ES15.7)') ie, DGINP
            end if
         else                                           ! LLY Lloyd
            t_lloyd%dgref(:,:,:,ie_num) = DGINP(:,:,:)  ! LLY Lloyd
         endif                                          ! LLY Lloyd
         if(t_lloyd%g0tr_to_file) then                  ! LLY Lloyd
            write (682,FMT='(2E24.16)') LLY_G0TR_IE     ! LLY Lloyd
            if (test('fileverb')) then
              write(682682,'(i9,200000ES15.7)') ie_num, LLY_G0TR_IE
            end if
         else                                           ! LLY Lloyd
            t_lloyd%g0tr(ie_num) = LLY_G0TR_IE          ! LLY Lloyd
         end if                                         ! LLY Lloyd
      end if ! LLY.NE.0                                 ! LLY Lloyd
      !
      if (TEST('flow    ').and.(t_inc%i_write>0)) write (1337,FMT=*)'G(n,lm,n,lm) (Ref) o.k.'
   end do ! IE
   !----------------------------------------------------------------------------
   if (t_tgmat%gref_to_file) then
      CLOSE (68)
   end if

   if (LLY.NE.0 .and. t_lloyd%dgref_to_file) close(681)
   if (LLY.NE.0 .and. t_lloyd%g0tr_to_file) close(682)
   i_all=-product(shape(GINP))*kind(GINP)
   deallocate(GINP,stat=i_stat)
   call memocc(i_stat,i_all,'GINP','TBREF')
   if(LLY.NE.0) then
      i_all=-product(shape(DGINP))*kind(DGINP)
      deallocate(DGINP,stat=i_stat)
      call memocc(i_stat,i_all,'DGINP','TBREF')

      i_all=-product(shape(LLY_G0TR))*kind(LLY_G0TR)
      deallocate(LLY_G0TR,stat=i_stat)
      call memocc(i_stat,i_all,'LLY_G0TR','TBREF')
   end if

end subroutine TBREF
