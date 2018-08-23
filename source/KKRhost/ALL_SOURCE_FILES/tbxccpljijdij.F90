!-------------------------------------------------------------------------------
! MODULE: mod_tbxccpljijdij
!> @author Bernd Zimmermann
!-------------------------------------------------------------------------------
module mod_tbxccpljijdij
      Use mod_datatypes, Only: dp
   implicit none

   private
   public :: tbxccpljijdij

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: tbxccpljijdij
   !> @author Bernd Zimmermann
   !----------------------------------------------------------------------------
   subroutine tbxccpljijdij( naezd,natypd,lmmaxd,lmgf0d,natomimpd,iemxd,   & !dimensions
      thetas, phis,                                                        &
      natomimp,atomimp,nofgijD,iqat,rclsimp,                               & !imp-cluser
      ijtabcalc,ijtabcalc_I,ijtabsh,ijtabsym,                              & !shells
      ielast,ez,wez,npol,                                                  & !energies
      dsymll,                                                              & !symmetries
      noq,itoq, ncpa)                                                        !CPA
      !   ********************************************************************
      !   * Bernd: add description here
      !   ********************************************************************

#ifdef CPP_MPI
      use mpi
      use mod_types, only: t_mpi_c_grid
#endif
      use mod_types, only: t_tgmat, t_dtmatJij, t_cpa
      use mod_mympi, only: myrank, master
      use mod_version_info
      use mod_md5sums
      use Constants
      use Profiling
      use mod_cmatmul
      use mod_test

      implicit none

      ! parameters
      integer, parameter :: nalpha = 2

      ! various
      integer, intent(in) :: lmmaxd, lmgf0d, naezd, natypd
      real (kind=dp), dimension(natypd), intent(in) :: phis
      real (kind=dp), dimension(natypd), intent(in) :: thetas

      ! Energy contour
      integer, intent(in) :: ielast, iemxd, npol
      complex (kind=dp), dimension(iemxd), intent(in) :: ez
      complex (kind=dp), dimension(iemxd), intent(in) :: wez
      integer :: ie, ie_end, ie_num
#ifdef CPP_MPI
      integer :: ie_start
#endif

      ! Shell indexing
      integer, intent(in) :: nofgijD
      integer, intent(in) :: natomimp
      integer, intent(in) :: natomimpd
      integer, dimension(natypd), intent(in)    :: iqat
      integer, dimension(natomimpd), intent(in) :: atomimp
      integer, dimension(nofgijD), intent(in)   :: ijtabsh
      integer, dimension(nofgijD), intent(in)   :: ijtabsym
      integer, dimension(nofgijD), intent(in)   :: ijtabcalc
      integer, dimension(nofgijD), intent(in)   :: ijtabcalc_I
      real (kind=dp), dimension(3,natomimpd), intent(in) :: rclsimp

      ! Symmetries
      complex (kind=dp), dimension(lmmaxd,lmmaxd,nsymaxd), intent(in) :: dsymll

      ! CPA arrays
      integer, intent(in) :: ncpa
      integer, dimension(naezd), intent(in) :: noq
      integer, dimension(natypd,naezd), intent(in) :: itoq

      !Green functions and GF-work arrays of size LMSxLMS
      complex (kind=dp), dimension(:,:), allocatable :: w1
      complex (kind=dp), dimension(:,:), allocatable :: w2
      complex (kind=dp), dimension(:,:), allocatable :: w3
      complex (kind=dp), dimension(:,:), allocatable :: Tik
      complex (kind=dp), dimension(:,:), allocatable :: Tjl
      complex (kind=dp), dimension(:,:), allocatable :: gij
      complex (kind=dp), dimension(:,:), allocatable :: gji
      complex (kind=dp), dimension(:,:), allocatable :: gij_proj
      complex (kind=dp), dimension(:,:), allocatable :: gji_proj
      !
      complex (kind=dp), dimension(:,:), allocatable :: jxcijint
      complex (kind=dp), dimension(:,:,:,:), allocatable :: Jijmat
#ifdef CPP_MPI
      complex (kind=dp), dimension(:,:,:,:), allocatable :: Jijmat_tmp
      integer :: ierr,i_all
#endif
      integer, dimension(:,:), allocatable :: indxarr

      integer :: i1,i2,nn,kk,ish,isym,iq,jq,iJ1,iI1,jt,it,kalpha,lalpha
      integer :: irec,lm1,lm2,lm3,istore,ncount
      integer :: nstore,i_stat
      real (kind=dp) :: rsh
      complex (kind=dp) :: csum
      integer, dimension(2*natomimp) :: isort, istoretmp
      real (kind=dp), dimension(3) :: rdiff
      real (kind=dp), dimension(2*natomimp) :: dists
      complex (kind=dp), dimension(4) :: jtmp
      real (kind=dp), parameter :: fpi= 4.0d0*PI
      character(len=1) :: cnt
      character(len=22) :: fmt2
      character(len=13) :: jfnam
      character(len=35) :: jfnam2

      logical :: test
      external :: test

      allocate(w1(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(w1))*kind(w1),'w1','tbxccpljijdij')
      allocate(w2(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(w2))*kind(w2),'w2','tbxccpljijdij')
      allocate(w3(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(w3))*kind(w3),'w3','tbxccpljijdij')
      allocate(gij(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(gij))*kind(gij),'gij','tbxccpljijdij')
      allocate(gji(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(gji))*kind(gji),'gji','tbxccpljijdij')
      allocate(Tik(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(Tik))*kind(Tik),'Tik','tbxccpljijdij')
      allocate(Tjl(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(Tjl))*kind(Tjl),'Tjl','tbxccpljijdij')
      allocate(gij_proj(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(gij_proj))*kind(gij_proj),'gij_proj','tbxccpljijdij')
      allocate(gji_proj(lmmaxd,lmmaxd),stat=i_stat)
      call memocc(i_stat,product(shape(gji_proj))*kind(gji_proj),'gji_proj','tbxccpljijdij')

      !get nstore = dimension of storage array
      nstore=0
      do i2=1,natomimp
         do i1=1,natomimp
            nn = (i1-1)*natomimp+i2 !key for the pair (i,j)
            kk = (i2-1)*natomimp+i1 !key for the pair (j,i)
            if(ijtabcalc_I(nn)==1 .and. i1/=i2)then
               !cross-check if also the block g(ji) was calculated
               if(ijtabcalc(kk)/=1) then
                  write(1337,'(A,I0,",",I0)') 'Gji not calculated for (i,j)=(',i1,i2,') calculated'
                  stop 'Blocks not consistent in TBXCCPLJIJDIJ'
               end if !ijtabcalc(kk)/=1

               iq = atomimp(i1)
               jq = atomimp(i2)
               ! --------------------------------------------- loop over occupants
               do iJ1 = 1,noq(jq)
                  jt = itoq(ij1,jq)
                  do iI1 = 1,noq(iq)
                     it = itoq(ii1,iq)
                     nstore=nstore+1
                  end do !iI1
               end do !iJ1
            end if !ijtabcalc(nn)==1
         end do !i2
      end do !i1

      allocate(Jijmat(nalpha,nalpha,nstore,ielast),stat=i_stat)
      call memocc(i_stat,product(shape(Jijmat))*kind(Jijmat),'Jijmat','tbxccpljijdij')
      allocate(indxarr(4,nstore),stat=i_stat)
      call memocc(i_stat,product(shape(indxarr))*kind(indxarr),'indxarr','tbxccpljijdij')
      Jijmat(:,:,:,:) = czero

      !   write(*,'(A)') ' k, l,it,jt,i1,i2,ie'
      ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES
#ifdef CPP_MPI
      ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
      ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

      do ie_num=1,ie_end
         IE = ie_start+ie_num
#else
      do IE = 1,IELAST
         ie_num = ie
         ie_end = ielast
#endif

         istore=0
         ! PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP Atom pairs
         !loop over pairs (i,j)=(i1,i2)
         do i2=1,natomimp
            do i1=1,natomimp
               !get the pointer to the green-function blocks
               nn = (i1-1)*natomimp+i2 !key for the pair (i,j)
               kk = (i2-1)*natomimp+i1 !key for the pair (j,i)

               !check if block g(ij) was calculated
               if(ijtabcalc_I(nn)==1 .and. i1/=i2)then

                  !=================================================================!
                  ! read in the g(ij) block and transform with appropriate symmetry !
                  !     using the appropiate rotation (similar to ROTGLL routine)   !
                  !=================================================================!
                  ish  = ijtabsh(nn)
                  isym = ijtabsym(nn)

                  irec = ie_num + ie_end*(ish-1)
                  w1   = t_tgmat%gmat(:,:,irec)
                  !----------------------------------------------------------------+
                  ! for REL CASE look if it is a unitary / ANTI - unitary rotation |
                  !----------------------------------------------------------------+
                  ! Comment by Bernd: REL case needs special care and testing. Disabled
                  ! for now.
                  cnt = 'N'
                  !       IF ( .NOT.SYMUNITARY(ISYM) ) CNT = 'T'

                  CALL ZGEMM('C',cnt,lmmaxd,lmmaxd,lmmaxd,cone,&
                  & dsymll(1,1,isym),lmmaxd,w1,      &
                  & lmmaxd,czero,w2,lmmaxd           )

                  CALL ZGEMM('N','N',lmmaxd,lmmaxd,lmmaxd,cone,&
                  & w2,lmmaxd,dsymll(1,1,isym),      &
                  & lmmaxd,czero,w1,lmmaxd           )

                  gij(:,:) = w1(:,:)

                  !=================================================================!
                  ! read in the g(ji) block and transform with appropriate symmetry !
                  !=================================================================!
                  ish  = ijtabsh(kk)
                  isym = ijtabsym(kk)

                  irec = ie_num + ie_end*(ish-1)
                  w1   = t_tgmat%gmat(:,:,irec)
                  !----------------------------------------------------------------!
                  ! for REL CASE look if it is a unitary / ANTI - unitary rotation !
                  !----------------------------------------------------------------!
                  ! Comment by Bernd: REL case needs special care and testing. Disabled
                  ! for now.
                  cnt = 'N'
                  !       IF ( .NOT.SYMUNITARY(ISYM) ) CNT = 'T'

                  call ZGEMM('C',cnt,lmmaxd,lmmaxd,lmmaxd,cone,&
                  & dsymll(1,1,isym),lmmaxd,w1,      &
                  & lmmaxd,czero,w2,lmmaxd           )

                  call ZGEMM('N','N',lmmaxd,lmmaxd,lmmaxd,cone,&
                  & w2,lmmaxd,dsymll(1,1,isym),      &
                  & lmmaxd,czero,w1,lmmaxd           )

                  gji(:,:) = w1(:,:)

                  !-----------------------------------------------------------------!
                  ! ==> calculate the exchange coupling constant J_matrix via       !
                  !     Eq.9 of Ebert+Mankovsky (modified for G instead of tau):    !
                  !      J_ij^{kl} ~ Trace [ dt_i/de_k * Gij * dt_j/de_l * Gji ]    !
                  !   in case of alloy system: perform projection on atom types     !
                  !-----------------------------------------------------------------!
                  iq = atomimp(i1)
                  jq = atomimp(i2)
                  ! --------------------------------------------- loop over occupants
                  do iJ1 = 1,noq(jq)
                     jt = itoq(ij1,jq)
                     do iI1 = 1,noq(iq)
                        it = itoq(ii1,iq)

                        istore=istore+1

                        if(ncpa==0)then
                           gij_proj = gij
                           gji_proj = gji
                        else !ncpa==0
                           irec = ie_num
                           !bernd:       ! ! --> Giq,jq is projected on it,jt ==> Git,jt
                           call CMATMUL(lmmaxd,lmmaxd,gij,t_cpa%dtilts(1,1,jt,irec),w2)
                           call CMATMUL(lmmaxd,lmmaxd,t_cpa%dmatts(1,1,it,irec),w2,gij_proj)
                           ! --> Gjq,iq is projected on jt,it ==> Gjt,it
                           call CMATMUL(lmmaxd,lmmaxd,gji,t_cpa%dtilts(1,1,it,irec),w2)
                           call CMATMUL(lmmaxd,lmmaxd,t_cpa%dmatts(1,1,jt,irec),w2,gji_proj)
                        end if !ncpa==0
                        !             write(*,'(5I4)') it,jt,i1,i2,ie
                        do lalpha=1,nalpha
                           ! get Tjl, rotate to global frame and perform Tjl * Gji
                           Tjl = t_dtmatJij(jt)%dtmat_xyz(:,:,lalpha,ie_num)
                           call rotatematrix(Tjl,thetas(jt),phis(jt),lmgf0d,0)
                           call CMATMUL(lmmaxd,lmmaxd,Tjl,gji_proj,w2)

                           do kalpha=1,nalpha
                              ! get Tik, rotate to global frame and perform Tik * Gij
                              Tik = t_dtmatJij(it)%dtmat_xyz(:,:,kalpha,ie_num)
                              call rotatematrix(Tik,thetas(it),phis(it),lmgf0d,0)
                              call CMATMUL(lmmaxd,lmmaxd,Tik,gij_proj,w1)
                              !perform multiplication Tik * Gij * Tjl * Gji
                              call CMATMUL(lmmaxd,lmmaxd,w1,w2,w3)
                              !calculate Tr[ ... ]
                              csum = czero
                              do lm1=1,lmmaxd
                                 csum = csum + w3(lm1,lm1)
                              end do !lm1
                              !store to array
                              Jijmat(kalpha,lalpha,istore,ie) = csum
                           end do !lalpha
                        end do !kalpha
                     end do !iI1
                  end do !iJ1
               end if !ijtabcalc(nn)==1
            end do !i2
         end do !i1
         ! PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP Atom pairs
      end do !ie
      ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES

      ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !   MPI Communication
#ifdef CPP_MPI
      allocate(Jijmat_tmp(nalpha,nalpha,nstore,ielast),stat=i_stat)
      call memocc(i_stat,product(shape(Jijmat_tmp))*kind(Jijmat_tmp),'Jijmat_tmp','tbxccpljijdij')
      Jijmat_tmp(:,:,:,:) = czero

      lm1 = nalpha*nalpha*nstore*ielast
      call MPI_REDUCE( Jijmat,Jijmat_tmp,lm1,MPI_DOUBLE_COMPLEX, &
         MPI_SUM,master,t_mpi_c_grid%myMPI_comm_at,ierr)
      if(ierr/=MPI_SUCCESS)then
         write(*,*) 'Problem in MPI_REDUCE(Jijmat)'
         stop 'TBXCCPLJIJDIJ'
      end if

      Jijmat = Jijmat_tmp
      i_all=-product(shape(Jijmat_tmp))*kind(Jijmat_tmp)
      deallocate(Jijmat_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'Jijmat_tmp','tbxccpljijdij')

#endif
      ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      !   Perform energy integration and output
      if(myrank==master)then
         allocate(jxcijint(4,nstore),stat=i_stat)
         call memocc(i_stat,product(shape(jxcijint))*kind(jxcijint),'jxcijint','tbxccpljijdij')

         jxcijint = czero

         istore = 0
         do i2=1,natomimp
            do i1=1,natomimp
               nn = (i1-1)*natomimp+i2 !key for the pair (i,j)
               kk = (i2-1)*natomimp+i1 !key for the pair (j,i)
               if(ijtabcalc_I(nn)==1 .and. i1/=i2)then

                  iq = atomimp(i1)
                  jq = atomimp(i2)
                  !-------------------------------------------------------------
                  ! Loop over occupants
                  !-------------------------------------------------------------
                  do iJ1 = 1,noq(jq)
                     jt = itoq(ij1,jq)
                     do iI1 = 1,noq(iq)
                        it = itoq(ii1,iq)
                        !count one up
                        istore=istore+1
                        indxarr(1,istore) = it
                        indxarr(2,istore) = jt
                        indxarr(3,istore) = i1
                        indxarr(4,istore) = i2
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! perform substraction instead of addition
                        ! because WGTE ~ -1/pi (WGTE = WEZ(IE)/NSPIN)
                        ! Write out energy-resorved integrand and integral
                        ! Phivos Mavropoulos 24.10.2012
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if(npol==0 .or. test('Jijenerg'))then
                           fmt2 = '(4(A,I5.5))'
                           write (jfnam2,fmt2) 'Jij_enrg.',it,'.',jt,'.',i1,'.',i2
                           open (499,file=jfnam2,status='unknown')
                           call version_print_header(499,'; '//md5sum_potential//'; '//md5sum_shapefun)
                           write(499,fmt='(a)') &
                           & '# Energy Re,Im ; j(E) Re,Im ; J(E) Re,Im ;&
                           & d(E) Re,Im ; D(E) Re,Im ;&
                           & s(E) Re,Im ; S(E) Re,Im ;&
                           & a(E) Re,Im ; A(E) Re,Im '
                           write(499,fmt='(4(a,i5))') &
                           &  '# IT=',IT,' JT=',JT,' ISITE=',i1,' JSITE=',i2
                           write(499,fmt='(a,i6)') '#ENERGIES: ',ielast
                        endif !test('Jijenerg')
                        !            Jijmat_real = 0d0
                        do ie=1,ielast
                           !                 Jijmat_real(:,:)=Jijmat_real(:,:)+aimag(wez(ie)*Jijmat(:,:,istore,ie))/2d0
                           jtmp(1) = (Jijmat(1,1,istore,ie)+Jijmat(2,2,istore,ie))/2d0
                           !jtmp(2) = (Jijmat(2,1,istore,ie)-Jijmat(1,2,istore,ie))/2d0 !<- old defiition (until Apr.2017) which assumed +Dij.(SixSj)
                           jtmp(2) = (Jijmat(1,2,istore,ie)-Jijmat(2,1,istore,ie))/2d0  !<- changed to -Dij.(SixSj), to be consistent with KKRwiki
                           jtmp(3) = (Jijmat(1,1,istore,ie)-Jijmat(2,2,istore,ie))/2d0
                           jtmp(4) = (Jijmat(2,1,istore,ie)+Jijmat(1,2,istore,ie))/2d0

                           jxcijint(:,istore) = jxcijint(:,istore) - wez(ie)*jtmp/4d0
                           !factor 2 for NSPIN, another factor of 2 to be consistent
                           !with definition in tbxccpljij (different from impurity program)
                           if(npol==0 .or. test('Jijenerg'))then
                              write (499,fmt='(18e12.4)') &
                              & ez(ie),jtmp(1)/fpi,jxcijint(1,istore),&
                              &        jtmp(2)/fpi,jxcijint(2,istore),&
                              &        jtmp(3)/fpi,jxcijint(3,istore),&
                              &        jtmp(4)/fpi,jxcijint(4,istore)
                           endif !test('Jijenerg')
                        end do !ie

                        if(npol==0 .or. test('Jijenerg')) close(499)
                        !            write(34536254,'(A,2I3)') '# coupling between',i1, i2
                        !            write(34536254,'(A)') '# gen. Jij matrix is'
                        !            write(34536254,'(3ES24.15E3)') Jijmat_real(1,1),Jijmat_real(1,2),Jijmat_real(1,3)
                        !            write(34536254,'(3ES24.15E3)') Jijmat_real(2,1),Jijmat_real(2,2),Jijmat_real(2,3)
                        !            write(34536254,'(3ES24.15E3)') Jijmat_real(3,1),Jijmat_real(3,2),Jijmat_real(3,3)
                        !            write(34536254,'(A)') '################################################################'
                     end do !iI1
                  end do !iJ1
               end if !ijtabcalc(nn)==1
            end do !i2
         end do !i1

         do lm1=1,natypd
            if(count(indxarr(1,:)==lm1)>0)then
               loop: do lm3=1,nstore
                  if(indxarr(1,lm3)==lm1)then
                     i1 = indxarr(3,lm3)
                     exit loop
                  end if
               end do loop !lm3

               write(jfnam,'(A,I5.5)') 'Jij.atom', lm1
               open(49,file=jfnam,form='formatted',action='write')
               call version_print_header(49,'; '//md5sum_potential//'; '//md5sum_shapefun)
               WRITE (49,99009) lm1,IQAT(lm1),i1

               do lm2=1,natypd
                  !determine the pairs (ij) that match
                  ncount = 0
                  dists = 1d38
                  istoretmp = 0
                  do istore = 1,nstore
                     if(indxarr(1,istore)==lm1 .and. indxarr(2,istore)==lm2)then
                        i1 = indxarr(3,istore)
                        i2 = indxarr(4,istore)

                        ncount = ncount+1
                        rdiff = rclsimp(:,i2)-rclsimp(:,i1)
                        rsh   = sqrt(sum(rdiff**2))
                        dists(ncount) = rsh
                        istoretmp(ncount) = istore
                     end if
                  end do !istore

                  !sort according to distance
                  if(ncount>0)then
                     call bubblesort(ncount, dists(1:ncount), isort(1:ncount))

                     !output
                     do lm3=1,ncount
                        istore = istoretmp(isort(lm3))
                        i1 = indxarr(3,istore)
                        i2 = indxarr(4,istore)

                        rdiff = rclsimp(:,i2)-rclsimp(:,i1)
                        rsh   = sqrt(sum(rdiff**2))
                        write(49,99010) rsh, aimag(jxcijint(1,istore)),&
                        & aimag(jxcijint(2,istore)),&
                        & aimag(jxcijint(3,istore)),&
                        & aimag(jxcijint(4,istore)),&
                        & rdiff, lm2, i2
                     end do

                     if(lm2/=natypd) write(49,'(A)') '#&'
                  end if
               end do !lm2
               close(49)
            end if !count>0
         end do !i1
      end if !myrank==master
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

#ifdef CPP_MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      99009 FORMAT("# off-diagonal exchange coupling constants ",/,&
      &     "# for atom IT = ",I5," on site IQ = ",I5," impurity site = ",I5,/,&
      &     "# R_IQ,JQ      J_IT,JT     D_IT,JT   S_IT,JT     A_IT,JT       Rvec    JT",/,&
      &     "# ( ALAT )       ( Ry )",/,"#      ")
      99010 FORMAT(F12.8,4E16.8,2X,3F12.8,2X,2I5)

      ! stop 'test stop'
   end subroutine tbxccpljijdij


   subroutine bubblesort(n, Xin, Iout)

      implicit none
      integer,          intent(in)  :: n
      real (kind=dp), intent(in)  :: Xin(n)
      integer,          intent(out) :: Iout(n)

      logical          :: swaped
      integer          :: ii, ntmp, itmp

      swaped = .true.
      ntmp   = n

      do ii=1,n
         Iout(ii) = ii
      end do

      do while ( swaped .and. ntmp>1 )
         swaped = .false.
         do ii=1,ntmp-1
            if( Xin(Iout(ii))>Xin(Iout(ii+1)) ) then
               itmp       = Iout(ii)
               Iout(ii)   = Iout(ii+1)
               Iout(ii+1) = itmp
               swaped = .true.
            end if
         end do !ii
         ntmp = ntmp-1
      end do !while

   end subroutine bubblesort

end module mod_tbxccpljijdij
