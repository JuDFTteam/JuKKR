module mod_rhoqtools

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rhoq_write_kmesh(nofks, nxyz, volbz, bzkp, volcub, recbv, bravais)

      implicit none
   
      integer, intent(in) :: nofks, nxyz(3)
      double precision, intent(in) :: volbz, bzkp(3,nofks), volcub(nofks), recbv(3,3), bravais(3,3)
      ! local
      integer :: ks, i
   
      ! write out kpoints
      open(8888, file='kpts.txt', form='formatted')
      write(8888,'(3I9)') NOFKS, NXYZ(1), NXYZ(2)
      write(8888,'(E16.7)') VOLBZ
      do ks=1,nofks
         write(8888,'(4E16.7)') (BZKP(i,KS), i=1,3), VOLCUB(KS)
      end do
      write(8888,'(100E16.7)') RECBV(1:3,1:3),BRAVAIS(1:3,1:3)
      close(8888)

   end subroutine rhoq_write_kmesh

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rhoq_read_mu0_scoef(iatomimp, mu, nscoef, imin)

#ifdef CPP_MPI
      use mpi
#endif
      use mod_mympi, only: myrank, master
      implicit none

      integer, intent(out) :: mu, nscoef, imin
      integer, allocatable, intent(inout) :: iatomimp(:)
      ! local
      integer :: i1
#ifdef CPP_MPI
      integer :: ierr
#endif

      ! read in mu0
      if(myrank==master) then
         open(8888,file='mu0',form='formatted')
         read(8888,*) mu, nscoef
         allocate(iatomimp(nscoef))
         do i1=1,nscoef
            read(8888,*) iatomimp(i1)
         end do
         close(8888)
      end if
      
#ifdef CPP_MPI
      ! communicate mu and nscoef
      call MPI_Bcast(mu, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'Error Bcast mu0'
      call MPI_Bcast(nscoef,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      if(ierr/=MPI_SUCCESS) stop 'Error Bcast nscoef'
      if(.not.allocated(iatomimp)) allocate(iatomimp(nscoef))
      call MPI_Bcast(iatomimp, nscoef, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'Error Bcast iatomimp'
#endif
      
      !find imin
      imin = 100000
      do i1=1,nscoef
         if(iatomimp(i1)<imin) imin = iatomimp(i1)
      end do
      nscoef = nscoef-1

   end subroutine rhoq_read_mu0_scoef

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rhoq_find_kmask(nofks, k_end, bzkp, kmask, rhoq_kmask)

#ifdef CPP_MPI
      use mpi
#endif
      use mod_mympi, only: myrank, master

      implicit none

      integer, intent(in) :: nofks
      integer, intent(out) :: k_end
      integer, allocatable, intent(out) :: kmask(:)
      double precision, intent(in) :: bzkp(3, nofks)
      double precision, allocatable, intent(out) :: rhoq_kmask(:,:)
      ! local
      integer :: i, j, kpt, kmask_mode, k_start
      logical :: kmask_info
      double precision :: k_mask_bounds(4), recbv(3,3), kp(3)
#ifdef CPP_MPI
      integer :: ierr
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (myrank==master) then
         !read recbv  
         open(8888, file='kpts.txt', form='formatted')
         read(8888,*)
         read(8888,*)
         do kpt=1,nofks
            read(8888,*) 
         end do
         read(8888,'(100E16.7)') RECBV(1:3,1:3)
         close(8888)
         
         allocate(kmask(nofks))
        
         ! determine kmask parameters
         inquire(file='kmask_info.txt', exist=kmask_info)
         if (kmask_info) then
            write(*,*) "found 'kmask_info.txt' file, start reading..."
            open(8888, file='kmask_info.txt', form='formatted')
            read(8888,*) kmask_mode
            if (kmask_mode==1) then ! read R0=(x,y), then R1 and R2 (outer and inner radius around R0)
               write(*,*) "kmask_mode is 1: spherical region"
               read(8888,*) k_mask_bounds(1), k_mask_bounds(2)
               write(*,*) "R0=", k_mask_bounds(1), k_mask_bounds(2)
               read(8888,*) k_mask_bounds(3)
               write(*,*) "R1=", k_mask_bounds(3)
               read(8888,*) k_mask_bounds(4)
               write(*,*) "R2=", k_mask_bounds(4)
            elseif (kmask_mode==2) then ! read xmin, xmax, ymin, ymax of kmask_box
               write(*,*) "kmask_mode is 2: box"
               read(8888,*) k_mask_bounds(1)
               read(8888,*) k_mask_bounds(2)
               read(8888,*) k_mask_bounds(3)
               read(8888,*) k_mask_bounds(4)
               write(*,*) "xmin=", k_mask_bounds(1)
               write(*,*) "xmax=", k_mask_bounds(2)
               write(*,*) "ymin=", k_mask_bounds(3)
               write(*,*) "ymax=", k_mask_bounds(4)
            end if ! kmask_mode 1 or 2
            ! close kmask_info.txt
            close(8888)
         else
            kmask_mode = 0
         end if !kmask_info.txt file found
        
         if (kmask_mode==3) then ! read kmask from file
            write(*,*) "kmask_mode is 3: read 'kpts_mask.txt' file"
            open(8888, file='kpts_mask.txt', form='formatted')
         end if
         
         ! use these as counters
         k_start = 1
         k_end = 0
         
         ! find kmask and number points in box
         do kpt=1,nofks
            !findig kmask
            ! default is take all
            kmask(kpt) = 1
            if (kmask_mode==1) then ! sph mode
               kmask(kpt)= 0
               do i=-1,1,1
                  do j=-1,1,1
                     kp(1:3) = bzkp(1:3,kpt)+i*recbv(1:3,1)+j*recbv(1:3,2)
                     ! first shift kpt to be centered around R0
                     kp(1) = kp(1) - k_mask_bounds(1)
                     kp(2) = kp(2) - k_mask_bounds(2)
                     ! then apply rules concerning inner and outer radius
                     if(dsqrt(kp(1)**2+kp(2)**2)<k_mask_bounds(3)) kmask(kpt)= 1
                     if(dsqrt(kp(1)**2+kp(2)**2)<k_mask_bounds(4)) kmask(kpt)= 0
                  end do ! j
               end do ! i
            elseif (kmask_mode==2) then ! box mode
               do i=-1,1,1
                  do j=-1,1,1
                     kp(1:3) = bzkp(1:3,kpt)+i*recbv(1:3,1)+j*recbv(1:3,2)
                     if (kp(1)<k_mask_bounds(1)) kmask(kpt)= 0
                     if (kp(1)>k_mask_bounds(2)) kmask(kpt)= 0
                     if (kp(2)<k_mask_bounds(3)) kmask(kpt)= 0
                     if (kp(2)>k_mask_bounds(4)) kmask(kpt)= 0
                  end do ! j
               end do ! i
            elseif (kmask_mode==3) then ! read kmask from file
               read(8888, *) kmask(kpt)
            end if ! kmask_mode
            ! count number of kpts in reduced part
            if(kmask(kpt)>0) k_end = k_end+1
         end do ! kpt loop
        
         ! close kmask file
         if (kmask_mode==3)then ! read kmask from file
            close(8888)
         end if
         
         ! fill rhoq_kmask (on reduced set of kpts)
         allocate(rhoq_kmask(1:5,k_end))
         do kpt=1, nofks
            if(kmask(kpt)>0) then
              rhoq_kmask(1:3,k_start) = bzkp(1:3,kpt)
              rhoq_kmask(4,k_start) = dfloat(kpt)
              rhoq_kmask(5,k_start) = dfloat(kmask(kpt))
              k_start = k_start+1
            end if
         end do
        
         write(*,*) 'found ', k_end,  'kpoints'
        
         open(8888, file='rhoq_kmask.test', form='formatted')
         do kpt=1, k_end
            write(8888, '(5F14.7)') rhoq_kmask(1:3,kpt), rhoq_kmask(4,kpt), rhoq_kmask(5,kpt)
         end do
         close(8888)
        
      end if !(myrank==master)

#if defined(CPP_MPI)
      ! communicate kmask stuff from master to all others
      call MPI_BCAST(k_end, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      if (myrank .ne. master) then
        allocate(kmask(nofks))
        allocate(rhoq_kmask(1:5,k_end))
      end if
      call MPI_BCAST(kmask, nofks, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(rhoq_kmask, 5*k_end, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
#endif

   end subroutine rhoq_find_kmask

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rhoq_saveG(nscoef, rhoq_kmask, kpt, nofks, k_end, kp, i, j, mu, imin, iatomimp, lmmaxd, G)

#ifdef CPP_HYBRID
      use omp_lib
#endif

      implicit none

      integer, intent(in) :: i, j, mu, imin, lmmaxd, nscoef, nofks, k_end, kpt
      integer, intent(in) :: iatomimp(nscoef)
      double precision, intent(in) :: rhoq_kmask(5, k_end), kp(3)
      double complex, intent(in) :: G(lmmaxd, lmmaxd)
      ! local
      integer :: ix, jx, lm1, irec 

#ifdef CPP_HYBRID
     !$omp critical
#endif
     irec = (nscoef*2)*(int(rhoq_kmask(4,kpt))-1)
     if( ((i==mu) .and. any(j==iatomimp(1:nscoef))) ) then
       ix=0
       jx=0
       lm1 = 1
       do while (ix==0 .and. lm1<=nscoef)
         if(iatomimp(lm1)==j) ix = j - imin + 1
         lm1 = lm1 + 1
       end do
       irec = irec + nscoef + ix
       write(9889,rec=irec) KP(1:2), G(1:LMMAXD,1:LMMAXD)! * 
!                                     rhoq_kmask(5,kpt)
     end if
     
     irec = (nscoef*2)*(int(rhoq_kmask(4,kpt))-1)
     if( ((j==mu) .and. any(i==iatomimp(1:nscoef))) ) then
       ix=0
       jx=0
       lm1 = 1
       do while (jx==0 .and. lm1<=nscoef+1)
         if(iatomimp(lm1)==i) jx = i - imin + 1
         lm1 = lm1 + 1
       end do
       irec = irec + jx
       write(9889,rec=irec) kp(1:2), g(1:LMMAXD,1:LMMAXD)! * !*ETAIKR(ISYM,NS) not this factor because it is dealt with explicitly in rhoq module
!                                     rhoq_kmask(5,kpt)
     end if ! i==mu ...
#ifdef CPP_HYBRID
     !$omp end critical
#endif

   end subroutine rhoq_saveG

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rhoq_write_tau0(nofks, nshell, nsh1, nsh2, nsymat, nscoef, mu, iatomimp, kmask, lmmaxd, bzkp, imin)

   use mod_mympi, only: myrank, master

   implicit none

   double complex, parameter :: CZERO = (0.0d0, 0.0d0)

   integer, intent(in) :: nofks, nshell, nsymat, nscoef, mu, lmmaxd, imin
   integer, intent(in) :: nsh1(nshell), nsh2(nshell), iatomimp(nscoef)
   double precision, intent(in) :: bzkp(3, nofks)
   integer, allocatable, intent(inout) :: kmask(:)
   !local
   integer :: kpt, ns, i, j, isym, irec, ix, jx, lm1
   double precision :: kp(3)
   double complex :: G(lmmaxd, lmmaxd)

   if(myrank==master) then
      write(*,*)                      ! status bar
      write(*,*) 'rhoq: write-out loop'
      write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
      write(*,FMT=190) !beginning of statusbar
      
      ! write out fort.998899, fort.998888
      ! ======================================================================
      do kpt=1, nofks
        DO NS = 1,NSHELL
          I = NSH1(NS)
          J = NSH2(NS)
          DO ISYM = 1,NSYMAT
            irec = (nscoef*2)*(kpt-1)
            if( ((i==mu) .and. any(j==iatomimp(1:nscoef))) ) then
              ix=0
              jx=0
              lm1 = 1
              do while (ix==0 .and. lm1<=nscoef)
                if(iatomimp(lm1)==j) ix = j - imin + 1
                lm1 = lm1 + 1
              end do
              irec = irec + nscoef + ix 
              if(kmask(kpt)>0) then
                 read(9889,rec=irec) KP(1:2), G(1:LMMAXD,1:LMMAXD)
              else
                 KP(1:3) = bzkp(1:3,kpt)
                 G(1:LMMAXD,1:LMMAXD) = CZERO
              end if
              write(998899,'(10000ES15.7)') KP(1:2), G(1:LMMAXD,1:LMMAXD) * dfloat(kmask(kpt))
            end if
            irec = (nscoef*2)*(kpt-1)
            if( ((j==mu) .and. any(i==iatomimp(1:nscoef))) ) then
              ix=0
              jx=0
              lm1 = 1
              do while (jx==0 .and. lm1<=nscoef+1)
                if(iatomimp(lm1)==i) jx = i - imin + 1
                lm1 = lm1 + 1
              end do
              irec = irec + jx
              if(kmask(kpt)>0) then
                 read(9889,rec=irec) KP(1:2), G(1:LMMAXD,1:LMMAXD)
              else
                 KP(1:3) = bzkp(1:3,kpt)
                 G(1:LMMAXD,1:LMMAXD) = CZERO
              end if
              write(998888,'(10000ES15.7)') KP(1:2), G(1:LMMAXD,1:LMMAXD) * dfloat(kmask(kpt))
            end if ! iii==mu ...
          END DO ! ISYM = 1,NSYMAT
        END DO !NS = 1,NSHELL
        
        if(nofks>=50) then
          if(mod(KPT-0,nofks/50)==0) write(6,FMT=200)
        else
          write(6,FMT=200)
        end if
        
      end do !kpt=1,nofks
      ! ======================================================================
      write(6,*)                      ! status bar
   end if !myrank==master
   deallocate(kmask)

190   FORMAT('                 |'$)      ! status bar
200   FORMAT('|'$)                       ! status bar

   end subroutine rhoq_write_tau0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rhoq_save_rmesh(natyp, irmd, ipand, irmin, irws, ipan, rmesh, ntcell, ircut, r_log, npan_log, npan_eq)

   implicit none

   integer, intent(in) :: natyp, irmd, ipand, npan_log, npan_eq
   integer, intent(in) :: irmin(natyp), irws(natyp), ipan(natyp), ntcell(natyp), ircut(0:ipand,natyp)
   double precision, intent(in) :: r_log
   double precision, intent(in) :: rmesh(irmd, natyp)
   !local
   integer :: i1

   !read mu_0
   open(9999, file='mu0', form='formatted')
   read(9999,*) i1
   close(9999)
   ! write out corresponding mesh information
   open(9999, file='rmesh_mu0.txt', form='formatted')
   write(9999,'(A)') '# mu_0, IRMIN, IRWS, IPAN'
   write(9999,'(4I9)') i1,IRMIN(i1), IRWS(i1), IPAN(i1)
   write(9999,'(A)') '# Rmesh(1:IRWS)'
   write(9999,'(1000E22.15)') Rmesh(1:IRWS(i1),i1)
   write(9999,'(A)') '# NTCELL(1:NATYP)'
   write(9999,'(1000I9)') NTCELL(1:NATYP)
   write(9999,'(A)') '# IRCUT(1:IPAN)'
   write(9999,'(1000I9)') IRCUT(1:IPAN(i1),i1)
   write(9999,'(A)') '# R_LOG, NPAN_LOG, NPAN_EQ'
   write(9999,'(E22.15,2I9)') R_LOG, NPAN_LOG, NPAN_EQ
   close(9999)            


   end subroutine rhoq_save_rmesh 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rhoq_save_refpot(ielast, i1, nref, natyp, refpot, wlength, lmmaxd, ie, trefll)

   implicit none

   integer, intent(in) :: ielast, i1, nref, natyp, wlength, lmmaxd, ie
   integer, intent(in) :: refpot(natyp)
   double complex, intent(in) :: trefll(lmmaxd,lmmaxd,natyp)
   !local
   integer :: irec

   if(i1==1) then
      open(99991, file='refinfo')
      write(99991,'(2I9)') NREF, NATYP
      write(99991, '(1000I9)') refpot(1:NATYP)
      close(99991)
      open(99992, file='tref', access='direct', recl=WLENGTH*4*lmmaxd**2, form='unformatted')
   end if
   irec = ie + ielast*(i1-1)
   write(99992, rec=irec) trefll(:,:,i1)

   end subroutine rhoq_save_refpot


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_rhoqtools
