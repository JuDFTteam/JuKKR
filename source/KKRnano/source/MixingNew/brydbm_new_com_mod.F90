!> Contains the classic Broyden and generalised Anderson mixing
!> routines + helpers.
!>
!> Modified to remove all atom dependencies.

      MODULE BRYDBM_new_com_mod
      implicit none
      private
      public :: brydbm_new_com, brysh1_new, brysh2_new, brysh3_new

      contains

!*********************************************************************
!> broyden mixing.
!>
!> \verbatim
!>    imix :
!>      3      broyden's            f i r s t  m e t h o d
!>      4      broyden's          s e c o n d  m e t h o d
!>      5      anderson's     g e n e r a l i z e d   m e t h o d
!>
!>    implemented here according to notes of s.b.
!>    broyden's iteration scheme following the papers of :
!>    srivastava, j. phys. , 17 (1984) , pp l317
!>    c.g. broyden in math.comput., 19 , pp 577, 1965
!>    c.g. broyden in ibid, 21 ,pp 368 ,1967
!>    the method has been generalized to include a metric.the
!>    definition of the necessary inner products are similar to the
!>    discription given in the notes of m.weinert.the algorithm
!>    discribed in the paper srivastava  has been simplified
!>    ( see notes of s.b.)
!>    the files ui,vi are stored on high speed ssd memory.
!>    broyden's update treats charge and spin on the same footing
!>                 s. bluegel , kfa , may 1987
!>    the anderson method (d.g. anderson, j. acm 12, 547 (1964)) has
!>    been generalized and reformulated as an improvement of broyden's
!>    second method. successive linesearch is replaced by successive
!>    search on hyperplanes. ( see notes of s.b. )
!>                 s. bluegel , issp , july 1989
!>
!>    modified for non spherical potential
!>                 b. drittler , aug. 1988
!>
!>    parallelized over the number of atoms
!>                 a. thiess , jun. 2008
!>
!>    made independent of atoms e.r.
!>
!> \endverbatim
!*********************************************************************

!*********************************************************************
      subroutine brydbm_new_com(visp,v,vins,lmpot,r,drdi,alpha,irc,irmin,nspin,imix,iter,ui2,vi2,wit,sm1s,fm1s,mylrank,communicator, itdbryd, irmd, irnsd, nspind)
      implicit none
      include 'mpif.h'
      external :: MPI_Allreduce

      integer, intent(in) :: itdbryd !< array dimension for broyden mixing history
      integer, intent(in) :: irmd !< number of radial mesh points - total
      integer, intent(in) :: irnsd !< number of radial mesh points of non-spherical part
      integer, intent(in) :: nspind !< number of spin directions
      integer, intent(in) :: imix
      integer, intent(in) :: lmpot, nspin

!      parameter (lmpotd= (lpotd+1)**2)
!      parameter (irmind=irmd-irnsd)
!      parameter (ntird=(irmd+(irnsd+1)*(lmpotd-1))*nspind)
      double precision, intent(in) :: drdi(irmd), r(irmd)
      double precision, intent(inout) :: v(irmd,lmpot,*)
!     double precision vins(irmind:irmd,lmpot,*)
      double precision, intent(in) :: vins(irmd-irnsd:irmd,lmpot,*)
      double precision, intent(in) :: visp(irmd,*)
      integer, intent(in) :: irc, irmin
      integer, intent(in) :: mylrank
      integer, intent(in) :: communicator
      double precision, intent(in) :: alpha
      double precision, intent(inout) :: sm1s(:), fm1s(:) !> dim((irmd+(irnsd+1)*(lmpot-1))*nspind)
      double precision, intent(inout) :: ui2(:,:), vi2(:,:) !> dim((irmd+(irnsd+1)*(lmpot-1))*nspind,2:itdbryd)
      double precision, intent(inout) :: wit(2:itdbryd)
      
!     ..
!     .. local scalars ..
      double precision :: rmixiv,vmdeno,vmnorm,volinv,smnorm_local,smnorm_global,ddot_local,ddot_global,cmm_local,cmm_global 
      integer :: ij, imap, ir, isp, it, lm, mit, iter

      integer :: ntird
      integer :: ierr, ist
!     ..
!     .. external functions ..
      double precision, external :: ddot
!     ..
      double precision, parameter :: zero=0.d0, one=1.d0
      integer, parameter :: itdthd=40
      integer, parameter :: MEM_ALIGN=1
!     ..
!     .. local arrays ..
!     .. these arrays are automatic (dynamically allocated)
!        fortran arrays
      double precision :: am(2:itdbryd), bm(2:itdbryd), am_local(2:itdbryd), bm_local(2:itdbryd)

!     the following arrays have dimension ntird
!     parameter (ntird=(irmd+(irnsd+1)*(lmpotd-1))*nspind)

      double precision, allocatable, dimension(:) :: fm, fm1, g, sm, sm1, work, vi3, ui3 !> dim((irmd+(irnsd+1)*(lmpot-1))*nspind)

      ntird = (irmd + (irnsd + 1)*(lmpot - 1))*nspind
      if (size(sm1s) < ntird) stop __LINE__
      if (size(fm1s) < ntird) stop __LINE__
      ! introduce memory alingment here
      ntird = ((ntird - 1)/MEM_ALIGN + 1)*MEM_ALIGN

      if (mylrank == 0) write(*,*) 'iteration: ',iter

      if (itdbryd > itdthd .or. itdthd > 200) stop 'itdbryd'
      
      allocate(fm(ntird), fm1(ntird), g(ntird), sm(ntird), sm1(ntird), work(ntird), vi3(ntird), ui3(ntird), stat=ist)
      if (ist /= 0) stop "Broyden new com: allocation failed!"
      
      if (imix <= 2 .or. imix > 5) stop 'imixd'

      mit = iter
      do while (mit > itdbryd)
        mit = mit - itdbryd
      enddo ! while

!     for the moment the broyden mixing is started from scratch if itdbryd 
!     is exceeded - in principle also possible to make the transition
!     continous                                         2008/09/16     


      if (mylrank == 0) then
        if (imix == 3) write (6,fmt='('' broyden"s 1st method used '')')
        if (imix == 4) write (6,fmt='('' broyden"s 2nd method used '')')
        if (imix == 5) write (6,fmt='('' gen. anderson method used '')')
        write(6,*) 'ITERATION: ',mit
      endif

!
      rmixiv = one/alpha
!
!---->  the following block is activated only one iteration before
!        broyden iteration scheme is used
!---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
!                   metric  g := r*r*drdi
!---->  map data of all muffin-tin spheres into one single vector
!
!
      imap = brysh3_new(sm1,visp,vins,irmin,irc,nspin,lmpot,irmd,irnsd)

      imap = brysh1_new(fm1,v,irmin,irc,nspin,lmpot,irmd)
      if (imap > ntird) stop 'nirdbry'

      do ij = 1,imap
        fm1(ij) = rmixiv*(fm1(ij) - sm1(ij))
      enddo ! ij

      ! determine diagonal metric g(:)
      ij = 0
      do isp = 1,nspin

        volinv = 3.d0/(r(irc)**3)
        do ir = 1,irc
          ij = ij + 1
          g(ij) = volinv*r(ir)*r(ir)*drdi(ir)
        enddo ! ir

        do lm = 2,lmpot
          do ir = irmin,irc
            ij = ij + 1
            g(ij) = volinv*r(ir)*r(ir)*drdi(ir)
          enddo ! ir
        enddo ! lm

      enddo ! isp

!=========================================================================
!=========================================================================
!=========================================================================
!=====  for mit gt 1 activ  ==============================================
!=========================================================================
!=========================================================================
!=========================================================================

      if (mit > 1) then

        do ij = 1,imap
          sm1(ij) = sm1s(ij)
          fm1(ij) = fm1s(ij)
        enddo ! ij
!
!----> map rho(m) of all mt-spheres into one single vector
!
        imap = brysh3_new(sm,visp,vins,irmin,irc,nspin,lmpot,irmd,irnsd)

!
!----> map f[m] = f(m) - rho(m) = f(rho(m)) - rho(m) of all mt-spheres
!      into one single vector
!
        imap = brysh1_new(fm,v,irmin,irc,nspin,lmpot,irmd)

        do ij = 1,imap
          fm(ij) = rmixiv* (fm(ij)-sm(ij))
        enddo ! ij
!
!----> calculate  sm = rho(m) - rho(m-1)
!----> calculate dfm = f[m] - f[m-1]
!
        do ij = 1,imap
          sm1(ij) = sm(ij) - sm1(ij)
          fm1(ij) = fm(ij) - fm1(ij)
        enddo ! ij

!
!----> loop to generate u[m] = u(ij,mit)
!
        do ij = 1,imap
          ui3(ij) = alpha*fm1(ij) + sm1(ij)
        enddo ! ij


        do it = 2,mit - 1
          do ij = 1,imap
            work(ij) = vi2(ij,it)
          enddo ! ij

          am_local(it) = ddot(imap,fm1,1,work,1)
        enddo

        call MPI_Allreduce(am_local,am,(itdbryd-1),MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

        do it = 2,mit - 1
          do ij = 1,imap
            work(ij) = ui2(ij,it)
          enddo ! ij
          call daxpy(imap,-am(it),work,1,ui3,1)
        enddo ! it

!
!----> print amj = the importance of the history of ui
!
!        write (ipf,fmt='(5x,'' amj , ---> j=2,'',i3,/,(9x,1p,7d10.2))') mit - 1, (am(it),it=2,mit-1)
!
!

!=========================================================================
!=========================================================================
!=========== broyd 1 =====================================================
!=========================================================================
!=========================================================================

        if (imix == 3) then
!-------->     b r o y d e n ' s   f i r s t   m e t h o d
!
!----> calculate dsmnorm
!
          smnorm_local = zero
          do ij = 1,imap
            smnorm_local = smnorm_local + sm1(ij)*g(ij)*sm1(ij)
          enddo ! ij

        call MPI_Allreduce(smnorm_local,smnorm_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)


!
!----> convolute dsm with the metric g
!
          do ij = 1,imap
            sm1(ij) = g(ij)*sm1(ij)
          enddo ! ij
!
!----> loop to generate v[m] = v(ij,mit)
!
          do ij = 1,imap
            vi3(ij) = alpha*sm1(ij)
          enddo ! ij


          do it = 2,mit - 1
            do ij = 1,imap
              work(ij) = ui2(ij,it)
            enddo ! ij
            bm_local(it) = ddot(imap,sm1,1,work,1)
          enddo

        call MPI_Allreduce(bm_local,bm,(itdbryd-1),MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

          do it = 2,mit - 1
            do ij = 1,imap
              work(ij) = vi2(ij,it)
            enddo ! ij
            call daxpy(imap,-bm(it),work,1,vi3,1)
          enddo ! it


!
!----> complete the evaluation of v[m]
!
        ddot_local = ddot(imap,sm1,1,ui3,1)

        call MPI_Allreduce(ddot_local,ddot_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

        vmdeno = ddot_global - smnorm_global

          if (abs(vmdeno) < 1.d-70) stop 'bry0sn'

          call dscal(imap,one/vmdeno,vi3,1)
!
!----> print bmj = the importance of the history of vi
!
!          write (ipf,fmt='(5x,'' bmj , ---> j=2,'',i3,/,(9x,1p,7d10.2))') mit - 1, (bm(it),it=2,mit-1)
!


!=========================================================================
!=========================================================================
!=========== broyd 2 =====================================================
!=========================================================================
!=========================================================================

        else if (imix == 4) then

!-------->     b r o y d e n ' s   s e c o n d    m e t h o d
!----> calculate v[m] ; convoluted with the metric g

          do ij = 1,imap
            vi3(ij) = g(ij)*fm1(ij)
          enddo ! ij

!----> calculate #vm# and normalize v[m]

          ddot_local = ddot(imap,vi3,1,fm1,1)

          call MPI_Allreduce(ddot_local,ddot_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

          vmnorm = ddot_global

          call dscal(imap,one/vmnorm,vi3,1)

!=========================================================================
!=========================================================================
!=========== anderson ====================================================
!=========================================================================
!=========================================================================

        else if (imix == 5) then

!-------->     g e n e r a l i z e d   a n d e r s o n   m e t h o d
!
!----> calculate v[m] ; convoluted with the metric g
!
          do ij = 1,imap
            vi3(ij) = g(ij)*fm1(ij)
          enddo ! ij
  
          do it = 2,mit - 1
            call daxpy(imap,-am(it)*wit(it),vi2,1,vi3,1)
          enddo ! it

!----> complete the evaluation of v[m]

          ddot_local = ddot(imap,fm1,1,vi3,1)

          call MPI_Allreduce(ddot_local,ddot_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

          vmdeno = ddot_global

          if (abs(vmdeno) < 1.d-70) stop 'bry1sn'

          call dscal(imap,one/vmdeno,vi3,1)

!----> save wit(mit) for next iteration

          wit(mit) = vmdeno

        endif

!=========================================================================
!=========================================================================
!============ endmixing, now output =====================================
!=========================================================================
!=========================================================================

!
!----> write u3(ij) and v3(ij) on disk
!
        do ij = 1,imap
          ui2(ij,mit) = ui3(ij)
          vi2(ij,mit) = vi3(ij)
        enddo ! ij

!
!----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
!
        do ij = 1,imap
          fm1(ij) = fm(ij)
          sm1(ij) = sm(ij)
        enddo ! ij
!
!----> calculate cmm
!
        cmm_local = ddot(imap,fm,1,vi3,1)

        call MPI_Allreduce(cmm_local,cmm_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

!           write (ipf,fmt='(5x,'' cmm = '',1p,d12.4)') cmm
!
!----> update rho(m+1)
!

        call daxpy(imap,one-cmm_global,ui3,1,sm,1)

!
!----> map solution back into each mt-sphere
!
        imap = brysh2_new(sm,v,irmin,irc,nspin,lmpot,irmd)

      endif

!=========================================================================
!=========================================================================
!=========================================================================
!=====  for mit gt 1 activ  ==============================================
!=========================================================================
!=========================================================================
!=========================================================================

!      mit = mit + 1

    do ij = 1,imap
      sm1s(ij) = sm1(ij)
      fm1s(ij) = fm1(ij)
    enddo ! ij

    deallocate(fm, fm1, g, sm, sm1, work, vi3, ui3, stat=ist) ! ignore status
  endsubroutine ! brydbm_new_com

! ************************************************************************
  integer function brysh1_new(y,x,irmin,irc,nspin,lmpot,irmd) result(imap)
!*********************************************************************
!     shifts the density or potential of all mt-cell into one single
!     vector and projects out the coulomb part only.
!                                    s. bluegel , kfa , 1987
!     modified for parallelization
!                                    a. thiess, jun 2008
!     make it independent of atoms
!                                    e. rabel
! ------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: irmd, lmpot, nspin, irc, irmin
    double precision, intent(in) :: x(irmd,lmpot,*)
    double precision, intent(out) :: y(*)

!     .. local scalars ..
    integer :: ir, is, lm

    imap = 0
    do is = 1,nspin

      do ir = 1,irc
        imap = imap + 1
        y(imap) = x(ir,1,is)
      enddo ! ir

      do lm = 2,lmpot
        do ir = irmin,irc
          imap = imap + 1
          y(imap) = x(ir,lm,is)
        enddo ! ir
      enddo ! lm

    enddo ! is

  endfunction ! brysh1_new

! ************************************************************************
  integer function brysh2_new(y,x,irmin,irc,nspin,lmpot,irmd) result(imap)
!*********************************************************************
!     maps the density or potential back from one single vector into
!     the proper bins of each single mt-cell . the magnetization
!     density is also added in.
!                                    s. bluegel , kfa , 1987
! ------------------------------------------------------------------------
    implicit none

    double precision, intent(in) :: y(*)
    double precision, intent(out) :: x(irmd,lmpot,*)
    integer, intent(in) :: irmd, lmpot, nspin, irc, irmin
    
!     .. local scalars ..
    integer :: ir, is, lm

    imap = 0
    do is = 1,nspin

      do ir = 1,irc
        imap = imap + 1
        x(ir,1,is) = y(imap)
      enddo ! ir

      do lm = 2,lmpot
        do ir = irmin,irc
          imap = imap + 1
          x(ir,lm,is) = y(imap)
        enddo ! ir
      enddo ! lm
          
    enddo ! is

  endfunction ! brysh2_new

!>    throws stuff from x and z into y.
!>    z after x
!>    y = (x(:), z(:))  - why so complicated

!>    @param y single vector containing all lm,spin components
!>    @param x spherical potential
!>    @param z non-spherical potential
! ************************************************************************
  integer function brysh3_new(y,x,z,irmin,irc,nspin,lmpot,irmd,irnsd) result(imap)
!*********************************************************************
!     shifts the density or potential of all mt-cell into one single
!     vector and projects out the coulomb part only.
!
!                                    s. bluegel , kfa , 1987
!     modified for parallelization
!                                    a. thiess , jun 2008
! ------------------------------------------------------------------------
    implicit none
    double precision, intent(out) :: y(*)
    double precision, intent(in) :: x(irmd,*), z(irmd-irnsd:irmd,lmpot,*)
    integer, intent(in) :: irmd, irnsd, lmpot, nspin, irc, irmin

!     .. local scalars ..
    integer :: ir, is, lm
      
    imap = 0
    do is = 1,nspin

      do ir = 1,irc
        imap = imap + 1
        y(imap) = x(ir,is)
      enddo ! ir

      do lm = 2,lmpot
        do ir = irmin,irc
          imap = imap + 1
          y(imap) = z(ir,lm,is)
        enddo ! ir
      enddo ! lm

    enddo ! is

  endfunction ! brysh3_new

endmodule ! brydbm_new_com_mod
