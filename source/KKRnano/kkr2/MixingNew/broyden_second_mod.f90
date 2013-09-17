!> Extracted Broyden's 2nd method from classic code.

module broyden_second_mod

contains

    !*********************************************************************
    !> Broyden mixing.
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

  subroutine calc_metric(g, lmpot,r,drdi,irc,irmin,nspin,imap)
    implicit none

    double precision, intent(out) :: g(imap)
    integer :: lmpot, imap, irc, irmin, nspin
    double precision :: r(:), drdi(:)

    integer ij, ir, lm
    integer isp
    double precision volinv

    ij = 0
    do 60 isp = 1,nspin

      volinv = 3.0d0/(r(irc)**3)
      do 20 ir = 1,irc
        ij = ij + 1
        g(ij) = volinv*r(ir)*r(ir)*drdi(ir)
20    continue
      !
      if (lmpot.gt.1) then

        do 40 lm = 2,lmpot
          do 30 ir = irmin,irc
            ij = ij + 1
            g(ij) = volinv*r(ir)*r(ir)*drdi(ir)
30        continue
40      continue
      end if

50  continue

60 continue
 end subroutine

 !*********************************************************************
 ! PARAMETERS STILL UNCLEAR!!!
 ! arrays destroyed on output, result in sm
 ! sm = V_in on input, new V_in on output???
 ! sm1 = V_in
 ! fm1 = V_in + alpha * V_out
 subroutine broyden_second(sm, alpha, sm1, fm1, ui2,vi2,sm1s,fm1s, g, &
 communicator, itdbryd, imap, mit)

   implicit none

   include 'mpif.h'

   !     array dimension for broyden mixing history
   integer itdbryd

   double precision, parameter :: ONE = 1.0d0

   !     ..
   !     .. Local Scalars ..
   double precision rmixiv,vmnorm, &
   ddot_local, &
   ddot_global,cmm_local,cmm_global
   integer ij,imap,it,mit

   integer ierr
   !     ..
   !     .. External Functions ..
   double precision ddot
   external ddot
   !     ..
   !     .. Local Arrays ..
   !     .. these arrays are automatic (dynamically allocated)
   !        Fortran arrays
   double precision am(2:itdbryd), am_local(2:itdbryd)

   !     the following arrays have dimension NTIRD
   !     PARAMETER (NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND)
   double precision fm  (imap)
   double precision fm1 (imap)
   double precision g   (imap)
   double precision sm  (imap)
   double precision sm1 (imap)
   double precision sm1s(imap)
   double precision fm1s(imap)
   double precision work(imap)
   double precision vi3 (imap)


   double precision ui3(imap)

   double precision ui2(imap, 2:itdbryd)
   double precision vi2(imap, 2:itdbryd)

   !     .. Scalar Arguments ..
   double precision alpha

   integer communicator
   !
   external mpi_allreduce

   rmixiv = one/alpha
   !
   !---->  the following block is activated only one iteration before
   !        broyden iteration scheme is used
   !---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
   !                   metric  g := r*r*drdi

   ! from simple mixed V_out -> reconstruct unmixed V_out ! OMG!
   do 10 ij = 1,imap
     fm1(ij) = rmixiv* (fm1(ij)-sm1(ij))
10 continue
    !

   !=====  For MIT GT 1 activ  ==============================================

   if (mit.gt.1) then

     do ij = 1,imap
       sm1(ij)=sm1s(ij)
       fm1(ij)=fm1s(ij)
     enddo

     ! fm = 1/alpha * (V_out[V_in] - V_in)
     do 70 ij = 1,imap
       fm(ij) = rmixiv* (fm(ij)-sm(ij))
70   continue
     !
     !----> calculate  sm = rho(m) - rho(m-1)
     !----> calculate dfm = f[m] - f[m-1]
     !
     do 80 ij = 1,imap
       sm1(ij) = sm(ij) - sm1(ij)
       fm1(ij) = fm(ij) - fm1(ij)
80   continue

     ! sm1 = \Delta V_in
     ! fm1 = \Delta V_out

     !
     !----> loop to generate u[m] = u(ij,mit)
     !
     ! = alpha*\Delta V_out + \Delta V_in
     do 90 ij = 1,imap
       ui3(ij) = alpha*fm1(ij) + sm1(ij)
90   continue


     do it = 2,mit - 1
       do ij = 1,imap
         work(ij)=vi2(ij,it)
       enddo

       am_local(it) = ddot(imap,fm1,1,work,1)
     enddo

     call mpi_allreduce(am_local,am,(itdbryd-1), &
     mpi_double_precision,mpi_sum,communicator,ierr)

     ! active only from 3rd iteration on
     ! ui3 = \sum_i^j -am_i (alpha*\Delta V_out + \Delta V_in)_i
     ! with
     ! am_i = (\Delta V_out)_i^T . G . (\Delta V_out)_j ???
     do it = 2,mit - 1
       do ij = 1,imap
         work(ij)=ui2(ij,it)
       enddo
       call daxpy(imap,-am(it),work,1,ui3,1)
     enddo

     !-------->     b r o y d e n ' s   s e c o n d    m e t h o d
     !----> calculate v[m] ; convoluted with the metric g

     ! vi3 = G . \Delta V_out
     do 150 ij = 1,imap
       vi3(ij) = g(ij)*fm1(ij)
150  continue

     !----> calculate #vm# and normalize v[m]

     ! = (\Delta V_out)^T . G . (\Delta V_out) ! using diagonal metric matrix G
     ddot_local = ddot(imap,vi3,1,fm1,1)

     call mpi_allreduce(ddot_local,ddot_global,1, &
     mpi_double_precision,mpi_sum,communicator,ierr)

     vmnorm = ddot_global

     ! normalize (\Delta V_out)
     call dscal(imap,one/vmnorm,vi3,1)

     !============ END MIXING, NOW OUTPUT =====================================

     !
     !----> write u3(ij) and v3(ij) on disk
     !
     do ij = 1,imap
       ui2(ij,mit)=ui3(ij)
       vi2(ij,mit)=vi3(ij)
     enddo

     !
     !----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
     !
     do 180 ij = 1,imap
       fm1(ij) = fm(ij)
       sm1(ij) = sm(ij)
180  continue
     !
     !----> calculate cmm
     !

     ! cmm = 1/alpha * (V_out[V_in] - V_in] . (\Delta V_out / ||\Delta V_out||)
     cmm_local = ddot(imap,fm,1,vi3,1)

     call mpi_allreduce(cmm_local,cmm_global,1, &
     mpi_double_precision,mpi_sum,communicator,ierr)

     !----> update rho(m+1)
     !
     ! V_in_new = (1 - cmm) * (alpha*\Delta V_out + \Delta V_in) + V_in
     call daxpy(imap,one-cmm_global,ui3,1,sm,1)

   end if
   !=====  For MIT GT 1 activ  ==============================================

   !      MIT = MIT + 1

   do ij = 1,imap
     sm1s(ij)=sm1(ij)
     fm1s(ij)=fm1(ij)
   enddo

   return

 end subroutine

end module
