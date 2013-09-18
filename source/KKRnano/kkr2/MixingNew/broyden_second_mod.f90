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

  subroutine calc_metric(g_metric, lmpot,r,drdi,irc,irmin,nspin,imap)
    implicit none

    double precision, intent(out) :: g_metric(imap)
    integer :: lmpot, imap, irc, irmin, nspin
    double precision :: r(:), drdi(:)

    integer ij, ir, lm
    integer isp
    double precision volinv

    ij = 0
    do isp = 1,nspin

      volinv = 3.0d0/(r(irc)**3)
      do ir = 1,irc
        ij = ij + 1
        g_metric(ij) = volinv*r(ir)*r(ir)*drdi(ir)
      end do
      !
      if (lmpot.gt.1) then

        do lm = 2,lmpot
          do ir = irmin,irc
            ij = ij + 1
            g_metric(ij) = volinv*r(ir)*r(ir)*drdi(ir)
          end do
        end do
      end if

   end do
 end subroutine

 !*********************************************************************
 ! PARAMETERS STILL UNCLEAR!!!
 ! arrays destroyed on output, result in sm
 ! sm = input: V_in from current iteration, new V_in on output
 ! fm = input: simple mixed pot. from current iteration
 !      output: changed on output - output needed for next iteration
 ! sm1 = it>1: work array
 ! fm1 = it>1: work array
 ! g_metric ... diagonal of metric matrix
subroutine broyden_second(sm, fm, sm1, fm1, ui2,vi2, g_metric, alpha, &
                          communicator, itdbryd, imap, mit)

  implicit none

  include 'mpif.h'

  ! Parameters
  double precision, parameter :: ONE =1.0d0

  ! Arguments
  double precision, intent(inout), dimension(imap) :: sm
  double precision, intent(inout), dimension(imap) :: fm
  double precision, intent(inout), dimension(imap) :: sm1
  double precision, intent(inout), dimension(imap) :: fm1
  double precision, intent(inout), dimension(imap, 2:itdbryd) :: ui2
  double precision, intent(inout), dimension(imap, 2:itdbryd) :: vi2

  double precision, intent(in) :: alpha
  double precision, intent(in), dimension(imap) :: g_metric
  integer, intent(in) :: communicator
  integer, intent(in) :: itdbryd
  integer, intent(in) :: imap
  integer, intent(in) :: mit

  ! Local variables of broyden_second
  double precision :: rmixiv
  double precision :: vmnorm
  double precision :: ddot_local
  double precision :: ddot_global
  double precision :: cmm_local
  double precision :: cmm_global
  integer :: ij
  integer :: it
  integer :: ierr
  double precision, dimension(2:itdbryd) :: am
  double precision, dimension(2:itdbryd) :: am_local
  double precision, dimension(imap) :: work
  double precision, dimension(imap) :: vi3
  double precision, dimension(imap) :: ui3

   !    .. External Functions ..
   double precision, external :: ddot
   external MPI_Allreduce

   rmixiv = one/alpha
   !
   !---->  the following block is activated only one iteration before
   !        broyden iteration scheme is used
   !---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
   !                   metric  g := r*r*drdi

   ! from simple mixed V_out -> reconstruct unmixed V_out ! OMG!
   ! fm = 1/alpha * (V_out[V_in] - V_in)
   do ij = 1,imap
     fm(ij) = rmixiv* (fm(ij)-sm(ij))
   end do

   !=====  For MIT GT 1 activ  ==============================================

   if (mit.gt.1) then

     !----> calculate  sm = rho(m) - rho(m-1)
     !----> calculate dfm = f[m] - f[m-1]
     !
     do ij = 1,imap
       sm1(ij) = sm(ij) - sm1(ij)
       fm1(ij) = fm(ij) - fm1(ij)
     end do

     ! sm1 = \Delta V_in
     ! fm1 = \Delta V_out

     !----> loop to generate u[m] = u(ij,mit)
     !
     ! See Srivastava (16).
     ! = alpha*\Delta V_out + \Delta V_in
     do ij = 1,imap
       ui3(ij) = alpha*fm1(ij) + sm1(ij)
     end do

     do it = 2,mit - 1
       do ij = 1,imap
         work(ij)=vi2(ij,it)
       enddo
       am_local(it) = ddot(imap,fm1,1,work,1)
     enddo

     call MPI_Allreduce(am_local,am,(itdbryd-1), &
     MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

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
     do ij = 1,imap
       vi3(ij) = g_metric(ij)*fm1(ij)
     end do

     !----> calculate #vm# and normalize v[m]

     ! = (\Delta V_out)^T . G . (\Delta V_out) ! using diagonal metric matrix G
     ddot_local = ddot(imap,vi3,1,fm1,1)

     call MPI_Allreduce(ddot_local,ddot_global,1, &
     MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

     vmnorm = ddot_global

     ! normalize (\Delta V_out)
     ! v^T(i) in Srivastava paper, equ. (16)
     call dscal(imap,one/vmnorm,vi3,1)

     !============ END MIXING, NOW OUTPUT =====================================

     !----> store u3(ij) and v3(ij) for following iterations
     do ij = 1,imap
       ui2(ij,mit)=ui3(ij)
       vi2(ij,mit)=vi3(ij)
     enddo

     !----> calculate cmm
     !

     ! See coefficients c_mi in Srivastava - why only c_mm???
     ! c_mi should change from iteration to iteration???
     ! is only 1 term of sum used???
     ! Is there a relation between cmm and am(i) ???
     ! cmm = v^T(i) . G . F(m)   i=m ???
     ! cmm = 1/alpha * (V_out[V_in] - V_in] . (\Delta V_out / ||\Delta V_out||)
     cmm_local = ddot(imap,fm,1,vi3,1)

     call MPI_Allreduce(cmm_local,cmm_global,1, &
     MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

     !----> update rho(m+1)
     !
     ! V_in_new = (1 - cmm) * (alpha*\Delta V_out + \Delta V_in) + V_in
     call daxpy(imap,one-cmm_global,ui3,1,sm,1)

   end if
   !=====  For MIT GT 1 activ  ==============================================

   !      MIT = MIT + 1
   !----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
   do ij = 1,imap
     fm1(ij) = fm(ij)
     sm1(ij) = sm(ij)
   end do

 end subroutine

end module
