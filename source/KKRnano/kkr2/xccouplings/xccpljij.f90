!     TODO: BUG - energy resolved Jijs - does not work as expected
!     - there exists some confusion with writing to the Eij files

!=======================================================================

!>    This routine is to be included in a loop over energy points
!>    @param IER energy point index
!>    @param WGTE energy weight factor/nspin
!>    @param IXCP cluster data?
!>    @param RXCCLS cluster data?
!>    @param GMATXIJ off-diagonal Greens function elements?
!>    @param DTIXIJ  \Delta T_up - \Delta T_down
!>    @param ERESJIJ true: write energy resolved Jij
!>    @param JXCIJINT integrated Jij(E), for IER=1 they are initialised
!>           then pass old JXCIJINT for next energy point -
!>           therefore integration is achieved: in,out
subroutine XCCPLJIJ_START( &
I1,IER,WGTE, &
RXIJ,NXIJ,IXCP,RXCCLS, &
GMATXIJ,DTIXIJ, &
communicator, &
JXCIJINT,ERESJIJ, &
naez, lmmaxd, nxijd, nspind)
  !   ********************************************************************
  !   *                                                                  *
  !   *  calculates the site-off diagonal  XC-coupling parameters  J_ij  *
  !   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
  !   *                                                                  *
  !   *  adopted for TB-KKR code from Munich SPR-KKR package Sep 2004    *
  !   *  adopted for KKRnano, Jun 2009                                   *
  !   ********************************************************************

  implicit none

  INCLUDE 'mpif.h'

  integer, intent(in) :: naez
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: nxijd
  integer, intent(in) :: nspind

  !     ..
  !     .. Scalar arguments
  double complex, intent(in) :: WGTE
  integer, intent(in) :: I1
  integer, intent(in) :: IER
  integer, intent(in) :: NXIJ
  logical, intent(in) :: ERESJIJ

  double complex :: GMATXIJ(lmmaxd,lmmaxd,NXIJD,NSPIND)
  double complex :: DTIXIJ(lmmaxd,lmmaxd)
  double precision :: RXIJ(NXIJD)
  double precision :: RXCCLS(3,NXIJD)
  integer :: IXCP(NXIJD)
  integer, intent(in) :: communicator

  double complex :: JXCIJINT(NXIJD)

  !     .. Parameters
  double complex :: CONE
  double complex :: CZERO
  parameter        ( CONE  = (1D0,0D0) )
  parameter        ( CZERO = (0D0,0D0) )

  !     ..
  !     .. Local scalars
  integer :: XIJ
  integer :: ISPIN
  integer :: LM1
  integer :: LM2
  integer :: D1
  integer :: D10
  integer :: D100
  integer :: D1000
  double complex :: CSUM
  double complex :: JSCAL
  double precision :: JOUT
  character(len=12)::FNAME
  !     ..
  !     .. Local arrays
  integer :: OFF(3)

  double complex :: GMIJ_down(lmmaxd,lmmaxd)
  double complex :: GMJI_up(lmmaxd,lmmaxd)
  double complex :: W1(lmmaxd,lmmaxd)
  double complex :: W2(lmmaxd,lmmaxd)
  double complex :: W3(lmmaxd,lmmaxd)

  !     large local array
  double complex, allocatable, dimension(:,:,:) :: DTNXIJ_ALL

  !     .. MPI ..
  integer :: IERR

  !     ..
  !     .. Intrinsic Functions ..
  intrinsic        SQRT
  !     ..
  !     .. External Subroutines ..
  external         CINIT,CMATMUL,ZCOPY
  !     ..

  integer :: memory_stat
  logical :: memory_fail

  memory_fail = .false.

  JSCAL = CONE/4D0

  ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
  ! ==>                   IE.EQ.1 -- initialisation step --            <==
  ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

  !     Allocate array

  allocate(DTNXIJ_ALL(LMMAXD,LMMAXD,NAEZ), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail .eqv. .true.) then
    write(*,*) "XCCPLJIJ: FATAL Error, failure to allocate memory."
    write(*,*) "       Probably out of memory."
    stop
  end if

  if ( IER==1 ) then

    JXCIJINT = CZERO

    if (ERESJIJ) then
!      D1 = mod(I1,10)
!      D10 = int( (mod(I1,100) + 0.5)/10 )
!      D100 = int( (mod(I1,1000) + 0.5)/100 )
!      D1000 = int( (mod(I1,10000) + 0.5)/1000 )
!
!      OFF(1) = iachar('1')-1
!      OFF(2) = iachar('1')-1
!      OFF(3) = iachar('1')-1
!
!      if ( D10>=10 ) OFF(1) = iachar('7')
!      if ( D100>=100 ) OFF(2) = iachar('7')
!      if ( D1000>=1000 ) OFF(3) = iachar('7')
!
!      FNAME='Eij.' &
!      //achar(D1000+OFF(3)) &
!      //achar(D100+OFF(2)) &
!      //achar(D10+OFF(1)) &
!      //achar(D1+iachar('1')-1) &
!      //'.dat'
!
!      open(75,file=FNAME,form='formatted')
    endif
  endif

  ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
  ! ==>                   INITIALISATION END                           <==
  ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  ! ==>  get DTJXIJ = T(UP) - T(DOWN) for all atoms


  call MPI_ALLGATHER(DTIXIJ,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX, &
  DTNXIJ_ALL,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX, &
  communicator,IERR)


  do XIJ = 2, NXIJ ! loop XIJ = 1, NXIJ(I1)

    ! ==>  get the off-diagonal Green function matrix Gij(UP) and Gji(DOWN)

    do ISPIN = 1,2
      if ( ISPIN==1 ) then
        call ZCOPY(LMMAXD*LMMAXD,GMATXIJ(1,1,XIJ,ISPIN), &
        1,GMIJ_down,1)
      else
        do LM2 = 1,LMMAXD
          do LM1 = 1,LMMAXD

            ! -> use Gji = Gij^T

            GMJI_up(LM1,LM2) = GMATXIJ(LM2,LM1,XIJ,ISPIN)

          enddo
        enddo
      endif
    enddo

    ! ----------------------------------------------------------------------

    ! ==> calculate the exchange coupling constant J_ij via Eq. (19)
    !     modified for G instead of tau:
    !          J_ij ~ Trace [ (t_i(D)-t_i(U)) * Gij(U)
    !                       * (t_j(D)-t_j(U)) * Gji(D)]

    ! -------------------------------------------------- loop over occupants
    ! --> Delta_j * Gjt,it

    call CMATMUL(LMMAXD,LMMAXD,DTNXIJ_ALL(1,1,IXCP(XIJ)),GMJI_up,W2)

    ! --> Delta_i * Git,jt

    call CMATMUL(LMMAXD,LMMAXD,DTIXIJ,GMIJ_down,W3)

    ! --> Delta_i * Git,jt * Delta_j * Gjt,it

    call CMATMUL(LMMAXD,LMMAXD,W3,W2,W1)

    CSUM = CZERO
    do LM1 = 1,LMMAXD
      CSUM = CSUM + W1(LM1,LM1)
    enddo

    JOUT = -DIMAG(WGTE*CSUM*JSCAL)

    if (ERESJIJ) then
    !  write(75,73002) &
    !  IER,XIJ,RXIJ(XIJ),JOUT, &
    !  RXCCLS(1,XIJ),RXCCLS(2,XIJ),RXCCLS(3,XIJ),IXCP(XIJ)
    endif

    JXCIJINT(XIJ) = JXCIJINT(XIJ) - WGTE*CSUM

  !                  -------> perform substraction instead of addition
  !                           because WGTE ~ -1/pi
  ! ======================================================================
  enddo             ! loop XIJ = 1, NXIJ
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  deallocate(DTNXIJ_ALL)

  !if (ERESJIJ) close(75) ! WRONG: close only after last energy point!!!

73002 format(I3,1X,I3,3X,F9.5,6X,D9.3,6X,3(1X,F7.4),I5)
end

!-----------------------------------------------------------------------
subroutine writeJiJs(I1, &
                     RXIJ,NXIJ,IXCP,RXCCLS, &
                     JXCIJINT, nxijd)
  implicit none

  integer :: I1
  integer :: NXIJ
  integer :: nxijd
  double precision :: RXIJ(NXIJD)
  double precision :: RXCCLS(3,NXIJD)
  integer :: IXCP(NXIJD)
  double complex :: JXCIJINT(NXIJD)

  !     local variables
  double complex :: CONE
  parameter        ( CONE  = (1D0,0D0) )

  integer :: XIJ
  integer :: D1
  integer :: D10
  integer :: D100
  integer :: D1000

  integer :: OFF(3)

  double complex :: JSCAL
  character(len=12) :: FNAME


  if (nxij > nxijd) then
    write(*,*) "writeJijs: nxij > nxijd"
    stop
  endif

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  JSCAL = CONE/4D0
  ! .. ...........................................................
  ! write Jij's to file Jij.I1.dat
  ! ..
  D1 = mod(I1,10)
  D10 = int( (mod(I1,100) + 0.5)/10 )
  D100 = int( (mod(I1,1000) + 0.5)/100 )
  D1000 = int( (mod(I1,10000) + 0.5)/1000 )

  OFF(1) = iachar('1')-1
  OFF(2) = iachar('1')-1
  OFF(3) = iachar('1')-1

  if ( D10>=10 ) OFF(1) = iachar('7')
  if ( D100>=100 ) OFF(2) = iachar('7')
  if ( D1000>=1000 ) OFF(3) = iachar('7')

  FNAME='Jij.' &
  //achar(D1000+OFF(3)) &
  //achar(D100+OFF(2)) &
  //achar(D10+OFF(1)) &
  //achar(D1+iachar('1')-1) &
  //'.dat'

  open(73,file=FNAME,form='formatted')

  write(73,73000) I1

  do XIJ = 2, NXIJ

    JXCIJINT(XIJ) = JSCAL*JXCIJINT(XIJ)
    write(73,73001) &
    XIJ,RXIJ(XIJ),DIMAG(JXCIJINT(XIJ)), &
    RXCCLS(1,XIJ),RXCCLS(2,XIJ),RXCCLS(3,XIJ),IXCP(XIJ)

  enddo

  close(73)

73000 format("# off-diagonal exchange-coupling constants Jij ",/, &
  "# for atom i = ",I1,/, &
  "# j    R_ij( ALAT )   J_ij( Ry )      RXCCLS      ", &
  "             IXCP")
73001 format(I3,3X,F9.5,6X,D9.3,6X,3(1X,F7.4),I5)

end




subroutine CMATMUL(N,M,A,B,C)
  !   ********************************************************************
  !   *                                                                  *
  !   *   perform  the matrix-matrix operation           C = A * B       *
  !   *                                                                  *
  !   *   A,B,C   complex  SQUARE  N x N - matrices                      *
  !   *   N       dimension of A, B and C                                *
  !   *   M       array size of A, B, C with M >= N                      *
  !   *                                                                  *
  !   ********************************************************************
  implicit double complex(a-h,o-z)

  ! PARAMETER definitions

  double complex :: C0
  parameter (C0=(0.0D0,0.0D0))

  ! Dummy arguments

  integer :: M
  integer :: N
  double complex :: A(M,M)
  double complex :: B(M,M)
  double complex :: C(M,M)

  ! Local variables

  double complex :: BLJ
  integer :: I
  integer :: J
  integer :: L

  do J = 1,N
    do I = 1,N
      C(I,J) = C0
    end do
  end do

  do J = 1,N
    do L = 1,N
      BLJ = B(L,J)
      if ( BLJ/=C0 ) then
        do I = 1,N
          C(I,J) = C(I,J) + A(I,L)*BLJ
        end do
      end if
    end do
  end do

end
