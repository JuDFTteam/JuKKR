#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!> Collection of the Jij calculation routines - have fun.
module jij_calc_mod

use JijData_mod

! Let this point to the JijData workspace!
type(JijData), save, pointer :: global_jij_data => null()

CONTAINS

! ************************************************************************
subroutine clsjij( &
     i1,naez,rr,nr,rbasis,rcut,nsymat, &
     isymindex, &
     ixcp,nxcp,nxij,rxij,rxccls,zkrxij, &
     !                       new parameters after inc.p removal
     nrd, nxijd)
  ! ************************************************************************
  ! This subroutine is used to create the clusters around each atom 
  ! where Jij's are calculated
  !
  ! called by: main2
  !
  ! STRATEGY :
  ! Calculate the cluster of each atom by the lattice
  ! parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
  ! this routine is executed only for the assigned atom I1 of the 
  ! responsible processor
  !                                                          A.Thiess 7/2009
  ! ************************************************************************
  implicit none

  integer nrd
  integer nxijd

  !     .. array arguments
  !
  double precision rbasis(3,naez),   & ! pos. of basis atoms in EZ 
       rr(3,0:nrd),      & ! set of lattice vectors
       rxij(nxijd),      & ! interatomic distance Ri-Rj
       zkrxij(48,3,nxijd)  ! enters the exp-factor of G in kkrmat01
  integer          ixcp(nxijd),      & ! index to atom in elem/cell at site in cluster
       nxcp(nxijd),      & ! index to bravais lattice  at site in cluster
       isymindex(48)
  !
  !
  !     .. scalar arguments
  !
  double precision rcut
  integer          nxij, &          ! number of atoms in cluster 
       naez, &          ! number of atoms in EZ 
       nr, &            ! number of lattice vectors RR 
       i1, &            ! processor calling this routine deals with atom I1
       nsymat
  !
  !     .. local arrays
  !
  double precision rxccls(3,nxijd), & ! real space pos of atom in cluster
       tmp(3), &
       ircls(3,nxijd),rsort(nxijd), &
       rmat(64,3,3)
  integer          iixcp(nxijd),inxcp(nxijd),isort(nxijd)
  !
  !     .. local scalars
  !
  double precision epsshl,rcut2,rtmp
  integer          iaez,ib,id,ir,iv,ix,pos
  character*10     rotname(64)
  !
  !
  external         xsort
  intrinsic        sqrt
  !
  data             epsshl   / 1.0d-4 /
  !
  ! ------------------------------------------------------------------------
  ! This is generating the clusters which have a distance smaller
  ! than RCUT and RCUTXY in plane .
  ! The cluster atoms are ordered with radious and then z>y>x 
  ! The ordering allows an easy comparison of clusters.
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  rcut2   = (rcut+epsshl)*(rcut+epsshl)
  !======================================================================
  ! loop in all atoms begin
  !======================================================================

  nxij = 0           ! counter for atoms in cluster
  do iaez = 1,naez  ! loop in all atoms
     do ir = 0, nr    ! loop in all bravais vectors    
        do iv=1,3
           tmp(iv) = rr(iv,ir)+rbasis(iv,iaez)-rbasis(iv,i1)
        enddo
        rtmp   =  tmp(3)**2 + tmp(1)**2+tmp(2)**2
        if (rtmp.le.rcut2)  then
           nxij = nxij + 1
           if (nxij.gt.nxijd) then 
              write (6,*) &
                   ' ERROR: Dimension NXIJD in inc.cls too small', &
                   nxij, nxijd
              stop '   < CLSJIJ >'
           endif
           !
           ixcp(nxij) = iaez       ! store the atom in elem cell
           nxcp(nxij) = ir         ! store the bravais vector
           !
           do iv=1,3
              rxccls(iv,nxij) = tmp(iv)
           enddo
        endif
     enddo              ! IR loop in bravais

  enddo                 ! IAEZ loop in NAEZ
  !
  !     sort the atoms of the cluster in increasing order. First by distance
  !     Then by z then by y
  !
  do ix=1,nxij
     rsort(ix) = sqrt(rxccls(1,ix)**2+ &
          rxccls(2,ix)**2+ &
          rxccls(3,ix)**2)
     rsort(ix) = 100000000.d0*rsort(ix)+ &
          10000.d0*rxccls(3,ix)+ &
          10.d0*rxccls(2,ix)+ &
          0.1d0*rxccls(1,ix)
  enddo
  !
  call xsort(rsort,isort,nxij,pos)
  !
  !     Rearange exchange ia with ib
  ! MAP temporarily to another array
  !
  do ix=1,nxij
     do iv=1,3
        ircls(iv,ix)    = rxccls(iv,ix)
     enddo
     iixcp(ix) = ixcp(ix)
     inxcp(ix) = nxcp(ix)
  enddo
  !
  ! Now use correct order
  !
  call pointgrp(rmat,rotname)
  !
  do ix =1,nxij
     ib = isort(ix)
     do iv=1,3
        rxccls(iv,ix) = ircls(iv,ib)
     enddo
     ixcp(ix) = iixcp(ib)
     nxcp(ix) = inxcp(ib) 
     rxij(ix) = &
          sqrt(rxccls(1,ix)**2+rxccls(2,ix)**2+rxccls(3,ix)**2) ! store interatomic distance
     !
     do id = 1,nsymat
        do iv = 1,3
           zkrxij(id,iv,ix) = rmat(isymindex(id),iv,1)*rxccls(1,ix) + &
                rmat(isymindex(id),iv,2)*rxccls(2,ix) + &
                rmat(isymindex(id),iv,3)*rxccls(3,ix) - &
                rbasis(iv,ixcp(ix)) + &
                rbasis(iv,i1)    !ART
        enddo
     enddo
     !
  enddo
  !
end subroutine clsjij


!------------------------------------------------------------------------------
!     WARNING: Assumes one-to-one mapping of atom to rank in
!              'communicator'
!     NOT COMPATIBLE WITH CUTOFFS
!
!     GSXIJ has to be set to (0,0) before k-loop!
!     call this routine for every k-point
!     output: only GSXIJ
subroutine kkrjij( &
     kpoint,k_weight, &
     nsymat,naez,i3, &
     nxij,ixcp,zkrxij, &
     gllke1, &
     gsxij, &
     communicator, &
     lmmaxd, nxijd)

  implicit none
  ! ------------------------------------------------------------------------
  ! a) performs scattering of the off-diagonal GLLKE-elements
  !    required to calculate Jij's
  ! b) multiplication with appropriate exp-factor, k-dependent
  !                                                       A. Thiess Sep'09
  !    moderately successful "simplification": Elias Rabel, 2012
  ! ------------------------------------------------------------------------
  include 'mpif.h'

  integer communicator
  integer comm_size
  integer lmmaxd
  integer nxijd

  integer          nsymaxd
  parameter        (nsymaxd=48)

  double complex   czero
  parameter        (czero=(0.0d0,0.0d0))
  double complex   cione
  parameter        (cione  = ( 0.0d0,-1.0d0))
  !     ..
  !     .. GLOBAL SCALAR ARGUMENTS ..
  integer          naez, &         ! number of atoms per unit cell
  nsymat, &      ! number of active symmetries  
  nxij           ! number of atoms in Jij-cluster
  !     ..

  !     .. ARRAY ARGUMENTS ..
  double complex   gsxij(lmmaxd,lmmaxd,nsymaxd,nxijd)

  double precision kpoint(3), &   
  k_weight, &
  zkrxij(48,3,nxijd)                ! position of atoms and sites
                                    ! connected by symmetry
  integer          ixcp(nxijd)      ! corresp. to Jij no. XIJ on the
                                    ! real space lattice

  !double complex   gllke1(naez * lmmaxd, lmmaxd)
  double complex   gllke1(:, :)

  !     .. LOCAL ARRAYS ..
  !     .. Fortran 90 automatic arrays
  double complex   gsend(lmmaxd, lmmaxd)
  double complex   grecv(lmmaxd, lmmaxd)
  double complex   gxij (lmmaxd, lmmaxd, nxijd)  ! large?
  double complex   ekrxij(48, nxijd)

  integer          ixcps(nxijd)       ! copydummy for IXCP used MPI

  !     .. INTRINSIC FUNCTIONS ..
  intrinsic atan,exp


  !     .. LOCAL SCALARS ..
  integer          i3,i5,iv,isym,lm1,lm2,lm,ilm,xij
  double complex   carg

  !     .. MPI ..
  integer, dimension(mpi_status_size) :: status

  !     .. N-MPI
  integer          mapblock,ierr

  call MPI_Comm_Size(communicator, comm_size, ierr)

  CHECKASSERT(naez == comm_size)
  CHECKASSERT(size(GLLKE1, 1) == naez * lmmaxd)
  CHECKASSERT(size(GLLKE1, 2) == lmmaxd)

  ! ------------------------------------------------------------------------

  ! ================================================================
  !       XCCPL communicate off-diagonal elements
  !       Note: O(N**2) scaling!
  ! ================================================================

  do i5 = 1, naez
     ! .....
     do xij = 1, nxijd
        ! ....
        ixcps(xij) = 0
        if (xij.le.nxij) ixcps(xij) = ixcp(xij)
        ! ....
     enddo
     !
     !         broadcast IXCPS from processor of atom I5 to all
     !
     call mpi_bcast(ixcps,nxijd,mpi_integer, &
          mapblock(i5,1,naez,1,0,comm_size-1), &
          communicator,ierr)
     !
     call mpi_barrier(communicator,ierr)
     !
     do xij = 2, nxij
        ! ....
        !           if local GLL(k,E) is required just copy
        !
        if (ixcps(xij).eq.i5) then
           ! ...
           if (i3.eq.i5) then
              ! ..
              do lm = 1, lmmaxd
                 ! .
                 ilm = lmmaxd*(i3-1) + 1
                 call zcopy(lmmaxd,gllke1(ilm,lm),1, &
                      gxij(1,lm,xij),1)
                 ! .
              enddo
              ! ..
           endif
           ! ...
           !           required GLL(k,E) send by IXCPS(XIJ).EQ.I3
           !           and recieved by I5.EQ.I3
           !
        else
           ! ...
           if (ixcps(xij).eq.i3) then
              ! ..
              do lm = 1, lmmaxd
                 ! .
                 ilm = lmmaxd*(i5-1) + 1
                 call zcopy(lmmaxd,gllke1(ilm,lm),1, &
                      gsend(1,lm),1)
                 ! .
              enddo
              !
              call mpi_send(gsend,lmmaxd*lmmaxd,mpi_double_complex, &
                   mapblock(i5,1,naez,1,0,comm_size-1), &
                   99,communicator,ierr)
              ! ..
           endif


           if (i5.eq.i3) then
              ! ..
              call mpi_recv(grecv,lmmaxd*lmmaxd,mpi_double_complex, &
                   mapblock(ixcps(xij),1,naez,1,0,comm_size-1), &
                   99,communicator,status,ierr)

              do lm = 1, lmmaxd
                 ! .
                 call zcopy(lmmaxd,grecv(1,lm),1, &
                      gxij(1,lm,xij),1)
                 ! .
              enddo
              ! ..
           endif

           call mpi_barrier(communicator,ierr)
           ! ...
        endif
        ! ....
     enddo
     ! .....
  enddo

  ! ================================================================
  !       XCCPL calculate exponential factor and multiply with GXIJ
  ! ================================================================

  do xij = 2, nxij
     do isym  = 1,nsymat

        carg = czero
        do iv = 1,3
           carg =  carg + zkrxij(isym,iv,xij)*kpoint(iv)
        enddo

        ekrxij(isym,xij) = &
             k_weight * exp(carg*(cione*8.d0*atan(1.d0)))

     enddo
  enddo

  do xij = 2, nxij
     do isym = 1,nsymat
        do lm2=1,lmmaxd
           do lm1=1,lmmaxd
              gsxij(lm1,lm2,isym,xij) = &
                   gsxij(lm1,lm2,isym,xij) + &
                   ekrxij(isym,xij) * gxij(lm1,lm2,xij)
           enddo
        enddo
     enddo             ! ISYM = 1,NSYMAT
  enddo

  ! ================================================================
  !       XCCPL calculated integrated over KPT GSXIJ for off-diag-
  !       onal elements
  !       now back to KLOOPZ1 ..
  ! ================================================================

  return

end subroutine kkrjij

!------------------------------------------------------------------------------
subroutine symjij( &
     alat,tauvbz, &
     nsymat,dsymll, &
     nxij,ixcp, &
     tmatll,mssq, &
     gsxij, &
     gmatxij, &
     !                       new input parameters after inc.p removal
     naez, lmmaxd, nxijd)
  ! =======================================================================
  !
  !  a) symmetrize GSXIJ
  !
  !  b) GMATXIJ = - Delta_t^(-1) * GLL * Delta_t^(-1)
  !                                                     A.Thiess Sep'09
  ! =======================================================================
  !
  implicit none

  integer naez
  integer lmmaxd
  integer nxijd

  integer          nsymaxd
  parameter        (nsymaxd=48)

  double complex   cone,czero
  parameter        (cone  = ( 1.0d0,0.0d0))
  parameter        (czero  = ( 0.0d0,0.0d0))
  !     ..
  !     .. Scalar Arguments ..
  !     ..
  double complex   tauvbz
  double precision alat
  integer          nsymat,xij,nxij
  !     ..
  !     .. Array Arguments ..
  !     ..

  double complex   gmatxij(lmmaxd,lmmaxd,nxijd)
  double complex   gsxij  (lmmaxd,lmmaxd,nsymaxd,nxijd)
  double complex   tmatll (lmmaxd,lmmaxd,naez)
  double complex   mssq   (lmmaxd,lmmaxd)
  double complex   dsymll (lmmaxd,lmmaxd,nsymaxd)

  integer          ixcp(nxijd)
  !     ..
  !     .. Local Scalars ..
  !     ..
  double precision rfctor
  integer          iu,lm1,lm2,info
  !     ..
  !     .. Local Arrays ..
  !     ..

  double complex       gll(lmmaxd,lmmaxd)
  double complex   mssxcpl(lmmaxd,lmmaxd)
  double complex       tpg(lmmaxd,lmmaxd)
  double complex        xc(lmmaxd,lmmaxd)
  double complex        w1(lmmaxd,lmmaxd)
  integer             ipvt(lmmaxd)
  !     ..
  external zcopy,zaxpy,zgetrf,zgetrs,zgetri,zgemm,zscal
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,dble

  !     ..
  !     ! = ALAT/(2*PI)
  rfctor = alat/(8.d0*atan(1.0d0))



  !================================
  do xij = 2, nxij
     !================================
     do lm2 = 1,lmmaxd
        do lm1 = 1,lmmaxd
           mssxcpl(lm1,lm2) = tmatll(lm1,lm2,ixcp(xij))
        end do
     end do
     !
     ! ---> inversion 
     !
     call zgetrf(lmmaxd,lmmaxd,mssxcpl,lmmaxd,ipvt,info)
     call zgetri(lmmaxd,mssxcpl,lmmaxd,ipvt,w1, &
          lmmaxd*lmmaxd,info)
     !
     !
     !-------------------------------------------------------- SYMMETRISE GLL
     !
     do iu = 1,nsymat
        !
        ! --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^T)
        !
        if ( iu.eq.1 ) then
           !
           ! --->    ull(1) is equal to unity matrix
           !
           call zcopy(lmmaxd*lmmaxd,gsxij(1,1,1,xij),1,gll,1)
           call zscal(lmmaxd*lmmaxd,tauvbz,gll,1)
           !
        else
           !
           ! --->      tpg = tauvbz * DLL * GS
           !
           call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,tauvbz, &
                dsymll(1,1,iu),lmmaxd,gsxij(1,1,iu,xij),lmmaxd, &
                czero,tpg,lmmaxd)
           !
           ! --->    GLL = GLL + TPG * DLL(i)^T
           !
           call zgemm('N','C',lmmaxd,lmmaxd,lmmaxd,cone,tpg,lmmaxd, &
                dsymll(1,1,iu),lmmaxd,cone,gll,lmmaxd)
        end if
        !
     end do

     !-------------------------------------------------------- IU = 1,NSYMAT
     !
     ! In case of more than one atom per unit cell different MSSQ's have to
     ! be prepared
     !
     ! --->  XC = Delta_t^(-1) * GLL
     !
     call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,mssq, &
          lmmaxd,gll,lmmaxd,czero,xc,lmmaxd)
     !
     ! --->  GLL = - Delta_t^(-1) * GLL * Delta_t^(-1)
     !
     call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,-cone,xc,lmmaxd, &
          mssxcpl,lmmaxd,-czero,gll,lmmaxd)
     !
     ! --->  GMATXIJ = GMATLL = GLL/RFCTOR
     !
     do lm1 = 1,lmmaxd
        do lm2 = 1,lmmaxd
           gmatxij(lm2,lm1,xij) = gll(lm2,lm1)/rfctor
        end do
     end do

     !================================
  enddo
  !================================
  !xccpl
  !
  return

end subroutine symjij

!------------------------------------------------------------------------------

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
  external         ZCOPY
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
end subroutine

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

end subroutine


!------------------------------------------------------------------------------
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

end subroutine


!------------------------------------------------------------------------------
subroutine pointgrp(rotmat,rotname)
  ! **********************************************
  ! This subroutine defines the rotation matrices for
  ! all the 32 point groups and names them after
  ! J.F. Cornwell (Group Theory??) second edition
  ! Appendix D, p 324-325
  !
  ! *********************************************
  implicit none
  integer i,j,i1,is
  double precision rotmat(64,3,3)
  double precision rthree,half
  character*10 rotname(64)

  rthree = sqrt(3.d0)/2.d0
  half = 0.5d0
  ! set to zero
  do i1=1,64
     do i=1,3
        do j=1,3
           rotmat(i1,i,j) = 0.d0
        end do
     end do
  end do
  !
  rotmat(1,1,1) =  1.d0
  rotmat(1,2,2) =  1.d0
  rotmat(1,3,3) =  1.d0
  rotname(1) = 'E'
  !
  rotmat(2,1,2) =  1.d0
  rotmat(2,2,3) = -1.d0
  rotmat(2,3,1) = -1.d0
  rotname(2) = 'C3alfa'
  !
  rotmat(3,1,2) = -1.d0
  rotmat(3,2,3) = -1.d0
  rotmat(3,3,1) =  1.d0
  rotname(3) = 'C3beta '
  !
  rotmat(4,1,2) = -1.d0
  rotmat(4,2,3) =  1.d0
  rotmat(4,3,1) = -1.d0
  rotname(4) = 'C3gamma'
  !
  rotmat(5,1,2) = 1.d0
  rotmat(5,2,3) = 1.d0
  rotmat(5,3,1) = 1.d0
  rotname(5) = 'C3delta '
  !
  rotmat(6,1,3) = -1.d0
  rotmat(6,2,1) =  1.d0
  rotmat(6,3,2) = -1.d0
  rotname(6) = 'C3alfa-1'
  !
  rotmat(7,1,3) =  1.d0
  rotmat(7,2,1) = -1.d0
  rotmat(7,3,2) = -1.d0
  rotname(7) = 'C3beta-1 '
  !
  rotmat(8,1,3) = -1.d0
  rotmat(8,2,1) = -1.d0
  rotmat(8,3,2) =  1.d0
  rotname(8) = 'C3gamma-1'
  !
  rotmat(9,1,3) =  1.d0
  rotmat(9,2,1) =  1.d0
  rotmat(9,3,2) =  1.d0
  rotname(9) = 'C3delta-1'
  !
  rotmat(10,1,1) =  1.d0
  rotmat(10,2,2) = -1.d0
  rotmat(10,3,3) = -1.d0
  rotname(10) = 'C2x'
  !
  rotmat(11,1,1) = -1.d0
  rotmat(11,2,2) =  1.d0
  rotmat(11,3,3) = -1.d0
  rotname(11) = 'C2y'
  !
  rotmat(12,1,1) = -1.d0
  rotmat(12,2,2) = -1.d0
  rotmat(12,3,3) =  1.d0
  rotname(12) = 'C2z'
  !
  rotmat(13,1,1) =  1.d0
  rotmat(13,2,3) =  1.d0
  rotmat(13,3,2) = -1.d0
  rotname(13) = 'C4x'
  !
  rotmat(14,1,3) = -1.d0
  rotmat(14,2,2) =  1.d0
  rotmat(14,3,1) =  1.d0
  rotname(14) = 'C4y '
  !
  rotmat(15,1,2) =  1.d0
  rotmat(15,2,1) = -1.d0
  rotmat(15,3,3) =  1.d0
  rotname(15) = 'C4z'
  !
  rotmat(16,1,1) =  1.d0
  rotmat(16,2,3) = -1.d0
  rotmat(16,3,2) =  1.d0
  rotname(16) = 'C4x-1 '
  !
  rotmat(17,1,3) =  1.d0
  rotmat(17,2,2) =  1.d0
  rotmat(17,3,1) = -1.d0
  rotname(17) = 'C4y-1'
  !
  rotmat(18,1,2) = -1.d0
  rotmat(18,2,1) =  1.d0
  rotmat(18,3,3) =  1.d0
  rotname(18) = 'C4z-1'
  !
  rotmat(19,1,2) =  1.d0
  rotmat(19,2,1) =  1.d0
  rotmat(19,3,3) = -1.d0
  rotname(19) = 'C2a'
  !
  rotmat(20,1,2) = -1.d0
  rotmat(20,2,1) = -1.d0
  rotmat(20,3,3) = -1.d0
  rotname(20) = 'C2b'
  !
  rotmat(21,1,3) =  1.d0
  rotmat(21,2,2) = -1.d0
  rotmat(21,3,1) =  1.d0
  rotname(21) = 'C2c'
  !
  rotmat(22,1,3) = -1.d0
  rotmat(22,2,2) = -1.d0
  rotmat(22,3,1) = -1.d0
  rotname(22) = 'C2d'
  !
  rotmat(23,1,1) = -1.d0
  rotmat(23,2,3) =  1.d0
  rotmat(23,3,2) =  1.d0
  rotname(23) = 'C2e'
  !
  rotmat(24,1,1) = -1.d0
  rotmat(24,2,3) = -1.d0
  rotmat(24,3,2) = -1.d0
  rotname(24) = 'C2f'
  do i1=1,24
     do i=1,3
        do j=1,3
           rotmat(i1+24,i,j) = -rotmat(i1,i,j)
        end do
     end do
     rotname(i1+24) = 'I'//rotname(i1)
  end do
  !
  !
  !*********************************************
  ! Trigonal and hexagonal groups
  !*********************************************
  !
  rotmat(49,1,1) = -half
  rotmat(49,1,2) =  rthree
  rotmat(49,2,1) = -rthree
  rotmat(49,2,2) = -half
  rotmat(49,3,3) =  1.d0
  rotname(49) = 'C3z'
  !
  rotmat(50,1,1) = -half
  rotmat(50,1,2) = -rthree
  rotmat(50,2,1) =  rthree
  rotmat(50,2,2) = -half
  rotmat(50,3,3) =  1.d0
  rotname(50) = 'C3z-1'
  !
  rotmat(51,1,1) =  half
  rotmat(51,1,2) =  rthree
  rotmat(51,2,1) = -rthree
  rotmat(51,2,2) =  half
  rotmat(51,3,3) =  1.d0
  rotname(51) = 'C6z'
  !
  rotmat(52,1,1) =  half
  rotmat(52,1,2) = -rthree
  rotmat(52,2,1) =  rthree
  rotmat(52,2,2) =  half
  rotmat(52,3,3) =  1.d0
  rotname(52) = 'C6z-1'
  !
  rotmat(53,1,1) = -half
  rotmat(53,1,2) =  rthree
  rotmat(53,2,1) =  rthree
  rotmat(53,2,2) =  half
  rotmat(53,3,3) = -1.d0
  rotname(53) = 'C2A'
  !
  rotmat(54,1,1) = -half
  rotmat(54,1,2) = -rthree
  rotmat(54,2,1) = -rthree
  rotmat(54,2,2) =  half
  rotmat(54,3,3) = -1.d0
  rotname(54) = 'C2B'
  !
  rotmat(55,1,1) =  half
  rotmat(55,1,2) = -rthree
  rotmat(55,2,1) = -rthree
  rotmat(55,2,2) = -half
  rotmat(55,3,3) = -1.d0
  rotname(55) = 'C2C'
  !
  rotmat(56,1,1) =  half
  rotmat(56,1,2) =  rthree
  rotmat(56,2,1) =  rthree
  rotmat(56,2,2) = -half
  rotmat(56,3,3) = -1.d0
  rotname(56) = 'C2D'
  do is=1,8
     do i=1,3
        do j=1,3
           rotmat(56+is,i,j) = -rotmat(48+is,i,j)
        end do
     end do
     rotname(56+is) = 'I'//rotname(48+is)
  end do

  !ccccccccccccccccccccccccccccccccccccccccccccccccc
end subroutine pointgrp

end module
