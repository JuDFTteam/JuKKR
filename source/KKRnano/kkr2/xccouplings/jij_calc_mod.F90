#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!> Collection of the Jij calculation routines - have fun.
module jij_calc_mod
  use JijData_mod, only: JijData
  implicit none
  private
  
  public :: XCCPLJIJ_START, writeJiJs
  public :: kkrjij, symjij, clsjij
  
  ! Let this point to the JijData workspace!
  type(JijData), pointer, public :: global_jij_data => null()

  integer, parameter :: nsymaxd=48
  double complex, parameter :: cone=(1.d0,0.d0), czero=(0.d0,0.d0)
  
  contains

! ************************************************************************
subroutine clsjij(i1, naez, rr, nr, rbasis, rcut, nsymat, isymindex, ixcp, nxcp, nxij, rxij, rxccls, zkrxij, nrd, nxijd)
  use Sorting_mod, only: dsort
  use Symmetry_mod, only: pointgrp
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
  integer nrd
  integer nxijd
  double precision rbasis(3,naez),   & ! pos. of basis atoms in EZ 
       rr(3,0:nrd),      & ! set of lattice vectors
       rxij(nxijd),      & ! interatomic distance Ri-Rj
       zkrxij(48,3,nxijd)  ! enters the exp-factor of G in kkrmat01
  integer          ixcp(nxijd),      & ! index to atom in elem/cell at site in cluster
       nxcp(nxijd),      & ! index to bravais lattice  at site in cluster
       isymindex(48)
  double precision rcut
  integer          nxij, &          ! number of atoms in cluster 
       naez, &          ! number of atoms in EZ 
       nr, &            ! number of lattice vectors RR 
       i1, &            ! processor calling this routine deals with atom I1
       nsymat

  !     .. locals
  double precision rxccls(3,nxijd), & ! real space pos of atom in cluster
       tmp(3), &
       ircls(3,nxijd),rsort(nxijd), &
       rmat(64,3,3)
  integer :: iixcp(nxijd),inxcp(nxijd),isort(nxijd)
  double precision :: rcut2, rtmp
  integer :: iaez, ib, id, ir, iv, ix
  character*10 :: rotname(64)
  double precision, parameter :: epsshl = 1.d-4

  ! ------------------------------------------------------------------------
  ! This is generating the clusters which have a distance smaller
  ! than RCUT and RCUTXY in plane .
  ! The cluster atoms are ordered with radious and then z>y>x 
  ! The ordering allows an easy comparison of clusters.
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  rcut2 = (rcut + epsshl)**2
  !======================================================================
  ! loop in all atoms begin
  !======================================================================

  nxij = 0           ! counter for atoms in cluster
  do iaez = 1, naez  ! loop in all atoms
     do ir = 0, nr    ! loop in all bravais vectors    
        tmp(1:3) = rr(1:3,ir) + rbasis(1:3,iaez) - rbasis(1:3,i1)
        rtmp = sum(tmp(1:3)**2)
        if (rtmp <= rcut2)  then
           nxij = nxij + 1
           if (nxij > nxijd) stop '  Dimension error NXIJD < CLSJIJ >'
           !
           ixcp(nxij) = iaez       ! store the atom in elem cell
           nxcp(nxij) = ir         ! store the bravais vector
           rxccls(1:3,nxij) = tmp(1:3)
        endif
     enddo ! IR loop in bravais
  enddo ! IAEZ loop in NAEZ
  !
  !     sort the atoms of the cluster in increasing order. First by distance
  !     Then by z then by y
  do ix = 1, nxij
    rsort(ix) = sqrt(sum(rxccls(1:3,ix)**2))
!   rsort(ix) = 1.d8*rsort(ix) + 1.d4*rxccls(3,ix) + 1.d1*rxccls(2,ix) + 1.d-1*rxccls(1,ix)
    rsort(ix) = 1.d9*rsort(ix) + 1.d5*rxccls(3,ix) + 1.d2*rxccls(2,ix) + rxccls(1,ix) ! different weights than in Voronoi code
  enddo ! ix

  call dsort(rsort, isort, nxij)
  !
  !     Rearange exchange ia with ib
  ! MAP temporarily to another array
  !
  do ix = 1, nxij
    ircls(1:3,ix) = rxccls(1:3,ix)
    iixcp(ix) = ixcp(ix)
    inxcp(ix) = nxcp(ix)
  enddo ! ix
  
  ! Now use correct order
  call pointgrp(rmat, rotname) ! warning memory order for rmat has changed: (3,3,64)

  do ix = 1, nxij
     ib = isort(ix)
     rxccls(1:3,ix) = ircls(1:3,ib)
     ixcp(ix) = iixcp(ib)
     nxcp(ix) = inxcp(ib) 
     rxij(ix) = sqrt(sum(rxccls(1:3,ix)**2)) ! store interatomic distance

     do id = 1, nsymat
       do iv = 1, 3
         zkrxij(id,iv,ix) = dot_product(rmat(iv,1:3,isymindex(id)), rxccls(1:3,ix)) - rbasis(iv,ixcp(ix)) + rbasis(iv,i1) !ART
       enddo ! iv
     enddo ! id
  enddo ! ix

endsubroutine ! clsjij


!------------------------------------------------------------------------------
!     WARNING: Assumes one-to-one mapping of atom to rank in
!              'communicator'
!     NOT COMPATIBLE WITH CUTOFFS
!
!     GSXIJ has to be set to (0,0) before k-loop!
!     call this routine for every k-point
!     output: only GSXIJ
subroutine kkrjij(&
     kpoint,k_weight, &
     nsymat,naez,i3, &
     nxij,ixcp,zkrxij, &
     gllke1, &
     gsxij, &
     communicator, &
     lmmaxd, nxijd)
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

  !     ..
  !     .. GLOBAL SCALAR ARGUMENTS ..
  integer          naez, &         ! number of atoms per unit cell
  nsymat, &      ! number of active symmetries  
  nxij           ! number of atoms in Jij-cluster
  !     ..

  !     .. ARRAY ARGUMENTS ..
  double complex   gsxij(lmmaxd,lmmaxd,nsymaxd,nxijd)

  double precision kpoint(3), k_weight, &
  zkrxij(48,3,nxijd)                ! position of atoms and sites connected by symmetry
  integer          ixcp(nxijd)      ! corresp. to Jij no. XIJ on the real space lattice

  !double complex   gllke1(naez*lmmaxd,lmmaxd)
  double complex :: gllke1(:,:)

  
  double complex, parameter :: cione=(0.d0,-1.d0)
  !     .. LOCAL ARRAYS ..
  !     .. Fortran 90 automatic arrays
  double complex   gsend(lmmaxd, lmmaxd)
  double complex   grecv(lmmaxd, lmmaxd)
  double complex   gxij (lmmaxd, lmmaxd, nxijd)  ! large?
  double complex   ekrxij(48,nxijd)

  integer          ixcps(nxijd)       ! copydummy for IXCP used MPI

  !     .. INTRINSIC FUNCTIONS ..
  intrinsic atan,exp


  !     .. LOCAL SCALARS ..
  integer          i3,i5,isym,lm,ilm,xij,ierr
  double complex   carg, mtpii
  integer :: status(MPI_status_size)
  integer, external :: mapblock

  call MPI_Comm_Size(communicator, comm_size, ierr)

  CHECKASSERT(naez == comm_size)
  CHECKASSERT(size(GLLKE1, 1) == naez * lmmaxd)
  CHECKASSERT(size(GLLKE1, 2) == lmmaxd)

  ! ------------------------------------------------------------------------

  ! ================================================================
  !       XCCPL communicate off-diagonal elements
  !       Note: O(N**2) scaling!
  ! ================================================================

  mtpii = dcmplx(0.d0,-8.d0*atan(1.d0)) ! = -2*pi*i
  
  do i5 = 1, naez
     
    ixcps(:) = 0
    ixcps(1:nxij) = ixcp(1:nxij)
    
    !
    !         broadcast IXCPS from processor of atom I5 to all
    !
    call MPI_bcast(ixcps,nxijd,MPI_integer, mapblock(i5,1,naez,1,0,comm_size-1), communicator,ierr)
    !
    call MPI_barrier(communicator,ierr)
    !
    do xij = 2, nxij
    ! ....
    !           if local GLL(k,E) is required just copy
    !
      if (ixcps(xij) == i5) then
        ! ...
        if (i3 == i5) then
          do lm = 1, lmmaxd
            ilm = lmmaxd*(i3-1) + 1
            call zcopy(lmmaxd,gllke1(ilm,lm),1,gxij(1,lm,xij),1)
          enddo ! lm
        endif ! i3 == i5
          
          !           required GLL(k,E) sendby IXCPS(XIJ) == I3
          !           and recieved by I5 == I3
      else
      
        if (ixcps(xij) == i3) then
          do lm = 1, lmmaxd
            ilm = lmmaxd*(i5-1) + 1
            call zcopy(lmmaxd,gllke1(ilm,lm),1,gsend(1,lm),1)
          enddo ! lm

          call MPI_send(gsend,lmmaxd*lmmaxd,MPI_double_complex, mapblock(i5,1,naez,1,0,comm_size-1), 99,communicator,ierr)
        endif


        if (i5 == i3) then
          call MPI_recv(grecv,lmmaxd*lmmaxd,MPI_double_complex, mapblock(ixcps(xij),1,naez,1,0,comm_size-1), 99,communicator,status,ierr)

          do lm = 1, lmmaxd
            call zcopy(lmmaxd,grecv(1,lm),1,gxij(1,lm,xij),1)
          enddo ! lm
        endif

        call MPI_barrier(communicator, ierr)

      endif
    enddo ! xij
  enddo ! i5

  ! ================================================================
  !       XCCPL calculate exponential factor and multiply with GXIJ
  ! ================================================================

  do xij = 2, nxij
    do isym = 1, nsymat
      carg = dot_product(zkrxij(isym,1:3,xij), kpoint(1:3))
      ekrxij(isym,xij) = k_weight*exp(mtpii*carg)
      gsxij(:,:,isym,xij) = gsxij(:,:,isym,xij) + ekrxij(isym,xij) * gxij(:,:,xij)
    enddo ! ISYM = 1, NSYMAT
  enddo ! xij

  ! XCCPL calculated integrated over KPT GSXIJ for off-diagonal elements, now back to KLOOPZ1 ..
endsubroutine kkrjij

!------------------------------------------------------------------------------
subroutine symjij(&
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
  integer naez
  integer lmmaxd
  integer nxijd
  double precision tauvbz
  double precision alat
  integer          nsymat,xij,nxij
  double complex   gmatxij(lmmaxd,lmmaxd,nxijd)
  double complex   gsxij  (lmmaxd,lmmaxd,nsymaxd,nxijd)
  double complex   tmatll (lmmaxd,lmmaxd,naez)
  double complex   mssq   (lmmaxd,lmmaxd)
  double complex   dsymll (lmmaxd,lmmaxd,nsymaxd)
  integer          ixcp(nxijd)

  !     .. Locals

  double precision rfctor
  integer          iu,info
  double complex       gll(lmmaxd,lmmaxd)
  double complex   mssxcpl(lmmaxd,lmmaxd)
  double complex       tpg(lmmaxd,lmmaxd)
  double complex        xc(lmmaxd,lmmaxd)
  double complex        w1(lmmaxd,lmmaxd)
  integer             ipvt(lmmaxd)
  
  external :: zcopy,zaxpy,zgetrf,zgetrs,zgetri,zgemm,zscal ! BLAS & LAPACK

  rfctor = alat/(8.d0*atan(1.d0)) ! = ALAT/(2*PI)

  do xij = 2, nxij
  
     mssxcpl(:,:) = tmatll(:,:,ixcp(xij))

     ! inversion 
     call zgetrf(lmmaxd,lmmaxd,mssxcpl,lmmaxd,ipvt,info)
     call zgetri(lmmaxd,mssxcpl,lmmaxd,ipvt,w1,lmmaxd*lmmaxd,info)

     !-------------------------------------------------------- SYMMETRISE GLL
     ! --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^T)
        
     ! assume that for iu == 1 we have th unit matrix ==> ull(1) is equal to unity matrix
      call zcopy(lmmaxd*lmmaxd,gsxij(1,1,1,xij),1,gll,1)
      call zscal(lmmaxd*lmmaxd,cmplx(tauvbz, 0.d0),gll,1)
     
    do iu = 2, nsymat
      ! --->      tpg = tauvbz * DLL * GS
      call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cmplx(tauvbz, 0.d0),dsymll(1,1,iu),lmmaxd,gsxij(1,1,iu,xij),lmmaxd,czero,tpg,lmmaxd)
      !
      ! --->    GLL = GLL + TPG * DLL(i)^T
      call zgemm('N','C',lmmaxd,lmmaxd,lmmaxd,cone,tpg,lmmaxd,dsymll(1,1,iu),lmmaxd,cone,gll,lmmaxd)
    enddo ! iu

    !-------------------------------------------------------- IU = 1,NSYMAT
    !
    ! In case of more than one atom per unit cell different MSSQ's have to
    ! be prepared
    !
    ! --->  XC = Delta_t^(-1) * GLL
    call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,mssq,lmmaxd,gll,lmmaxd,czero,xc,lmmaxd)
    !
    ! --->  GLL = - Delta_t^(-1) * GLL * Delta_t^(-1)
    call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,-cone,xc,lmmaxd,mssxcpl,lmmaxd,-czero,gll,lmmaxd)
    !
    ! --->  GMATXIJ = GMATLL = GLL/RFCTOR
    gmatxij(:,:,xij) = gll(:,:)/rfctor
     
  enddo ! xij
  
endsubroutine symjij

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
subroutine XCCPLJIJ_START(i1, ier, wgte, rxij,nxij,ixcp,rxccls, gmatxij,dtixij, communicator, jxcijint,eresjij, naez, lmmaxd, nxijd, nspind) ! todo: remove i1, rxij, rxccls
  !   ********************************************************************
  !   *                                                                  *
  !   *  calculates the site-off diagonal  XC-coupling parameters  J_ij  *
  !   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
  !   *                                                                  *
  !   *  adopted for TB-KKR code from Munich SPR-KKR package Sep 2004    *
  !   *  adopted for KKRnano, Jun 2009                                   *
  !   ********************************************************************

  include 'mpif.h'

  integer, intent(in) :: naez
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: nxijd
  integer, intent(in) :: nspind

  !     ..
  !     .. scalar arguments
  double complex, intent(in) :: wgte
  integer, intent(in) :: i1
  integer, intent(in) :: ier
  integer, intent(in) :: nxij
  logical, intent(in) :: eresjij

  double complex :: gmatxij(lmmaxd,lmmaxd,nxijd,nspind)
  double complex :: dtixij(lmmaxd,lmmaxd)
  double precision :: rxij(nxijd)
  double precision :: rxccls(3,nxijd)
  integer :: ixcp(nxijd)
  integer, intent(in) :: communicator

  double complex :: jxcijint(nxijd)


  external :: zcopy
  !     .. locals
  integer :: xij
  integer :: ispin
  integer :: lm1
  double complex :: csum
  double complex :: jscal
  double precision :: jout
! character(len=16) :: filename
  double complex :: gmij_down(lmmaxd,lmmaxd)
  double complex :: gmji_up(lmmaxd,lmmaxd)
  double complex :: w1(lmmaxd,lmmaxd)
  double complex :: w2(lmmaxd,lmmaxd)
  double complex :: w3(lmmaxd,lmmaxd)
  double complex, allocatable :: dtnxij_all(:,:,:)
  integer :: memory_stat, ierr
  
  jscal = 0.25d0*cone

  ! ==>                   ie == 1 -- initialisation step --            <==

  allocate(dtnxij_all(lmmaxd,lmmaxd,naez), stat=memory_stat)
  if (memory_stat /= 0) stop 'xccpljij: fatal error, failure to allocate memory.'

  if (ier == 1) then

    jxcijint = czero

    if (eresjij) then
!      write(unit=filename, fmt="(a,i4.4,a)") 'eij.',i1,'.dat'
!      open(75, file=filename, form='formatted', action='write')
    endif
  endif

  ! ==>                   initialisation end                          <==

  ! ==>  get dtjxij = t(up) - t(down) for all atoms


  call MPI_allgather(dtixij,lmmaxd*lmmaxd,MPI_double_complex, dtnxij_all,lmmaxd*lmmaxd,MPI_double_complex, communicator,ierr)


  do xij = 2, nxij ! loop xij = 1, nxij(i1)

    ! ==>  get the off-diagonal green function matrix gij(up) and gji(down)

    do ispin = 1, 2
      if (ispin == 1) then
        call zcopy(lmmaxd*lmmaxd,gmatxij(1,1,xij,ispin),1,gmij_down,1)
      else
        ! -> use gji = gij^t
        gmji_up(:,:) = transpose(gmatxij(:,:,xij,ispin))
      endif
    enddo ! ispin

    ! ----------------------------------------------------------------------

    ! ==> calculate the exchange coupling constant j_ij via eq. (19)
    !     modified for g instead of tau:
    !          j_ij ~ trace [ (t_i(d)-t_i(u)) * gij(u)
    !                       * (t_j(d)-t_j(u)) * gji(d)]

    ! -------------------------------------------------- loop over occupants
    ! --> delta_j * gjt,it

    call cmatmul(lmmaxd,lmmaxd,dtnxij_all(1,1,ixcp(xij)),gmji_up,w2)

    ! --> delta_i * git,jt

    call cmatmul(lmmaxd,lmmaxd,dtixij,gmij_down,w3)

    ! --> delta_i * git,jt * delta_j * gjt,it

    call cmatmul(lmmaxd,lmmaxd,w3,w2,w1)

    csum = czero
    do lm1 = 1, lmmaxd
      csum = csum + w1(lm1,lm1)
    enddo ! lm1

    jout = -dimag(wgte*csum*jscal)

!   if (eresjij) write(75, fmt='(i3,1x,i3,3x,f9.5,6x,d9.3,6x,3(1x,f7.4),i5)') ier,xij,rxij(xij),jout,rxccls(1:3,xij),ixcp(xij)

    jxcijint(xij) = jxcijint(xij) - wgte*csum

    ! -------> perform substraction instead of addition because wgte ~ -1/pi
  enddo ! xij = 1, nxij

  deallocate(dtnxij_all, stat=memory_stat)

  ! if (eresjij) close(75) ! wrong: close only after last energy point!!!

endsubroutine

!-----------------------------------------------------------------------
subroutine writeJiJs(i1, rxij, nxij, ixcp, rxccls, jxcijint, nxijd)
  integer, intent(in) :: i1, nxij, nxijd
  double precision, intent(in)  :: rxij(nxijd)
  double precision, intent(in)  :: rxccls(3,nxijd)
  integer, intent(in)           :: ixcp(nxijd)
  double complex, intent(inout) :: jxcijint(nxijd)

  ! local variables
  integer :: xij, ios
  character(len=16) :: filename

  if (nxij > nxijd) stop "writeJijs: nxij > nxijd"

  ! write Jij's to file Jij.I1.dat
  write(unit=filename, fmt="(a,i4.4,a)", iostat=ios) 'Jij.',I1,'.dat' 
  open (73, file=filename, form='formatted', action='write', iostat=ios)
  write(73, fmt='(a)') "# off-diagonal exchange-coupling constants Jij "
  write(73, fmt='(a,i0)') "# for atom i = ",I1
  write(73, fmt='(a)') "# j    R_ij(ALAT)   J_ij(Ry)      RXCCLS                   IXCP"
  do xij = 2, nxij
    jxcijint(xij) = 0.25d0*jxcijint(xij)
    write(73, fmt='(i3,3x,f9.5,6x,d9.3,6x,3(1x,f7.4),i5)') xij,rxij(xij),dimag(jxcijint(xij)),rxccls(1:3,xij),ixcp(xij)
  enddo ! xij
  close(73, iostat=ios)
endsubroutine


!------------------------------------------------------------------------------
subroutine cmatmul(n,m,a,b,c)
  integer, intent(in) :: n, m
  double complex, intent(in) :: a(m,m), b(m,m)
  double complex, intent(out) :: c(m,m)
  !   ********************************************************************
  !   *   perform  the matrix-matrix operation           c = a * b       *
  !   *                                                                  *
  !   *   a,b,c   complex  square  n x n - matrices                      *
  !   *   n       dimension of a, b and c                                *
  !   *   m       array size of a, b, c with m >= n                      *
  !   ********************************************************************
  integer :: j, l

  c(1:n,1:n) = czero
  do j = 1, n
    do l = 1, n
      if (b(l,j) /= czero) c(1:n,j) = c(1:n,j) + a(1:n,l)*b(l,j)
    enddo ! l
  enddo ! j

endsubroutine

#if 0
c ************************************************************************
      SUBROUTINE CLSJIJ0(NAEZ, RR, NR, RBASIS, RCUT, JIJ, NRD, NXIJD)
c ************************************************************************
c This subroutine is used check the cluster size around each atom 
c where Jij's are calculated
c all arguments are input arguments
c
c called by: main0
c                                                          A.Thiess 7/2009
c ************************************************************************
      IMPLICIT NONE
c
c      INCLUDE 'inc.p'
c
c
c     .. array arguments
c
      DOUBLE PRECISION, INTENT(IN) :: RBASIS(3,NAEZ) ! positions in EZ
      DOUBLE PRECISION, INTENT(IN) :: RR(3,0:NRD) ! lattice vectors
c
c
c     .. scalar arguments
c
      INTEGER, INTENT(IN) :: NRD
      INTEGER, INTENT(IN) :: NXIJD
      DOUBLE PRECISION, INTENT(IN) :: RCUT
      INTEGER, INTENT(IN) :: NAEZ ! number of atoms in EZ
      INTEGER, INTENT(IN) :: NR   ! number of lattice vectors RR
      LOGICAL          JIJ
c
c     .. local arrays
c
      DOUBLE PRECISION TMP(3)
c
c     .. local scalars
c
      DOUBLE PRECISION EPSSHL,RCUT2,RTMP
      INTEGER          IAEZ,IR,NXIJ,I1
c
c
      DATA             EPSSHL   / 1.0D-4 /
c
c ------------------------------------------------------------------------
c This is generating the clusters which have a distance smaller
c than RCUT2.
c ------------------------------------------------------------------------
       IF (JIJ) THEN
C
       RCUT2 = (RCUT + EPSSHL)**2
C======================================================================
C loop in all atoms begin
C======================================================================
       DO I1 = 1, NAEZ ! loop over all Jij-centers
C
         NXIJ = 0           ! counter for atoms in cluster
         DO IAEZ = 1, NAEZ   ! loop in all atoms
           DO IR = 0, NR    ! loop in all bravais vectors    
             TMP(1:3) = RR(1:3,IR) + RBASIS(1:3,IAEZ) - RBASIS(1:3,I1)
             RTMP = TMP(1)**2 + TMP(2)**2 + TMP(3)**2
             IF (RTMP <= RCUT2)  THEN
               NXIJ = NXIJ + 1
               IF (NXIJ > NXIJD) THEN 
                 WRITE (6,*) 
     &           ' ERROR: Dimension NXIJD in inc.cls too small',
     &           NXIJ, NXIJD
                 STOP '   < CLSJIJ >'
               ENDIF
             ENDIF
           ENDDO              ! IR loop in bravais
         ENDDO                ! IAEZ loop in NAEZ
C
       ENDDO ! loop over all Jij-centers

      WRITE(6,'(79(1H=),/,15X,A)') 
     &     'CLSJIJ0: checking Jij-cluster size ........ OK'
      WRITE (6,'(79(1H=),/)')
C======================================================================
C======================================================================
       ENDIF
C
      END
#endif

endmodule
