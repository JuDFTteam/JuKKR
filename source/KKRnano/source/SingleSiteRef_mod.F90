module SingleSiteRef_mod
!-------------------------------------------------------------------------------
!> Summary: Construction of the single site reference system
!> Author: Rudolf Zeller, Elias Rabel, Alexander R Thiess, Paul F Baumeister, Marcel Bornemann
!> Category: KKRnano, single-site, solver, reference-system
!-------------------------------------------------------------------------------
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: gref, tref

  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0), ci=(0.d0,1.d0), cmone=(-1.d0,0.d0)

  contains
  
!-----------------------------------------------------------------------
!>    Calculate reference system's T-matrix.
!>    @param     e complex energy
!>    @param     vref repulsive reference potential field strength
!>    @param     lmax angular momentum cutoff
!>    @param     rMTref repulsive reference potential muffin-tin radius
!>    @param     trefLL reference system t-matrix
!>    @param     dtrefLL energy derivative of reference system t-matrix only calculated if Lly=1
!>    @param     Lly do lloyd's formula calculation 0=no/1=yes
  subroutine tref(e, vref, lmax, rMTref, trefLL, dtrefLL, derive)
    use SingleSiteHelpers_mod, only: bessel
    integer, intent(in) :: lmax
    double complex, intent(in) :: e
    double precision, intent(in) :: rMTref, vref
    double complex, intent(out) :: trefLL(0:lmax) ! storage scheme exploits that tref is m-degenerate and diagonal
    double complex, intent(out), optional :: dtrefLL(0:lmax)
    logical, intent(in) :: derive

    integer :: l
    double complex :: a1, b1, da1, db1, sqEmV, sqE
    double complex, dimension(0:lmax+1) :: Bjw1, Bjw2, Byw1, Byw2, Hws1, Hws2, dBjw1, dBjw2, dHws1

    !-----------------------------------------------------------------------
    !---- t-matrix and derivative of t-matrix of the reference system

    !     the analytical formula for the derivative of spherical Bel
    !     functions is used:

    !     d                     l+1
    !     --  j (x) = j   (x) - --- j (x)
    !     dx   l       l-1       x   l

    !     d
    !     --  j (x) = - j (x)
    !     dx   0         1

    !     which for x = sqrt(E0)*r leads to

    !      d          r*r             (l+1)
    !     --- j (x) = --- ( j   (x) - ----- j (x) )
    !     dE0  l      2 x    l-1        x    l

    !      d            r*r
    !     --- j (x) = - --- j (x)
    !     dE0  0        2 x  1

    !-----------------------------------------------------------------------

    sqE = sqrt(e)
    sqEmV = sqrt(e - vref)
    
    a1 = rMTref*sqE
    b1 = rMTref*sqEmV
    call bessel(Bjw1, Byw1, Hws1, a1, lmax+1)
    call bessel(Bjw2, Byw2, Hws2, b1, lmax+1)
    
    if (derive) then
    
      if (present(dtrefLL)) then
      
        dtrefLL = zero ! clear
        
        ! derivative of the J0(x/a) = -J1(x/a)/a 
        dBjw1(0) = -Bjw1(1)/a1
        dBjw2(0) = -Bjw2(1)/b1
        dHws1(0) = -Hws1(1)/a1

        do l = 1, lmax+1
          ! recursion formulae
          dBjw1(l) = (Bjw1(l-1) - (l+1)*Bjw1(l)/a1)/a1
          dBjw2(l) = (Bjw2(l-1) - (l+1)*Bjw2(l)/b1)/b1
          dHws1(l) = (Hws1(l-1) - (l+1)*Hws1(l)/a1)/a1
        enddo ! l

        ! scale
        dBjw1(:) = dBjw1(:)*0.5d0*rMTref**2
        dBjw2(:) = dBjw2(:)*0.5d0*rMTref**2
        dHws1(:) = dHws1(:)*0.5d0*rMTref**2

        do l = 0, lmax
          a1 = sqE*Bjw1(l+1)*Bjw2(l) - sqEmV*Bjw1(l)*Bjw2(l+1)

          da1 = 0.5d0/sqE*Bjw1(l+1)*Bjw2(l) - 0.5d0/sqEmV*Bjw1(l)*Bjw2(l+1) + sqE*dBjw1(l+1)*Bjw2(l) &
                    - sqEmV*dBjw1(l)*Bjw2(l+1) + sqE*Bjw1(l+1)*dBjw2(l) - sqEmV*Bjw1(l)*dBjw2(l+1)

          b1 = sqE*Hws1(l+1)*Bjw2(l) - sqEmV*Hws1(l)*Bjw2(l+1)

          db1 = 0.5d0/sqE*Hws1(l+1)*Bjw2(l) - 0.5d0/sqEmV*Hws1(l)*Bjw2(l+1) + sqE*dHws1(l+1)*Bjw2(l) &
                    - sqEmV*dHws1(l)*Bjw2(l+1) + sqE*Hws1(l+1)*dBjw2(l) - sqEmV*Hws1(l)*dBjw2(l+1)

          dtrefLL(l) = 0.5d0/sqE**3*a1/b1 - 1.d0/sqE*(da1/b1 - a1*db1/b1**2) ! add lm-diagonal, m-degenerate part
        enddo ! l
      else
        warn(6, "in the calculation of tref, the derivative was required but no array dtrefLL was passed!")
      endif ! present dtrefLL

    endif ! derive for Lloyd's formula
    
    ! regular reference t matrix
    do l = 0, lmax
      a1 = sqE*Bjw1(l+1)*Bjw2(l) - sqEmV*Bjw1(l)*Bjw2(l+1)
      b1 = sqE*Hws1(l+1)*Bjw2(l) - sqEmV*Hws1(l)*Bjw2(l+1)

      trefLL(l) = -1.d0/sqE*a1/b1 ! add lm-diagonal, m-degenerate part
    enddo ! l

  endsubroutine ! tref


!------------------------------------------------------------------------------
  subroutine gref(e, alat, iend, cleb, rcls, icleb, loflm, nacls, tref_ell, dtref_ell, grefn, Lly_g0tr, lmax, ncleb, Lly)
    integer, intent(in) :: lmax, nacls, Lly
    double complex, intent(in) :: e
    double precision, intent(in) :: alat
    integer, intent(in) :: iend
    double precision, intent(in) :: rcls(3,nacls)
    integer, intent(in) :: ncleb
    double precision, intent(in) :: cleb(ncleb,2)
    integer, intent(in) :: icleb(ncleb,3)
    integer, intent(in) :: loflm((2*lmax+1)**2)
    double complex, intent(in)  :: tref_ell(0:lmax,nacls), dtref_ell(0:lmax,nacls)
    double complex, intent(out) :: grefn(:,:,0:,:) ! grefn((lmax+1)**2,(lmax+1)**2,0:1,nacls)
    double complex, intent(out) :: Lly_g0tr

    integer :: iacls, lm, lmmaxd, ist
    double complex, allocatable :: ginp(:,:,:,:)
    
    lmmaxd = (lmax + 1)**2
    
    assert(all(shape(grefn) >= [lmmaxd, lmmaxd, 1+Lly, nacls])) ! prepare for deferred shape interface
    
    allocate(ginp(lmmaxd,nacls,lmmaxd,0:1))
    
    call gll95(e, ncleb, cleb(:,2), icleb, loflm, iend, tref_ell, dtref_ell, rcls, nacls, alat, ginp, Lly_g0tr, lmax, Lly)

    do iacls = 1, nacls
      do lm = 1, lmmaxd
        ! ginp has dim(lmmaxd,nacls,lmmaxd,0:1), so we interchange dimensions here
        grefn(:,lm,0,iacls) = ginp(:,iacls,lm,0) ! value
        if (Lly > 0) & 
        grefn(:,lm,Lly,iacls) = ginp(:,iacls,lm,1) ! energy derivative
      enddo ! lm
    enddo ! iacls
    
    deallocate(ginp, stat=ist)

  endsubroutine ! gref
  
  

!**********************************************************************
!> @param e complex energy
!> @param cleb array of gaunt coefficients
!> @param icleb index array for gaunt coefficients
!> @param loflm array that maps (lm)-index to l-index
!> @param iend ???
!> @param trefLL reference t-matrix
!> @param dtrefLL derivative of reference t-matrix
!> @param ratom real space positions of atoms in ref. cluster
!> @param nacls number of atoms in reference cluster
!> @param alat length of unit vector in bohr
!> @param gref0 todo
!> @param dgref0 ??? energy derivative of green's function
!> @param Lly_g0tr trace(m^-1 dm/de) with m = (1 - g0 \delta t_ref)
!> @param lmax angular momentum cutoff (it would be better to rewrite routine to pass lmmaxd)
!> @param ncleb number of gaunt coefficients in cleb
!> @param Lly do lloyd's formula calculations 1=yes/0=no
  subroutine gll95(e, ncleb, cleb, icleb, loflm, iend, tref_ell, dtref_ell, ratom, nacls, alat, gref0, Lly_g0tr, lmax, Lly)
  ! **********************************************************************
  !
  !     solution of the dyson equation for a cluster of potentials
  !     (trefLL) centered at positions ratom in free space,
  !
  ! ----------------------------------------------------------------------
    integer, intent(in) :: lmax, nacls, Lly
    double complex, intent(in) :: e
    double precision, intent(in) :: alat
    double complex, intent(in) :: tref_ell(0:,1:) ! (0:lmax,nacls) ! trefLL((lmax+1)**2,(lmax+1)**2,nacls)
    double complex, intent(in) :: dtref_ell(0:,1:) ! dtrefLL((lmax+1)**2,(lmax+1)**2,nacls)
    integer, intent(in) :: ncleb
    double precision, intent(in) :: cleb(ncleb)
    integer, intent(in) :: icleb(ncleb,3)
    integer, intent(in) :: loflm(:)
    integer, intent(in) :: iend
    double precision, intent(in) :: ratom(3,*) ! first dim: 3
    double complex, intent(out) :: Lly_g0tr
    double complex, intent(out) :: gref0((lmax+1)**2,nacls,(lmax+1)**2,0:Lly) !> dim(lmmaxd,nacls,lmmaxd,0:Lly)

    external :: zcopy, zgemm, zgetrs ! from BLAS
    double complex, allocatable :: gref(:,:,:,:,:), gtref(:,:,:), trefLL(:,:)
    double complex, allocatable :: dgtde0(:,:,:,:), dgtde(:,:,:) ! Lly
    integer, allocatable :: ipvt(:)
    integer :: lmmaxd, ngd, iacls, info, ist

    lmmaxd = (lmax + 1)**2
    
    assert(all(shape(gref0) == [lmmaxd, nacls, lmmaxd, 1+Lly]))
    
    allocate(gref(lmmaxd,nacls,lmmaxd,nacls,0:Lly), trefLL(lmmaxd,lmmaxd), stat=ist)
    if (ist /= 0) die_here("gll95: fatal error, failure to allocate memory, probably out of memory.")   

    call calcFreeGreens(gref(:,:,:,:,0), e, lmax, lmmaxd, nacls, ratom, alat, cleb, icleb, ncleb, iend, loflm, derive=.false.)

    ngd = lmmaxd*nacls ! combine the 1st and 2nd dimension to the assumed leading dimension for BLAS and LAPACK calls
    
    ! construct right hand side of linear equation system for grefsy
    ! the first lmmaxd columns of gref are copied into gref0
    ! gref0 then contains g0^{(1)n'}_{ll'}, the free space structural
    ! green's function for the central cluster atom (e.r.)
    ! --------------------------------------------------------------

    if (Lly > 0) then

      call calcFreeGreens(gref(:,:,:,:,Lly), e, lmax, lmmaxd, nacls, ratom, alat, cleb, icleb, ncleb, iend, loflm, derive=.true.)
      
      allocate(dgtde0(lmmaxd,nacls,lmmaxd,nacls), dgtde(lmmaxd,nacls,lmmaxd), stat=ist)
      if (ist /= 0) die_here("gll95: fatal error, failure to allocate memory, probably out of memory.")

      ! independent iacls, private(trefLL)
      do iacls = 1, nacls
        
        trefLL(:,:) = unfold_m_deg_diag_rep(lmax, tref_ell(0:,iacls))
        ! dgtde0(:,: , :,iacls) = - dg_0/de * \delta t_ref
        call zgemm('n','n',ngd,lmmaxd,lmmaxd,cmone, gref(1,1,1,iacls,Lly),ngd, trefLL,lmmaxd, zero, dgtde0(1,1,1,iacls),ngd)
        ! dgtde0(:,: , :,iacls) =  dgtde0(:,: , :,iacls)  - g_0 * d(\delta t_ref)/de
        !                       = -dg_0/de * \delta t_ref - g_0 * d(\delta t_ref)/de
#define dtrefLL trefLL
        dtrefLL(:,:) = unfold_m_deg_diag_rep(lmax, dtref_ell(0:,iacls))
        call zgemm('n','n',ngd,lmmaxd,lmmaxd,cmone, gref(1,1,1,iacls,0),ngd, dtrefLL,lmmaxd, cone, dgtde0(1,1,1,iacls),ngd)
#undef  dtrefLL

        if (iacls == 1) dgtde(:,:,:) = dgtde0(:,:,:,1) ! copy slice #1
      enddo ! iacls

    else
      allocate(dgtde(1,1,1), stat=ist)
    endif ! Lly > 0

    gref0(:,:,:,0) = gref(:,:,:,1,0) ! should be equivalent to call zcopy(ngd*lmmaxd,gref,1,gref0,1)

    allocate(gtref(lmmaxd,nacls,lmmaxd), ipvt(lmmaxd*nacls), stat=ist)
    if (ist /= 0) die_here("gll95: fatal error, failure to allocate memory, probably out of memory.")   
    do iacls = 1, nacls
      ! -g_ref * \delta t_ref  -- stored in gtref
      trefLL(:,:) = unfold_m_deg_diag_rep(lmax, tref_ell(0:,iacls))
      call zgemm('n','n',ngd,lmmaxd,lmmaxd,cmone, gref(1,1,1,iacls,0),ngd, trefLL,lmmaxd, zero, gtref,ngd)
      gref(:,:,:,iacls,0) = gtref(:,:,:) ! copy slice back into gref, here, the minus sign comes in
    enddo ! iacls
    deallocate(gtref, stat=ist)

    ! solve Dyson-equation for the reference system
    ! solves (1 - g0 \delta t) g_ref = g0 for g_ref.
    call grefsy(gref(:,:,:,:,0), gref0(:,:,:,0), ipvt, dgtde, Lly_g0tr, nacls, lmmaxd, Lly)

    if (Lly > 0) then
    
      call zgemm('n','n',ngd,lmmaxd,ngd,cmone, dgtde0,ngd, gref0(1,1,1,0),ngd, cone, gref(1,1,1,1,Lly),ngd)
      call zgetrs('n',ngd,lmmaxd, gref(1,1,1,1,0),ngd, ipvt, gref(1,1,1,1,Lly),ngd, info)
      gref0(:,:,:,Lly) = gref(:,:,:,1,Lly) ! copy slice #1
      
    endif ! Lly > 0

    deallocate(gref, dgtde, dgtde0, trefLL, stat=ist) ! ignore status
  endsubroutine ! gll95

  
  function unfold_m_deg_diag_rep(lmax, tref_ell) result(tref)
    integer, intent(in) :: lmax
    double complex, intent(in) :: tref_ell(0:) ! dim(0:lmax) if we passed an assumed shape array, an array temporary would be created at runtime 
    double complex :: tref((lmax+1)**2,(lmax+1)**2) ! result
    integer :: l, m, lm
    tref = zero
    do l = 0, lmax
      do m = -l, l
        lm = l*l + l + m + 1
        tref(lm,lm) = tref_ell(l) ! stored m-degenerate and diagonal
      enddo ! m
    enddo ! l
  endfunction ! unfold
  
  
  
 !> calculates the free-space green-function or the derivative of the free-space green-function
 !> @param[out] greenfree   the free space green-function
 !> @param      nacls       number of atoms in reference cluster
 !> @param[in]  derivative  .false. = calculate free-space-greens function
 !>                         .true.  = calculate derivative of free-space-greens function
  subroutine calcFreeGreens(greenfree, e, lmax, lmmaxd, nacls, ratom, alat, cleb, icleb, ncleb, iend, loflm, derive)
    double complex, intent(out) :: greenfree(:,:,:,:) !> dim(lmmaxd,nacls,lmmaxd,nacls)
    double precision, intent(in) :: alat
    double complex, intent(in) :: e
    integer, intent(in) :: ncleb
    integer, intent(in) :: icleb(ncleb,3)
    double precision, intent(in) :: cleb(:)
    integer, intent(in) :: iend
    integer, intent(in) :: lmax
    integer, intent(in) :: lmmaxd, nacls
    integer, intent(in) :: loflm(:)
    double precision, intent(in) :: ratom(3,*)
    logical, intent(in) :: derive

    ! local
    double precision :: rdiff(3)
    integer :: iacls, jacls
    double complex :: gll(lmmaxd,lmmaxd) ! automatic array
    
    assert(all(shape(greenfree) == [lmmaxd, nacls, lmmaxd, nacls]))
    
    !
    ! ---> construct free green's function
    ! the free space structural green's function g0 is a matrix of dimension lmmaxd x lmmaxd
    ! (for a certain pair of reference cluster atoms n and n')
    ! it is calculated in routine gfree and stored in gll
    !
    ! then the matrix gref^{nn'}_{ll'} is constructed
    ! the blocks n /= n' contain g0,
    ! the other n=n' blocks are set to zero (green's function is not defined for r=r')
    ! (e.r.)
    
    do iacls = 1, nacls
      do jacls = 1, nacls
      
        if (jacls /= iacls) then
          rdiff(1:3) = (ratom(1:3,iacls) - ratom(1:3,jacls))*alat ! difference vector

          call gfree(rdiff, e, gll, cleb, icleb, loflm, iend, lmax, ncleb, derive) ! fills gll with the Green function or the derivative of the Green function 
          greenfree(:,jacls,:,iacls) = gll(:,:)

        else  ! jacls /= iacls
        
          greenfree(:,jacls,:,iacls) = zero
        
        endif ! jacls /= iacls

      enddo ! jacls
    enddo ! iacls
    
  endsubroutine ! calcFreeGreens

  

!>    Solves (1 - g0 \Delta t) G_ref = g0 for G_ref (Full inversion).
!
!> @param  GTMAT     on input it has to contain (-1) * g0 * \Delta t
!>                   dimension (LMMAXD*NACLSD) x (LMMAXD*NACLSD)
!>                   on output: LU-factorisation of the input matrix
!> @param  GMAT  input: g0 free-space Green's function for the central
!>                      reference cluster atom: g0^{(1)N'}_{LL'}
!>               on output it contains G_ref
!> @param  IPVT         integer work array of dimension (LMMAXD*NACLSD)
!> @param  NDIM         NDIM = #cluster atoms * maximal LM
!>                      NDIM <= (LMMAXD*NACLSD)
!>
!>    \verbatim
!>                 NDIM
!>         +---------+----+
!>         |         |    |
!>         | contents|    |
!>         |         |    |      (GTMAT)
!>         |         |    |
!>    NDIM +---------+    |
!>         |    (uninit.) |
!>         +--------------+ (LMMAXD*NACLSD)
!
!>    \endverbatim
!
!> @param   NACLSD    MAXIMAL number of cluster atoms
!>
!>    @author: ???, commented by E. Rabel, Nov 2011

  subroutine grefsy(gtmat, gmat, ipvt, dgtde, Lly_g0tr, nacls, lmmaxd, Lly)
  !
  !---> solve the Dyson equation to get reference green function
  !
    integer, intent(in) :: nacls, lmmaxd, Lly ! lloyd's formula switch 0 (inactive)/ 1 (active)
    double complex, intent(inout) :: gtmat(:,:,:,:) ! (lmmaxd,nacls,lmmaxd,nacls) !> contains - g0 * \delta t on entry
    double complex, intent(inout) :: gmat(:,:,:) ! (lmmaxd,nacls,lmmaxd)
    double complex, intent(inout) :: dgtde(:,:,:) ! (lmmaxd,nacls,lmmaxd)
    integer,        intent(out) :: ipvt(:) ! (lmmaxd*nacls)
    double complex, intent(out) :: Lly_g0tr

    external :: zgetrf, zgetrs ! from LAPACK: L-U-factorisation and linear system solver
    integer :: info, ngd, lm, iacls

    assert(all(shape(gtmat) == [lmmaxd, nacls, lmmaxd, nacls]))
    assert(all(shape(gmat)  == [lmmaxd, nacls, lmmaxd]))

    do iacls = 1, nacls
      do lm = 1, lmmaxd
        gtmat(lm,iacls,lm,iacls) = gtmat(lm,iacls,lm,iacls) + cone ! add unity matrix
      enddo ! lm
    enddo ! iacls
    ! now gtmat == 1 - g0 * \delta t
    
    ngd = lmmaxd*nacls ! combine the first two dims of the arrays

    call zgetrf(ngd,ngd, gtmat,ngd, ipvt, info) ! L-U-factorisation
    
    call zgetrs('n',ngd,lmmaxd, gtmat,ngd, ipvt, gmat,ngd, info) ! solve linear system

    Lly_g0tr = zero
    if (Lly > 0) then ! Lloyd

      assert(all(shape(dgtde) == [lmmaxd, nacls, lmmaxd]))

      call zgetrs('n',ngd,lmmaxd, gtmat,ngd, ipvt, dgtde,ngd, info) ! solve linear system

      do lm = 1, lmmaxd
        Lly_g0tr = Lly_g0tr - dgtde(lm,1,lm)
      enddo ! lm

    endif ! Lly > 0 ! Lloyd

  endsubroutine ! grefsy

  
  subroutine gfree(rdiff, e, gmLL, cleb, icleb, loflm, iend, lmax, ncleb, derive)
    use Constants_mod, only: Pi
    use SingleSiteHelpers_mod, only: BesHan
    use Harmonics_mod, only: Ymy

    double precision, intent(in) :: rdiff(3)
    integer, intent(in) :: lmax, ncleb, iend
    double complex, intent(in) :: e
    double complex, intent(inout) :: gmLL((lmax+1)**2,(lmax+1)**2)
    double precision, intent(in) :: cleb(ncleb)
    integer, intent(in) :: icleb(ncleb,3), loflm(*)
    logical, intent(in) :: derive
    
    double precision :: fpi, rfpi
    integer :: ifac, j, lm1, lm2, lm3, lp1
    double complex :: bl(2*lmax+1), hl(2*lmax+1), nl(2*lmax+1), dhl(2*lmax+1) ! only used for the derivative
    double complex :: hyl((2*lmax+1)**2), sqE
    double precision :: yl((2*lmax+1)**2), rabs
    integer :: lmmax, lmx2sq

    lmmax = (lmax + 1)**2
    lmx2sq = (2*lmax + 1)**2

    fpi = 4.d0*pi
    rfpi = sqrt(fpi)
    sqE = sqrt(e)
!
!---- calculation of free electron green's function :  g(m)ll'(e)
!
    call Ymy(rdiff(1), rdiff(2), rdiff(3), rabs, yl, 2*lmax)
    call BesHan(hl, bl, nl, sqE*rabs, 2*lmax)

    if (derive) then
!-----------------------------------------------------------------------
!---- derivative of free electron green function matrix elements
!
!     the analytical formula for the derivative of spherical hankel
!     functions is used:
!
!     d                     l+1
!     --  h (x) = h   (x) - --- h (x)   
!     dx   l       l-1       x   l
!
!     which for x = sqrt(e)*r leads to
!
!     d                       r           rl
!     --- ( sqrt(e) h (x) ) = - h   (x) - -- h (x) )
!     de             l        2  l-1      2x  l
!
!-----------------------------------------------------------------------
    
      dhl(1) = 0.5d0*ci*rabs*hl(1) ! start recursion here
      do lp1 = 2, 2*lmax + 1
        dhl(lp1) = 0.5d0*(rabs*hl(lp1-1) - (lp1-1)*hl(lp1)/sqE)
      enddo ! lp1

      do lm1 = 1, lmx2sq
        hyl(lm1) = -fpi*ci*yl(lm1)*dhl(loflm(lm1)+1)
      enddo ! lm1
      
    else  ! derive
    
      do lm1 = 1, lmx2sq
        hyl(lm1) = -fpi*ci*sqE*yl(lm1)*hl(loflm(lm1)+1)
      enddo ! lm1

    endif ! derive
    
    do lm1 = 1, lmmax
      gmLL(lm1,lm1) = hyl(1)/rfpi
      do lm2 = 1, lm1-1
        gmLL(lm1,lm2) = zero
      enddo ! lm2
    enddo ! lm1

    do j = 1, iend
      lm1 = icleb(j,1)
      lm2 = icleb(j,2)
      lm3 = icleb(j,3)
      gmLL(lm1,lm2) = gmLL(lm1,lm2) + cleb(j)*hyl(lm3)
    enddo ! j

    do lm1 = 1, lmmax
      do lm2 = 1, lm1 - 1
        ifac = (-1)**(loflm(lm1) + loflm(lm2))
        gmLL(lm2,lm1) = ifac*gmLL(lm1,lm2)
      enddo ! lm2
    enddo ! lm1
    
  endsubroutine ! gfree

endmodule ! SingleSiteRef_mod
