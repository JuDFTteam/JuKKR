module SingleSiteRef_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: gll95, gref, tref

  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0), ci=(0.d0,1.d0)

  contains

!**********************************************************************
!> @param e complex energy
!> @param cleb array of gaunt coefficients
!> @param icleb index array for gaunt coefficients
!> @param loflm array that maps (lm)-index to l-index
!> @param iend ???
!> @param trefLL reference t-matrix
!> @param dtrefLL derivative of reference t-matrix
!> @param ratom real space positions of atoms in ref. cluster
!> @param natom number of atoms in reference cluster
!> @param alat length of unit vector in bohr
!> @param gref0 todo
!> @param dgdeout ??? energy derivative of green's function
!> @param Lly_g0tr trace(m^-1 dm/de) with m = (1 - g0 \delta t_ref)
!> @param lmaxd angular momentum cutoff (it would be better to rewrite routine to pass lmmaxd)
!> @param naclsd dimension array: maximal number of atoms in reference clusters
!> @param ncleb number of gaunt coefficients in cleb
!> @param Lly do lloyd's formula calculations 1=yes/0=no
  subroutine gll95(e,cleb,icleb,loflm,iend,trefLL,dtrefLL, ratom,natom,alat,gref0,dgdeout, Lly_g0tr, lmaxd, naclsd, ncleb, Lly)
  ! **********************************************************************
  !
  !     solution of the dyson equation for a cluster of potentials
  !     (trefLL) centered at positions ratom in free space,
  !
  ! ----------------------------------------------------------------------
    integer, intent(in) :: lmaxd, naclsd, ncleb, Lly
    double complex, intent(in) :: e
    double complex, intent(out) :: Lly_g0tr
    double precision, intent(in) :: alat
    integer, intent(in) :: natom
    double complex, intent(inout) :: gref0(naclsd*(lmaxd+1)**2,(lmaxd+1)**2)
    double complex, intent(in) :: trefLL((lmaxd+1)**2,(lmaxd+1)**2, naclsd)
    double complex, intent(in) :: dtrefLL((lmaxd+1)**2,(lmaxd+1)**2, naclsd)
    double complex, intent(out) :: dgdeout(Lly*(naclsd*(lmaxd+1)**2-1)+1,(lmaxd+1)**2)
    double precision, intent(in) :: cleb(ncleb)
    integer, intent(in) :: icleb(ncleb,3)
    integer, intent(in) :: loflm(:)
    integer, intent(in) :: iend
    double precision, intent(in) :: ratom(3,*)  ! first dim: 3

    external :: zcopy, zgemm ! from BLAS
    integer :: info
    integer :: ipvt(naclsd*(lmaxd+1)**2)
    integer :: n2, ndim, site_lm_index2
    integer :: memory_stat, memory_fail
  ! ---------------------------------------------------------------------
  !     the following arrays can be very large (> 60 mb in typical cases)
  !     therefore they are allocated on the heap
  ! ---------------------------------------------------------------------
    double complex, allocatable :: gref(:,:), gtref(:,:), dgtde(:,:)
    double complex, allocatable :: dgtde0(:,:), dgde(:,:)
    integer :: lmgf0d, lmmaxd, ngd

    lmmaxd = (lmaxd+1)**2
    lmgf0d = lmmaxd
    ngd = lmmaxd*naclsd

    !     allocate arrays
    memory_stat = 0
    memory_fail = .false.

    allocate(gref(ngd,ngd), stat=memory_stat) ; memory_fail = memory_fail + memory_stat**2
    allocate(gtref(ngd,lmmaxd), stat=memory_stat) ; memory_fail = memory_fail + memory_stat**2
    allocate(dgtde(Lly*(lmmaxd*naclsd-1)+1,lmmaxd), stat=memory_stat) ; memory_fail = memory_fail + memory_stat**2  !fix: workaround, need to pass it to grefsy

    if (Lly == 1) then
      !allocate(dgtde(ngd,lmmaxd), stat=memory_stat) ; memory_fail = memory_fail + memory_stat**2 ! wrong: need to allocate in any case because it has to be passed to grefsy
      allocate(dgtde0(ngd,ngd), stat=memory_stat) ; memory_fail = memory_fail + memory_stat**2
      allocate(dgde(ngd,ngd), stat=memory_stat) ; memory_fail = memory_fail + memory_stat**2
    end if

    if (memory_fail /= 0) die_here("gll95: fatal error, failure to allocate memory, probably out of memory.")   

    ndim = lmmaxd*natom ! ndim can be smaller than ngd=lmmaxd*naclsd
    call calcFreeGreens(gref, e, lmmaxd, natom, ratom, alat, cleb, icleb, ncleb, iend, loflm, derivative=.false.)

    if (Lly == 1) then
      call calcFreeGreens(dgde, e, lmmaxd, natom, ratom, alat, cleb, icleb, ncleb, iend, loflm, derivative=.true.)
    endif ! Lly == 1

  ! construct right hand side of linear equation system for grefsy
  ! the first lmmaxd columns of gref are copied into gref0
  ! gref0 then contains g0^{(1)n'}_{ll'}, the free space structural
  ! green's function for the central cluster atom (e.r.)
  ! --------------------------------------------------------------
    call zcopy(ngd*lmmaxd,gref,1,gref0,1)
  ! --------------------------------------------------------------

    if (Lly == 1) then
      do n2 = 1, natom
        site_lm_index2 = (n2-1)*lmmaxd + 1
        ! -dg_0/de * \delta t_ref    -- stored in gtref
        call zgemm('n','n',ndim,lmmaxd,lmmaxd,-cone,dgde(1,site_lm_index2),ngd,trefLL(1,1,n2), lmmaxd,zero,gtref,ngd)
        !   - g_0 * d(\delta t_ref)/de + gtref  -- stored again in gtref
        ! = -dg_0/de * \delta t_ref - g_0 * d(\delta t_ref)/de
        call zgemm('n','n',ndim,lmmaxd,lmmaxd,-cone,gref(1,site_lm_index2),ngd,dtrefLL(1,1,n2), lmmaxd,cone,gtref,ngd)
        ! copy gtref to dgtde0 - gtref is reused
        call zcopy(ngd*lmmaxd,gtref,1,dgtde0(1,site_lm_index2),1)
      enddo ! n2
    endif ! (Lly==1)

    do n2 = 1, natom
      site_lm_index2 = (n2-1)*lmmaxd + 1
      ! -g_ref * \delta t_ref  -- stored in gtref
      call zgemm('n','n',ndim,lmmaxd,lmmaxd,-cone,gref(1,site_lm_index2),ngd,trefLL(1,1,n2), lmmaxd,zero,gtref,ngd)
      call zcopy(ngd*lmmaxd,gtref,1,gref(1,site_lm_index2),1)
    enddo ! n2

    if (Lly == 1) then
      dgtde(1:ngd,1:lmmaxd) = dgtde0(1:ngd,1:lmmaxd)
    endif ! Lly == 1
      
    ! solve dyson-equation for reference system
    ! solves (1 - g0 \delta t) g_ref = g0 for g_ref.
    Lly_g0tr = zero
    call grefsy(gref,gref0,ipvt,ndim,dgtde,Lly_g0tr,naclsd, lmmaxd, Lly)

    if (Lly == 1) then
      call zgemm('n','n',ndim,lmmaxd,ndim,-cone,dgtde0,ngd,gref0,ngd,cone,dgde,ngd)
      call zgetrs('n',ndim,lmmaxd,gref,ngd,ipvt,dgde,ngd,info)
      dgdeout(1:ngd,1:lmmaxd) = dgde(1:ngd,1:lmmaxd)
    endif ! Lly == 1

    deallocate(gref, gtref, dgtde)
    if (Lly == 1) deallocate(dgtde0, dgde)

  endsubroutine gll95

 !> calculates the free-space green-function or the derivative of the free-space green-function
 !> @param[out] greenfree   the free space green-function
 !> @param      natom       number of atoms in reference cluster
 !> @param[in]  derivative  .false. = calculate free-space-greens function
 !>                         .true.  = calculate derivative of free-space-greens function
  subroutine calcFreeGreens(greenfree, energy, lmmaxd, natom, ratom, alat, cleb, icleb, ncleb, iend, loflm, derivative)
    use kkr_helpers_mod, only: lmmaxtolmax
    double complex, intent(out) :: greenfree(:,:)
    integer, intent(in) :: ncleb
    double precision, intent(in) :: alat
    double precision, intent(in) :: cleb(:)
    double complex, intent(in) :: energy
    integer, intent(in) :: icleb(ncleb,3)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: loflm(:)
    integer, intent(in) :: natom, iend
    double precision, intent(in) :: ratom(3,*)
    logical, intent(in) :: derivative

    ! local
    double precision :: rdiff(3)
    integer :: lm1, lm2, lmaxd
    double complex :: gll(lmmaxd,lmmaxd) ! automatic array
    integer :: n1, n2, site_lm_index1, site_lm_index2

    lmaxd = lmmaxtolmax(lmmaxd)
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
    do n1 = 1, natom
      do n2 = 1, natom
        !            rdiff(i) = (ratom(i,n1) - ratom(i,n2))*alat
        !           changed p.z. 4.7.97
        rdiff(1:3) = (ratom(1:3,n2) - ratom(1:3,n1))*alat

        if (n1 /= n2) then

          call gfree(rdiff,energy,gll,cleb,icleb,loflm,iend, lmaxd, ncleb, derivative) ! fills gll with the Green function or the derivative of the Green function 
          
        else  ! n1 /= n2
        
          gll(:,:) = zero
        
        endif ! n1 /= n2

        do lm2 = 1, lmmaxd
          site_lm_index2 = (n2-1)*lmmaxd + lm2
          do lm1 = 1, lmmaxd
            site_lm_index1 = (n1-1)*lmmaxd + lm1
            
            greenfree(site_lm_index1,site_lm_index2) = gll(lm1,lm2)
            
          enddo ! lm1
        enddo ! lm2
        
      enddo ! n2
    enddo ! n1
    
  endsubroutine calcFreeGreens

  
  
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
    double complex, intent(out) ::  trefLL((lmax+1)**2,(lmax+1)**2)
    double complex, intent(out), optional :: dtrefLL((lmax+1)**2,(lmax+1)**2)
    logical, intent(in) :: derive

    integer :: m, l, lm
    double complex :: a1, b1, da1, db1, sqEmV, sqE
    double complex, dimension(0:lmax+1) :: Bjw1, Bjw2, Byw1, Byw2, Hws1, Hws2, dBjw1, dBjw2, dHws1
    integer :: lmaxd, lmgf0d

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

    lmaxd = lmax
    lmgf0d = (lmaxd+1)**2

    sqE = sqrt(e)
    sqEmV = sqrt(e - vref)
    
    a1 = rMTref*sqE
    b1 = rMTref*sqEmV
    call bessel(Bjw1, Byw1, Hws1, a1, lmaxd+1)
    call bessel(Bjw2, Byw2, Hws2, b1, lmaxd+1)

    ! trefLL could be initalized to zero here !!! PFB
    
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

          do m = -l, l
            lm = l*l + l + m + 1
            dtrefLL(lm,lm) = 0.5d0/sqE**3*a1/b1 - 1.d0/sqE*(da1/b1 - a1*db1/b1**2) ! add lm-diagonal, m-degenerate part
          enddo ! m

        enddo ! l
      else
        warn(6, "in the calculation of tref, the derivative was required but no array dtrefLL was passed!")
      endif ! present dtrefLL

    endif ! derive for Lloyd's formula
    
    ! regular reference t matrix
    do l = 0, lmax
      a1 = sqE*Bjw1(l+1)*Bjw2(l) - sqEmV*Bjw1(l)*Bjw2(l+1)
      b1 = sqE*Hws1(l+1)*Bjw2(l) - sqEmV*Hws1(l)*Bjw2(l+1)

      do m = -l, l
        lm = l*l + l + m + 1
        trefLL(lm,lm) = -1.d0/sqE*a1/b1 ! add lm-diagonal, m-degenerate part
      enddo ! m

    enddo ! l

  endsubroutine ! tref


!------------------------------------------------------------------------------
  subroutine gref(e,alatc,iend, cleb,rcls,icleb,loflm,nacls, trefLL,dtrefLL,grefn,dgrefn, Lly_g0tr, lmaxd, naclsd, ncleb, Lly)
    integer, intent(in) :: lmaxd, naclsd, ncleb, Lly
    double complex, intent(in) :: e
    double precision, intent(in) :: alatc
    integer, intent(in) :: iend
    double precision, intent(in) :: cleb(ncleb,2), rcls(3,naclsd)
    integer, intent(in) :: icleb(ncleb,3)
    integer, intent(in) :: loflm((2*lmaxd+1)**2)
    integer, intent(in) :: nacls
    double complex, intent(out) :: Lly_g0tr
    double complex, intent(in) :: trefLL((lmaxd+1)**2,(lmaxd+1)**2, naclsd), dtrefLL((lmaxd+1)**2,(lmaxd+1)**2, naclsd)
    double complex, intent(out) :: dgrefn((lmaxd+1)**2,(lmaxd+1)**2,naclsd), grefn((lmaxd+1)**2,(lmaxd+1)**2,naclsd)
    
    integer :: ig, ig1, lm, lm2, lmgf0d, lmmaxd
    double complex ::  ginp(naclsd*(lmaxd+1)**2,(lmaxd+1)**2)
    double complex :: dginp(naclsd*(lmaxd+1)**2,(lmaxd+1)**2)

    lmmaxd = (lmaxd+1)**2
    lmgf0d = (lmaxd+1)**2

    ! attention in this subroutine i3 labels the fixed atom - i1 is a variable !

    call gll95(e,cleb(1,2),icleb,loflm,iend,trefLL,dtrefLL, rcls,nacls,alatc,ginp,dginp, Lly_g0tr, lmaxd, naclsd, ncleb, Lly)

    do ig = 1, naclsd
      do lm = 1, lmgf0d
        do lm2 = 1, lmgf0d
          ig1 = (ig-1)*lmgf0d+lm2
           grefn(lm2,lm,ig) =  ginp(ig1,lm)
          dgrefn(lm2,lm,ig) = dginp(ig1,lm)
        enddo ! lm2
      enddo ! lm
    enddo ! ig

  endsubroutine ! gref
  

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

  subroutine grefsy(gtmat,gmat,ipvt,ndim,dgtde, Lly_g0tr, naclsd, lmmaxd, Lly)
!
!---> solve the dyson equation to get reference green function
!
    integer, intent(in) :: ndim, naclsd, lmmaxd, Lly ! lloyd's formula switch 0 (inactive)/ 1 (active)
    double complex, intent(out) :: Lly_g0tr
    double complex, intent(inout) :: gmat(lmmaxd*naclsd,lmmaxd)
    double complex, intent(inout) :: gtmat(lmmaxd*naclsd,lmmaxd*naclsd)
    double complex, intent(inout) :: dgtde(Lly*(lmmaxd*naclsd-1)+1,lmmaxd)
    
    external :: zgetrf, zgetrs ! from BLAS
    integer :: ipvt(lmmaxd*naclsd)
    integer :: i, info, lmgf0d, ngd

    lmgf0d = lmmaxd
    ngd = lmmaxd*naclsd

    do i = 1, ndim
      gtmat(i,i) = cone + gtmat(i,i) ! gtmat= 1 - g0 * \delta t
    enddo ! i
!
!---> solve the system of linear equations
!
    call zgetrf(ndim,ndim,gtmat,ngd,ipvt,info) ! L-U-factorisation
    call zgetrs('n',ndim,lmgf0d,gtmat,ngd,ipvt,gmat,ngd,info)

! ..  lloyd
    if (Lly == 1) then

      call zgetrs('n',ndim,lmgf0d,gtmat,ngd,ipvt,dgtde,ngd,info)

      Lly_g0tr = zero

      do i = 1, lmgf0d
        Lly_g0tr = Lly_g0tr - dgtde(i,i)
      enddo ! i

    else  ! Lly == 1

      Lly_g0tr = zero

    endif ! Lly == 1
! .. lloyd

  endsubroutine ! grefsy
  
  
  
  subroutine gfree(rdiff,e0,gmll,cleb,icleb,loflm,iend,lmax, ncleb, derivative)
    use SingleSiteHelpers_mod, only: beshan
    use Harmonics_mod, only: ymy
    
    integer, intent(in) :: lmax, ncleb, iend
    double complex, intent(in) :: e0
    double complex, intent(inout) :: gmll((lmax+1)**2,(lmax+1)**2)
    double precision, intent(in) :: cleb(ncleb), rdiff(3)
    integer, intent(in) :: icleb(ncleb,3), loflm(*)
    logical, intent(in) :: derivative
    
    double precision fpi,pi,rfpi
    integer ifac,j,lm1,lm2,lm3,lp1
    double complex bl(lmax*2+1),hl(lmax*2+1)
    double complex hyl((lmax*2+1)**2),nl(lmax*2+1)
    double complex :: dhl(lmax*2+1) ! only used for the derivative
    double precision :: yl((lmax*2+1)**2), rabs
    integer :: lmgf0d, lmx2sq

    lmgf0d = (lmax+1)**2
    lmx2sq = (lmax*2+1)**2

    pi = 4.d0*atan(1.d0)
    fpi = 4.d0*pi
    rfpi = sqrt(fpi)
!
!---- calculation of free electron green's function :  g(m)ll'(e0)
!
    call ymy(rdiff(1),rdiff(2),rdiff(3),rabs,yl,lmax*2)
    call beshan(hl,bl,nl,sqrt(e0)*rabs,lmax*2)
    
    if (derivative) then
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
!     which for x = sqrt(e0)*r leads to
!
!      d                       r           rl
!     --- ( sqrt(e0) h (x) ) = - h   (x) - -- h (x) )
!     de0             l        2  l-1      2x  l
!
!-----------------------------------------------------------------------
    
      dhl(1) = 0.5d0*ci*rabs*hl(1) ! start recursion here
      do lp1 = 2, lmax*2+1
        dhl(lp1) = 0.5d0*(rabs*hl(lp1-1)-(lp1-1)*hl(lp1)/sqrt(e0))
      enddo ! lp1
      
      do lm1 = 1, lmx2sq
        hyl(lm1) = -fpi*ci*yl(lm1)*dhl(loflm(lm1)+1)
      enddo ! lm1
      
    else  ! derivative
    
      do lm1 = 1, lmx2sq
        hyl(lm1) = -fpi*ci*sqrt(e0)*yl(lm1)*hl(loflm(lm1)+1)
      enddo ! lm1
    
    endif ! derivative
    
    do lm1 = 1, lmgf0d
      gmll(lm1,lm1) = hyl(1)/rfpi
      do lm2 = 1, lm1-1
        gmll(lm1,lm2) = zero
      enddo ! lm2
    enddo ! lm1
    
    do j = 1, iend
      lm1 = icleb(j,1)
      lm2 = icleb(j,2)
      lm3 = icleb(j,3)
      gmll(lm1,lm2) = gmll(lm1,lm2) + cleb(j)*hyl(lm3)
    enddo ! j
  
    do lm1 = 1, lmgf0d
      do lm2 = 1, lm1-1
        ifac = (-1)**(loflm(lm1)+loflm(lm2))
        gmll(lm2,lm1) = ifac*gmll(lm1,lm2)
      enddo ! lm2
    enddo ! lm1
    
  endsubroutine gfree

endmodule SingleSiteRef_mod
