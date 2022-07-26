      MODULE MOD_RLLSLL
        CONTAINS
      SUBROUTINE RLLSLL(RPANBOUND,RMESH,VLL,RLL,SLL,TLLP, &
                        NCHEB,NPAN,LMSIZE,LMSIZE2,NRMAX, &
                        nvec,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR, &
                        cmoderll,cmodesll,cmodetest,idotime)
! ************************************************************************
! radial wave functions by the integral equation method of
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! which has been extended for KKR using non-sperical potentials.
! Further information can be found in 
!
! David Bauer, 
! "Development of a relativistic full-potential first-principles multiple scattering 
! Green function method applied to complex magnetic textures of nano structures 
! at surfaces", PhD Thesis, 2014
!
! http://darwin.bth.rwth-aachen.de/opus3/volltexte/2014/4925/
!
!
!
! ************************************************************************
! This routine solves the following two equations:
!
! ULL(r) = J(r) - PRE * J(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * ULL(r') ) 
!               + PRE * H(r) * int_0^r( dr' r'^2 J2(r') * op(V(r')) * ULL(r') )
!
! SLL(r) = H(r) - PRE * H(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * RLL(r') ) 
!               + PRE * J(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * SLL(r') )
!
! where the integral int_0^r() runs from 0 to r
! ************************************************************************
! Potential matrix : VLL(LMSIZE*NVEC,LMSIZE*NVEC)
! LMSIZE = LMMAX (number of LM components) x Number of spin components
! LMSIZE2 = NVEC* LMSIZE 
! NVEC is 2 for a spinor and 1 in case of a non-rel. calculation
! 
! ************************************************************************
! Green function prefacor PRE=GMATPREFACTOR (scalar value)
! tipically \kappa for non-relativistic and M_0 \kappa for SRA 
! 
! ************************************************************************


! ************************************************************************
! The discretization of the Lippmann-Schwinger equation results in a matrix
! equation which is solved in this routine. Further information is given
! in section 5.2.3, page 90 of Bauer, PhD 
!
! Source terms : 
!   right solution:  J, H  (nvec*lmsize,lmsize) or (lmsize,nvec*lmsize)
!    left solution:  J2,H2 (lmsize,nvec*lmsize) or (nvec*lmsize,lmsize)
!
! Example:
! The source term J is for LMSIZE=3 and NVEC=2 given by:
! J =      / jlk(jlk_index(1))                                          \
!          |       0            jlk(jlk_index(2))                       |
!          |       0                   0            jlk(jlk_index(3))   |
!          | jlk(jlk_index(4))                                          |
!          |       0            jlk(jlk_index(5))                       |
!          \       0                   0            jlk(jlk_index(6))   /
!
! first 3 rows are for the large and the last 3 rows for the small component
! ************************************************************************
! Operator op() can be chosen to be a unity or a transpose operation
!     The unity operation is used to calculate the right solution
!     The transpose operation is used to calculate the left solution
! ************************************************************************
! RMESH      - radial mesh
! RPANBOUND  - panel bounds RPANBOUND(0) left  panel border of panel 1
!                           RPANBOUND(1) right panel border of panel 1
! NCHEB      - highes chebyshev polynomial
!              number of points per panel = NCHEB + 1
! NPAN       - number of panels
! LMSIZE     - number of colums for the source matrix J etc...
! LMSIZE2    - number of rows   for the source matrix J etc...
! NRMAX      - total number of radial points (NPAN*(NCHEB+1))
! NVEC       - number of LMSIZE*LMSIZE blocks in J (LMSIZE2=NVEC*LMSIZE)
! ************************************************************************
use mod_timing                            ! timing routine
use mod_beshank                           ! calculates bessel and hankel func.
use mod_chebint                           ! chebyshev integration routines
use mod_config, only: config_testflag     ! reads if testflags are present
use mod_rllslltools, only: inverse        ! inverse matrix routine
use mod_physic_params,only: cvlight       ! speed of light
use sourceterms                           
use mod_chebyshev
implicit none
      integer :: ncheb                               ! number of chebyshev nodes
      integer :: npan                                ! number of panels
      integer :: lmsize                              ! lm-components * nspin 
      integer :: lmsize2                             ! lmsize * nvec
      integer :: nvec                                ! spinor integer
                                                     ! nvec=1 non-rel, nvec=2 for sra and dirac
      integer :: nrmax                               ! total number of rad. mesh points

      double complex,parameter:: ci= (0.0d0,1.0d0), &! complex i
                                 cone=(1.0d0,0.0d0),&!         1
                                 czero=(0.0d0,0.0d0) !         0
      ! running indices
      integer ivec, ivec2                            
      integer l1,l2,lm1,lm2,lm3
      integer info,icheb3,icheb2,icheb,icheb1,ipan,mn,nm,mn2,nplm

      ! source terms
      double complex :: gmatprefactor               ! prefactor of green function
                                                    ! non-rel: = kappa = sqrt e
      double complex :: hlk(:,:), jlk(:,:), &       ! right sol. source terms
                        hlk2(:,:), jlk2(:,:)        ! left sol. source terms
                                                    ! (tipically bessel and hankel fn)

      integer jlk_index(:)                          ! mapping array l = jlk_index(lm)
                                                    ! in: lm-index
                                                    ! corresponding l-index used hlk,..
                                                    ! hlk(l) = jlk_index(lm)

      character(len=1) :: cmoderll,cmodesll,cmodetest  ! These define the op(V(r)) in the eqs. above
                                                       ! (comment in the beginning of this subroutine)
                                                       ! cmoderll ="1" : op( )=identity       for reg. solution
                                                       ! cmoderll ="T" : op( )=transpose in L for reg. solution
                                                       ! cmodesll: same for irregular

      double complex ::  sll(lmsize2,lmsize,nrmax), &  ! irr. volterra sol.
                         rll(lmsize2,lmsize,nrmax), &  ! reg. fredholm sol.
                         tllp(lmsize,lmsize), &        ! t-matrix
                         vll(lmsize*nvec,lmsize*nvec,nrmax) ! potential term in 5.7 
                                                       ! bauer, phd
      double complex,allocatable ::  ull(:,:,:)        ! reg. volterra sol.

      double complex,allocatable ::  &
                     work(:,:), &
                     work2(:,:), &
                     allp(:,:,:),bllp(:,:,:), &                  ! eq. 5.9, 5.10 for reg. sol
                     cllp(:,:,:),dllp(:,:,:), &                  ! same for the irr. sol
                     slv(:,:,:,:),srv(:,:,:,:), &                ! a in eq 5.68
                     slv1(:,:,:,:),srv1(:,:,:,:), &              !****************
                     slv2(:,:,:,:),srv2(:,:,:,:), &              ! used for sra trick
                     slv3(:,:,:,:),srv3(:,:,:,:), &              !****************
                     mrnvy(:,:,:),mrnvz(:,:,:), &                ! ***************
                     mrjvy(:,:,:),mrjvz(:,:,:), &                !    eq. 5.19-5.22
                     mihvy(:,:,:),mihvz(:,:,:), &                !
                     mijvy(:,:,:),mijvz(:,:,:), &                ! ***************
                     yill(:,:,:),zill(:,:,:), &                  ! source terms  (i:irreg., r: regular)
                     yrll(:,:,:),zrll(:,:,:),yrlltmp(:,:,:), &   ! 
                     yill1(:,:,:),zill1(:,:,:), &                ! source terms in case of sratrick
                     yrll1(:,:,:),zrll1(:,:,:), &
                     yill2(:,:,:),zill2(:,:,:), &
                     yrll2(:,:,:),zrll2(:,:,:), &
                     vhlr(:,:,:), &                               ! vhlr = h * v (regular sol.) 
                     vjlr(:,:,:), &                               ! vjlr = j * v (regular sol.)
                     vhli(:,:,:), &                               ! vhli = h * v (irregular sol.)
                     vjli(:,:,:)                                  ! vjli = j * v (irregular sol.)
      double complex,allocatable :: yif(:,:,:,:), &               ! source terms (different array
                     yrf(:,:,:,:), &                              !               ordering)
                     zif(:,:,:,:), &
                     zrf(:,:,:,:)
      ! chebyshev arrays
      double complex zslc1sum(0:ncheb)
      double precision c1(0:ncheb,0:ncheb),rpanbound(0:npan)
      double precision cslc1(0:ncheb,0:ncheb), & ! Integration matrix from left ( C*S_L*C^-1 in eq. 5.53)
                       csrc1(0:ncheb,0:ncheb), & ! Same from right ( C*S_R*C^-1 in eq. 5.54)
                       tau(0:ncheb,0:npan), &    ! Radial mesh point
                       slc1sum(0:ncheb),rmesh(nrmax)

      integer ipiv(0:ncheb,lmsize2)
      integer,allocatable :: ipiv2(:)
      logical test
      integer :: ierror,use_sratrick
      integer :: idotime
      integer,parameter  :: directsolv=1

      external zgetrf,zgetrs
      intrinsic abs,atan,cos,dimag,exp,max,min,sin,sqrt

! ***********************************************************************
!                                  SRA trick
! ***********************************************************************
! on page 68 of Bauer, PhD, a method is described how to speed up the 
! calculations in case of the SRA. A similar approach is implemented 
! here by using Eq. 4.132 and substituting DV from 4.133, and discretising
! the radial mesh of the Lippmann-Schwinger eq. according to 5.68. 
! The Lippmann-Schwinger equation leads to a matrix inversion 
! problem. The matrix M which needs to be inverted has a special form
! if the SRA approximation is used:
! 
! matrix A ( C 1)     (same as in eq. 5.68)
!          ( B 0)
! (C, B are matricies here)
!
! inverse of A is   (C^-1    0 )
!                   (-B C^-1 1 )
! Thus, it is sufficient to only inverse the matrix C which saves computational
! time. This is refered to as the SRA trick.
! ***********************************************************************
! in future implementation equation 4.134 is supposed to be 
! implemented which should lead to an additional speed-up.
! ***********************************************************************

if ( .not. config_testflag('sph') .or. lmsize==1 ) then
  use_sratrick=0
elseif ( config_testflag('sph') ) then
  use_sratrick=1
else
  stop '[rllsll] use_sratrick error'
end if

if (idotime==1) call timing_start('rllsll')


allocate ( ull(lmsize2,lmsize,nrmax) )

if ( use_sratrick==0 ) then
  allocate (slv(0:ncheb,lmsize2,0:ncheb,lmsize2),srv(0:ncheb,lmsize2,0:ncheb,lmsize2) )
elseif ( use_sratrick==1 ) then
  allocate (work2((ncheb+1)*lmsize,(ncheb+1)*lmsize), ipiv2((ncheb+1)*lmsize))
  allocate (slv1(0:ncheb,lmsize,0:ncheb,lmsize),srv1(0:ncheb,lmsize,0:ncheb,lmsize), &
            slv2(0:ncheb,lmsize,0:ncheb,lmsize),srv2(0:ncheb,lmsize,0:ncheb,lmsize), &
            slv3(0:ncheb,lmsize,0:ncheb,lmsize),srv3(0:ncheb,lmsize,0:ncheb,lmsize) )
  allocate (yill1(0:ncheb,lmsize,lmsize),zill1(0:ncheb,lmsize,lmsize), &
            yrll1(0:ncheb,lmsize,lmsize),zrll1(0:ncheb,lmsize,lmsize), &
            yill2(0:ncheb,lmsize,lmsize),zill2(0:ncheb,lmsize,lmsize), &
            yrll2(0:ncheb,lmsize,lmsize),zrll2(0:ncheb,lmsize,lmsize),&
            yrlltmp(0:ncheb,lmsize,lmsize)  )
else
  stop '[rllsll] error with testflag sph'
end if

allocate( work(lmsize,lmsize),&
          allp(lmsize,lmsize,0:npan),bllp(lmsize,lmsize,0:npan),&
          cllp(lmsize,lmsize,0:npan),dllp(lmsize,lmsize,0:npan),&
          mrnvy(lmsize,lmsize,npan),mrnvz(lmsize,lmsize,npan),&
          mrjvy(lmsize,lmsize,npan),mrjvz(lmsize,lmsize,npan),&
          mihvy(lmsize,lmsize,npan),mihvz(lmsize,lmsize,npan),&
          mijvy(lmsize,lmsize,npan),mijvz(lmsize,lmsize,npan),&
          yill(0:ncheb,lmsize2,lmsize),zill(0:ncheb,lmsize2,lmsize),&
          yrll(0:ncheb,lmsize2,lmsize),zrll(0:ncheb,lmsize2,lmsize),&
          vjlr(lmsize,lmsize2,0:ncheb),vhlr(lmsize,lmsize2,0:ncheb),&
          vjli(lmsize,lmsize2,0:ncheb),vhli(lmsize,lmsize2,0:ncheb))

yrll=(0.0d0,0.0d0)
zill=(0.0d0,0.0d0)
yrll=(0.0d0,0.0d0)
zill=(0.0d0,0.0d0)

allocate( yif(lmsize2,lmsize,0:ncheb,npan),&
          yrf(lmsize2,lmsize,0:ncheb,npan),&
          zif(lmsize2,lmsize,0:ncheb,npan),&
          zrf(lmsize2,lmsize,0:ncheb,npan) )

do ipan = 1,npan
  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
    tau(icheb,ipan) = rmesh(mn)
  end do
end do

call chebint(cslc1,csrc1,slc1sum,c1,ncheb)

if (idotime==1) call timing_start('local')

do ipan = 1,npan

  if (idotime==1) call timing_start('local1')
  
  vhlr=czero
  vjlr=czero
  vhli=czero
  vjli=czero

  if (use_sratrick==0) then

    yrll=czero
    zrll=czero
    yill=czero
    zill=czero
  else
    yrll1=czero
    zrll1=czero
    yill1=czero
    zill1=czero
    yrll2=czero
    zrll2=czero
    yill2=czero
    zill2=czero
  end if

!---------------------------------------------------------------------
! 1. prepare VJLR, VNL, VHLR, which appear in the integrands
! TAU(K,IPAN) is used instead of TAU(K,IPAN)**2, which directly gives
! RLL(r) and SLL(r) multiplied with r. TAU is the radial mesh.
!
! 2. prepare the source terms YR, ZR, YI, ZI 
! because of the conventions used by
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! a factor sqrt(E) is included in the source terms
! this factor is removed by the definition of ZSLC1SUM given below
!
!vjlr = \kappa * J * V = \kappa * r * j *V 
!vhlr = \kappa * H * V = \kappa * r * h *V 
!
! i.e. prepare terms kappa*J*DV, kappa*H*DV appearing in 5.11, 5.12.

  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
    if     (cmoderll=='1') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
              vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do 
    elseif (cmoderll=='T') then ! transposed matrix (might not be needed anymore)
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
              vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmoderll=='0') then ! as a test option
              vjlr(:,:,icheb) = czero
              vhlr(:,:,icheb) = czero
    else
      stop'[rllsll] mode not known'
    end if

    if     (cmodesll=='1') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjli(lm1,lm2+lmsize*(ivec2-1),icheb) = vjli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
              vhli(lm1,lm2+lmsize*(ivec2-1),icheb) = vhli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmodesll=='T') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjli(lm1,lm2+lmsize*(ivec2-1),icheb) = vjli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
              vhli(lm1,lm2+lmsize*(ivec2-1),icheb) = vhli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmodesll=='0') then
              vjli(:,:,icheb) = czero
              vhli(:,:,icheb) = czero
    else
      stop'[rllsll] mode not known'
    end if

    ! calculation of the J (and H) matrix according to equation 5.69 (2nd eq.)
    if ( use_sratrick==0 ) then
      do ivec=1,nvec ! index for large/small component
        do lm1 = 1,lmsize
          l1 = jlk_index( lm1+lmsize*(ivec-1) )
          yrll(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn) 
          zrll(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn) 
          yill(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn)
          zill(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn)
        end do
      end do !ivec=1,nvec
    elseif ( use_sratrick==1 ) then
      do lm1 = 1,lmsize
        l1 = jlk_index( lm1+lmsize*(1-1) )
        l2 = jlk_index( lm1+lmsize*(2-1) )
        yrll1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn) 
        zrll1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn) 
        yill1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn)
        zill1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn)
        yrll2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l2,mn) 
        zrll2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l2,mn) 
        yill2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l2,mn)
        zill2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l2,mn)
      end do
    end if
  end do ! icheb

  ! calculation of A in 5.68
  if ( use_sratrick==0 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize2
          do ivec=1,nvec
            do lm3 = 1,lmsize
              lm1=lm3+(ivec-1)*lmsize
              l1 = jlk_index(lm1)
              slv(icheb,lm1,icheb2,lm2) = &
            ( tau(icheb,ipan)*jlk(l1,mn)*cslc1(icheb,icheb2)*vhlr(lm3,lm2,icheb2) &
              -tau(icheb,ipan)*hlk(l1,mn)*cslc1(icheb,icheb2)*vjlr(lm3,lm2,icheb2))&
            *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
              srv(icheb,lm1,icheb2,lm2) = &
            (-tau(icheb,ipan)*jlk(l1,mn)*csrc1(icheb,icheb2)*vhli(lm3,lm2,icheb2) &
              +tau(icheb,ipan)*hlk(l1,mn)*csrc1(icheb,icheb2)*vjli(lm3,lm2,icheb2)) &
                *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0
            end do
          end do
        end do
      end do
    end do
    do lm1 = 1,lmsize2
      do icheb = 0,ncheb
        slv(icheb,lm1,icheb,lm1) = slv(icheb,lm1,icheb,lm1) + 1.d0
        srv(icheb,lm1,icheb,lm1) = srv(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do
  elseif  ( use_sratrick==1 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize
          do ivec=1,1
            do lm3 = 1,lmsize
              lm1=lm3+(ivec-1)*lmsize
              l1 = jlk_index(lm1)

              ! this is the block to be inverted in SRAtrick. (named C in comment above):
              slv1(icheb,lm1,icheb2,lm2) = &      
            ( tau(icheb,ipan)*jlk(l1,mn)*cslc1(icheb,icheb2)*vhlr(lm3,lm2,icheb2) &
              -tau(icheb,ipan)*hlk(l1,mn)*cslc1(icheb,icheb2)*vjlr(lm3,lm2,icheb2))&
            *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0

              srv1(icheb,lm1,icheb2,lm2) = &
            (-tau(icheb,ipan)*jlk(l1,mn)*csrc1(icheb,icheb2)*vhli(lm3,lm2,icheb2) &
              +tau(icheb,ipan)*hlk(l1,mn)*csrc1(icheb,icheb2)*vjli(lm3,lm2,icheb2)) &
                *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0

            end do
          end do
        end do
      end do
    end do
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize
          do ivec=2,2
            do lm3 = 1,lmsize
              lm1=lm3+(ivec-1)*lmsize
              l1 = jlk_index(lm1)

              slv2(icheb,lm3,icheb2,lm2) = &
            ( tau(icheb,ipan)*jlk(l1,mn)*cslc1(icheb,icheb2)*vhlr(lm3,lm2,icheb2) &
              -tau(icheb,ipan)*hlk(l1,mn)*cslc1(icheb,icheb2)*vjlr(lm3,lm2,icheb2))&
            *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0

              srv2(icheb,lm3,icheb2,lm2) = &
            (-tau(icheb,ipan)*jlk(l1,mn)*csrc1(icheb,icheb2)*vhli(lm3,lm2,icheb2) &
              +tau(icheb,ipan)*hlk(l1,mn)*csrc1(icheb,icheb2)*vjli(lm3,lm2,icheb2)) &
                *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0

            end do
          end do
        end do
      end do
    end do
    do lm1 = 1,lmsize
      do icheb = 0,ncheb
        slv1(icheb,lm1,icheb,lm1) = slv1(icheb,lm1,icheb,lm1) + 1.d0
        srv1(icheb,lm1,icheb,lm1) = srv1(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do

  else
    stop '[rllsll] error in inversion'
  end if

  if (idotime==1) call timing_pause('local1')
  if (idotime==1) call timing_start('local2')

!-------------------------------------------------------
! determine the local solutions
! solve the equations SLV*YRLL=S and SLV*ZRLL=C 
!                 and SRV*YILL=C and SRV*ZILL=S
! i.e., solve system A*U=J, see eq. 5.68.

  if ( use_sratrick==0 ) then
    nplm = (ncheb+1)*lmsize2

    if (cmoderll/='0') then
      if (idotime==1) call timing_start('inversion')
      call zgetrf(nplm,nplm,slv,nplm,ipiv,info)
      if (idotime==1) call timing_stop('inversion','test')
      if (info/=0) stop'rllsll: zgetrf'
      call zgetrs('n',nplm,lmsize,slv,nplm,ipiv,yrll,nplm,info)
      call zgetrs('n',nplm,lmsize,slv,nplm,ipiv,zrll,nplm,info)
    end if
    if (cmodesll/='0') then
      if (directsolv==1) then
        call zgetrf(nplm,nplm,srv,nplm,ipiv,info)
        if (info/=0) stop'rllsll: zgetrf'
        call zgetrs('n',nplm,lmsize,srv,nplm,ipiv,yill,nplm,info)
        call zgetrs('n',nplm,lmsize,srv,nplm,ipiv,zill,nplm,info)
      else
        call iterativesol (ncheb,lmsize2,lmsize,srv,yill)
        call iterativesol (ncheb,lmsize2,lmsize,srv,zill)
      end if
    end if
  elseif ( use_sratrick==1 ) then
    nplm = (ncheb+1)*lmsize
    call inverse(nplm,slv1,work2,ipiv2)
    call inverse(nplm,srv1,work2,ipiv2)

    call zgemm('n','n',nplm,lmsize,nplm,cone,slv1, &
               nplm,yrll1,nplm,czero,yrlltmp,nplm)
    yrll1=yrlltmp
    call zgemm('n','n',nplm,lmsize,nplm,-cone,slv2, &
        nplm,yrll1,nplm,cone,yrll2,nplm)

    call zgemm('n','n',nplm,lmsize,nplm,cone,slv1, &
        nplm,zrll1,nplm,czero,yrlltmp,nplm)
    zrll1=yrlltmp
    call zgemm('n','n',nplm,lmsize,nplm,-cone,slv2, &
        nplm,zrll1,nplm,cone,zrll2,nplm)

    call zgemm('n','n',nplm,lmsize,nplm,cone,srv1, &
        nplm,yill1,nplm,czero,yrlltmp,nplm)
    yill1=yrlltmp
    call zgemm('n','n',nplm,lmsize,nplm,-cone,srv2, &
        nplm,yill1,nplm,cone,yill2,nplm)

    call zgemm('n','n',nplm,lmsize,nplm,cone,srv1, &
        nplm,zill1,nplm,czero,yrlltmp,nplm)
    zill1=yrlltmp
    call zgemm('n','n',nplm,lmsize,nplm,-cone,srv2, &
        nplm,zill1,nplm,cone,zill2,nplm)

  else
    stop '[rllsll] error in inversion'
  end if

  ! Reorient indices for later use
  if ( use_sratrick==0 ) then
    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize2
          yrf(lm1,lm2,icheb,ipan) = yrll(icheb,lm1,lm2)
          zrf(lm1,lm2,icheb,ipan) = zrll(icheb,lm1,lm2)
          yif(lm1,lm2,icheb,ipan) = yill(icheb,lm1,lm2)
          zif(lm1,lm2,icheb,ipan) = zill(icheb,lm1,lm2)
        end do
      end do
    end do

  elseif ( use_sratrick==1 ) then

    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrf(lm1,lm2,icheb,ipan) = yrll1(icheb,lm1,lm2)
          zrf(lm1,lm2,icheb,ipan) = zrll1(icheb,lm1,lm2)
          yif(lm1,lm2,icheb,ipan) = yill1(icheb,lm1,lm2)
          zif(lm1,lm2,icheb,ipan) = zill1(icheb,lm1,lm2)
        end do
      end do
    end do

    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrf(lm1+lmsize,lm2,icheb,ipan) = yrll2(icheb,lm1,lm2)
          zrf(lm1+lmsize,lm2,icheb,ipan) = zrll2(icheb,lm1,lm2)
          yif(lm1+lmsize,lm2,icheb,ipan) = yill2(icheb,lm1,lm2)
          zif(lm1+lmsize,lm2,icheb,ipan) = zill2(icheb,lm1,lm2)
        end do
      end do
    end do

  end if

  if (idotime==1) call timing_pause('local2')
  if (idotime==1) call timing_start('local3')

  ! Calculation of eq. 5.19-5.22

  do icheb = 0,ncheb
    zslc1sum(icheb) = slc1sum(icheb) * (rpanbound(ipan)-rpanbound(ipan-1))/ (2.d0)
  end do
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhlr(1,1,0), &
        lmsize,yrf(1,1,0,ipan),lmsize2,czero,mrnvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjlr(1,1,0), &
        lmsize,yrf(1,1,0,ipan),lmsize2,czero,mrjvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhlr(1,1,0), &
        lmsize,zrf(1,1,0,ipan),lmsize2,czero,mrnvz(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjlr(1,1,0), &
        lmsize,zrf(1,1,0,ipan),lmsize2,czero,mrjvz(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhli(1,1,0), &
        lmsize,yif(1,1,0,ipan),lmsize2,czero,mihvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjli(1,1,0), &
        lmsize,yif(1,1,0,ipan),lmsize2,czero,mijvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhli(1,1,0), &
        lmsize,zif(1,1,0,ipan),lmsize2,czero,mihvz(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjli(1,1,0), &
        lmsize,zif(1,1,0,ipan),lmsize2,czero,mijvz(1,1,ipan),lmsize)
  do icheb = 1,ncheb
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhlr(1,1,icheb), &
          lmsize,yrf(1,1,icheb,ipan),lmsize2,cone,mrnvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjlr(1,1,icheb), &
          lmsize,yrf(1,1,icheb,ipan),lmsize2,cone,mrjvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhlr(1,1,icheb), &
          lmsize,zrf(1,1,icheb,ipan),lmsize2,cone,mrnvz(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjlr(1,1,icheb), &
          lmsize,zrf(1,1,icheb,ipan),lmsize2,cone,mrjvz(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhli(1,1,icheb), &
          lmsize,yif(1,1,icheb,ipan),lmsize2,cone,mihvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjli(1,1,icheb), &
          lmsize,yif(1,1,icheb,ipan),lmsize2,cone,mijvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhli(1,1,icheb), &
          lmsize,zif(1,1,icheb,ipan),lmsize2,cone,mihvz(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjli(1,1,icheb), &
          lmsize,zif(1,1,icheb,ipan),lmsize2,cone,mijvz(1,1,ipan),lmsize)
  end do
  if (idotime==1) call timing_pause('local3')

end do !ipan
! end the big loop over the subintervals



if (idotime==1) call timing_stop('local')
if (idotime==1) call timing_start('afterlocal')

! ***********************************************************************
! calculate A(M), B(M), C(M), D(M)
! according to 5.17-5.18 (regular solution) of Bauer PhD
! C,D are calculated accordingly for the irregular solution
! (starting condition: A(0) = 1, B(0) = 0, C(MMAX) = 0 and D(MMAX) = 1)
! ***********************************************************************

! regular 
do lm2 = 1,lmsize
  do lm1 = 1,lmsize
    bllp(lm1,lm2,0) = czero
    allp(lm1,lm2,0) = czero
  end do
end do

do lm1 = 1,lmsize
  allp(lm1,lm1,0) = cone
end do

do ipan = 1,npan
  call zcopy(lmsize*lmsize,allp(1,1,ipan-1),1,allp(1,1,ipan),1)
  call zcopy(lmsize*lmsize,bllp(1,1,ipan-1),1,bllp(1,1,ipan),1)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvy(1,1,ipan), &
          lmsize,allp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvz(1,1,ipan), &
          lmsize,bllp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mrjvy(1,1,ipan), &
          lmsize,allp(1,1,ipan-1),lmsize,cone,bllp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mrjvz(1,1,ipan), &
          lMSIZE,BLLP(1,1,IPAN-1),LMSIZE,CONE,BLLP(1,1,IPAN),LMSIZE)
end do

! irregular 
do lm2 = 1,lmsize
  do lm1 = 1,lmsize
    dllp(lm1,lm2,npan) = 0.d0
    cllp(lm1,lm2,npan) = 0.d0
  end do
end do
do lm1 = 1,lmsize
  dllp(lm1,lm1,npan) = 1.d0
end do
do ipan = npan,1,-1
  call zcopy(lmsize*lmsize,cllp(1,1,ipan),1,cllp(1,1,ipan-1),1)
  call zcopy(lmsize*lmsize,dllp(1,1,ipan),1,dllp(1,1,ipan-1),1)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mihvz(1,1,ipan), &
          lmsize,cllp(1,1,ipan),lmsize,cone,cllp(1,1,ipan-1),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mihvy(1,1,ipan), &
          lmsize,dllp(1,1,ipan),lmsize,cone,cllp(1,1,ipan-1),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mijvz(1,1,ipan), &
          lmsize,cllp(1,1,ipan),lmsize,cone,dllp(1,1,ipan-1),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mijvy(1,1,ipan), &
          lmsize,dllp(1,1,ipan),lmsize,cone,dllp(1,1,ipan-1),lmsize)
end do

! ***********************************************************************
! determine the regular solution ull by using 5.14
! and the irregular solution sll accordingly
! ***********************************************************************
do ipan = 1,npan
  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,yrf(1,1,icheb,ipan), &
          lmsize2,allp(1,1,ipan-1),lmsize,czero,ull(1,1,mn),lmsize2)
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,zrf(1,1,icheb,ipan), &
          lmsize2,bllp(1,1,ipan-1),lmsize,cone,ull(1,1,mn),lmsize2)
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,zif(1,1,icheb,ipan), &
          lmsize2,cllp(1,1,ipan),lmsize,czero,sll(1,1,mn),lmsize2)
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,yif(1,1,icheb,ipan), &
          lmsize2,dllp(1,1,ipan),lmsize,cone,sll(1,1,mn),lmsize2)
  end do
end do

if (idotime==1) call timing_stop('afterlocal')
if (idotime==1) call timing_start('endstuff')

! ***********************************************************************
! next part converts the volterra solution u of equation (5.7) to
! the fredholm solution r by employing eq. 4.122 and 4.120 of bauer, phd
! and the t-matrix is calculated
! ***********************************************************************

call zgetrf(lmsize,lmsize,allp(1,1,npan),lmsize,ipiv,info)                     !invert alpha
call zgetri(lmsize,allp(1,1,npan),lmsize,ipiv,work,lmsize*lmsize,info)         !invert alpha -> transformation matrix rll=alpha^-1*rll
! calculation of the t-matrix 
call zgemm('n','n',lmsize,lmsize,lmsize,cone/gmatprefactor,bllp(1,1,npan), &   ! calc t-matrix tll = bll*alpha^-1 
            lmsize,allp(1,1,npan),lmsize,czero,tllp,lmsize)

do nm = 1,nrmax
call zgemm('n','n',lmsize2,lmsize,lmsize,cone,ull(1,1,nm), &
            lmsize2,allp(1,1,npan),lmsize,czero,rll(1,1,nm),lmsize2)
end do

if (idotime==1) call timing_stop('endstuff')
if (idotime==1) call timing_start('checknan')
if (idotime==1) call timing_stop('checknan')
if (idotime==1) call timing_stop('local1')
if (idotime==1) call timing_stop('local2')
if (idotime==1) call timing_stop('local3')
if (idotime==1) call timing_stop('rllsll')

end subroutine

! subroutine inverse(nmat,mat)
! !interface
! integer        :: nmat
! double complex :: mat(nmat,nmat)
! double complex :: work3(nmat,nmat)
! !local
! integer        :: IPIV3(nmat)
! integer        :: info

! call ZGETRF( nmat, nmat, mat, nmat, IPIV3, INFO )
! if (info/=0) stop '[inverse] error INFO' 
! call ZGETRI( nmat, mat, nmat, IPIV3, WORK3, nmat*nmat, INFO )
! if (info/=0) stop '[inverse] error INFO' 
! end subroutine inverse


subroutine iterativesol (NCHEB,LMSIZE2,LMSIZE,MMAT,BMAT)
integer :: NCHEB
integer :: LMSIZE,LMSIZE2
double complex :: MMAT(0:NCHEB,LMSIZE2,0:NCHEB,LMSIZE2)
double complex :: BMAT(0:NCHEB,LMSIZE2,LMSIZE)
double complex :: XMAT(0:NCHEB,LMSIZE2,LMSIZE)
!########################################################
! solves the system of linear equations
! MMAT*XMAT = BMAT
!########################################################

NPLM = (NCHEB+1)*LMSIZE2
CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,CONE,SRV, &
    NPLM,ZILL,NPLM,CZERO,OUT,NPLM)



end subroutine iterativesol



END MODULE MOD_RLLSLL


