  ! **********************************************************************
  ! *                                                                    *
  !>*  calculation of lattice sums for l .le. 2*lpot :                   *
  !>*                                                                    *
  !>*                   ylm( q(i) - q(j) - rm )                          *
  !>*        sum      ===========================                        *
  !>*                 | q(i) - q(j) - rm |**(l+1)                        *
  !>*                                                                    *
  !>*         - summed over all lattice vectors rm  -                    *
  !>*                                                                    *
  !>*  ylm       : real spherical harmic to given l,m                    *
  !>*  q(i),q(j) : basis vectors of the unit cell                        *
  !>*                                                                    *
  !>*  in the case of i = j, rm = 0 is omitted.                          *
  !>*                                                                    *
  !>*  the ewald method is used to perform the lattice summations        *
  !>*  the splitting parameter lamda is set equal sqrt(pi)/alat          *
  !>*  (alat is the lattice constant) .                                  *
  !>*                                                                    *
  !>*  if the contribution of the last shell of the direct and the       *
  !>*  reciprocal lattice is greater than 1.0e-8 a message is written    *
  !>*                                                                    *
  !>*                                    b.drittler may 1989             *
  !>*                                                                    *
  !>*  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
  !>*  one not being used (see also lattice3d)     v.popescu May 2004    *
  ! *                                                                    *
  ! **********************************************************************

! OpenMP parallelised, needs threadsafe erfcex, gamfc and ymy E.R.

subroutine strmat(alat,lpot,naez,ngmax,nrmax,nsg,nsr,nshlg,nshlr, &
     gn,rm,qi0,smat,vol,lassld,lmxspd,naezd,i1)

  implicit none
  ! Parameters
  double complex, parameter :: CI =(0d0,1d0)
  double precision, parameter :: BOUND =1d-8

  ! Arguments
  double precision :: alat
  double precision :: vol
  integer :: lpot
  integer :: naez
  integer :: ngmax
  integer :: nrmax
  integer :: nshlg
  integer :: nshlr
  integer :: lassld
  integer :: lmxspd
  integer :: naezd
  double precision, dimension(3,*) :: gn
  double precision, dimension(3,*) :: qi0
  double precision, dimension(3,*) :: rm
  double precision, dimension(lmxspd,*) :: smat
  integer, dimension(*) :: nsg
  integer, dimension(*) :: nsr
  integer :: i1

 !local variables of strmat
 double complex :: bfac
 double precision :: alpha
 double precision :: beta
 double precision :: dq1
 double precision :: dq2
 double precision :: dq3
 double precision :: dqdotg
 double precision :: expbsq
 double precision :: fpi
 double precision :: g1
 double precision :: g2
 double precision :: g3
 double precision :: ga
 double precision :: lamda
 double precision :: pi
 double precision :: r
 double precision :: r1
 double precision :: r2
 double precision :: r3
 double precision :: rfac
 double precision :: s
 integer :: i
 integer :: i2
 integer :: it
 integer :: l
 integer :: lm
 integer :: lmx
 integer :: lmxsp
 integer :: m
 integer :: nge
 integer :: ngs
 integer :: nre
 integer :: nrs
 integer :: nstart

 double complex, dimension(lmxspd) :: stest
 double precision, dimension(0:lassld) :: g
 double precision, dimension(lmxspd) :: ylm

  !     ..................................................................

  lmx = 2*lpot
  lmxsp = (lmx+1)*(lmx+1)
  pi = 4.0d0*atan(1.0d0)
  fpi = 4.0d0*pi

  ! --> choose proper splitting parameter

  lamda = sqrt(pi)/alat

  ! **********************************************************************
  !$omp parallel do private(I2,DQ1,DQ2,DQ3,STEST,LM,NSTART,IT, &
  !$omp                     NRS,NGS,NRE,NGE,I,R1,R2,R3, &
  !$omp                     YLM,R,ALPHA,G,RFAC,L,M, &
  !$omp                     G1,G2,G3,GA,BETA,EXPBSQ,DQDOTG,BFAC,S)
  do i2 = 1,naez
     !======================================================================
     dq1 = (qi0(1,i1) - qi0(1,i2)) * alat
     dq2 = (qi0(2,i1) - qi0(2,i2)) * alat
     dq3 = (qi0(3,i1) - qi0(3,i2)) * alat

     stest(1) = -sqrt(fpi)/vol/(4d0*lamda*lamda)
     do lm = 2,lmxsp
        stest(lm) = 0.0d0
     end do

     ! --> exclude the origine and add correction if i1.eq.i2

     if ( i1.eq.i2 ) then
        stest(1) = stest(1) - lamda/pi
        nstart = 2
     else
        nstart = 1
     end if
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     ! --> loop first over n-1 shells of real and reciprocal lattice - then
     !     add the contribution of the last shells to see convergence

     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     do it = 1,2
        if ( it.eq.1 ) then
           nrs = nstart
           ngs = 2
           nre = nrmax - nsr(nshlr)
           nge = ngmax - nsg(nshlg)
        else
           nrs = nre + 1
           ngs = nge + 1
           nre = nrmax
           nge = ngmax
        end if

        ! --> sum over real lattice

        ! ---------------------------------------------------------------------
        do i = nrs,nre
           r1 = dq1 - rm(1,i)
           r2 = dq2 - rm(2,i)
           r3 = dq3 - rm(3,i)

           call ymy(r1,r2,r3,r,ylm,lmx)
           alpha = lamda*r
           call gamfc(alpha,g,lmx,r)

           do l = 0,lmx
              rfac = g(l)/sqrt(pi)
              do m = -l,l
                 lm = l*(l+1) + m + 1
                 stest(lm) = stest(lm) + ylm(lm)*rfac
              end do
           end do

        end do
        ! ---------------------------------------------------------------------

        ! --> sum over reciprocal lattice

        ! ---------------------------------------------------------------------
        do i = ngs,nge
           g1 = gn(1,i)
           g2 = gn(2,i)
           g3 = gn(3,i)

           call ymy(g1,g2,g3,ga,ylm,lmx)
           beta = ga/lamda
           expbsq = exp(beta*beta/4.0d0)
           dqdotg = dq1*g1 + dq2*g2 + dq3*g3

           bfac = fpi*exp(CI*dqdotg)/(ga*ga*expbsq*vol)

           do l = 0,lmx
              do m = -l,l
                 lm = l*(l+1) + m + 1
                 stest(lm) = stest(lm) + ylm(lm)*bfac
              end do
              bfac = bfac*ga/CI/dble(2*l+1)
           end do
        end do
        ! ---------------------------------------------------------------------
        if ( it.eq.1 ) then
           do lm = 1,lmxsp
              if ( abs(dimag(stest(lm))).gt.BOUND ) then
                 write (6,*) ' ERROR: Imaginary contribution', &
                      ' to REAL lattice sum'
                 stop
              end if
              smat(lm,i2) = dble(stest(lm))
              stest(lm) = 0.0d0
           end do
        else

           ! --> test convergence

           do lm = 1,lmxsp
              s = dble(stest(lm))
              smat(lm,i2) = smat(lm,i2) + s
              !IF (2.LT.1.AND. ABS(S).GT.BOUND ) WRITE (6,FMT=99001) I1,I2, &
              !LM,ABS(S)
           end do
        end if
        ! ---------------------------------------------------------------------
     end do
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do   ! I2
  !$omp end parallel do
  ! **********************************************************************

99001 format (5x,'WARNING : Convergence of SMAT(',i2,',',i2,') ', &
       ' for LMXSP =',i3,' is ',1p,d8.2,' > 1D-8',/,15x, &
       'You should use more lattice vectors (RMAX/GMAX)')
end subroutine strmat
