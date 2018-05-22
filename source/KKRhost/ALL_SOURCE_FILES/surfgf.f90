!-------------------------------------------------------------------------------
! SUBROUTINE: SURFGF
!> @brief Solve surface green's function: \f$ f(x)=ml\left(m0-x\right)^{\left(-1\right)*mr} \f$
!> @details method: decimation technique
!> input:  ml,m0,mr - complex rectangular matrices
!> output: x        - result, matrix of same type as before
!> @note NEW VERSION (speeded up) by V.Bellini (march,1999)
!>
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine SURFGF(NDIM,ML,M0,MR,X,ITERMAX,ERRMAX,ICHCK,LMMAXD)

   use Constants
   use Profiling
   use global_variables

   implicit none

   !-------------------------------------------------------------------------------
   ! For KREL = 1 (relativistic mode)
   !
   !  NPOTD = 2 * NATYPD
   !  LMMAXD = 2 * (LMAXD+1)^2
   !  NSPIND = 1
   !  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green
   !          function, set up in the spin-independent non-relativstic
   !          (l,m_l)-representation
   !
   !-------------------------------------------------------------------------------
   ! .. Input variables
   integer, intent(in) :: NDIM
   integer, intent(in) :: ICHCK
   integer, intent(in) :: LMMAXD       !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: ITERMAX
   double precision, intent(in) :: ERRMAX
   ! .. Input arrays
   double complex, dimension(NDIM,NDIM), intent(in) :: M0
   double complex, dimension(NDIM,NDIM), intent(in) :: ML
   double complex, dimension(NDIM,NDIM), intent(in) :: MR
   ! .. Output variables
   double complex, dimension(NDIM,NDIM), intent(out) :: X
   ! .. Local Scalars
   integer :: I,INFO,ITER,J,N
   double precision :: ERR,SUM,XIM,XRE
   ! .. Local Arrays
   integer, dimension(NDIM) :: IPVT
   double complex, dimension(NDIM,NDIM) :: AA,ALFA,BB,BETA,CC,CUNIT,EPS,TEMPIN
   double complex, dimension(NDIM,NDIM) :: TEMPOUT,Y1,Y2
   ! .. External Subroutines ..
   external CINIT,ZAXPY,ZCOPY,ZGEMM,ZGETRF,ZGETRS
   ! .. Intrinsic Functions ..
   intrinsic :: DBLE,DIMAG
   !     ..

   call CINIT(NDIM*NDIM,CUNIT)
   do N = 1,NDIM
      CUNIT(N,N) = CONE
   end do

   call ZCOPY(NDIM*NDIM,M0,1,EPS,1)
   call ZCOPY(NDIM*NDIM,ML,1,ALFA,1)
   call ZCOPY(NDIM*NDIM,MR,1,BETA,1)
   call ZCOPY(NDIM*NDIM,M0,1,X,1)

   ITER = 1
   10 continue

   call ZCOPY(NDIM*NDIM,EPS,1,Y1,1)
   call ZCOPY(NDIM*NDIM,Y1,1,TEMPIN,1)
   call ZGETRF(NDIM,NDIM,TEMPIN,NDIM,IPVT,INFO)

   !     aa = eps^-1 * alfa
   call ZCOPY(NDIM*NDIM,ALFA,1,TEMPOUT,1)
   call ZGETRS('N',NDIM,NDIM,TEMPIN,NDIM,IPVT,TEMPOUT,NDIM,INFO)
   call ZCOPY(NDIM*NDIM,TEMPOUT,1,AA,1)

   !     bb = eps^-1 * beta

   call ZCOPY(NDIM*NDIM,BETA,1,TEMPOUT,1)
   call ZGETRS('N',NDIM,NDIM,TEMPIN,NDIM,IPVT,TEMPOUT,NDIM,INFO)
   call ZCOPY(NDIM*NDIM,TEMPOUT,1,BB,1)

   !     alfa_new = alfa * aa

   call ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,ALFA,NDIM,AA,NDIM,CZERO,Y1,NDIM)

   !     beta_new = beta * bb

   call ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,BETA,NDIM,BB,NDIM,CZERO,Y2,NDIM)

   !     cc = - alfa * bb

   call ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,ALFA,NDIM,BB,NDIM,CZERO,CC,NDIM)

   !     x_new = x + cc

   call ZAXPY(NDIM*NDIM,CONE,CC,1,X,1)

   !     cc = eps + cc

   call ZAXPY(NDIM*NDIM,CONE,CC,1,EPS,1)

   !     eps_new = cc - beta * aa

   call ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,BETA,NDIM,AA,NDIM,CONE,EPS,NDIM)

   call ZCOPY(NDIM*NDIM,Y1,1,ALFA,1)
   call ZCOPY(NDIM*NDIM,Y2,1,BETA,1)

   SUM = 0.d0
   do I = 1,NDIM
      do J = 1,NDIM
         XRE = DBLE(ALFA(I,J))
         XIM = DIMAG(ALFA(I,J))
         SUM = SUM + XRE*XRE + XIM*XIM
      end do
   end do

   ERR = SQRT(SUM)
   if (ERR.LT.ERRMAX .OR. ITER.GT.ITERMAX) go to 20
   ITER = ITER + 1
   go to 10

   20 continue

   call ZCOPY(NDIM*NDIM,X,1,TEMPIN,1)
   call ZCOPY(NDIM*NDIM,CUNIT,1,TEMPOUT,1)
   call ZGETRF(NDIM,NDIM,TEMPIN,NDIM,IPVT,INFO)
   call ZGETRS('N',NDIM,NDIM,TEMPIN,NDIM,IPVT,TEMPOUT,NDIM,INFO)
   call ZCOPY(NDIM*NDIM,TEMPOUT,1,X,1)

   call ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,X,NDIM,MR,NDIM,CZERO,TEMPIN,NDIM)
   call ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,ML,NDIM,TEMPIN,NDIM,CZERO,X,NDIM)

   if (ITER.GT.ITERMAX) then
      write (6,FMT='('' itermax too small.  iter='',i3)') ITER
      write(6,'('' Surfgf:  iter='',i4,''  error='',d14.7)') iter,err
   end if
   if (ICHCK.EQ.0) return
   !      write(6,'('' Surfgf:  iter='',i4,''  error='',d12.7)') iter,err
   !      write(6,'(/'' X matrix'')')
   !      write(6,*)
   !
   return
end subroutine SURFGF
