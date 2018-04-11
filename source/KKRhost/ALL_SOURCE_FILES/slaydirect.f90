!-------------------------------------------------------------------------------
! SUBROUTINE: SLAYDIRECT
!> @brief This subroutine returns
!> \f$ SUM_L^{i_1,i_2}=\sum_n \frac{Y_L\left(R^{i_1}-R_{n}^{i_2}\right)}{\left|R^{i_1}-R^{i_2}_n \right|^{l+1}}\f$
!
!> @details  \f$i_1\f$, \f$i_2\f$ are plane indeces and the sum is done in the
!> \f$i_2\f$ plane
!>
!> In case \f$i_1=i_2\f$ the \f$r=0\f$ term is excluded!
!>
!> Error check is performed!
!> @note Jonathan Chico Apr. 2018: Removed inc.p dependencies, rewrote to Fortran90
!> and renamed the SUM array to SUML to avoid confusion with intrinsic SUM function
!-------------------------------------------------------------------------------
subroutine SLAYDIRECT(LPOT,VEC1,VEC2,ALAT,BR,SUML)

   use Constants
   use global_variables

   implicit none

   ! .. Parameters
   integer :: LMPOTD
   integer :: L2POTD
   integer :: LM2POTD
   parameter (L2POTD=2*LPOTD)
   parameter (LMPOTD=(LPOTD+1)**2)
   parameter (LM2POTD=(2*LPOTD+1)**2)
   ! .. Input variables
   double precision, intent(in) :: ALAT   !< Lattice constant in a.u.
   double precision, dimension(3), intent(in) :: vec1
   double precision, dimension(3), intent(in) :: vec2
   double precision, dimension(3,3), intent(in) :: BR
   ! .. Output variables
   double precision, dimension(LM2POTD), intent(out) :: SUML
   ! .. Local variables
   integer :: I1,I2
   integer :: L,lm,M,K,LPOT,n
   integer :: I,J,lm2pot,N1,MAXN,IA,N2
   integer :: ipar
   double precision :: ZZZ,R,CC,ZZ,CR,X0,Y0,R2
   double precision :: zoffx,zoffy,rmax
   double precision :: r0,a1,a2,b1,b2,rtest, ccc
   logical :: ltest,TEST
   double precision, dimension(LM2POTD) :: ylm
   double precision, dimension(LM2POTD) :: sumt
   !----------------------------------------------------------------------------
   do lm=1,lm2potd
      suml(lm)  = 0.d0
      sumT(lm) = 0.d0
   end do
   ! Run this sub only in the test option "electro"
   if (.NOT.TEST('electro ')) return

   ZOFFX = (VEC2(1) - VEC1(1))*ALAT
   ZOFFY = (VEC2(2) - VEC1(2))*ALAT
   ZZ    = (VEC2(3) - VEC1(3))*ALAT
   ZZZ   = ZZ*ZZ
   MAXN =130
   IA = 0
   ipar = 0
   ltest=.true.
   a1 = br(1,1)
   a2 = br(2,1)
   b1 = br(1,2)
   b2 = br(2,2)
   rmax = MAXN* MIN(sqrt(a1*a1+a2*a2),sqrt(b1*b1+b2*b2))
   rmax = rmax*1.00001d0
   rtest = (MAXN-20)*MIN(sqrt(a1*a1+a2*a2),sqrt(b1*b1+b2*b2))
   rtest = rtest*1.00001d0

   do n1=-MAXN,MAXN
      n2=-MAXN
      R2=1e9
      R0=0.d0
      do while(n2<=MAXN.and.(r2.LE.1.D-6).and.(sqrt(r0).gt.rmax))
         n2=n2+1
         !do 10 n2=-MAXN,MAXN
         x0 = zoffx - (a1*N1+b1*N2)
         y0 = zoffy - (a2*N1+b2*N2)
         R2 = x0*x0 + y0*y0 + zzz
         r0 = (a1*N1+b1*N2)**2+(a2*N1+b2*N2)**2
         !IF (r2.LE.1.D-6) GO TO 10
         !if (sqrt(r0).gt.rmax) goto 10
         ipar = ipar +n1+n2
         IA = IA + 1
         call YMY(x0,y0,zz,R,YLM,L2POTD)
         CR = 1.d0/r
         do l=1,2*lpot
            CC = CR**(L+1)
            do m=-l,l
               lm =  l*(l+1)+m+1
               CCC = CC*YLM(LM)
               suml(lm) = suml(lm)+ccc
            end do
         end do
         !     test the convergence
         if (sqrt(r0).gt.rtest) then
            do l=1,2*lpot
               CC = CR**(L+1)
               do m=-l,l
                  lm =  l*(l+1)+m+1
                  CCC = CC*YLM(LM)
                  sumt(lm) = sumt(lm)+ccc
               end do
            end do
         end if
      enddo ! Do while
   end do ! all atoms in plane

   if (ipar.ne.0) write(1337,*) 'SUMLAYER  Asymetric sum '
   LM2POT=(2*LPOT+1)**2
   ! Now test the sum
   IF (TEST('electro '))THEN
      do lm=2 ,lm2pot ! test only l=3 terms!
         if (abs(suml(lm)).gt.1.d-8) then ! test
            write(1337,100) vec1(1),vec1(2),vec2(1),vec2(2),lm,suml(lm),sumt(lm)
         end if
      end do
   end if
   100 format('SUMLAYER ',4F8.4,' LM=',I3,' SUM=',2D12.5)
   101 format(I5,5F12.6)
   !         do lm=2,lm2pot
   !            if (dabs(sum(lm)).gt.1.d-8) then
   !               write(6,4005) i1,i2,lm,sum(lm)
   !            end if
   !         end do
   !     -------------------------------------------------------------
   4005 format(1x,' I1, I2 LM, DIRECT SUM=',3I4,(D18.9,D18.9))

end subroutine SLAYDIRECT
