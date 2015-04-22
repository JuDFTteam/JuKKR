      SUBROUTINE SLAYDIRECT(LPOT,VEC1,VEC2,ALAT,BR,
     &                      SUM)
      implicit none
c -----------------------------------------------------------------
c   This subroutine returns 
c
c                          i1   i2
c                     Y  (R  - R   )
c        i1,i2  ---    L        n    
c     SUM    =  \    -----------------
cc       L      /       |  i1   i2 | l+1 
cc              ---     | R  - R   |
c                n      |       n  |          
c
c
c
c
c   i1,i2 are plane indeces and the sum is done in the 
c        i2 plane  
c   In case i1=i2 the r=0 term is excluded! 
c   
c   Error check is performed!
c -----------------------------------------------------------------  
c     .. parameters ..
      include 'inc.p'
      integer LMPOTD,L2POTD,LM2POTD
      PARAMETER (LMPOTD=(LPOTD+1)**2,LM2POTD=(2*LPOTD+1)**2)
      PARAMETER (L2POTD=2*LPOTD)
      !integer nmaxd,ishld
      !PARAMETER (NMAXD=19500,ISHLD=1500)      

c     ..
c     .. array arguments ..
      
c
      double precision BR(3,3)
      double precision vec1(3),vec2(3)
      double precision sum(LM2POTD),sumt(lm2potd)
c---  > integer variables
c   
      INTEGER I1,I2
      INTEGER L,lm,M,K,LPOT,n
      INTEGER I,J,lm2pot,N1,MAXN,IA,N2
      integer ipar
c---  > real variables
c
      DOUBLE PRECISION ALAT
      DOUBLE PRECISION ZZZ,R,CC,ZZ,CR,X0,Y0,R2
      double precision zoffx,zoffy,rmax
      double precision ccc,ylm(LM2POTD),r0,a1,a2,b1,b2,rtest
      logical ltest,TEST
C ------------------------------------------------------------
         do lm=1,lm2potd
           sum(lm)  = 0.d0
           sumT(lm) = 0.d0
         end do
c Run this sub only in the test option "electro"       
         IF (.NOT.TEST('electro ')) RETURN
       
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
            do 10 n2=-MAXN,MAXN
               x0 = zoffx - (a1*N1+b1*N2) 
               y0 = zoffy - (a2*N1+b2*N2) 
               R2 = x0*x0 + y0*y0 + zzz
               r0 = (a1*N1+b1*N2)**2+(a2*N1+b2*N2)**2
               IF (r2.LE.1.D-6) GO TO 10
               if (sqrt(r0).gt.rmax) goto 10
               ipar = ipar +n1+n2
               IA = IA + 1 
               CALL YMY(x0,y0,zz,R,YLM,L2POTD)
               CR = 1.d0/r
               do l=1,2*lpot
                  CC = CR**(L+1)  
                  do m=-l,l
                     lm =  l*(l+1)+m+1
                     CCC = CC*YLM(LM)
                     sum(lm) = sum(lm)+ccc
                  end do
               end do
c     test the convergence
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
 10         CONTINUE
            
         end do                 ! all atoms in plane 

         if (ipar.ne.0) write(6,*) 'SUMLAYER  Asymetric sum '             
         LM2POT=(2*Lpot+1)**2
c   Now test the sum
         IF (TEST('electro '))THEN
         do lm=2 ,lm2pot ! test only l=3 terms!
            if (dabs(sum(lm)).gt.1.d-8) then ! test 
        write(6,100) vec1(1),vec1(2),vec2(1),vec2(2),lm,sum(lm),sumt(lm)
            end if 
         end do
         END IF         
 100     FORMAT('SUMLAYER ',4F8.4,' LM=',I3,' SUM=',2D12.5)
 101     format(I5,5F12.6)
c         do lm=2,lm2pot
c            if (dabs(sum(lm)).gt.1.d-8) then
c               write(6,4005) i1,i2,lm,sum(lm)
c            end if
c         end do
c     -------------------------------------------------------------
         
 4005    FORMAT(1x,' I1, I2 LM, DIRECT SUM=',3I4,(D18.9,D18.9))
         END  
