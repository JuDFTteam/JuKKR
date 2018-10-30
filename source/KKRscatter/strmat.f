      subroutine strmat(alat,lpot,natom,ngmax,nrmax,nsg,nsr,nshlg,nshlr,
     +                  gn,rm,qi0,smat,vol)
      implicit none
c-----------------------------------------------------------------------
c
c     calculation of lattice sums for l .le. 2*lpot :
c
c                      ylm( q(i) - q(j) - rm )
c           sum      ===========================
c                    | q(i) - q(j) - rm |**(l+1)
c
c            - summed over all lattice vectors rm  -
c
c     ylm       : real spherical harmic to given l,m
c     q(i),q(j) : basis vectors of the unit cell
c
c     in the case of i = j is rm = 0 omitted .
c
c     the ewald method is used to perform the lattice summations
c     the splitting parameter lamda is set equal sqrt(pi)/alat
c     (alat is the lattice constant) .
c
c     if the contribution of the last shell of the direct and the
c     reciprocal lattice is greater than 1.0e-8 a message is written
c
c                                     b.drittler may 1989
c
c-----------------------------------------------------------------------
c     .. parameters ..
c      integer lmaxd,lpotd
c      PARAMETER (LMAXD=16,LPOTD=16)
      include 'inc.p'
      integer l2maxd,l2mmxd
      parameter (l2maxd=2*LPOTD,l2mmxd= (l2maxd+1)**2)
c     ..
c     .. scalar arguments ..
      double precision alat,vol
      integer lpot,natom,ngmax,nrmax,nshlg,nshlr
c     ..
c     .. array arguments ..
      double precision  gn(4,*),qi0(3,*),qi(3,NAEZD+NEMBD),rm(4,*),
     +     smat(NAEZD,NAEZD,l2mmxd)
      integer nsg(*),nsr(*)
c     ..
c     .. local scalars ..
      double complex bfac,cfac,ci
      double precision alpha,beta,bound,dq1,dq2,dq3,dqdotg,expbsq,fpi,
     +                 g1,g2,g3,ga,
     +                 lamda,pi,r,r1,r2,r3,rfac,s,sgm
      integer i,i1,i2,it,l,lm,lmax,lmmax,m,nge,ngs,nre,nrs,nstart
      logical TEST
c     ..
c     .. local arrays ..
      double complex stest(l2mmxd)
      double precision g(0:l2maxd),ylm(l2mmxd)
c     ..
c     .. external subroutines ..
      external gamfc,ymy
      external test      
c     ..
c     .. intrinsic functions ..
      intrinsic abs,aimag,atan,exp,real,sqrt
c     ..
c     .. save statement ..
      save ci,bound
c     ..
c     .. data statements ..
      data ci/ (0.0d0,1.0d0)/,bound/1.0d-7/
c     ..
      lmax = 2*lpot
      lmmax = (lmax+1)* (lmax+1)
      pi = 4.0d0*datan(1.0d0)
      fpi = 4.0d0*pi
      write(6,*) '>>>>>>>>> STRMAT'
c
c---> choose proper splitting parameter
c
      lamda = sqrt(pi)/alat
c
c---> loop over atoms per unit cell
c
      write(6,*) '>>>>>> In STRMAT'
c scale basis atoms with alat (this might change in future versions)
      do i2=1,natom
        do i1 =1,3
          qi(i1,i2) = qi0(i1,i2)*alat
        end do
       ! write(6,*) (qi(i1,i2),i1=1,3)
      end do 
ccccccccccccccccccc
      do 10 i1 = 1,natom
         do 20 i2 = 1,natom
c
            dq1 = qi(1,i1) - qi(1,i2)
            dq2 = qi(2,i1) - qi(2,i2)
            dq3 = qi(3,i1) - qi(3,i2)
            
c
            stest(1) = -sqrt(fpi)/vol/ (4.0*lamda*lamda)
            do 30 lm = 2,lmmax
               stest(lm) = 0.0d0
   30       continue
c
c---> exclude the origine and add correction if i1.eq.i2
c
            if (i1.eq.i2) then
               stest(1) = stest(1) - lamda/pi
               nstart = 2
 
            else
               nstart = 1
            end if
c
c---> loop first over n-1 shells of real and recipro. lattice - then
c      add the contribution of the last shells to see convergence
c
            do 40 it = 1,2
               if (it.eq.1) then
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
c
c---> sum over real lattice
c
               do 50 i = nrs,nre
                  r1 = dq1 - rm(1,i)
                  r2 = dq2 - rm(2,i)
                  r3 = dq3 - rm(3,i)
                  call ymy(r1,r2,r3,r,ylm,lmax)
c
                  alpha = lamda*r
c
                  call gamfc(alpha,g,lmax,r)
c
                  do 60 l = 0,lmax
c
                     rfac = g(l)/sqrt(pi)
                     do 70 m = -l,l
                        lm = l* (l+1) + m + 1
                        stest(lm) = stest(lm) + ylm(lm)*rfac
   70                continue
   60             continue
 
   50          continue
c
c---> sum over reciprocal lattice
c
               do 80 i = ngs,nge
                  g1 = gn(1,i)
                  g2 = gn(2,i)
                  g3 = gn(3,i)
                  call ymy(g1,g2,g3,ga,ylm,lmax)
c
                  beta = ga/lamda
                  expbsq = exp(beta*beta/4.0d0)
                  dqdotg = dq1*g1 + dq2*g2 + dq3*g3
c
                  bfac = fpi*exp(ci*dqdotg)/ (ga*ga*expbsq*vol)
                  do 90 l = 0,lmax
c
                     do 100 m = -l,l
                        lm = l* (l+1) + m + 1
                        stest(lm) = stest(lm) + ylm(lm)*bfac
  100                continue
                     bfac = bfac*ga/ci/real(2*l+1)
   90             continue
   80          continue
c
               if (it.eq.1) then
                  do 110 lm = 1,lmmax
                     if (abs(dimag(stest(lm))).gt.bound) then
                        go to 120
 
                     else
c                        smat(i1,i2,lm) = real(stest(lm))
c linux change
                        smat(i1,i2,lm) = stest(lm)
                        stest(lm) = 0.0d0
                     end if
 
  110             continue
 
               else
c
c---> test convergence
c
                  do 130 lm = 1,lmmax
c                     s = real(stest(lm))
c   changed for linux 
                     s = stest(lm)
                     smat(i1,i2,lm) = smat(i1,i2,lm) + s
                     if (abs(s).gt.bound) write (6,fmt=9000) i1,i2,lm,
     +                   abs(s)
  130             continue
               end if
 
   40       continue
c
            if (TEST('madelsum')) THEN 
            do 140 lm = 1,lmmax
               if (abs(smat(i1,i2,lm)).gt.1.0d-8) write (*,fmt=*) i1,i2,
     +             lm,smat(i1,i2,lm)
  140       continue
            end if
c
   20    continue
   10 continue
      return
 
  120 stop ' imaginary contribution to real lattice sum '
 
 
 9000 format (1x,' convergence of smat(',i4,i4,i4,') : ',d10.5,
     +       ' is less than 1.0d-8 - use more lattice vectors ')
 
      end




