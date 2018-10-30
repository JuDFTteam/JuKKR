      subroutine findreals(darry,narry,yes)

      implicit none
      integer narry
      double precision darry(narry)
      logical yes

      integer i1,i2,nmax
      parameter (nmax=500)
      integer idone(nmax)

      integer nsqr,nmul
      parameter(nsqr=4,nmul=5)
      integer isqr(nsqr),imul(nmul)
      integer div,divmax
      parameter(divmax=6)
      double precision tol
      parameter (tol = 1d-4)
      double precision dsq,x,xn
      
      data isqr / 2, 3, 5, 7 /
      data imul / 3, 7, 11, 13, 17 /

      do i1=1,nmax
         idone(i1) = 0
      end do
      yes = .false.
      do i1 = 1,nsqr
         do div = 1, divmax
            dsq = dsqrt(dble(div*div*isqr(i1)))
            do i2 = 1, narry
               if (idone(i2).eq.0) then
                  x = darry(i2) * dsq
                  xn=dnint(x)
                  if (dabs(x-xn)/dsq.lt.tol.and.xn.ne.0.d0) then
C                     write(6,300) dabs(darry(i2)),isqr(i1)
C     &                    ,iabs(idnint(xn)),iabs(isqr(i1)*div)
                     darry(i2)=xn/dsq
                     idone(i2)=1
                     yes = .true.
                  endif
               end if
            enddo
         end do
      end do
C
      do i1=1,nmul
         do div=1,divmax
            dsq = dble(div*imul(i1))
            do i2=1,narry
               if (idone(i2).eq.0) then
                  x = darry(i2) * dsq
                  xn = dnint(x)
                  if (dabs(x-xn)/dsq.lt.tol.and.xn.ne.0.d0) then
C                     write(6,301) dabs(darry(i2)),imul(i1)
C     &                    ,iabs(idnint(xn)),div
                     darry(i2)=xn/dsq
                     idone(i2)=1
                     yes = .true.
                  endif
               end if
            enddo
         end do
      end do
                     

 300  format(' FINDSQ: recognize ',f11.8,' as dsqrt(',i1,')*',
     .     i3,'/',i3)
 301  format(' FINDSQ: recognize ',f11.8,' as 1/',i2,' * ',
     .     i2,'/',i1)
            
      return
      end

      

