      program complexdos
      implicit none
c Principle of DOS here: Two-point contour integration
c for DOS in the middle of the two points. The input DOS
c and energy must be complex. Parameter deltae should be
c of the order of magnitude of eim. 
c
c    
c      <-2*deltae->   _
c           /\        |     DOS=(n(1)+n(2))/2 + (n(1)-n(2))*eim/deltae
c          /  \       |
c        (1)  (2)   2*i*eim=2*i*pi*Kb*Tk
c        /      \     |
c       /        \    |
c------------------------ (Real E axis)
      integer*4 IEMAXD,LMAXD
      parameter(IEMAXD=1000,LMAXD=10)
      integer*4 NPOT,IEMAX,LMAX
      real*8 eim,deltae,Tk,Kb,pi
c Total dos stored in DOS(LMAX+1,IE)
      complex*16 DOS(0:LMAXD+1,IEMAXD),EZ(IEMAXD)
      complex*16 DOSNEW(0:LMAXD+1)
      real*8 TEMP(2*LMAXD+6)
      real*8 ENEW
      integer*4 IE,II,LL,IHEADER
      character*80 TEXT

c If only Tk is known, use Bolzmann constant, pi to find eim:
c     Kb=0.6333659D-5
c     pi=3.14159265358979312d0
c     eim=pi*Kb*Tk
      OPEN (49,FILE='complex.dos',FORM='formatted')
      OPEN (50,FILE='new.dos',FORM='formatted')
      READ (49,*) NPOT
      READ (49,*) IEMAX
      READ (49,*) LMAX
      if (IEMAX.gt.IEMAXD) stop 'IEMAX.gt.IEMAXD'
      if (LMAX.gt.LMAXD) stop 'LMAX.gt.LMAXD'


      do II=1,NPOT
c Read header:
         do IHEADER=1,9
            read(49,9002) TEXT
            write(50,9002) TEXT
         enddo
c Read dos: (total dos stored at DOS(LMAX+1,IE))
         do IE=1,IEMAX
            read(49,*) EZ(IE),(DOS(LL,IE),LL=0,LMAX+1)
         enddo
C
c Compute and write out corrected dos at new (middle) energy points:
C
         do IE=1,IEMAX-1
            DELTAE = DREAL(EZ(IE+1) - EZ(IE))
            if (DELTAE.GT.1D-4) then
               EIM = DIMAG(EZ(IE))
               ENEW = 0.5d0*DREAL(EZ(IE+1) + EZ(IE)) ! Real quantity

               do LL=0,LMAX+1
                  DOSNEW(LL) = 0.5d0*(DOS(LL,IE)+DOS(LL,IE+1))
     &                 + (DOS(LL,IE)-DOS(LL,IE+1))
     &                 *DCMPLX(0.d0,EIM)/DELTAE
               enddo
         
               write(50,9001) ENEW,(DIMAG(DOSNEW(LL)),LL=0,LMAX+1)
            end if

         enddo                  ! IE=1,IEMAX-1

         if (II.ne.NPOT) then 
            read(49,*) TEXT
            write(50,9002) ' '
         endif

      enddo                     ! II=1,NPOT
      
      close(49)
      close(50)
      
      stop 'telos'
 9000 format(10e12.4)
 9001 format(10e14.6)
 9002 format(A80)
 9003 format(8(2e12.4))
      end

