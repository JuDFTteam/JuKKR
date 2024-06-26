      program shiftpot
      implicit none
c Shift a potential by a constant
      character*100 text1
      character*100 pot
      integer iline,ir,irm,inew,ifp,ios,ncore
      real*8 vpot(10000),vshift


      write(*,*) 'Filename of original potential (max 100 char)'
      read(*,*) pot
      open(12,file=pot,form='formatted',iostat=ios)
      if (ios.ne.0) then
         write(6,*) ' I cannot find file >',pot
         stop
      end if
      write(*,*) 'Type in 0 for ASA, 1 for full-potential'
      read(*,*) ifp
      write(*,*) 'vshift in Ryd?'
      read(*,*) vshift
      if(ifp.ne.0) then
         write(*,*) 'Attention: Convolution with shapes not supported.'
         write(*,*) 'Full shift from ir=1 to irm by',vshift
      endif


      open(11,file='shifted.pot')

      do iline = 1,4
      read(12,fmt='(a100)') text1
      write(11,fmt='(a100)') text1
      enddo

      read(12,*) irm
      write(11,fmt='(i3)') irm

      read(12,fmt='(a100)') text1
      write(11,fmt='(a100)') text1

      read(12,fmt='(i2,a100)') ncore,text1
      write(11,fmt='(i2,a100)') ncore,text1

      do iline = 1,ncore
         read(12,fmt='(a100)') text1
         write(11,fmt='(a100)') text1
      enddo

      if (ifp.ne.0) then
         read(12,fmt='(a100)') text1
         write(11,fmt='(a100)') text1
      endif

      read(12,fmt='(1p,4d20.13)') (vpot(ir),ir=1,irm)
      write(11,fmt='(1p,4d20.13)') (vpot(ir)+vshift,ir=1,irm)

      iline = 0

      do while (ios.eq.0)

         read(12,fmt='(a100)',END=100) text1
         write(11,fmt='(a100)') text1
         
      end do 
      
 100  write(6,*) ' **** Output written in shifted.pot'

      end
