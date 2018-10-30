program scaledos
! Change dos units:
! mode=1 : From n(E) to n(E-EF)
! mode=2 : From n(E) to n(E-EF) and from ryd to eV
implicit none
integer nlm ! How many lm components
integer ne  ! How many energies
real*8 ev
parameter(ev=13.6058)
integer mode ! see above
real*8 ef    ! Fermi level
real*8 energy
integer ie,ilm,ios
character*100 filename
real*8,allocatable :: dos(:)

open(11,file='scaled.dos')

10 write(*,*) "Original filename? (max. 100 characters)"
read(*,*) filename
open(10,file=filename,form='formatted',iostat=ios)
if (ios.ne.0) then
   write(6,*) ' I cannot find file >',filename
   goto 10
endif

20 write(*,*) "Operation mode?"
write(*,*) "1: From n(E) to n(E-EF)"
write(*,*) "2: From n(E) to 0.5*n(E-EF)"
write(*,*) "3: From n(E) to n(E-EF) and from ryd to eV"
write(*,*) "4: From n(E) to 0.5*n(E-EF) and from ryd to eV"
read(*,*) mode
if (mode.ne.1.and.mode.ne.2.and.mode.ne.3.and.mode.ne.4) then
   write(*,*) "No such mode:",mode
   goto 20
endif

write(*,*) "Fermi level (ryd)?"
read(*,*) ef

write(*,*) "How many DOS components per energy?"
read(*,*) nlm
allocate( dos(nlm) )

ne=0
ios=0
do while (ios.eq.0)
   read(10,*,END=100) energy,(dos(ilm),ilm=1,nlm)
   ne = ne + 1
   energy = energy - ef
   if (mode.eq.3) then
      energy = energy * ev
      dos(1:nlm) = dos(1:nlm) / ev
   endif
   if (mode.eq.2.or.mode.eq.4) dos(1:nlm) = dos(1:nlm) / 2.d0

   write(11,fmt='(50e16.8)') energy, (dos(ilm),ilm=1,nlm)
enddo

100 write(*,*) 'End of file; energies:',ne

write(*,*) 'Output written in file scaled.dos'


end
