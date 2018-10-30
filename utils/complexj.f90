program complexj
implicit none
! Extrapolate the integrand of the Heisenberg exchange coupling to the real axis
! run: complexj.exe < infile > outfile

integer ne ! Number of energies
complex*16,allocatable:: ez(:),jz(:),djz(:) ! Energies, exchange
real*8 jreal  ! interpolated value of exchange at real axis
real*8 er,ei,jr,ji
complex*16 de  ! distance between energy points
complex*16 eps ! Imaginary part of energy
real*8 pi
integer ie
character*80 text1
character*11 text2

read(*,*) text1
read(*,*) text1
read(*,*) text2,ne

allocate(ez(ne),jz(ne),djz(ne))

do ie = 1,ne
   read(*,*) er,ei,jr,ji
   ez(ie) = dcmplx(er,ei)
   jz(ie) = dcmplx(jr,ji)
enddo

de = ez(2) - ez(1)
eps = dcmplx(0.d0,dimag(ez(1)))

! Derivative
djz(1) = (jz(2) - jz(1))/de
djz(ne) = (jz(ne) - jz(ne-1))/de
do ie = 2,ne - 1
   djz(ie) = (jz(ie+1) - jz(ie-1))/de
enddo

pi = 4.d0*datan(1.d0)
do ie = 1,ne
jreal = dimag( jz(ie) - eps*djz(ie) )    ! Taylor to the real axis
write(*,fmt='(2e12.4)') dreal(ez(ie)),jreal
enddo

end

