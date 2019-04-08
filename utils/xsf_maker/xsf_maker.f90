program xsf_maker

implicit none

double precision, parameter   :: a0 = 5.2917721067e-1

logical, parameter :: lprint = .false.
character(len=100) :: filename
logical            :: langstrom
double precision   :: alat, bravais(3,3), rshift(3)
integer            :: natom
integer,          allocatable :: zatom(:)
double precision, allocatable :: rbasis(:,:)
integer          :: i
double precision :: veclen

!-----------------------------------------------------------------------
! output file name
read(*,*) ! comment line
read(*,*) filename
filename = trim(adjustl(filename))//'.xsf'
if (lprint) write(*,'("filename=",(a))') filename
!-----------------------------------------------------------------------
read(*,*) ! comment line
read(*,*) alat, langstrom
if (langstrom) alat = alat*a0
if (lprint) write(*,'("alat=",f16.8)') alat
!-----------------------------------------------------------------------
read(*,*) ! comment line
do i=1,3
  read(*,*) bravais(1:3,i)
  bravais(1:3,i) = bravais(1:3,i)*alat
  veclen = sqrt(sum(bravais(1:3,i)**2))
  if (veclen < 1.d-6) stop 'bravais vector has zero length!'
  if (lprint) write(*,'("bravais=",3f16.8)') bravais(1:3,i)
end do
!-----------------------------------------------------------------------
read(*,*) ! comment line
read(*,*) rshift(1:3)
if (lprint) write(*,'("rshift=",3f16.8)') rshift(1:3)
!-----------------------------------------------------------------------
read(*,*) ! comment line
read(*,*) natom
if (lprint) write(*,'("natom=",i8)') natom
allocate(zatom(natom),rbasis(3,natom))
read(*,*) ! comment line
do i=1,natom
  read(*,*) zatom(i), rbasis(1:3,i)
  rbasis(1:3,i) = (rbasis(1:3,i) + rshift(1:3))*alat
  if (lprint) write(*,'("zatom=",i4,"  rbasis=",3f16.8)')  zatom(i), rbasis(1:3,i)
end do
!-----------------------------------------------------------------------
! write xsf file
open(file=filename,unit=10,status='replace')
write(10,'("CRYSTAL")')
write(10,'("PRIMVEC")')
do i=1,3
  write(10,'(3f12.8)') bravais(1:3,i)
end do
write(10,'("CONVVEC")')
do i=1,3
  write(10,'(3f12.8)') bravais(1:3,i)
end do
write(10,'("PRIMCOORD")')
write(10,'(2i8)') natom, 1
do i=1,natom
  write(10,'(i4,3f12.8)') zatom(i), rbasis(1:3,i)
end do
close(10)
!-----------------------------------------------------------------------
deallocate(zatom,rbasis)
! Done!
end program xsf_maker
