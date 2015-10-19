! module mod_clustcomp
! contains
!   include 'clustcomp.f90'
! ! end module mod_clustcomp
! ! 
! 
! ! module mod_dsort
! ! contains
!   include 'dsort.f90'
! ! end module mod_dsort
! 
! ! module mod_clsgen_imp
! ! contains
!   include 'clsgen_imp.f90'
! ! end module mod_clsgen_imp
! 
! 
! 
include  'nrtype'
! include 'mathtools.f90'
program test
 use nrtype
!  use mod_mathtools
integer,parameter           :: n=3
complex(kind=dpc)           :: A(3,3)
complex(kind=dpc)           :: b(3,1)

A=0
A(1,1)=2
A(2,2)=3
A(3,3)=4
b=0
b(3,1)=1
! b(2,2)=1
! b(3,3)=1
write(*,*) A(:,1:2)
! A=(/ (/ (1,1),(2,2),(3,3) /), (/ (4,4),(5,5),(6,6) /) , (/ (7,7),(8,8),(9,9) /) /)
! b=(/ (/ 1,0,0 /), (/0,1,0/),(/ 0,0,1 /) /)
! do i = 1, n
! write(*,*) b(i,:)
! end do


! call linearsolve(A,b)

! do i = 1, n
! write(*,*) b(i,:)
! end do
! use mod_dsort
! use mod_clustcomp
! use mod_clsgen
! implicit none
! 
!   integer,parameter  ::  ntotatom=10
!   real*8             ::  rcls(3,ntotatom,ntotatom)
!   real*8             ::  rcut
!   real*8             ::  rr(3,ntotatom)
!   integer            ::  nr
!   real*8             ::  rmt(ntotatom)
!   integer            :: atom(ntotatom,ntotatom)
!   integer            :: naclsd,nclsd
!   real*8             ::  z(ntotatom)
!   integer             ::  cls(ntotatom)
! 
!   integer             :: iatom,natom
!   integer             :: nacls(ntotatom)
! 
! naclsd=ntotatom
! nclsd=ntotatom
! nr=ntotatom
! 
! rcut=2
! 
! do iatom=1,ntotatom
! 
! rr(:,iatom)=(/ 1.0D0*iatom, 0, 0 /)
! z(iatom)=1
! rmt(iatom)=2.0D0
! natom=ntotatom-1
! end do
! rmt(3)=2.1D0
! write(*,*) 'rr',rr
! write(*,*) 'rmt',rmt
! 
! 
! !       SUBROUTINE CLSGEN(NATOM,RR,NR,Z,RMT,CLS,NACLS, &
! !                          ATOM,RCLS, RCUT,NACLSD,NCLSD)
! 
!   call          CLSGEN(NATOM,RR,NR,RMT,CLS,NACLS,ATOM,RCLS, RCUT,NACLSD,NCLSD)
! !  call         CLSGEN99(RR,NR,RMT,Z,CLS,NACLS,ATOM,RCLS, RCUT,NACLSD,NCLSD)
! 


end program test