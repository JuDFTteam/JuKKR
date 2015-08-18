subroutine writehoststructure(bravais,nrbasis,rbasis,NAEZD,NEMBD)
implicit none
! INCLUDE 'inc.p'
!interface
double precision     :: bravais(3,3)
integer              :: nrbasis
double precision     :: RBASIS(3,NAEZD+NEMBD)
integer              :: NAEZD
integer              :: NEMBD

!local
integer              :: iatom
real(kind=8)         :: wght
CHARACTER*256        :: UIO
integer              :: ier
integer              :: itemp1(12),it

open(unit=3463453,file='kkrflex_hoststructure.dat')

write(3463453,'(100A)') '[bravais]'
write(3463453,'(100A)') '#   x   y   z - component'
do iatom=1,3
  write(3463453,'(3F)') BRAVAIS(:,iatom)
end do !nbasis
write(3463453,'(100A)') '[basis]'
write(3463453,*) nrbasis
write(3463453,'(100A)') '#   x   y   z    weight'
do iatom=1,nrbasis
  CALL IoInput('ATOMINFO        ',UIO,iatom+1,7,IER)
  IF (IER.EQ.0) THEN 
     READ (UNIT=UIO,FMT=*)  (itemp1(it),it=1,12),WGHT
  ELSE
     wght = 1.d0
  ENDIF
  write(3463453,'(4F)') RBASIS(:,iatom),wght
end do !nbasis

close(3463453)



end subroutine writehoststructure
