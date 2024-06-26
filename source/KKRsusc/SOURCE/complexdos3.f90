module mod_complexdos3
! Interpolates the complex density of states to the real axes
! does the same thing as compexdos3.f from the Juelich-Muenchen 
! source
! Principle of DOS here: Two-point contour integration
! for DOS in the middle of the two points. The input DOS
! and energy must be complex. Parameter deltae should be
! of the order of magnitude of eim. 
!
!    
!      <-2*deltae->   _
!           /\        |     DOS=(n(1)+n(2))/2 + (n(1)-n(2))*eim/deltae
!          /  \       |
!        (1)  (2)   2*i*eim=2*i*pi*Kb*Tk
!        /      \     |
!       /        \    |
!------------------------ (Real E axis)
contains

subroutine complexdos3(lmax,iemax,iatom,ispin,nspin,dos,doslm,ez)
implicit none
integer       :: lmax
integer       :: iemax
integer       :: iatom
integer       :: ispin
integer       :: nspin
double complex :: dos(0:lmax+2,nspin,iemax)
double complex :: doslm((lmax+1)**2,nspin,iemax)
double complex :: ez(iemax)
!local
double precision :: eim,deltae,enew !tk,kb,pi
integer          :: ie,lval,lmmax
double complex :: dosnew(0:lmax+2)
double complex :: doslmnew(1:(lmax+1)**2)
integer,parameter        :: pi=3.1415926535897931d0

lmmax=(lmax+1)**2

open(unit=30, &
     file="out_ldos.interpol.atom="//char(48+iatom/10)//char(48+mod(iatom,10))//"_spin"//char(48+ispin)//".dat")

open(unit=31, &
     file="out_lmdos.interpol.atom="//char(48+iatom/10)//char(48+mod(iatom,10))//"_spin"//char(48+ispin)//".dat")


do ie=2,iemax-1
  deltae = dreal(ez(ie+1) - ez(ie))
  eim = dimag(ez(ie))
  enew = dreal(ez(ie)) ! real quantity

  do lval=0,lmax+2
      dosnew(lval) = dos(lval,ispin,ie) &
        + 0.5d0*(dos(lval,ispin,ie-1)-dos(lval,ispin,ie+1))*dcmplx(0.d0,eim)/deltae
  enddo

  do lval=1,lmmax
      doslmnew(lval) = doslm(lval,ispin,ie) &
        + 0.5d0*(doslm(lval,ispin,ie-1)-doslm(lval,ispin,ie+1))*dcmplx(0.d0,eim)/deltae
  enddo


  write(30,'(300g24.16)') enew,-dimag(dosnew(lmax+2))/pi,    & ! e, l-summed dos,
                        (-dimag(dosnew(lval))/pi,lval=0,lmax+1)       ! l-resolved dos
  write(31,'(300g24.16)') enew,-dimag(doslmnew(lmax+2))/pi,    & ! e, l-summed dos,
                        (-dimag(doslmnew(lval))/pi,lval=1,lmmax)       ! l-resolved dos

!   write(50,9001) enew,dimag(dosnew(lmax+1)) &
!                ,(dimag(dosnew(lval)),lval=0,lmax)

enddo ! ie=2,iemax-1

close(30)
close(31)
end subroutine complexdos3


end module mod_complexdos3
