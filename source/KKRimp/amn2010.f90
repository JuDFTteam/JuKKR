!------------------------------------------------------------------------------------
!> Summary: Module handling the structure constants for the intersite potential
!> Author:
!> For details see vinters2010; also used to construct the intercell potential of the host
!------------------------------------------------------------------------------------
module mod_amn2010
REAL*8,allocatable                           ::     CLEB(:,:)
INTEGER,allocatable                          ::     ICLEB(:,:)
integer                                      ::     iend,ncleb
integer                                      ::     first=1,amnlmaxpot

contains

!-------------------------------------------------------------------------------
!> Summary: Structure constants for the intersite potential
!> Author:
!> Category: electrostatics, potential, KKRimp
!> Deprecated: False 
!> For details see vinters2010
!-------------------------------------------------------------------------------
subroutine amn2010(jatom,amat,bmat,alat,gauntcoeff,lpotmax,ntotatom,ratom)
use type_gauntcoeff
use mod_ymy
use mod_amngaunt
implicit none
! integer lmmaxd,l3d,lm3d
! PARAMETER (LMMAXD= (LPOTD+1)**2,L3D=2*LPOTD,LM3D= (L3D+1)**2)

!interface
integer                               ::     jatom
real*8                                ::     amat(ntotatom,(lpotmax+1)**2,(lpotmax+1)**2)
real*8                                ::     bmat(ntotatom,(lpotmax+1)**2)
real*8                                ::     alat
type(gauntcoeff_type),intent(in)      ::     gauntcoeff
integer                               ::     lpotmax
integer ntotatom
real*8                                ::     ratom(3,ntotatom)
!local 
integer                               ::     iatom
integer                               ::     lmpotmax
real*8                                ::     epi,fpi
real*8                                ::     rlength,r1,r2,r3
real*8                                ::     dfac(0:lpotmax,0:lpotmax)
integer                               ::     lx, ly,ival
integer                               ::     lm1,lm3,lm4,l1,l2
real*8                                ::     y((2*lpotmax+1)**2)

! pi    = 3.14159265359D0
fpi   = 4.0d0*pi
epi   = 8.0d0*pi
lmpotmax = (lpotmax+1)**2
amat     = 0.0d0
bmat     = 0.0d0

ncleb=(2*lpotmax+1)**2 * (lpotmax+1)**2 
IF (NCLEB.LT.gauntcoeff%IEND) THEN
  STOP 'Error NCLEB < IEND'
END IF

if (first==1) then
  amnlmaxpot=lpotmax
  allocate(CLEB(NCLEB,2),ICLEB(NCLEB,4))

!   write(*,*) ncleb
!   stop
    cleb=0.0D0
     icleb = 0.0D0
!   write(20000,*) LpotMAX,CLEB,ICLEB,IEND,gauntcoeff%Wg,gauntcoeff%YRG
  CALL AMNGAUNT(LpotMAX,CLEB,ICLEB,IEND,gauntcoeff%Wg,gauntcoeff%YRG,2*LpotMAX,2*LpotMAX,ncleb,(2*LpotMAX+1)**2)
!   write(20001,*) LpotMAX,CLEB,ICLEB,IEND,gauntcoeff%Wg,gauntcoeff%YRG
!   write(*,*) CLEB
!   write(*,*) ICLEB
!   stop
  first=0
else 
  if (amnlmaxpot/=lpotmax) then
    write(*,*) '[AMN2010] Warning: Lpotmax for amngaunt has changed'
    deallocate(CLEB,ICLEB)
    allocate(CLEB(NCLEB,2),ICLEB(NCLEB,4))
    cleb=0.0D0
    icleb = 0.0D0

    CALL AMNGAUNT(LpotMAX,CLEB,ICLEB,IEND,gauntcoeff%Wg,gauntcoeff%YRG,2*LpotMAX,2*LpotMAX,ncleb,(2*LpotMAX+1)**2)
    amnlmaxpot=lpotmax
  end if
end if


!
!--->calculate:                  (2*(l+l')-1)!!
!                dfac(l,l')= ----------------------
!                            (2*l+1)!! * (2*l'+1)!!
dfac(0,0) = 1.d0
do lx = 1,lpotmax
  dfac(lx,0) = dfac(lx-1,0)*real(2*lx-1)/real(2*lx+1)
  dfac(0,lx) = dfac(lx,0)
  do ly = 1,lx
    dfac(lx,ly) = dfac(lx,ly-1)*real(2* (lx+ly)-1)/real(2*ly+1)
    dfac(ly,lx) = dfac(lx,ly)
  end do !ly
end do !lx

do iatom=1,ntotatom
  if (jatom.ne.iatom) then
    r1 = ratom(1,jatom) - ratom(1,iatom)
    r2 = ratom(2,jatom) - ratom(2,iatom)
    r3 = ratom(3,jatom) - ratom(3,iatom)
    call ymy(r1,r2,r3,rlength,y,2*lpotmax)
    rlength = rlength*alat
    do lm1 = 1,lmpotmax
!       write(*,*) lm1
      l1 = gauntcoeff%loflm(lm1)
      bmat(iatom,lm1) = bmat(iatom,lm1) - epi/real(2*l1+1)*y(lm1)/ (rlength** (l1+1))
    end do
!        write(*,*) iend
!       stop
    do ival = 1,iend
      lm1 = icleb(ival,1)
      lm4 = icleb(ival,2)
      lm3 = icleb(ival,3)
!       write(*,*) lm1,lm4
      l1  = gauntcoeff%loflm(lm1)
      l2  = gauntcoeff%loflm(lm4)
!       write(*,*) lm1,lm4,lm3,l1,l2
      amat(iatom,lm1,lm4) = amat(iatom,lm1,lm4) + epi*fpi*dfac(l1,l2)*y(lm3) * cleb(ival,1)/(rlength** (l1+l2+1))
!       write(40000,*) lm1,lm4,lm3,l1,l2
!       write(40000,*) epi,fpi,dfac(l1,l2),y(lm3),cleb(ival,1),(rlength** (l1+l2+1)),amat(iatom,lm1,lm4)
    end do
  else
!      do lm1=1,(lpotmax+1)**2
!        bmat(iatom,lm1)=0.0d0
!        do lm3=1,(lpotmax+1)**2
!          amat(iatom,lm1,lm3)=0.0d0
!        end do !lm3=1,(4*lmax+1)**2
!      end do !lm1=1,(4*lmax+1)**2
  end if
end do !iatom=1,ntotatom
! amat=0.0D0
! bmat=0.0D0
! stop
! write(*,*) amat
end subroutine amn2010
end module mod_amn2010
