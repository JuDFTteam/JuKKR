PROGRAM complexdos
IMPLICIT NONE
! Principle of DOS here: Two-point contour integration
! for DOS in the middle of the two points. The input DOS
! and energy must be complex. Parameter deltae should be
! of the order of magnitude of eim.


!      <-2*deltae->   _
!           /\        |     DOS=(n(1)+n(2))/2 + (n(1)-n(2))*eim/deltae
!          /  \       |
!        (1)  (2)   2*i*eim=2*i*pi*Kb*Tk
!        /      \     |
!       /        \    |
!------------------------ (Real E axis)
INTEGER*4 iemaxd,lmaxd
PARAMETER(iemaxd=1000,lmaxd=10)
INTEGER*4 npot,iemax,lmax
REAL*8 eim,deltae,tk,kb,pi
! Total dos stored in DOS(LMAX+1,IE)
COMPLEX*16 dos(0:lmaxd+1,iemaxd),ez(iemaxd)
COMPLEX*16 dosnew(0:lmaxd+1)
REAL*8 temp(2*lmaxd+6)
REAL*8 enew,ef,ev
INTEGER*4 ie,ii,ll,iheader
CHARACTER (LEN=80) :: text
CHARACTER (LEN=45) :: text1
CHARACTER (LEN=20) :: text2

ev = 13.6058

! If only Tk is known, use Bolzmann constant, pi to find eim:
!     Kb=0.6333659D-5
!     pi=3.14159265358979312d0
!     eim=pi*Kb*Tk
OPEN (49,FILE='complex.dos',FORM='formatted')
OPEN (50,FILE='new3.dos',FORM='formatted')
OPEN (51,FILE='new3_eV_EF.dos',FORM='formatted')
READ (49,*) text2 !dummy readin of header, may be replaced later
READ (49,*) npot
READ (49,*) iemax
READ (49,*) lmax
IF (iemax > iemaxd) STOP 'IEMAX.gt.IEMAXD'
IF (lmax > lmaxd) STOP 'LMAX.gt.LMAXD'


DO ii=1,npot
  WRITE(*,*) 'Reading potential',ii
! Read header:
  DO iheader=1,3
    READ(49,9002) text
    WRITE(50,9002) text
    WRITE(51,9002) text
    IF (iheader == 1) WRITE(*,9002) text
  END DO
  READ(49,FMT='(A45,F10.6,A20)') text1,ef,text2
  WRITE(50,FMT='(A45,F10.6,A20)') text1,ef,text2
  WRITE(51,FMT='(A45,F10.6,A20)') text1,ef,text2
  DO iheader=5,9
    READ(49,9002) text
    WRITE(50,9002) text
    WRITE(51,9002) text
  END DO
! Read dos: (total dos stored at DOS(LMAX+1,IE))
  DO ie=1,iemax
    READ(49,*) ez(ie),(dos(ll,ie),ll=0,lmax+1)
  END DO
! Compute and write out corrected dos at new (middle) energy points:
  DO ie=2,iemax-1
    deltae = dreal(ez(ie+1) - ez(ie))
    eim = DIMAG(ez(ie))
    enew = dreal(ez(ie)) ! Real quantity
    
    DO ll=0,lmax+1
      dosnew(ll) = dos(ll,ie)  &
          + 0.5D0*(dos(ll,ie-1)-dos(ll,ie+1))*DCMPLX(0.d0,eim)/deltae
    END DO
    
    WRITE(50,9001) enew,DIMAG(dosnew(lmax+1)) ,(DIMAG(dosnew(ll)),ll=0,lmax)
    WRITE(51,9001) (enew-ef)*ev,DIMAG(dosnew(lmax+1))/ev  &
        ,(DIMAG(dosnew(ll))/ev,ll=0,lmax)
    
  END DO ! IE=2,IEMAX-1
  
  IF (ii /= npot) THEN
    READ(49,*)
    WRITE(50,9002) ' '
    WRITE(51,9002) ' '
  END IF
  
END DO ! II=1,NPOT

CLOSE(49)
CLOSE(50)
CLOSE(51)

STOP 'telos'
9000 FORMAT(10E12.4)
9001 FORMAT(10E14.6)
9002 FORMAT(a80)
9003 FORMAT(8(2E12.4))
END PROGRAM complexdos

