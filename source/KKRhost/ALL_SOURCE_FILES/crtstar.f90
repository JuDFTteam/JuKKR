! 20.07.96 ***************************************************************
SUBROUTINE crtstar(ratom,nshell,nd,irot,isymindex,rrot)
! ************************************************************************
!  THE SYMMETRY OPERATIONS OF THE SYMMETRY GROUP ARE APPLIED TO THE
!  INPUT VECTOR RATOM
! ------------------------------------------------------------------------
implicit none
INTEGER IROT,NSHELL
DOUBLE PRECISION ND(64,3,*),RATOM(3,*),RROT(48,3,*)
INTEGER ISYMINDEX(*)

INTEGER I,ID,NS,K,J,ISYM
LOGICAL TEST
EXTERNAL TEST
! ------------------------------------------------------------------------

DO  ns = 1,nshell
  DO  id = 1,irot
    isym = isymindex(id)
    DO  i = 1,3
      rrot(id,i,ns) = nd(isym,i,1)*ratom(1,ns) + nd(isym,i,2)*ratom(2,ns) +  &
          nd(isym,i,3)*ratom(3,ns)
    END DO
  END DO
END DO

IF (test('ND      ')) WRITE(1337,FMT='((I3,3(/,3f6.2)))')  &
    (k,((nd(k,i,j),j=1,3), i=1,3),k=1,irot)

RETURN

END SUBROUTINE crtstar
