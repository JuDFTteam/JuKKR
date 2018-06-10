SUBROUTINE epathtb(ez,df,efermi,npnt,iesemicore,idosemicore,  &
        ebotval,emuval,tkval,npolval,n1val,n2val,n3val,  &
        ebotsem,emusem,tksem,npolsem,n1sem,n2sem,n3sem,  &
        iemxd)
! **********************************************************************
! *                                                                    *
! * Generating the energy mesh.                                        *
! *                                                                    *
! * Calls the routine EMESHT once for the valence contour and once for *
! * the semicore contour.                                              *
! * In the semicore range, -NPOLSEM is used to create a rectangular    *
! * contour.                                                           *
! *              ph. mavropoulos, v.popescu Juelich/Munich 2004        *
! *                                                                    *
! **********************************************************************
      use mod_types, only: t_inc
      IMPLICIT NONE
      INTEGER IEMXD
      DOUBLE COMPLEX EZ(*),DF(*),EZSEMI(IEMXD),DFSEMI(IEMXD)
      DOUBLE COMPLEX EZVAL(IEMXD),DFVAL(IEMXD)
      DOUBLE PRECISION EBOTSEM,EMUSEM,TKSEM,EBOTVAL,EMUVAL,TKVAL,EFERMI
      INTEGER NPOLSEM,N1SEM,N2SEM,N3SEM
      INTEGER NPOLVAL, N1VAL, N2VAL, N3VAL
      INTEGER IESEMICORE,NPNT,NPNTSEMI,NPNTVAL
      INTEGER IE,JE
      INTEGER IDOSEMICORE


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
IF(t_inc%i_write>0) THEN
  WRITE (1337,*)
  WRITE (1337,'(79(1H=))')
  WRITE (1337,'(20X,A)') 'EPATHTB: generates a complex E contour'
  WRITE (1337,'(79(1H=))')
  WRITE (1337,*)
END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

iesemicore = 0
IF ( idosemicore == 1 ) THEN
  IF(t_inc%i_write>0) WRITE(1337,99001) 'semi-core contour'
  CALL emesht(ezsemi,dfsemi,npntsemi,ebotsem,emusem,efermi,  &
      tksem,-npolsem,n1sem,n2sem,n3sem,iemxd)
  iesemicore = npntsemi
  IF(t_inc%i_write>0) WRITE(1337,99001) 'valence contour'
END IF
CALL emesht(ezval,dfval,npntval,ebotval,emuval,efermi,tkval,  &
    npolval,n1val,n2val,n3val,iemxd)

npnt = iesemicore + npntval

DO ie = 1,iesemicore
  ez(ie) = ezsemi(ie)
  df(ie) = dfsemi(ie)
END DO

DO ie = iesemicore+1,npnt
  je = ie - iesemicore
  ez(ie) = ezval(je)
  df(ie) = dfval(je)
END DO

99001 FORMAT(7X,'* ',a,/,7X,20(1H-),/)
END SUBROUTINE epathtb
