      SUBROUTINE WRITESHAPE(NPAN,MESHN,NM,XRN,DRN,NFU,LMIFUN,THETAS,
     &                      ISHAPE)
      use mod_version_info, only: serialnr
      IMPLICIT NONE
c#@# KKRtags: VORONOI input-output shape-functions
c Write out shape function into unit 15.
      include 'inc.geometry'
      INTEGER IBMAXD
      PARAMETER (IBMAXD=(LMAXD1+1)*(LMAXD1+1))
c Input:
      INTEGER NPAN,MESHN,NFU,ISHAPE
      INTEGER NM(NPAND),LMIFUN(IBMAXD)
      REAL*8  DRN(IRID),THETAS(IRID,IBMAXD),XRN(IRID)

c Local:
      INTEGER IPAN1,IR,IFUN
      CHARACTER(len=150) ::char1
      
      ! write serial number after Shape number
      char1 = ';     # serial: ' // trim(serialnr)
      
      WRITE(15,FMT=9020) NPAN,MESHN,ISHAPE,trim(char1)
      WRITE (15,FMT=9000) (NM(IPAN1),IPAN1=1,NPAN)
      WRITE (15,FMT=9010) (XRN(IR),DRN(IR),IR=1,MESHN)
      WRITE (15,FMT=9000) NFU
      DO IFUN = 1,NFU
         WRITE (15,FMT=9000) LMIFUN(IFUN)
         WRITE (15,FMT=9010) (THETAS(IR,IFUN),IR=1,MESHN)
      END DO

      RETURN

 9000 FORMAT (16I5)
 9010 FORMAT (4D20.12)
 9020 FORMAT (2I5,' NPAN,MESHN;  Shape number',I6,A)

      END
