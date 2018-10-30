      SUBROUTINE READ_DELTA_TMAT(DTMATLM,DELTAMAT,LMCL,NCL)

      implicit none

      INTEGER,INTENT(IN)     ::  LMCL,NCL
      COMPLEX,INTENT(OUT)    ::  DTMATLM(LMCL,LMCL),
     +                           DELTAMAT(LMCL,LMCL)
      DOUBLE PRECISION       ::  RCL(3,NCL)

      integer                ::  LM1,LM2,LM1D,LM2D,NCLT,STATUS_OPEN

      DTMATLM=0.d0
      DELTAMAT=0.d0
      RCL=0d0

      open (unit=56, file="DTMTRX",form="formatted",
     +              action="read",iostat=status_open)

      if (status_open /= 0) then
         stop 'DELTA_TMAT does not exist'
      end if

      READ(56, "(I5)") NCLT
      IF (NCLT .NE. NCL) STOP "Number of cluster atoms"
      DO LM1=1,NCL
        READ(56,"(3e17.9)") (RCL(LM2,LM1),LM2=1,3)
      END DO

      DO LM1=1,LMCL
        DO LM2=1,LMCL
          READ(56,"((2I5),(4e17.9))") LM2D,LM1D,DTMATLM(LM2,LM1),
     +                                       DELTAMAT(LM2,LM1)
        END DO
      END DO

      CLOSE(56)

      END SUBROUTINE READ_DELTA_TMAT
