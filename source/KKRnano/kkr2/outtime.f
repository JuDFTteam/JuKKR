      SUBROUTINE OUTTIME(RANK,NAME,TIME_I,ITER)

c**********************************************************************
c     every time this subroutine is called it produces one line of
c     output on the file with the unit-number 2, which consitst
c     of up to 59 characters (from the variable name) and the time.
c     here the time must be in seconds. in the output the time is
c     given in seconds, but also in hours minutes and seconds.
c                                              p.kurz   8.2.96
c**********************************************************************

      IMPLICIT NONE

C     .. Scalar Arguments ..
      REAL          TIME_I
      INTEGER       ITER,RANK
      CHARACTER*(*) NAME
C     ..
C     .. Local Scalars ..
      REAL          TIME_F,REST,SECONDS
      INTEGER       ihours,iminutes
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC     real,int
C     ..
C     only one processor proceeds ...
C
      IF (RANK.EQ.0) THEN
C
      CALL CPU_TIME(TIME_F)
C
C     calculate time in hours, minutes and seconds
C
      REST = TIME_F - TIME_I

      ihours = int(rest/3600.0)
      rest = rest - real(ihours)*3600

      iminutes = int(rest/60.0)
      seconds = rest - real(iminutes)*60

c     output of the results

      WRITE (2,FMT=8000) 
     + iter,name,(TIME_F-TIME_I),ihours,iminutes,seconds
      WRITE (6,FMT=8000) 
     + iter,name,(TIME_F-TIME_I),ihours,iminutes,seconds
C
      ENDIF
C
 8000 FORMAT ('ITER: ',i3,2x,a,t42,f9.2,
     +        ' sec = ',i3,' h ',i2,' min ',f5.2,' sec')
C
      END SUBROUTINE OUTTIME
