c ************************************************************************
      SUBROUTINE SNAME(NAME,NEW,BAND)
C ************************************************************************
C     .. scalar arguments
c
      INTEGER BAND
      CHARACTER*40 NAME,NEW

c     .. locals
c
      INTEGER I,L,LO
      CHARACTER*1  CH(50),POI
      CHARACTER*10  S

      INTEGER LENGTH
      EXTERNAL LENGTH

C ------------------------------------------------------------------------
      POI='.'
      if (band.lt.0) then
        lo=log(real(-band))/log(10.0d0)+1
      else if (band.eq.0) then 
        lo=0
      else 
        LO=LOG(REAL(BAND))/LOG(10.0D0)
      end if

c      write(6,*) 'LO ',lo

      READ(NAME,fmt='(255a1)')(CH(I),I=1,40)
      L = LENGTH(CH,40)
c      write(6,*) 'L  ',l

c      write(6,*) 'CH ',(CH(I),I=1,25)

      WRITE(S,FMT='(I10)') BAND
c      write(6,*) 'S  ',s

      READ(S,FMT='(255A1)') (CH(I),I=L+1,L+10)
c      write(6,*) 'CH ',(CH(I),I=L+1,L+10)

      WRITE(NEW,FMT='(255A1)') 
     +     (CH(I),I=1,L),POI,(CH(I),I=L+10-LO,L+10)

      RETURN
      END
