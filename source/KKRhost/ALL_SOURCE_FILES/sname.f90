! ************************************************************************
SUBROUTINE sname(NAME,NEW,band)
! ************************************************************************
!.. scalar arguments
INTEGER BAND
CHARACTER*40 NAME,NEW

!.. locals
INTEGER I,L,LO
CHARACTER*1  CH(50),POI
CHARACTER*10  S

INTEGER LENGTH
EXTERNAL LENGTH
! ------------------------------------------------------------------------
poi='.'
IF (band < 0) THEN
  lo=LOG(REAL(-band))/LOG(10.0D0)+1
ELSE IF (band == 0) THEN
  lo=0
ELSE
  lo=LOG(REAL(band))/LOG(10.0D0)
END IF

!      write(6,*) 'LO ',lo

READ(NAME,FMT='(255a1)')(ch(i),i=1,40)
l = length(ch,40)
!      write(6,*) 'L  ',l

!      write(6,*) 'CH ',(CH(I),I=1,25)

WRITE(s,FMT='(I10)') band
!      write(6,*) 'S  ',s

READ(s,FMT='(255A1)') (ch(i),i=l+1,l+10)
!      write(6,*) 'CH ',(CH(I),I=L+1,L+10)

WRITE(NEW,FMT='(255A1)') (ch(i),i=1,l),poi,(ch(i),i=l+10-lo,l+10)

RETURN
END SUBROUTINE sname
