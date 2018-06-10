SUBROUTINE cmatstr(str,lstr,a,n,m,mlin,mcol,ijq,tolp,k_fmt_fil)
!   ********************************************************************
!   *                                                                  *
!   *   writes structure of COMPLEX   NxN   matrix   A                 *
!   *                                                                  *
!   *   M           is the actual array - size used for   A            *
!   *   MLIN/COL    MODE for line and column indexing                  *
!   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
!   *   TOL         tolerance for difference                           *
!   *   IJQ         if IJQ > 1000    pick  IQ-JQ-block matrix          *
!   *               assuming  IJQ = IQ*1000 + JQ                       *
!   *               else: no IQ-JQ-indexing                            *
!   *   K_FMT_FIL   output channel                                     *
!   *               a negative sign suppresses table at the end        *
!   *                                                                  *
!   *   any changes should be done in RMATSTR as well !!!!!!!!!!!!!!!  *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! PARAMETER definitions
DOUBLE COMPLEX CI
PARAMETER (CI=(0.0D0,1.0D0))

! Dummy arguments
INTEGER IJQ,K_FMT_FIL,LSTR,M,MCOL,MLIN,N
CHARACTER (len=LSTR) :: STR
DOUBLE PRECISION TOLP
DOUBLE COMPLEX A(M,M)

! Local variables
DOUBLE COMPLEX B(N,N),CA,CB,ARG,DTAB(0:N*N)
CHARACTER CHAR
LOGICAL SAME,SMALL
CHARACTER (len=1) :: CTAB(0:N*N),VZ(-1:+1)
DOUBLE PRECISION DBLE
CHARACTER (len=150) ::FMT1,FMT2,FMT3,FMT4
INTEGER I,I1,IC0,ID,IL,ILSEP(20),IPT(218),IQ,ISL,IW(M),J, &
        J0,JP,JQ,K,L3,LF,MM,N1,N2,N3,NC,ND,NFIL,NK,NM,NM1,NM2,NM3, &
        NNON0,NSL
INTEGER ICHAR,ISIGN,NINT
DOUBLE PRECISION TOL

DATA VZ/'-',' ',' '/
small(arg) = ABS(arg*tol) < 1.0D0

same(ca,cb) = small(1.0D0-ca/cb)

nfil = ABS(k_fmt_fil)

tol = 1.0D0/tolp

!----------------------------------------------- set block indices IQ JQ

IF ( ijq > 1000 ) THEN
  iq = ijq/1000
  jq = ijq - iq*1000
  IF ( iq*n > m .OR. iq*n > m ) THEN
    WRITE (1337,99002) ijq,iq,jq,iq*n,jq*n,n,m
    RETURN
  END IF
ELSE
  iq = 1
  jq = 1
END IF

!----------------------------------------------------- copy matrix block

j0 = n*(jq-1)
DO j = 1,n
  i1 = n*(iq-1)+1
  jp = j0 + j
  CALL zcopy(n,a(i1,jp),1,b(1,j),1)
END DO

!------------------------------------------------ set up character table

nc = 0
DO i = 1,26
  nc = nc + 1
  ipt(nc) = 62 + i
END DO
DO i = 1,8
  nc = nc + 1
  ipt(nc) = 96 + i
END DO
DO i = 10,26
  nc = nc + 1
  ipt(nc) = 96 + i
END DO
DO i = 191,218
  nc = nc + 1
  ipt(nc) = i
END DO
DO i = 35,38
  nc = nc + 1
  ipt(nc) = i
END DO
DO i = 40,42
  nc = nc + 1
  ipt(nc) = i
END DO
DO i = 91,93
  nc = nc + 1
  ipt(nc) = i
END DO

!---------------------------------------------------------------- header
ic0 = ICHAR('0')
n3 = n/100
n2 = n/10 - n3*10
n1 = n - n2*10 - n3*100

IF ( n <= 18 ) THEN
  fmt1 = '(8X,I3,''|'','
  fmt2 = '( 9X,''--|'','
  fmt3 = '( 9X,'' #|'','
  fmt4 = '( 9X,''  |'','
ELSE
  fmt1 = '(   I4,''|'','
  fmt2 = '( 2X,''--|'','
  fmt3 = '( 2X,'' #|'','
  fmt4 = '( 2X,''  |'','
END IF

lf = 11
l3 = 11
IF ( mcol == 0 ) THEN
  fmt1 = fmt1(1:lf)//CHAR(ic0+n3)//CHAR(ic0+n2)//CHAR(ic0+n1)  &
      //'( 2A1),''|'',I3)'
  fmt2 = fmt2(1:lf)//CHAR(ic0+n3)//CHAR(ic0+n2)//CHAR(ic0+n1)  &
      //'(''--''),''|'',I3)'
  fmt3 = fmt3(1:lf)//'60(2X,I2))'
  fmt4 = fmt4(1:lf)//'60(I2,2X))'
  lf = 21
ELSE
  IF ( mcol == 1 ) THEN
    nk = nint(SQRT(DBLE(n)))
  ELSE IF ( mcol == 2 ) THEN
    nk = nint(SQRT(DBLE(n/2)))
  ELSE IF ( mcol == 3 ) THEN
    nk = 2*nint(SQRT(DBLE(n/2))) - 1
  END IF
  DO k = 1,nk
    IF ( mcol <= 2 ) THEN
      nm = 2*k - 1
    ELSE
      nm = 2*((k+1)/2)
    END IF
    nm2 = nm/10
    nm1 = nm - nm2*10
    nm3 = nm/2
    fmt1 = fmt1(1:lf)//CHAR(ic0+nm2)//CHAR(ic0+nm1) //'( 2A1),''|'','
    fmt2 = fmt2(1:lf)//CHAR(ic0+nm2)//CHAR(ic0+nm1) //'(''--''),''|'','
    
    IF ( mcol <= 2 ) THEN
      DO mm = 1,nm
        IF ( MOD(mm,2) == MOD(k,2) ) THEN
          fmt3 = fmt3(1:l3)//'2X,'
          fmt4 = fmt4(1:l3)//'I2,'
        ELSE
          fmt3 = fmt3(1:l3)//'I2,'
          fmt4 = fmt4(1:l3)//'2X,'
        END IF
        l3 = l3 + 3
      END DO
      fmt3 = fmt3(1:l3)//'''|'','
      fmt4 = fmt4(1:l3)//'''|'','
      l3 = l3 + 4
    ELSE
      fmt3 = fmt3(1:lf)//CHAR(ic0+nm3)//'(2X,I2),''|'','
      fmt4 = fmt4(1:lf)//CHAR(ic0+nm3)//'(I2,2X),''|'','
      l3 = l3 + 13
    END IF
    lf = lf + 13
  END DO
  IF ( mcol == 2 ) THEN
    fmt1 = fmt1(1:lf)//fmt1(12:lf)
    fmt2 = fmt2(1:lf)//fmt2(12:lf)
    fmt3 = fmt3(1:l3)//fmt3(12:l3)
    fmt4 = fmt4(1:l3)//fmt4(12:l3)
    lf = 2*lf - 11
    l3 = 2*l3 - 11
  END IF
  fmt1 = fmt1(1:lf)//'I3)'
  fmt2 = fmt2(1:lf)//'I3)'
  fmt3 = fmt3(1:l3)//'I3)'
  fmt4 = fmt4(1:l3)//'I3)'
END IF
IF ( mlin == 0 ) THEN
  nsl = 1
  ilsep(1) = n
ELSE IF ( mlin == 1 ) THEN
  nsl = nint(SQRT(DBLE(n)))
  DO il = 1,nsl
    ilsep(il) = il**2
  END DO
ELSE IF ( mlin == 2 ) THEN
  nsl = nint(SQRT(DBLE(n/2)))
  DO il = 1,nsl
    ilsep(il) = il**2
  END DO
  DO il = 1,nsl
    ilsep(nsl+il) = ilsep(nsl) + il**2
  END DO
  nsl = 2*nsl
ELSE IF ( mlin == 3 ) THEN
  nsl = 2*nint(SQRT(DBLE(n/2))) - 1
  ilsep(1) = 2
  DO k = 2,nsl
    ilsep(k) = ilsep(k-1) + 2*((k+1)/2)
  END DO
END IF


WRITE (nfil,99001) str(1:lstr)
IF ( ijq > 1000 ) WRITE (nfil,99003) iq,jq
WRITE (nfil,fmt3) (i,i=2,n,2)
WRITE (nfil,fmt4) (i,i=1,n,2)
WRITE (nfil,FMT=fmt2)
!------------------------------------------------------------ header end
nnon0 = 0
nd = 0
ctab(0) = ' '
dtab(0) = 9999D0

DO i = 1,n
  DO j = 1,n
    IF ( .NOT.small(b(i,j)) ) THEN
      nnon0 = nnon0 + 1
      DO id = 1,nd
        IF ( same(b(i,j),+dtab(id)) ) THEN
          iw(j) = +id
          GO TO 50
        END IF
        IF ( same(b(i,j),-dtab(id)) ) THEN
          iw(j) = -id
          GO TO 50
        END IF
      END DO
!----------------------------------------------------------- new element
      nd = nd + 1
      iw(j) = nd
      dtab(nd) = b(i,j)
      IF ( ABS(dtab(nd)-1.0D0)*tol < 1.0D0 ) THEN
        ctab(nd) = '1'
      ELSE IF ( ABS(dtab(nd)+1.0D0)*tol < 1.0D0 ) THEN
        dtab(nd) = +1.0D0
        ctab(nd) = '1'
        iw(j) = -nd
      ELSE IF ( ABS(dtab(nd)-ci)*tol < 1.0D0 ) THEN
        ctab(nd) = 'i'
      ELSE IF ( ABS(dtab(nd)+ci)*tol < 1.0D0 ) THEN
        dtab(nd) = +ci
        ctab(nd) = 'i'
        iw(j) = -nd
      ELSE
        ctab(nd) = CHAR(ipt(1+MOD((nd+1),nc)))
      END IF
    ELSE
      iw(j) = 0
    END IF
  50      END DO
!------------------------------------------------------------ write line
  WRITE (nfil,FMT=fmt1) i, (vz(ISIGN(1,iw(j))),ctab(ABS(iw(j))),j=1,  &
      n),i
  
  DO isl = 1,nsl
    IF ( i == ilsep(isl) ) WRITE (nfil,FMT=fmt2)
  END DO
END DO

!------------------------------------------------------------------ foot

WRITE (nfil,fmt4) (i,i=1,n,2)
WRITE (nfil,fmt3) (i,i=2,n,2)

IF ( k_fmt_fil > 0 ) THEN
  WRITE (nfil,99004) (id,ctab(id),dtab(id),id=1,nd)
  WRITE (nfil,99005) nnon0,tolp,n*n - nnon0,tolp
ELSE
  WRITE (nfil,*) ' '
END IF

99001 FORMAT (/,8X,a,/)
99002 FORMAT (/,1X,79('*'),/,10X,'inconsistent call of <CMATSTR>',/,10X,  &
    'argument IJQ =',i8,'  implies IQ=',i3,'   JQ=',i3,/,10X,  &
    'IQ*N=',i6,' > M   or   JQ*N=',i6,' > M   for N =',i4,  &
    ' M=',i4,/,1X,79('*'),/)
99003 FORMAT (8X,'IQ-JQ-block  for  IQ = ',i3,'   JQ = ',i3,/)
99004 FORMAT (/,8X,'symbols used:',/,(8X,i3,3X,a1,2X,2F20.12))
99005 FORMAT (/,8X,i5,' elements   >',1PE9.1,/,  &
    8X,i5,' elements   <',1PE9.1,/)
END SUBROUTINE cmatstr
