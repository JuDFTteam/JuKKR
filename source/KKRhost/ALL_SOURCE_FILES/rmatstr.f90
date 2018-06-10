SUBROUTINE rmatstr(str,lstr,a,n,m,mlin,mcol,tolp,nfil)
!   ********************************************************************
!   *                                                                  *
!   *   writes structure of REAL      NxN   matrix   A                 *
!   *                                                                  *
!   *   M           is the actual array - size used for   A            *
!   *   MLIN/COL    MODE for line and comlun indexing                  *
!   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
!   *   TOL         tolerance for difference                           *
!   *                                                                  *
!   *                                                         03/01/96 *
!   ********************************************************************
IMPLICIT NONE

INTEGER, PARAMETER :: NDIFMAX =  250
INTEGER, intent(in) :: lstr, m, n, mlin, mcol, nfil
double precision, intent(in) :: tolp

CHARACTER (len=LSTR) :: STR
CHARACTER (len=150) :: FMT1,FMT2,FMT3,FMT4
CHARACTER (len=1) :: CTAB(0:NDIFMAX), VZ(-1:+1)
INTEGER    IW(M), ILSEP(20)
INTEGER :: ICAUC , ICZUC, NATOZ, ICALC, ICILC, K, NK, NM, NM2, &
           NM1, NM3, IC0, N3, N2, N1, LF, L3, MM, NSL, IL, I, &
           NNON0, ND, NALF, J, ID, ICND, ISL
double precision :: tol
DOUBLE PRECISION     A(M,M), CNUM,CA,CB,DTAB(0:NDIFMAX)  
LOGICAL  CSMALL, CSAME

SAVE VZ
DATA VZ / '-', ' ', ' ' /

csmall(cnum) =   ABS(cnum*tol) < 1.0D0

csame(ca,cb) =   csmall(1.0D0 - ca/cb)

tol = 1.0D0 / tolp

icauc = ICHAR('A')
iczuc = ICHAR('Z')
natoz = iczuc - icauc + 1
icalc = ICHAR('a')
icilc = ICHAR('i')

!---------------------------------------------------------------- header
ic0 = ICHAR('0')
n3=n/100
n2=n/10 - n3*10
n1=n    - n2*10 - n3*100

fmt1='(8X,I3,''|'','
fmt2='( 9X,''--|'','
fmt3='( 9X,'' #|'','
fmt4='( 9X,''  |'','

lf  = 11
l3  = 11
IF( mcol == 0 ) THEN
  fmt1=fmt1(1:lf)//CHAR(ic0+n3)//CHAR(ic0+n2)  &
      //CHAR(ic0+n1)//'( 2A1),''|'',I3)'
  fmt2=fmt2(1:lf)//CHAR(ic0+n3)//CHAR(ic0+n2)  &
      //CHAR(ic0+n1)//'(''--''),''|'',I3)'
  fmt3=fmt3(1:lf)//'60(2X,I2))'
  fmt4=fmt4(1:lf)//'60(I2,2X))'
  lf = 21
ELSE
  IF( mcol == 1 ) THEN
    nk = nint(SQRT( DBLE(n) ))
  ELSE IF( mcol == 2 ) THEN
    nk = nint(SQRT( DBLE(n/2) ))
  ELSE IF( mcol == 3 ) THEN
    nk = 2*nint(SQRT( DBLE(n/2) )) - 1
  END IF
  DO k=1,nk
    IF( mcol <= 2 ) THEN
      nm = 2*k-1
    ELSE
      nm = 2*((k+1)/2)
    END IF
    nm2=nm/10
    nm1=nm - nm2*10
    nm3=nm/2
    fmt1=fmt1(1:lf)//CHAR(ic0+nm2) //CHAR(ic0+nm1)//'( 2A1),''|'','
    fmt2=fmt2(1:lf)//CHAR(ic0+nm2) //CHAR(ic0+nm1)//'(''--''),''|'','
    
    IF( mcol <= 2 ) THEN
      DO mm=1,nm
        IF( MOD(mm,2) == MOD(k,2) ) THEN
          fmt3=fmt3(1:l3)//'2X,'
          fmt4=fmt4(1:l3)//'I2,'
        ELSE
          fmt3=fmt3(1:l3)//'I2,'
          fmt4=fmt4(1:l3)//'2X,'
        END IF
        l3=l3+3
      END DO
      fmt3=fmt3(1:l3)//'''|'','
      fmt4=fmt4(1:l3)//'''|'','
      l3=l3+4
    ELSE
      fmt3=fmt3(1:lf)//CHAR(ic0+nm3)//'(2X,I2),''|'','
      fmt4=fmt4(1:lf)//CHAR(ic0+nm3)//'(I2,2X),''|'','
      l3=l3+13
    END IF
    lf  = lf + 13
  END DO
  IF( mcol == 2 ) THEN
    fmt1=fmt1(1:lf)//fmt1(12:lf)
    fmt2=fmt2(1:lf)//fmt2(12:lf)
    fmt3=fmt3(1:l3)//fmt4(12:l3)
    fmt4=fmt4(1:l3)//fmt3(12:l3)
    lf  = 2*lf - 11
  END IF
  fmt1=fmt1(1:lf)//'I3)'
  fmt2=fmt2(1:lf)//'I3)'
  fmt3=fmt3(1:l3)//'I3)'
  fmt4=fmt4(1:l3)//'I3)'
END IF
IF( mlin == 0 ) THEN
  nsl = 1
  ilsep(1) = n
ELSE IF( mlin == 1 ) THEN
  nsl = nint(SQRT( DBLE(n) ))
  DO il=1,nsl
    ilsep(il) = il**2
  END DO
ELSE IF( mlin == 2 ) THEN
  nsl = nint(SQRT( DBLE(n/2) ))
  DO il=1,nsl
    ilsep(il) = il**2
  END DO
  DO il=1,nsl
    ilsep(nsl+il) = ilsep(nsl) + il**2
  END DO
  nsl = 2*nsl
ELSE IF( mlin == 3 ) THEN
  nsl = 2*nint(SQRT( DBLE(n/2) )) - 1
  ilsep(1) = 2
  DO k=2,nsl
    ilsep(k) = ilsep(k-1) + 2*((k+1)/2)
  END DO
END IF


WRITE(nfil,9000) str(1:lstr)
WRITE(nfil,fmt3) (i,i=2,n,2)
WRITE(nfil,fmt4) (i,i=1,n,2)
WRITE(nfil,FMT=fmt2)
!------------------------------------------------------------ header end
nnon0   = 0
nd      = 0
nalf    = 0
ctab(0) = ' '
dtab(0) = 9999D0

DO  i=1,n
  DO  j=1,n
    IF( .NOT. csmall( a(i,j) ) ) THEN
      nnon0 = nnon0 + 1
      DO  id=1,nd
        IF( csame(a(i,j),+dtab(id)) ) THEN
          iw(j) = + id
          GO TO 40
        END IF
        IF( csame(a(i,j),-dtab(id)) ) THEN
          iw(j) = - id
          GO TO 40
        END IF
      END DO
!----------------------------------------------------------- new element
      nd = nd + 1
      IF( nd > ndifmax ) THEN
        WRITE(nfil,'('' TROUBLE IN <RMATSTR> !!!!'',/,  &
            '' ND > ARRAY SIZE NDIFMAX='',I3)') ndifmax
        STOP
      END IF
      iw(j)    = nd
      dtab(nd) = a(i,j)
      IF( ABS(dtab(nd)-1.0D0)*tol < 1.0D0 ) THEN
        ctab(nd) = '1'
      ELSE IF( ABS(dtab(nd)+1.0D0)*tol < 1.0D0 ) THEN
        dtab(nd) = +1.0D0
        ctab(nd) = '1'
        iw(j) = - nd
      ELSE
        nalf = nalf + 1
        IF( nalf <= natoz ) THEN
          ctab(nd) = CHAR( icauc + nalf - 1 )
        ELSE
          icnd = icalc + nalf-natoz - 1
          IF( icnd < icilc ) THEN
            ctab(nd) = CHAR(icnd)
          ELSE
            ctab(nd) = CHAR(icnd+1)
          END IF
        END IF
      END IF
      40          CONTINUE
    ELSE
      iw(j) = 0
    END IF
  END DO
!------------------------------------------------------------ write line
  WRITE(nfil,FMT=fmt1) i,(vz(ISIGN(1,iw(j))),ctab(ABS(iw(j))),j=1,n),i
  
  
  DO isl=1,nsl
    IF( i == ilsep(isl) ) WRITE(nfil,FMT=fmt2)
  END DO
END DO

!------------------------------------------------------------------ foot

WRITE(nfil,fmt4) (i,i=1,n,2)
WRITE(nfil,fmt3) (i,i=2,n,2)

WRITE(nfil,9030) (id,ctab(id),dtab(id),id=1,nd)
WRITE(nfil,9040)  nnon0,n*n-nnon0

9000 FORMAT(/,8X,a,/)
9030 FORMAT(/,8X,'symbols used:',/,(8X,i3,3X,a1,2X, f20.12) )
9040 FORMAT(/,8X,'elements <> 0:',i4,/, 8X,'elements  = 0:',i4)
RETURN
END SUBROUTINE rmatstr
