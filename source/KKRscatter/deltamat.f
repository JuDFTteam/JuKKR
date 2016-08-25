C ************************************************************************
      SUBROUTINE DELTAMAT(TMATLL,DELTALL,E)
C ************************************************************************
      include 'inc.p'
      INTEGER LMAXSQ
      PARAMETER(LMAXSQ = (LMAXD+1)**2)
      DOUBLE COMPLEX E
      DOUBLE COMPLEX TMATLL(LMAXSQ,*),
     +               DELTALL(LMAXSQ,*)
C     .. local scalars
      INTEGER I,L
C     .. External Functions ..
      DOUBLE COMPLEX EXPIDL
      LOGICAL TEST
      EXTERNAL CINIT,EXPIDL,TEST
      DOUBLE PRECISION PI
      PARAMETER (PI= 3.14159265358979312d0)
      INTEGER LF(144)
      DATA LF/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11/
c ------------------------------------------------------------------------
      CALL CINIT(LMAXSQ**2,DELTALL)
      DO 10 I=1,LMAXSQ
        DELTALL(I,I) = EXPIDL(TMATLL(I,I),E,LF(I))
 10   CONTINUE
      IF (TEST('deltall ')) THEN
        write(6,*) 'DELTALL :'
        DO L=1,LMAXD+1
          I=L**2
          WRITE(6,FMT='(1p,2(d14.3,d11.3),0p,f14.6)') 
     +         DELTALL(I,I),
     +         -(1.d0,0.d0)/sqrt(e)*dimag(DELTALL(I,I))*DELTALL(I,I),
     +         DIMAG(1.D0/PI*ZLOG(DELTALL(I,I)))
        END DO
      END IF  
      RETURN
      END
