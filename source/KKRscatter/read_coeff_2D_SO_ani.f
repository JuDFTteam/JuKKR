      SUBROUTINE READ_COEFF_2D_SO_ANI(COEFF_0,BZKP,LMMAXSO,NKPOID,NAEZ,
     +                            ENTG,NSPOH,ITHETA,IPHI)

      implicit none 

      INTEGER,INTENT(OUT)  ::  ENTG(NKPOID)
      INTEGER,INTENT(IN)  ::  NKPOID,LMMAXSO,NAEZ,NSPOH
      COMPLEX,INTENT(OUT) ::  COEFF_0(LMMAXSO,NAEZ,4,NKPOID)
      DOUBLE PRECISION,INTENT(OUT) ::  BZKP(3,NKPOID)

c  local variables
      INTEGER             ::  NK,LM1,STATUS_READ,K(NKPOID),
     +                        STATUS_OPEN,J,I1,ENT,EN,KSUM
      INTEGER             ::  ITHETA,IPHI
      CHARACTER*9         ::  FILENAME
      COEFF_0=(0d0,0d0)
      WRITE(6,*) "open coefficients",IPHI,ITHETA
c       WRITE(FILENAME,"((A5),(I3.3))") 'fort.',200+(ITHETA-1)*10+IPHI
       WRITE(FILENAME,"((A5),(I3.3))") 'coeff',(ITHETA-1)*10+IPHI
       WRITE(6,*) FILENAME

      IF (NSPOH == 1 ) THEN
c        open(unit=318, file="coefficients_0" , form="formatted",
c     +              action="read",iostat=status_open)
        open(unit=318, file="coefficients_0" , form="unformatted")
      ELSE
c        open(unit=318, file=FILENAME , form="formatted",
c     +              action="read",iostat=status_open)
        open(unit=318, file=FILENAME , form="unformatted")
      END IF

c      READ (318,"(I9)") KSUM
      READ (318) KSUM

          DO NK=1,KSUM
c            READ(318,"((I9),(3e17.9))",iostat=status_read)
c     +              K(NK),(BZKP(J,NK),J=1,3)
c            READ(318,"(I5)") ENTG(NK)
            READ(318) K(NK),(BZKP(J,NK),J=1,3)
            READ(318) ENTG(NK)
c             IF (status_read /=0 ) THEN
c               STOP "problem to read coefficients_0"
c             END IF
            DO ENT=1,ENTG(NK)
c              READ (318,"(I5)") EN
              READ (318) EN
              IF (ENT .NE. EN) stop "problem when reading coefficients"
              DO I1=1,NAEZ
                DO LM1=1,LMMAXSO
c                  READ(318,"(2(E20.12,2x))")  
c     +                 COEFF_0(LM1,I1,ENT,NK)
                  READ(318) COEFF_0(LM1,I1,ENT,NK)
                END DO
              END DO
            END DO
          ENDDO
         CLOSE(318)

      END SUBROUTINE READ_COEFF_2D_SO_ANI
