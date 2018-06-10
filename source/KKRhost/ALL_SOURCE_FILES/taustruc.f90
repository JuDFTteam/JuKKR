SUBROUTINE taustruct(drot,nsym,symunitary,nkm,nq,nqmax,nkmmax, iprint,irel)
!   ********************************************************************
!   *                                                                  *
!   *   find the structure of the site-diagonal TAU - matrices  TAUQ   *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! PARAMETER definitions
double COMPLEX C0
PARAMETER (C0=(0.0D0,0.0D0))

! Dummy arguments
INTEGER IPRINT,IREL,NKM,NKMMAX,NQ,NQMAX,NSYM
double COMPLEX DROT(NKMMAX,NKMMAX,48)
LOGICAL SYMUNITARY(48)

! Local variables
double precision ABST,X
DOUBLE PRECISION DBLE
INTEGER I,I0,IMWEIGHT,IQ,ISYM,IW,IWR,J,K,L,LIN,NELMT,NKMQ(NQMAX), &
        NKMTOP,NLIN,NON0(NQMAX)
double COMPLEX ST(NKMMAX,NKMMAX),TAUK(NKMMAX,NKMMAX,NQMAX)

DO iq = 1,nq
  nkmq(iq) = nkm
END DO

imweight = 0
nelmt = 0
nlin = 0
iw = 6

DO iq = 1,nq
  non0(iq) = 0
  nkmtop = nkmq(iq)
  
  IF ( iprint > 0 ) WRITE (1337,99004) iq
  DO i = 1,nkmtop
    DO j = 1,nkmtop
      st(i,j) = 0.0D0
      
      CALL cinit(nkmmax*nkmmax*nqmax,tauk)
      
      DO isym = 1,nsym
        i0 = iq
        
        IF ( symunitary(isym) ) THEN
          DO l = 1,nkmtop
            DO k = 1,nkmtop
              tauk(k,l,i0) = tauk(k,l,i0) + drot(i,k,isym)  &
                  *DCONJG(drot(j,l,isym))
            END DO
          END DO
        ELSE
          DO l = 1,nkmtop
            DO k = 1,nkmtop
              tauk(l,k,i0) = tauk(l,k,i0) + drot(i,k,isym)  &
                  *DCONJG(drot(j,l,isym))
            END DO
          END DO
        END IF
      END DO
      
      lin = 0
      iwr = 0
      DO k = 1,nkmq(iq)
        DO l = 1,nkmq(iq)
          abst = ABS(tauk(k,l,iq))
          st(i,j) = st(i,j) + abst
          IF ( abst > 1D-8 ) THEN
            IF ( DIMAG(tauk(k,l,iq)) > 1D-5 ) THEN
              IF ( iprint > 0 ) WRITE (1337,*) ' Im(Weight) > 1D-5 ',i,j,k,l
              imweight = 1
            END IF
            x = dreal(tauk(k,l,iq))/DBLE(nsym)
            
            IF ( iprint > 1 ) THEN
              IF ( iwr == 0 ) THEN
                iwr = 1
                WRITE (iw,99002) i,j,iq,x,k + (iq-1)*nkm, l + (iq-1)*nkm
              ELSE
                WRITE (iw,99003) x,k + (iq-1)*nkm, l + (iq-1)*nkm
              END IF
            END IF
            lin = lin + 1
          END IF
          
        END DO
      END DO
      
      IF ( lin > 0 ) THEN
        nlin = nlin + lin
        nelmt = nelmt + 1
        non0(iq) = non0(iq) + 1
      END IF
      
      IF ( ABS(st(i,j)) > 1D-5 ) st(i,j) = 2
      
    END DO
  END DO
  
  IF ( iprint > 1 ) CALL cmatstr('TAU-MAT',7,st,nkmtop,nkmmax,  &
      irel,irel,0,1D-8,6)
END DO

WRITE (1337,99005) nelmt,(non0(iq),iq=1,nq)
WRITE (1337,99006) nlin

IF ( imweight /= 0 ) WRITE (1337,99007)

!-----------------------------------------------------------------------

RETURN
99002 FORMAT ('     TAUQ(',i2,',',i2,',',i2,') =  ','   ',f10.4,' * <',  &
    i3,'|T(K)|',i3,'>')
99003 FORMAT (23X,' + ',f10.4,' * <',i3,'|T(K)|',i3,'>')
99004 FORMAT (//,  &
    ' ==========================================================='  &
    ,/,'   structure of  TAU-matrix   INT <i|t(k)|j>     IQ=',i3, /,  &
    ' ===========================================================' ,/)
99005 FORMAT (/,5X,'non-0 TAU-elements          ',i5,'   Q:',80I4)
99006 FORMAT (5X,'terms to sum up             ',i5,/)
99007 FORMAT (/,5X,50('#'),/,5X,'WARNING: complex TAU weights found',/,  &
    5X,'this may occur for rotated magnetic moments',/,5X,  &
    'relevant only for tetrahedron BZ-integration',/,5X, 50('#'),/)
END SUBROUTINE taustruct
