SUBROUTINE initabjij(iprint,naez,natyp,natomimp,nofgij,nqcalc,  &
    nsmax,nshell,iqcalc,atomimp,ish,jsh,  &
    ijtabcalc,ijtabsh,ijtabsym,nijcalc,kijsh, nijmax,nshell0,nsheld)
!   ********************************************************************
!   *  subroutine called by < TBXCCPLJIJ > to set up some auxiliary    *
!   *  arrays allowing the indexing of shells, sites, atomic types     *
!   ********************************************************************

IMPLICIT NONE

! Arguments
INTEGER IPRINT,NAEZ,NATOMIMP,NATYP,NIJMAX,NOFGIJ,NQCALC,NSHELD, &
        NSHELL0,NSMAX
INTEGER ATOMIMP(*),IJTABCALC(*),IJTABSH(*),IJTABSYM(*),IQCALC(*), &
        ISH(NSHELD,*),JSH(NSHELD,*),KIJSH(NIJMAX,NSHELL0), &
        NIJCALC(NSHELL0),NSHELL(0:NSHELD)

! Locals
INTEGER I1,IA,IDONE(NAEZ),IQTOJQ(NIJMAX),J1,JA,LM1,LM2,NS
INTEGER NIDONE

! ======================================================================
DO ns = nsmax + 1,nshell(0)
  DO i1 = 1,nijmax
    iqtojq(i1) = 0
  END DO
! ----------------------------------------------------------------------
  DO i1 = 1,nshell(ns)
    ia = atomimp(ish(ns,i1))
    ja = 0
    DO j1 = 1,nijcalc(ns)
      IF ( ia == iqtojq(j1) ) THEN
        ja = 1
        GO TO 20
      END IF
    END DO
    20         CONTINUE
    IF ( ja == 0 ) THEN
      nijcalc(ns) = nijcalc(ns) + 1
      IF ( nijcalc(ns) > nijmax ) THEN
        WRITE (6,99001) 'local','NIJMAX',nijcalc(ns)
        STOP '       in < TBXCCPLJIJ > '
      END IF
      iqtojq(nijcalc(ns)) = ia
      kijsh(nijcalc(ns),ns) = i1
    END IF
  END DO
END DO
! ======================================================================
IF ( iprint <= 0 ) RETURN
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
nidone = 0
DO ia = 1,naez
  idone(ia) = 0
  DO i1 = 1,nqcalc
    IF ( iqcalc(i1) == ia ) THEN
      nidone = nidone + 1
      idone(nidone) = ia
    END IF
  END DO
END DO

lm2 = MIN(25,natomimp)
WRITE (1337,99002) naez,natyp,natomimp,nofgij,nshell(0),lm2
DO i1 = 1,nidone
  DO ia = 1,natomimp
    IF (atomimp(ia) == idone(i1)) THEN
      lm1 = (ia-1)*natomimp
      WRITE (1337,99003) ia,(ijtabcalc(lm1+ja),ja=1,lm2)
      GO TO 100
    END IF
  END DO
  100     CONTINUE
END DO
WRITE (1337,99004) lm2
DO i1 = 1,nidone
  DO ia = 1,natomimp
    IF (atomimp(ia) == idone(i1)) THEN
      lm1 = (ia-1)*natomimp
      WRITE (1337,99003) ia,(ijtabsh(lm1+ja),ja=1,lm2)
      GO TO 110
    END IF
  END DO
  110     CONTINUE
END DO
WRITE (1337,99005) lm2
DO i1 = 1,nidone
  DO ia = 1,natomimp
    IF (atomimp(ia) == idone(i1)) THEN
      lm1 = (ia-1)*natomimp
      WRITE (1337,99003) ia,(ijtabsym(lm1+ja),ja=1,lm2)
      GO TO 120
    END IF
  END DO
  120     CONTINUE
END DO
lm2 = 0
DO ns = nsmax + 1,nshell(0)
  lm2 = MAX(lm2,nijcalc(ns))
END DO
lm2 = MIN(5,lm2)
WRITE (1337,99006)
DO ns = nsmax + 1,nshell(0)
  WRITE (1337,99007) ns,(ish(ns,kijsh(i1,ns)),  &
      jsh(ns,kijsh(i1,ns)),i1=1,MIN(nijcalc(ns),lm2))
  WRITE (1337,99008) (atomimp(ish(ns,kijsh(i1,ns))),  &
      atomimp(jsh(ns,kijsh(i1,ns))), ijtabsym((ish(ns,kijsh(i1,ns))-1)  &
      *natomimp+jsh(ns,kijsh(i1,ns))),i1=1, MIN(nijcalc(ns),lm2))
END DO
!ccc      WRITE (6,99009) MIN(NATYP,25)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

99001 FORMAT (6X,'Dimension ERROR: please increase the ',a,' parameter',  &
    /,6X,a,' to a value >=',i5,/)
99002 FORMAT (8X,60('-'),/8X,'Data used for J_ij calculation:',//,10X,  &
    'Number of sites/types i        (NAEZ/NATYP) :',2(1X,i3),  &
    /,10X,'Number of atoms in the cluster   (NATOMIMP) :',1X,  &
    i3,/,10X,'Number of ij pairs                 (NOFGIJ) :', 1X,i3,/,10X,  &
    'Number of representative pairs     (NSHELL) :',1X,i3,//,  &
    10X,'ij-pairs calculation table ( 1 = calculated )',/,10X,  &
    'IA   JA = 1 ..',i3)
99003 FORMAT (10X,i3,3X,25(i3))
99004 FORMAT (/,10X,'ij-shells table ',/,10X,'IA   JA = 1 ..',i3)
99005 FORMAT (/,10X,'ij-symmetries table ',/,10X,'IA   JA = 1 ..',i3)
99006 FORMAT (/,10X,'effectively calculated pairs/shells',/,10X,  &
    'SHELL   (IAT,JAT) ',/,10X, 'SHELL   (IQ,JQ - ISYM) ')
99007 FORMAT (10X,i4,3X,5(i3,',',i3,5X))
99008 FORMAT (10X,4X,3X,5(i3,',',i3,' - ',i2))
99009 FORMAT (/,10X,'effectively calculated type-type pairs (shells)',  &
    /,10X,'IT   JT = 1 ..',i3)
END SUBROUTINE initabjij
