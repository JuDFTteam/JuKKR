MODULE MOD_ETOTB1
  CONTAINS

  SUBROUTINE ETOTB1(EFERMI,LMAXATOM,ENERGYPARTS,CORESTATE,NSPIN,NATOM,LPOTD)
    !     calculate the total energy of the cluster .
    !     gather all energy-parts which are calculated in different
    !     subroutines .
    !     since the program uses group theory only shell-indices
    !     are used instead of atom-indices .
    !
    !                               b.drittler   may 1987
    !
    !     modified for supercells with nshell(i) atoms of type i in the
    !     unit cell
    !                               p.zahn       oct. 95
    !
    !     adopted for more atoms per site (CPA) v.popescu feb. 02
    !-----------------------------------------------------------------------
    USE TYPE_ENERGYPARTS
    USE TYPE_CORESTATE
    IMPLICIT NONE

    ! PARAMETER definitions
    INTEGER :: LMAXATOM(NATOM)
    TYPE(ENERGYPARTS_TYPE) ENERGYPARTS
    TYPE(CORESTATE_TYPE) CORESTATE(NATOM)
    INTEGER :: LPOTD
    ! Dummy arguments
    DOUBLE PRECISION EFERMI
    INTEGER KPRE,NATOM,NSPIN
    ! Local variables
    DOUBLE PRECISION BANDESUM,BANDET,ECOUS,EDC,EFCTOR,ET,ETOT,EXCS
    DOUBLE PRECISION ETOTLDAU
    INTEGER IATOM,IPOT,IS,ISPIN,L
    CHARACTER (len=4) TEXTL(0:6)
    CHARACTER (len=5) TEXTNS
    CHARACTER (len=13) TEXTS(3)
    !     .. Data statements ..
    DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
    DATA TEXTS/' spin down   ',' spin  up    ',' paramagnetic'/
    DATA TEXTNS/' ns ='/
    ! ------------------------------------------------------------------------
    EFCTOR = 1.0D0/13.6058D0

    ETOT = 0.0D0
    BANDESUM = 0.0D0
    ETOTLDAU = 0.0D0

    KPRE=1
    IF ( KPRE.EQ.1 ) WRITE (1337,FMT=99001)

    !---> loop over host atoms
    DO IATOM = 1,NATOM

       IF ( KPRE.EQ.1 ) WRITE (1337,FMT=99002) IATOM

       EDC = 0.0D0
       ET = 0.0D0
       BANDET = 0.0D0
       ECOUS = 0.0D0
       EXCS = 0.0D0

       IS = 0
       IF ( NSPIN.EQ.1 ) IS = IS + 2
       DO ISPIN = 1,NSPIN
          IS = IS + 1
          IPOT = (IATOM-1)*NSPIN + ISPIN

          IF ( KPRE.EQ.1 ) THEN
           WRITE (1337,FMT=99003) TEXTS(IS)
           WRITE (1337,FMT=99004) (TEXTL(L),ENERGYPARTS%ESPC(L,ISPIN,IATOM),L=0,corestate(IATOM)%LCOREMAX)
           WRITE (1337,FMT=99005) (TEXTL(L),ENERGYPARTS%ESPV(L,ISPIN,IATOM),L=0,LMAXATOM(IATOM))
           WRITE (1337,FMT=99006) TEXTNS,ENERGYPARTS%ESPV(LMAXATOM(IATOM)+1,ISPIN,IATOM)
          END IF

          DO L = 0,corestate(IATOM)%LCOREMAX
             ET = ET + ENERGYPARTS%ESPC(L,ISPIN,IATOM)
          END DO

          DO L = 0,LMAXATOM(IATOM)
             BANDET = BANDET + ENERGYPARTS%ESPV(L,ISPIN,IATOM)
             ET = ET + ENERGYPARTS%ESPV(L,ISPIN,IATOM)
          END DO
          BANDET = BANDET + ENERGYPARTS%ESPV(LMAXATOM(IATOM)+1,ISPIN,IATOM)
          ET = ET + ENERGYPARTS%ESPV(LMAXATOM(IATOM)+1,ISPIN,IATOM)
       END DO

    ! -> LDA+U
    !          ET = ET + EU(IATOM)
    !          BANDET = BANDET + EU(IATOM)
    !          IF ( KPRE.EQ.1 .AND. IDOLDAU.EQ.1 .AND. LOPT(IATOM).GE.0 )
    !        WRITE(6,99019) EU(IATOM)

    ! --->  sum up Coulomb and Ex.-Corel. contribution
       DO L = 0,2*LMAXATOM(IATOM)
          ECOUS = ECOUS + ENERGYPARTS%ECOU(L,IATOM)
          EXCS = EXCS + ENERGYPARTS%EXC(L,IATOM)
       END DO

       IF ( KPRE.EQ.1 ) THEN
          WRITE (1337,FMT=99007) ET
          WRITE (1337,FMT=99008) BANDET
          WRITE (1337,FMT=99009) (L, ENERGYPARTS%ECOU(L,IATOM),L=0,2*LMAXATOM(IATOM))
          WRITE (1337,FMT=99010)
          WRITE (1337,FMT=99018) ECOUS
          WRITE (1337,FMT=99011) (L,ENERGYPARTS%EXC(L,IATOM), L=0,2*LMAXATOM(IATOM))
          WRITE (1337,FMT=99010)
          WRITE (1337,FMT=99017) EXCS
          WRITE (1337,FMT=99015) ENERGYPARTS%EPOTIN(IATOM)
       END IF

       WRITE (22349378,*) IATOM,ET/EFCTOR
       WRITE (22349379,*) IATOM,BANDET/EFCTOR


    !          IF ( .NOT.(TEST('NoMadel ')) ) THEN
          ET = ET + ECOUS + EXCS
          EDC = EDC + ECOUS + EXCS


          ET = ET + ENERGYPARTS%EPOTIN(IATOM) !- EDCLDAU(IATOM)
          EDC = EDC + ENERGYPARTS%EPOTIN(IATOM) !- EDCLDAU(IATOM)

          IF ( KPRE.EQ.1 ) THEN
    !                IF ( IDOLDAU.EQ.1 .AND. LOPT(IATOM).GE.0 )
    !              WRITE(6,99020) -EDCLDAU(IATOM)
             WRITE (1337,FMT=99016) EDC
          END IF

    !          END IF

       IF ( NATOM.GT.1 ) THEN
          WRITE (1337,FMT=99012) IATOM,ET
    !             IF ( KPRE.EQ.1 .AND. IDOLDAU.EQ.1 .AND. LOPT(IATOM).GE.0 )
    !           WRITE(6,99021) EU(IATOM) - EDCLDAU(IATOM)
    !             WRITE (6,FMT=99022)
       END IF

       ETOT = ETOT + ET !*DBLE(IATOM) !*CONC(IATOM)
       BANDESUM = BANDESUM + BANDET !*DBLE(IATOM) !*CONC(IATOM)

    END DO                        ! IATOM = 1,NATOM

    WRITE (1337,FMT=99013) BANDESUM
    WRITE (1337,FMT=99014) ETOT,ETOT/EFCTOR
    !       WRITE (*,FMT=99013) BANDESUM
    WRITE (*,FMT=99014) ETOT,ETOT/EFCTOR
    WRITE (22349375,*) ETOT/EFCTOR
    WRITE (22349376,*) BANDESUM/EFCTOR

    RETURN

    99001 FORMAT (32('='),' TOTAL ENERGIES ',31('='),/)
    99002 FORMAT (3X,'Total energies atom ',i3,/,3x,23('-'))
    99003 FORMAT (5X,'single particle energies ',a13)
    99004 FORMAT (7X,'  core   contribution : ',2(a4,f15.8),/,        (31X,2(a4,f15.8)))
    99005 FORMAT (7X,'valence  contribution : ',2(a4,f15.8),/,        (31X,2(a4,f15.8)))
    99006 FORMAT (7x,'                        ',a4,f15.8)
    99007 FORMAT (5X,68('-'),/,5X,        'total contribution of the single particle energies :',1X,        f15.8)
    99008 FORMAT (5X,'                              band energy per atom :',        1X,f15.10,/)
    99009 FORMAT (5X,'coulomb  contribution : ',2(i3,1X,f15.8),/,        (29X,2(i3,1X,f15.8)))
    99010 FORMAT (5X,68('-'))

    99011 FORMAT (5X,'ex.-cor. contribution : ',2(i3,1X,f15.8),/,        (29X,2(i3,1X,f15.8)))
    99012 FORMAT (/,3X,'Total contribution of atom',i3,' =',f15.8)
    99013 FORMAT (5X,'                              sum of band energies :',        1X,F15.10,/,3X,70('-'))
    99014 FORMAT (/,3X,70('+'),/,15x,'TOTAL ENERGY in ryd. : ',f17.8,/15x,        '                 eV  : ',F17.8,/,3x,70('+'))
    99015 FORMAT (5X,'eff. pot. contribution     : ',f15.8)
    99016 FORMAT (5X,'total double counting contribution                 :',        1X,f15.8)
    99017 FORMAT (5X,'tot. ex.-cor. contribution : ',f15.8,/)
    99018 FORMAT (5X,'tot. coulomb contribution : ',f15.8,/)
    99019 FORMAT (/,5X,        'LDA+U correction to the single particle energy     :',        F16.8)
    99020 FORMAT (/,5X,        'LDA+U double counting contribution                 :',        F16.8)
    99021 FORMAT (3X,'   including LDA+U correction :',F15.8)
    99022 FORMAT (3X,70('-'))
  END SUBROUTINE
END MODULE MOD_ETOTB1
