       SUBROUTINE TESTDIM(NSPIN,NAEZ,NEMB,NATYP,LMAX,IRM,INS,INSREF,
     &                    NREF,IRNS,NCLS,NLAYER,
     &                    KREL,LMAXD,NSPIND,NAEZD,NATYPD,NREFD,NCLSD,
     &                    NEMBD,NPRINCD,KNOSPH,IRMD,IRNSD,KORBIT)
       IMPLICIT NONE
       INTEGER KREL,LMAXD,NSPIND,NAEZD,NATYPD,NREFD,NCLSD,NEMBD,NPRINCD
       INTEGER KNOSPH,IRMD,IRNSD,KORBIT
C
       INTEGER NSPIN,NAEZ,NEMB,NATYP,LMAX,IRM,INS,INSREF,NCLS
       INTEGER NREF,NLAYER,STOP_MARK
       INTEGER IRNS(*)
       INTEGER I,J
       LOGICAL TEST,OPT
       EXTERNAL TEST,OPT
c
c ---> dimension tests
c
      write(6,2050) 
c
      stop_mark=0
C------------ although NSPIND is fixed to 1 in REL mode,
C             NSPIN should be used as 1 or 2 at this stage
C             to indicate a non- or spin-polarised potential
C             that has to be read in.
C
      if((nspin.gt.nspind).and.(krel.eq.0)) then
        write(6,*) 'Please, change the parameter nspind in',
     *   ' inc.p to',nspin
        stop_mark=1
      endif
      if(naez.ne.naezd) then
        write(6,*) 'Please, change the parameter naezd (',naezd,
     +       ') in inc.p to',naez
        stop_mark=1
      endif
      if(nemb.gt.nembd) then
        write(6,*) 'Please, change the parameter nembd (',nembd,
     +       ') in inc.p to',nemb
        stop_mark=1
      endif
      IF ( .not. OPT('VIRATOMS') ) THEN
        if(natyp.ne.natypd .and. .not.opt('no test ')) then
          write(6,*) 'Please, change the parameter natypd in',
     *   ' inc.p to',natyp
          stop_mark=1
        endif
      ELSE
        if(natyp+1.ne.natypd .and. .not.opt('no test ')) then
          write(6,*) 'Please, change the parameter natypd in',
     *   ' inc.p to for virtual atom option',natyp+1
          stop_mark=1
        endif
      END IF! ( .not. OPT('VIRATOMS') ) THEN
      if(lmax.ne.lmaxd) then
        write(6,*) 'Please, change the parameter lmaxd (',lmaxd,
     +       ') in inc.p to',lmax
        stop_mark=1
      endif
      if(irm.gt.irmd) then
        write(6,*) 'Please, change the parameter irmd (',irmd,
     +       ') in inc.p to',irm
        stop_mark=1
      endif
      if(max(ins,insref).gt.knosph) then
        write(6,*) 'Please, change the parameter insd in',
     *   ' inc.p to',max(ins,insref)
        stop_mark=1
      endif
      if(nref.gt.nrefd) then
        write(6,*) 'Please, change the parameter nrefd in',
     *   ' inc.p to',nref
        stop_mark=1
      endif

      J=1
      DO I=1,NATYP 
        J=MAX(J,IRNS(I)) 
      ENDDO
      IF (INS.EQ.0 .AND. J.GT.1) THEN
        WRITE(6,*) 
     +       'IRNS(*) is set to 1 in case of ',
     +       'spherical potential treatment.'
        DO I=1,NATYP 
          IRNS(I) = 1
        ENDDO
        J = 1
      END IF
c
      if(j.gt.irnsd) then
        write(6,*) 'Please, change the parameter irnsd in',
     *   ' inc.p to',j
        stop_mark=1
      endif
c
      if(ncls.gt.nclsd) then
c        write(6,*) 'Please, change the parameter nclsd in',
c     *   ' inc.p to',ncls
c        stop_mark=1
      endif
c
c      if(ispard.lt.0) then
c        write(6,*) 'Please, check the parameter ifulld and islabd in',
c     *   ' inc.p.'
c        stop_mark=1
c      endif
c
      IF ( .not. OPT('VIRATOMS') ) THEN
        if(nref.gt.natyp) then
          write(6,*) 'There are some inconsistencies
     &               in the input file./',
     &          ' nref(=',nref,') is greater than natyp (=',natyp,').'
          stop_mark=1
        endif
      END IF

      IF((KREL.EQ.1).AND.(KORBIT.EQ.1)) then
        WRITE(6,*) 'Full relativistic for ASA and new SO solver',
     &   'KREL',KREL,'KORBIT',KORBIT
        stop_mark=1
      ENDIF

      IF(.NOT.OPT('NEWSOSOL').AND.KORBIT.EQ.1) then
        WRITE(6,*) 
     &        'Option NEWSOSOL not found, change KORBIT in inc.p from',
     &   KORBIT,'to 0' 
        stop_mark=1
      ENDIF
      
      IF(OPT('NEWSOSOL').AND.KORBIT.EQ.0) then
        WRITE(6,*) 'Using option NEWSOSOL, change KORBIT in inc.p from',
     &   KORBIT,'to 1' 
        stop_mark=1
      ENDIF
C ------------------------------------------------------------------------
c
c ---> dimension tests for matrix inversion
c
c      IF (OPT('full inv') .NEQV. (IFULLD.EQ.1)) THEN
c        write(6,*) 'Please, change the parameter ifulld in',
c     *       ' inc.p to',1-ifulld
c        stop_mark=1
c      END IF

c      IF (OPT('SPARSE  ') .NEQV. (ISPARD.EQ.1)) THEN
c        write(6,*) 'Please, change the parameter ispard in',
c     *       ' inc.p to',1-ispard
c        stop_mark=1
c      END IF

c      IF (OPT('SPARSE  ') .AND. (NPRINCD.NE.1) ) THEN
c        write(6,*) 'Please, change the parameter nprincd in',
c     *       ' inc.p to 1.'
c        stop_mark=1
c      END IF

c      IF (OPT('wfct    ') .AND. .NOT.OPT('full inv') ) THEN
c        write(6,*) 'Please, use option ''full inv'' for correct ',
c     +             'determination of eigenvectors.'
c        stop_mark=1
c      END IF

c
c ---> OPT 'WIRE' is only useful with OPT 'full inv' or 
c      OPT 'SPARSE  ' because of sparsity of
c      the KKR matrix ( not tridiagonal like for 2D and 3D systems)
c
      IF (OPT('WIRE    ') .AND. 
     +     .NOT. (OPT('full inv') .or. OPT('SPARSE  ') )) THEN
        write(6,*) 'Use option ''full inv'' or ''SPARSE  '' ',
     +       'for WIRE calculation.'
        stop_mark=1
      END IF
c
      IF (OPT('COMPLEX ') .and.
     +     .not.( OPT('EigenV  ') .or. 
     +            OPT('wfct    ') .or. 
     +            OPT('iso surf')     )   ) THEN
        write(6,*) 
     +     'Use option ''COMPLEX '' only for eigenvalue determination.'
        stop_mark=1
      END IF
c
Cccc      IF (OPT('WIRE    ')) THEN
Cccc        LATT = 11
Cccc        write(6,*) 
Cccc     +       'WIRE calculation. LATT is set to ',LATT,'.'
Cccc      END IF
c
      IF (TEST('CONT    ')) THEN
        NEMB = 0
        write(6,*) 
     +       'No usage of embedding points. NEMB is set to ',NEMB,'.'
      END IF
c
      IF (.NOT.OPT('full inv').AND. .NOT.OPT('SPARSE  ')) THEN
c
c --->  constants for O(N) algorithm for matrix inversion
c
        NLAYER=NAEZ/NPRINCD
        WRITE(6,2020) NPRINCD,NLAYER
        WRITE(6,2112)
        IF (NLAYER*NPRINCD.NE.NAEZ) THEN
          write(6,*) 'NLAYER*NPRINCD ( = ',NLAYER*NPRINCD,
     +         ').NE.NAEZ ( = ',NAEZ,')'
          stop_mark=1
        END IF

      END IF
c ------------------------------------------------------------------------
c
C ---> STOP IF A DIMENSION ERROR OCCURED
c
      if (stop_mark.gt.0) STOP 'STOP : Dimension Error.'
 2020 FORMAT(' NPRINCD  NLAYER'/,2I8)
 2112 format( 2(7(1H-),1H+) ,63(1H-))  
 2050 FORMAT(' Dimension and Input Data CHECK')    
      RETURN
      END
