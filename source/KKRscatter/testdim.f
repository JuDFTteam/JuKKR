       SUBROUTINE TESTDIM(nspin,nspo,nspoh,naez,nemb,natyp,lmax,irm,ins,
     &                    insref,nref,M2,IRNS,ncls,nlayer,LATT)
       implicit none
       include 'inc.p'
       include 'inc.cls'
       INTEGER nspin,nspo,naez,nemb,nz,natyp,lmax,irm,ins,insref,ncls
       INTEGER M2,nspoh
       INTEGER nref,nlayer,stop_mark
       INTEGER IRNS(*)
       INTEGER I,J,LATT
       LOGICAL TEST,OPT
       EXTERNAL TEST,OPT
c
c ---> dimension tests
c
      write(6,2050) 
c
      stop_mark=0
      if(nspin.gt.nspind) then
        write(6,*) 'Please, change the parameter nspind in',
     *   ' inc.p to',nspin
        stop_mark=1
      endif
      if(nspo.ne.nspod.and.nspoh.ne.nspod) then
        write(6,*) 'Please, change the parameter nspod (',nspod,
     +       ') in inc.p to',nspo
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
      if(natyp.ne.natypd .and. .not.opt('no test ')) then
        write(6,*) 'Please, change the parameter natypd in',
     *   ' inc.p to',natyp
        stop_mark=1
      endif
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
      if(max(ins,insref).gt.insd) then
        write(6,*) 'Please, change the parameter insd in',
     *   ' inc.p to',max(ins,insref)
        stop_mark=1
      endif
      if(nref.ne.nrefd) then
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
      if(ncls.ne.nclsd) then
        write(6,*) 'Please, change the parameter nclsd in',
     *   ' inc.cls to',ncls
        stop_mark=1
      endif
c
c      if(ispard.lt.0) then
c        write(6,*) 'Please, check the parameter ifulld and islabd in',
c     *   ' inc.p.'
c        stop_mark=1
c      endif
c
      if(nref.gt.natyp) then
        write(6,*) 'There are some inconsistencies in the input file./',
     *   ' nref(=',nref,') is greater than natyp (=',natyp,').'
        stop_mark=1
      endif
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
      IF (OPT('WIRE    ')) THEN
        LATT = 11
        write(6,*) 
     +       'WIRE calculation. LATT is set to ',LATT,'.'
      END IF
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







