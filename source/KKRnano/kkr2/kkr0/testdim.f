       SUBROUTINE TESTDIM(nspin,naez,lmax,
     &                    nref,
     &                    LMAXD,NREFD,NSPIND)
       implicit none

       INTEGER LMAXD
       INTEGER NREFD
       INTEGER NSPIND

       INTEGER nspin,naez,lmax
       INTEGER nref,stop_mark
c
c ---> dimension tests
c
      write(6,2050) 
c
      stop_mark=0
      if(nspin.ne.nspind) then
        write(6,*) 'Please, change the parameter nspind in',
     *   ' inc.p to',nspin
        stop_mark=1
      endif
c
      if(lmax.ne.lmaxd) then
        write(6,*) 'Please, change the parameter lmaxd (',lmaxd,
     +       ') in inc.p to',lmax
        stop_mark=1
      endif
c
      if(nref.ne.nrefd) then
        write(6,*) 'Please, change the parameter nrefd in',
     *   ' inc.p to',nref
        stop_mark=1
      endif

c
      if(nref.gt.NAEZ) then
        write(6,*) 'There are some inconsistencies in the input file./',
     *   ' nref(=',nref,') is greater than NAEZ (=',NAEZ,').'
        stop_mark=1
      endif
C-----------------------------------------------------------------------
C
C ---> STOP IF A DIMENSION ERROR OCCURED
C
      if (stop_mark.gt.0) STOP 'STOP : Dimension Error.'
 2050 FORMAT(' Dimension and Input Data CHECK')
      RETURN
      END
