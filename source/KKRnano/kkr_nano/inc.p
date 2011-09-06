C ======================================================================
C              parameters file for the host TBKKR package
C                                          last update:   06/07/2009
C ======================================================================
C 
C Description of parameters follows below
C
C a) parameters that have to be set by hand before compilation
C ============================================================
C
C PARAMETER   MEANING (settings)
C ----------------------------------------------------------------------
C KREL         switch for non-relativistic (KREL=0) or relativistic 
C              (KREL=1) program. Attention: several other parameters
C              depend explicitly on KREL, they are set automatically
C KNOSPH       switch for spherical (KNOSPH=0) or non-spherical 
C              (KNOSPH=1) program. Same obs. as for KREL applies.
C KSP          switch for spin-polarised (KSP=1) or paramagnetic 
C              (KSP=0) calculation. It determines the value of NSPIND 
C LMAXD        cut-off for the angular momentum expansion
C IEMXD        dimension for energy-dependent arrays 
C IRMD         number of radial mesh points in (0,...,RWS)
C IRNSD        number of radial mesh points in (RMT,...,RWS)
C NRD          number of real space 
C KPOIBZ       and reciprocal space vectors
C NMAXD,ISHLD  paremeters for the Ewald summations
C NTREFD       parameter in broyden subroutine
C              MUST BE 0 for the host program 
C NAEZD        number of different atomic sites in the unit cell
C NREFD        number of reference potentials
C NXIJD        number of Jij taken maximally into account
C NCELLD       number of cells (shapes) in non-spherical part
C IPAND        number of panels in non-spherical part
C NFUND,IRID
C NGSHD        shape functions parameters in non-spherical part 
C IGUESSD      if 1 arrays for initial guess of the QMR-iteration 
C              are allocated ( otherwise 0 )   
C
C
C b) parameters derived from the values of the previous ones
C ============================================================
C
C PARAMETER   MEANING (settings)
C ----------------------------------------------------------------------
C NSPIND      number of spin directions, depending on KSP and KREL
C LPOTD       highest expansion in potential 
C NCLEB       number of Clebsch-Gordon coefficients
C NSPOTD      number of potentials for storing non-sph. potentials
C NTPERD      parameter in broyden subroutines
C NGUESSD     dimension of the array used for initial guess  
C
C =====================================================================
C                           fixed parameters
C =====================================================================
      INTEGER KREL,KNOSPH,KSP,LMAXD,IEMXD,IRMD,IRNSD,
     &        NRD,KPOIBZ,NMAXD,ISHLD,NTREFD,NAEZD,
     &        NREFD,NCELLD,IPAND,NFUND,IRID,
     &        NGSHD,LLY,ITDBRYD,BCPD,IGUESSD,
     &        XDIM,YDIM,ZDIM,NATBLD,
     &        EKMD,LMPID,SMPID,EMPID,NTHRDS,NXIJD,
     &        TRC,NATRCD,NUTRCD
C ---------------------------------------------------------------------
C     general settings
      PARAMETER ( KREL = 0 )
      PARAMETER ( KNOSPH = 1 )
      PARAMETER ( KSP = 0 )
      PARAMETER ( LMAXD = 2 )
      PARAMETER ( IEMXD = 45)
      PARAMETER ( IRMD = 484, IRNSD = 208 )
      PARAMETER ( NRD = 6000, KPOIBZ = 820 )
      PARAMETER ( NMAXD = 20000, ISHLD= 500 )
      PARAMETER ( NTREFD = 0 )             ! must be 0 for host program
      PARAMETER ( LLY = 0 )              ! LLY = 1(0), (in)active
C ---------------------------------------------------------------------
C     structure-dependent
      PARAMETER ( NAEZD = 4 )
      PARAMETER ( NREFD =  1 )
      PARAMETER ( NXIJD = 100 )
      PARAMETER ( TRC    = 0   )
      PARAMETER ( NATRCD = 341 )
      PARAMETER ( NUTRCD = 341 )
C ---------------------------------------------------------------------
C     non-spherical potential 
      PARAMETER ( NCELLD = 1, IPAND = 5 )
      PARAMETER ( NFUND = 15, IRID = 135, NGSHD = 13079)
C
C ---------------------------------------------------------------------
C     mixing - parallelized  
      PARAMETER ( ITDBRYD = 30 )           ! choose equal ITDBRY in inp
C
C ---------------------------------------------------------------------
C     PRECONDITIONING
      PARAMETER ( BCPD = 0)
      PARAMETER ( IGUESSD = 1)
      PARAMETER ( EKMD =  2000 )
      PARAMETER ( NATBLD = 4 )                     ! number of atoms per block
      PARAMETER ( XDIM = 1 , YDIM = 1 , ZDIM = 1 ) ! number of blocks in x,y,z direction
C ---------------------------------------------------------------------
C     L-parallelization
      PARAMETER ( LMPID = 1 )
      PARAMETER ( SMPID = 1 )
      PARAMETER ( EMPID = 2 )
      PARAMETER ( NTHRDS = 1 )
C =====================================================================
C                         derived parameters
C =====================================================================
      INTEGER NSPIND,LPOTD,NCLEB,NSATYPD,NSPOTD,NTPERD,NGUESSD,
     &        IELLYD,LRECTRC,NFLEXD
C
      PARAMETER ( NSPIND = KREL + (1-KREL)*(KSP+1) )
      PARAMETER ( LPOTD = 2*LMAXD )
      PARAMETER ( NCLEB = (LMAXD*2+1)**2 * (LMAXD+1)**2 )
      PARAMETER ( NSATYPD = (NAEZD-1)*KNOSPH+1,
     &            NSPOTD = (2*KREL + (1-KREL)*NSPIND) * NSATYPD )
      PARAMETER ( NTPERD=NAEZD-NTREFD )
      PARAMETER ( NGUESSD = 1 + IGUESSD * ( NAEZD * (LMAXD+1)**2 - 1 ) )
      PARAMETER ( LRECTRC   = 4*(NATRCD*3+2) )
      PARAMETER ( NFLEXD =  NAEZD + TRC*(-NAEZD + NUTRCD) )
      PARAMETER ( IELLYD = IEMXD * ( 1 + 2*LLY ) )
C =====================================================================
