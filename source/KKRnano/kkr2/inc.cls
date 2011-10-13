C ======================================================================
C   parameters file for the host TBKKR package -- TB-clusters settings
C   and for the TRUNCATION (TRC) - clusters
C                                          last update:   25/11/2008
C ======================================================================
C 
C     Description of parameters follows below (have to be set by hand
C     before compilation)
C ======================================================================
C
C PARAMETER   MEANING (settings)
C ----------------------------------------------------------------------
C NCLSD        maximum number of different TB-clusters
C NACLSD       maximum number of atoms in a TB-cluster
C ======================================================================
      INTEGER NCLSD, NACLSD
      INTEGER NCLSD_TRC, NACLSD_TRC
C
      PARAMETER ( NCLSD = 1 )
      PARAMETER ( NACLSD =51 )
C ======================================================================
