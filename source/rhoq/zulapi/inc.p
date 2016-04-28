c number of different atoms in unitcell
      integer natypd
      parameter (NATYPD =26)
c number of atoms in unitcell
      integer naezd
      parameter (NAEZD =NATYPD)
c highest valence orbital quantum number     
      integer lmaxd  
      parameter (LMAXD = 2)
c maximum number of reference potentials
      integer nrefd
      parameter (NREFD = 1)
c number of spin directions (nspind= 1 non-spin-polarized
c                                    2  spin-polarized    )
      integer NSPIND
      parameter (NSPIND = 2)
c spin-orbit coupling yes/no  (nspo= 1 no spin-orbit coupl
c                                    2 with spin-orbit coupl)
      integer NSPOD
      parameter (NSPOD = 2)
c combined parameter for SOC and NSPIND 
      integer NSPD                   
      parameter (NSPD=MAX(NSPOD,NSPIND))
c number of layers in one principal layer
      integer NPRINCD
      parameter (NPRINCD = 2)
c parameter for non-/spherical potentials (0/1)
      integer INSD
      parameter (INSD = 1)
c dimension of KKR-matrix for full matrix inversion
      integer ndimgk
      parameter (NDIMGK = NAEZD)
c number of principal layers
      integer nlayerd
      parameter (NLAYERD = (NAEZD/NPRINCD))
c modified number of atoms for storing non-sph. radial functions
      integer nsatypd
      parameter (NSATYPD = (NATYPD-1)*INSD+1)
c modified number of potentials for storing non-sph. potentials
      integer nspotd
      parameter (NSPOTD = NSPIND*NSATYPD)
c number of shells (for storing GF in GS (*,*,*,nsheld))
      integer nsheld
      parameter (NSHELD = 1500)
c number of symmetries of the system optimizing reduces 
c memory size significantly
       integer NSYMD
       parameter (NSYMD=6 )
c highest potential orbital quantum number     
      integer lpotd
      parameter (LPOTD = 2*LMAXD)
c number of clebsch gordon coeffizients
      integer ncleb
      parameter (NCLEB = ((lmaxd*2+1)**2)*((lmaxd+1)**2) )
c number of points in complex energy plane for integration
      integer iemxd
      parameter (IEMXD = 3  )
c number of r points in (0  , rmt (? rws))
      integer irmd
c      parameter (IRMD = 353)
      parameter (IRMD = 600)
c number of r points in ( rmt,rws)
      integer irnsd
c      parameter (IRNSD =1)
      parameter (IRNSD =500)    
c     ..
      integer ijd
      parameter(ijd=434)
c     ..   
      integer nfund,irid,ngshd,ngfd
c     parameter (NFUND = 35,IRID = 200,NGSHD = 3000,NGFD= 1)
!      parameter (NFUND = 35,IRID = 145,NGSHD = 1500,NGFD= 1)
      parameter (NFUND = 100,IRID = 200,NGSHD = 1500,NGFD= 1)
c      parameter (NFUND = 49,IRID = 135,NGSHD = 3500,NGFD= 1)
c ..  number of cells (shapes), panels in non-spherical part
      integer ncelld,ipand
      parameter (NCELLD = 24,IPAND = 24)
c number of real space and reciprocal space vectors
      integer nrd,kpoibz
c      parameter (NRD = 10000,KPOIBZ =3000000)
      parameter (NRD = 10000,KPOIBZ =300000)
c
c number of principal layers for SPARSE matrix inversion
c now put to naezd, Better NOT change this
c 
      integer nlspd
      parameter (NLSPD = NAEZD)
c
c nauxspd: nonzero blocks in G(L,L,k,E) for SPARSE matrix inversion 
c You can put to naezd*nprincd*3 for the case of principal layer
c here I disconect the sparce matrix from the principal layer
c 
      integer nauxspd
      parameter (NAUXSPD =50*NAEZD)
c parameter nembd for embeding positions around the real atomic pos.
      integer nembd
      parameter (NEMBD = 20 )
c Parameter for Impurity Calculations
      integer natomimpd
      PARAMETER (NATOMIMPD=300)

