C **********************************************************************
      SUBROUTINE read_shape(THETAS_out,lmsp_out,ifunm_out,IRID,IRMD,
     &                      NFUND,LPOT,mu0,IPAND,NTCELL,NATYP)
C **********************************************************************
c   reads the input potentials
c
c    units :       ry - units for energy
c                  the lattice constant and all other lengths
c                                           given in bohr units
c                  the planck constant h/2pi=1
c                  the electron charge e=sqrt(2)
c                  the electron mass m=1/2
c                  the speed of light c = 2/alpha = 274.0720442
c                      with the fine structure constant alpha
c
c    remember that the input potentials do not include the electro-
c             static contribution of the nucleus of the cell itself
c             this has to be added explicitly !
c
c   as input is used: lmax=maximum angular momentum
c                    nbeg .. nend=number of different atoms
c
c
c     in case of shape corrections this routine  reads from unit 19
c     a suitable radial  mesh 'xrn',its derivate 'drn' and the shape
c     functions 'thetas' .          thus, the region from the muffin
c     tin to the circumscribed  sphere radii is divided  into  'npan'
c     pannels, each one containing 'nm(ipan)' points in order to take
c     care of the  discontinuities of the shape-function  derivative.
c     at the output one obtains :
c            llmsp (icell,ifun)       = integer array giving the com-
c                                       posite  index  lm=l*(l+1)+m+1
c                                       of the ifun-th shape function
c            lmsp  (icell,lm)         = (0,1)  if the lm-th component
c                                       is vanishing or not
c            nfu   (icell)            = number  of   shape   function
c                                       components in cell 'icell'
c
c
c     modified for bandstructure code
c
c                                 b.drittler nov. 1989
C **********************************************************************
C     .. Parameters ..
      IMPLICIT NONE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION, intent(out) :: THETAS_out(IRID,NFUND)
      INTEGER :: IFUNM_OUT((2*LPOT+1)**2), LMSP_OUT((2*LPOT+1)**2), 
     &           NTCELL(NATYP)
C     ..
C     .. Scalar Arguments ..
      INTEGER, intent(in) :: IRID, IRMD, NFUND, LPOT, mu0, IPAND, NATYP
C     ..
C     .. Local Scalars ..
      INTEGER :: NCELL, IPAN1, IR, ICELL, IFUN, LM, N, NFUN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION, allocatable :: DRN(:,:),SCALE(:),U(:),XRN(:,:),
     &                                 THETAS(:,:,:)
      INTEGER, allocatable :: MESHN(:),NM(:,:),NPAN(:),NFU(:),IFUNM(:),
     &                        LMSP(:)
C     ..
c-----------------------------------------------------------------------

      OPEN (19,FILE='shapefun',STATUS='old',FORM='formatted')
c
c---> read radial mesh information of the shape functions and
c     shape functions THETAS in the first iteration - if needed
c
        READ (19,FMT=9000) NCELL
        WRITE (1337,FMT=*) '  ncell : ',NCELL
        
      ! allocate arrays depending on NCELL  
      ALLOCATE( DRN(IRID,NCELL),SCALE(NCELL),U(IRMD),
     +                 XRN(IRID,NCELL) )
      ALLOCATE( MESHN(NCELL),NM(IPAND,NCELL),NPAN(NCELL), NFU(NCELL) )
      ALLOCATE( THETAS(IRID,NFUND,NCELL) )
      ALLOCATE( IFUNM((2*LPOT+1)**2), LMSP((2*LPOT+1)**2) )
c
        READ (19,FMT=9010) (SCALE(ICELL),ICELL=1,NCELL)
        DO ICELL = 1,NCELL
          READ (19,FMT=9000) NPAN(ICELL),MESHN(ICELL)
c
          IF(NPAN(ICELL)+1.GT.IPAND) THEN
            WRITE(6,*) 'Please, change the parameter ipand (',IPAND,
     +           ') in inc.p to',NPAN(ICELL)+1
            STOP 'STARTB - IPAND'
          ENDIF
c
          IF(MESHN(ICELL).GT.IRID) THEN
            WRITE(6,*) 'Please, change the parameter irid (',IRID,
     +           ') in inc.p to',MESHN(ICELL)
            STOP 'STARTB - IRID'
          ENDIF
c
          READ (19,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
          READ (19,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1,
     +      MESHN(ICELL))

          READ (19,FMT=9000) NFU(ICELL)
          NFUN = NFU(ICELL)
          WRITE (1337,FMT=*) '  nfun  : ',NFUN,NFUND
c
          IF(NFUN.GT.NFUND) THEN
            WRITE(6,*) 'Please, change the parameter nfund (',NFUND,
     +           ') in inc.p to',NFUN
            STOP 'STARTB - NFUND'
          ENDIF
c
          DO LM = 1,(2*LPOT+1)**2
            LMSP(LM) = 0
          END DO ! LM

          DO IFUN = 1,NFUN
            READ (19,FMT=9000) LM
!             LLMSP(ICELL,IFUN) = LM
            LMSP(LM) = 1
            IFUNM(LM) = IFUN
            READ (19,FMT=9010) (THETAS(N,IFUN,ICELL),N=1,MESHN(ICELL))
            ! save only the shapefunction that is needed
!             write(*,*) MU0,ICELL,ifun,thetas(50,ifun,icell)
            if(NTCELL(mu0)==ICELL) THETAS_OUT(1:MESHN(ICELL),IFUN) = 
     &                                THETAS(1:MESHN(ICELL),IFUN,ICELL)
            if(NTCELL(mu0)==ICELL) lmsp_out(lm) = lmsp(lm)
            if(NTCELL(mu0)==ICELL) ifunm_out(lm) = ifunm(lm)
          END DO !IFUN=1,NFUN

        END DO !ICELL=1,NCELL
        
!         stop
   
        DEALLOCATE( DRN, SCALE, U, XRN, MESHN, NM, NPAN, THETAS, 
     &              lmsp, ifunm )
        
        close(19)
      
c-----------------------------------------------------------------------

      RETURN

 9000 FORMAT (16i5)
 9010 FORMAT (4d20.12)
 
      END                           ! STARTB1
