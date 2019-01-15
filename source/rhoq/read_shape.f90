! **********************************************************************
      SUBROUTINE read_shape(THETAS_out,lmsp_out,ifunm_out,IRID,IRMD, &
     &                      NFUND,LPOT,mu0,IPAND,NTCELL,NATYP)
! **********************************************************************
!   reads the input potentials
!
!    units :       ry - units for energy
!                  the lattice constant and all other lengths
!                                           given in bohr units
!                  the planck constant h/2pi=1
!                  the electron charge e=sqrt(2)
!                  the electron mass m=1/2
!                  the speed of light ! = 2/alpha = 274.0720442
!                      with the fine structure constant alpha
!
!    remember that the input potentials do not include the electro-
!             static contribution of the nucleus of the cell itself
!             this has to be added explicitly !
!
!   as input is used: lmax=maximum angular momentum
!                    nbeg .. nend=number of different atoms
!
!
!     in case of shape corrections this routine  reads from unit 19
!     a suitable radial  mesh 'xrn',its derivate 'drn' and the shape
!     functions 'thetas' .          thus, the region from the muffin
!     tin to the circumscribed  sphere radii is divided  into  'npan'
!     pannels, each one containing 'nm(ipan)' points in order to take
!     care of the  discontinuities of the shape-function  derivative.
!     at the output one obtains :
!            llmsp (icell,ifun)       = integer array giving the com-
!                                       posite  index  lm=l*(l+1)+m+1
!                                       of the ifun-th shape function
!            lmsp  (icell,lm)         = (0,1)  if the lm-th component
!                                       is vanishing or not
!            nfu   (icell)            = number  of   shape   function
!                                       components in cell 'icell'
!
!
!     modified for bandstructure code
!
!                                 b.drittler nov. 1989
! **********************************************************************
      IMPLICIT NONE
!     .. Array Arguments ..
      DOUBLE PRECISION, intent(out) :: THETAS_out(IRID,NFUND)
      INTEGER :: IFUNM_OUT((2*LPOT+1)**2), LMSP_OUT((2*LPOT+1)**2),NTCELL(NATYP)
!     .. Scalar Arguments ..
      INTEGER, intent(in) :: IRID, IRMD, NFUND, LPOT, mu0, IPAND, NATYP
!     .. Local Scalars ..
      INTEGER :: NCELL, IPAN1, IR, ICELL, IFUN, LM, N, NFUN
!     .. Local Arrays ..
      DOUBLE PRECISION, allocatable :: DRN(:,:),SCALE(:),U(:),XRN(:,:),THETAS(:,:,:)
      INTEGER, allocatable :: MESHN(:),NM(:,:),NPAN(:),NFU(:),IFUNM(:),LMSP(:)
!-----------------------------------------------------------------------

      OPEN (19,FILE='shapefun_mu0',STATUS='old',FORM='formatted')

!---> read radial mesh information of the shape functions and
!     shape functions THETAS in the first iteration - if needed

        READ (19,FMT=9000) NCELL
        WRITE (1337,FMT=*) '  ncell : ',NCELL
        
      ! allocate arrays depending on NCELL  
      ALLOCATE( DRN(IRID,NCELL),SCALE(NCELL),U(IRMD), XRN(IRID,NCELL) )
      ALLOCATE( MESHN(NCELL),NM(IPAND,NCELL),NPAN(NCELL), NFU(NCELL) )
      ALLOCATE( THETAS(IRID,NFUND,NCELL) )
      ALLOCATE( IFUNM((2*LPOT+1)**2), LMSP((2*LPOT+1)**2) )

        READ (19,FMT=9010) (SCALE(ICELL),ICELL=1,NCELL)
        DO ICELL = 1,NCELL
          READ (19,FMT=9000) NPAN(ICELL),MESHN(ICELL)

          IF(NPAN(ICELL)+1.GT.IPAND) THEN
            WRITE(6,*) 'Please, change the parameter ipand (',IPAND,') in inc.p to',NPAN(ICELL)+1
            STOP 'STARTB - IPAND'
          ENDIF

          IF(MESHN(ICELL).GT.IRID) THEN
            WRITE(6,*) 'Please, change the parameter irid (',IRID,') in inc.p to',MESHN(ICELL)
            STOP 'STARTB - IRID'
          ENDIF

          READ (19,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
          READ (19,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1, MESHN(ICELL))

          READ (19,FMT=9000) NFU(ICELL)
          NFUN = NFU(ICELL)
          WRITE (1337,FMT=*) '  nfun  : ',NFUN,NFUND

          IF(NFUN.GT.NFUND) THEN
            WRITE(6,*) 'Please, change the parameter nfund (',NFUND,') in inc.p to',NFUN
            STOP 'STARTB - NFUND'
          ENDIF

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
            if(NTCELL(mu0)==ICELL) THETAS_OUT(1:MESHN(ICELL),IFUN) = THETAS(1:MESHN(ICELL),IFUN,ICELL)
            if(NTCELL(mu0)==ICELL) lmsp_out(lm) = lmsp(lm)
            if(NTCELL(mu0)==ICELL) ifunm_out(lm) = ifunm(lm)
          END DO !IFUN=1,NFUN

        END DO !ICELL=1,NCELL
        
        DEALLOCATE( DRN, SCALE, U, XRN, MESHN, NM, NPAN, THETAS, lmsp, ifunm )
        
        close(19)
      
!-----------------------------------------------------------------------

      RETURN

 9000 FORMAT (16i5)
 9010 FORMAT (4d20.12)
 
      END                           ! STARTB1
