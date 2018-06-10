SUBROUTINE rholm(den,df,gmat,nsra,rho2ns,drdi,ipan,ircut,pz,fz,  &
        qz,sz,cleb,icleb,iend,jend,ekl)
!-----------------------------------------------------------------------
!     calculate in the paramagnetic case (nspin=1) :
!         the valence charge density times r**2 from the greensfunction
!     calculate in the spin-polarized case (nspin=2) :
!         the valence charge density times r**2 and the valence spin
!         density times r**2 from the greensfunction ,
!         ( convention spin density :=
!                            density(spin up)-density(spin down) )
!     calculate the valence density of states , in the spin-polarized
!      case spin dependent ; splitted into its l-contributions .

!     in this subroutine an implicit energy-spin integration is  done :
!        this subroutine is called for each energy and spin value
!        and n(r,e) times df (the energy weight) is calculated .

!     recognize that the density of states is always complex also in
!      the case of "real-energy-integation" (ief>0) since in that case
!      the energy integration is done parallel to the real energy axis
!      but not on the real energy axis .
!      in the paramagnetic case only rho2ns(irmd,lmxtsq,natypd,1)
!      is used containing  the charge density times r**2 .
!      in the spin-polarized case rho2ns(...,1) contains the charge
!      density times r**2 and rho2ns(...,2) the spin density times
!      r**2 .

!     the charge density is expanded in spherical harmonics :

!             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)

!          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
!                                                           unit sphere)
!     in the case of spin-polarization :
!       the spin density is developed in spherical harmonics :

!            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)

!         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
!                                                           unit sphere)
!     n(r,e) is developed in

!        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }

!     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
!     has to be used to calculate the lm-contribution of the charge
!     density .
!             (see notes by b.drittler)

!     attention : the gaunt coeffients are stored in an index array
!                 (see subroutine gaunt)
!                 the structure part of the greens-function (gmat) is
!                 symmetric in its lm-indices , therefore only one
!                 half of the matrix is calculated in the subroutine
!                 for the back-symmetrisation . the gaunt coeffients
!                 are symmetric too (since the are calculated for
!                 real spherical harmonics) . that is why the lm2-
!                 loop only goes up to lm1 and the summands are
!                 multiplied by a factor of 2 in the case of lm1
!                 not equal to lm2 .

!                               b.drittler   may 1987
!                                   changed  dec 1988
!-----------------------------------------------------------------------
!     .. Parameters ..
INCLUDE 'inc.p'

! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
      INTEGER LMMAXD
      parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1= LMAXD+1)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO= (0.0D0,0.0D0))
!..
!.. Scalar Arguments ..
      DOUBLE COMPLEX DF
      INTEGER IEND,IPAN,NSRA
!..
!.. Array Arguments ..
DOUBLE COMPLEX DEN(0:LMAXD1),EKL(0:LMAXD),FZ(IRMD,0:LMAXD), &
               GMAT(LMMAXD,LMMAXD),PZ(IRMD,0:LMAXD), &
               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD)
DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD)
INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAND),JEND(LMPOTD,0:LMAXD,0:LMAXD)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX FFZ,GMATL,PPZ
      DOUBLE PRECISION C0LL,FACSYM,PI
      INTEGER I,J,J0,J1,L,L1,L2,LM3,LM3MAX,LN1,LN2,LNE,LNS
!..
!.. Local Arrays ..
      DOUBLE COMPLEX DENR(IRMD),WR(IRMD,0:LMAXD,0:LMAXD)
!..
!.. External Subroutines ..
      EXTERNAL CSIMPK
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,DIMAG,SQRT
!..
!.. Save statement ..
      SAVE
!..
pi = 4.0D0*ATAN(1.0D0)
c0ll = 1.0D0/SQRT(4.0D0*pi)


lm3max = icleb(iend,3)


!---> set up of wr(ir,l1,l2) = pz(ir,l1)*pz(ir,l2)

IF (nsra == 2) THEN
  DO  l1 = 0,lmaxd
    DO  l2 = 0,l1
      DO  i = 2,ircut(1)
        wr(i,l1,l2) = pz(i,l1)*pz(i,l2) + fz(i,l1)*fz(i,l2)
      END DO
    END DO
  END DO
  
ELSE
  
  DO  l1 = 0,lmaxd
    DO  l2 = 0,l1
      DO  i = 2,ircut(1)
        wr(i,l1,l2) = pz(i,l1)*pz(i,l2)
      END DO
    END DO
  END DO
  
END IF

!---> first calculate only the spherically symmetric contribution

DO  l = 0,lmaxd
  gmatl = czero
  lns = l*l + 1
  lne = lns + 2*l
  DO  ln1 = lns,lne
    gmatl = gmatl + gmat(ln1,ln1)
  END DO
  
!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
  
  denr(1) = czero
  IF (nsra == 2) THEN
    DO  i = 2,ircut(1)
      ppz = pz(i,l)
      ffz = fz(i,l)
      denr(i) = ppz* (gmatl*ppz+ekl(l)*qz(i,l)) +  &
          ffz* (gmatl*ffz+ekl(l)*sz(i,l))
      rho2ns(i,1) = rho2ns(i,1) + c0ll*DIMAG(df*denr(i))
    END DO
    
  ELSE
    
    DO  i = 2,ircut(1)
      ppz = pz(i,l)
      denr(i) = ppz* (gmatl*ppz+ekl(l)*qz(i,l))
      rho2ns(i,1) = rho2ns(i,1) + c0ll*DIMAG(df*denr(i))
    END DO
  END IF
  
  
!---> calculate density of states
  
  CALL csimpk(denr,den(l),ipan,ircut,drdi)
END DO
den(lmaxd1) = 0.0D0

!---> calculate the non spherically symmetric contribution
!        to speed up the pointer jend generated in gaunt is used
!        remember that the wavefunctions are l and not lm dependent

j0 = 1

DO  i = 1,ircut(1)
  denr(i) = 0.0D0
END DO
DO  lm3 = 2,lm3max
  DO  l1 = 0,lmaxd
    DO  l2 = 0,l1
      
      j1 = jend(lm3,l1,l2)
      
      IF (j1 /= 0) THEN
        
        gmatl = czero
        
!---> sum over m1,m2 for fixed lm3,l1,l2
        
        DO  j = j0,j1
          facsym = 2.0D0
          ln1 = icleb(j,1)
          ln2 = icleb(j,2)
          IF (ln1 == ln2) facsym = 1.0D0
          gmatl = gmatl + facsym*cleb(j)*df*gmat(ln2,ln1)
        END DO
        
        j0 = j1 + 1
        
        DO  i = 2,ircut(1)
          rho2ns(i,lm3) = rho2ns(i,lm3) + DIMAG(gmatl*wr(i,l1,l2))
        END DO
        
      END IF
      
    END DO
    
  END DO
  
END DO

END SUBROUTINE rholm
