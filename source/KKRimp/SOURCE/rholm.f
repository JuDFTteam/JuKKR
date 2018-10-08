      MODULE MOD_RHOLM
      CONTAINS
!-------------------------------------------------------------------------
!> Summary: Driver for valence charge density for spherical potential
!> Category: physical-observables, KKRimp
!>
!> calculate in the paramagnetic case (nspin=1) :
!>     the valence charge density times r**2 from the greensfunction
!> calculate in the spin-polarized case (nspin=2) :
!>     the valence charge density times r**2 and the valence spin
!>     density times r**2 from the greensfunction ,
!>     ( convention spin density :=
!>                        density(spin up)-density(spin down) )
!> calculate the valence density of states , in the spin-polarized
!>  case spin dependent ; splitted into its l-contributions .
!>
!>   in this subroutine an implicit energy-spin integration is  done :
!>    this subroutine is called for each energy and spin value
!>    and n(r,e) times df (the energy weight) is calculated .
!>
!>  recognize that the density of states is always complex also in
!>  the case of "real-energy-integation" (ief>0) since in that case
!>  the energy integration is done parallel to the real energy axis
!>  but not on the real energy axis .
!>  in the paramagnetic case only rho2ns(irmd,lmxtsq,natypd,1)
!>  is used containing  the charge density times r**2 .
!>  in the spin-polarized case rho2ns(...,1) contains the charge
!>  density times r**2 and rho2ns(...,2) the spin density times
!>  r**2 .
!>
!>  the charge density is expanded in spherical harmonics :
!>
!>           rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)
!>
!>        rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
!>                                                       unit sphere)
!> in the case of spin-polarization :
!>   the spin density is developed in spherical harmonics :
!>
!>          sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)
!>
!>       sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
!>                                                       unit sphere)
!> n(r,e) is developed in
!>
!>      n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }
!>
!>   therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
!> has to be used to calculate the lm-contribution of the charge
!> density .
!>         (see notes by b.drittler)
!>
!>   attention : the gaunt coeffients are stored in an index array
!>             (see subroutine gaunt)
!>             the structure part of the greens-function (gmat) is
!>             symmetric in its lm-indices , therefore only one
!>             half of the matrix is calculated in the subroutine
!>             for the back-symmetrisation . the gaunt coeffients
!>             are symmetric too (since the are calculated for
!>             real spherical harmonics) . that is why the lm2-
!>             loop only goes up to lm1 and the summands are
!>             multiplied by a factor of 2 in the case of lm1
!>             not equal to lm2 .
!>
!>                             b.drittler   may 1987
!>                               changed  dec 1988
!>
!> For KREL = 1 (relativistic mode)
!>  NPOTD = 2 * NATYPD             
!>  LMMAXD = 2 * (LMAXD+1)^2       
!>  NSPIND = 1                     
      SUBROUTINE RHOLM(DEN,DF,GMAT,NSRA,RHO2NS,DRDI,IPAN,IRCUT,PZ,FZ,
     +                   QZ,SZ,CLEB,ICLEB,IEND,JEND,EKL,
     +                   IRMD,NCLEB,LMAXD,LMMAXD,LMPOTD)

      USE MOD_CSIMPK
      IMPLICIT NONE

      INTEGER IRMD,NCLEB
      INTEGER LMAXD,LMMAXD
!       parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
!       INTEGER LMAXD1
!       PARAMETER (LMAXD1= LMAXD+1)
      INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO= (0.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF
      INTEGER IEND,IPAN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DEN(0:LMAXD+1),EKL(0:LMAXD),FZ(IRMD,0:LMAXD),
     +               GMAT(LMMAXD,LMMAXD),PZ(IRMD,0:LMAXD),
     +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD)
      DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD)
      INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAN),JEND(LMPOTD,0:LMAXD,0:LMAXD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX FFZ,GMATL,PPZ
      DOUBLE PRECISION C0LL,FACSYM,PI
      INTEGER I,J,J0,J1,L,L1,L2,LM3,LM3MAX,LN1,LN2,LNE,LNS
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX DENR(IRMD),WR(IRMD,0:LMAXD,0:LMAXD)
C     ..
C     .. External Subroutines ..
!       EXTERNAL CSIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DIMAG,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
      PI = 4.0D0*ATAN(1.0D0)
      C0LL = 1.0D0/SQRT(4.0D0*PI)
c

      LM3MAX = ICLEB(IEND,3)

c
c---> set up of wr(ir,l1,l2) = pz(ir,l1)*pz(ir,l2)
c
      IF (NSRA.EQ.2) THEN
        DO 30 L1 = 0,LMAXD
          DO 20 L2 = 0,L1
            DO 10 I = 2,IRCUT(1)
              WR(I,L1,L2) = PZ(I,L1)*PZ(I,L2) + FZ(I,L1)*FZ(I,L2)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE

      ELSE

        DO 60 L1 = 0,LMAXD
          DO 50 L2 = 0,L1
            DO 40 I = 2,IRCUT(1)
              WR(I,L1,L2) = PZ(I,L1)*PZ(I,L2)
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE

      END IF
c
c---> first calculate only the spherically symmetric contribution
c
      DO 100 L = 0,LMAXD
        GMATL = CZERO
        LNS = L*L + 1
        LNE = LNS + 2*L
        DO 70 LN1 = LNS,LNE
          GMATL = GMATL + GMAT(LN1,LN1)
   70   CONTINUE
c
c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
        DENR(1) = CZERO
        IF (NSRA.EQ.2) THEN
          DO 80 I = 2,IRCUT(1)
            PPZ = PZ(I,L)
            FFZ = FZ(I,L)
            DENR(I) = PPZ* (GMATL*PPZ+EKL(L)*QZ(I,L)) +
     +                FFZ* (GMATL*FFZ+EKL(L)*SZ(I,L))
            RHO2NS(I,1) = RHO2NS(I,1) + C0LL*DIMAG(DF*DENR(I))
   80     CONTINUE

        ELSE

          DO 90 I = 2,IRCUT(1)
            PPZ = PZ(I,L)
            DENR(I) = PPZ* (GMATL*PPZ+EKL(L)*QZ(I,L))
            RHO2NS(I,1) = RHO2NS(I,1) + C0LL*DIMAG(DF*DENR(I))
   90     CONTINUE
        END IF
c
c
c---> calculate density of states
c
        CALL CSIMPK(DENR,DEN(L),IPAN,IRCUT,DRDI)
  100 CONTINUE
      DEN(LMAXD+1) = 0.0D0
c
c---> calculate the non spherically symmetric contribution
c        to speed up the pointer jend generated in gaunt is used
c        remember that the wavefunctions are l and not lm dependent
c
      J0 = 1
c
      DO 110 I = 1,IRCUT(1)
        DENR(I) = 0.0D0
  110 CONTINUE
      DO 160 LM3 = 2,LM3MAX
        DO 150 L1 = 0,LMAXD
          DO 140 L2 = 0,L1
c
            J1 = JEND(LM3,L1,L2)
c
            IF (J1.NE.0) THEN
c
              GMATL = CZERO
c
c---> sum over m1,m2 for fixed lm3,l1,l2
c
              DO 120 J = J0,J1
                FACSYM = 2.0D0
                LN1 = ICLEB(J,1)
                LN2 = ICLEB(J,2)
                IF (LN1.EQ.LN2) FACSYM = 1.0D0
                GMATL = GMATL + FACSYM*CLEB(J)*DF*GMAT(LN2,LN1)
  120         CONTINUE
c
              J0 = J1 + 1
c
              DO 130 I = 2,IRCUT(1)
                RHO2NS(I,LM3) = RHO2NS(I,LM3) + DIMAG(GMATL*WR(I,L1,L2))
  130         CONTINUE

            END IF

  140     CONTINUE

  150   CONTINUE

  160 CONTINUE

      END SUBROUTINE RHOLM
      END MODULE MOD_RHOLM
