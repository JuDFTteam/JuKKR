!-----------------------------------------------------------------------
!>    Calculate reference system's T-matrix.
!>    @param     E complex energy
!>    @param     VREF repulsive reference potential field strength
!>    @param     LMAX angular momentum cutoff
!>    @param     RMTREF repulsive reference potential muffin-tin radius
!>    @param     TREFLL reference system T-matrix
!>    @param     DTREFLL energy derivative of reference system T-matrix
!>               only calculated if LLY=1
!>    @param     LLY do Lloyd's formula calculation 0=no/1=yes
subroutine TREF(E,VREF,LMAX,RMTREF,TREFLL,DTREFLL, &
LLY) ! new input parameter after inc.p removal

  implicit none

  integer LLY
  !     ..
  !     .. Scalar Arguments ..
  integer LMAX
  !     ..
  !     .. Array Arguments ..
  double precision RMTREF,VREF
  !     ..
  !     .. Local Scalars ..
  double complex A1,B1,DA1,DB1,E
  integer I,L,LM1
  !     ..
  !     .. Local Arrays ..
  double complex BESSJW1(0:LMAX+1),BESSJW2(0:LMAX+1), &
  BESSYW1(0:LMAX+1),BESSYW2(0:LMAX+1), &
  HANKWS1(0:LMAX+1),HANKWS2(0:LMAX+1)

  double complex TREFLL((LMAX+1)**2, (LMAX+1)**2)

  double complex DBESSJW1(0:LMAX+1),DBESSJW2(0:LMAX+1), &
  DHANKWS1(0:LMAX+1)

  double complex DTREFLL((LMAX+1)**2, (LMAX+1)**2)

  !     .. External Subroutines ..
  external BESSEL
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic SQRT
  !     ..
  !     .. Save statement ..
  !     SAVE
  integer LMAXD
  integer LMGF0D

  !-----------------------------------------------------------------------
  !---- t-matrix and derivative of t-matrix of the reference system

  !     the analytical formula for the derivative of spherical Bessel
  !     functions is used:

  !     d                     l+1
  !     --  j (x) = j   (x) - --- j (x)
  !     dx   l       l-1       x   l

  !     d
  !     --  j (x) = - j (x)
  !     dx   0         1

  !     which for x = sqrt(E0)*r leads to

  !      d          r*r             (l+1)
  !     --- j (x) = --- ( j   (x) - ----- j (x) )
  !     dE0  l      2 x    l-1        x    l

  !      d            r*r
  !     --- j (x) = - --- j (x)
  !     dE0  0        2 x  1

  !-----------------------------------------------------------------------

  LMAXD = LMAX
  LMGF0D= (LMAXD+1)**2

  A1 = SQRT(E)*RMTREF
  B1 = SQRT(E-VREF)*RMTREF
  call BESSEL(BESSJW1,BESSYW1,HANKWS1,A1,LMAXD+1)
  call BESSEL(BESSJW2,BESSYW2,HANKWS2,B1,LMAXD+1)

  if(LLY.eq.1) then

    DBESSJW1(0) = - BESSJW1(1)/A1
    DBESSJW2(0) = - BESSJW2(1)/B1
    DHANKWS1(0) = - HANKWS1(1)/A1

    do L = 1,LMAX + 1
      DBESSJW1(L) = (BESSJW1(L-1) - (L+1)*BESSJW1(L)/A1)/A1
      DBESSJW2(L) = (BESSJW2(L-1) - (L+1)*BESSJW2(L)/B1)/B1
      DHANKWS1(L) = (HANKWS1(L-1) - (L+1)*HANKWS1(L)/A1)/A1
    end do

    do L = 0,LMAX + 1
      DBESSJW1(L) = 0.5D0*DBESSJW1(L)*RMTREF**2
      DBESSJW2(L) = 0.5D0*DBESSJW2(L)*RMTREF**2
      DHANKWS1(L) = 0.5D0*DHANKWS1(L)*RMTREF**2
    end do

  endif

  do L = 0,LMAX
    A1 = SQRT(E)*BESSJW1(L+1)*BESSJW2(L) - &
    SQRT(E-VREF)*BESSJW1(L)*BESSJW2(L+1)

    B1 = SQRT(E)*HANKWS1(L+1)*BESSJW2(L) - &
    SQRT(E-VREF)*HANKWS1(L)*BESSJW2(L+1)

    do I = -L,L
      LM1 = L* (L+1) + I + 1
      TREFLL(LM1,LM1) = -1.D0/SQRT(E)*A1/B1
    end do

  end do


  call CINIT(LMGF0D*LMGF0D,DTREFLL)

  if(LLY.eq.1) then

    do L = 0,LMAX
      A1 = SQRT(E)*BESSJW1(L+1)*BESSJW2(L) - &
      SQRT(E-VREF)*BESSJW1(L)*BESSJW2(L+1)

      DA1 = 0.5D0/SQRT(E)*BESSJW1(L+1)*BESSJW2(L) - &
      0.5D0/SQRT(E-VREF)*BESSJW1(L)*BESSJW2(L+1) + &
      SQRT(E)*DBESSJW1(L+1)*BESSJW2(L) - &
      SQRT(E-VREF)*DBESSJW1(L)*BESSJW2(L+1) + &
      SQRT(E)*BESSJW1(L+1)*DBESSJW2(L) - &
      SQRT(E-VREF)*BESSJW1(L)*DBESSJW2(L+1)

      B1 = SQRT(E)*HANKWS1(L+1)*BESSJW2(L) - &
      SQRT(E-VREF)*HANKWS1(L)*BESSJW2(L+1)

      DB1 = 0.5D0/SQRT(E)*HANKWS1(L+1)*BESSJW2(L) - &
      0.5D0/SQRT(E-VREF)*HANKWS1(L)*BESSJW2(L+1) + &
      SQRT(E)*DHANKWS1(L+1)*BESSJW2(L) - &
      SQRT(E-VREF)*DHANKWS1(L)*BESSJW2(L+1) + &
      SQRT(E)*HANKWS1(L+1)*DBESSJW2(L) - &
      SQRT(E-VREF)*HANKWS1(L)*DBESSJW2(L+1)

      do I = -L,L
        LM1 = L* (L+1) + I + 1
        DTREFLL(LM1,LM1) = 0.5D0/SQRT(E)**3*A1/B1 &
        - 1.D0/SQRT(E)*(DA1/B1-A1*DB1/B1**2)
      end do

    end do
  endif

end


!------------------------------------------------------------------------------
subroutine GREF(E,ALATC,IEND,NCLS,NAEZ, &
CLEB,RCLS,ATOM,CLS,ICLEB,LOFLM,NACLS, &
TREFLL,DTREFLL,GREFN,DGREFN, &
LLY_G0TR, &
lmaxd, naclsd, ncleb, nclsd, &
LLY)

  use SingleSiteRef_mod
  implicit none

  !     .. Parameters ..

  integer  lmaxd
  integer  naclsd
  integer  ncleb
  integer  nclsd
  integer  LLY

  !     ..
  !     .. Scalar Arguments ..
  double precision ALATC
  integer          IEND,NCLS,NAEZ
  !     ..
  !     .. Array Arguments ..
  double precision CLEB(NCLEB,2),RCLS(3,NACLSD,NCLSD)
  integer          ATOM(NACLSD,NAEZ),CLS(NAEZ),ICLEB(NCLEB,3)

  integer          LOFLM((2*LMAXD+1)**2)
  integer          NACLS(NCLSD)
  double complex   LLY_G0TR(NCLSD)

  double complex   TREFLL((LMAXD+1)**2,(LMAXD+1)**2), &
  DTREFLL((LMAXD+1)**2,(LMAXD+1)**2)
  double complex   DGREFN((LMAXD+1)**2,(LMAXD+1)**2,NACLSD,NCLSD), &
  GREFN((LMAXD+1)**2,(LMAXD+1)**2,NACLSD,NCLSD)

  !     ..
  !     .. Local Scalars ..
  double complex   E
  integer          I1,IC,ICLS,IG,IG1,LM1,LM2
  !     ..

  !     Local arrays
  !     Fortran 90 automatic arrays - medium size

  double complex   DGINP(NACLSD*(LMAXD+1)**2,(LMAXD+1)**2)
  double complex   GINP(NACLSD*(LMAXD+1)**2,(LMAXD+1)**2)

  integer          LMGF0D
  integer          LMMAXD

  LMMAXD= (LMAXD+1)**2
  LMGF0D= (LMAXD+1)**2

  ! attention in this subroutine I3 labels the fixed atom - I1 is a variable !

  !=====================================================================
  do ICLS = 1,NCLS
  !=====================================================================

      I1 = 1
      IC = 0
      ! find atom index corresponding to cluster index
      do while (IC.eq.0 .and. I1.le.NAEZ)
        if (CLS(I1).eq.ICLS) IC = I1
        I1 = I1 + 1
      end do
      if (IC.eq.0) stop 'Error in CLS(*) array in tbref'

      call GLL95(E,CLEB(1,2),ICLEB,LOFLM,IEND,TREFLL,DTREFLL, &
      ATOM(1,IC),RCLS(1,1,ICLS),NACLS(ICLS), &
      ALATC,GINP,DGINP, &
      LLY_G0TR(ICLS), &
      lmaxd, naclsd, ncleb, LLY )

      do IG=1,NACLSD
        do LM2=1,LMGF0D
          IG1 = (IG-1)*LMGF0D + LM2
          do LM1=1,LMGF0D
            GREFN(LM2,LM1,IG,ICLS)=GINP(IG1,LM1)
            DGREFN(LM2,LM1,IG,ICLS)=DGINP(IG1,LM1)
          enddo
        enddo
      enddo

  !=====================================================================
  end do
  !=====================================================================

end subroutine
