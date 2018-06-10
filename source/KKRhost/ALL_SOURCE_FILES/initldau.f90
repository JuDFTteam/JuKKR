SUBROUTINE initldau(nsra,ntldau,itldau,lopt,ueff,jeff,erefldau,  &
    visp,nspin,r,drdi,z,ipan,ircut,phi,uldau)

!   *******************************************************************
!   *  Calculates ULDAU                                               *
!   *  fivos will add some comments later                             *
!   *                                                                 *
!   *  Munich, March 2003, h.ebert, v.popescu and ph.mavropoulos      *
!   *                                                                 *
!   *  It is already later, but no comments have been added yet.      *
!   *  But wait up, they're coming...    Munich, Feb.2004, Phivos     *
!   *                                                                 *
!   *******************************************************************

use mod_types, only: t_inc
IMPLICIT NONE
INCLUDE 'inc.p'

! PARAMETER definitions
INTEGER NPOTD,MMAXD
PARAMETER (NPOTD=(2*KREL+(1-KREL)*NSPIND)*NATYPD)
PARAMETER (MMAXD=2*LMAXD+1)

! Dummy arguments
INTEGER NTLDAU,NSPIN,NSRA
DOUBLE PRECISION DRDI(IRMD,NATYPD),R(IRMD,NATYPD), &
                 VISP(IRMD,NPOTD),Z(NATYPD)
DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD)
! DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
DOUBLE COMPLEX PHI(IRMD,NATYPD)
INTEGER ITLDAU(NATYPD),LOPT(NATYPD)
INTEGER IPAN(NATYPD),IRCUT(0:IPAND,NATYPD)

! Local variables
DOUBLE PRECISION AA(MMAXD,MMAXD,MMAXD,MMAXD,0:2*LMAXD), &
                 EREFLDAU(NATYPD),FACT(0:100),FCLMB(0:2*LMAXD+1), &
                 G12,G34,JEFF(NATYPD),PI,RLOP, &
                 RPW(IRMD,2*LMAXD+1),SCL,SG(IRMD),SL(IRMD), &
                 SUM,SUMFCLMB,TG(IRMD),TL(IRMD),UEFF(NATYPD), &
                 WGTFCLMB(0:2*LMAXD+1),WIG3J,WINT(IRMD),W2(IRMD)
DOUBLE PRECISION CGCRAC,GAUNTC
DOUBLE PRECISION ATAN,DBLE
INTEGER I1,IM1,IM2,IM3,IM4,IPAN1,IR,IRC1,IT,KK, &
        L1,LF,LFMAX,LL,M1,M2,M3,M4
INTEGER NINT


!      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

pi = 4.d0*ATAN(1.0D0)
fact(0) = 1.0D0
DO i1 = 1,100
  fact(i1) = fact(i1-1)*DBLE(i1) ! factorial
END DO
IF(t_inc%i_write>0) WRITE (1337,'(/,79(1H=),/,22X,A,/,79(1H=))')  &
    'LDA+U:  INITIALISE Coulomb matrix U'

! -> Calculate test functions Phi. Phi is already normalised to
!    int phi**2 dr =1, thus it also contains factor r.

! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
!                                                        Loop over atoms
!                                         which need LDA+U ( LOPT >= 0 )
DO it = 1,ntldau
  i1 = itldau(it)
  IF ( lopt(i1)+1 == 0 ) STOP ' this atom should be LDA+U'
  CALL phicalc(i1,lopt(i1),visp,ipan,ircut,r,drdi,z,  &
      erefldau(i1),phi(1,i1),nspin,nsra)
END DO
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
IF(t_inc%i_write>0) WRITE (1337,'(6X,43(1H-),/,6X,A,/,6X,43(1H-))')  &
    'Slater integrals F^n'
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
DO it = 1,ntldau
  i1 = itldau(it)
  ipan1 = ipan(i1)
  irc1 = ircut(ipan1,i1)
  
! -> define r**l in array rpw:
  
  lfmax = 2*lopt(i1)
  DO ir = 2,irmd
    rpw(ir,1) = r(ir,i1)
    DO l1 = 2,2*lmaxd + 1
      rpw(ir,l1) = r(ir,i1)*rpw(ir,l1-1)
    END DO
  END DO
  
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  
! 1.  Calculate slater integrals FCLMB
!     (here using only the large component of sra wavefunction.
!     Whoever wants can also do the small component.)
  
! 1a. Calculate slater integrals
  
  rlop = DBLE(lopt(i1))
  sumfclmb = 0.d0
! ----------------------------------------------------------------------
  DO lf = 2,lfmax,2
    tl(1) = 0.0D0
    tg(1) = 0.0D0
    
! Note on integrand:
! Integrals are up to IRC1 because we integrate in sphere,
! without thetas.
! In case of cell integration, from IRCUT(1)+1 to IRC1 a convolution
! with thetas and gaunts is needed:
! Int dr R_l(r)**2 Sum_L' Gaunt_{lm,lm,l'm'}*thetas_{l'm'}(r).
! But then the result is m-dependent. Here, the easy way is used!
    
    DO ir = 2,irc1
      wint(ir) = dreal( DCONJG(phi(ir,i1)) * phi(ir,i1) )
      w2(ir) = 2.d0 * drdi(ir,i1) * wint(ir)
      tl(ir) = w2(ir) * rpw(ir,lf)
      tg(ir) = w2(ir) / rpw(ir,lf+1)
    END DO
    
    CALL soutk(tl,sl,ipan(i1),ircut(0,i1))
    CALL soutk(tg,sg,ipan(i1),ircut(0,i1))
    
    sl(1) = 0.0D0
    DO ir = 2,irc1
      sl(ir) = sl(ir)/rpw(ir,lf+1) + (sg(irc1)-sg(ir))*rpw(ir,lf)
    END DO
    
    sg(1) = 0.0D0
    
! See Note on integrand above.
    
    DO ir = 2,irc1
      sg(ir) = wint(ir) * sl(ir)
    END DO
    
    CALL simpk(sg,fclmb(lf),ipan1,ircut(0,i1),drdi(1,i1))
    
    wig3j = (-1)**nint(rlop) * (1D0/SQRT(2D0*rlop+1D0))  &
        * cgcrac(fact,rlop,DBLE(lf),rlop,0D0,0D0,0D0)
    
    wgtfclmb(lf) = ((2*rlop+1)/(2*rlop))*wig3j**2
    sumfclmb = sumfclmb + wgtfclmb(lf)*fclmb(lf)
  END DO
! ----------------------------------------------------------------------
  
! 1b.   Normalise slater integrals FCLMB
  
  scl = jeff(i1)/sumfclmb
  fclmb(0) = ueff(i1)
  
  DO lf = 2,lfmax,2
    fclmb(lf) = scl*fclmb(lf)
  END DO
  
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  
! ======================================================================
  
! 2.   Calculate coefficient matrix AA = a(m1,m2,m3,m4)
  
  CALL rinit(mmaxd*mmaxd*mmaxd*mmaxd*(2*lmaxd+1),aa)
  ll = lopt(i1)
  DO lf = 0,lfmax,2
    DO m3 = -ll,ll
      im3 = ll + m3 + 1
      DO m2 = -ll,ll
        im2 = ll + m2 + 1
        DO m1 = -ll,ll
          im1 = ll + m1 + 1
          m4 = m1 - m2 + m3
          IF ( -ll <= m4 .AND. m4 <= ll ) THEN
            im4 = ll + m4 + 1
            sum = 0.d0
            
            DO kk = -lf,lf
              g12 = gauntc(fact,ll,m1,lf,kk,ll,m2)
              g34 = gauntc(fact,ll,m3,lf,-kk,ll,m4)
              sum = sum + g12*g34*(-1)**ABS(kk)
            END DO
            
            aa(im1,im2,im3,im4,lf) = sum*4.d0*pi/(2.d0*DBLE(lf)+1.d0)
          END IF
        END DO
      END DO
    END DO
  END DO
! ======================================================================
  
! UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
! 3.  Calculate ULDAU
  
  CALL rinit(mmaxd*mmaxd*mmaxd*mmaxd,uldau(1,1,1,1,i1))
  
  DO lf = 0,lfmax,2
    DO im4 = 1,2*ll + 1
      DO im3 = 1,2*ll + 1
        DO im2 = 1,2*ll + 1
          DO im1 = 1,2*ll + 1
            uldau(im1,im2,im3,im4,i1) = uldau(im1,im2,im3,im4,i1)  &
                + aa(im1,im2,im3,im4,lf)*fclmb(lf)
          END DO
        END DO
      END DO
    END DO
  END DO
  
! UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
  
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  IF(t_inc%i_write>0) THEN
    WRITE (1337,'(/,8X,A,I3,/)') 'ATOM: ',i1
    WRITE (1337,'(12X,A,F8.4,A)') 'LDA+U reference energy :',  &
        erefldau(i1),' Ry'
    WRITE (1337,'(12X,A,2F8.4,A,/)') 'Ueff and Jeff = ',ueff(i1),  &
        jeff(i1),' Ry'
    WRITE (1337,'(12X,"Scaling factor for F^n :",F10.6,/)') scl
    WRITE (1337,'(12X,"  n   F^n calculated   F^n scaled ")')
    DO lf = 2,lfmax,2
      WRITE (1337,'(12X,I3,2(2X,F12.8," Ry"))') lf,fclmb(lf)/scl,fclmb(lf)
    END DO
    IF ( it < ntldau ) WRITE(1337,'(8X,58(1H~))')
  END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  
END DO
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

IF(t_inc%i_write>0) THEN
  WRITE (1337,'(/,6X,60(1H-),/,6X,A,/,6X,60(1H-))')  &
      'Coulomb matrix U(m1,m1,m3,m3)'
  DO it = 1,ntldau
    i1 = itldau(it)
    ll = lopt(i1)
    ll = MIN(3,ll)
    WRITE(1337,'(/,8X,"ATOM :",I3,/)') i1
    
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    DO im1 = 1,2*ll + 1
      WRITE (1337,99001) (uldau(im1,im1,im3,im3,i1),im3=1,2*ll+1)
    END DO
    IF ( it < ntldau ) WRITE(1337,'(8X,58(1H~))')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
  END DO
  WRITE (1337,'(/,6X,60(1H-),/)')
END IF
99001 FORMAT(6X,7F10.6)
END SUBROUTINE initldau

!*==cgcrac.f    processed by SPAG 6.05Rc at 15:27 on  7 Mar 2003

FUNCTION cgcrac(fact,j1,j2,j3,m1,m2,m3)
!   ********************************************************************
!   *                                                                  *
!   *     CLEBSCH GORDAN COEFFICIENTS FOR ARBITRARY                    *
!   *     QUANTUM NUMBERS  J1,J2 ...                                   *
!   *     ACCORDING TO THE FORMULA OF   RACAH                          *
!   *     SEE: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM          *
!   *          EQUATION (3.19)                                         *
!   *          EDMONDS EQ. (3.6.11) PAGE 45                            *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

!Dummy arguments
DOUBLE PRECISION J1,J2,J3,M1,M2,M3
DOUBLE PRECISION CGCRAC
DOUBLE PRECISION FACT(0:100)

! Local variables
DOUBLE PRECISION DSQRT
INTEGER J,N,N1,N2,N3,N4,N5,NBOT,NTOP
INTEGER NINT
DOUBLE PRECISION RFACT
DOUBLE PRECISION S,SUM,VF,X,Y

! INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
rfact(x) = fact(nint(x))


cgcrac = 0.0D0
IF ( ABS(m3-(m1+m2)) > 1.0D-6 ) RETURN
IF ( ABS(j1-j2) > j3 ) RETURN
IF ( (j1+j2) < j3 ) RETURN
IF ( ABS(m1) > (j1+1.0D-6) ) RETURN
IF ( ABS(m2) > (j2+1.0D-6) ) RETURN
IF ( ABS(m3) > (j3+1.0D-6) ) RETURN

DO j = ABS(nint(2*(j1-j2))),nint(2*(j1+j2)),2
  IF ( j == nint(2*j3) ) GO TO 100
END DO
RETURN


100  CONTINUE
x = (2.0D0*j3+1.0D0)*rfact(j1+j2-j3)*rfact(j1-j2+j3)  &
    *rfact(-j1+j2+j3)*rfact(j1+m1)*rfact(j1-m1)*rfact(j2+m2)  &
    *rfact(j2-m2)*rfact(j3+m3)*rfact(j3-m3)

y = rfact(j1+j2+j3+1)

vf = DSQRT(x/y)


n1 = nint(j1+j2-j3)
n2 = nint(j1-m1)
n3 = nint(j2+m2)
n4 = nint(j3-j2+m1)
n5 = nint(j3-j1-m2)
ntop = MIN(n1,n2,n3)
nbot = MAX(0,-n4,-n5)

n = nbot + 1
IF ( n == (2*(n/2)) ) THEN
  s = +1.0D0
ELSE
  s = -1.0D0
END IF
sum = 0.0D0

DO n = nbot,ntop
  s = -s
  y = fact(n)*fact(n1-n)*fact(n2-n)*fact(n3-n)*fact(n4+n) *fact(n5+n)
  sum = sum + (s/y)
END DO
cgcrac = vf*sum
END FUNCTION cgcrac
!*==gauntc.f    processed by SPAG 6.05Rc at 15:27 on  7 Mar 2003

FUNCTION gauntc(fact,l1,m1,l2,m2,l3,m3)
!   ********************************************************************
!   *                                                                  *
!   *     GAUNT COEFFICIENTS for complex spherical harmonics  Y[l,m]   *
!   *                                                                  *
!   *            G = INT dr^  Y[l1,m1]* Y[l2,m2] Y[l3,m3]              *
!   *                                                                  *
!   * see: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM  Eq. (4.34)  *
!   *                                                                  *
!   * 26/01/95  HE                                                     *
!   ********************************************************************
IMPLICIT NONE

! PARAMETER definitions
DOUBLE PRECISION PI
PARAMETER (PI=3.141592653589793238462643D0)

! Dummy arguments
INTEGER L1,L2,L3,M1,M2,M3
DOUBLE PRECISION FACT(0:100)
DOUBLE PRECISION GAUNTC

! Local variables
DOUBLE PRECISION CGCRAC
DOUBLE PRECISION DBLE
DOUBLE PRECISION G

IF ( (l1 < 0) .OR. (l2 < 0) .OR. (l3 < 0) ) THEN
  g = 0.0D0
ELSE
  g = (DBLE(2*l2+1)*DBLE(2*l3+1)/(4.0D0*pi*DBLE(2*l1+1)))  &
      **0.5D0*cgcrac(fact,DBLE(l3),DBLE(l2),DBLE(l1),DBLE(m3),  &
      DBLE(m2),DBLE(m1))*cgcrac(fact,DBLE(l3),DBLE(l2),DBLE(l1),  &
      0.0D0,0.0D0,0.0D0)
END IF
gauntc = g
END FUNCTION gauntc
