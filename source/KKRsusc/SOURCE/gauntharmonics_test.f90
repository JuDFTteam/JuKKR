module mod_gauntharmonics_test
  use nrtype
  use mod_config, only: config_testflag
  use type_gauntcoeff
  type(GAUNTCOEFF_TYPE),allocatable        :: gauntcoeff(:)
!---------------------------------------------------------------------
! old comment out
!---------------------------------------------------------------------
!   type                               :: harmonicstype
!      real(kind=dp),allocatable       :: wg(:)
!      real(kind=dp),allocatable       :: yrg(:,:,:)
!   end type harmonicstype
!   type(HARMONICSTYPE),protected    :: harmonics

contains
!---------------------------------------------------------------------
!-- SUBROUTINE : gauntharmonics_set  
!--   sets calculates the gaunt coefficients and sph. harmonics
!---------------------------------------------------------------------
subroutine gauntharmonics_set_test()!,lmaxatom)
  implicit none
!local variables
  real(kind=dp),allocatable               :: wg(:),yrg(:,:,:)
  integer,allocatable                     :: icleb(:,:)                !: pointer array
  integer,allocatable                     :: loflm(:)                  !: l of lm=(l,m) (gaunt)
  real(kind=dp),allocatable               :: cleb(:,:)                 !: gaunt coefficients (gaunt)
  integer                                 :: iend                      !: number of nonzero gaunt coeffizients
  integer,allocatable                     :: jend(:,:,:)               !: pointer array for icleb()
  integer                                 :: ncleb 
  integer                                 :: lval
  integer                                 :: lmaxbounds(2)
  integer                                 :: lmaxd,lpotd,lmmaxd,lmpotd,lm2d
  integer                                 :: iatom

lmaxbounds=0
! call GAUNTHARMONICS_getlmaxbounds(lmaxatom,natom,lmaxbounds)
lmaxbounds(1)=2
lmaxbounds(2)=4
allocate( gauntcoeff( lmaxbounds(1) : lmaxbounds(2) ) )

write(1337,*) '------------------------------------------------------'
write(1337,*) '------- MODULE: GAUNTHAROMICS           --------------'
write(1337,*) '------------------------------------------------------'
write(1337,*) ' min. LMAX VALUE = ',lmaxbounds(1)
write(1337,*) ' max. LMAX VALUE = ',lmaxbounds(2)
write(1337,*) ' creating GAUNTCOEFF in this range'


do lval = lmaxbounds(1), lmaxbounds(2)
  lmaxd  = lval
  lpotd  = 2*lmaxd
  lmmaxd = (lmaxd+1)**2
  lmpotd = (lpotd+1)**2
  lm2d   = (2*lmaxd+1)**2
  ncleb  = 100*(lmaxd*2+1)**2 * (lmaxd+1)**2

  allocate( icleb(ncleb,4), cleb(ncleb,2),&
            jend(lmpotd,0:lmaxd,0:lmaxd), loflm(lm2d),&
            wg(4*lmaxd), yrg(4*lmaxd,0:4*lmaxd,0:4*lmaxd) )

  icleb=0;cleb=0;jend=0;loflm=0;wg=0;yrg=0

  call gauntharmonics_gaunt2(wg,yrg,4*lmaxd)

  call gauntharmonics_gaunt1(lmaxd,lpotd,wg,yrg,cleb,loflm,&
                             icleb,iend,jend, ncleb,lmaxd,lmmaxd,lmpotd)

  allocate( GAUNTcoeff(lval)%icleb(ncleb,4), GAUNTcoeff(lval)%cleb(ncleb,2), &
            GAUNTcoeff(lval)%jend(lmpotd,0:lmaxd,0:lmaxd), GAUNTcoeff(lval)%loflm(lm2d), &
            GAUNTcoeff(lval)%wg(4*lmaxd), GAUNTcoeff(lval)%yrg(4*lmaxd,0:4*lmaxd,0:4*lmaxd) )

  GAUNTcoeff(lval) % cleb = cleb
  GAUNTcoeff(lval) % loflm = loflm
  GAUNTcoeff(lval) % icleb = icleb
  GAUNTcoeff(lval) % iend  = iend
  GAUNTcoeff(lval) % jend = jend
  GAUNTcoeff(lval) % ncleb = ncleb
  GAUNTcoeff(lval) % wg    = wg
  GAUNTcoeff(lval)%yrg      = yrg
  deallocate( icleb, cleb, jend, loflm, wg, yrg ) 

end do !lval

!---------------------------------------------------------------------
! if a testflag is set write out the gaunt coefficient data
!---------------------------------------------------------------------
! if (config_testflag('writegaunt')) then
!   open(unit=435234,file='test_gaunt')
!   write(*,*)   'Gaunt coefficients'
!   write(*,*) 'GAUNTcoeff%ncleb'
!   write(*,*) ncleb
!   write(*,*) 'GAUNTcoeff%cleb(:,1)'
!   write(*,*) cleb(:,1)
!   write(*,*) 'GAUNTcoeff%cleb(:,2)'
!   write(*,*) cleb(:,2)
!   write(*,*) 'loflm'
!   write(*,*) loflm
!   write(*,*) 'icleb'
!   write(*,*) icleb
!   write(*,*) 'iend'
!   write(*,*) iend
!   write(*,*) 'jend'
!   write(*,*) jend
!   write(*,*) 'ncleb'
!   write(*,*) ncleb
!   close(435234)
! end if

if (config_testflag('writefilegaunt')) then
  open(unit=3453453,file='gaunt.dat') 
  write(3453453,*) 'WG'
  write(3453453,*) WG
  write(3453453,*) 'YRG'
  write(3453453,*) YRG
  write(3453453,*) 'CLEB'
  write(3453453,*) CLEB
  write(3453453,*) 'LOFLM'
  write(3453453,*) LOFLM
  write(3453453,*) 'ICLEB'
  write(3453453,*) ICLEB
  write(3453453,*) 'IEND'
  write(3453453,*) IEND
  write(3453453,*) 'JEND'
  write(3453453,*) JEND
  write(3453453,*) 'NCLEB'
  write(3453453,*) NCLEB
  close(3453453)
end if

end subroutine gauntharmonics_set_test

!---------------------------------------------------------------------
!-- SUBROUTINE : GAUNTHARMONICS_getlmaxbounds  
!---------------------------------------------------------------------
subroutine GAUNTHARMONICS_getlmaxbounds(lmaxatom,natom,lmaxbounds)
  implicit none
  integer, intent(in)       :: lmaxatom(natom)
  integer, intent(in)       :: natom
  integer, intent(inout)    :: lmaxbounds(2)
  integer                   :: iatom
lmaxbounds=0
do iatom=1,natom
  if (lmaxbounds(1)==0 .and. lmaxbounds(2)==0) lmaxbounds=lmaxatom(iatom)
  if (lmaxatom(iatom)<lmaxbounds(1)) lmaxbounds(1)=lmaxatom(iatom)
  if (lmaxatom(iatom)>lmaxbounds(2)) lmaxbounds(2)=lmaxatom(iatom)
end do !iatom
end subroutine GAUNTHARMONICS_getlmaxbounds


!---------------------------------------------------------------------
!-- SUBROUTINE : GAUNTHARMONICS_GAUNT1  
!---------------------------------------------------------------------
      SUBROUTINE GAUNTHARMONICS_GAUNT1(LMAX,LPOT,W,YR,CLEB,LOFLM,ICLEB,IEND,JEND, &
                       NCLEB,LMAXD,LMGF0D,LMPOTD)
! ************************************************************************
!
!   - fills the array cleb with the gaunt coeffients ,i.e.
!      the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
!      but only for lm2.le.lm1 and lm3>1
!   - calculate the pointer array jend  to project the indices
!      array cleb with the same lm3,l1,l2 values - because of
!      the special ordering of array cleb only the last index
!      has to be determined .
!     (the parameter n has to be chosen that l1+l2+l3 .lt. 2*n)
!     using gaussian quadrature as given by
!     m. abramowitz and i.a. stegun, handbook of mathematical functions,
!     nbs applied mathematics series 55 (1968), pages 887 and 916
!     m. weinert and e. wimmer
!     northwestern university march 1980
!
!     an index array -icleb- is used to save storage place .
!     fills the array loflm which is used to determine the
!     l-value of a given lm-value .
!     this subroutine has to be called only once !
!
!                               b.drittler   november 1987
!
!     modified gaunt coefficients are als calculated defined by
!     the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
!-----------------------------------------------------------------------
!
!---> attention : ncleb is an empirical factor - it has to be optimized
!
!     .. Parameters ..
!cccccccc      include 'inc.p'
!
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************
!
!     ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
!     ..
      INTEGER LMPOTD,LMGF0D,LMAXD,NCLEB
!     ..
!     .. Scalar Arguments ..
      INTEGER IEND,LMAX,LPOT
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION CLEB(NCLEB,2),W(*), &
                       YR(4*LMAXD,0:4*LMAXD,0:4*LMAXD)
      INTEGER ICLEB(NCLEB,4),JEND(LMPOTD,0:LMAXD,0:LMAXD),LOFLM(*)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION CLECG,FACTOR,FCI,S
      INTEGER I,J,L,L1,L1P,L2,L2P,L3,LM1,LM2,LM3,LM3P,LMPOT,M,M1,M1A, &
              M1S,M2,M2A,M2S,M3,M3A,M3S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,MOD,REAL,SIGN
!     ..
!     .. External Subroutines ..
!       EXTERNAL RCSTOP
!     ..
      I = 1
      DO L = 0,2*LMAX
        DO M = -L,L
          LOFLM(I) = L
          I = I + 1
        END DO
      END DO
!
      IF (LPOT.EQ.0) THEN
        IEND = 1
        ICLEB(1,1) = (LMAX+1)**2
        ICLEB(1,3) = 1
      END IF
!
      IF (LPOT.NE.0) THEN
!
!---> set up of the gaunt coefficients with an index field
!
        I = 1
        DO L3 = 0,LPOT
          DO M3 = -L3,L3
!
            DO L1 = 0,LMAX
              DO L2 = 0,L1
!
!                 IF (MOD((L1+L2+L3),2).NE.1 .AND. (L1+L2-L3).GE.0 .AND. &
!                     (L1-L2+L3).GE.0 .AND. (L2-L1+L3).GE.0) THEN
                IF (1 &
                      ) THEN
!
                  FCI = DBLE(CI** (L2-L1+L3))
                  DO M1 = -L1,L1
                    DO M2 = -L2,L2
!
!---> store only gaunt coeffients for lm2.le.lm1
!
                      LM1 = L1* (L1+1) + M1 + 1
                      LM2 = L2* (L2+1) + M2 + 1
                      IF (1) THEN
!
                        M1S = SIGN(1,M1)
                        M2S = SIGN(1,M2)
                        M3S = SIGN(1,M3)
!
!                         IF (M1S*M2S*M3S.GE.0) THEN
                        IF (1) THEN
!
                          M1A = ABS(M1)
                          M2A = ABS(M2)
                          M3A = ABS(M3)
!
                          FACTOR = 0.0
!
                          IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR + &
                              REAL(3*M3S+SIGN(1,-M3))/8.0d0
                          IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR + &
                              REAL(M1S)/4.0d0
                          IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR + &
                              REAL(M2S)/4.0d0
!
                          IF (1) THEN
!
                            IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR. &
                                M1S*M3S.NE.1) FACTOR = -FACTOR
!
                            S = 0.0
                            DO J = 1,4*LMAXD
                              S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)* &
                                  YR(J,L3,M3A)
                            END DO !J
                            CLECG = S*FACTOR
                            IF (1) THEN
                              CLEB(I,1) = CLECG
                              CLEB(I,2) = FCI*CLECG
                              ICLEB(I,1) = LM1
                              ICLEB(I,2) = LM2
                              ICLEB(I,3) = L3* (L3+1) + M3 + 1
                              ICLEB(I,4) = LM2*LMGF0D - &
                                           (LM2*LM2-LM2)/2 + LM1 - &
                                           LMGF0D
!                               if (ICLEB(I,1)==ICLEB(I,2)) write(*,*) ICLEB(I,1),ICLEB(I,2),ICLEB(I,3),CLECG
                              if (1==ICLEB(I,3) .and. abs(CLECG)>10e-10) write(*,*) ICLEB(I,1),ICLEB(I,2),ICLEB(I,3),CLECG
!                               if (1) write(*,*) ICLEB(I,1),ICLEB(I,2),ICLEB(I,3),CLECG
                              I = I + 1
                            END IF

                          END IF

                        END IF

                      END IF

                    END DO !M2
                  END DO !M1
                END IF

              END DO !L2
            END DO !L1
          END DO !M3
        END DO !L3
        IEND = I - 1
        IF (NCLEB.LT.IEND) THEN
          WRITE (6,FMT=9000) NCLEB,IEND
          STOP'[GAUNT.f] stop'

        ELSE
!
!---> set up of the pointer array jend,use explicitly
!     the ordering of the gaunt coeffients
!
          LMPOT = (LPOT+1)* (LPOT+1)
          DO L1 = 0,LMAX
            DO L2 = 0,L1
              DO LM3 = 2,LMPOT
                JEND(LM3,L1,L2) = 0
              END DO !LM3
            END DO !L2
          END DO !L1
!
          LM3 = ICLEB(1,3)
          L1 = LOFLM(ICLEB(1,1))
          L2 = LOFLM(ICLEB(1,2))
!
          DO J = 2,IEND
            LM3P = ICLEB(J,3)
            L1P = LOFLM(ICLEB(J,1))
            L2P = LOFLM(ICLEB(J,2))
!
            IF (LM3.NE.LM3P .OR. L1.NE.L1P .OR. L2.NE.L2P) THEN
              JEND(LM3,L1,L2) = J - 1
              LM3 = LM3P
              L1 = L1P
              L2 = L2P
            END IF

          END DO !J
          JEND(LM3,L1,L2) = IEND
!
!
        END IF

      END IF
!
!

 9000 FORMAT (13x,'error stop in gaunt : dimension of NCLEB = ',i10, &
             ' too small ',/, &
             13x,'change NCLEB to ',i6)
      END SUBROUTINE GAUNTHARMONICS_GAUNT1





!---------------------------------------------------------------------
!-- SUBROUTINE : GAUNTHARMONICS_GAUNT2  
!---------------------------------------------------------------------
      SUBROUTINE GAUNTHARMONICS_GAUNT2(W,YR,N)
! ************************************************************************
!     sets up values needed for gaunt
!        m. weinert  january 1982
!
!     changed for calculating with real spherical harmonics
!                                           b.drittler  july 1987
!
!     W(N)        integration weights on 4*LMAXD points in the intervall
!                 (-1,0) (from routine GRULE)
!
!     YR(N,L,M)   spherical harmonics on 4*LMAXD points to angular 
!                 momentum indices (l,m) scaled with a factor 
!                 of RF=(4*pi)**(1/3)
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!     .. Arguments
      INTEGER N
      DOUBLE PRECISION W(*),YR(N,0:N,0:N)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION A,CD,CTH,FAC,FPI,RF,STH,T
      INTEGER K,L,LOMAX,M
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION P(0:N+1,0:N),X(N)
!     ..
!     .. External Subroutines ..
!       EXTERNAL GRULE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
!     ..
      FPI = 16D0*ATAN(1D0)
      RF = FPI**(1D0/3D0)
      LOMAX = N
!
!--->    obtain gauss-legendre points and weights
!
      CALL GAUNTHARMONICS_GRULE(2*N,X,W)
!
!--->    generate associated legendre functions for m.ge.0
!
      DO K = 1,N
         CTH = X(K)
         STH = SQRT(1.D0-CTH*CTH)
         FAC = 1.D0
!
!--->    loop over m values
!
         DO M = 0,LOMAX
            FAC = - DBLE(2*M-1)*FAC
            P(M,M) = FAC
            P(M+1,M) = DBLE(2*M+1)*CTH*FAC
!
!--->    recurse upward in l
!
            DO L = M + 2,LOMAX
               P(L,M) = ( DBLE(2*L-1)*CTH*P(L-1,M) &
                        - DBLE(L+M-1)    *P(L-2,M) ) / DBLE(L-M)
            END DO
!
            FAC = FAC*STH
         END DO
!
!--->    multiply in the normalization factors
!
         DO L = 0,LOMAX
            A = RF*SQRT((2*L+1)/FPI)
            CD = 1.D0
            YR(K,L,0) = A*P(L,0)
!
            DO M = 1,L
               T = DBLE( (L+1-M)* (L+M))
               CD = CD/T
               YR(K,L,M) = A*SQRT(2.D0*CD)*P(L,M)
            END DO
         END DO
      END DO
      END SUBROUTINE GAUNTHARMONICS_GAUNT2




!---------------------------------------------------------------------
!-- SUBROUTINE : GAUNTHARMONICS_GRULE  
!---------------------------------------------------------------------
      SUBROUTINE GAUNTHARMONICS_GRULE(N,X,W)
!
!***********************************************************************
!
!     determines the (n+1)/2 nonnegative points x(i) and
!     the corresponding weights w(i) of the n-point
!     gauss-legendre integration rule, normalized to the
!     interval [-1,1]. the x(i) appear in descending order.
!
!     this routine is from 'methods of numerical integration',
!     p.j. davis and p. rabinowitz, page 369.
!
!***********************************************************************
!
!
!     .. Scalar Arguments ..
      INTEGER N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION W(*),X(*)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION D1,D2PN,D3PN,D4PN,DEN,DP,DPN,E1,FX,H,P,PI,PK, &
                       PKM1,PKP1,T,T1,U,V,X0
      INTEGER I,IT,K,M
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC COS,ATAN
!     ..
      PI = 4.D0*ATAN(1.D0)
      M  = (N+1)/2
      E1 = N* (N+1)
      DO I = 1,M
        T = (4*I-1)*PI/ (4*N+2)
        X0 = (1.0d0- (1.0d0-1.0d0/N)/ (8.0d0*N*N))*COS(T)
!
!--->    iterate on the value  (m.w. jan. 1982)
!
        DO IT = 1,2
          PKM1 = 1.
          PK = X0
          DO K = 2,N
            T1 = X0*PK
            PKP1 = T1 - PKM1 - (T1-PKM1)/K + T1
            PKM1 = PK
            PK = PKP1
          END DO !K
          DEN = 1. - X0*X0
          D1 = N* (PKM1-X0*PK)
          DPN = D1/DEN
          D2PN = (2.*X0*DPN-E1*PK)/DEN
          D3PN = (4.*X0*D2PN+ (2.-E1)*DPN)/DEN
          D4PN = (6.*X0*D3PN+ (6.-E1)*D2PN)/DEN
          U = PK/DPN
          V = D2PN/DPN
          H = -U* (1.+.5*U* (V+U* (V*V-U*D3PN/ (3.*DPN))))
          P = PK + H* (DPN+.5*H* (D2PN+H/3.* (D3PN+.25*H*D4PN)))
          DP = DPN + H* (D2PN+.5*H* (D3PN+H*D4PN/3.))
          H = H - P/DP
          X0 = X0 + H
        END DO !IT
        X(I) = X0
        FX = D1 - H*E1* (PK+.5*H* (DPN+H/3.* (D2PN+.25*H* (D3PN+ &
             .2*H*D4PN))))
        W(I) = 2.* (1.-X(I)*X(I))/ (FX*FX)
      END DO !I
      IF (M+M.GT.N) X(M) = 0.
      END SUBROUTINE GAUNTHARMONICS_GRULE




end module mod_gauntharmonics_test