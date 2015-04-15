      SUBROUTINE RLLSLL(RPANBOUND,RMESH,VLL,RLL,SLL,TLLP, &
                        NCHEB,NPAN,LMSIZE,LMSIZE2,LBESSEL,NRMAX,NRMAXD, &
                        NVEC,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR, &
                        CMODERLL,CMODESLL,CMODETEST,USE_SRATRICK1 )
! ************************************************************************
! radial wave functions by the integral equation method of
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! ************************************************************************
! This routine solves the following two equations:
!
! RLL(r) = J(r) - PRE * J(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * RLL(r') ) 
!               + PRE * H(r) * int_0^r( dr' r'^2 J2(r') * op(V(r')) * RLL(r') )
!
! SLL(r) = H(r) - PRE * H(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * RLL(r') ) 
!               + PRE * J(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * SLL(r') )
!
! where the integral int_0^r() runs from 0 to r
! ************************************************************************
! Potential matrix : VLL(LMSIZE*NVEC,LMSIZE*NVEC)
! ************************************************************************
! Green function prefacor PRE=GMATPREFACTOR (scalar value)
! ************************************************************************
! Source terms : J, H  (nvec*lmsize,lmsize) or (lmsize,nvec*lmsize)
!                J2,H2 (lmsize,nvec*lmsize) or (nvec*lmsize,lmsize) 
!
!
! The source term J is for LMSIZE=3 and NVEC=2 given by:
! J =      / jlk(jlk_index(1))                                          \
!          |       0            jlk(jlk_index(2))                       |
!          |       0                   0            jlk(jlk_index(3))   |
!          | jlk(jlk_index(4))                                          |
!          |       0            jlk(jlk_index(5))                       |
!          \       0                   0            jlk(jlk_index(6))   /
! ************************************************************************
! Operator op() can be chosen to be a unity or a transpose operation
!     The unity operation is used to calculate the right solution
!     The transpose operation is used to calculate the left solution
! ************************************************************************
! RMESH      - radial mesh
! RPANBOUND  - panel bounds RPANBOUND(0) left  panel border of panel 1
!                           RPANBOUND(1) right panel border of panel 1
! NCHEB      - highes chebyshev polynomial
!              number of points per panel = NCHEB + 1
! NPAN       - number of panels
! LMSIZE     - number of colums for the source matrix J etc...
! LMSIZE2    - number of rows   for the source matrix J etc...
! NRMAX      - total number of radial points (NPAN*(NCHEB+1))
! NVEC       - number of LMSIZE*LMSIZE blocks in J (LMSIZE2=NVEC*LMSIZE)
! ************************************************************************
IMPLICIT NONE
      INTEGER NCHEB,NPAN,LMSIZE,LMSIZE2,NRMAX,NRMAXD,LBESSEL
      INTEGER IVEC, IVEC2, NVEC
      DOUBLE COMPLEX CI,CONE,CZERO
      PARAMETER (CI= (0.0D0,1.0D0),CONE=(1.0D0,0.0D0))
      PARAMETER (CZERO=(0.0D0,0.0D0))
!     ..
!     .. Local Scalars ..
      DOUBLE COMPLEX GMATPREFACTOR
      DOUBLE PRECISION PLLM
      INTEGER INFO,ICHEB3,ICHEB2,ICHEB,ICHEB1,IPAN,MN,NM,MN2,NPLM
      INTEGER L1,L2,LM1,LM2,LM3
!     ..
!     .. Local Arrays ..


      DOUBLE COMPLEX :: HLK(LBESSEL,NRMAX), &
                        JLK(LBESSEL,NRMAX), &
                        HLK2(LBESSEL,NRMAX), &
                        JLK2(LBESSEL,NRMAX) 

      character(len=1) :: CMODERLL,CMODESLL,CMODETEST


      DOUBLE COMPLEX &
                     SLL(LMSIZE2,LMSIZE,NRMAX), &
                     RLL(LMSIZE2,LMSIZE,NRMAX), TLLP(LMSIZE,LMSIZE), &
                     VLL(LMSIZE*NVEC,LMSIZE*NVEC,NRMAX)


      DOUBLE COMPLEX,allocatable ::  ULL(:,:,:)


      DOUBLE COMPLEX,allocatable ::  &
                     WORK(:,:), &
                     WORK2(:,:), &
                     ALLP(:,:,:),BLLP(:,:,:), &
                     CLLP(:,:,:),DLLP(:,:,:), &
                     SLV(:,:,:,:),SRV(:,:,:,:), &
                     SLV1(:,:,:,:),SRV1(:,:,:,:), &
                     SLV2(:,:,:,:),SRV2(:,:,:,:), &
                     SLV3(:,:,:,:),SRV3(:,:,:,:), &
                     MRNVY(:,:,:),MRNVZ(:,:,:), &
                     MRJVY(:,:,:),MRJVZ(:,:,:), &
                     MIHVY(:,:,:),MIHVZ(:,:,:), &
                     MIJVY(:,:,:),MIJVZ(:,:,:), &
                     YILL(:,:,:),ZILL(:,:,:), &
                     YRLL(:,:,:),ZRLL(:,:,:),YRLLTMP(:,:,:), &
                     YILL1(:,:,:),ZILL1(:,:,:), &
                     YRLL1(:,:,:),ZRLL1(:,:,:), &
                     YILL2(:,:,:),ZILL2(:,:,:), &
                     YRLL2(:,:,:),ZRLL2(:,:,:), &
                     VHLR(:,:,:), &
                     VJLR(:,:,:), &
                     VHLI(:,:,:), &
                     VJLI(:,:,:)
      DOUBLE COMPLEX,ALLOCATABLE :: YIF(:,:,:,:), &
                     YRF(:,:,:,:), &
                     ZIF(:,:,:,:), &
                     ZRF(:,:,:,:)
      DOUBLE COMPLEX ZSLC1SUM(0:NCHEB)
      DOUBLE PRECISION C1(0:NCHEB,0:NCHEB),RPANBOUND(0:NPAN)
      DOUBLE PRECISION CSLC1(0:NCHEB,0:NCHEB),CSRC1(0:NCHEB,0:NCHEB), &
                       TAU(0:NCHEB,0:NPAN),CDDRC1(0:NCHEB,0:NCHEB), &
                       SLC1SUM(0:NCHEB),RMESH(NRMAXD)
      INTEGER JLK_INDEX(LMSIZE2),IPIV(0:NCHEB,LMSIZE2)
      INTEGER,ALLOCATABLE :: IPIV2(:)
      LOGICAL TEST
      INTEGER :: IERROR,USE_SRATRICK,USE_SRATRICK1
      INTEGER,PARAMETER  :: DIRECTSOLV=1


!     .. External Subroutines ..
      EXTERNAL ZGETRF,ZGETRS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,DIMAG,EXP,MAX,MIN,SIN,SQRT
!     ..
! use pherical potential as references or not

      ALLOCATE ( ULL(LMSIZE2,LMSIZE,NRMAX) )

    IF (LMSIZE.EQ.1) THEN
     USE_SRATRICK=0
    ELSE
     USE_SRATRICK=USE_SRATRICK1
    ENDIF
    IF ( USE_SRATRICK==0 ) then

      ALLOCATE (SLV(0:NCHEB,LMSIZE2,0:NCHEB,LMSIZE2),SRV(0:NCHEB,LMSIZE2,0:NCHEB,LMSIZE2) )

    ELSEIF ( USE_SRATRICK==1 ) then

      ALLOCATE (WORK2((NCHEB+1)*LMSIZE,(NCHEB+1)*LMSIZE), IPIV2((NCHEB+1)*LMSIZE))

      ALLOCATE (SLV1(0:NCHEB,LMSIZE,0:NCHEB,LMSIZE),SRV1(0:NCHEB,LMSIZE,0:NCHEB,LMSIZE), &
                SLV2(0:NCHEB,LMSIZE,0:NCHEB,LMSIZE),SRV2(0:NCHEB,LMSIZE,0:NCHEB,LMSIZE), &
                SLV3(0:NCHEB,LMSIZE,0:NCHEB,LMSIZE),SRV3(0:NCHEB,LMSIZE,0:NCHEB,LMSIZE) )

     ALLOCATE (      YILL1(0:NCHEB,LMSIZE,LMSIZE),ZILL1(0:NCHEB,LMSIZE,LMSIZE), &
                     YRLL1(0:NCHEB,LMSIZE,LMSIZE),ZRLL1(0:NCHEB,LMSIZE,LMSIZE), &
                     YILL2(0:NCHEB,LMSIZE,LMSIZE),ZILL2(0:NCHEB,LMSIZE,LMSIZE), &
                     YRLL2(0:NCHEB,LMSIZE,LMSIZE),ZRLL2(0:NCHEB,LMSIZE,LMSIZE),&
                     YRLLTMP(0:NCHEB,LMSIZE,LMSIZE)  )



   ENDIF

      ALLOCATE(  &
                WORK(LMSIZE,LMSIZE),&
                ALLP(LMSIZE,LMSIZE,0:NPAN),BLLP(LMSIZE,LMSIZE,0:NPAN),&
                CLLP(LMSIZE,LMSIZE,0:NPAN),DLLP(LMSIZE,LMSIZE,0:NPAN),&
                MRNVY(LMSIZE,LMSIZE,NPAN),MRNVZ(LMSIZE,LMSIZE,NPAN),&
                MRJVY(LMSIZE,LMSIZE,NPAN),MRJVZ(LMSIZE,LMSIZE,NPAN),&
                MIHVY(LMSIZE,LMSIZE,NPAN),MIHVZ(LMSIZE,LMSIZE,NPAN),&
                MIJVY(LMSIZE,LMSIZE,NPAN),MIJVZ(LMSIZE,LMSIZE,NPAN),&
                YILL(0:NCHEB,LMSIZE2,LMSIZE),ZILL(0:NCHEB,LMSIZE2,LMSIZE),&
                YRLL(0:NCHEB,LMSIZE2,LMSIZE),ZRLL(0:NCHEB,LMSIZE2,LMSIZE),&
                VJLR(LMSIZE,LMSIZE2,0:NCHEB),VHLR(LMSIZE,LMSIZE2,0:NCHEB),&
                VJLI(LMSIZE,LMSIZE2,0:NCHEB),VHLI(LMSIZE,LMSIZE2,0:NCHEB))
      YRLL=(0.0D0,0.0D0)
      ZILL=(0.0D0,0.0D0)
      YRLL=(0.0D0,0.0D0)
      ZILL=(0.0D0,0.0D0)

      ALLOCATE(&
                     YIF(LMSIZE2,LMSIZE,0:NCHEB,NPAN),&
                     YRF(LMSIZE2,LMSIZE,0:NCHEB,NPAN),&
                     ZIF(LMSIZE2,LMSIZE,0:NCHEB,NPAN),&
                     ZRF(LMSIZE2,LMSIZE,0:NCHEB,NPAN) )

      DO IPAN = 1,NPAN
        DO ICHEB = 0,NCHEB
          MN = IPAN*NCHEB + IPAN - ICHEB
          TAU(ICHEB,IPAN) = RMESH(MN)
        END DO
      END DO

      CALL GETCLAMBDACINV(NCHEB,CDDRC1)

      CALL CHEBINT(CSLC1,CSRC1,SLC1SUM,C1,NCHEB)

! loop over subintervals
      DO IPAN = 1,NPAN

! initialization

        VHLR=CZERO
        VJLR=CZERO
        VHLI=CZERO
        VJLI=CZERO

        IF (USE_SRATRICK==0) then

          YRLL=CZERO
          ZRLL=CZERO
          YILL=CZERO
          ZILL=CZERO
        ELSE
          YRLL1=CZERO
          ZRLL1=CZERO
          YILL1=CZERO
          ZILL1=CZERO
          YRLL2=CZERO
          ZRLL2=CZERO
          YILL2=CZERO
          ZILL2=CZERO
        ENDIF

!---------------------------------------------------------------------
! 1. prepare VJLR, VNL, VHLR, which appear in the integrands
! TAU(K,IPAN) is used instead of TAU(K,IPAN)**2, which directly gives
! RLL(r) and SLL(r) multiplied with r
!
! 2. prepare the source terms YR, ZR, YI, ZI
! because of the conventions used by
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! a factor sqrt(E) is included in the source terms
! this factor is removed by the definition of ZSLC1SUM given below
!


        DO ICHEB = 0,NCHEB
          MN = IPAN*NCHEB + IPAN - ICHEB

          IF     (CMODERLL=='1') THEN
            DO IVEC2=1,NVEC
              DO LM2 = 1,LMSIZE
                DO IVEC=1,NVEC
                  DO LM1 = 1,LMSIZE
                    L1 = JLK_INDEX( LM1+LMSIZE*(IVEC-1) )
                    VJLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VJLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*JLK2(L1,MN)*VLL(LM1+LMSIZE*(IVEC-1),LM2+LMSIZE*(IVEC2-1),MN)
                    VHLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VHLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*HLK2(L1,MN)*VLL(LM1+LMSIZE*(IVEC-1),LM2+LMSIZE*(IVEC2-1),MN)
                  END DO
                END DO
              END DO
            END DO !NVEC
          ELSEIF (CMODERLL=='T') THEN
            DO IVEC2=1,NVEC
              DO LM2 = 1,LMSIZE
                DO IVEC=1,NVEC
                  DO LM1 = 1,LMSIZE
                    L1 = JLK_INDEX( LM1+LMSIZE*(IVEC-1) )
                    VJLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VJLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*JLK2(L1,MN)*VLL(LM2+LMSIZE*(IVEC-1),LM1+LMSIZE*(IVEC2-1),MN)
                    VHLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VHLR(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*HLK2(L1,MN)*VLL(LM2+LMSIZE*(IVEC-1),LM1+LMSIZE*(IVEC2-1),MN)
                  END DO
                END DO
              END DO
            END DO !NVEC
          ELSEIF (CMODERLL=='0') THEN
                    VJLR(:,:,ICHEB) = CZERO
                    VHLR(:,:,ICHEB) = CZERO
          ELSE
            STOP'[RLLSLL] MODE NOT KNOWN'
          END IF
          IF     (CMODESLL=='1') THEN
            DO IVEC2=1,NVEC
              DO LM2 = 1,LMSIZE
                DO IVEC=1,NVEC
                  DO LM1 = 1,LMSIZE
                    L1 = JLK_INDEX( LM1+LMSIZE*(IVEC-1) )
                    VJLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VJLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*JLK2(L1,MN)*VLL(LM1+LMSIZE*(IVEC-1),LM2+LMSIZE*(IVEC2-1),MN)
                    VHLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VHLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*HLK2(L1,MN)*VLL(LM1+LMSIZE*(IVEC-1),LM2+LMSIZE*(IVEC2-1),MN)
                  END DO
                END DO
              END DO
            END DO !NVEC
          ELSEIF (CMODESLL=='T') THEN
            DO IVEC2=1,NVEC
              DO LM2 = 1,LMSIZE
                DO IVEC=1,NVEC
                  DO LM1 = 1,LMSIZE
                    L1 = jlk_index( LM1+LMSIZE*(IVEC-1) )
                    VJLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VJLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*JLK2(L1,MN)*VLL(LM2+LMSIZE*(IVEC-1),LM1+LMSIZE*(IVEC2-1),MN)
                    VHLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) = VHLI(LM1,LM2+LMSIZE*(IVEC2-1),ICHEB) + &
                        GMATPREFACTOR*TAU(ICHEB,IPAN)*HLK2(L1,MN)*VLL(LM2+LMSIZE*(IVEC-1),LM1+LMSIZE*(IVEC2-1),MN)
                  END DO
                END DO
              END DO
            END DO !NVEC
          ELSEIF (CMODESLL=='0') THEN
                    VJLI(:,:,ICHEB) = CZERO
                    VHLI(:,:,ICHEB) = CZERO
          ELSE
            STOP'[RLLSLL] MODE NOT KNOWN'
          END IF
    


          IF ( USE_SRATRICK==0 ) then
            DO IVEC=1,NVEC
              DO LM1 = 1,LMSIZE
                L1 = JLK_INDEX( LM1+LMSIZE*(IVEC-1) )
                YRLL(ICHEB,LM1+LMSIZE*(IVEC-1),LM1) =  TAU(ICHEB,IPAN)*JLK(L1,MN) 
                ZRLL(ICHEB,LM1+LMSIZE*(IVEC-1),LM1) =  TAU(ICHEB,IPAN)*HLK(L1,MN) 
                YILL(ICHEB,LM1+LMSIZE*(IVEC-1),LM1) =  TAU(ICHEB,IPAN)*HLK(L1,MN)
                ZILL(ICHEB,LM1+LMSIZE*(IVEC-1),LM1) =  TAU(ICHEB,IPAN)*JLK(L1,MN)
              END DO
            END DO !IVEC=1,NVEC

          ELSEIF ( USE_SRATRICK==1 ) then

            DO LM1 = 1,LMSIZE
              L1 = JLK_INDEX( LM1+LMSIZE*(1-1) )
              L2 = JLK_INDEX( LM1+LMSIZE*(2-1) )
              YRLL1(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*JLK(L1,MN) 
              ZRLL1(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*HLK(L1,MN) 
              YILL1(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*HLK(L1,MN)
              ZILL1(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*JLK(L1,MN)
              YRLL2(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*JLK(L2,MN) 
              ZRLL2(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*HLK(L2,MN) 
              YILL2(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*HLK(L2,MN)
              ZILL2(ICHEB,LM1+LMSIZE*(1-1),LM1) =  TAU(ICHEB,IPAN)*JLK(L2,MN)
            END DO
          ENDIF

        END DO ! ICHEB

! determine the matrices in equations (4.5a) and (4.5b)
        IF ( USE_SRATRICK==0 ) THEN
          DO ICHEB2 = 0,NCHEB
            DO ICHEB = 0,NCHEB
              MN = IPAN*NCHEB + IPAN - ICHEB
              DO LM2 = 1,LMSIZE2
                DO IVEC=1,NVEC
                  DO LM3 = 1,LMSIZE
                    LM1=LM3+(IVEC-1)*LMSIZE
                    L1 = JLK_INDEX(LM1)
                    SLV(ICHEB,LM1,ICHEB2,LM2) = &
                  ( TAU(ICHEB,IPAN)*JLK(L1,MN)*CSLC1(ICHEB,ICHEB2)*VHLR(LM3,LM2,ICHEB2) &
                    -TAU(ICHEB,IPAN)*HLK(L1,MN)*CSLC1(ICHEB,ICHEB2)*VJLR(LM3,LM2,ICHEB2))&
                  *(RPANBOUND(IPAN)-RPANBOUND(IPAN-1))/ 2.D0
                    SRV(ICHEB,LM1,ICHEB2,LM2) = &
                  (-TAU(ICHEB,IPAN)*JLK(L1,MN)*CSRC1(ICHEB,ICHEB2)*VHLI(LM3,LM2,ICHEB2) &
                    +TAU(ICHEB,IPAN)*HLK(L1,MN)*CSRC1(ICHEB,ICHEB2)*VJLI(LM3,LM2,ICHEB2)) &
                      *(RPANBOUND(IPAN)-RPANBOUND(IPAN-1))/ 2.D0
                  END DO
                END DO
              END DO
            END DO
          END DO
          DO LM1 = 1,LMSIZE2
            DO ICHEB = 0,NCHEB
              SLV(ICHEB,LM1,ICHEB,LM1) = SLV(ICHEB,LM1,ICHEB,LM1) + 1.D0
              SRV(ICHEB,LM1,ICHEB,LM1) = SRV(ICHEB,LM1,ICHEB,LM1) + 1.D0
            END DO
          END DO

        ELSEIF  ( USE_SRATRICK==1 ) then

          DO ICHEB2 = 0,NCHEB
            DO ICHEB = 0,NCHEB
              MN = IPAN*NCHEB + IPAN - ICHEB
              DO LM2 = 1,LMSIZE
                DO IVEC=1,1
                  DO LM3 = 1,LMSIZE
                    LM1=LM3+(IVEC-1)*LMSIZE
                    L1 = jlk_index(LM1)

                    SLV1(ICHEB,LM1,ICHEB2,LM2) = &
                  ( TAU(ICHEB,IPAN)*JLK(L1,MN)*CSLC1(ICHEB,ICHEB2)*VHLR(LM3,LM2,ICHEB2) &
                   -TAU(ICHEB,IPAN)*HLK(L1,MN)*CSLC1(ICHEB,ICHEB2)*VJLR(LM3,LM2,ICHEB2))&
                  *(RPANBOUND(IPAN)-RPANBOUND(IPAN-1))/ 2.D0

                    SRV1(ICHEB,LM1,ICHEB2,LM2) = &
                  (-TAU(ICHEB,IPAN)*JLK(L1,MN)*CSRC1(ICHEB,ICHEB2)*VHLI(LM3,LM2,ICHEB2) &
                   +TAU(ICHEB,IPAN)*HLK(L1,MN)*CSRC1(ICHEB,ICHEB2)*VJLI(LM3,LM2,ICHEB2)) &
                      *(RPANBOUND(IPAN)-RPANBOUND(IPAN-1))/ 2.D0

                  END DO
                END DO
              END DO
            END DO
          END DO

          DO ICHEB2 = 0,NCHEB
            DO ICHEB = 0,NCHEB
              MN = IPAN*NCHEB + IPAN - ICHEB
              DO LM2 = 1,LMSIZE
                DO IVEC=2,2
                  DO LM3 = 1,LMSIZE
                    LM1=LM3+(IVEC-1)*LMSIZE
                    L1 = JLK_INDEX(LM1)

                    SLV2(ICHEB,LM3,ICHEB2,LM2) = &
                  ( TAU(ICHEB,IPAN)*JLK(L1,MN)*CSLC1(ICHEB,ICHEB2)*VHLR(LM3,LM2,ICHEB2) &
                   -TAU(ICHEB,IPAN)*HLK(L1,MN)*CSLC1(ICHEB,ICHEB2)*VJLR(LM3,LM2,ICHEB2))&
                  *(RPANBOUND(IPAN)-RPANBOUND(IPAN-1))/ 2.D0

                    SRV2(ICHEB,LM3,ICHEB2,LM2) = &
                  (-TAU(ICHEB,IPAN)*JLK(L1,MN)*CSRC1(ICHEB,ICHEB2)*VHLI(LM3,LM2,ICHEB2) &
                   +TAU(ICHEB,IPAN)*HLK(L1,MN)*CSRC1(ICHEB,ICHEB2)*VJLI(LM3,LM2,ICHEB2)) &
                      *(RPANBOUND(IPAN)-RPANBOUND(IPAN-1))/ 2.D0

                  END DO
                END DO
              END DO
            END DO
          END DO


          DO LM1 = 1,LMSIZE
            DO ICHEB = 0,NCHEB
              SLV1(ICHEB,LM1,ICHEB,LM1) = SLV1(ICHEB,LM1,ICHEB,LM1) + 1.D0
              SRV1(ICHEB,LM1,ICHEB,LM1) = SRV1(ICHEB,LM1,ICHEB,LM1) + 1.D0
!               SLV2(ICHEB,LM1,ICHEB,LM1) = SLV2(ICHEB,LM1,ICHEB,LM1) + 1.D0
!               SRV2(ICHEB,LM1,ICHEB,LM1) = SRV2(ICHEB,LM1,ICHEB,LM1) + 1.D0
            END DO
          END DO
        ELSE
          stop '[rllsll] error in inversion'
        END IF

      IF ( USE_SRATRICK==0 ) then
        NPLM = (NCHEB+1)*LMSIZE2

          IF (CMODERLL/='0') THEN
            CALL ZGETRF(NPLM,NPLM,SLV,NPLM,IPIV,INFO)
            IF (INFO/=0) STOP'RLLSLL: ZGETRF'
            CALL ZGETRS('N',NPLM,LMSIZE,SLV,NPLM,IPIV,YRLL,NPLM,INFO)
            CALL ZGETRS('N',NPLM,LMSIZE,SLV,NPLM,IPIV,ZRLL,NPLM,INFO)
          END IF
          IF (CMODESLL/='0') THEN

            IF (DIRECTSOLV==1) THEN
             CALL ZGETRF(NPLM,NPLM,SRV,NPLM,IPIV,INFO)
             IF (INFO/=0) STOP'RLLSLL: ZGETRF'
             CALL ZGETRS('N',NPLM,LMSIZE,SRV,NPLM,IPIV,YILL,NPLM,INFO)
             CALL ZGETRS('N',NPLM,LMSIZE,SRV,NPLM,IPIV,ZILL,NPLM,INFO)
           ELSE 
             CALL ITERATIVESOL (NCHEB,LMSIZE2,LMSIZE,SRV,YILL)
             CALL ITERATIVESOL (NCHEB,LMSIZE2,LMSIZE,SRV,ZILL)
           ENDIF

         END IF

      ELSEIF ( USE_SRATRICK==1 ) THEN

        NPLM = (NCHEB+1)*LMSIZE

         CALL INVERSE(NPLM,SLV1)

         CALL INVERSE(NPLM,SRV1)


         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,CONE,SLV1, &
              NPLM,YRLL1,NPLM,CZERO,YRLLTMP,NPLM)
         YRLL1=YRLLTMP
         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,-CONE,SLV2, &
              NPLM,YRLL1,NPLM,CONE,YRLL2,NPLM)

         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,CONE,SLV1, &
              NPLM,ZRLL1,NPLM,CZERO,YRLLTMP,NPLM)
         ZRLL1=YRLLTMP
         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,-CONE,SLV2, &
              NPLM,ZRLL1,NPLM,CONE,ZRLL2,NPLM)

         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,CONE,SRV1, &
              NPLM,YILL1,NPLM,CZERO,YRLLTMP,NPLM)
         YILL1=YRLLTMP
         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,-CONE,SRV2, &
              NPLM,YILL1,NPLM,CONE,YILL2,NPLM)

         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,CONE,SRV1, &
              NPLM,ZILL1,NPLM,CZERO,YRLLTMP,NPLM)
         ZILL1=YRLLTMP
         CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,-CONE,SRV2, &
              NPLM,ZILL1,NPLM,CONE,ZILL2,NPLM)

        ELSE
          stop '[rllsll] error in inversion'
        END IF


      IF ( USE_SRATRICK==0 ) THEN

        DO ICHEB = 0,NCHEB
          DO LM2 = 1,LMSIZE
            DO LM1 = 1,LMSIZE2
              YRF(LM1,LM2,ICHEB,IPAN) = YRLL(ICHEB,LM1,LM2)
              ZRF(LM1,LM2,ICHEB,IPAN) = ZRLL(ICHEB,LM1,LM2)
              YIF(LM1,LM2,ICHEB,IPAN) = YILL(ICHEB,LM1,LM2)
              ZIF(LM1,LM2,ICHEB,IPAN) = ZILL(ICHEB,LM1,LM2)
            END DO
          END DO
        END DO

      ELSEIF ( USE_SRATRICK==1 ) then

        DO ICHEB = 0,NCHEB
          DO LM2 = 1,LMSIZE
            DO LM1 = 1,LMSIZE
              YRF(LM1,LM2,ICHEB,IPAN) = YRLL1(ICHEB,LM1,LM2)
              ZRF(LM1,LM2,ICHEB,IPAN) = ZRLL1(ICHEB,LM1,LM2)
              YIF(LM1,LM2,ICHEB,IPAN) = YILL1(ICHEB,LM1,LM2)
              ZIF(LM1,LM2,ICHEB,IPAN) = ZILL1(ICHEB,LM1,LM2)
            END DO
          END DO
        END DO

        DO ICHEB = 0,NCHEB
          DO LM2 = 1,LMSIZE
            DO LM1 = 1,LMSIZE
              YRF(LM1+LMSIZE,LM2,ICHEB,IPAN) = YRLL2(ICHEB,LM1,LM2)
              ZRF(LM1+LMSIZE,LM2,ICHEB,IPAN) = ZRLL2(ICHEB,LM1,LM2)
              YIF(LM1+LMSIZE,LM2,ICHEB,IPAN) = YILL2(ICHEB,LM1,LM2)
              ZIF(LM1+LMSIZE,LM2,ICHEB,IPAN) = ZILL2(ICHEB,LM1,LM2)
            END DO
          END DO
        END DO

      ENDIF

! (3.8), and (3.9) with the equations on page 143.

      DO ICHEB = 0,NCHEB
        ZSLC1SUM(ICHEB) = SLC1SUM(ICHEB) * (RPANBOUND(IPAN)-RPANBOUND(IPAN-1))/ (2.D0)
      END DO
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VHLR(1,1,0), &
              LMSIZE,YRF(1,1,0,IPAN),LMSIZE2,CZERO,MRNVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VJLR(1,1,0), &
              LMSIZE,YRF(1,1,0,IPAN),LMSIZE2,CZERO,MRJVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VHLR(1,1,0), &
              LMSIZE,ZRF(1,1,0,IPAN),LMSIZE2,CZERO,MRNVZ(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VJLR(1,1,0), &
              LMSIZE,ZRF(1,1,0,IPAN),LMSIZE2,CZERO,MRJVZ(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VHLI(1,1,0), &
              LMSIZE,YIF(1,1,0,IPAN),LMSIZE2,CZERO,MIHVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VJLI(1,1,0), &
              LMSIZE,YIF(1,1,0,IPAN),LMSIZE2,CZERO,MIJVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VHLI(1,1,0), &
              LMSIZE,ZIF(1,1,0,IPAN),LMSIZE2,CZERO,MIHVZ(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(0),VJLI(1,1,0), &
              LMSIZE,ZIF(1,1,0,IPAN),LMSIZE2,CZERO,MIJVZ(1,1,IPAN),LMSIZE)
        DO ICHEB = 1,NCHEB
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VHLR(1,1,ICHEB), &
              LMSIZE,YRF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MRNVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VJLR(1,1,ICHEB), &
              LMSIZE,YRF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MRJVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VHLR(1,1,ICHEB), &
              LMSIZE,ZRF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MRNVZ(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VJLR(1,1,ICHEB), &
              LMSIZE,ZRF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MRJVZ(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VHLI(1,1,ICHEB), &
              LMSIZE,YIF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MIHVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VJLI(1,1,ICHEB), &
              LMSIZE,YIF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MIJVY(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VHLI(1,1,ICHEB), &
              LMSIZE,ZIF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MIHVZ(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE2,ZSLC1SUM(ICHEB),VJLI(1,1,ICHEB), &
              LMSIZE,ZIF(1,1,ICHEB,IPAN),LMSIZE2,CONE,MIJVZ(1,1,IPAN),LMSIZE)
        END DO

      END DO ! IPAN

! end of loop over the subintervals

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculate A(M), B(M), C(M), D(M) for m from 1 to MMAX
! starting from A(0) = 1, B(0) = 0, C(MMAX) = 0 and D(MMAX) = 1

      DO LM2 = 1,LMSIZE
        DO LM1 = 1,LMSIZE
          BLLP(LM1,LM2,0) = 0.D0
          ALLP(LM1,LM2,0) = 0.D0
        END DO
      END DO
      DO LM1 = 1,LMSIZE
        ALLP(LM1,LM1,0) = 1.D0
      END DO
      DO IPAN = 1,NPAN
        CALL ZCOPY(LMSIZE*LMSIZE,ALLP(1,1,IPAN-1),1,ALLP(1,1,IPAN),1)
        CALL ZCOPY(LMSIZE*LMSIZE,BLLP(1,1,IPAN-1),1,BLLP(1,1,IPAN),1)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,-CONE,MRNVY(1,1,IPAN), &
                LMSIZE,ALLP(1,1,IPAN-1),LMSIZE,CONE,ALLP(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,-CONE,MRNVZ(1,1,IPAN), &
                LMSIZE,BLLP(1,1,IPAN-1),LMSIZE,CONE,ALLP(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE, CONE,MRJVY(1,1,IPAN), &
                LMSIZE,ALLP(1,1,IPAN-1),LMSIZE,CONE,BLLP(1,1,IPAN),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE, CONE,MRJVZ(1,1,IPAN), &
                LMSIZE,BLLP(1,1,IPAN-1),LMSIZE,CONE,BLLP(1,1,IPAN),LMSIZE)
      END DO

      DO LM2 = 1,LMSIZE
        DO LM1 = 1,LMSIZE
          DLLP(LM1,LM2,NPAN) = 0.D0
          CLLP(LM1,LM2,NPAN) = 0.D0
        END DO
      END DO

      DO LM1 = 1,LMSIZE
        DLLP(LM1,LM1,NPAN) = 1.D0
      END DO

      DO IPAN = NPAN,1,-1
        CALL ZCOPY(LMSIZE*LMSIZE,CLLP(1,1,IPAN),1,CLLP(1,1,IPAN-1),1)
        CALL ZCOPY(LMSIZE*LMSIZE,DLLP(1,1,IPAN),1,DLLP(1,1,IPAN-1),1)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE, CONE,MIHVZ(1,1,IPAN), &
                 LMSIZE,CLLP(1,1,IPAN),LMSIZE,CONE,CLLP(1,1,IPAN-1),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE, CONE,MIHVY(1,1,IPAN), &
                 LMSIZE,DLLP(1,1,IPAN),LMSIZE,CONE,CLLP(1,1,IPAN-1),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,-CONE,MIJVZ(1,1,IPAN), &
                 LMSIZE,CLLP(1,1,IPAN),LMSIZE,CONE,DLLP(1,1,IPAN-1),LMSIZE)
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,-CONE,MIJVY(1,1,IPAN), &
                 LMSIZE,DLLP(1,1,IPAN),LMSIZE,CONE,DLLP(1,1,IPAN-1),LMSIZE)
      END DO

!---------------------------------------------------------------------
! determine the regular solution ULL and the irregular solution SLL

      DO IPAN = 1,NPAN
        DO ICHEB = 0,NCHEB
          MN = IPAN*NCHEB + IPAN - ICHEB
        CALL ZGEMM('N','N',LMSIZE2,LMSIZE,LMSIZE,CONE,YRF(1,1,ICHEB,IPAN), &
                 LMSIZE2,ALLP(1,1,IPAN-1),LMSIZE,CZERO,ULL(1,1,MN),LMSIZE2)
        CALL ZGEMM('N','N',LMSIZE2,LMSIZE,LMSIZE,CONE,ZRF(1,1,ICHEB,IPAN), &
                 LMSIZE2,BLLP(1,1,IPAN-1),LMSIZE,CONE,ULL(1,1,MN),LMSIZE2)
        CALL ZGEMM('N','N',LMSIZE2,LMSIZE,LMSIZE,CONE,ZIF(1,1,ICHEB,IPAN), &
                 LMSIZE2,CLLP(1,1,IPAN),LMSIZE,CZERO,SLL(1,1,MN),LMSIZE2)
        CALL ZGEMM('N','N',LMSIZE2,LMSIZE,LMSIZE,CONE,YIF(1,1,ICHEB,IPAN), &
                 LMSIZE2,DLLP(1,1,IPAN),LMSIZE,CONE,SLL(1,1,MN),LMSIZE2)
        END DO
      END DO

!-----------------------------------------------------------------
! replace regular wave function in the first subinterval by a
! linear function times r**(l+1)
! replace irregular wave function in the first subinterval by a
! linear function divided by r**l
!       IF(config_testflag('wforigin')) THEN
!       DO LM2 =1,LMSIZE
!       DO LM1 =1,LMSIZE
!       PLLM = 0.5D0*(LOFLM(LM1)+LOFLM(LM2))
!       PLLM = LOFLM(LM1)
!       F1U = ULL(LM1,LM2,2*N+2-N)/TAU(N,2)**(PLLM+1.D0)
!       F2U = ULL(LM1,LM2,2*N+2-N/2)/TAU(N/2,2)**(PLLM+1.D0)
!       AU = (F1U*TAU(N/2,2)-F2U*TAU(N,2))/(TAU(N/2,2)-TAU(N,2))
!       BU = (F1U-F2U)/(TAU(N,2)-TAU(N/2,2))
!       F1S = SLL(LM1,LM2,2*N+2-N)*TAU(N,2)**PLLM
!       F2S = SLL(LM1,LM2,2*N+2-N/2)*TAU(N/2,2)**PLLM
!       AS = (F1S*TAU(N/2,2)-F2S*TAU(N,2))/(TAU(N/2,2)-TAU(N,2))
!       BS = (F1S-F2S)/(TAU(N,2)-TAU(N/2,2))
!       DO K = 0,N
!         MN = N + 1 - K
!         ULL(LM1,LM2,MN) = (AU+BU*TAU(K,1))*TAU(K,1)**(PLLM+1.D0)
!         SLL(LM1,LM2,MN) = (AS+BS*TAU(K,1))/TAU(K,1)**PLLM
!       END DO
!       END DO
!       END DO
!       END IF
! C-----------------------------------------------------------------
! transform from Volterra solution to Fredholm solution
! calculate alpha and t matrices
!       CALL ZAXPY(LMSIZE*LMSIZE,CI,BLLP(1,1,MMAX),1,ALLP(1,1,MMAX),1)       ! calculate the transformation matrix alpha
                                                                          ! assuming A is calculated with a neuman function
                                                                          ! n=h+ij (?)
!David
!       DO NM = 1,NRMAX
!       CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,CONE,SLL(1,1,NM),
!      +            LMSIZE,ALLP(1,1,MMAX),LMSIZE,CZERO,HLL(1,1,NM),LMSIZE)
!       END DO
! end David

      CALL ZGETRF(LMSIZE,LMSIZE,ALLP(1,1,NPAN),LMSIZE,IPIV,INFO)              !invert alpha
      CALL ZGETRI(LMSIZE,ALLP(1,1,NPAN),LMSIZE,IPIV,WORK,LMSIZE*LMSIZE,INFO)   !invert alpha -> transformation matrix RLL=alpha^-1*RLL
      CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,CONE/GMATPREFACTOR,BLLP(1,1,NPAN), &      ! calc t-matrix TLL = BLL*alpha^-1 
                  LMSIZE,ALLP(1,1,NPAN),LMSIZE,CZERO,TLLP,LMSIZE)

      DO NM = 1,NRMAX
      CALL ZGEMM('N','N',LMSIZE2,LMSIZE,LMSIZE,CONE,ULL(1,1,NM), &
                  LMSIZE2,ALLP(1,1,NPAN),LMSIZE,CZERO,RLL(1,1,NM),LMSIZE2)
      END DO
      DEALLOCATE(ULL)
      IF (USE_SRATRICK.EQ.0) THEN
       DEALLOCATE(SRV,SLV)
      ELSE
       DEALLOCATE(WORK2,IPIV2)
       DEALLOCATE(SLV1,SRV1,SLV2,SRV2,SLV3,SRV3)
       DEALLOCATE(YILL1,ZILL1,YRLL1,ZRLL1,YILL2,ZILL2,YRLL2,ZRLL2,YRLLTMP)
      ENDIF
      DEALLOCATE(WORK,ALLP,BLLP,CLLP,DLLP,MRNVY,MRNVZ,MRJVY,MRJVZ, &
                 MIHVY,MIHVZ,MIJVY,MIJVZ,YILL,ZILL,YRLL,ZRLL,&
                 VJLR,VHLR,VJLI,VHLI)
      DEALLOCATE(YIF,YRF,ZIF,ZRF)
     
      
!      RETURN
      END SUBROUTINE


subroutine inverse(nmat,mat)
!interface
integer        :: nmat
double complex :: mat(nmat,nmat)
double complex :: work(nmat,nmat)
!local
integer        :: IPIV(nmat)
integer        :: info

call ZGETRF( nmat, nmat, mat, nmat, IPIV, INFO )
if (info/=0) stop '[inverse] error INFO' 
call ZGETRI( nmat, mat, nmat, IPIV, WORK, nmat*nmat, INFO )
if (info/=0) stop '[inverse] error INFO' 
end subroutine inverse


subroutine iterativesol (NCHEB,LMSIZE2,LMSIZE,MMAT,BMAT)
integer :: NCHEB
integer :: LMSIZE,LMSIZE2
double complex :: MMAT(0:NCHEB,LMSIZE2,0:NCHEB,LMSIZE2)
double complex :: BMAT(0:NCHEB,LMSIZE2,LMSIZE)
double complex :: XMAT(0:NCHEB,LMSIZE2,LMSIZE)
!########################################################
! should solve the system of linear equations
! MMAT*XMAT = BMAT
!########################################################

NPLM = (NCHEB+1)*LMSIZE2
CALL ZGEMM('N','N',NPLM,LMSIZE,NPLM,CONE,SRV, &
    NPLM,ZILL,NPLM,CZERO,OUT,NPLM)



end subroutine iterativesol





