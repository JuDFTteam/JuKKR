C***********************************************************************
C
      SUBROUTINE MMINVMOD(GLLH1,X2,TMATLL,NUMN0,INDN0,N2B,
     +                    IAT,SCITER,ITCOUNT,
     +                    GLLHBLCK,BCP,IGUESS,CNVFAC,
     +                    TOL)
C
C
      IMPLICIT NONE
C
C Parameters ..
      include 'inc.p'
      include 'inc.cls'
C      
      INTEGER            LMMAXD
      PARAMETER         (LMMAXD= (LMAXD+1)**2)
      INTEGER            NDIM,NAEZ
      PARAMETER         (NAEZ=NAEZD,NDIM=NAEZD*LMMAXD)  
      INTEGER            NGTBD
      PARAMETER         (NGTBD = NACLSD*LMMAXD)
      INTEGER            NBLCKD
      PARAMETER         (NBLCKD= XDIM * YDIM * ZDIM) 
      INTEGER            NLEN
      PARAMETER         (NLEN=NAEZD*LMMAXD)      
      DOUBLE COMPLEX     CONE, CZERO
      PARAMETER         (CONE = (1.0D0,0.0D0),CZERO = (0.0D0,0.0D0))
C ..
C
C global scalars ..
      DOUBLE PRECISION   CNVFAC,TOL
      INTEGER            IAT,SCITER,ITCOUNT,BCP,IGUESS
      LOGICAL            QMRABS
C ..
C
C global arrays ..
      DOUBLE COMPLEX     TMATLL(LMMAXD,LMMAXD,NAEZD),
     +                   X2(LMMAXD*NAEZD,LMMAXD),
     +                   GLLH1(LMMAXD,NGTBD,NAEZD),
     +                   GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD)
      INTEGER            NUMN0(NAEZD),
     +                   INDN0(NAEZD,NACLSD)
C ..
C
C local scalars ..
      DOUBLE COMPLEX     ZTMP
      DOUBLE PRECISION   DTMP,RESNAV,TOLAV
      INTEGER            NLIM,INIT,LM2,LM1,IL1,I,DIMVEC,IT,PROBE,
     +                   ITPROBE
      LOGICAL            EXITIT,CNVCONST
C ..
C
C local arrays ..
      DOUBLE COMPLEX     B(NDIM,LMMAXD),
     +                   VECS(NDIM,LMMAXD,9),
     +                   RHO(LMMAXD),ETA(LMMAXD),BETA(LMMAXD),
     +                   ALPHA(LMMAXD),
     +                   DUMMY(LMMAXD*NAEZD,LMMAXD)
      DOUBLE PRECISION   R0(LMMAXD),RESN(LMMAXD),VAR(LMMAXD),
     +                   TAU(LMMAXD),COS(LMMAXD),N2B(LMMAXD),
     +                   TOLB(LMMAXD),
     +                   HSTR(10)
      INTEGER            ITER(LMMAXD),HSTI(3)
      LOGICAL            DONE(LMMAXD)
C ..
C
C external ..
      EXTERNAL           DZNRM2,ZDOTU,ZRANDN
      DOUBLE COMPLEX     ZDOTU
      DOUBLE PRECISION   DZNRM2
C
C=======================================================================
C INITIALIZATION I
C=======================================================================
C
      QMRABS = .TRUE.           ! QMRABS tolerance for residual norm is defined globally
      CNVCONST = .TRUE.         ! CNVFAC is constant for all sc-steps
C
      IF (CNVCONST) THEN
        CNVFAC  = 1000
      ENDIF
C
      DO I=1,3
        HSTI(I) = 1             ! HSTI and HSTR store the history of QMR-convergency,
      ENDDO                     ! which is the basis to identify stagnation
      DO I=1,10
        HSTR(I) = 99+I
      ENDDO
C
      ITCOUNT = 0
      ITPROBE = 0
C
C=======================================================================
C INITIALIZATION II
C=======================================================================
C
      CALL CINIT(NDIM*LMMAXD,B)
C
      DO LM1=1,LMMAXD
        IL1=LMMAXD*(IAT-1)+LM1
        DO LM2=1,LMMAXD
          B(IL1,LM2)=-TMATLL(LM1,LM2,IAT)
        ENDDO
      ENDDO
C
      IF (IGUESS.EQ.1.AND.SCITER.GT.1) THEN
      DO I=1,NAEZD
      DO LM1=1,LMMAXD
        IL1=LMMAXD*(I-1)+LM1
        DO LM2=1,LMMAXD
          B(IL1,LM2)=-TMATLL(LM1,LM2,I)
        ENDDO
      ENDDO
      ENDDO
      ENDIF
C
      DO I = 1,NDIM
        DO LM2=1,LMMAXD
          VECS(I,LM2,2) = - B(I,LM2)
          VECS(I,LM2,1) = CZERO
          DO DIMVEC=3,9
            VECS(I,LM2,DIMVEC) = CZERO
          ENDDO
        ENDDO
      ENDDO
C


      NLIM = 2000
      INIT = 0
C
C--------------
      DO LM2=1,LMMAXD
C--------------
C
C        N2B(LM2)   = DZNRM2(NDIM,B(1,LM2),1)           ! norm of right-hand-sight
        IF (QMRABS) THEN
          TOLB(LM2)  = TOL
        ELSE
          TOLB(LM2)  = MAX(TOL*N2B(LM2),1.0D-10)          ! adapted relative tolerance
        ENDIF
C
        TOLAV      = TOLAV  + TOLB(LM2)/LMMAXD
C
        DONE(LM2) = .FALSE.
        ITER(LM2) = 0
        PROBE     = 1
C
        CALL ZAXPBY(NLEN,VECS(1,LM2,5),
     +              CONE,VECS(1,LM2,2),CZERO,VECS(1,LM2,5))
        CALL ZAXPBY(NLEN,VECS(1,LM2,1),
     +              CZERO,VECS(1,LM2,1),CZERO,VECS(1,LM2,1))
C
C--------------
      ENDDO
C--------------
C
C
C--------------
      DO LM2=1,LMMAXD
C--------------        
C        
C        R0(LM2) = DZNRM2(NLEN,VECS(1,LM2,5),1)
        R0(LM2) = N2B(LM2)
C
C     Check whether the auxiliary vector must be supplied.
C
        IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,LM2,3),1)
C
C     Initialize the variables.
C
        RESN(LM2) = 1.0D0
        RHO(LM2)  = CONE
        VAR(LM2)  = 0.0D0
        ETA(LM2)  = CZERO
        TAU(LM2)  = R0(LM2) * R0(LM2)
        CALL ZAXPBY(NLEN,VECS(1,LM2,8),
     +              CZERO,VECS(1,LM2,8),CZERO,VECS(1,LM2,8))
        CALL ZAXPBY(NLEN,VECS(1,LM2,4),
     +              CZERO,VECS(1,LM2,4),CZERO,VECS(1,LM2,4))
        CALL ZAXPBY(NLEN,VECS(1,LM2,6),
     +              CZERO,VECS(1,LM2,6),CZERO,VECS(1,LM2,6))      
C
C--------------
      ENDDO
C--------------
C 
C============================================================================
C============================================================================
C ITERATION
C
      DO IT=1, NLIM
C      
C============================================================================
C============================================================================
C
C--------------
        DO LM2=1,LMMAXD
        IF (.NOT.DONE(LM2)) THEN
C--------------
C
        ZTMP      = ZDOTU(NLEN,VECS(1,LM2,3),1,VECS(1,LM2,5),1)
        BETA(LM2) = ZTMP / RHO(LM2)
        RHO(LM2)  = ZTMP
        CALL ZAXPBY(NLEN,VECS(1,LM2,4),
     +              BETA(LM2),VECS(1,LM2,4),CONE,VECS(1,LM2,8)) 
        CALL ZAXPBY(NLEN,VECS(1,LM2,6),
     +              CONE,VECS(1,LM2,5),BETA(LM2),VECS(1,LM2,6))
C
C--------------
        ENDIF
        ENDDO
C--------------
C
        IF (BCP.EQ.1) THEN
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
C            WRITE(*,*) IT,LM2,DZNRM2(NDIM,DUMMY(1,LM2),1),
C     +                 BETA(LM2)
          ENDDO
C
          CALL APPBLCKCIRC(DUMMY,GLLHBLCK,
     &                     naez,lmaxd,nthrds,natbld,xdim,ydim,zdim)
C
        ELSE
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
          ENDDO
C
        ENDIF
C
        CALL SPRSZMM(
     >               IAT,GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE,
     >               CONE,CZERO,
     <               VECS(1,1,9))
C
C     VECS(:,:,6) input vector to be multiplied by A = GLLH1
C     VECS(:,:,9) result
C
C--------------
        DO LM2=1,LMMAXD
        IF (.NOT.DONE(LM2)) THEN
C--------------
C
        CALL ZAXPBY(NLEN,VECS(1,LM2,4),
     +              BETA(LM2),VECS(1,LM2,4),CONE,VECS(1,LM2,9))
C
        ZTMP = ZDOTU(NLEN,VECS(1,LM2,3),1,VECS(1,LM2,4),1)
C
        ALPHA(LM2) = RHO(LM2) / ZTMP

        ZTMP       = VAR(LM2) * ETA(LM2) / ALPHA(LM2)
        CALL ZAXPBY (NLEN,VECS(1,LM2,7),
     +               CONE,VECS(1,LM2,6),ZTMP,VECS(1,LM2,7))
        CALL ZAXPBY (NLEN,VECS(1,LM2,5),
     +               CONE,VECS(1,LM2,5),-ALPHA(LM2),VECS(1,LM2,9))
C
        DTMP = DZNRM2(NLEN,VECS(1,LM2,5),1)
        DTMP = DTMP * DTMP
        VAR(LM2)  = DTMP / TAU(LM2)
        COS(LM2)  = 1.0D0 / ( 1.0D0 + VAR(LM2) )
        TAU(LM2)  = DTMP * COS(LM2)
        ETA(LM2)  = ALPHA(LM2) * COS(LM2)
C
        CALL ZAXPBY(NLEN,VECS(1,LM2,1),
     +              CONE,VECS(1,LM2,1),ETA(LM2),VECS(1,LM2,7))
C
C--------------
        ENDIF
        ENDDO
C--------------
C      
C--------------
        DO LM2=1,LMMAXD
        IF (.NOT.DONE(LM2)) THEN
C--------------
C
        CALL ZAXPBY(NLEN,VECS(1,LM2,6),
     +              CONE,VECS(1,LM2,6),-ALPHA(LM2),VECS(1,LM2,4))
        ZTMP = VAR(LM2) * COS(LM2)
        CALL ZAXPBY(NLEN,VECS(1,LM2,7),
     +              CONE,VECS(1,LM2,6),ZTMP,VECS(1,LM2,7))
C
C--------------
        ENDIF
        ENDDO
C--------------
C
        EXITIT = .TRUE.
        DO LM2 = 1, LMMAXD
          IF (.NOT.DONE(LM2)) EXITIT = .FALSE.
        ENDDO
        IF (EXITIT) GO TO 67
C
C
        IF (BCP.EQ.1) THEN
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
          ENDDO
C
          CALL APPBLCKCIRC(DUMMY,GLLHBLCK,
     &                     naez,lmaxd,nthrds,natbld,xdim,ydim,zdim)
C
        ELSE
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
          ENDDO
C
        ENDIF
C
        CALL SPRSZMM(                               
     >               IAT,GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE,
     >               CONE,CZERO,                    
     <               VECS(1,1,8))                        
C
C     VECS(:,:,6) input vector to be multiplied by A = GLLH1
C     VECS(:,:,8) result
C
C--------------
        DO LM2=1,LMMAXD
        IF (.NOT.DONE(LM2)) THEN
C--------------
C
        CALL ZAXPBY(NLEN,VECS(1,LM2,5),
     +              CONE,VECS(1,LM2,5),-ALPHA(LM2),VECS(1,LM2,8))
C
        DTMP = DZNRM2(NLEN,VECS(1,LM2,5),1)
        DTMP = DTMP * DTMP
        VAR(LM2)  = DTMP / TAU(LM2)
        COS(LM2)  = 1.0D0 / ( 1.0D0 + VAR(LM2) )
        TAU(LM2)  = DTMP * COS(LM2)
        ETA(LM2)  = ALPHA(LM2) * COS(LM2)
C
        CALL ZAXPBY(NLEN,VECS(1,LM2,1),
     +              CONE,VECS(1,LM2,1),ETA(LM2),VECS(1,LM2,7))
C
C--------------
        ENDIF
        ENDDO
C--------------
C
C
C >>>>>>>>>>>
        IF (MOD(IT,PROBE).EQ.0) THEN
C
C in case of right-preconditioning
C                  -1
C         r = A * M  * y - b      otherwise   r = A * y - b
C                  2
C >>>>>>>>>>>
C has to be performed ..
C
        IF (BCP.EQ.1) THEN
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
          ENDDO
C
          CALL APPBLCKCIRC(DUMMY,GLLHBLCK,
     &                     naez,lmaxd,nthrds,natbld,xdim,ydim,zdim)
C
        ELSE
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
          ENDDO
C
        ENDIF
C
        CALL SPRSZMM(                               
     >               IAT,GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE,
     >               CONE,CZERO,                    
     <               VECS(1,1,9))
C
C     VECS(:,:,1) input vector to be multiplied by A = GLLH1
C     VECS(:,:,9) result
C
C--------------
        ITPROBE = ITPROBE + 1
        HSTI(MOD(ITPROBE,3)+1) = IT
C
        RESNAV = 0.0D0

        DO LM2=1,LMMAXD
        IF (.NOT.DONE(LM2)) THEN
C--------------
C
        CALL ZAXPBY(NLEN,VECS(1,LM2,9),
     +              CONE,VECS(1,LM2,2),-CONE,VECS(1,LM2,9))
        RESN(LM2) = DZNRM2(NLEN,VECS(1,LM2,9),1) / R0(LM2)
C
        IF (RESN(LM2).LE.TOLB(LM2)) THEN 
          DONE(LM2) = .TRUE.
          ITER(LM2) = 2*IT
        ENDIF
        RESNAV = RESNAV + RESN(LM2)/LMMAXD
C
C--------------
        ENDIF
        ENDDO
C--------------
C
        HSTR(MOD(ITPROBE,10)+1) = RESNAV
C
C check if QMR stagnated .. if yes leave cycle ..
        IF (MAX(HSTR(1),HSTR(2),HSTR(3),HSTR(4),HSTR(5),
     +          HSTR(6),HSTR(7),HSTR(8),HSTR(9),HSTR(10))
     +  .EQ.MIN(HSTR(1),HSTR(2),HSTR(3),HSTR(4),HSTR(5),
     +          HSTR(6),HSTR(7),HSTR(8),HSTR(9),HSTR(10))) THEN
          DO LM2=1,LMMAXD
          IF (.NOT.DONE(LM2)) THEN
            DONE(LM2) = .TRUE.
            ITER(LM2) = 2*IT
            WRITE(6,*) 'stgn.atm. ',IAT,LM2,' w. ',RESN(LM2),TOLB(LM2)
          ENDIF
          ENDDO
        ENDIF
C ..
C
        PROBE  = MAX(1,INT( LOG( RESNAV/(TOLAV*CNVFAC) ) ) )
C
        EXITIT = .TRUE.
        DO LM2 = 1, LMMAXD
          IF (.NOT.DONE(LM2)) EXITIT = .FALSE.
        ENDDO
        IF (EXITIT) GO TO 67
C
C <<<<<<<<<<<<
        ENDIF
C <<<<<<<<<<<<
C
C
C============================================================================
C============================================================================
      ENDDO
C
C ITERATION
C============================================================================
C============================================================================
C
   67 CONTINUE
C     Done.
C
C >>>>>>>>>>>
C
C in case of right-preconditioning
C              -1
C         x = M  * y
C              2
C has to be performed ..
C
        IF (BCP.EQ.1) THEN
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
          ENDDO
C
          CALL APPBLCKCIRC(DUMMY,GLLHBLCK,
     &                     naez,lmaxd,nthrds,natbld,xdim,ydim,zdim)
C
          DO LM2 = 1,LMMAXD
            CALL ZCOPY(NDIM,DUMMY(1,LM2),1,VECS(1,LM2,1),1)
          ENDDO
C
        ENDIF
C
C <<<<<<<<<<
C
      CALL CINIT(NDIM*LMMAXD,X2)
C
      DO LM2 = 1,LMMAXD
        CALL ZCOPY(NDIM,VECS(1,LM2,1),1,X2(1,LM2),1)
        ITCOUNT = ITCOUNT + ITER(LM2)
      ENDDO
C
C
      ITPROBE = MAX(HSTI(1),HSTI(2),HSTI(3))
     &        - MIN(HSTI(1),HSTI(2),HSTI(3))
      IF (CNVFAC.LT.1D+8.AND.CNVFAC.GT.1D-8) THEN
        CNVFAC = CNVFAC * 10.0**(MIN(ITPROBE,8)-5)
      ELSEIF (CNVFAC.GE.1D+8.AND.(MIN(ITPROBE,8)-5).LT.0) THEN
        CNVFAC = CNVFAC * 10.0**(MIN(ITPROBE,8)-5)
      ELSEIF (CNVFAC.LE.1D-8.AND.(MIN(ITPROBE,8)-5).GT.0) THEN
        CNVFAC = CNVFAC * 10.0**(MIN(ITPROBE,8)-5)
      ENDIF
C
      RETURN
      END
