      SUBROUTINE ROTGLL(GMATLL,NATOMIMP,IJTABSYM,IJTABSH,
     +                  DSYMLL,SYMUNITARY,IGF,RC,CREL,RREL,
     +                  KREL,LMMAXD,IREC)
C **********************************************************************
C *                                                                    *
C *   it calculates all the elements of the Green Function of          *
C *   the impurity cluster using the GF calculated for the             *
C *   representative pairs.                                            *
C *   the representative pair and the symmetry operation D are given   *
C *   through the arrays IJTABSH and IJTABSYM set up in < SHELLGEN2K > *
C *                                                                    *
C *     _ _                                                            *
C *     n n'                      n n'                                 *
C *     m m'               T      m m'                                 *
C *    G    (E) = SUM    D     * G    (E) * D                          *
C *     L L'      L1 L2   L L1    L1L2      L2 L'                      *
C *                                                                    *
C *                   _                 _                              *
C *              n    n            n'   n'                             *
C *   where   D R  = R     and  D R  = R                               *
C *              m    m            m    m                              *
C *                                                                    *
C **********************************************************************
      use mod_mympi, only: myrank, master
      
      IMPLICIT NONE
C     ..
C     .. Parameter definitions
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO= (0.0D0,0.0D0),CONE= (1.D0,0.D0))
C     ..
C     .. Scalar arguments
      INTEGER NGCLUSD,NGCLUS,LMMAXD,IREC
      INTEGER IGF,KREL,NATOMIMP
C     ..
C     .. Array arguments
      INTEGER IJTABSYM(*),IJTABSH(*)
      DOUBLE COMPLEX GMATLL(LMMAXD,LMMAXD,*),
     &               DSYMLL(LMMAXD,LMMAXD,*)
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RC(LMMAXD,LMMAXD),
     &               RREL(LMMAXD,LMMAXD)
      LOGICAL SYMUNITARY(*)
C     ..
C     .. Local arrays
      DOUBLE COMPLEX GLL(:,:,:,:),TPG(:,:)
      COMPLEX*8 GCLUST(:)
      ALLOCATABLE GLL,TPG,GCLUST
CF77--------------------------------------------------------------------
Cccc      DOUBLE COMPLEX GLL(LMMAXD,LMMAXD,NATOMIMP,NATOMIMP),
Cccc     &               TPG(LMMAXD,LMMAXD)
Cccc      COMPLEX*8 GCLUST(NGCLUSD*NGCLUSD)
CF77--------------------------------------------------------------------
CF90--------------------------------------------------------------------
Cccc      DOUBLE COMPLEX GLL(:,:,:,:),TPG(:,:)
Cccc      COMPLEX*8 GCLUST(:)
Cccc      ALLOCATABLE GLL,TPG,GCLUST
CF90--------------------------------------------------------------------
C     ..
C     .. Local scalars
      INTEGER ILIN,IQ,ICALL,ISH,ISYM,JQ
      INTEGER LM1,LM2,NLIN
      CHARACTER*1 CNT
      CHARACTER*4 STR4I,STR4J
      CHARACTER*18 STR18
C     ..
C     .. External Subroutines
      EXTERNAL CHANGEREP,CMATSTR,ZGEMM,OPT
C     ..
C     .. External Functions 
      LOGICAL TEST,OPT
      EXTERNAL TEST
C     ..
C     .. Data statement
      DATA ICALL / 1 /
C     ..
C     .. Save statement
      SAVE ICALL
C     ..
CF90--------------------------------------------------------------------
      ALLOCATE (GLL(LMMAXD,LMMAXD,NATOMIMP,NATOMIMP),
     &          TPG(LMMAXD,LMMAXD),STAT=LM1)
      IF ( LM1.NE.0 ) THEN
         WRITE(6,99001) ' GLL/TPG'
         STOP '           < ROTGLL > '
      END IF
99001 FORMAT(6X,"ERROR: failed to allocate array(s) :",A,/)
CF90--------------------------------------------------------------------
C **********************************************************************
      IF ( ICALL.EQ.1) THEN
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
         if(myrank==master) WRITE (6,'(79(1H=))')
         if(myrank==master) WRITE (6,'(6X,2A)') 
     &     'ROTGLL : Expand GF for all pairs by rotation',
     &     ' and write out (all E-points)'      
         if(myrank==master) WRITE (6,'(79(1H=))')
         if(myrank==master) WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      END IF
C***********************************************************************
      DO IQ = 1,NATOMIMP
         DO JQ = 1,NATOMIMP
C-----------------------------------------------------------------------
            ILIN = (IQ-1)*NATOMIMP + JQ
            ISH = IJTABSH(ILIN)
            ISYM = IJTABSYM(ILIN)
C-----------------------------------------------------------------------
C    for REL CASE look if it is a unitary / ANTI - unitary rotation
C-----------------------------------------------------------------------
            CNT = 'N'
            IF ( .NOT.SYMUNITARY(ISYM) ) CNT = 'T'
C     
            CALL ZGEMM('C',CNT,LMMAXD,LMMAXD,LMMAXD,CONE,
     +           DSYMLL(1,1,ISYM),LMMAXD,GMATLL(1,1,ISH),
     +           LMMAXD,CZERO,TPG,LMMAXD)   
C                     
            CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,
     +           TPG,LMMAXD,DSYMLL(1,1,ISYM),
     +           LMMAXD,CZERO,GLL(1,1,IQ,JQ),LMMAXD)
C-----------------------------------------------------------------------
         END DO
      END DO
C***********************************************************************
C
C     visualise Gij
C
      IF ( TEST('Gmatij  ') ) THEN
         WRITE (6,'(/,4X,70(1H+),/,4X,A,I4)')
     &                     'cluster G_ij matrices for i,j = 1,',NATOMIMP
C
         DO IQ = 1,NATOMIMP
            WRITE(6,'(/,8X,66(1H=))')
            DO JQ = 1,NATOMIMP
C
               WRITE(STR4I,'(I4)') IQ
               WRITE(STR4J,'(I4)') JQ
               STR18 = '   i ='//STR4I(1:4)//' j ='//STR4J(1:4)
               IF (KREL.EQ.0) THEN 
                  CALL CMATSTR(STR18,18,GLL(1,1,IQ,JQ),LMMAXD,LMMAXD,
     &                         0,0,0,1d-8,6)
               ELSE
                  CALL CHANGEREP(GLL(1,1,IQ,JQ),'REL>RLM',TPG,LMMAXD,
     &                           LMMAXD,RC,CREL,RREL,STR18,18)
               END IF
               IF (JQ.LT.NATOMIMP) WRITE(6,'(/,9X,65(1H-))')
            END DO
            IF (IQ.EQ.NATOMIMP) WRITE(6,'(/,8X,66(1H=),/)')
         END DO
         WRITE(6,'(4X,70(1H+))')
      END IF
C     
C***********************************************************************
C
C --> output of GF
C
      NGCLUS=LMMAXD*NATOMIMP
      IF ( IGF.NE.0 ) THEN
         ICALL = ICALL + 1
         NLIN = 0
C
CF90--------------------------------------------------------------------
         ALLOCATE (GCLUST(NGCLUS*NGCLUS),STAT=LM1)
         IF ( LM1.NE.0 ) THEN
            WRITE(6,99001) ' GCLUST'
            STOP '           < ROTGLL > '
         END IF
CF90--------------------------------------------------------------------
         DO  JQ=1,NATOMIMP
            DO LM2=1,LMMAXD
               DO  IQ=1,NATOMIMP
                  DO LM1=1,LMMAXD
                     NLIN = NLIN + 1
                     IF ( NLIN.GT.NGCLUS*NGCLUS ) 
     &                    STOP "<ROTGLL>: NLIN.GT.(NATOMIMP*LMMAXD)**2"
                     GCLUST(NLIN) = GLL(LM1,LM2,IQ,JQ)
                  END DO
               END DO
            END DO
         END DO
C
!           write(*,*) 'icall ',icall
!          IF (ICALL==2) then
!           DO LM2=0,NATOMIMP*LMMAXD-1
!             WRITE(8888,'(50000F)') GCLUST(LM2*NATOMIMP*LMMAXD+1:
!      &                          (LM2+1)*NATOMIMP*LMMAXD)
!           END DO !LM2=0,NATOMIMP*LMMAXD-1
!          END IF !ICALL=1 then



         IF ( ( OPT('KKRFLEX ') ) ) THEN
#ifdef CPP_MPI
           IREC = IREC
#else
           IREC = ICALL
#endif
           WRITE(888,REC=IREC) GCLUST
           IF ( ( OPT('GPLAIN  ') ) ) THEN
             WRITE(8888,'(50000E)') GCLUST
           END IF
         END IF

         IF ( .not. OPT('KKRFLEX ') ) THEN
           WRITE(88,REC=ICALL) GCLUST
         END IF

CF90--------------------------------------------------------------------
         DEALLOCATE (GCLUST)
CF90--------------------------------------------------------------------
      ENDIF
C***********************************************************************
CF90--------------------------------------------------------------------
      DEALLOCATE (GLL,TPG)
CF90--------------------------------------------------------------------
      RETURN
      END                       ! SUBROUTINE ROTGLL
