c **********************************************************************
      SUBROUTINE GLL13(EZ,CLEB,ICLEB,LOFLM,IEND,TMATLL,DTMATLL,ATOM,
     +                 REFPOT,RATOM,NATOM,TOLRDIF,ALAT,OUT_WR,GREF0,
     &                 DGDEOUT,NACLSMAX,LLY_G0TR,LLY)
c **********************************************************************
c
c     solution of the DYSON equation for a cluster of potentials
c     (TMATLL) centered at positions RATOM in free space,
c
c     (modified version of GLL91 by P. Zahn, Sept. 95)
!     (modified by Phivos Mavropoulos to apply Lloyds formula
!      ported from KKRnano, Oct. 2013)
c ----------------------------------------------------------------------
#ifdef CPP_hybrid
      use omp_lib
#endif
      use mod_types, only: t_inc
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER LMAX,NATOMD
      PARAMETER (LMAX=LMAXD,NATOMD=NACLSD)
      INTEGER LMGF0D,NGD
      PARAMETER (LMGF0D= (LMAX+1)**2,NGD=LMGF0D*NATOMD)
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EZ
      DOUBLE PRECISION ALAT,TOLRDIF ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
      INTEGER IEND,NATOM,OUT_WR,NACLSMAX
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX GREF0(NACLSMAX*LMGF0D,LMGF0D),
     &               TMATLL(LMGF0D,LMGF0D,*)
      DOUBLE PRECISION CLEB(*),RATOM(3,*)
      INTEGER ATOM(*),ICLEB(NCLEB,4),LOFLM(*),REFPOT(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,LM1,LM2,M,N,N1,N2,NDIM,NLM1,NLM2,INFO,NGD1
C     ..
C     .. Local Arrays ..
      INTEGER IPVT(:)
      DOUBLE PRECISION RDIFF(3),ABSRDIFF
      DOUBLE COMPLEX DTMATLL(LMGF0D,LMGF0D,*) ! Derivative of ref.-sys t-matrix
      DOUBLE COMPLEX GREF(:,:),GLL(:,:),GTREF(:,:)
      DOUBLE COMPLEX DGTDE(:,:),DGTDE0(:,:) ! LLY (1-gt)^-1 * d(1-gt)/dE (after grefsy13)
      DOUBLE COMPLEX DGLLDE(:,:),DGDE(:,:),
     &               DGDEOUT(NACLSMAX*LMGF0D,LMGF0D)
      DOUBLE COMPLEX LLY_G0TR   ! LLY Trace of  DTGLL for Lloyds formula
      INTEGER LLY     ! LLY =0 : no Lloyd's formula; <>0: use Lloyd's formula
      ALLOCATABLE GREF,GLL,GTREF,DGLLDE,DGTDE,DGTDE0,DGDE,IPVT
C     ..
C     .. External Subroutines ..
      EXTERNAL GFREE13,GREFSY13,ZCOPY,ZGEMM
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE
#ifdef CPP_hybrid
      INTEGER thread_id
#endif
C     ..
      NGD1 = NACLSMAX*LMGF0D
      NDIM = LMGF0D*NATOM

      ALLOCATE( GTREF(NGD1,LMGF0D),DGTDE(NGD1,LMGF0D), stat=lm1)
      IF(LM1/=0) stop 'Error allocating gtref etc. <GLL95>'
      GTREF(:,:) = CZERO
      DGTDE(:,:) = CZERO
      
      ALLOCATE(GREF(NGD1,NGD1),IPVT(NGD1),STAT=LM1)
      IF(LM1/=0) stop 'Error allocating gref etc. <GLL95>'
      GREF(:,:) = CZERO
      IPVT(:) = 0
      IF (LLY.NE.0) THEN
         ALLOCATE ( DGTDE0(NGD1,NGD1) , DGDE(NGD1,NGD1), STAT=LM1 )
         IF(LM1/=0) stop 'Error allocating dgtde0 etc. <GLL95>'
         DGTDE0(:,:) = CZERO
         DGTDE(:,:) = CZERO
      ENDIF
99001 FORMAT(6X,"ERROR: failed to allocate array(s) :",A,/)
      IF (TEST('flow    ').and.(t_inc%i_write>0)) 
     &    WRITE (1337,FMT=*) '>>> GLL95'


#ifdef CPP_hybrid
!$omp parallel default(shared)
!$omp& private(n1, n2, rdiff, absrdiff, lm2, lm1, nlm2, nlm1, GLL)
!$omp& private(thread_id, gtref, dgtde, DGLLDE)
      thread_id = omp_get_thread_num()
#endif
      ! allocate here, inside omp parallel region
      ALLOCATE(GLL(LMGF0D,LMGF0D),DGLLDE(LMGF0D,LMGF0D), stat=lm1)
      IF(LM1/=0) stop 'Error allocating gll etc. <GLL95>'
      GLL(:,:) = CZERO
      DGLLDE(:,:) = CZERO

c
c ---> construct free Green's function
c
# ifdef CPP_hybrid
!$omp do
#endif
      DO N1 = 1,NATOM
        DO N2 = 1,NATOM
           RDIFF(1:3) = - (RATOM(1:3,N1)-RATOM(1:3,N2))*ALAT
           ABSRDIFF=SQRT(RDIFF(1)**2+RDIFF(2)**2+RDIFF(3)**2)
C
           IF (N1.NE.N2 .and. (ABSRDIFF.GT.TOLRDIF) ) THEN
              CALL GFREE13(RDIFF,EZ,GLL,DGLLDE,CLEB,ICLEB,LOFLM,IEND)
              DO LM2 = 1,LMGF0D
                 NLM2 = (N2-1)*LMGF0D + LM2
                 DO LM1 = 1,LMGF0D
                    NLM1 = (N1-1)*LMGF0D + LM1
                    GREF(NLM1,NLM2) = GLL(LM1,LM2)
                    IF (LLY.NE.0) DGDE(NLM1,NLM2) = DGLLDE(LM1,LM2)
                 END DO
              END DO
           ELSE
              DO LM2 = 1,LMGF0D
                 NLM2 = (N2-1)*LMGF0D + LM2
                 DO LM1 = 1,LMGF0D
                    NLM1 = (N1-1)*LMGF0D + LM1
                    GREF(NLM1,NLM2) = CZERO
                    IF (LLY.NE.0) DGDE(NLM1,NLM2) = CZERO
                 END DO
              END DO
           END IF

        END DO
      END DO
#ifdef CPP_hybrid
!$omp end do
      ! deallocate in omp parallel region
      DEALLOCATE( GLL, DGLLDE, STAT=LM1)
      IF ( LM1.NE.0 ) STOP ' [gll13] dealloc'
!$omp end parallel
#endif
      IF (TEST('flow    ').and.(t_inc%i_write>0)) 
     &      WRITE (1337,FMT=*) 'GFREE o.k.'
c ----------------------------------------------------------------------

      ! GREF0 = g:= gfree
      CALL ZCOPY(NGD1*LMGF0D,GREF,1,GREF0,1)

! ----------------------------------------------------------------------
! LLY Lloyd 
! Prepare source term -dg/dE * t - g * dt/dE 
      IF (LLY.NE.0) THEN

         DO N2 = 1,NATOM
            NLM2 = (N2-1)*LMGF0D + 1
            ! GTREF = -DGDE*t = -dg/dE * t
            CALL ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,DGDE(1,NLM2),
     +             NGD1,TMATLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D,
     +             CZERO,GTREF,NGD1)
            ! GTREF = GTREF - GREF*DTMATLL = -dg/dE * t - g * dt/dE  (here, GREF=g:=gfree)
            CALL ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,GREF(1,NLM2),
     +             NGD1,DTMATLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D,
     +             CONE,GTREF,NGD1)
            CALL ZCOPY(NGD1*LMGF0D,GTREF,1,DGTDE0(1,NLM2),1)
         END DO
         DO N2 = 1,LMGF0D
            DO N1 = 1,NGD1 
               DGTDE(N1,N2) = DGTDE0(N1,N2) 
            ENDDO
         ENDDO
         ! Now DGTDE = DGTDE0 = -dg/dE * t - g * dt/dE 
         ! (DGTDE is reduced matrix; DGTDE0 is full matrix)

      END IF ! (LLY.NE.0)
! LLY Lloyd
! ----------------------------------------------------------------------
      DO N2 = 1,NATOM
         NLM2 = (N2-1)*LMGF0D + 1
         ! GTREF =  -g*t
         CALL ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,GREF(1,NLM2),
     +              NGD1,TMATLL(1,1,REFPOT(ABS(ATOM(N2)))),
     +              LMGF0D,CZERO,GTREF,NGD1)
         CALL ZCOPY(NGD1*LMGF0D,GTREF,1,GREF(1,NLM2),1)
         ! Now GREF =  -g*t
         IF (TEST('REFPOT  ').and.(t_inc%i_write>0)) WRITE (1337,FMT=*)
     &                         N2,REFPOT(ABS(ATOM(N2)))
      END DO

      IF (TEST('WAIT    ')) WRITE (6,FMT=*) 'Input I'
      IF (TEST('WAIT    ')) READ (5,FMT=*) I

      CALL GREFSY13(GREF,GREF0,DGTDE,LLY_G0TR,IPVT,
     &              NDIM,LLY,LMGF0D,NGD1)
      ! Now GREF contains LU(1-gt) (full matrix NGD1xNGD1)
      ! DGTDE contains (1-gt)^-1 * d(1-gt)/dE (Thiess PhD Eq.5.28)
      ! between atoms 0 and n (n running through the cluster)
      ! and LLY_G0TR contains -Trace[ (1-gt)^-1 * d(1-gt)/dE ] (Thiess PhD Eq.5.38)


! ----------------------------------------------------------------------
! LLY Lloyd 
      DGDEOUT(:,:) = CZERO
      IF (LLY.NE.0) THEN

         ! Prepare dg/de + (dg/dE * t + g * dt/dE)*Gref (Thiess PhD Eq.5.42)

         ! DGDE = DGDE - DGTDE0*GREF0 = dg/de + (dg/dE * t + g * dt/dE)*Gref
         CALL ZGEMM('N','N',NDIM,LMGF0D,NDIM,-CONE,DGTDE0,NGD1,
     +           GREF0,NGD1,CONE,DGDE,NGD1)

         ! Solve linear system: (remember GREF contains LU(1-gt))
         ! (1-gt)*DGDE = dg/de + (dg/dE * t + g * dt/dE)*Gref
         CALL ZGETRS('N',NDIM,LMGF0D,GREF,NGD1,IPVT,DGDE,NGD1,INFO)
         ! Result is DGDE = dGref/dE

         DO N2 = 1,LMGF0D
            DO N1 = 1,NGD1 
               DGDEOUT(N1,N2) = DGDE(N1,N2)
            ENDDO
         ENDDO

      ENDIF
! LLY Lloyd
! ----------------------------------------------------------------------


      IF (TEST('flow    ').and.(t_inc%i_write>0)) 
     &      WRITE (1337,FMT=*) 'GREFSY o.k.'

      IF (OUT_WR.GT.0) WRITE (OUT_WR) ((GREF0(N,M),M=1,LMGF0D),N=1,NGD1)

      ! deallocate arrays
      DEALLOCATE( GTREF, DGTDE, GREF, IPVT, STAT=LM1)
      IF ( LM1.NE.0 ) STOP ' [gll13] dealloc'
      IF (LLY.NE.0) THEN
         DEALLOCATE ( DGTDE0 , DGDE, STAT=LM1 )
         IF ( LM1.NE.0 ) STOP ' [gll13] dealloc'
      END IF
      END
