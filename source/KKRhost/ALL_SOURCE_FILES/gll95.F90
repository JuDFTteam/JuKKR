module mod_gll13

contains

! **********************************************************************
SUBROUTINE gll13(ez,cleb,icleb,loflm,iend,tmatll,dtmatll,atom,  &
    refpot,ratom,natom,tolrdif,alat,out_wr,gref0, dgdeout,naclsmax,lly_g0tr,lly)
! **********************************************************************

!     solution of the DYSON equation for a cluster of potentials
!     (TMATLL) centered at positions RATOM in free space,

!     (modified version of GLL91 by P. Zahn, Sept. 95)
!     (modified by Phivos Mavropoulos to apply Lloyds formula
!      ported from KKRnano, Oct. 2013)
! ----------------------------------------------------------------------
#ifdef CPP_HYBRID
use omp_lib
#endif
use mod_types, only: t_inc
  use global_variables
   Use mod_datatypes, Only: dp
   use mod_gfree13
   use mod_grefsy13
      IMPLICIT NONE
!.. Parameters ..
      complex (kind=dp) CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
!..
!.. Scalar Arguments ..
      complex (kind=dp) EZ
      real (kind=dp) ALAT,TOLRDIF ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
      INTEGER IEND,NATOM,OUT_WR,NACLSMAX
!..
!.. Array Arguments ..
      complex (kind=dp) GREF0(NACLSMAX*LMGF0D,LMGF0D), &
                     TMATLL(LMGF0D,LMGF0D,nrefd)
      real (kind=dp) CLEB(ncleb),RATOM(3,naclsd)
      INTEGER ATOM(naclsd),ICLEB(NCLEB,4),LOFLM(lm2d),REFPOT(naezd+nembd)
!..
!.. Local Scalars ..
      INTEGER I,LM1,LM2,M,N,N1,N2,NDIM,NLM1,NLM2,INFO,NGD1
!..
!.. Local Arrays ..
      INTEGER IPVT(:)
      real (kind=dp) RDIFF(3),ABSRDIFF
      complex (kind=dp) DTMATLL(LMGF0D,LMGF0D,nrefd) ! Derivative of ref.-sys t-matrix
      complex (kind=dp) GREF(:,:),GLL(:,:),GTREF(:,:)
      complex (kind=dp) DGTDE(:,:),DGTDE0(:,:) ! LLY (1-gt)^-1 * d(1-gt)/dE (after grefsy13)
      complex (kind=dp) DGLLDE(:,:),DGDE(:,:), &
                     DGDEOUT(NACLSMAX*LMGF0D,LMGF0D)
      complex (kind=dp) LLY_G0TR   ! LLY Trace of  DTGLL for Lloyds formula
      INTEGER LLY     ! LLY =0 : no Lloyd's formula; <>0: use Lloyd's formula
      ALLOCATABLE GREF,GLL,GTREF,DGLLDE,DGTDE,DGTDE0,DGDE,IPVT
!..
!.. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DBLE
#ifdef CPP_HYBRID
      INTEGER thread_id
#endif

!     ..
ngd1 = naclsmax*lmgf0d
ndim = lmgf0d*natom

allocate( gtref(ngd1,lmgf0d),dgtde(ngd1,lmgf0d), stat=lm1)
IF(lm1/=0) STOP 'Error allocating gtref etc. <GLL95>'
gtref(:,:) = czero
dgtde(:,:) = czero

allocate(gref(ngd1,ngd1),ipvt(ngd1),stat=lm1)
IF(lm1/=0) STOP 'Error allocating gref etc. <GLL95>'
gref(:,:) = czero
ipvt(:) = 0
IF (lly /= 0) THEN
  allocate ( dgtde0(ngd1,ngd1) , dgde(ngd1,ngd1), stat=lm1 )
  IF(lm1/=0) STOP 'Error allocating dgtde0 etc. <GLL95>'
  dgtde0(:,:) = czero
  dgtde(:,:) = czero
endif
99001 FORMAT(6X,"ERROR: failed to allocate array(s) :",a,/)
IF (test('flow    ').AND.(t_inc%i_write>0)) WRITE (1337,FMT=*) '>>> GLL95'


#ifdef CPP_HYBRID
!$omp parallel default(shared) &
!$omp& private(n1, n2, rdiff, absrdiff, lm2, lm1, nlm2, nlm1, GLL) &
!$omp& private(thread_id, gtref, dgtde, DGLLDE)
thread_id = omp_get_thread_num()
#endif
! allocate here, inside omp parallel region
allocate(gll(lmgf0d,lmgf0d),dgllde(lmgf0d,lmgf0d), stat=lm1)
IF(lm1/=0) STOP 'Error allocating gll etc. <GLL95>'
gll(:,:) = czero
dgllde(:,:) = czero


! ---> construct free Green's function

#ifdef CPP_HYBRID
!$omp do
#endif
DO n1 = 1,natom
  DO n2 = 1,natom
    rdiff(1:3) = - (ratom(1:3,n1)-ratom(1:3,n2))*alat
    absrdiff=SQRT(rdiff(1)**2+rdiff(2)**2+rdiff(3)**2)
    
    IF (n1 /= n2 .AND. (absrdiff > tolrdif) ) THEN
      CALL gfree13(rdiff,ez,gll,dgllde,cleb,icleb,loflm,iend)
      DO lm2 = 1,lmgf0d
        nlm2 = (n2-1)*lmgf0d + lm2
        DO lm1 = 1,lmgf0d
          nlm1 = (n1-1)*lmgf0d + lm1
          gref(nlm1,nlm2) = gll(lm1,lm2)
          IF (lly /= 0) dgde(nlm1,nlm2) = dgllde(lm1,lm2)
        END DO
      END DO
    ELSE
      DO lm2 = 1,lmgf0d
        nlm2 = (n2-1)*lmgf0d + lm2
        DO lm1 = 1,lmgf0d
          nlm1 = (n1-1)*lmgf0d + lm1
          gref(nlm1,nlm2) = czero
          IF (lly /= 0) dgde(nlm1,nlm2) = czero
        END DO
      END DO
    endif
    
  END DO
END DO
#ifdef CPP_HYBRID
!$omp end do
! deallocate in omp parallel region
deallocate( gll, dgllde, stat=lm1)
IF ( lm1 /= 0 ) STOP ' [gll13] dealloc'
!$omp end parallel
#endif
IF (test('flow    ').AND.(t_inc%i_write>0)) WRITE (1337,FMT=*) 'GFREE o.k.'
! ----------------------------------------------------------------------

! GREF0 = g:= gfree
CALL zcopy(ngd1*lmgf0d,gref,1,gref0,1)

! ----------------------------------------------------------------------
! LLY Lloyd
! Prepare source term -dg/dE * t - g * dt/dE
IF (lly /= 0) THEN
  
  DO n2 = 1,natom
    nlm2 = (n2-1)*lmgf0d + 1
! GTREF = -DGDE*t = -dg/dE * t
    CALL zgemm('N','N',ndim,lmgf0d,lmgf0d,-cone,dgde(1,nlm2),  &
        ngd1,tmatll(1,1,refpot(ABS(atom(n2)))),lmgf0d, czero,gtref,ngd1)
! GTREF = GTREF - GREF*DTMATLL = -dg/dE * t - g * dt/dE  (here, GREF=g:=gfree)
    CALL zgemm('N','N',ndim,lmgf0d,lmgf0d,-cone,gref(1,nlm2),  &
        ngd1,dtmatll(1,1,refpot(ABS(atom(n2)))),lmgf0d, cone,gtref,ngd1)
    CALL zcopy(ngd1*lmgf0d,gtref,1,dgtde0(1,nlm2),1)
  END DO
  DO n2 = 1,lmgf0d
    DO n1 = 1,ngd1
      dgtde(n1,n2) = dgtde0(n1,n2)
    END DO
  END DO
! Now DGTDE = DGTDE0 = -dg/dE * t - g * dt/dE
! (DGTDE is reduced matrix; DGTDE0 is full matrix)
  
endif ! (LLY.NE.0)
! LLY Lloyd
! ----------------------------------------------------------------------
DO n2 = 1,natom
  nlm2 = (n2-1)*lmgf0d + 1
! GTREF =  -g*t
  CALL zgemm('N','N',ndim,lmgf0d,lmgf0d,-cone,gref(1,nlm2),  &
      ngd1,tmatll(1,1,refpot(ABS(atom(n2)))), lmgf0d,czero,gtref,ngd1)
  CALL zcopy(ngd1*lmgf0d,gtref,1,gref(1,nlm2),1)
! Now GREF =  -g*t
  IF (test('REFPOT  ').AND.(t_inc%i_write>0)) WRITE (1337,FMT=*)  &
      n2,refpot(ABS(atom(n2)))
END DO

IF (test('WAIT    ')) WRITE (6,FMT=*) 'Input I'
IF (test('WAIT    ')) READ (5,FMT=*) i

CALL grefsy13(gref,gref0,dgtde,lly_g0tr,ipvt, ndim,lly,lmgf0d,ngd1)
! Now GREF contains LU(1-gt) (full matrix NGD1xNGD1)
! DGTDE contains (1-gt)^-1 * d(1-gt)/dE (Thiess PhD Eq.5.28)
! between atoms 0 and n (n running through the cluster)
! and LLY_G0TR contains -Trace[ (1-gt)^-1 * d(1-gt)/dE ] (Thiess PhD Eq.5.38)


! ----------------------------------------------------------------------
! LLY Lloyd
dgdeout(:,:) = czero
IF (lly /= 0) THEN
  
! Prepare dg/de + (dg/dE * t + g * dt/dE)*Gref (Thiess PhD Eq.5.42)
  
! DGDE = DGDE - DGTDE0*GREF0 = dg/de + (dg/dE * t + g * dt/dE)*Gref
  CALL zgemm('N','N',ndim,lmgf0d,ndim,-cone,dgtde0,ngd1,  &
      gref0,ngd1,cone,dgde,ngd1)
  
! Solve linear system: (remember GREF contains LU(1-gt))
! (1-gt)*DGDE = dg/de + (dg/dE * t + g * dt/dE)*Gref
  CALL zgetrs('N',ndim,lmgf0d,gref,ngd1,ipvt,dgde,ngd1,info)
! Result is DGDE = dGref/dE
  
  DO n2 = 1,lmgf0d
    DO n1 = 1,ngd1
      dgdeout(n1,n2) = dgde(n1,n2)
    END DO
  END DO
  
endif
! LLY Lloyd
! ----------------------------------------------------------------------


IF (test('flow    ').AND.(t_inc%i_write>0)) WRITE (1337,FMT=*) 'GREFSY o.k.'

IF (out_wr > 0) WRITE (out_wr) ((gref0(n,m),m=1,lmgf0d),n=1,ngd1)

! deallocate arrays
deallocate( gtref, dgtde, gref, ipvt, stat=lm1)
IF ( lm1 /= 0 ) STOP ' [gll13] dealloc'
IF (lly /= 0) THEN
  deallocate ( dgtde0 , dgde, stat=lm1 )
  IF ( lm1 /= 0 ) STOP ' [gll13] dealloc'
endif
END SUBROUTINE gll13

end module mod_gll13
