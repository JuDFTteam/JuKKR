
SUBROUTINE calctmat(icst,ins,ielast, nsra,ispin,nspin,i1,ez,  &
    drdi,rmesh,vins,visp,zat,irmin,ipan, &               ! Added IRMIN 1.7.2014  &
    ircut,cleb,loflm,icleb,iend,solver,soctl,ctl,  &
    vtrel,btrel,rmrel,drdirel,r2drdirel, zrel,jwsrel,idoldau,lopt,wldau,  &
    lly,deltae) ! LLY

! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *  LDA+U implementation     Mar. 2002-Dec.2004                      *
! *                           ph.mavropoulos, h. ebert, v. popescu    *
! * Notes:                                                            *
! *  average WLDAU for spherical wavefunctions:                       *
! *  The spherical part of the d or f wavefunction is found by adding *
! *  the average interaction potential WLDAUAV to the spherical       *
! *  potential. Then the non-spherical parts are found by using only  *
! *  the deviation of WLDAU from the average. This speeds up the      *
! *  convergence of the Born series. See also subroutines             *
! *  regsol, pnstmat and pnsqns                                       *
! *                                                                   *
! *********************************************************************
#ifdef CPP_MPI
use mpi
#endif
use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
use mod_mympi, only: distribute_linear_on_tasks, mpiadapt
use mod_timing
#endif
use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid,init_tgmat,  &
    t_lloyd,init_tlloyd
use mod_DataTypes

IMPLICIT NONE

!.. Parameters ..
INCLUDE 'inc.p'
INTEGER LMMAXD
PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
INTEGER LMPOTD
PARAMETER (LMPOTD= (LPOTD+1)**2)
INTEGER IRMIND
PARAMETER (IRMIND=IRMD-IRNSD)
INTEGER MMAXD
PARAMETER ( MMAXD = 2*LMAXD+1 )
INTEGER LM2D
PARAMETER (LM2D= (2*LMAXD+1)**2)
real (kind=dp) CVLIGHT
PARAMETER (CVLIGHT=274.0720442D0)
complex (kind=dp) CZERO,CONE
PARAMETER (CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))
!..
!.. Scalar Arguments ..
real (kind=dp) ZAT
INTEGER I1,ICST,IELAST,IEND,INS,IPAN,ISPIN,NSPIN,NSRA
INTEGER IDOLDAU,JWSREL,LOPT,ZREL
INTEGER LLY ! LLY <> 0 for applying Lloyds formula
INTEGER SIGNDE,IDERIV
complex (kind=dp) DELTAE  ! Energy difference for numerical derivative
!..
!.. Array Arguments ..
complex (kind=dp) EZ(IEMXD), &
               TMAT0(LMMAXD,LMMAXD), &
               TMATLL(LMMAXD,LMMAXD), &
               DTMATLL(LMMAXD,LMMAXD) ! LLY t-matrix derivative
real (kind=dp) CLEB(NCLEB,2),DRDI(IRMD),RMESH(IRMD), &
                 VINS(IRMIND:IRMD,LMPOTD), &
                 VISP(IRMD)
CHARACTER*10 SOLVER
real (kind=dp) SOCTL(KREL*LMAXD+1)
real (kind=dp) CTL(KREL*LMAXD+1)
INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAND), &
        LOFLM(LM2D)
real (kind=dp) VTREL(IRMD*KREL+(1-KREL))
real (kind=dp) BTREL(IRMD*KREL+(1-KREL))
real (kind=dp) DRDIREL(IRMD*KREL+(1-KREL)), &
                 R2DRDIREL(IRMD*KREL+(1-KREL)), &
                 RMREL(IRMD*KREL+(1-KREL))
!..
!.. Local Scalars ..
complex (kind=dp) ERYD,EK
INTEGER IE,IREC,LM1,LM2,LMHI,LMLO,M1,MMAX,IRMIN
real (kind=dp) WLDAUAV
!.. this routine does not need irregular wavefunctions
LOGICAL LIRRSOL
PARAMETER ( LIRRSOL = .TRUE. )
!..
!.. Local Arrays ..
complex (kind=dp) ALPHA(0:LMAXD), &
               FZ(IRMD,0:LMAXD), &
               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD), &
               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),TMAT(0:LMAXD), &
               ALPHALL(LMMAXD,LMMAXD),DALPHALL(LMMAXD,LMMAXD), &  ! LLY alpha matrix and derivative
               ALPHA0(LMMAXD,LMMAXD),AUX(LMMAXD,LMMAXD), &         ! LLY alpha-matrix
               TRALPHA,TRALPHA1  ! LLY Trace of alpha^-1 * d alpha/dE for Lloyds formula
real (kind=dp) WLDAU(MMAXD,MMAXD,NSPIND)
real (kind=dp) RS(IRMD,0:LMAXD),SL(0:LMAXD)
real (kind=dp) CUTOFF(IRMD)
INTEGER IPIV(LMMAXD) ! LLY
CHARACTER*9 TXTS(2)
#ifdef CPP_MPI
integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
#endif
integer :: ie_end, ie_num, ie_start
!..
!.. External Functions ..
LOGICAL TEST
EXTERNAL TEST
!..
!.. External Subroutines ..
EXTERNAL CRADWF,PNSTMAT,WFMESH,CMATSTR,DRVRELTMAT,ZGEINV1,ZGEMM
!..
!.. Data Statements
DATA TXTS /'spin   UP','spin DOWN'/
!..................................................................

! save timing information to array, needed to tackle load imbalance adaptively



! ==LDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAU

IF ( idoldau == 1 ) THEN
  wldauav = 0.d0
  lmlo = lopt*lopt + 1
  lmhi = (lopt+1)*(lopt+1)
  mmax = lmhi - lmlo + 1
  DO m1 = 1,mmax
    wldauav = wldauav + wldau(m1,m1,ispin)
  END DO
  wldauav = wldauav/DBLE(mmax)
  
  
! -> Note: Application if WLDAU makes the potential discontinuous.
!    A cutoff can be used if necessary to make the potential continuous
!    for example (array bounds should be adjusted):
  
!ccc            CUTOFF(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) *
!ccc     &                   ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
!ccc            CUTOFF(IR) = 1D0/CUTOFF(IR)
  
  DO m1 = 1,irmd
    cutoff(m1) = 1.d0
  END DO
endif

! ==LDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAU

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
IF(myrank==master .AND. t_inc%i_write>0) WRITE(1337,*) 'atom: ',i1
#ifdef CPP_MPI
CALL distribute_linear_on_tasks(t_mpi_c_grid%nranks_at,  &
    t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at+(i1-1), ! print this info only for first atom at master  &
    master,ielast,ntot_pt,ioff_pt,.true.)

ie_start = ioff_pt(t_mpi_c_grid%myrank_at)
ie_end   = ntot_pt(t_mpi_c_grid%myrank_at)

t_mpi_c_grid%ntot2=ie_end !t_mpi_c_grid%dims(1)
IF(.NOT. (allocated(t_mpi_c_grid%ntot_pt2) .OR.  &
    allocated(t_mpi_c_grid%ioff_pt2)))  &
    allocate(t_mpi_c_grid%ntot_pt2(0:t_mpi_c_grid%nranks_at-1),  &
    t_mpi_c_grid%ioff_pt2(0:t_mpi_c_grid%nranks_at-1))
t_mpi_c_grid%ntot_pt2 = ntot_pt
t_mpi_c_grid%ioff_pt2 = ioff_pt
! now initialize arrays for tmat, gmat, and gref
CALL init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
IF(lly /= 0) CALL init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)
#else
! now initialize arrays for tmat, gmat, and gref
CALL init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
IF(lly /= 0) CALL init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)

ie_start = 0
ie_end = ielast
#endif

DO ie_num=1,ie_end
  
  ie = ie_num+ie_start
  
#ifdef CPP_MPI
!start timing measurement for this pair of ie and i1, needed for MPIadapt
  IF(mpiadapt) CALL timing_start('time_1a_ieiatom')
#endif

IF(t_inc%i_write>0) WRITE(1337,*)'CALCTMAT: IE=',ie,' ATOM:',i1

tmatll(:,:) = czero
dtmatll(:,:) = czero

alphall(:,:) = czero
dalphall(:,:) = czero

! In case of Lloyds formula the derivative of t is needed.
! Then calculate t at E+dE, E-dE and average for t, subtract for dt/dE
ideriv = 0                         ! LLY
IF (lly /= 0) ideriv = 1           ! LLY
DO signde = -ideriv,ideriv,2       ! LLY
  pz(:,:) = czero
  qz(:,:) = czero
  fz(:,:) = czero
  sz(:,:) = czero
  pns(:,:,irmind:irmd,:) = czero
  tmat0(:,:) = czero
  alpha0(:,:) = czero
  
  eryd = ez(ie) + signde * deltae / 2.d0 ! LLY
  IF(t_inc%i_write>0) WRITE(1337,*) 'energy:',ie,'',eryd
  
!=======================================================================
! non/scalar-relativistic OR relativistic
  
  IF ( krel == 0 ) THEN
    CALL wfmesh(eryd,ek,cvlight,nsra,zat,rmesh,sl,rs, ircut(ipan),irmd,lmaxd)
    CALL cradwf(eryd,ek,nsra,alpha,ipan,ircut,cvlight,rs,sl,  &
        pz,fz,qz,sz,tmat,visp,drdi,rmesh,zat,lirrsol, idoldau,lopt,wldauav,cutoff)
    
    DO lm1 = 1,lmmaxd                      ! LLY
      alpha0(lm1,lm1) = alpha(loflm(lm1)) ! LLY spherical (diag.) part of alpha
    END DO                                 ! LLY
!-----------------------------------------------------------------------
! spherical/non-spherical
    
    IF ( ins == 0 ) THEN
      DO lm1 = 1,lmmaxd
        tmat0(lm1,lm1) = tmat(loflm(lm1))
      END DO
    ELSE
      CALL pnstmat(drdi,ek,icst,pz,qz,fz,sz,pns,tmat0,  &
          vins,irmin,ipan,ircut,nsra,cleb,icleb,iend,loflm, &     ! Added IRMIN 1.7.2014  &
          tmat,lmaxd,idoldau,lopt,lmlo,lmhi, wldau(1,1,ispin),wldauav,cutoff,  &
          alpha0)                 ! LLY In goes diag. alpha, out comes full alpha
    endif
    
!-----------------------------------------------------------------------
!=======================================================================
!commented this out because compiling with debug options did not work
!due to an erro in the interface to drvreltmat > talk to long,
!philipp 20_08_2015
  ELSE
    STOP 'WARNING check DRVRELTMAT interface in code!!!'
!                CALL DRVRELTMAT(ERYD,TMAT0,VTREL,BTREL,RMREL,
!      &                      DRDIREL,R2DRDIREL,ZREL,JWSREL,SOLVER,SOCTL,
!      &                      CTL,LMMAXD,LMAXD,IRMD)
  endif
  
! non/scalar-relativistic OR relativistic
!=======================================================================
  
! In case of derivative calculation, average the t-matrix
! of the two energies and calculate the difference
  tmatll(:,:) = tmatll(:,:)  + tmat0(:,:)          ! LLY
  alphall(:,:) =  alphall(:,:)  + alpha0(:,:)      ! LLY
  IF (lly /= 0) THEN
    dtmatll(:,:) =dtmatll(:,:) + signde * tmat0(:,:)      ! LLY
    dalphall(:,:) = dalphall(:,:) + signde * alpha0(:,:)  ! LLY
  endif
  
END DO  ! SIGNDE = -IDERIV,IDERIV,2      ! LLY

! Average values of t-matrix and alpha at E+DE and E-DE
tmatll(:,:) = tmatll(:,:) / real(1+ideriv, kind=dp)      !/ 1 or 2 (IDERIV=1 or 2)! LLY
alphall(:,:) = alphall(:,:) / real(1+ideriv, kind=dp)    !/ 1 or 2 ! LLY

IF (lly /= 0) THEN
! Construct derivative of t-matrix and alpha
  dtmatll(:,:) = dtmatll(:,:) / deltae
  dalphall(:,:) = dalphall(:,:) / deltae
! LLY Calculate Trace [alpha^-1 * d alpha/dE] for Lloyd's formula
  alpha0(:,:) = czero
  aux(:,:) = czero
! ALPHA0 = ALPHALL^-1
  CALL zgeinv1(alphall,alpha0,aux,ipiv,lmmaxd)
! AUX = ALPHALL^-1 * DALPHALL
  CALL zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,alpha0,lmmaxd,  &
      dalphall,lmmaxd,czero,aux,lmmaxd)
! Trace of AUX is Trace of [alpha^-1 * d alpha/dE]
  tralpha = czero
  DO lm1 = 1,lmmaxd
    tralpha = tralpha + aux(lm1,lm1)
  END DO
  
  IF (lly == -1) THEN ! Calculate Tr(t^-1 dt/dE) instead of Tr(alpha^-1 * d alpha/de)
    alpha0(:,:) = tmatll(:,:)
    tmat0(:,:) = czero
    aux(:,:) = czero
! ALPHA0 = TMATLL^-1
    CALL zgeinv1(alpha0,tmat0,aux,ipiv,lmmaxd)
! Build trace
    tralpha1 = czero
    DO lm1 = 1,lmmaxd
      DO lm2 = 1,lmmaxd
        tralpha1 = tralpha1 + tmat0(lm1,lm2) * dtmatll(lm2,lm1)
      END DO
    END DO
    tralpha1 = czero
    DO lm1 = 1,lmmaxd
      tralpha1 = tralpha1 + dtmatll(lm1,lm1)/tmatll(lm1,lm1)
    END DO
    
  endif
  
endif ! (LLY.NE.0)


tmat0(:,:) = tmatll(:,:)
irec = ie + ielast*(ispin-1) + ielast*nspin* (i1-1)
IF (t_tgmat%tmat_to_file) THEN
  WRITE(69,REC=irec) tmat0
ELSE
#ifdef CPP_MPI
  irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*  &
      (i1-t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie)-1)
#else
  irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i1-1)
#endif
t_tgmat%tmat(:,:,irec) = tmat0
endif
IF (lly /= 0) THEN                                          ! LLY
  IF(t_lloyd%dtmat_to_file) THEN
    irec = ie + ielast*(ispin-1) + ielast*nspin* (i1-1)
    WRITE(691,REC=irec) dtmatll                                ! LLY
  ELSE
    irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i1-1)
    
    t_lloyd%dtmat(:,:,irec) = dtmatll(:,:)
  endif ! t_lloyd%dtmat_to_file
  IF(t_lloyd%tralpha_to_file) THEN
    irec = ie + ielast*(ispin-1) + ielast*nspin* (i1-1)
    WRITE(692,REC=irec) tralpha                              ! LLY
  ELSE
    irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i1-1)
    t_lloyd%tralpha(irec) = tralpha
  endif
endif                                                       ! LLY

! ----------------------------------------------------------------------
IF ( test('tmat    ') .AND. (t_inc%i_write>0)) THEN
  WRITE (1337,*)
  WRITE (1337,99001, advance='no') '-----> t matrix for atom: ',i1
  IF ( krel == 0 ) WRITE (1337,99002, advance='no') txts(ispin)
  WRITE (1337,99003) ', energy: ',eryd
  CALL cmatstr(' ',1,tmat0,lmmaxd,lmmaxd, 2*krel+1,2*krel+1,0,1D-8,6)
  WRITE (1337,*)
endif
! ----------------------------------------------------------------------

#ifdef CPP_MPI
!stop timing measurement for this pair of ie and i1, needed for MPIadapt
IF(mpiadapt) CALL timing_stop('time_1a_ieiatom',save_out=timings_1a(ie, i1))
#endif

END DO ! IE = 1,IELAST
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

99001 FORMAT (a,i3)
99002 FORMAT (', ',a)
99003 FORMAT (a,2F10.6)

RETURN
END SUBROUTINE calctmat
