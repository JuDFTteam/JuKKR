!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module NonCollinearMagnetism_mod
!-------------------------------------------------------------------------------
!> Summary: Single site solver based on direct inversion supporting non-collinear magnetism
!> Author: Marcel Bornemann, David S G Bauer
!> Category: KKRnano, single-site, solver
!>
!> ToDo: adopt coding style to real F90
!-------------------------------------------------------------------------------
use RadialMeshData_mod!, only:
use ChebMeshData_mod!, only:
implicit none
private

public :: tmat_newsolver
public :: rhovalnew
public :: rotatematrix

  contains

SUBROUTINE drvbastrans(rc,crel,rrel,srrel,nrrel,irrel,  &
    nlmax,nkmmax,nmuemax,nkmpmax,nkmax,linmax)
!   ********************************************************************
!   *                                                                  *
!   *                                                                  *
!   ********************************************************************
IMPLICIT REAL*8(a-h,o-z)

COMPLEX*16, INTENT(IN OUT)               :: rc(nkmmax,nkmmax)
COMPLEX*16, INTENT(IN OUT)               :: crel(nkmmax,nkmmax)
COMPLEX*16, INTENT(IN OUT)               :: rrel(nkmmax,nkmmax)
COMPLEX*16, INTENT(IN OUT)               :: srrel(2,2,nkmmax)
INTEGER, INTENT(IN OUT)                  :: nrrel(2,nkmmax)
INTEGER, INTENT(IN OUT)                  :: irrel(2,2,nkmmax)
INTEGER, INTENT(IN)                      :: nlmax
INTEGER, INTENT(IN)                      :: nkmmax
INTEGER, INTENT(IN)                      :: nmuemax
INTEGER, INTENT(IN)                      :: nkmpmax
INTEGER, INTENT(IN)                      :: nkmax
INTEGER, INTENT(IN)                      :: linmax

!*** Start of declarations rewritten by SPAG

! Local variables

REAL*8 cgc(nkmpmax,2)
INTEGER :: i,ikm1lin(linmax),ikm2lin(linmax),il,imue,iprint,  &
    kaptab(nmuemax),ltab(nmuemax),mmax,nmuetab(nmuemax), nsollm(nlmax,nmuemax)

!*** End of declarations rewritten by SPAG

IF (nkmmax /= 2*nlmax**2) STOP ' Check NLMAX,NKMMAX in < DRVBASTRANS > '
IF (nmuemax /= 2*nlmax) STOP ' Check NLMAX,NMUEMAX in < DRVBASTRANS > '
IF (nkmpmax /= (nkmmax+2*nlmax))  &
    STOP ' Check NLMAX,NKMMAX,NKMPMAX in < DRVBASTRANS > '
IF (nkmax /= 2*nlmax-1) STOP ' Check NLMAX,NKMAX in < DRVBASTRANS > '
IF (linmax /= (2*nlmax*(2*nlmax-1)))  &
    STOP ' Check NLMAX,LINMAX in < DRVBASTRANS > '

iprint = 0

DO i = 1,nmuemax
  ltab(i) = i/2
  IF ( 2*ltab(i) == i ) THEN
    kaptab(i) = ltab(i)
  ELSE
    kaptab(i) = -ltab(i) - 1
  END IF
  nmuetab(i) = 2*ABS(kaptab(i))
END DO

DO il = 1,nlmax
  mmax = 2*il
  DO imue = 1,mmax
    IF ( (imue == 1) .OR. (imue == mmax) ) THEN
      nsollm(il,imue) = 1
    ELSE
      nsollm(il,imue) = 2
    END IF
  END DO
END DO

CALL ikmlin(iprint,nsollm,ikm1lin,ikm2lin,nlmax,nmuemax,linmax, nlmax)

CALL calccgc(ltab,kaptab,nmuetab,cgc,nkmax,nmuemax,nkmpmax)

! ---------------------------- now calculate the transformation matrices

CALL strsmat(nlmax-1,cgc,srrel,nrrel,irrel,nkmmax,nkmpmax)

CALL bastrmat(nlmax-1,cgc,rc,crel,rrel,nkmmax,nkmpmax)

RETURN
END SUBROUTINE drvbastrans

SUBROUTINE changerep(a,mode,b,n,m,rc,crel,rrel,text,ltext)
!   ********************************************************************
!   *                                                                  *
!   *   change the representation of matrix A and store in B           *
!   *   according to MODE:                                             *
!   *                                                                  *
!   *   RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)        *
!   *   REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.          *
!   *   CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)      *
!   *   REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.        *
!   *   RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm. *
!   *   CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm. *
!   *                                                                  *
!   *   the non-relat. representations include the  spin index         *
!   *                                                                  *
!   *   for LTEXT > 0 the new matrix  B  is printed                    *
!   *                                                                  *
!   ********************************************************************
IMPLICIT REAL*8(a-h,o-z)


COMPLEX*16, INTENT(IN OUT)               :: a(m,m)
CHARACTER (LEN=7), INTENT(IN)            :: mode
COMPLEX*16, INTENT(IN OUT)               :: b(m,m)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: m
COMPLEX*16, INTENT(IN OUT)               :: rc(m,m)
COMPLEX*16, INTENT(IN OUT)               :: crel(m,m)
COMPLEX*16, INTENT(IN OUT)               :: rrel(m,m)
CHARACTER (LEN=*), INTENT(IN)        :: text
INTEGER, INTENT(IN)                      :: ltext

!*** Start of declarations rewritten by SPAG

! PARAMETER definitions

COMPLEX*16, PARAMETER :: c1=(1.0D0,0.0D0)
COMPLEX*16, PARAMETER :: c0=(0.0D0,0.0D0)

! Dummy arguments






! Local variables

INTEGER :: key
COMPLEX*16 w1(m,m)

!*** End of declarations rewritten by SPAG

!---------------------- transform MAT from (kappa,mue) to REAL (l,ml,ms)
IF ( mode == 'REL>RLM' ) THEN
  CALL zgemm('N','N',n,n,n,c1,rrel,m,a,m,c0,w1,m)
  CALL zgemm('N','C',n,n,n,c1,w1,m,rrel,m,c0,b,m)
  key = 2
ELSE IF ( mode == 'RLM>REL' ) THEN
  CALL zgemm('C','N',n,n,n,c1,rrel,m,a,m,c0,w1,m)
  CALL zgemm('N','N',n,n,n,c1,w1,m,rrel,m,c0,b,m)
  key = 3
ELSE IF ( mode == 'REL>CLM' ) THEN
  CALL zgemm('N','N',n,n,n,c1,crel,m,a,m,c0,w1,m)
  CALL zgemm('N','C',n,n,n,c1,w1,m,crel,m,c0,b,m)
  key = 2
ELSE IF ( mode == 'CLM>REL' ) THEN
  CALL zgemm('C','N',n,n,n,c1,crel,m,a,m,c0,w1,m)
  CALL zgemm('N','N',n,n,n,c1,w1,m,crel,m,c0,b,m)
  key = 3
ELSE IF ( mode == 'CLM>RLM' ) THEN
  CALL zgemm('N','N',n,n,n,c1,rc,m,a,m,c0,w1,m)
  CALL zgemm('N','C',n,n,n,c1,w1,m,rc,m,c0,b,m)
  key = 2
ELSE IF ( mode == 'RLM>CLM' ) THEN
  CALL zgemm('C','N',n,n,n,c1,rc,m,a,m,c0,w1,m)
  CALL zgemm('N','N',n,n,n,c1,w1,m,rc,m,c0,b,m)
  key = 2
ELSE
  WRITE (*,*) ' MODE = ',mode
  STOP 'in <ROTATE>  MODE not allowed'
END IF

IF ( ltext > 0 ) CALL cmatstr(text,ltext,b,n,m,key,key,0,1D-8,6)
!     IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-12,6)
END SUBROUTINE changerep

SUBROUTINE create_newmesh_new(r_log,npan_log,npan_eq,ncheb,  &
                              npan_lognew,npan_eqnew,  &
                              npan_tot,rnew,rpan_intervall,ipan_intervall,  &
                              thetasnew,thetas,nfu,radial_mesh) ! new parameters
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-03-29  Time: 17:20:39

!use DimParams_mod
!use InputParams_mod
use read_formatted_shapefun_mod, only: ShapefunFile

!INTEGER, INTENT(IN)                      :: nspin
DOUBLE PRECISION, INTENT(IN)             :: r_log 
INTEGER, INTENT(IN)                      :: npan_log
INTEGER, INTENT(IN)                      :: npan_eq
INTEGER, INTENT(IN)                      :: ncheb
INTEGER, INTENT(OUT)                     :: npan_lognew
INTEGER, INTENT(OUT)                     :: npan_eqnew
INTEGER, INTENT(OUT)                     :: npan_tot
DOUBLE PRECISION, INTENT(OUT)            :: rnew(:)
!DOUBLE PRECISION, INTENT(IN)            :: thetas(irid,nfund,ncelld)
DOUBLE PRECISION, INTENT(OUT)            :: rpan_intervall(0:)
INTEGER, INTENT(OUT)                     :: ipan_intervall(0:)
DOUBLE PRECISION, INTENT(OUT)            :: thetasnew(:,:)
DOUBLE PRECISION, INTENT(IN)             :: thetas(:,:)
INTEGER, INTENT(IN)                      :: nfu
type(RadialMeshData), INTENT(IN)         :: radial_mesh
!type(ShapefunFile), intent(in)           :: sfile
!type(InputParams), intent(in)            :: params
!type(DimParams), intent(in)              :: dims

!IMPLICIT NONE
!INCLUDE 'inc.p'


!INTEGER, PARAMETER :: lmmaxd= (lmaxd+1)**2

!INTEGER, PARAMETER :: lmpotd= (lpotd+1)**2

!INTEGER, PARAMETER :: irmind= irmd-irnsd
!INTEGER :: npan_inst
!INTEGER :: npan_lognew(dims%naez)



DOUBLE PRECISION, PARAMETER :: fac=2D0
INTEGER :: ipotm,ir2,ip,  &
    ishift,ilogpanshift,ilinpanshift,npan_logtemp,npan_inst,imin,imax,iminnew,imaxnew,lm1
DOUBLE PRECISION :: rmin,rmax,rval
!type(RadialMeshData) :: radial_mesh

! set up radial mesh
!call createRadialMeshData(radial_mesh, dims%irmd, sfile%mesh(i1)%npan+1)

! local arrays to be filled with data from 'sfile'
! ---------------------------------------------
!DOUBLE PRECISION, allocatable       :: r(:,:)
!INTEGER, allocatable                :: irmin(:)
!INTEGER, allocatable                :: irws(:)
!INTEGER, allocatable                :: ipan(:)
!INTEGER, allocatable                :: ircut(:,:) 
!DOUBLE PRECISION, allocatable       :: thetasin(:,:,:)

!allocate(r(size(sfile%mesh(i1)%r),dims%naez))
!allocate(irmin(dims%naez))
!allocate(irws(dims%naez))
!allocate(ipan(dims%naez))
!allocate(ircut(size(sfile%mesh(i1)%ircut),dims%naez))
!allocate(thetasin(dims%irid,sfile%shapes(1)%nfu,dims%naez)) ! #shapefunctions = #atoms in KKRnano, can nfu be different for different shapes?!

!DO i1 = 1,dims%naez
  
   !call createRadialMeshData(meshdata, irmd, ipand)
   !call initRadialMesh(meshdata, alat, sfile%mesh(cell_index)%xrn, &
   !                    sfile%mesh(cell_index)%drn, sfile%mesh(cell_index)%nm, irmd-irid, irns)
!   r(:,i1)   = sfile%mesh(i1)%r
!   irmin(i1) = sfile%mesh(i1)%irmin
!   irws(i1)  = sfile%mesh(i1)%irws
!   ipan(i1)  = sfile%mesh(i1)%ipan
!   ircut(:,i1)  = sfile%mesh(i1)%ircut
!   DO lm1=1,sfile%shapes(i1)%nfu
!    thetasin(:,lm1,i1)=sfile%shapes(i1)%thetas(:,lm1) ! #shapefunctions = #atoms in KKRnano 
!   ENDDO
!ENDDO
! ----------------------------------------------


thetasnew=0D0
ipotm=0

!DO i1 = 1,dims%naez ! JM: NATYPD instead of NAEZ
 
  ! set up radial mesh for atom i1
!  call createRadialMeshData(radial_mesh, dims%irmd, sfile%mesh(i1)%npan+1)
!  call initRadialMesh(radial_mesh, params%alat, sfile%mesh(i1)%xrn, &
!                      sfile%mesh(i1)%drn, sfile%mesh(i1)%nm, &
!                      dims%irmd-dims%irid, dims%irnsd)
 
!  ipot=dims%nspind*(i1-1)+1
  npan_inst = radial_mesh%ipan-1
  npan_tot  = npan_log+npan_eq+npan_inst
  
  
! log panel
  rmin=radial_mesh%r(2)
  rmax=r_log
  rval=0D0
  ishift=0
  IF (r_log > radial_mesh%r(radial_mesh%irmin)) THEN
    ilogpanshift=1
    ilinpanshift=0
  ELSE
    ilogpanshift=0
    ilinpanshift=1
  END IF
  
  IF (ilinpanshift == 1) THEN
    STOP 'non-spherical part of the potential needs to be inside the log panel'
  END IF
  
  DO ip=0,npan_log-ilogpanshift
    rval=(fac**ip-1D0)/(fac**(npan_log-ilogpanshift)-1D0)
    rpan_intervall(ip+ishift)= rmin+rval*(rmax-rmin)
    ipan_intervall(ip+ishift)= (ip+ishift)*(ncheb+1)
    IF (ishift == 0.AND. rpan_intervall(ip) > radial_mesh%r(radial_mesh%irmin)) THEN
      ishift=1
      npan_logtemp=ip
      rpan_intervall(ip+1)=rpan_intervall(ip)
      ipan_intervall(ip+1)=(ip+ishift)*(ncheb+1)
      rpan_intervall(ip)=radial_mesh%r(radial_mesh%irmin)
      ipan_intervall(ip)=ip*(ncheb+1)
    END IF
  END DO ! NPAN_LOG
  
! equivalent panel
  ishift=0
  rmin=r_log
  rmax=radial_mesh%r(radial_mesh%ircut(1))
  DO ip=0,npan_eq-ilinpanshift
    rpan_intervall(ip+ishift+npan_log)=rmin+ip*(rmax-rmin)/  &
        (npan_eq-ilinpanshift)
    ipan_intervall(ip+ishift+npan_log)=(npan_log+ip+ishift)* (ncheb+1)
  END DO ! NPAN_EQ
  
! intersection zone
  DO ip=1,npan_inst
    rpan_intervall(npan_log+npan_eq+ip)=radial_mesh%r(radial_mesh%ircut(ip+1))
    ipan_intervall(npan_log+npan_eq+ip)=(npan_log+npan_eq+ip)* (ncheb+1)
  END DO ! NPAN_INST
  
  npan_eqnew=npan_eq+npan_log-npan_logtemp
  npan_lognew=npan_logtemp
  
  CALL chebmesh(npan_tot,ncheb,rpan_intervall(0:), rnew(1:))
  
! interpolate shape function THETAS to new shape function THETASNEW
! save THETAS to THETASIN
!  icell = i1 ! #shapefunctions = #atoms in KKRnano, JM: ICELL = NTCELL(I1)
  DO lm1=1,nfu
!    thetasin(:,lm1,icell)=sfile%shapes(icell)%thetas(:,lm1) ! get thetas for specific atom and lm1
    ir2=0
    DO ip=npan_lognew+npan_eqnew+1,npan_tot
      ir2=ir2+1
      imin=radial_mesh%ircut(ir2)+1
      imax=radial_mesh%ircut(ir2+1)
      iminnew=ipan_intervall(ip-1)+1
      imaxnew=ipan_intervall(ip)
      CALL interpolspline(radial_mesh%r(imin:imax),rnew(iminnew:imaxnew),  &
           thetas(imin-radial_mesh%ircut(1):imax-radial_mesh%ircut(1),lm1),  &
           thetasnew(iminnew:imaxnew,lm1), imax-imin+1,imaxnew-iminnew+1)
    END DO
  END DO

!call destroyRadialMeshData(radial_mesh)
!END DO ! I1
END SUBROUTINE create_newmesh_new


SUBROUTINE chebmesh(npan,ncheb,ri,ro)

INTEGER, INTENT(IN)                      :: npan
INTEGER, INTENT(IN)                      :: ncheb
DOUBLE PRECISION, INTENT(IN)             :: ri(0:npan)
DOUBLE PRECISION, INTENT(OUT)            :: ro(npan*(ncheb+1))
!IMPLICIT NONE



INTEGER :: i,k,ik
DOUBLE PRECISION :: tau,pi

pi=4D0*DATAN(1D0)
DO i=1,npan
  DO k=0,ncheb
    ik=i*ncheb+i-k
    tau=DCOS(((2*k+1)*pi)/(2*(ncheb+1)))
    tau=0.5D0*((ri(i)-ri(i-1))*tau+ri(i)+ri(i-1))
    ro(ik)=tau
  END DO
END DO
END SUBROUTINE chebmesh


SUBROUTINE bastrmat(lmax,cgc,rc,crel,rrel,nkmmax,nkmpmax)
!   ********************************************************************
!   *                                                                  *
!   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
!   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
!   *                                                                  *
!   *    this is a special version of <STRSMAT> passing the            *
!   *    full BASis TRansformation MATrices  RC, CREL and RREL         *
!   *                                                                  *
!   * 13/01/98  HE                                                     *
!   ********************************************************************

IMPLICIT REAL*8(a-h,o-z)

INTEGER, INTENT(IN)                      :: lmax
REAL*8, INTENT(IN)                       :: cgc(nkmpmax,2)
COMPLEX*16, INTENT(OUT)                  :: rc(nkmmax,nkmmax)
COMPLEX*16, INTENT(OUT)                  :: crel(nkmmax,nkmmax)
COMPLEX*16, INTENT(IN OUT)               :: rrel(nkmmax,nkmmax)
INTEGER, INTENT(IN)                  :: nkmmax
INTEGER, INTENT(IN)                  :: nkmpmax

!*** Start of declarations rewritten by SPAG

! PARAMETER definitions

COMPLEX*16, PARAMETER :: ci=(0.0D0,1.0D0)
COMPLEX*16, PARAMETER :: c1=(1.0D0,0.0D0)
COMPLEX*16, PARAMETER :: c0=(0.0D0,0.0D0)

! Local variables

INTEGER :: i,ikm,j,jp05,k,l,lm,lnr,m,muem05,muep05,nk,nkm,nlm
REAL*8 w

!*** End of declarations rewritten by SPAG

nk = 2*(lmax+1) + 1
nlm = (lmax+1)**2
nkm = 2*nlm
!     ===================================================
!     INDEXING:
!     IKM  = L*2*(J+1/2) + J + MUE + 1
!     LM   = L*(L+1)     +     M   + 1
!     ===================================================

! ----------------------------------------------------------------------
! CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
! ----------------------------------------------------------------------
CALL cinit(nkmmax*nkmmax,crel)

lm = 0
DO lnr = 0,lmax
  DO m = -lnr,lnr
    lm = lm + 1
    
    ikm = 0
    DO k = 1,nk
      l = k/2
      IF ( 2*l == k ) THEN
        jp05 = l
      ELSE
        jp05 = l + 1
      END IF
      
      DO muem05 = -jp05,(jp05-1)
        muep05 = muem05 + 1
        ikm = ikm + 1
        
        IF ( l == lnr ) THEN
          IF ( muep05 == m ) crel(lm,ikm) = cgc(ikm,1)
          IF ( muem05 == m ) crel(lm+nlm,ikm) = cgc(ikm,2)
        END IF
        
      END DO
    END DO
    
  END DO
END DO

! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
CALL cinit(nkmmax*nkmmax,rc)

w = 1.0D0/SQRT(2.0D0)

DO l = 0,lmax
  DO m = -l,l
    i = l*(l+1) + m + 1
    j = l*(l+1) - m + 1
    
    IF ( m < 0 ) THEN
      rc(i,i) = -ci*w
      rc(j,i) = w
      rc(i+nlm,i+nlm) = -ci*w
      rc(j+nlm,i+nlm) = w
    END IF
    IF ( m == 0 ) THEN
      rc(i,i) = c1
      rc(i+nlm,i+nlm) = c1
    END IF
    IF ( m > 0 ) THEN
      rc(i,i) = w*(-1.0D0)**m
      rc(j,i) = ci*w*(-1.0D0)**m
      rc(i+nlm,i+nlm) = w*(-1.0D0)**m
      rc(j+nlm,i+nlm) = ci*w*(-1.0D0)**m
    END IF
  END DO
END DO

! ----------------------------------------------------------------------
! RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
! ----------------------------------------------------------------------

CALL zgemm('N','N',nkm,nkm,nkm,c1,rc,nkmmax,crel,nkmmax,c0,rrel, nkmmax)

END SUBROUTINE bastrmat

SUBROUTINE calccgc(ltab,kaptab,nmuetab,cgc,nkmax,nmuemax,nkmpmax)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:05:10
 
!   ********************************************************************
!   *                                                                  *
!   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
!   *                                                                  *
!   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
!   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
!   *   IS= 1/2  SPIN DOWN/UP                                          *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ltab(nmuemax)
INTEGER, INTENT(IN)                      :: kaptab(nmuemax)
INTEGER, INTENT(IN)                      :: nmuetab(nmuemax)
REAL*8, INTENT(OUT)                      :: cgc(nkmpmax,2)
INTEGER, INTENT(IN)                      :: nkmax
INTEGER, INTENT(IN)                      :: nmuemax
INTEGER, INTENT(IN)                      :: nkmpmax


! Local variables

INTEGER :: ikm,k,kappa,m
REAL*8 j,l,mue,twolp1

ikm = 0
DO k = 1,(nkmax+1)
  l = ltab(k)
  kappa = kaptab(k)
  j = ABS(kappa) - 0.5D0
  mue = -j - 1.0D0
  twolp1 = 2.0D0*l + 1.0D0
  
  IF ( kappa < 0 ) THEN
    
!     J = L + 1/2
    DO m = 1,nmuetab(k)
      
      mue = mue + 1.0D0
      ikm = ikm + 1
      cgc(ikm,1) = DSQRT((l-mue+0.5D0)/twolp1)
      cgc(ikm,2) = DSQRT((l+mue+0.5D0)/twolp1)
    END DO
  ELSE
!     J = L - 1/2
    DO m = 1,nmuetab(k)
      
      mue = mue + 1.0D0
      ikm = ikm + 1
      cgc(ikm,1) = DSQRT((l+mue+0.5D0)/twolp1)
      cgc(ikm,2) = -DSQRT((l-mue+0.5D0)/twolp1)
      
    END DO
  END IF
  
  
END DO

END SUBROUTINE calccgc

!*==cmatstr.f    processed by SPAG 6.05Rc at 15:50 on 12 Oct 2002
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:05:17

SUBROUTINE cmatstr(str,lstr,a,n,m,mlin,mcol,ijq,tolp,k_fmt_fil)
!   ********************************************************************
!   *                                                                  *
!   *   writes structure of COMPLEX   NxN   matrix   A                 *
!   *                                                                  *
!   *   M           is the actual array - size used for   A            *
!   *   MLIN/COL    MODE for line and column indexing                  *
!   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
!   *   TOL         tolerance for difference                           *
!   *   IJQ         if IJQ > 1000    pick  IQ-JQ-block matrix          *
!   *               assuming  IJQ = IQ*1000 + JQ                       *
!   *               else: no IQ-JQ-indexing                            *
!   *   K_FMT_FIL   output channel                                     *
!   *               a negative sign suppresses table at the end        *
!   *                                                                  *
!   *   any changes should be done in RMATSTR as well !!!!!!!!!!!!!!!  *
!   *                                                                  *
!   ********************************************************************

IMPLICIT COMPLEX*16(a-h,o-z)

CHARACTER (LEN=*), INTENT(IN)            :: str
INTEGER, INTENT(IN)                      :: lstr
COMPLEX*16, INTENT(IN OUT)               :: a(m,m)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: mlin
INTEGER, INTENT(IN)                      :: mcol
INTEGER, INTENT(IN)                      :: ijq
REAL*8, INTENT(IN)                       :: tolp
INTEGER, INTENT(IN)                      :: k_fmt_fil

!*** Start of declarations rewritten by SPAG

! PARAMETER definitions

COMPLEX*16, PARAMETER :: ci=(0.0D0,1.0D0)

! Local variables

COMPLEX*16 b(n,n),ca,cb,arg,dtab(0:n*n)
CHARACTER (LEN=1) :: CHAR
LOGICAL :: same,small
CHARACTER (LEN=1) :: ctab(0:n*n),vz(-1:+1)
DOUBLE PRECISION :: DBLE
CHARACTER (LEN=150) :: fmt1,fmt2,fmt3,fmt4
INTEGER :: i,i1,ic0,id,il,ilsep(20),ipt(218),iq,isl,iw(m),j,  &
    j0,jp,jq,k,l3,lf,mm,n1,n2,n3,nc,nd,nfil,nk,nm,nm1,nm2,nm3, nnon0,nsl
INTEGER :: ICHAR,ISIGN,nint
REAL*8 tol

!*** End of declarations rewritten by SPAG

DATA vz/'-',' ',' '/

small(arg) = ABS(arg*tol) < 1.0D0

same(ca,cb) = small(1.0D0-ca/cb)

nfil = ABS(k_fmt_fil)

tol = 1.0D0/tolp

!----------------------------------------------- set block indices IQ JQ

IF ( ijq > 1000 ) THEN
  iq = ijq/1000
  jq = ijq - iq*1000
  IF ( iq*n > m .OR. iq*n > m ) THEN
    WRITE (6,99002) ijq,iq,jq,iq*n,jq*n,n,m
    RETURN
  END IF
ELSE
  iq = 1
  jq = 1
END IF

!----------------------------------------------------- copy matrix block

j0 = n*(jq-1)
DO j = 1,n
  i1 = n*(iq-1)+1
  jp = j0 + j
  CALL zcopy(n,a(i1,jp),1,b(1,j),1)
END DO

!------------------------------------------------ set up character table

nc = 0
DO i = 1,26
  nc = nc + 1
  ipt(nc) = 62 + i
END DO
DO i = 1,8
  nc = nc + 1
  ipt(nc) = 96 + i
END DO
DO i = 10,26
  nc = nc + 1
  ipt(nc) = 96 + i
END DO
DO i = 191,218
  nc = nc + 1
  ipt(nc) = i
END DO
DO i = 35,38
  nc = nc + 1
  ipt(nc) = i
END DO
DO i = 40,42
  nc = nc + 1
  ipt(nc) = i
END DO
DO i = 91,93
  nc = nc + 1
  ipt(nc) = i
END DO

!---------------------------------------------------------------- header
ic0 = ICHAR('0')
n3 = n/100
n2 = n/10 - n3*10
n1 = n - n2*10 - n3*100

IF ( n <= 18 ) THEN
  fmt1 = '(8X,I3,''|'','
  fmt2 = '( 9X,''--|'','
  fmt3 = '( 9X,'' #|'','
  fmt4 = '( 9X,''  |'','
ELSE
  fmt1 = '(   I4,''|'','
  fmt2 = '( 2X,''--|'','
  fmt3 = '( 2X,'' #|'','
  fmt4 = '( 2X,''  |'','
END IF

lf = 11
l3 = 11
IF ( mcol == 0 ) THEN
  fmt1 = fmt1(1:lf)//CHAR(ic0+n3)//CHAR(ic0+n2)//CHAR(ic0+n1)  &
      //'( 2A1),''|'',I3)'
  fmt2 = fmt2(1:lf)//CHAR(ic0+n3)//CHAR(ic0+n2)//CHAR(ic0+n1)  &
      //'(''--''),''|'',I3)'
  fmt3 = fmt3(1:lf)//'60(2X,I2))'
  fmt4 = fmt4(1:lf)//'60(I2,2X))'
  lf = 21
ELSE
  IF ( mcol == 1 ) THEN
    nk = nint(SQRT(DBLE(n)))
  ELSE IF ( mcol == 2 ) THEN
    nk = nint(SQRT(DBLE(n/2)))
  ELSE IF ( mcol == 3 ) THEN
    nk = 2*nint(SQRT(DBLE(n/2))) - 1
  END IF
  DO k = 1,nk
    IF ( mcol <= 2 ) THEN
      nm = 2*k - 1
    ELSE
      nm = 2*((k+1)/2)
    END IF
    nm2 = nm/10
    nm1 = nm - nm2*10
    nm3 = nm/2
    fmt1 = fmt1(1:lf)//CHAR(ic0+nm2)//CHAR(ic0+nm1) //'( 2A1),''|'','
    fmt2 = fmt2(1:lf)//CHAR(ic0+nm2)//CHAR(ic0+nm1) //'(''--''),''|'','
    
    IF ( mcol <= 2 ) THEN
      DO mm = 1,nm
        IF ( MOD(mm,2) == MOD(k,2) ) THEN
          fmt3 = fmt3(1:l3)//'2X,'
          fmt4 = fmt4(1:l3)//'I2,'
        ELSE
          fmt3 = fmt3(1:l3)//'I2,'
          fmt4 = fmt4(1:l3)//'2X,'
        END IF
        l3 = l3 + 3
      END DO
      fmt3 = fmt3(1:l3)//'''|'','
      fmt4 = fmt4(1:l3)//'''|'','
      l3 = l3 + 4
    ELSE
      fmt3 = fmt3(1:lf)//CHAR(ic0+nm3)//'(2X,I2),''|'','
      fmt4 = fmt4(1:lf)//CHAR(ic0+nm3)//'(I2,2X),''|'','
      l3 = l3 + 13
    END IF
    lf = lf + 13
  END DO
  IF ( mcol == 2 ) THEN
    fmt1 = fmt1(1:lf)//fmt1(12:lf)
    fmt2 = fmt2(1:lf)//fmt2(12:lf)
    fmt3 = fmt3(1:l3)//fmt3(12:l3)
    fmt4 = fmt4(1:l3)//fmt4(12:l3)
    lf = 2*lf - 11
    l3 = 2*l3 - 11
  END IF
  fmt1 = fmt1(1:lf)//'I3)'
  fmt2 = fmt2(1:lf)//'I3)'
  fmt3 = fmt3(1:l3)//'I3)'
  fmt4 = fmt4(1:l3)//'I3)'
END IF
IF ( mlin == 0 ) THEN
  nsl = 1
  ilsep(1) = n
ELSE IF ( mlin == 1 ) THEN
  nsl = nint(SQRT(DBLE(n)))
  DO il = 1,nsl
    ilsep(il) = il**2
  END DO
ELSE IF ( mlin == 2 ) THEN
  nsl = nint(SQRT(DBLE(n/2)))
  DO il = 1,nsl
    ilsep(il) = il**2
  END DO
  DO il = 1,nsl
    ilsep(nsl+il) = ilsep(nsl) + il**2
  END DO
  nsl = 2*nsl
ELSE IF ( mlin == 3 ) THEN
  nsl = 2*nint(SQRT(DBLE(n/2))) - 1
  ilsep(1) = 2
  DO k = 2,nsl
    ilsep(k) = ilsep(k-1) + 2*((k+1)/2)
  END DO
END IF


WRITE (nfil,99001) str(1:lstr)
IF ( ijq > 1000 ) WRITE (nfil,99003) iq,jq
WRITE (nfil,fmt3) (i,i=2,n,2)
WRITE (nfil,fmt4) (i,i=1,n,2)
WRITE (nfil,FMT=fmt2)
!------------------------------------------------------------ header end
nnon0 = 0
nd = 0
ctab(0) = ' '
dtab(0) = 9999D0

DO i = 1,n
  DO j = 1,n
    IF ( .NOT.small(b(i,j)) ) THEN
      nnon0 = nnon0 + 1
      DO id = 1,nd
        IF ( same(b(i,j),+dtab(id)) ) THEN
          iw(j) = +id
          GO TO 50
        END IF
        IF ( same(b(i,j),-dtab(id)) ) THEN
          iw(j) = -id
          GO TO 50
        END IF
      END DO
!----------------------------------------------------------- new element
      nd = nd + 1
      iw(j) = nd
      dtab(nd) = b(i,j)
      IF ( ABS(dtab(nd)-1.0D0)*tol < 1.0D0 ) THEN
        ctab(nd) = '1'
      ELSE IF ( ABS(dtab(nd)+1.0D0)*tol < 1.0D0 ) THEN
        dtab(nd) = +1.0D0
        ctab(nd) = '1'
        iw(j) = -nd
      ELSE IF ( ABS(dtab(nd)-ci)*tol < 1.0D0 ) THEN
        ctab(nd) = 'i'
      ELSE IF ( ABS(dtab(nd)+ci)*tol < 1.0D0 ) THEN
        dtab(nd) = +ci
        ctab(nd) = 'i'
        iw(j) = -nd
      ELSE
        ctab(nd) = CHAR(ipt(1+MOD((nd+1),nc)))
      END IF
    ELSE
      iw(j) = 0
    END IF
  50      END DO
!------------------------------------------------------------ write line
  WRITE (nfil,FMT=fmt1) i, (vz(ISIGN(1,iw(j))),ctab(ABS(iw(j))),j=1,  &
      n),i
  
  DO isl = 1,nsl
    IF ( i == ilsep(isl) ) WRITE (nfil,FMT=fmt2)
  END DO
END DO

!------------------------------------------------------------------ foot

WRITE (nfil,fmt4) (i,i=1,n,2)
WRITE (nfil,fmt3) (i,i=2,n,2)

IF ( k_fmt_fil > 0 ) THEN
  WRITE (nfil,99004) (id,ctab(id),dtab(id),id=1,nd)
  WRITE (nfil,99005) nnon0,tolp,n*n - nnon0,tolp
ELSE
  WRITE (nfil,*) ' '
END IF

99001 FORMAT (/,8X,a,/)
99002 FORMAT (/,1X,79('*'),/,10X,'inconsistent call of <CMATSTR>',/,10X,  &
    'argument IJQ =',i8,'  implies IQ=',i3,'   JQ=',i3,/,10X,  &
    'IQ*N=',i6,' > M   or   JQ*N=',i6,' > M   for N =',i4,  &
    ' M=',i4,/,1X,79('*'),/)
99003 FORMAT (8X,'IQ-JQ-block  for  IQ = ',i3,'   JQ = ',i3,/)
99004 FORMAT (/,8X,'symbols used:',/,(8X,i3,3X,a1,2X,2F20.12))
99005 FORMAT (/,8X,i5,' elements   >',1PE9.1,/,  &
    8X,i5,' elements   <',1PE9.1,/)
END SUBROUTINE cmatstr

FUNCTION ikapmue(kappa,muem05)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:21:58

!   ********************************************************************
!   *                                                                  *
!   *  INDEXING OF MATRIX-ELEMENTS:                                    *
!   *                                                                  *
!   *  I = 2*L*(J+1/2) + J + MUE + 1                                   *
!   *                                                                  *
!   ********************************************************************


IMPLICIT NONE 

INTEGER, INTENT(IN)                      :: kappa
INTEGER, INTENT(IN)                      :: muem05


! Dummy arguments


INTEGER :: ikapmue

! Local variables

INTEGER :: IABS 
INTEGER :: jp05,l

jp05 = IABS(kappa)

IF ( kappa < 0 ) THEN 
  l = -kappa - 1
ELSE
  l = kappa
END IF

ikapmue = 2*l*jp05 + jp05 + muem05 + 1

END FUNCTION ikapmue


SUBROUTINE ikmlin(iprint,nsollm,ikm1lin,ikm2lin,nlmax,nmuemax,  &
        linmax,nl)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:05:20

!   ********************************************************************
!   *                                                                  *
!   * SETUP TABLE OF INDICES    IKM(INT)                               *
!   *                                                                  *
!   *  IKM IS STANDARD INDEX IN  (KAPPA,MUE)-REPRESENTATION            *
!   *  IKM = 2*L*(J+1/2) + J + MUE + 1                                 *
!   *                                                                  *
!   *  INT NUMBERS LINEARLY ONLY NON-VANISHING ELEMENTS OF M-SS        *
!   *  USED TO CALCULATE DOS ...                                       *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: iprint
INTEGER, INTENT(IN)                      :: nsollm(nlmax,nmuemax)
INTEGER, INTENT(OUT)                     :: ikm1lin(linmax)
INTEGER, INTENT(OUT)                     :: ikm2lin(linmax)
INTEGER, INTENT(IN)                      :: nlmax
INTEGER, INTENT(IN)                      :: nmuemax
INTEGER, INTENT(IN)                      :: linmax
INTEGER, INTENT(IN)                      :: nl


! Dummy arguments




! Local variables

INTEGER :: i,il,imue,k1,k2,kap(2),l,lin,muem05,nsol
!INTEGER :: ikapmue

lin = 0

DO il = 1,nl
  l = il - 1
  muem05 = -il - 1
  kap(1) = -l - 1
  kap(2) = +l
  
  DO imue = 1,2*il
    muem05 = muem05 + 1
    nsol = nsollm(il,imue)
    
    DO k2 = 1,nsol
      DO k1 = 1,nsol
        lin = lin + 1
        ikm1lin(lin) = ikapmue(kap(k1),muem05)
        ikm2lin(lin) = ikapmue(kap(k2),muem05)
      END DO
    END DO
    
  END DO
END DO

IF ( iprint < 2 ) RETURN
WRITE (6,FMT='('' INT='',I3,''  IKM=('',I3,'','',I3,'')'')')  &
    (i,ikm1lin(i),ikm2lin(i),i=1,lin)
END SUBROUTINE ikmlin

!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:26:23

SUBROUTINE splint_real(xa,ya,y2a,n,x,y,yderiv)

IMPLICIT NONE

REAL*8, INTENT(IN)                       :: xa(*)
REAL*8, INTENT(IN)                       :: ya(*)
REAL*8, INTENT(IN OUT)                   :: y2a(*)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN)                       :: x
REAL*8, INTENT(OUT)                      :: y
REAL*8, INTENT(OUT)                      :: yderiv


! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the xai's in order), and given the array y2a(1:n), which
! is the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y and the derivative yderiv.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
INTEGER :: k,khi,klo
REAL*8         a,b,h
! We will  nd the right place in the table by means of bisection.
! This is optimal if sequential calls to this routine are at random
! values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the
! next call.
klo=1
khi=n
1    IF (khi-klo > 1) THEN
  k=(khi+klo)/2
  IF(xa(k) > x)THEN
    khi=k
  ELSE
    klo=k
  END IF
  GO TO 1
END IF
! klo and khi now bracket the input value of x.
h=xa(khi)-xa(klo)
! The xa's must be distinct.
IF (h == 0.d0) STOP 'bad xa input in splint' !! used to be PAUSE
! Cubic spline polynomial is now evaluated.
a = (xa(khi)-x)/h
b = (x-xa(klo))/h
y = a*ya(klo) + b*ya(khi) +  &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (h**2)/6.d0
yderiv = (ya(khi)-ya(klo))/h -  &
    ((3.d0*a*a-1.d0)*y2a(klo) - (3.d0*b*b-1.d0)*y2a(khi))*h/6.d0

RETURN
END SUBROUTINE splint_real


SUBROUTINE interpolspline(rmesh,rmeshnew,vpot,vpotnew,  &
        nrmax,nrmaxnew)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:05:24
 
IMPLICIT NONE
!interface
INTEGER :: nrmax
INTEGER :: nrmaxnew
DOUBLE PRECISION :: rmesh(nrmax)
DOUBLE PRECISION :: rmeshnew(nrmaxnew)
DOUBLE PRECISION :: vpot(nrmax)
DOUBLE PRECISION :: vpotnew(nrmaxnew)
!local
DOUBLE PRECISION :: maxa
DOUBLE PRECISION :: spline(nrmax)
DOUBLE PRECISION :: parsum, parsumderiv,r0
INTEGER :: ir
maxa = 1.d35
CALL spline_real(nrmax,rmesh,vpot,nrmax,maxa,maxa,spline)
!           CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)

DO ir = 1,nrmaxnew
  r0 = rmeshnew(ir)
  CALL splint_real(rmesh,vpot,spline,nrmax,r0,parsum,parsumderiv)
  vpotnew(ir) = parsum
END DO
END SUBROUTINE  interpolspline

SUBROUTINE strsmat(lmax,cgc,srrel,nrrel,irrel,nkmmax,nkmpmax)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:05:34
 
!   ********************************************************************
!   *                                                                  *
!   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
!   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
!   *                                                                  *
!   *    ONLY THE NON-0 ELEMENTS OF THE MATRIX ARE STORED              *
!   *                                                                  *
!   * 25/10/95  HE  proper convention of trans. matrix introduced      *
!   ********************************************************************

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: lmax
REAL*8, INTENT(IN)                       :: cgc(nkmpmax,2)
COMPLEX*16, INTENT(OUT)                  :: srrel(2,2,nkmmax)
INTEGER, INTENT(OUT)                     :: nrrel(2,nkmmax)
INTEGER, INTENT(OUT)                     :: irrel(2,2,nkmmax)
INTEGER, INTENT(IN)                  :: nkmmax
INTEGER, INTENT(IN)                  :: nkmpmax

! PARAMETER definitions

COMPLEX*16, PARAMETER :: ci=(0.0D0,1.0D0)
COMPLEX*16, PARAMETER :: c1=(1.0D0,0.0D0)
COMPLEX*16, PARAMETER :: c0=(0.0D0,0.0D0)

! Dummy arguments






! Local variables

COMPLEX*16 crel(nkmmax,nkmmax),rc(nkmmax,nkmmax), rrel(nkmmax,nkmmax)
INTEGER :: i,ikm,j,jp05,k,l,lam,lm,lnr,lr,m,muem05,muep05,nk,nkm,nlm, ns1,ns2
REAL*8 w

nk = 2*(lmax+1) + 1
nlm = (lmax+1)**2
nkm = 2*nlm
!     ===================================================
!     INDEXING:
!     IKM  = L*2*(J+1/2) + J + MUE + 1
!     LM   = L*(L+1)     +     M   + 1
!     ===================================================

! ----------------------------------------------------------------------
! CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
! ----------------------------------------------------------------------
CALL cinit(nkmmax*nkmmax,crel)

lm = 0
DO lnr = 0,lmax
  DO m = -lnr,lnr
    lm = lm + 1
    
    ikm = 0
    DO k = 1,nk
      l = k/2
      IF ( 2*l == k ) THEN
        jp05 = l
      ELSE
        jp05 = l + 1
      END IF
      
      DO muem05 = -jp05,(jp05-1)
        muep05 = muem05 + 1
        ikm = ikm + 1
        
        IF ( l == lnr ) THEN
          IF ( muep05 == m ) crel(lm,ikm) = cgc(ikm,1)
          IF ( muem05 == m ) crel(lm+nlm,ikm) = cgc(ikm,2)
        END IF
        
      END DO
    END DO
    
  END DO
END DO

! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
CALL cinit(nkmmax*nkmmax,rc)

w = 1.0D0/SQRT(2.0D0)

DO l = 0,lmax
  DO m = -l,l
    i = l*(l+1) + m + 1
    j = l*(l+1) - m + 1
    
    IF ( m < 0 ) THEN
      rc(i,i) = -ci*w
      rc(j,i) = w
      rc(i+nlm,i+nlm) = -ci*w
      rc(j+nlm,i+nlm) = w
    END IF
    IF ( m == 0 ) THEN
      rc(i,i) = c1
      rc(i+nlm,i+nlm) = c1
    END IF
    IF ( m > 0 ) THEN
      rc(i,i) = w*(-1.0D0)**m
      rc(j,i) = ci*w*(-1.0D0)**m
      rc(i+nlm,i+nlm) = w*(-1.0D0)**m
      rc(j+nlm,i+nlm) = ci*w*(-1.0D0)**m
    END IF
  END DO
END DO

! ----------------------------------------------------------------------
! RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
! ----------------------------------------------------------------------
CALL zgemm('N','N',nkm,nkm,nkm,c1,rc,nkmmax,crel,nkmmax,c0,rrel, nkmmax)

!     ---------------------------------------------------
!     store the elements of  RREL
!     ---------------------------------------------------
DO lam = 1,nkm
  ns1 = 0
  ns2 = 0
  
  DO lr = 1,2*nlm
    IF ( CDABS(rrel(lr,lam)) > 1D-6 ) THEN
      IF ( lr <= nlm ) THEN
        ns1 = ns1 + 1
        IF ( ns1 > 2 ) STOP ' IN <STRSMAT>   NS1 > 2'
        srrel(ns1,1,lam) = rrel(lr,lam)
        irrel(ns1,1,lam) = lr
      ELSE
        ns2 = ns2 + 1
        IF ( ns2 > 2 ) STOP ' IN <STRSMAT>   NS2 > 2'
        srrel(ns2,2,lam) = rrel(lr,lam)
        irrel(ns2,2,lam) = lr - nlm
      END IF
    END IF
  END DO
  
  nrrel(1,lam) = ns1
  nrrel(2,lam) = ns2
END DO

END SUBROUTINE strsmat


!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-01  Time: 12:22:03

SUBROUTINE spline_real(nmax,x,y,n,yp1,ypn,y2)

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nmax
REAL*8, INTENT(IN)                       :: x(nmax)
REAL*8, INTENT(IN)                       :: y(nmax)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN OUT)                   :: yp1
REAL*8, INTENT(IN OUT)                   :: ypn
REAL*8, INTENT(OUT)                      :: y2(nmax)


!      REAL*8           x(NMAX)
!      COMPLEX          yp1,ypn,y(NMAX),y2(NMAX)
! Given arrays x(1:n) and  y(1:n) containing a tabulated function,
! i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn
! for the 1rst derivative of the interpolating function at points
! 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating
! function at the tabulated points xi.
! If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
! signaled to set the corresponding boundary condition for a natural
! spline, with zero second derivative on that boundary.
! Parameter: NMAX is the largest anticipated value of n.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
INTEGER :: i,k
REAL*8          p,qn,sig,un,u(nmax)
!      COMPLEX          p,qn,sig,un,u(NMAX)

IF (n > nmax) STOP 'SPLINE: n > NMAX.'
IF (ABS(yp1) > 0.99D30) THEN
! The lower boundary condition is set either to be "natural"
  y2(1) = 0.d0
  u(1) = 0.d0
ELSE
! or else to have a specified first derivative.
  y2(1) = -0.5D0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
END IF

DO i = 2,n-1
! This is the decomposition loop of the tridiagonal algorithm. y2 and u
! are used for temporary storage of the decomposed factors.
  sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
  p = sig * y2(i-1) + 2.d0
  y2(i) = (sig-1.d0)/p
  u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))  &
      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1)) / p
END DO

IF (ABS(ypn) > 0.99D30) THEN
! The upper boundary condition is set either to be "natural"
  qn = 0.d0
  un = 0.d0
ELSE
! or else to have a specified 1rst derivative.
  qn = 0.5D0
  un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
END IF
y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0)
DO k = n-1,1,-1
! This is the backsubstitution loop of the tridiagonal algorithm.
  y2(k)=y2(k)*y2(k+1)+u(k)
END DO

RETURN
END SUBROUTINE spline_real


SUBROUTINE tmat_newsolver(ie,nspin,lmax,zat,socscale,  &
        ez,nsra,cleb,icleb,iend,ncheb,npan_tot,  &
        rpan_intervall,ipan_intervall,  &
        rnew,vinsnew,theta,phi,ipot,  &
       ! lly,        &
        lmpotd,irmd_new,TmatN,soc) ! new input parameters
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-18  Time: 14:58:02
 
#ifdef cpp_omp
!use omp_lib        ! necessary for omp functions
#endif
#ifdef cpp_mpi
!use mpi
#endif
!use mod_mympi, only: myrank, nranks, master
#ifdef cpp_mpi
!                   & ,distribute_linear_on_tasks
#endif
!use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid,init_tgmat,  &
!    t_lloyd,init_tlloyd

!!use JijDij_mod, only: type_dtmatJijDij, init_t_dtmatJij_at, calc_dtmatJij

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ie
!INTEGER, INTENT(IN)                      :: ielast
INTEGER, INTENT(IN)                      :: nspin
INTEGER, INTENT(IN)                      :: lmax
!DOUBLE PRECISION, INTENT(IN)             :: rmesh(:)
DOUBLE PRECISION, INTENT(IN)             :: zat
DOUBLE PRECISION, INTENT(IN)             :: socscale
DOUBLE COMPLEX, INTENT(IN)               :: ez(:)
INTEGER, INTENT(IN)                      :: nsra
DOUBLE PRECISION, INTENT(IN)             :: cleb(:)
INTEGER, INTENT(IN)                      :: icleb(:,:)
INTEGER, INTENT(IN)                      :: iend
INTEGER, INTENT(IN)                      :: ncheb
INTEGER, INTENT(IN)                      :: npan_tot
DOUBLE PRECISION, INTENT(IN)             :: rpan_intervall(0:)
INTEGER, INTENT(IN)                      :: ipan_intervall(0:)
DOUBLE PRECISION, INTENT(IN)             :: rnew(:)
DOUBLE PRECISION, INTENT(IN)             :: vinsnew(:,:,:)
DOUBLE PRECISION, INTENT(IN)             :: theta
DOUBLE PRECISION, INTENT(IN)             :: phi
!INTEGER, INTENT(IN)                      :: i1
INTEGER, INTENT(IN)                      :: ipot
!INTEGER, INTENT(IN)                      :: lly
!DOUBLE COMPLEX, INTENT(IN)               :: deltae
INTEGER, INTENT(IN)                      :: lmpotd
!INTEGER, INTENT(IN)                      :: lmaxd
INTEGER, INTENT(IN)                      :: irmd_new
DOUBLE COMPLEX, INTENT(OUT)              :: TmatN(:,:)
LOGICAL, INTENT(IN)                      :: soc
!INCLUDE 'inc.p'


INTEGER :: lmmaxd
INTEGER :: lmmaxso
INTEGER :: nrmaxd

DOUBLE COMPLEX eryd

DOUBLE PRECISION, PARAMETER :: cvlight=274.0720442D0
DOUBLE COMPLEX, PARAMETER :: czero=(0D0,0D0)
DOUBLE COMPLEX, PARAMETER :: cone=(1D0,0D0)


DOUBLE COMPLEX, allocatable ::  tmatll(:,:)
INTEGER :: ir,use_sratrick,nvec,lm1,irmdnew
DOUBLE COMPLEX gmatprefactor
DOUBLE PRECISION, allocatable :: vins(:,:,:)
DOUBLE COMPLEX, allocatable :: vnspll0(:,:,:),vnspll1(:,:,:,:), vnspll(:,:,:,:)
DOUBLE COMPLEX, allocatable :: hlk(:,:,:),jlk(:,:,:), hlk2(:,:,:),jlk2(:,:,:)
DOUBLE COMPLEX, allocatable :: rll(:,:,:,:)
!DOUBLE COMPLEX, allocatable :: rllleft(:,:,:,:),sllleft(:,:,:,:) ! neded for D_ij calculation
DOUBLE COMPLEX, allocatable :: tmatsph(:,:)! TMAT_OUT(:,:), tmat_out necessary for parallel ie loop
DOUBLE COMPLEX, allocatable :: dtmatll(:,:),tmat0(:,:) ! LLY
DOUBLE COMPLEX, allocatable :: alphall(:,:),dalphall(:,:),alpha0(:,:),aux(:,:)         ! LLY
!DOUBLE COMPLEX, allocatable :: alphasph(:)!, DTMAT_OUT(:,:,:), ! LLY
INTEGER, allocatable        :: jlk_index(:)
! LLoyd:
!INTEGER :: ideriv,signde        ! LLY
!DOUBLE COMPLEX              :: tralpha            ! LLY
DOUBLE COMPLEX, allocatable :: ipiv(:)            ! LLY
!     .. OMP ..
INTEGER :: nth,ith         ! total number of threads and thread id

lmmaxd = (lmax+1)**2
lmmaxso=2*lmmaxd
nrmaxd=irmd_new

allocate(tmatll(lmmaxso,lmmaxso))
allocate(dtmatll(lmmaxso,lmmaxso))
allocate(tmat0(lmmaxso,lmmaxso))
allocate(alphall(lmmaxso,lmmaxso))
allocate(dalphall(lmmaxso,lmmaxso))
allocate(alpha0(lmmaxso,lmmaxso))
allocate(aux(lmmaxso,lmmaxso))
allocate(jlk_index(2*lmmaxso))
allocate(ipiv(lmmaxso))
!     .. OMP ..
! determine if omp parallelisation is used (compiled with -openmp flag and
! OMP_NUM_THREADS>1)
!$noomp parallel shared(nth,ith)
!$noomp single
nth = 1
ith = 0
!nth = omp_get_num_threads()
!$noomp end single
!$noomp end parallel
! write(*,*) 'nth =',nth

irmdnew= npan_tot*(ncheb+1)
allocate(vins(irmdnew,lmpotd,nspin))
vins=0D0
DO lm1=1,lmpotd
  DO ir=1,irmdnew
    vins(ir,lm1,1)=vinsnew(ir,lm1,ipot)
    vins(ir,lm1,nspin)=vinsnew(ir,lm1,ipot+nspin-1)
  END DO
END DO
!c set up the non-spherical ll' matrix for potential VLL'
 IF (NSRA.EQ.2) THEN
USE_SRATRICK=1
ELSEIF (NSRA.EQ.1) THEN
USE_SRATRICK=0
ENDIF
allocate(vnspll0(lmmaxso,lmmaxso,irmdnew))
allocate(vnspll1(lmmaxso,lmmaxso,irmdnew,0:nth-1))
vnspll0=czero
CALL vllmat(1,irmdnew,lmmaxd,lmmaxso,vnspll0,vins,  &
    cleb,icleb,iend,nspin,zat,rnew,use_sratrick)

! initial allocate
IF (nsra == 2) THEN
  allocate(vnspll(2*lmmaxso,2*lmmaxso,irmdnew,0:nth-1))
ELSE
  allocate(vnspll(lmmaxso,lmmaxso,irmdnew,0:nth-1))
END IF

allocate(hlk(1:4*(lmax+1),irmdnew,0:nth-1))
allocate(jlk(1:4*(lmax+1),irmdnew,0:nth-1))
allocate(hlk2(1:4*(lmax+1),irmdnew,0:nth-1))
allocate(jlk2(1:4*(lmax+1),irmdnew,0:nth-1))
allocate(tmatsph(2*(lmax+1),0:nth-1))
allocate(rll(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1))
!allocate(rllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1))
!allocate(sllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1))
!allocate(tmat_out(lmmaxso,lmmaxso))

! energy loop
!WRITE(6,*) 'atom: ',i1,' NSRA:',nsra

!!$noomp parallel do default(none)
!!$noomp& private(eryd,ie,i1,ir,irec,nvec,lm1,lm2,gmatprefactor)
!!$noomp& private(jlk_index,tmatll,ith)
!!$noomp& shared(nspin,nsra,lmax,iend,ipot,ielast,npan_tot,ncheb)
!!$noomp& shared(zat,socscale,ez,rmesh,cleb,rnew,nth)
!!$noomp& shared(rpan_intervall,vinsnew,ipan_intervall)
!!$noomp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)
!!$noomp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,tmat_out)
!!$noomp& shared(tmatsph)
!DO ie=1,ielast
! get current thread
!  IF (nth>=1) THEN
!    ith = omp_get_thread_num()
!  ELSE
    ith = 0
!  END IF
  eryd = ez(ie)
!!$noomp critical
!  WRITE(6,*) 'energy:',ie,'',eryd
!write(*,*) 'nested omp?',omp_get_nested()
!!$noomp end critical
  
! contruct the spin-orbit coupling hamiltonian and add to potential
  CALL spinorbit_ham(lmax,lmmaxd,vins,rnew,  &
      eryd,zat,cvlight,socscale,nspin,lmpotd,  &
      theta,phi,ipan_intervall,rpan_intervall, npan_tot,ncheb,irmdnew,nrmaxd,  &
      vnspll0(:,:,:),vnspll1(:,:,:,ith),'1',soc)
!c extend matrix for the SRA treatment
  vnspll(:,:,:,ith)=czero
  IF (nsra == 2) THEN
    IF (use_sratrick == 0) THEN
      CALL vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,  &
          lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'Ref=0')
    ELSE IF (use_sratrick == 1) THEN
      CALL vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,  &
          lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'Ref=Vsph')
    END IF
  ELSE
    vnspll(:,:,:,ith)=vnspll1(:,:,:,ith)
  END IF
  
!c calculate the source terms in the Lippmann-Schwinger equation
!c these are spherical hankel and bessel functions
  hlk(:,:,ith)=czero
  jlk(:,:,ith)=czero
  hlk2(:,:,ith)=czero
  jlk2(:,:,ith)=czero
  gmatprefactor=czero
  CALL rllsllsourceterms(nsra,nvec,eryd,rnew,irmdnew,nrmaxd,lmax,  &
      lmmaxso,1,jlk_index,hlk(:,:,ith),  &
      jlk(:,:,ith),hlk2(:,:,ith),jlk2(:,:,ith), gmatprefactor)
!c using spherical potential as reference
  IF (use_sratrick == 1) THEN
    CALL calcsph(nsra,irmdnew,nrmaxd,lmax,nspin,zat,cvlight,eryd,  &
        rnew,vins,ncheb,npan_tot,rpan_intervall,  &
        jlk_index,hlk(:,:,ith),jlk(:,:,ith),hlk2(:,:,ith),  &
        jlk2(:,:,ith),gmatprefactor,tmatsph(:,ith), use_sratrick)
  END IF
  
!c calculate the tmat and wavefunctions
  rll(:,:,:,ith)=czero
  
!c right solutions
  tmatll=czero
  CALL rll_only(rpan_intervall,rnew,vnspll(:,:,:,ith),  &
      rll(:,:,:,ith),tmatll(:,:),ncheb,  &
      npan_tot,lmmaxso,nvec*lmmaxso,4*(lmax+1),irmdnew,  &
      nsra,jlk_index,hlk(:,:,ith),jlk(:,:,ith),  &
      hlk2(:,:,ith),jlk2(:,:,ith),gmatprefactor, '1','1',use_sratrick)
!     &              ,ith) ! test fivos
!  IF (nsra == 2) THEN
!         RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
!     +            RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/C
!  END IF
!if(t_dtmatjij_at%calculate) then


  
 !for Jij-tensor calculation: allocate array to hold additional t-matrices
!  call init_t_dtmatJij_at(t_dtmatJij_at)
!  
!
!!       lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
!!       lllllllllll calculate the left-hand side solution lllllllllllllllllllllllllllllllllllll
!!       contruct the spin-orbit coupling hamiltonian and add to potential
!   call spinorbit_ham(lmax,lmmaxd,vins,rnew, &
!                      eryd,zat,cvlight,socscale,nsra,nspin,lmpotd, &
!                      theta,phi,ipan_intervall,rpan_intervall, &
!                      npan_tot,ncheb,irmdnew,nrmaxd, &
!                      vnspll0(:,:,:),vnspll1(:,:,:,ith), &
!                      'transpose',soc)
!
!!       extend matrix for the sra treatment
!   vnspll(:,:,:,ith)=czero
!   if (nsra.eq.2) then
!    if (use_sratrick.eq.0) then
!     call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew, &
!          lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'ref=0')
!    elseif (use_sratrick.eq.1) then
!     call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew, &
!         lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'ref=vsph')
!    endif
!   else
!    vnspll(:,:,:,ith)=vnspll1(:,:,:,ith)
!   endif
!
!!       calculate the source terms in the lippmann-schwinger equation
!!       these are spherical hankel and bessel functions
!   hlk(:,:,ith)=czero
!   jlk(:,:,ith)=czero
!   hlk2(:,:,ith)=czero
!   jlk2(:,:,ith)=czero
!   gmatprefactor=czero
!   jlk_index = 0
!   call rllsllsourceterms(nsra,nvec,eryd,rnew,irmdnew,nrmaxd,lmax, &
!                         lmmaxso,1,jlk_index,hlk(:,:,ith), &
!                         jlk(:,:,ith),hlk2(:,:,ith),jlk2(:,:,ith), &
!                         gmatprefactor)
!
!!       using spherical potential as reference
!!        notice that exchange the order of left and right hankel/bessel functions
!   if (use_sratrick.eq.1) then
!    tmatsph(:,ith)=czero
!    call calcsph(nsra,irmdnew,nrmaxd,lmax,nspin,zat,cvlight,eryd, &
!                lmpotd,lmmaxso,rnew,vins,ncheb,npan_tot,rpan_intervall, &
!                jlk_index,hlk2(:,:,ith),jlk2(:,:,ith),hlk(:,:,ith), &
!                jlk(:,:,ith),gmatprefactor,tmatsph(:,ith), &
!                use_sratrick)
!   endif
!   
!!       calculate the tmat and wavefunctions
!   rllleft(:,:,:,ith)=czero
!   sllleft(:,:,:,ith)=czero
!
!!       left solutions
!!        notice that exchange the order of left and right hankel/bessel functions
!   tmat0=czero
!   alpha0=czero ! lly
!   call rllsll(rpan_intervall,rnew,vnspll(:,:,:,ith), &
!               rllleft(:,:,:,ith),sllleft(:,:,:,ith),tmat0,ncheb, &
!               npan_tot,lmmaxso,nvec*lmmaxso,4*(lmax+1),irmdnew, &
!               nrmaxd,nsra,jlk_index,hlk2(:,:,ith),jlk2(:,:,ith), &
!               hlk(:,:,ith),jlk(:,:,ith),gmatprefactor, &
!               '1','1','0',use_sratrick)
!   if (nsra.eq.2) then
!    rllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)= &
!             rllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)/cvlight
!    sllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)= &
!             sllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)/cvlight
!   endif
!!       lllllllllll calculate the left-hand side solution lllllllllllllllllllllllllllllllllllll
!!       lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
!
!  call calc_dtmatjij(lmaxd,lmmaxd,lmmaxso,lmpotd,ntotd,nrmaxd, &
!         nsra,irmdnew,nspin,vins,rllleft(:,:,:,ith),rll(:,:,:,ith), &
!         rpan_intervall, &
!         ipan_intervall,npan_tot,ncheb,cleb,icleb,iend,ncleb,rnew, &
!         theta,phi,t_dtmatjij_at%dtmat_xyz(:,:,:,ie_num))

!  end if!t_dtmatjij_at%calculate
  
  
! add spherical contribution of tmatrix
  IF (use_sratrick == 1) THEN
    DO lm1=1,lmmaxso
      tmatll(lm1,lm1)=tmatll(lm1,lm1)+tmatsph(jlk_index(lm1),ith)
    END DO
  END IF
  TmatN(:,:) = tmatll(:,:)
!END DO ! IE loop
!!$noomp end parallel do

! serial write out after parallel calculation of tmat
!DO ie=1,ielast
!  irec = ie + ielast*(i1-1)
!  WRITE(69,REC=irec) tmat_out(:,:,ie)
!         write(696969,*) TMAT_out(:,:,ie)
!END DO

deallocate(vins)
deallocate(vnspll0)
deallocate(vnspll1)
deallocate(vnspll)
deallocate(hlk)
deallocate(jlk)
deallocate(hlk2)
deallocate(jlk2)
deallocate(tmatsph)
deallocate(rll)

END SUBROUTINE tmat_newsolver


! ************************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-18  Time: 14:28:39

SUBROUTINE vllmat(irmin,irc,lmmax,lmmaxso,vnspll0,vins,  &
    cleb,icleb,iend,nspin,z,rnew,use_sratrick)
! ************************************************************************
!     .. Parameters ..
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: irmin
!INTEGER, INTENT(IN)                      :: nrmaxd
INTEGER, INTENT(IN)                      :: irc
INTEGER, INTENT(IN)                      :: lmmax
INTEGER, INTENT(IN)                      :: lmmaxso
DOUBLE COMPLEX, INTENT(OUT)              :: vnspll0(:,:,irmin:)
DOUBLE PRECISION, INTENT(IN OUT)         :: vins(irmin:,:,:)
DOUBLE PRECISION, INTENT(IN)             :: cleb(:)
INTEGER, INTENT(IN)                      :: icleb(:,:)
INTEGER, INTENT(IN)                      :: iend
INTEGER, INTENT(IN)                      :: nspin
DOUBLE PRECISION, INTENT(IN)             :: z
DOUBLE PRECISION, INTENT(IN)             :: rnew(irmin:)
INTEGER, INTENT(IN OUT)                  :: use_sratrick
!INCLUDE 'inc.p'
!INTEGER :: lmpotd
!DOUBLE PRECISION, INTENT, PARAMETER :: lmpotd= (lpotd+1)**2
!     ..
!     .. Scalar Arguments ..

INTEGER :: isp
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION, allocatable :: vnspll(:,:,:,:)

!     ..
!     .. Local Scalars ..
INTEGER :: i,ir,j,lm1,lm2,lm3
!     ..

allocate(vnspll(lmmax,lmmax,irmin:irc,2))

DO isp=1,nspin
  DO  lm1 = 1,lmmax
    DO  lm2 = 1,lm1
      DO  ir = irmin,irc
        vnspll(lm1,lm2,ir,isp) = 0.0D0
      END DO
    END DO
  END DO
  
  DO  j = 1,iend
    lm1 = icleb(j,1)
    lm2 = icleb(j,2)
    lm3 = icleb(j,3)
    DO  i = irmin,irc
      vnspll(lm1,lm2,i,isp) = vnspll(lm1,lm2,i,isp) + cleb(j)*vins(i,lm3,isp)
    END DO
  END DO
  
!---> use symmetry of the gaunt coef.
  
  DO  lm1 = 1,lmmax
    DO  lm2 = 1,lm1 - 1
      DO  i = irmin,irc
        vnspll(lm2,lm1,i,isp) = vnspll(lm1,lm2,i,isp)
      END DO
    END DO
  END DO
  
  IF (use_sratrick == 0) THEN
    DO lm1=1,lmmax
      DO i=irmin,irc
        vnspll(lm1,lm1,i,isp)=vnspll(lm1,lm1,i,isp)+  &
            vins(i,1,isp)-2D0*z/rnew(i)
      END DO
    END DO
  END IF
  
END DO !NSPIN

! set vnspll as twice as large

vnspll0(1:lmmax,1:lmmax,irmin:irc)= vnspll(1:lmmax,1:lmmax,irmin:irc,1)

vnspll0(lmmax+1:lmmaxso,lmmax+1:lmmaxso,irmin:irc)=  &
    vnspll(1:lmmax,1:lmmax,irmin:irc,nspin)
END SUBROUTINE vllmat


SUBROUTINE spinorbit_ham(lmax,lmmaxd,vins,rnew,e,z,c,socscale,  &
        nspin,lmpotd,theta,phi,  &
        ipan_intervall,rpan_intervall,  &
        npan_tot,ncheb,irmdnew,nrmaxd,vnspll,vnspll1,  &
        mode,soc)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-18  Time: 14:28:35

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: lmax
INTEGER, INTENT(IN)                      :: lmmaxd
DOUBLE PRECISION, INTENT(IN)             :: vins(irmdnew,lmpotd,nspin)
DOUBLE PRECISION, INTENT(IN)             :: rnew(nrmaxd)
DOUBLE COMPLEX, INTENT(IN OUT)           :: e
DOUBLE PRECISION, INTENT(IN)             :: z
DOUBLE PRECISION, INTENT(IN)             :: c
DOUBLE PRECISION, INTENT(IN)             :: socscale
!INTEGER, INTENT(IN)                      :: nsra
INTEGER, INTENT(IN)                      :: nspin
INTEGER, INTENT(IN)                      :: lmpotd
DOUBLE PRECISION, INTENT(IN)             :: theta
DOUBLE PRECISION, INTENT(IN)             :: phi
INTEGER, INTENT(IN)                      :: ipan_intervall(0:)
DOUBLE PRECISION, INTENT(IN)             :: rpan_intervall(0:)
INTEGER, INTENT(IN)                      :: npan_tot
INTEGER, INTENT(IN)                      :: ncheb
INTEGER, INTENT(IN)                      :: irmdnew
INTEGER, INTENT(IN OUT)                  :: nrmaxd
DOUBLE COMPLEX, INTENT(IN)               :: vnspll(:,:,:)
DOUBLE COMPLEX, INTENT(OUT)              :: vnspll1(:,:,:)
CHARACTER(LEN=*), INTENT(IN)             :: mode
LOGICAL, INTENT(IN)                      :: soc !switches SOC on and off





DOUBLE PRECISION :: vr(irmdnew),dvdr(irmdnew)
DOUBLE PRECISION :: rmass(irmdnew),hsofac(irmdnew)
DOUBLE PRECISION :: rnucl,atn,widthfac
INTEGER :: ir,ip,lm1,lm2,ispin,irmin,irmax,ncoll
DOUBLE COMPLEX lsmh(2*lmmaxd,2*lmmaxd),temp
DOUBLE PRECISION :: clambdacinv(0:ncheb,0:ncheb)
!DOUBLE PRECISION :: matvec_dmdm
LOGICAL :: test,opt
EXTERNAL test,opt

vnspll1=(0D0,0D0)
vr=0D0
DO ispin=1,nspin
  DO ir=1,ipan_intervall(npan_tot)
    vr(ir)=vr(ir)+vins(ir,1,ispin)/nspin
  END DO
END DO
! derivative of potential
dvdr=0D0
CALL getclambdacinv(ncheb,clambdacinv)
DO ip=1,npan_tot
  irmin=ipan_intervall(ip-1)+1
  irmax=ipan_intervall(ip)
  widthfac= 2D0/(rpan_intervall(ip)-rpan_intervall(ip-1))
  CALL dgemv('N',ncheb+1,ncheb+1,1D0,clambdacinv,ncheb+1,  &
      vr(irmin:irmax),1,0D0,dvdr(irmin:irmax),1)
  dvdr(irmin:irmax)= dvdr(irmin:irmax)*widthfac
END DO
! core potential
IF (z > 24D0) THEN
  atn=-16.1532921+2.70335346*z
ELSE
  atn=0.03467714+2.04820786*z
END IF
rnucl=1.2D0/0.529177D0*atn**(1./3D0)*1.d-5

DO ir=1,ipan_intervall(npan_tot)
  IF (rnew(ir) <= rnucl) THEN
!        DVDR(IR)=DVDR(IR)+2d0*Z*RNEW(IR)/RNUCL**3d0
  ELSE
!        DVDR(IR)=DVDR(IR)+2d0*Z/RNEW(IR)**2d0
  END IF
  dvdr(ir)=dvdr(ir)+2D0*z/rnew(ir)**2D0
END DO
! contruct LS matrix

CALL spin_orbit_compl(lmax,lmmaxd,lsmh)

! roate LS matrix
ncoll=1
IF (ncoll == 1) THEN
  CALL rotatematrix(lsmh,theta,phi,lmmaxd,1)
END IF

IF (mode == 'transpose') THEN
  DO lm1=1,2*lmmaxd
    DO lm2=1,lm1-1
      temp=lsmh(lm2,lm1)
      lsmh(lm2,lm1)=lsmh(lm1,lm2)
      lsmh(lm1,lm2)=temp
    END DO
  END DO
ELSE IF (mode == '1') THEN
END IF
! contruct prefactor of spin-orbit hamiltonian

hsofac=0D0
DO ir=1,irmdnew
  rmass(ir)=0.5D0-0.5D0/c**2*((vr(ir)-REAL(e))-2D0*z/rnew(ir))
  IF (soc .eqv. .false. .OR. z < 1D-6) THEN
    hsofac(ir)=0D0
  ELSE
    hsofac(ir)=socscale/(2D0*rmass(ir)**2*c**2*rnew(ir))*dvdr(ir)
  END IF
  
! add to potential
  
  DO lm1=1,2*lmmaxd
    DO lm2=1,2*lmmaxd
      vnspll1(lm1,lm2,ir)=vnspll(lm1,lm2,ir)+hsofac(ir)*lsmh(lm1,lm2)
    END DO
  END DO
END DO
END SUBROUTINE spinorbit_ham


subroutine vllmatsra(vll0,vll,rmesh,lmsize,nrmax,nrmaxd,eryd,cvlight,lmax,lval_in,cmode)  
!************************************************************************************
! The perubation matrix for the SRA-equations are set up
!************************************************************************************
implicit none
!interface
  DOUBLE COMPLEX VLL(2*lmsize,2*lmsize,nrmax)
  DOUBLE COMPLEX VLL0(lmsize,lmsize,nrmax)
  double precision            :: rmesh(nrmaxd)
  double complex              :: eryd
  double precision            :: cvlight
  integer                     :: lmax,lval_in
  integer                     :: lmsize,nrmax,nrmaxd
  character(len=*)            :: cmode
!local
  integer                     :: ilm,lval,mval,ival,ir
  integer                     :: loflm(lmsize)
  double complex              :: Mass,Mass0
  double complex,parameter    :: cone=(1.0D0,0.0D0)
  double complex,parameter    :: czero=(0.0D0,0.0D0)


!************************************************************************************
! determine the bounds of the matricies to get the lm-expansion and the max. number
! of radial points
!************************************************************************************



!************************************************************************************
! calculate the index array to determine the L value of an LM index
! in case of spin-orbit coupling 2*(LMAX+1)**2 are used instead of (LMAX+1)**2
! the second half refers to the second spin and has the the same L value
!************************************************************************************
ilm=0

if (lmsize==1) then
  loflm(1)=lval_in
elseif ((lmax+1)**2 == lmsize) then
  do lval=0,lmax
    do mval = -lval,lval
      ilm=ilm+1
      loflm(ilm)=lval
    end do
  end do
elseif (2* (lmax+1)**2 ==lmsize ) then
  do ival=1,2
    do lval=0,lmax
      do mval = -lval,lval
        ilm=ilm+1
        loflm(ilm)=lval
      end do
    end do
  end do
else
  stop '[vllmatsra] error'
end if




vll=(0.0D0,0d0)




if     (cmode=='Ref=0') then
  vll(1:lmsize,1:lmsize,:)= vll0 !/cvlight

  do ir=1,nrmax
      do ival=1,lmsize  
        lval=loflm(ival)
        Mass =cone+(eryd-vll0(ival,ival,ir))/cvlight**2
        Mass0=cone+eryd/cvlight**2

  !************************************************************************************
  ! Conventional potential matrix
  !************************************************************************************

       vll(lmsize+ival,lmsize+ival,ir)= -vll0(ival,ival,ir)/cvlight**2 ! TEST 9/22/2011
       vll(ival,ival,ir)=vll(ival,ival,ir)+ (1.0D0/Mass-1.0D0/Mass0)*lval*(lval+1)/rmesh(ir)**2

  !************************************************************************************
  ! The pertubation matrix is changed in the following way
  !
  !     from  / V11  V12 \   to    / V21  V22 \
  !           \ V21  V22 /         \-V11 -V12 / 
  ! because of the convention used for the left solution
  !************************************************************************************
     end do !ival

  end do !ir
elseif     (cmode=='Ref=Vsph') then
 vll(lmsize+1:2*lmsize,1:lmsize,:)=vll0
endif


end subroutine vllmatsra


subroutine rllsllsourceterms(nsra,nvec,eryd,rmesh,nrmax,nrmaxd,lmax,lmsize,use_fullgmat,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR)
implicit none
! ************************************************************************
! calculates the source terms J,H and the left solution J2, H2 for:
! - non-relativistic
! - scalar-relativistic
! - full-relativistic
! calculations
! ************************************************************************
double complex,parameter   :: ci=(0.0d0,1.0d0)
double precision           :: cvlight
parameter (cvlight=274.0720442D0)
integer                    :: nsra,lmax,nrmax,nrmaxd,nvec
double complex             :: eryd
double precision           :: rmesh(nrmaxd)
integer                    :: jlk_index(2*lmsize)
integer                    :: l1,lm1,m1,ivec,ispinfullgmat,ir
integer                    :: use_fullgmat
integer                    :: lmsize

double complex             :: ek,ek2,gmatprefactor
double complex             :: hlk(1:4*(lmax+1),nrmax),jlk(1:4*(lmax+1),nrmax)
double complex             :: hlk2(1:4*(lmax+1),nrmax),jlk2(1:4*(lmax+1),nrmax)

if (nsra==2) then 
  nvec=2
elseif (nsra==1) then 
  nvec=1
end if


  lm1 = 1
  do ivec=1,nvec
    do ispinfullgmat=0,use_fullgmat
      do l1 = 0,lmax
        do m1 = -l1,l1
          jlk_index(lm1) = l1+(ivec-1)*(lmax+1)+1
          lm1 = lm1 + 1
        end do   
      end do  
    end do!ispinorbit=0,use_fullgmat
  end do !nvec


if (nsra==1) then 
  ek = sqrt(eryd)
  ek2 = sqrt(eryd)
elseif (nsra==2) then
  ek = sqrt(eryd+(eryd/cvlight)**2)
  ek2 = sqrt(eryd+(eryd/cvlight)**2) *(1.0d0+eryd/cvlight**2)
end if

                              

do ir = 1,nrmax

    call beshank(hlk(:,ir),jlk(:,ir),ek*rmesh(ir),lmax)
    if (nsra==2) then
      call beshank_smallcomp(hlk(:,ir),jlk(:,ir),&
                        ek*rmesh(ir),rmesh(ir),eryd,lmax)
    end if

    do l1 = 1,nvec*(lmax+1)
      hlk(l1,ir) = -ci*hlk(l1,ir)
    end do

    if (nsra==1) then
      do l1 = 1,nvec*(lmax+1)
        jlk2(l1,ir) = jlk(l1,ir)
        hlk2(l1,ir) = hlk(l1,ir)
      end do
    else if (nsra==2) then
    do l1 = 1,lmax+1
      jlk2(l1,ir) = jlk(l1,ir)
      hlk2(l1,ir) = hlk(l1,ir)
    end do
    do l1 = lmax+2,2*(lmax+1)
      jlk2(l1,ir) = -jlk(l1,ir)
      hlk2(l1,ir) = -hlk(l1,ir)
    end do
    end if

end do
gmatprefactor=ek2
end subroutine rllsllsourceterms


SUBROUTINE calcsph(nsra,irmdnew,nrmaxd,lmax,nspin,z,c,e,  &
        rnew,vins,ncheb,npan_tot,rpan_intervall,  &
        jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,tmat,  &
        use_sratrick)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-18  Time: 14:28:30

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nsra
INTEGER, INTENT(IN)                      :: irmdnew
INTEGER, INTENT(IN OUT)                  :: nrmaxd
INTEGER, INTENT(IN)                      :: lmax
INTEGER, INTENT(IN)                      :: nspin
DOUBLE PRECISION, INTENT(IN)             :: z
DOUBLE PRECISION, INTENT(IN)             :: c
DOUBLE COMPLEX, INTENT(OUT)              :: e
!INTEGER, INTENT(IN)                      :: lmpotd
!INTEGER, INTENT(IN OUT)                  :: lmmaxso
DOUBLE PRECISION, INTENT(IN)             :: rnew(:)
DOUBLE PRECISION, INTENT(IN)             :: vins(:,:,:)
INTEGER, INTENT(IN)                      :: ncheb
INTEGER, INTENT(IN)                      :: npan_tot
DOUBLE PRECISION, INTENT(IN)             :: rpan_intervall(0:)
INTEGER, INTENT(OUT)                     :: jlk_index(:)
DOUBLE COMPLEX, INTENT(IN OUT)           :: hlk(:,:)
DOUBLE COMPLEX, INTENT(IN OUT)           :: jlk(:,:)
DOUBLE COMPLEX, INTENT(IN OUT)           :: hlk2(:,:)
DOUBLE COMPLEX, INTENT(IN OUT)           :: jlk2(:,:)
DOUBLE COMPLEX, INTENT(IN OUT)           :: gmatprefactor
DOUBLE COMPLEX, INTENT(IN OUT)           :: tmat(:)
INTEGER, INTENT(IN OUT)                  :: use_sratrick
! construct wavefunctions for spherical potentials


! local
INTEGER :: lmsize,lmsize2,nvec
INTEGER :: ivec,lval,ir,ispin,lspin,lsra,i,l1,m1,lm1
INTEGER, allocatable :: jlk_indextemp(:)
DOUBLE COMPLEX, allocatable :: vll0(:,:,:)
DOUBLE COMPLEX, allocatable :: vll(:,:,:)
DOUBLE COMPLEX, allocatable :: rlltemp(:,:,:),slltemp(:,:,:),  &
    hlktemp(:,:),jlktemp(:,:), hlk2temp(:,:),jlk2temp(:,:),  &
    hlknew(:,:),jlknew(:,:)
DOUBLE COMPLEX, allocatable :: tmattemp(:,:)

lmsize=1
IF (nsra == 2) THEN
  lmsize2=2
  nvec=2
ELSE
  lmsize2=1
  nvec=1
END IF
allocate (rlltemp(lmsize2,lmsize,irmdnew))
allocate (slltemp(lmsize2,lmsize,irmdnew))
allocate (hlktemp(nvec,irmdnew))
allocate (jlktemp(nvec,irmdnew))
allocate (hlk2temp(nvec,irmdnew))
allocate (jlk2temp(nvec,irmdnew))
allocate (jlk_indextemp(lmsize2))
allocate (tmattemp(lmsize,lmsize))
allocate (hlknew(nvec*nspin*(lmax+1),irmdnew))
allocate (jlknew(nvec*nspin*(lmax+1),irmdnew))

DO ivec=1,nvec
  jlk_indextemp(ivec)=ivec
END DO
allocate(vll0(lmsize,lmsize,irmdnew))
IF (nsra == 2) THEN
  allocate(vll(2*lmsize,2*lmsize,irmdnew))
ELSE
  allocate(vll(lmsize,lmsize,irmdnew))
END IF
! spin loop
DO ispin=1,nspin
  
  lspin=(lmax+1)*(ispin-1)
  lsra=(lmax+1)*nvec
! each value of l, the Lippmann-Schwinger equation is solved using
! the free-potential wavefunctions and potentials corresponding to l-value
  DO lval=0,lmax
    
    DO ir=1,irmdnew
      vll0(lmsize,lmsize,ir)=vins(ir,1,ispin)-2D0*z/rnew(ir)
    END DO
    
    IF (nsra == 2) THEN
      CALL vllmatsra(vll0,vll,rnew,lmsize,irmdnew,nrmaxd,  &
          e,c,lmax,lval,'Ref=0')
    ELSE
      vll(:,:,:)=vll0(:,:,:)
    END IF
    
    jlktemp(1,:)=jlk(lval+1,:)
    hlktemp(1,:)=hlk(lval+1,:)
    jlk2temp(1,:)=jlk2(lval+1,:)
    hlk2temp(1,:)=hlk2(lval+1,:)
    IF (nsra == 2) THEN
      jlktemp(2,:)=jlk(lmax+lval+2,:)
      hlktemp(2,:)=hlk(lmax+lval+2,:)
      jlk2temp(2,:)=jlk2(lmax+lval+2,:)
      hlk2temp(2,:)=hlk2(lmax+lval+2,:)
    END IF
    CALL rllsll(rpan_intervall,rnew,vll,rlltemp,slltemp,tmattemp,  &
        ncheb,npan_tot,lmsize,lmsize2,nvec,irmdnew,nvec,  &
        jlk_indextemp,hlktemp,jlktemp,hlk2temp,jlk2temp,  &
        gmatprefactor,'1','1',use_sratrick)
    
    DO ir=1,irmdnew
      hlknew(lspin+lval+1,ir)=slltemp(1,1,ir)/rnew(ir)
      jlknew(lspin+lval+1,ir)=rlltemp(1,1,ir)/rnew(ir)
    END DO
    IF (nsra == 2) THEN
      DO ir=1,irmdnew
        hlknew(lspin+lsra+lval+1,ir)=slltemp(2,1,ir)/rnew(ir)
        jlknew(lspin+lsra+lval+1,ir)=rlltemp(2,1,ir)/rnew(ir)
      END DO
    END IF
    tmat(lspin+lval+1)=tmattemp(1,1)
  END DO ! LMAX
END DO ! NSPIN

lm1=1
DO ivec=1,nvec
  DO i=1,2
    DO l1=0,lmax
      DO m1=-l1,l1
        jlk_index(lm1)=l1+(ivec-1)*nspin*(lmax+1)+(i-1)*(lmax+1)+1
        lm1=lm1+1
      END DO
    END DO
  END DO
END DO
DO ir=1,irmdnew
  DO l1=1,nvec*(lmax+1)*nspin
    hlk(l1,ir)=hlknew(l1,ir)
    jlk(l1,ir)=jlknew(l1,ir)
  END DO
END DO
IF (nsra == 2) THEN
  DO ir=1,irmdnew
    DO l1=1,(lmax+1)*nspin
      hlk2(l1,ir)=-hlknew(l1+lmax+1,ir)
      jlk2(l1,ir)=-jlknew(l1+lmax+1,ir)
    END DO
    DO l1=nspin*(lmax+1)+1,nvec*(lmax+1)*nspin
      hlk2(l1,ir)=hlknew(l1-(lmax+1)*nspin,ir)
      jlk2(l1,ir)=jlknew(l1-(lmax+1)*nspin,ir)
    END DO
  END DO
ELSE
  DO ir=1,irmdnew
    DO l1=1,nvec*(lmax+1)*nspin
      hlk2(l1,ir)=-hlknew(l1,ir)
      jlk2(l1,ir)=-jlknew(l1,ir)
    END DO
  END DO
END IF

deallocate (rlltemp)
deallocate (slltemp)
deallocate (hlktemp)
deallocate (jlktemp)
deallocate (hlk2temp)
deallocate (jlk2temp)
deallocate (jlk_indextemp)
deallocate (tmattemp)
deallocate (hlknew)
deallocate (jlknew)
deallocate (vll0)
deallocate (vll)
END SUBROUTINE calcsph

#define hostcode ! comment this out to use the impurity code interface
! choose between interface for impurity and host code (different calling lists)
#ifndef hostcode
      MODULE MOD_RLL_ONLY
        CONTAINS
      SUBROUTINE RLL_ONLY(RPANBOUND,RMESH,VLL,RLL,TLLP, &
                        NCHEB,NPAN,LMSIZE,LMSIZE2,NRMAX, &
                        nvec,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR, &
                        cmoderll,cmodesll,cmodetest,idotime)
#else
      SUBROUTINE RLL_ONLY(RPANBOUND,RMESH,VLL,RLL,TLLP, &
                        NCHEB,NPAN,LMSIZE,LMSIZE2,LBESSEL,NRMAX, &
                        NVEC,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR, &
                        CMODERLL,CMODESLL,USE_SRATRICK1)   !   &
                       ! ALPHAGET) ! LLY
#endif
! ************************************************************************
! radial wave functions by the integral equation method of
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! which has been extended for KKR using non-sperical potentials.
! Further information can be found in 
!
! David Bauer, 
! "Development of a relativistic full-potential first-principles multiple scattering 
! Green function method applied to complex magnetic textures of nano structures 
! at surfaces", PhD Thesis, 2014
!
! http://darwin.bth.rwth-aachen.de/opus3/volltexte/2014/4925/
!
!
!
! ************************************************************************
! This routine solves the following two equations:
!
! ULL(r) = J(r) - PRE * J(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * ULL(r') ) 
!               + PRE * H(r) * int_0^r( dr' r'^2 J2(r') * op(V(r')) * ULL(r') )
!
! where the integral int_0^r() runs from 0 to r
! ************************************************************************
! Potential matrix : VLL(LMSIZE*NVEC,LMSIZE*NVEC)
! LMSIZE = LMMAX (number of LM components) x Number of spin components
! LMSIZE2 = NVEC* LMSIZE 
! NVEC is 2 for a spinor and 1 in case of a non-rel. calculation
! 
! ************************************************************************
! Green function prefacor PRE=GMATPREFACTOR (scalar value)
! tipically \kappa for non-relativistic and M_0 \kappa for SRA 
! 
! ************************************************************************


! ************************************************************************
! The discretization of the Lippmann-Schwinger equation results in a matrix
! equation which is solved in this routine. Further information is given
! in section 5.2.3, page 90 of Bauer, PhD 
!
! Source terms : 
!   right solution:  J, H  (nvec*lmsize,lmsize) or (lmsize,nvec*lmsize)
!    left solution:  J2,H2 (lmsize,nvec*lmsize) or (nvec*lmsize,lmsize)
!
! Example:
! The source term J is for LMSIZE=3 and NVEC=2 given by:
! J =      / jlk(jlk_index(1))                                          \
!          |       0            jlk(jlk_index(2))                       |
!          |       0                   0            jlk(jlk_index(3))   |
!          | jlk(jlk_index(4))                                          |
!          |       0            jlk(jlk_index(5))                       |
!          \       0                   0            jlk(jlk_index(6))   /
!
! first 3 rows are for the large and the last 3 rows for the small component
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
#ifndef hostcode
use mod_beshank                           ! calculates bessel and hankel func.
use mod_chebint                           ! chebyshev integration routines
use mod_config, only: config_testflag     ! reads if testflags are present
use mod_physic_params,only: cvlight       ! speed of light
use sourceterms                           
use mod_chebyshev
#endif
!use mod_timing                            ! timing routine
#ifdef CPP_hybrid
!use omp_lib ! omp functions
#endif
implicit none
      integer :: ncheb                               ! number of chebyshev nodes
      integer :: npan                                ! number of panels
      integer :: lmsize                              ! lm-components * nspin 
      integer :: lmsize2                             ! lmsize * nvec
      integer :: nvec                                ! spinor integer
                                                     ! nvec=1 non-rel, nvec=2 for sra and dirac
      integer :: nrmax                               ! total number of rad. mesh points
#ifdef hostcode
      integer :: LBESSEL, use_sratrick1      !  dimensions etc., needed only for host code interface
#endif

      double complex,parameter:: ci= (0.0d0,1.0d0), &! complex i
                                 cone=(1.0d0,0.0d0),&!         1
                                 czero=(0.0d0,0.0d0) !         0
      ! running indices
      integer ivec, ivec2                            
      integer l1,l2,lm1,lm2,lm3
      integer info,icheb2,icheb,ipan,mn,nm,nplm

      ! source terms
      double complex :: gmatprefactor               ! prefactor of green function
                                                    ! non-rel: = kappa = sqrt e
#ifndef hostcode
      double complex :: hlk(:,:), jlk(:,:), &       ! right sol. source terms
                        hlk2(:,:), jlk2(:,:)        ! left sol. source terms
                                                    ! (tipically bessel and hankel fn)
#else
      DOUBLE COMPLEX :: HLK(LBESSEL,NRMAX), &
                        JLK(LBESSEL,NRMAX), &
                        HLK2(LBESSEL,NRMAX), &
                        JLK2(LBESSEL,NRMAX) 
#endif

#ifndef hostcode
      integer jlk_index(:)                          ! mapping array l = jlk_index(lm)
                                                    ! in: lm-index
                                                    ! corresponding l-index used hlk,..
                                                    ! hlk(l) = jlk_index(lm)
#else
      INTEGER JLK_INDEX(2*LMSIZE)
#endif

      character(len=1) :: cmoderll,cmodesll,cmodetest  ! These define the op(V(r)) in the eqs. above
                                                       ! (comment in the beginning of this subroutine)
                                                       ! cmoderll ="1" : op( )=identity       for reg. solution
                                                       ! cmoderll ="T" : op( )=transpose in L for reg. solution
                                                       ! cmodesll: same for irregular

      double complex ::  rll(lmsize2,lmsize,nrmax), &  ! reg. fredholm sol.
                         tllp(lmsize,lmsize), &        ! t-matrix
                         vll(lmsize*nvec,lmsize*nvec,nrmax) ! potential term in 5.7 
                                                       ! bauer, phd
      double complex,allocatable ::  ull(:,:,:)        ! reg. volterra sol.

      double complex,allocatable ::  &
                     work(:,:), &
                     work2(:,:), &
                     allp(:,:,:),bllp(:,:,:), &                  ! eq. 5.9, 5.10 for reg. sol
                     slv(:,:,:,:),              &                ! a in eq 5.68
                     slv1(:,:,:,:), &                            !****************
!                    slv2(:,:,:,:), &                            ! used for sra trick
                     mrnvy(:,:,:),mrnvz(:,:,:), &                ! ***************
                     mrjvy(:,:,:),mrjvz(:,:,:), &                !    eq. 5.19-5.22
                     yrll(:,:,:),zrll(:,:,:), &                  ! 
                     yrll1(:,:,:),zrll1(:,:,:), &
                     yrll2(:,:,:),zrll2(:,:,:), &
                     vhlr(:,:,:), &                               ! vhlr = h * v (regular sol.) 
                     vjlr(:,:,:)                                  ! vjlr = j * v (regular sol.)
      double complex,allocatable ::  &
                     vhlr_yrll1(:,:), &                           !
                     vhlr_zrll1(:,:), &                           !
                     vjlr_yrll1(:,:), &                           !
                     vjlr_zrll1(:,:), &                           !
                     yrll1temp(:,:), &                            !
                     zrll1temp(:,:)                               !
      double complex,allocatable :: yrf(:,:,:,:), &               ! source terms (different array
                     zrf(:,:,:,:)                                 !               ordering)
      ! chebyshev arrays
      double complex zslc1sum(0:ncheb)
      double precision c1(0:ncheb,0:ncheb),rpanbound(0:npan)
      double precision cslc1(0:ncheb,0:ncheb), & ! Integration matrix from left ( C*S_L*C^-1 in eq. 5.53)
                       csrc1(0:ncheb,0:ncheb), & ! Same from right ( C*S_R*C^-1 in eq. 5.54)
                       tau(0:ncheb,0:npan), &    ! Radial mesh point
                       slc1sum(0:ncheb),rmesh(nrmax),taucslcr

      integer ipiv(0:ncheb,lmsize2)
      integer,allocatable :: ipiv2(:)
!      logical test
!      integer :: ierror
      integer :: use_sratrick
!      integer :: idotime
      integer,parameter  :: directsolv=1
#ifdef hostcode
!      DOUBLE COMPLEX ALPHAGET(LMSIZE,LMSIZE) ! LLY
#endif

#ifdef CPP_hybrid
!     openMP variable --sacin 23/04/2015
!      integer :: thread_id, number_of_openmp_threads,number_of_processor
#endif

      external zgetrf,zgetrs
      intrinsic abs,atan,cos,dimag,exp,max,min,sin,sqrt

! ***********************************************************************
!                                  SRA trick
! ***********************************************************************
! on page 68 of Bauer, PhD, a method is described how to speed up the 
! calculations in case of the SRA. A similar approach is implemented 
! here by using Eq. 4.132 and substituting DV from 4.133, and discretising
! the radial mesh of the Lippmann-Schwinger eq. according to 5.68. 
! The Lippmann-Schwinger equation leads to a matrix inversion 
! problem. The matrix M which needs to be inverted has a special form
! if the SRA approximation is used:
! 
! matrix A ( C 0)     (same as in eq. 5.68)
!          ( B 1)
! (C, B are matricies here)
!
! inverse of A is   (C^-1    0 )
!                   (-B C^-1 1 )
! Thus, it is sufficient to only inverse the matrix C which saves computational
! time. This is refered to as the SRA trick.
! ***********************************************************************
! in future implementation equation 4.134 is supposed to be 
! implemented which should lead to an additional speed-up.
! ***********************************************************************

#ifndef hostcode
if ( config_testflag('nosph') .or. lmsize==1 ) then
  use_sratrick=0
elseif ( .not. config_testflag('nosph') ) then
  use_sratrick=1
else
  stop '[rll] use_sratrick error'
end if
#else
if ( lmsize==1 ) then
  use_sratrick=0
else
  use_sratrick=use_sratrick1
end if
#endif

!#ifdef hostcode
!! turn timing output off if in the host code
!idotime = 0
!#endif
!if (idotime==1) call timing_start('rll')


do ipan = 1,npan
  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
    tau(icheb,ipan) = rmesh(mn)
  end do
end do

call chebint(cslc1,csrc1,slc1sum,c1,ncheb)



if(.not.allocated(ull)) allocate ( ull(lmsize2,lmsize,nrmax) )

if ( use_sratrick==0 ) then
  if(.not.allocated(slv)) allocate ( slv(0:ncheb,lmsize2,0:ncheb,lmsize2) )
elseif ( use_sratrick==1 ) then
  if(.not.allocated(work2)) allocate ( work2((ncheb+1)*lmsize,(ncheb+1)*lmsize), ipiv2((ncheb+1)*lmsize) )
  if(.not.allocated(slv1)) allocate ( slv1(0:ncheb,lmsize,0:ncheb,lmsize) )
! if(.not.allocated(slv2)) allocate ( slv2(0:ncheb,lmsize,0:ncheb,lmsize) )
  if(.not.allocated(yrll1)) allocate ( yrll1(0:ncheb,lmsize,lmsize), zrll1(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yrll2)) allocate ( yrll2(0:ncheb,lmsize,lmsize), zrll2(0:ncheb,lmsize,lmsize) )
else
  stop '[rll] error with testflag sph'
end if

if(.not.allocated(work)) allocate( work(lmsize,lmsize) )
if(.not.allocated(allp)) allocate( allp(lmsize,lmsize,0:npan), bllp(lmsize,lmsize,0:npan) )
if(.not.allocated(mrnvy)) allocate( mrnvy(lmsize,lmsize,npan), mrnvz(lmsize,lmsize,npan) )
if(.not.allocated(mrjvy)) allocate( mrjvy(lmsize,lmsize,npan), mrjvz(lmsize,lmsize,npan) )
if(.not.allocated(yrll)) allocate( yrll(0:ncheb,lmsize2,lmsize), zrll(0:ncheb,lmsize2,lmsize) )
if(.not.allocated(vjlr)) allocate( vjlr(lmsize,lmsize2,0:ncheb), vhlr(lmsize,lmsize2,0:ncheb) )

yrll=(0.0d0,0.0d0)
yrll=(0.0d0,0.0d0)

if(.not.allocated(yrf)) allocate( yrf(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(zrf)) allocate( zrf(lmsize2,lmsize,0:ncheb,npan) )



#ifdef CPP_hybrid
!call omp_set_num_threads(16)
!number_of_openmp_threads = omp_get_num_threads()
!write(*,*) 'number_of_openmp_threads: ', number_of_openmp_threads
!$NOOMP PARALLEL DEFAULT (PRIVATE) &
!$NOOMP&  SHARED(tau,npan,rpanbound,mrnvy,mrnvz,mrjvy,mrjvz,yrf, &
!$NOOMP&  zrf,nvec,lmsize,lmsize2,ncheb,jlk,jlk2,jlk_index,vll,gmatprefactor,hlk,hlk2,cslc1,csrc1,slc1sum, &
!$NOOMP&  cmoderll,cmodesll,cmodetest,use_sratrick, rmesh)

!thread_id = omp_get_thread_num()
#endif

if(.not.allocated(ull)) allocate ( ull(lmsize2,lmsize,nrmax) )

if ( use_sratrick==0 ) then
  if(.not.allocated(slv)) allocate ( slv(0:ncheb,lmsize2,0:ncheb,lmsize2) )
elseif ( use_sratrick==1 ) then
  if(.not.allocated(work2)) allocate ( work2((ncheb+1)*lmsize,(ncheb+1)*lmsize), ipiv2((ncheb+1)*lmsize) )
  if(.not.allocated(slv1)) allocate ( slv1(0:ncheb,lmsize,0:ncheb,lmsize) )
! if(.not.allocated(slv2)) allocate ( slv2(0:ncheb,lmsize,0:ncheb,lmsize) )
  if(.not.allocated(yrll1)) allocate ( yrll1(0:ncheb,lmsize,lmsize), zrll1(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yrll2)) allocate ( yrll2(0:ncheb,lmsize,lmsize), zrll2(0:ncheb,lmsize,lmsize) )
else
  stop '[rll] error with testflag sph'
end if

if(.not.allocated(work)) allocate( work(lmsize,lmsize) )
if(.not.allocated(allp)) allocate( allp(lmsize,lmsize,0:npan), bllp(lmsize,lmsize,0:npan) )
if(.not.allocated(mrnvy)) allocate( mrnvy(lmsize,lmsize,npan), mrnvz(lmsize,lmsize,npan) )
if(.not.allocated(mrjvy)) allocate( mrjvy(lmsize,lmsize,npan), mrjvz(lmsize,lmsize,npan) )
if(.not.allocated(yrll)) allocate( yrll(0:ncheb,lmsize2,lmsize), zrll(0:ncheb,lmsize2,lmsize) )
if(.not.allocated(vjlr)) allocate( vjlr(lmsize,lmsize2,0:ncheb), vhlr(lmsize,lmsize2,0:ncheb) )
if(.not.allocated(vjlr_yrll1)) allocate( vjlr_yrll1(lmsize,lmsize), vhlr_yrll1(lmsize,lmsize) )
if(.not.allocated(vjlr_zrll1)) allocate( vjlr_zrll1(lmsize,lmsize), vhlr_zrll1(lmsize,lmsize) )
if(.not.allocated(yrll1temp)) allocate( yrll1temp(lmsize,lmsize), zrll1temp(lmsize,lmsize) )

yrll=(0.0d0,0.0d0)
yrll=(0.0d0,0.0d0)

if(.not.allocated(yrf)) allocate( yrf(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(zrf)) allocate( zrf(lmsize2,lmsize,0:ncheb,npan) )

!if (idotime==1) call timing_start('local')

! loop over subintervals
#ifdef CPP_hybrid
! openMP pragmas added sachin, parallel region starts earlier to get allocations of arrays right
!$NOOMP DO
#endif
do ipan = 1,npan

!  if (idotime==1) call timing_start('local1')

  ! initialization
  
  vhlr=czero
  vjlr=czero

  if (use_sratrick==0) then

    yrll=czero
    zrll=czero
  else
    yrll1=czero
    zrll1=czero
    yrll2=czero
    zrll2=czero
  end if

!---------------------------------------------------------------------
! 1. prepare VJLR, VNL, VHLR, which appear in the integrands
! TAU(K,IPAN) is used instead of TAU(K,IPAN)**2, which directly gives
! RLL(r) multiplied with r. TAU is the radial mesh.
!
! 2. prepare the source terms YR, ZR
! because of the conventions used by
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! a factor sqrt(E) is included in the source terms
! this factor is removed by the definition of ZSLC1SUM given below
!
!vjlr = \kappa * J * V = \kappa * r * j *V 
!vhlr = \kappa * H * V = \kappa * r * h *V 
!
! i.e. prepare terms kappa*J*DV, kappa*H*DV appearing in 5.11, 5.12.

  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
    if     (cmoderll=='1') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
              vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do 
    elseif (cmoderll=='T') then ! transposed matrix (might not be needed anymore)
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
              vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmoderll=='0') then ! as a test option
              vjlr(:,:,icheb) = czero
              vhlr(:,:,icheb) = czero
    else
      stop '[rll] mode not known'
    end if

    ! calculation of the J (and H) matrix according to equation 5.69 (2nd eq.)
    if ( use_sratrick==0 ) then
      do ivec=1,nvec ! index for large/small component
        do lm1 = 1,lmsize
          l1 = jlk_index( lm1+lmsize*(ivec-1) )
          yrll(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn) 
          zrll(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn) 
        end do
      end do !ivec=1,nvec
    elseif ( use_sratrick==1 ) then
      do lm1 = 1,lmsize
        l1 = jlk_index( lm1+lmsize*(1-1) )
        l2 = jlk_index( lm1+lmsize*(2-1) )
        yrll1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn) 
        zrll1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn) 
        yrll2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l2,mn) 
        zrll2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l2,mn) 
      end do
    end if
  end do ! icheb

  ! calculation of A in 5.68
  if ( use_sratrick==0 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
         taucslcr = tau(icheb,ipan)*cslc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize2
          do ivec=1,nvec
            do lm3 = 1,lmsize
              lm1=lm3+(ivec-1)*lmsize
              l1 = jlk_index(lm1)
              slv(icheb,lm1,icheb2,lm2) = &
              taucslcr*(jlk(l1,mn)*vhlr(lm3,lm2,icheb2) &
                       -hlk(l1,mn)*vjlr(lm3,lm2,icheb2))
            end do
          end do
        end do
      end do
    end do
    do lm1 = 1,lmsize2
      do icheb = 0,ncheb
        slv(icheb,lm1,icheb,lm1) = slv(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do
  elseif  ( use_sratrick==1 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
         taucslcr = tau(icheb,ipan)*cslc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize
!          do ivec=1,1
            do lm3 = 1,lmsize
!             lm1=lm3+(ivec-1)*lmsize
              lm1=lm3
              l1 = jlk_index(lm1)

              ! this is the block to be inverted in SRAtrick. (named C in comment above):
              slv1(icheb,lm1,icheb2,lm2) = &
              taucslcr*(jlk(l1,mn)*vhlr(lm3,lm2,icheb2) &
                       -hlk(l1,mn)*vjlr(lm3,lm2,icheb2))

            end do
!         end do
        end do
      end do
    end do
!   do icheb2 = 0,ncheb
!     do icheb = 0,ncheb
!        taucslcr = tau(icheb,ipan)*cslc1(icheb,icheb2) &
!                   *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
!       mn = ipan*ncheb + ipan - icheb
!       do lm2 = 1,lmsize
!         do ivec=2,2
!           do lm3 = 1,lmsize
!             lm1=lm3+(ivec-1)*lmsize
!             lm1=lm3+lmsize
!             l1 = jlk_index(lm1)

!             slv2(icheb,lm3,icheb2,lm2) = &
!             taucslcr*(jlk(l1,mn)*vhlr(lm3,lm2,icheb2) &
!                      -hlk(l1,mn)*vjlr(lm3,lm2,icheb2))

!           end do
!         end do
!       end do
!     end do
!   end do
    do lm1 = 1,lmsize
      do icheb = 0,ncheb
        slv1(icheb,lm1,icheb,lm1) = slv1(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do

  else
    stop '[rll] error in inversion'
  end if

!  if (idotime==1) call timing_pause('local1')
!  if (idotime==1) call timing_start('local2')

!-------------------------------------------------------
! determine the local solutions
! solve the equations SLV*YRLL=S and SLV*ZRLL=C 
!                 and SRV*YILL=C and SRV*ZILL=S
! i.e., solve system A*U=J, see eq. 5.68.

  if ( use_sratrick==0 ) then
    nplm = (ncheb+1)*lmsize2

    if (cmoderll/='0') then
!      if (idotime==1) call timing_start('inversion')
      call zgetrf(nplm,nplm,slv,nplm,ipiv,info)
!      if (idotime==1) call timing_stop('inversion','test')
      if (info/=0) stop 'rll: zgetrf'
      call zgetrs('n',nplm,lmsize,slv,nplm,ipiv,yrll,nplm,info)
      call zgetrs('n',nplm,lmsize,slv,nplm,ipiv,zrll,nplm,info)
    end if
  elseif ( use_sratrick==1 ) then
    nplm = (ncheb+1)*lmsize

    call zgetrf(nplm,nplm,slv1,nplm,ipiv,info)
    call zgetrs('n',nplm,lmsize,slv1,nplm,ipiv,yrll1,nplm,info)
    call zgetrs('n',nplm,lmsize,slv1,nplm,ipiv,zrll1,nplm,info)

!   call zgemm('n','n',nplm,lmsize,nplm,-cone,slv2, &
!       nplm,yrll1,nplm,cone,yrll2,nplm)

!   call zgemm('n','n',nplm,lmsize,nplm,-cone,slv2, &
!       nplm,zrll1,nplm,cone,zrll2,nplm)

    do icheb2 = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrll1temp(lm1,lm2) = yrll1(icheb2,lm1,lm2)
          zrll1temp(lm1,lm2) = zrll1(icheb2,lm1,lm2)
        end do
      end do
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhlr(1,1,icheb2), &
        lmsize,yrll1temp,lmsize,czero,vhlr_yrll1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhlr(1,1,icheb2), &
        lmsize,zrll1temp,lmsize,czero,vhlr_zrll1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjlr(1,1,icheb2), &
        lmsize,yrll1temp,lmsize,czero,vjlr_yrll1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjlr(1,1,icheb2), &
        lmsize,zrll1temp,lmsize,czero,vjlr_zrll1,lmsize)

      do icheb = 0,ncheb
         taucslcr = - tau(icheb,ipan)*cslc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize
            do lm3 = 1,lmsize
              lm1=lm3+lmsize
              l1 = jlk_index(lm1)

              yrll2(icheb,lm3,lm2) = &
              yrll2(icheb,lm3,lm2) + &
              taucslcr*(jlk(l1,mn)*vhlr_yrll1(lm3,lm2) &
                       -hlk(l1,mn)*vjlr_yrll1(lm3,lm2))

              zrll2(icheb,lm3,lm2) = &
              zrll2(icheb,lm3,lm2) + &
              taucslcr*(jlk(l1,mn)*vhlr_zrll1(lm3,lm2) &
                       -hlk(l1,mn)*vjlr_zrll1(lm3,lm2))

            end do
        end do
      end do
    end do

  else
    stop '[rll] error in inversion'
  end if

  ! Reorient indices for later use
  if ( use_sratrick==0 ) then
    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize2
          yrf(lm1,lm2,icheb,ipan) = yrll(icheb,lm1,lm2)
          zrf(lm1,lm2,icheb,ipan) = zrll(icheb,lm1,lm2)
        end do
      end do
    end do

  elseif ( use_sratrick==1 ) then

    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrf(lm1,lm2,icheb,ipan) = yrll1(icheb,lm1,lm2)
          zrf(lm1,lm2,icheb,ipan) = zrll1(icheb,lm1,lm2)
        end do
      end do
    end do

    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrf(lm1+lmsize,lm2,icheb,ipan) = yrll2(icheb,lm1,lm2)
          zrf(lm1+lmsize,lm2,icheb,ipan) = zrll2(icheb,lm1,lm2)
        end do
      end do
    end do

  end if

!  if (idotime==1) call timing_pause('local2')
!  if (idotime==1) call timing_start('local3')

  ! Calculation of eq. 5.19-5.22

  do icheb = 0,ncheb
    zslc1sum(icheb) = slc1sum(icheb) * (rpanbound(ipan)-rpanbound(ipan-1))/ (2.d0)
  end do
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhlr(1,1,0), &
        lmsize,yrf(1,1,0,ipan),lmsize2,czero,mrnvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjlr(1,1,0), &
        lmsize,yrf(1,1,0,ipan),lmsize2,czero,mrjvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhlr(1,1,0), &
        lmsize,zrf(1,1,0,ipan),lmsize2,czero,mrnvz(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjlr(1,1,0), &
        lmsize,zrf(1,1,0,ipan),lmsize2,czero,mrjvz(1,1,ipan),lmsize)
  do icheb = 1,ncheb
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhlr(1,1,icheb), &
          lmsize,yrf(1,1,icheb,ipan),lmsize2,cone,mrnvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjlr(1,1,icheb), &
          lmsize,yrf(1,1,icheb,ipan),lmsize2,cone,mrjvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhlr(1,1,icheb), &
          lmsize,zrf(1,1,icheb,ipan),lmsize2,cone,mrnvz(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjlr(1,1,icheb), &
          lmsize,zrf(1,1,icheb,ipan),lmsize2,cone,mrjvz(1,1,ipan),lmsize)
  end do
!  if (idotime==1) call timing_pause('local3')

end do !ipan
#ifdef CPP_hybrid
!$NOOMP END DO
!$NOOMP END PARALLEL
#endif
! end the big loop over the subintervals



!if (idotime==1) call timing_stop('local')
!if (idotime==1) call timing_start('afterlocal')

! ***********************************************************************
! calculate A(M), B(M), C(M), D(M)
! according to 5.17-5.18 (regular solution) of Bauer PhD
! C,D are calculated accordingly for the irregular solution
! (starting condition: A(0) = 1, B(0) = 0, C(MMAX) = 0 and D(MMAX) = 1)
! ***********************************************************************

! regular 
do lm2 = 1,lmsize
  do lm1 = 1,lmsize
    bllp(lm1,lm2,0) = czero
    allp(lm1,lm2,0) = czero
  end do
end do

do lm1 = 1,lmsize
  allp(lm1,lm1,0) = cone
end do

do ipan = 1,npan
  call zcopy(lmsize*lmsize,allp(1,1,ipan-1),1,allp(1,1,ipan),1)
  call zcopy(lmsize*lmsize,bllp(1,1,ipan-1),1,bllp(1,1,ipan),1)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvy(1,1,ipan), &
          lmsize,allp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvz(1,1,ipan), &
          lmsize,bllp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mrjvy(1,1,ipan), &
          lmsize,allp(1,1,ipan-1),lmsize,cone,bllp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mrjvz(1,1,ipan), &
          lMSIZE,BLLP(1,1,IPAN-1),LMSIZE,CONE,BLLP(1,1,IPAN),LMSIZE)
end do

! ***********************************************************************
! determine the regular solution ull by using 5.14
! ***********************************************************************
do ipan = 1,npan
  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,yrf(1,1,icheb,ipan), &
          lmsize2,allp(1,1,ipan-1),lmsize,czero,ull(1,1,mn),lmsize2)
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,zrf(1,1,icheb,ipan), &
          lmsize2,bllp(1,1,ipan-1),lmsize,cone,ull(1,1,mn),lmsize2)
  end do
end do

!if (idotime==1) call timing_stop('afterlocal')
!if (idotime==1) call timing_start('endstuff')

! ***********************************************************************
! next part converts the volterra solution u of equation (5.7) to
! the fredholm solution r by employing eq. 4.122 and 4.120 of bauer, phd
! and the t-matrix is calculated
! ***********************************************************************

call zgetrf(lmsize,lmsize,allp(1,1,npan),lmsize,ipiv,info)                     !invert alpha
call zgetri(lmsize,allp(1,1,npan),lmsize,ipiv,work,lmsize*lmsize,info)         !invert alpha -> transformation matrix rll=alpha^-1*rll
#ifdef hostcode
! get alpha matrix
!      DO LM1=1,LMSIZE                          ! LLY
!       DO LM2=1,LMSIZE                         ! LLY
!        ALPHAGET(LM1,LM2)=ALLP(LM1,LM2,NPAN)   ! LLY
!       ENDDO                                   ! LLY
!      ENDDO                                    ! LLY
#endif
! calculation of the t-matrix 
call zgemm('n','n',lmsize,lmsize,lmsize,cone/gmatprefactor,bllp(1,1,npan), &   ! calc t-matrix tll = bll*alpha^-1 
            lmsize,allp(1,1,npan),lmsize,czero,tllp,lmsize)

do nm = 1,nrmax
call zgemm('n','n',lmsize2,lmsize,lmsize,cone,ull(1,1,nm), &
            lmsize2,allp(1,1,npan),lmsize,czero,rll(1,1,nm),lmsize2)
end do

!if (idotime==1) call timing_stop('endstuff')
!if (idotime==1) call timing_start('checknan')
!if (idotime==1) call timing_stop('checknan')
!if (idotime==1) call timing_stop('local1')
!if (idotime==1) call timing_stop('local2')
!if (idotime==1) call timing_stop('local3')
!if (idotime==1) call timing_stop('rll')

if ( use_sratrick==0 ) then
  if(allocated(slv)) deallocate ( slv )
elseif ( use_sratrick==1 ) then
  if(allocated(work2)) deallocate ( work2, ipiv2 )
  if(allocated(slv1)) deallocate ( slv1 )
! if(allocated(slv2)) deallocate ( slv2 )
  if(allocated(yrll1)) deallocate ( yrll1, zrll1 )
  if(allocated(yrll2)) deallocate ( yrll2, zrll2 )
end if

if(allocated(work)) deallocate( work )
if(allocated(allp)) deallocate( allp, bllp )
if(allocated(mrnvy)) deallocate( mrnvy, mrnvz )
if(allocated(mrjvy)) deallocate( mrjvy, mrjvz )
if(allocated(yrll)) deallocate( yrll, zrll )
if(allocated(vjlr)) deallocate( vjlr, vhlr )
if(allocated(vjlr_yrll1)) deallocate( vjlr_yrll1, vhlr_yrll1 )
if(allocated(vjlr_zrll1)) deallocate( vjlr_zrll1, vhlr_zrll1 )
if(allocated(yrll1temp)) deallocate( yrll1temp, zrll1temp )

if(allocated(yrf)) deallocate( yrf )
if(allocated(zrf)) deallocate( zrf )

end subroutine

#ifndef hostcode
END MODULE MOD_RLL_ONLY
#endif

#define hostcode ! comment this out to use the impurity code interface
! choose between interface for impurity and host code (different calling lists)
#ifndef hostcode
      MODULE MOD_RLLSLL
        CONTAINS
      SUBROUTINE RLLSLL(RPANBOUND,RMESH,VLL,RLL,SLL,TLLP, &
                        NCHEB,NPAN,LMSIZE,LMSIZE2,NRMAX, &
                        nvec,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR, &
                        cmoderll,cmodesll,cmodetest,idotime)
#else
      SUBROUTINE RLLSLL(RPANBOUND,RMESH,VLL,RLL,SLL,TLLP, &
                        NCHEB,NPAN,LMSIZE,LMSIZE2,LBESSEL,NRMAX, &
                        NVEC,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR, &
                        CMODERLL,CMODESLL,USE_SRATRICK1)   !   &
                       ! ALPHAGET) ! LLY
#endif
! ************************************************************************
! radial wave functions by the integral equation method of
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! which has been extended for KKR using non-sperical potentials.
! Further information can be found in 
!
! David Bauer, 
! "Development of a relativistic full-potential first-principles multiple scattering 
! Green function method applied to complex magnetic textures of nano structures 
! at surfaces", PhD Thesis, 2014
!
! http://darwin.bth.rwth-aachen.de/opus3/volltexte/2014/4925/
!
!
!
! ************************************************************************
! This routine solves the following two equations:
!
! ULL(r) = J(r) - PRE * J(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * ULL(r') ) 
!               + PRE * H(r) * int_0^r( dr' r'^2 J2(r') * op(V(r')) * ULL(r') )
!
! SLL(r) = H(r) - PRE * H(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * RLL(r') ) 
!               + PRE * J(r) * int_0^r( dr' r'^2 H2(r') * op(V(r')) * SLL(r') )
!
! where the integral int_0^r() runs from 0 to r
! ************************************************************************
! Potential matrix : VLL(LMSIZE*NVEC,LMSIZE*NVEC)
! LMSIZE = LMMAX (number of LM components) x Number of spin components
! LMSIZE2 = NVEC* LMSIZE 
! NVEC is 2 for a spinor and 1 in case of a non-rel. calculation
! 
! ************************************************************************
! Green function prefacor PRE=GMATPREFACTOR (scalar value)
! tipically \kappa for non-relativistic and M_0 \kappa for SRA 
! 
! ************************************************************************


! ************************************************************************
! The discretization of the Lippmann-Schwinger equation results in a matrix
! equation which is solved in this routine. Further information is given
! in section 5.2.3, page 90 of Bauer, PhD 
!
! Source terms : 
!   right solution:  J, H  (nvec*lmsize,lmsize) or (lmsize,nvec*lmsize)
!    left solution:  J2,H2 (lmsize,nvec*lmsize) or (nvec*lmsize,lmsize)
!
! Example:
! The source term J is for LMSIZE=3 and NVEC=2 given by:
! J =      / jlk(jlk_index(1))                                          \
!          |       0            jlk(jlk_index(2))                       |
!          |       0                   0            jlk(jlk_index(3))   |
!          | jlk(jlk_index(4))                                          |
!          |       0            jlk(jlk_index(5))                       |
!          \       0                   0            jlk(jlk_index(6))   /
!
! first 3 rows are for the large and the last 3 rows for the small component
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
#ifndef hostcode
use mod_beshank                           ! calculates bessel and hankel func.
use mod_chebint                           ! chebyshev integration routines
use mod_config, only: config_testflag     ! reads if testflags are present
use mod_physic_params,only: cvlight       ! speed of light
use sourceterms                           
use mod_chebyshev
#endif
!use mod_timing                            ! timing routine
#ifdef CPP_hybrid
!use omp_lib ! omp functions
#endif
implicit none
      integer :: ncheb                               ! number of chebyshev nodes
      integer :: npan                                ! number of panels
      integer :: lmsize                              ! lm-components * nspin 
      integer :: lmsize2                             ! lmsize * nvec
      integer :: nvec                                ! spinor integer
                                                     ! nvec=1 non-rel, nvec=2 for sra and dirac
      integer :: nrmax                               ! total number of rad. mesh points
#ifdef hostcode
      integer :: LBESSEL, use_sratrick1      !  dimensions etc., needed only for host code interface
#endif

      double complex,parameter:: ci= (0.0d0,1.0d0), &! complex i
                                 cone=(1.0d0,0.0d0),&!         1
                                 czero=(0.0d0,0.0d0) !         0
      ! running indices
      integer ivec, ivec2                            
      integer l1,l2,lm1,lm2,lm3
      integer info,icheb2,icheb,ipan,mn,nm,nplm

      ! source terms
      double complex :: gmatprefactor               ! prefactor of green function
                                                    ! non-rel: = kappa = sqrt e
#ifndef hostcode
      double complex :: hlk(:,:), jlk(:,:), &       ! right sol. source terms
                        hlk2(:,:), jlk2(:,:)        ! left sol. source terms
                                                    ! (tipically bessel and hankel fn)
#else
      DOUBLE COMPLEX :: HLK(LBESSEL,NRMAX), &
                        JLK(LBESSEL,NRMAX), &
                        HLK2(LBESSEL,NRMAX), &
                        JLK2(LBESSEL,NRMAX) 
#endif

#ifndef hostcode
      integer jlk_index(:)                          ! mapping array l = jlk_index(lm)
                                                    ! in: lm-index
                                                    ! corresponding l-index used hlk,..
                                                    ! hlk(l) = jlk_index(lm)
#else
      INTEGER JLK_INDEX(2*LMSIZE)
#endif

      character(len=1) :: cmoderll,cmodesll,cmodetest  ! These define the op(V(r)) in the eqs. above
                                                       ! (comment in the beginning of this subroutine)
                                                       ! cmoderll ="1" : op( )=identity       for reg. solution
                                                       ! cmoderll ="T" : op( )=transpose in L for reg. solution
                                                       ! cmodesll: same for irregular

      double complex ::  sll(lmsize2,lmsize,nrmax), &  ! irr. volterra sol.
                         rll(lmsize2,lmsize,nrmax), &  ! reg. fredholm sol.
                         tllp(lmsize,lmsize), &        ! t-matrix
                         vll(lmsize*nvec,lmsize*nvec,nrmax) ! potential term in 5.7 
                                                       ! bauer, phd
      double complex,allocatable ::  ull(:,:,:)        ! reg. volterra sol.

      double complex,allocatable ::  &
                     work(:,:), &
                     work2(:,:), &
                     allp(:,:,:),bllp(:,:,:), &                  ! eq. 5.9, 5.10 for reg. sol
                     cllp(:,:,:),dllp(:,:,:), &                  ! same for the irr. sol
                     slv(:,:,:,:),srv(:,:,:,:), &                ! a in eq 5.68
                     slv1(:,:,:,:),srv1(:,:,:,:), &              !****************
!                    slv2(:,:,:,:),srv2(:,:,:,:), &              ! used for sra trick
                     mrnvy(:,:,:),mrnvz(:,:,:), &                ! ***************
                     mrjvy(:,:,:),mrjvz(:,:,:), &                !    eq. 5.19-5.22
                     mihvy(:,:,:),mihvz(:,:,:), &                !
                     mijvy(:,:,:),mijvz(:,:,:), &                ! ***************
                     yill(:,:,:),zill(:,:,:), &                  ! source terms  (i:irreg., r: regular)
                     yrll(:,:,:),zrll(:,:,:), &                  ! 
                     yill1(:,:,:),zill1(:,:,:), &                ! source terms in case of sratrick
                     yrll1(:,:,:),zrll1(:,:,:), &
                     yill2(:,:,:),zill2(:,:,:), &
                     yrll2(:,:,:),zrll2(:,:,:), &
                     vhlr(:,:,:), &                               ! vhlr = h * v (regular sol.) 
                     vjlr(:,:,:), &                               ! vjlr = j * v (regular sol.)
                     vhli(:,:,:), &                               ! vhli = h * v (irregular sol.)
                     vjli(:,:,:)                                  ! vjli = j * v (irregular sol.)
      double complex,allocatable ::  &
                     vhlr_yrll1(:,:), &                           !
                     vhlr_zrll1(:,:), &                           !
                     vjlr_yrll1(:,:), &                           !
                     vjlr_zrll1(:,:), &                           !
                     yrll1temp(:,:), &                            !
                     zrll1temp(:,:), &                            !
                     yill1temp(:,:), &                            !
                     zill1temp(:,:), &                            ! 
                     vhli_yill1(:,:), &                           !
                     vhli_zill1(:,:), &                           !
                     vjli_yill1(:,:), &                           !
                     vjli_zill1(:,:)
      double complex,allocatable :: yif(:,:,:,:), &               ! source terms (different array
                     yrf(:,:,:,:), &                              !               ordering)
                     zif(:,:,:,:), &
                     zrf(:,:,:,:)
      ! chebyshev arrays
      double complex zslc1sum(0:ncheb)
      double precision c1(0:ncheb,0:ncheb),rpanbound(0:npan)
      double precision cslc1(0:ncheb,0:ncheb), & ! Integration matrix from left ( C*S_L*C^-1 in eq. 5.53)
                       csrc1(0:ncheb,0:ncheb), & ! Same from right ( C*S_R*C^-1 in eq. 5.54)
                       tau(0:ncheb,0:npan), &    ! Radial mesh point
                       slc1sum(0:ncheb),rmesh(nrmax),taucslcr,taucsrcr

      integer ipiv(0:ncheb,lmsize2)
      integer,allocatable :: ipiv2(:)
!      logical test
!      integer :: ierror
      integer :: use_sratrick
!      integer :: idotime
      integer,parameter  :: directsolv=1
#ifdef hostcode
!      DOUBLE COMPLEX ALPHAGET(LMSIZE,LMSIZE) ! LLY
#endif

#ifdef CPP_hybrid
!     openMP variable --sacin 23/04/2015
!      integer :: thread_id, number_of_openmp_threads,number_of_processor
#endif

      external zgetrf,zgetrs
      intrinsic abs,atan,cos,dimag,exp,max,min,sin,sqrt

! ***********************************************************************
!                                  SRA trick
! ***********************************************************************
! on page 68 of Bauer, PhD, a method is described how to speed up the 
! calculations in case of the SRA. A similar approach is implemented 
! here by using Eq. 4.132 and substituting DV from 4.133, and discretising
! the radial mesh of the Lippmann-Schwinger eq. according to 5.68. 
! The Lippmann-Schwinger equation leads to a matrix inversion 
! problem. The matrix M which needs to be inverted has a special form
! if the SRA approximation is used:
! 
! matrix A ( C 0)     (same as in eq. 5.68)
!          ( B 1)
! (C, B are matricies here)
!
! inverse of A is   (C^-1    0 )
!                   (-B C^-1 1 )
! Thus, it is sufficient to only inverse the matrix C which saves computational
! time. This is refered to as the SRA trick.
! ***********************************************************************
! in future implementation equation 4.134 is supposed to be 
! implemented which should lead to an additional speed-up.
! ***********************************************************************

#ifndef hostcode
if ( config_testflag('nosph') .or. lmsize==1 ) then
  use_sratrick=0
elseif ( .not. config_testflag('nosph') ) then
  use_sratrick=1
else
  stop '[rllsll] use_sratrick error'
end if
#else
if ( lmsize==1 ) then
  use_sratrick=0
else
  use_sratrick=use_sratrick1
end if
#endif

!#ifdef hostcode
!! turn timing output off if in the host code
!idotime = 0
!#endif
!if (idotime==1) call timing_start('rllsll')


do ipan = 1,npan
  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
    tau(icheb,ipan) = rmesh(mn)
  end do
end do

call chebint(cslc1,csrc1,slc1sum,c1,ncheb)



if(.not.allocated(ull)) allocate ( ull(lmsize2,lmsize,nrmax) )

if ( use_sratrick==0 ) then
  if(.not.allocated(slv)) allocate ( slv(0:ncheb,lmsize2,0:ncheb,lmsize2),srv(0:ncheb,lmsize2,0:ncheb,lmsize2) )
elseif ( use_sratrick==1 ) then
  if(.not.allocated(work2)) allocate ( work2((ncheb+1)*lmsize,(ncheb+1)*lmsize), ipiv2((ncheb+1)*lmsize) )
  if(.not.allocated(slv1)) allocate ( slv1(0:ncheb,lmsize,0:ncheb,lmsize), srv1(0:ncheb,lmsize,0:ncheb,lmsize) )
! if(.not.allocated(slv2)) allocate ( slv2(0:ncheb,lmsize,0:ncheb,lmsize), srv2(0:ncheb,lmsize,0:ncheb,lmsize) )
  if(.not.allocated(yill1)) allocate ( yill1(0:ncheb,lmsize,lmsize), zill1(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yrll1)) allocate ( yrll1(0:ncheb,lmsize,lmsize), zrll1(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yill2)) allocate ( yill2(0:ncheb,lmsize,lmsize), zill2(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yrll2)) allocate ( yrll2(0:ncheb,lmsize,lmsize), zrll2(0:ncheb,lmsize,lmsize) )
else
  stop '[rllsll] error with testflag sph'
end if

if(.not.allocated(work)) allocate( work(lmsize,lmsize) )
if(.not.allocated(allp)) allocate( allp(lmsize,lmsize,0:npan), bllp(lmsize,lmsize,0:npan) )
if(.not.allocated(cllp)) allocate( cllp(lmsize,lmsize,0:npan), dllp(lmsize,lmsize,0:npan) )
if(.not.allocated(mrnvy)) allocate( mrnvy(lmsize,lmsize,npan), mrnvz(lmsize,lmsize,npan) )
if(.not.allocated(mrjvy)) allocate( mrjvy(lmsize,lmsize,npan), mrjvz(lmsize,lmsize,npan) )
if(.not.allocated(mihvy)) allocate( mihvy(lmsize,lmsize,npan), mihvz(lmsize,lmsize,npan) )
if(.not.allocated(mijvy)) allocate( mijvy(lmsize,lmsize,npan), mijvz(lmsize,lmsize,npan) )
if(.not.allocated(yill)) allocate( yill(0:ncheb,lmsize2,lmsize), zill(0:ncheb,lmsize2,lmsize) )
if(.not.allocated(yrll)) allocate( yrll(0:ncheb,lmsize2,lmsize), zrll(0:ncheb,lmsize2,lmsize) )
if(.not.allocated(vjlr)) allocate( vjlr(lmsize,lmsize2,0:ncheb), vhlr(lmsize,lmsize2,0:ncheb) )
if(.not.allocated(vjli)) allocate( vjli(lmsize,lmsize2,0:ncheb), vhli(lmsize,lmsize2,0:ncheb) )
if(.not.allocated(vjlr_yrll1)) allocate( vjlr_yrll1(lmsize,lmsize), vhlr_yrll1(lmsize,lmsize) )
if(.not.allocated(vjlr_zrll1)) allocate( vjlr_zrll1(lmsize,lmsize), vhlr_zrll1(lmsize,lmsize) )
if(.not.allocated(yrll1temp)) allocate( yrll1temp(lmsize,lmsize), zrll1temp(lmsize,lmsize) )
if(.not.allocated(vjli_yill1)) allocate( vjli_yill1(lmsize,lmsize), vhli_yill1(lmsize,lmsize) )
if(.not.allocated(vjli_zill1)) allocate( vjli_zill1(lmsize,lmsize), vhli_zill1(lmsize,lmsize) )
if(.not.allocated(yill1temp)) allocate( yill1temp(lmsize,lmsize), zill1temp(lmsize,lmsize) )

yrll=(0.0d0,0.0d0)
zill=(0.0d0,0.0d0)
yrll=(0.0d0,0.0d0)
zill=(0.0d0,0.0d0)

if(.not.allocated(yif)) allocate( yif(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(yrf)) allocate( yrf(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(zif)) allocate( zif(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(zrf)) allocate( zrf(lmsize2,lmsize,0:ncheb,npan) )



#ifdef CPP_hybrid
!call omp_set_num_threads(16)
!number_of_openmp_threads = omp_get_num_threads()
!write(*,*) 'number_of_openmp_threads: ', number_of_openmp_threads
!$NOOMP PARALLEL DEFAULT (PRIVATE) &
!$NOOMP&  SHARED(tau,npan,rpanbound,mrnvy,mrnvz,mrjvy,mrjvz,mihvy,mihvz,mijvy,mijvz,yif,yrf, &
!$NOOMP&  zif,zrf,nvec,lmsize,lmsize2,ncheb,jlk,jlk2,jlk_index,vll,gmatprefactor,hlk,hlk2,cslc1,csrc1,slc1sum, &
!$NOOMP&  cmoderll,cmodesll,cmodetest,use_sratrick, rmesh)

!thread_id = omp_get_thread_num()
#endif

if(.not.allocated(ull)) allocate ( ull(lmsize2,lmsize,nrmax) )

if ( use_sratrick==0 ) then
  if(.not.allocated(slv)) allocate ( slv(0:ncheb,lmsize2,0:ncheb,lmsize2),srv(0:ncheb,lmsize2,0:ncheb,lmsize2) )
elseif ( use_sratrick==1 ) then
  if(.not.allocated(work2)) allocate ( work2((ncheb+1)*lmsize,(ncheb+1)*lmsize), ipiv2((ncheb+1)*lmsize) )
  if(.not.allocated(slv1)) allocate ( slv1(0:ncheb,lmsize,0:ncheb,lmsize), srv1(0:ncheb,lmsize,0:ncheb,lmsize) )
! if(.not.allocated(slv2)) allocate ( slv2(0:ncheb,lmsize,0:ncheb,lmsize), srv2(0:ncheb,lmsize,0:ncheb,lmsize) )
  if(.not.allocated(yill1)) allocate ( yill1(0:ncheb,lmsize,lmsize), zill1(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yrll1)) allocate ( yrll1(0:ncheb,lmsize,lmsize), zrll1(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yill2)) allocate ( yill2(0:ncheb,lmsize,lmsize), zill2(0:ncheb,lmsize,lmsize) )
  if(.not.allocated(yrll2)) allocate ( yrll2(0:ncheb,lmsize,lmsize), zrll2(0:ncheb,lmsize,lmsize) )
else
  stop '[rllsll] error with testflag sph'
end if

if(.not.allocated(work)) allocate( work(lmsize,lmsize) )
if(.not.allocated(allp)) allocate( allp(lmsize,lmsize,0:npan), bllp(lmsize,lmsize,0:npan) )
if(.not.allocated(cllp)) allocate( cllp(lmsize,lmsize,0:npan), dllp(lmsize,lmsize,0:npan) )
if(.not.allocated(mrnvy)) allocate( mrnvy(lmsize,lmsize,npan), mrnvz(lmsize,lmsize,npan) )
if(.not.allocated(mrjvy)) allocate( mrjvy(lmsize,lmsize,npan), mrjvz(lmsize,lmsize,npan) )
if(.not.allocated(mihvy)) allocate( mihvy(lmsize,lmsize,npan), mihvz(lmsize,lmsize,npan) )
if(.not.allocated(mijvy)) allocate( mijvy(lmsize,lmsize,npan), mijvz(lmsize,lmsize,npan) )
if(.not.allocated(yill)) allocate( yill(0:ncheb,lmsize2,lmsize), zill(0:ncheb,lmsize2,lmsize) )
if(.not.allocated(yrll)) allocate( yrll(0:ncheb,lmsize2,lmsize), zrll(0:ncheb,lmsize2,lmsize) )
if(.not.allocated(vjlr)) allocate( vjlr(lmsize,lmsize2,0:ncheb), vhlr(lmsize,lmsize2,0:ncheb) )
if(.not.allocated(vjli)) allocate( vjli(lmsize,lmsize2,0:ncheb), vhli(lmsize,lmsize2,0:ncheb) )

yrll=(0.0d0,0.0d0)
zill=(0.0d0,0.0d0)
yrll=(0.0d0,0.0d0)
zill=(0.0d0,0.0d0)

if(.not.allocated(yif)) allocate( yif(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(yrf)) allocate( yrf(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(zif)) allocate( zif(lmsize2,lmsize,0:ncheb,npan) )
if(.not.allocated(zrf)) allocate( zrf(lmsize2,lmsize,0:ncheb,npan) )

!if (idotime==1) call timing_start('local')

! loop over subintervals
#ifdef CPP_hybrid
! openMP pragmas added sachin, parallel region starts earlier to get allocations of arrays right
!$NOOMP DO
#endif
do ipan = 1,npan

!  if (idotime==1) call timing_start('local1')

  ! initialization
  
  vhlr=czero
  vjlr=czero
  vhli=czero
  vjli=czero

  if (use_sratrick==0) then

    yrll=czero
    zrll=czero
    yill=czero
    zill=czero
  else
    yrll1=czero
    zrll1=czero
    yill1=czero
    zill1=czero
    yrll2=czero
    zrll2=czero
    yill2=czero
    zill2=czero
  end if

!---------------------------------------------------------------------
! 1. prepare VJLR, VNL, VHLR, which appear in the integrands
! TAU(K,IPAN) is used instead of TAU(K,IPAN)**2, which directly gives
! RLL(r) and SLL(r) multiplied with r. TAU is the radial mesh.
!
! 2. prepare the source terms YR, ZR, YI, ZI 
! because of the conventions used by
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! a factor sqrt(E) is included in the source terms
! this factor is removed by the definition of ZSLC1SUM given below
!
!vjlr = \kappa * J * V = \kappa * r * j *V 
!vhlr = \kappa * H * V = \kappa * r * h *V 
!
! i.e. prepare terms kappa*J*DV, kappa*H*DV appearing in 5.11, 5.12.

  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
    if     (cmoderll=='1') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
              vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do 
    elseif (cmoderll=='T') then ! transposed matrix (might not be needed anymore)
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vjlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
              vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) = vhlr(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmoderll=='0') then ! as a test option
              vjlr(:,:,icheb) = czero
              vhlr(:,:,icheb) = czero
    else
      stop '[rllsll] mode not known'
    end if

    if     (cmodesll=='1') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjli(lm1,lm2+lmsize*(ivec2-1),icheb) = vjli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
              vhli(lm1,lm2+lmsize*(ivec2-1),icheb) = vhli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmodesll=='T') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjli(lm1,lm2+lmsize*(ivec2-1),icheb) = vjli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*jlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
              vhli(lm1,lm2+lmsize*(ivec2-1),icheb) = vhli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gmatprefactor*tau(icheb,ipan)*hlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmodesll=='0') then
              vjli(:,:,icheb) = czero
              vhli(:,:,icheb) = czero
    else
      stop '[rllsll] mode not known'
    end if

    ! calculation of the J (and H) matrix according to equation 5.69 (2nd eq.)
    if ( use_sratrick==0 ) then
      do ivec=1,nvec ! index for large/small component
        do lm1 = 1,lmsize
          l1 = jlk_index( lm1+lmsize*(ivec-1) )
          yrll(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn) 
          zrll(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn) 
          yill(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn)
          zill(icheb,lm1+lmsize*(ivec-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn)
        end do
      end do !ivec=1,nvec
    elseif ( use_sratrick==1 ) then
      do lm1 = 1,lmsize
        l1 = jlk_index( lm1+lmsize*(1-1) )
        l2 = jlk_index( lm1+lmsize*(2-1) )
        yrll1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn) 
        zrll1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn) 
        yill1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l1,mn)
        zill1(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l1,mn)
        yrll2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l2,mn) 
        zrll2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l2,mn) 
        yill2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*hlk(l2,mn)
        zill2(icheb,lm1+lmsize*(1-1),lm1) =  tau(icheb,ipan)*jlk(l2,mn)
      end do
    end if
  end do ! icheb

  ! calculation of A in 5.68
  if ( use_sratrick==0 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
         taucslcr = tau(icheb,ipan)*cslc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
         taucsrcr = tau(icheb,ipan)*csrc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0 
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize2
          do ivec=1,nvec
            do lm3 = 1,lmsize
              lm1=lm3+(ivec-1)*lmsize
              l1 = jlk_index(lm1)
              slv(icheb,lm1,icheb2,lm2) = &
              taucslcr*(jlk(l1,mn)*vhlr(lm3,lm2,icheb2) &
                       -hlk(l1,mn)*vjlr(lm3,lm2,icheb2))
              srv(icheb,lm1,icheb2,lm2) = &
              taucsrcr*(-jlk(l1,mn)*vhli(lm3,lm2,icheb2) &
                        +hlk(l1,mn)*vjli(lm3,lm2,icheb2))
!             slv(icheb,lm1,icheb2,lm2) = &
!           ( tau(icheb,ipan)*jlk(l1,mn)*cslc1(icheb,icheb2)*vhlr(lm3,lm2,icheb2) &
!             -tau(icheb,ipan)*hlk(l1,mn)*cslc1(icheb,icheb2)*vjlr(lm3,lm2,icheb2))&
!           *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
!             srv(icheb,lm1,icheb2,lm2) = &
!           (-tau(icheb,ipan)*jlk(l1,mn)*csrc1(icheb,icheb2)*vhli(lm3,lm2,icheb2) &
!             +tau(icheb,ipan)*hlk(l1,mn)*csrc1(icheb,icheb2)*vjli(lm3,lm2,icheb2)) &
!               *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0
            end do
          end do
        end do
      end do
    end do
    do lm1 = 1,lmsize2
      do icheb = 0,ncheb
        slv(icheb,lm1,icheb,lm1) = slv(icheb,lm1,icheb,lm1) + 1.d0
        srv(icheb,lm1,icheb,lm1) = srv(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do
  elseif  ( use_sratrick==1 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
         taucslcr = tau(icheb,ipan)*cslc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
         taucsrcr = tau(icheb,ipan)*csrc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0 
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize
!          do ivec=1,1
            do lm3 = 1,lmsize
!             lm1=lm3+(ivec-1)*lmsize
              lm1=lm3
              l1 = jlk_index(lm1)

              ! this is the block to be inverted in SRAtrick. (named C in comment above):
              slv1(icheb,lm1,icheb2,lm2) = &
              taucslcr*(jlk(l1,mn)*vhlr(lm3,lm2,icheb2) &
                       -hlk(l1,mn)*vjlr(lm3,lm2,icheb2))
              srv1(icheb,lm1,icheb2,lm2) = &
              taucsrcr*(-jlk(l1,mn)*vhli(lm3,lm2,icheb2) &
                        +hlk(l1,mn)*vjli(lm3,lm2,icheb2))

            end do
!         end do
        end do
      end do
    end do
!   do icheb2 = 0,ncheb
!     do icheb = 0,ncheb
!        taucslcr = tau(icheb,ipan)*cslc1(icheb,icheb2) &
!                   *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
!        taucsrcr = tau(icheb,ipan)*csrc1(icheb,icheb2) &
!                   *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0 
!       mn = ipan*ncheb + ipan - icheb
!       do lm2 = 1,lmsize
!         do ivec=2,2
!           do lm3 = 1,lmsize
!             lm1=lm3+(ivec-1)*lmsize
!             lm1=lm3+lmsize
!             l1 = jlk_index(lm1)

!             slv2(icheb,lm3,icheb2,lm2) = &
!             taucslcr*(jlk(l1,mn)*vhlr(lm3,lm2,icheb2) &
!                      -hlk(l1,mn)*vjlr(lm3,lm2,icheb2))
!             srv2(icheb,lm3,icheb2,lm2) = &
!             taucsrcr*(-jlk(l1,mn)*vhli(lm3,lm2,icheb2) &
!                       +hlk(l1,mn)*vjli(lm3,lm2,icheb2))

!           end do
!         end do
!       end do
!     end do
!   end do
    do lm1 = 1,lmsize
      do icheb = 0,ncheb
        slv1(icheb,lm1,icheb,lm1) = slv1(icheb,lm1,icheb,lm1) + 1.d0
        srv1(icheb,lm1,icheb,lm1) = srv1(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do

  else
    stop '[rllsll] error in inversion'
  end if

!  if (idotime==1) call timing_pause('local1')
!  if (idotime==1) call timing_start('local2')

!-------------------------------------------------------
! determine the local solutions
! solve the equations SLV*YRLL=S and SLV*ZRLL=C 
!                 and SRV*YILL=C and SRV*ZILL=S
! i.e., solve system A*U=J, see eq. 5.68.

  if ( use_sratrick==0 ) then
    nplm = (ncheb+1)*lmsize2

    if (cmoderll/='0') then
!      if (idotime==1) call timing_start('inversion')
      call zgetrf(nplm,nplm,slv,nplm,ipiv,info)
!      if (idotime==1) call timing_stop('inversion','test')
      if (info/=0) stop 'rllsll: zgetrf'
      call zgetrs('n',nplm,lmsize,slv,nplm,ipiv,yrll,nplm,info)
      call zgetrs('n',nplm,lmsize,slv,nplm,ipiv,zrll,nplm,info)
    end if
    if (cmodesll/='0') then
!      if (directsolv==1) then
        call zgetrf(nplm,nplm,srv,nplm,ipiv,info)
        if (info/=0) stop 'rllsll: zgetrf'
        call zgetrs('n',nplm,lmsize,srv,nplm,ipiv,yill,nplm,info)
        call zgetrs('n',nplm,lmsize,srv,nplm,ipiv,zill,nplm,info)
!      else
!        call iterativesol (ncheb,lmsize2,lmsize,srv,yill)
!        call iterativesol (ncheb,lmsize2,lmsize,srv,zill)
!      end if
    end if
  elseif ( use_sratrick==1 ) then
    nplm = (ncheb+1)*lmsize

    call zgetrf(nplm,nplm,slv1,nplm,ipiv,info)
    call zgetrs('n',nplm,lmsize,slv1,nplm,ipiv,yrll1,nplm,info)
    call zgetrs('n',nplm,lmsize,slv1,nplm,ipiv,zrll1,nplm,info)

!   call zgemm('n','n',nplm,lmsize,nplm,-cone,slv2, &
!       nplm,yrll1,nplm,cone,yrll2,nplm)

!   call zgemm('n','n',nplm,lmsize,nplm,-cone,slv2, &
!       nplm,zrll1,nplm,cone,zrll2,nplm)

    call zgetrf(nplm,nplm,srv1,nplm,ipiv,info)
    call zgetrs('n',nplm,lmsize,srv1,nplm,ipiv,yill1,nplm,info)
    call zgetrs('n',nplm,lmsize,srv1,nplm,ipiv,zill1,nplm,info)

!   call zgemm('n','n',nplm,lmsize,nplm,-cone,srv2, &
!       nplm,yill1,nplm,cone,yill2,nplm)

!   call zgemm('n','n',nplm,lmsize,nplm,-cone,srv2, &
!       nplm,zill1,nplm,cone,zill2,nplm)

    do icheb2 = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrll1temp(lm1,lm2) = yrll1(icheb2,lm1,lm2)
          zrll1temp(lm1,lm2) = zrll1(icheb2,lm1,lm2)
          yill1temp(lm1,lm2) = yill1(icheb2,lm1,lm2)
          zill1temp(lm1,lm2) = zill1(icheb2,lm1,lm2)
        end do
      end do
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhlr(1,1,icheb2), &
        lmsize,yrll1temp,lmsize,czero,vhlr_yrll1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhlr(1,1,icheb2), &
        lmsize,zrll1temp,lmsize,czero,vhlr_zrll1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjlr(1,1,icheb2), &
        lmsize,yrll1temp,lmsize,czero,vjlr_yrll1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjlr(1,1,icheb2), &
        lmsize,zrll1temp,lmsize,czero,vjlr_zrll1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhli(1,1,icheb2), &
        lmsize,yill1temp,lmsize,czero,vhli_yill1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhli(1,1,icheb2), &
        lmsize,zill1temp,lmsize,czero,vhli_zill1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjli(1,1,icheb2), &
        lmsize,yill1temp,lmsize,czero,vjli_yill1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjli(1,1,icheb2), &
        lmsize,zill1temp,lmsize,czero,vjli_zill1,lmsize)

      do icheb = 0,ncheb
         taucslcr = - tau(icheb,ipan)*cslc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0    ! *(b-a)/2 in eq. 5.53, 5.54
         taucsrcr =  tau(icheb,ipan)*csrc1(icheb,icheb2) &
                    *(rpanbound(ipan)-rpanbound(ipan-1))/ 2.d0 
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize
            do lm3 = 1,lmsize
              lm1=lm3+lmsize
              l1 = jlk_index(lm1)

              yrll2(icheb,lm3,lm2) = &
              yrll2(icheb,lm3,lm2) + &
              taucslcr*(jlk(l1,mn)*vhlr_yrll1(lm3,lm2) &
                       -hlk(l1,mn)*vjlr_yrll1(lm3,lm2))

              zrll2(icheb,lm3,lm2) = &
              zrll2(icheb,lm3,lm2) + &
              taucslcr*(jlk(l1,mn)*vhlr_zrll1(lm3,lm2) &
                       -hlk(l1,mn)*vjlr_zrll1(lm3,lm2))

              yill2(icheb,lm3,lm2) = &
              yill2(icheb,lm3,lm2) + &
              taucsrcr*(jlk(l1,mn)*vhli_yill1(lm3,lm2) &
                       -hlk(l1,mn)*vjli_yill1(lm3,lm2))

              zill2(icheb,lm3,lm2) = &
              zill2(icheb,lm3,lm2) + &
              taucsrcr*(jlk(l1,mn)*vhli_zill1(lm3,lm2) &
                       -hlk(l1,mn)*vjli_zill1(lm3,lm2))

            end do
        end do
      end do
    end do

  else
    stop '[rllsll] error in inversion'
  end if

  ! Reorient indices for later use
  if ( use_sratrick==0 ) then
    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize2
          yrf(lm1,lm2,icheb,ipan) = yrll(icheb,lm1,lm2)
          zrf(lm1,lm2,icheb,ipan) = zrll(icheb,lm1,lm2)
          yif(lm1,lm2,icheb,ipan) = yill(icheb,lm1,lm2)
          zif(lm1,lm2,icheb,ipan) = zill(icheb,lm1,lm2)
        end do
      end do
    end do

  elseif ( use_sratrick==1 ) then

    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrf(lm1,lm2,icheb,ipan) = yrll1(icheb,lm1,lm2)
          zrf(lm1,lm2,icheb,ipan) = zrll1(icheb,lm1,lm2)
          yif(lm1,lm2,icheb,ipan) = yill1(icheb,lm1,lm2)
          zif(lm1,lm2,icheb,ipan) = zill1(icheb,lm1,lm2)
        end do
      end do
    end do

    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yrf(lm1+lmsize,lm2,icheb,ipan) = yrll2(icheb,lm1,lm2)
          zrf(lm1+lmsize,lm2,icheb,ipan) = zrll2(icheb,lm1,lm2)
          yif(lm1+lmsize,lm2,icheb,ipan) = yill2(icheb,lm1,lm2)
          zif(lm1+lmsize,lm2,icheb,ipan) = zill2(icheb,lm1,lm2)
        end do
      end do
    end do

  end if

!  if (idotime==1) call timing_pause('local2')
!  if (idotime==1) call timing_start('local3')

  ! Calculation of eq. 5.19-5.22

  do icheb = 0,ncheb
    zslc1sum(icheb) = slc1sum(icheb) * (rpanbound(ipan)-rpanbound(ipan-1))/ (2.d0)
  end do
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhlr(1,1,0), &
        lmsize,yrf(1,1,0,ipan),lmsize2,czero,mrnvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjlr(1,1,0), &
        lmsize,yrf(1,1,0,ipan),lmsize2,czero,mrjvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhlr(1,1,0), &
        lmsize,zrf(1,1,0,ipan),lmsize2,czero,mrnvz(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjlr(1,1,0), &
        lmsize,zrf(1,1,0,ipan),lmsize2,czero,mrjvz(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhli(1,1,0), &
        lmsize,yif(1,1,0,ipan),lmsize2,czero,mihvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjli(1,1,0), &
        lmsize,yif(1,1,0,ipan),lmsize2,czero,mijvy(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhli(1,1,0), &
        lmsize,zif(1,1,0,ipan),lmsize2,czero,mihvz(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjli(1,1,0), &
        lmsize,zif(1,1,0,ipan),lmsize2,czero,mijvz(1,1,ipan),lmsize)
  do icheb = 1,ncheb
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhlr(1,1,icheb), &
          lmsize,yrf(1,1,icheb,ipan),lmsize2,cone,mrnvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjlr(1,1,icheb), &
          lmsize,yrf(1,1,icheb,ipan),lmsize2,cone,mrjvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhlr(1,1,icheb), &
          lmsize,zrf(1,1,icheb,ipan),lmsize2,cone,mrnvz(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjlr(1,1,icheb), &
          lmsize,zrf(1,1,icheb,ipan),lmsize2,cone,mrjvz(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhli(1,1,icheb), &
          lmsize,yif(1,1,icheb,ipan),lmsize2,cone,mihvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjli(1,1,icheb), &
          lmsize,yif(1,1,icheb,ipan),lmsize2,cone,mijvy(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhli(1,1,icheb), &
          lmsize,zif(1,1,icheb,ipan),lmsize2,cone,mihvz(1,1,ipan),lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjli(1,1,icheb), &
          lmsize,zif(1,1,icheb,ipan),lmsize2,cone,mijvz(1,1,ipan),lmsize)
  end do
!  if (idotime==1) call timing_pause('local3')

end do !ipan
#ifdef CPP_hybrid
!$NOOMP END DO
!$NOOMP END PARALLEL
#endif
! end the big loop over the subintervals



!if (idotime==1) call timing_stop('local')
!if (idotime==1) call timing_start('afterlocal')

! ***********************************************************************
! calculate A(M), B(M), C(M), D(M)
! according to 5.17-5.18 (regular solution) of Bauer PhD
! C,D are calculated accordingly for the irregular solution
! (starting condition: A(0) = 1, B(0) = 0, C(MMAX) = 0 and D(MMAX) = 1)
! ***********************************************************************

! regular 
do lm2 = 1,lmsize
  do lm1 = 1,lmsize
    bllp(lm1,lm2,0) = czero
    allp(lm1,lm2,0) = czero
  end do
end do

do lm1 = 1,lmsize
  allp(lm1,lm1,0) = cone
end do

do ipan = 1,npan
  call zcopy(lmsize*lmsize,allp(1,1,ipan-1),1,allp(1,1,ipan),1)
  call zcopy(lmsize*lmsize,bllp(1,1,ipan-1),1,bllp(1,1,ipan),1)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvy(1,1,ipan), &
          lmsize,allp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvz(1,1,ipan), &
          lmsize,bllp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mrjvy(1,1,ipan), &
          lmsize,allp(1,1,ipan-1),lmsize,cone,bllp(1,1,ipan),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mrjvz(1,1,ipan), &
          lMSIZE,BLLP(1,1,IPAN-1),LMSIZE,CONE,BLLP(1,1,IPAN),LMSIZE)
end do

! irregular 
do lm2 = 1,lmsize
  do lm1 = 1,lmsize
    dllp(lm1,lm2,npan) = 0.d0
    cllp(lm1,lm2,npan) = 0.d0
  end do
end do
do lm1 = 1,lmsize
  dllp(lm1,lm1,npan) = 1.d0
end do
do ipan = npan,1,-1
  call zcopy(lmsize*lmsize,cllp(1,1,ipan),1,cllp(1,1,ipan-1),1)
  call zcopy(lmsize*lmsize,dllp(1,1,ipan),1,dllp(1,1,ipan-1),1)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mihvz(1,1,ipan), &
          lmsize,cllp(1,1,ipan),lmsize,cone,cllp(1,1,ipan-1),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize, cone,mihvy(1,1,ipan), &
          lmsize,dllp(1,1,ipan),lmsize,cone,cllp(1,1,ipan-1),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mijvz(1,1,ipan), &
          lmsize,cllp(1,1,ipan),lmsize,cone,dllp(1,1,ipan-1),lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mijvy(1,1,ipan), &
          lmsize,dllp(1,1,ipan),lmsize,cone,dllp(1,1,ipan-1),lmsize)
end do

! ***********************************************************************
! determine the regular solution ull by using 5.14
! and the irregular solution sll accordingly
! ***********************************************************************
do ipan = 1,npan
  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,yrf(1,1,icheb,ipan), &
          lmsize2,allp(1,1,ipan-1),lmsize,czero,ull(1,1,mn),lmsize2)
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,zrf(1,1,icheb,ipan), &
          lmsize2,bllp(1,1,ipan-1),lmsize,cone,ull(1,1,mn),lmsize2)
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,zif(1,1,icheb,ipan), &
          lmsize2,cllp(1,1,ipan),lmsize,czero,sll(1,1,mn),lmsize2)
  call zgemm('n','n',lmsize2,lmsize,lmsize,cone,yif(1,1,icheb,ipan), &
          lmsize2,dllp(1,1,ipan),lmsize,cone,sll(1,1,mn),lmsize2)
  end do
end do

!if (idotime==1) call timing_stop('afterlocal')
!if (idotime==1) call timing_start('endstuff')

! ***********************************************************************
! next part converts the volterra solution u of equation (5.7) to
! the fredholm solution r by employing eq. 4.122 and 4.120 of bauer, phd
! and the t-matrix is calculated
! ***********************************************************************

call zgetrf(lmsize,lmsize,allp(1,1,npan),lmsize,ipiv,info)                     !invert alpha
call zgetri(lmsize,allp(1,1,npan),lmsize,ipiv,work,lmsize*lmsize,info)         !invert alpha -> transformation matrix rll=alpha^-1*rll
#ifdef hostcode
! get alpha matrix
!      DO LM1=1,LMSIZE                          ! LLY
!       DO LM2=1,LMSIZE                         ! LLY
!        ALPHAGET(LM1,LM2)=ALLP(LM1,LM2,NPAN)   ! LLY
!       ENDDO                                   ! LLY
!      ENDDO                                    ! LLY
#endif
! calculation of the t-matrix 
call zgemm('n','n',lmsize,lmsize,lmsize,cone/gmatprefactor,bllp(1,1,npan), &   ! calc t-matrix tll = bll*alpha^-1 
            lmsize,allp(1,1,npan),lmsize,czero,tllp,lmsize)

do nm = 1,nrmax
call zgemm('n','n',lmsize2,lmsize,lmsize,cone,ull(1,1,nm), &
            lmsize2,allp(1,1,npan),lmsize,czero,rll(1,1,nm),lmsize2)
end do

!if (idotime==1) call timing_stop('endstuff')
!if (idotime==1) call timing_start('checknan')
!if (idotime==1) call timing_stop('checknan')
!if (idotime==1) call timing_stop('local1')
!if (idotime==1) call timing_stop('local2')
!if (idotime==1) call timing_stop('local3')
!if (idotime==1) call timing_stop('rllsll')

if ( use_sratrick==0 ) then
  if(allocated(slv)) deallocate ( slv,srv )
elseif ( use_sratrick==1 ) then
  if(allocated(work2)) deallocate ( work2, ipiv2 )
  if(allocated(slv1)) deallocate ( slv1, srv1 )
! if(allocated(slv2)) deallocate ( slv2, srv2 )
  if(allocated(yill1)) deallocate ( yill1, zill1 )
  if(allocated(yrll1)) deallocate ( yrll1, zrll1 )
  if(allocated(yill2)) deallocate ( yill2, zill2 )
  if(allocated(yrll2)) deallocate ( yrll2, zrll2 )
end if

if(allocated(work)) deallocate( work )
if(allocated(allp)) deallocate( allp, bllp )
if(allocated(cllp)) deallocate( cllp, dllp )
if(allocated(mrnvy)) deallocate( mrnvy, mrnvz )
if(allocated(mrjvy)) deallocate( mrjvy, mrjvz )
if(allocated(mihvy)) deallocate( mihvy, mihvz )
if(allocated(mijvy)) deallocate( mijvy, mijvz )
if(allocated(yill)) deallocate( yill, zill )
if(allocated(yrll)) deallocate( yrll, zrll )
if(allocated(vjlr)) deallocate( vjlr, vhlr )
if(allocated(vjli)) deallocate( vjli, vhli )
if(allocated(vjlr_yrll1)) deallocate( vjlr_yrll1, vhlr_yrll1 )
if(allocated(vjli_yill1)) deallocate( vjli_yill1, vhli_yill1 )
if(allocated(vjlr_zrll1)) deallocate( vjlr_zrll1, vhlr_zrll1 )
if(allocated(vjli_zill1)) deallocate( vjli_zill1, vhli_zill1 )
if(allocated(yrll1temp)) deallocate( yrll1temp, zrll1temp )
if(allocated(yill1temp)) deallocate( yill1temp, zill1temp )

if(allocated(yif)) deallocate( yif )
if(allocated(yrf)) deallocate( yrf )
if(allocated(zif)) deallocate( zif )
if(allocated(zrf)) deallocate( zrf )

end subroutine

#ifndef hostcode
END MODULE MOD_RLLSLL
#endif

subroutine getCLambdaCinv(Ncheb,CLambdaCinv)
implicit none
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion
integer          :: Ncheb
double precision :: CLambdaCinv(0:Ncheb,0:Ncheb)
!local
double precision :: Lambda(0:Ncheb,0:Ncheb)
double precision :: Cmatrix(0:Ncheb,0:Ncheb)
double precision :: Cinvmatrix(0:Ncheb,0:Ncheb)
double precision :: temp1(0:Ncheb,0:Ncheb)
integer n
 Lambda=(0.0D0,0.0D0)
 Cmatrix=(0.0D0,0.0D0)
 Cinvmatrix=(0.0D0,0.0D0)
 Lambda=(0.0D0,0.0D0)
 temp1=(0.0D0,0.0D0)

call getLambda(Ncheb,Lambda)
call getCinvmatrix(Ncheb,Cinvmatrix)
call getCmatrix(Ncheb,Cmatrix)
n=Ncheb+1
 call dgemm('N','N',n,n,n,1d0,Lambda,n,Cinvmatrix,n,0d0,temp1,n)
 call dgemm('N','N',n,n,n,1d0,Cmatrix,n,temp1,n,0d0,CLambdaCinv,n)
! temp1=matmat_dmdm(Lambda,Cinvmatrix,Ncheb)
! CLambdaCinv=matmat_dmdm(Cmatrix,temp1,Ncheb)

end subroutine

subroutine rotatematrix(mat,theta,phi,lmmax,mode)
! rotates a matrix in the local frame pointing in
! the direction of phi and theta to the global frame
implicit none
!interface
double complex,intent(inout)    ::  mat(2*lmmax,2*lmmax)
double precision,intent(in)     :: phi
double precision,intent(in)     :: theta
integer                         :: lmmax
integer                         :: mode
!local
double complex   :: Umat(2*lmmax,2*lmmax)
double complex   :: Udeggamat(2*lmmax,2*lmmax)
double complex   :: mattemp(2*lmmax,2*lmmax)
!double precision :: matmat_zmzm

!***********************************************************************
! create the rotation matrix:
!     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
!  U= |                                                          |
!     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
!
!  Udegga = transpose(complex conjug ( U ) )
!***********************************************************************


call create_Umatrix(theta,phi,lmmax,Umat,Udeggamat)
!***********************************************************************
! calculate matrix in the global frame:
!
!  t_glob = U * t_loc * Udegga
!***********************************************************************


if (mode==0) then ! 'loc->glob'
  call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),mat,2*lmmax,Udeggamat,2*lmmax,(0d0,0d0),mattemp,2*lmmax)
  call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Umat,2*lmmax,mattemp,2*lmmax,(0d0,0d0),mat,2*lmmax)
elseif (mode==1) then !'glob->loc'
  call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),mat,2*lmmax,Umat,2*lmmax,(0d0,0d0),mattemp,2*lmmax)
  call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Udeggamat,2*lmmax,mattemp,2*lmmax,(0d0,0d0),mat,2*lmmax)
else
  stop '[rotatematrix] mode not known'
end if
!  writE(324,'(5000F)') tmat
! stop

end subroutine rotatematrix


SUBROUTINE spin_orbit_compl(lmax,lmmaxd,l_s)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: lmax
INTEGER, INTENT(IN)        :: lmmaxd
DOUBLE COMPLEX, INTENT(OUT):: l_s(:,:)
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics


!  local variableINTEGER    ::     i1,i2,i1l,rl,lm1,lm2
INTEGER    ::     rl,lm1,lm2
DOUBLE COMPLEX,allocatable     ::     ls_l(:,:)


!icompl=(0D0,1D0)


CALL cinit((2*lmmaxd)**2,l_s)

DO rl=0,lmax
  
  allocate(ls_l((2*rl+1)*2,(2*rl+1)*2))
  CALL cinit(((2*rl+1)*2)**2,ls_l)
  
  
  CALL spin_orbit_one_l(rl,ls_l)
  
  DO lm1=1,(2*rl+1)*2
    
    IF (lm1 <= 2*rl+1 ) THEN
      DO lm2=1,(2*rl+1)
        l_s(rl**2+lm1,rl**2+lm2)=0.5D0*ls_l(lm1,lm2)
      END DO
      DO lm2=(2*rl+1)+1,(2*rl+1)*2
        l_s(rl**2+lm1,lmmaxd+rl**2-(2*rl+1)+lm2)= 0.5D0*ls_l(lm1,lm2)
      END DO
    ELSE
      DO lm2=1,(2*rl+1)
        l_s(lmmaxd+rl**2-(2*rl+1)+lm1,rl**2+lm2)= 0.5D0*ls_l(lm1,lm2)
      END DO
      DO lm2=(2*rl+1)+1,(2*rl+1)*2
        l_s(lmmaxd+rl**2-(2*rl+1)+lm1,lmmaxd+rl**2-(2*rl+1)+lm2)=  &
            0.5D0*ls_l(lm1,lm2)
      END DO
    END IF
    
  END DO    !lm1
  
  deallocate(ls_l)
  
  
END DO     !rl=0,lmax


END SUBROUTINE spin_orbit_compl


SUBROUTINE beshank(hl,jl,z,lmax)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-19  Time: 12:22:05
 
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders lmin .le. l .le. lmax.
!  For |z| .lt. l+1 the taylor expansions of jl and nl are used.
!  For |z| .ge. l+1 the explicit expressions for hl(+), hl(-) are used.
!-----------------------------------------------------------------------
!     .. Parameters ..
DOUBLE COMPLEX ci
PARAMETER (ci= (0.0D0,1.0D0))
!     ..
!     .. Scalar Arguments ..
DOUBLE COMPLEX z
INTEGER :: lmax
!     ..
!     .. Array Arguments ..
DOUBLE COMPLEX hl(0:lmax),jl(0:lmax),nl(0:lmax)
!     ..
!     .. Local Scalars ..
DOUBLE COMPLEX termj,termn,z2,zj,zn
DOUBLE PRECISION :: rl,rn,rnm
INTEGER :: l,m,n
!     ..
!     .. Intrinsic Functions ..
INTRINSIC CDABS,EXP
!     ..
zj = 1.d0
zn = 1.d0
z2 = z*z
IF (CDABS(z) < lmax+1.d0) THEN
  DO  l = 0,lmax
    rl = l + l
    termj = -0.5D0/ (rl+3.d0)*z2
    termn = 0.5D0/ (rl-1.d0)*z2
    jl(l) = 1.d0
    nl(l) = 1.d0
    DO  n = 2,25
      jl(l) = jl(l) + termj
      nl(l) = nl(l) + termn
      rn = n + n
      termj = -termj/ (rl+rn+1.d0)/rn*z2
      termn = termn/ (rl-rn+1.d0)/rn*z2
    END DO
    jl(l) = jl(l)*zj
    nl(l) = -nl(l)*zn/z
    hl(l) = jl(l) + nl(l)*ci
    
    zj = zj*z/ (rl+3.d0)
    zn = zn/z* (rl+1.d0)
  END DO
END IF

DO  l = 0,lmax
  IF (CDABS(z) >= l+1.d0) THEN
    hl(l) = 0.d0
    nl(l) = 0.d0
    rnm = 1.d0
    DO  m = 0,l
      hl(l) = hl(l) + rnm/ (-ci* (z+z))**m
      nl(l) = nl(l) + rnm/ (ci* (z+z))**m
      rnm = rnm* (l*l+l-m*m-m)/ (m+1.d0)
    END DO
    hl(l) = hl(l)* (-ci)**l*EXP(ci*z)/ (ci*z)
    nl(l) = nl(l)*ci**l*EXP(-ci*z)/ (-ci*z)
    jl(l) = (hl(l)+nl(l))*0.5D0
    nl(l) = (hl(l)-jl(l))/ci
  END IF
END DO

RETURN

END SUBROUTINE

SUBROUTINE beshank_smallcomp(hl,jl,zval,tau,eryd,lmax)
IMPLICIT NONE
!-----------------------------------------------------------------------
!  takes the spherical bessel etc functions stored in an array up to LMAX
!  array entries from LMAX+1 to 2*LMAX are assumed to be empty
!  these values are filled with the potential-free solution of the
!  SRA-equations
!-----------------------------------------------------------------------
DOUBLE COMPLEX hl(0:2*(lmax+1)-1), jl(0:2*(lmax+1)-1),  &
    nl(0:2*(lmax+1)-1)
DOUBLE PRECISION :: cvlight
PARAMETER (cvlight=274.0720442D0)
DOUBLE COMPLEX zval
DOUBLE COMPLEX eryd
DOUBLE PRECISION :: tau
INTEGER :: lmax

!       DOUBLE PRECISION CVLIGHT
DOUBLE COMPLEX prefac
INTEGER :: il,il2


prefac = 1.0D0 / (1.0D0+eryd/cvlight**2) / tau !/cvlight  !last cvlight for small component test

il=0
il2=il+lmax+1
nl(il2)=prefac * (zval* (-nl(il+1)) )
jl(il2)=prefac * (zval* (-jl(il+1)) )
!       HL(IL2)=JL(IL2)+ CI*NL(IL2)
hl(il2)=prefac * (zval* (-hl(il+1)) )
!       write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
!       write(*,'(5000E)') tau,HL(0),JL(0)+ (0.0D0,1.0D0)*NL(0)

prefac = 1.0D0 / (1.0D0+eryd/cvlight**2) / tau !/cvlight !last cvlight for small component test

DO il=1,lmax
  il2=il+lmax+1
  nl(il2)=prefac * ( zval * nl(il-1)-(il+1)*nl(il) )
  jl(il2)=prefac * ( zval * jl(il-1)-(il+1)*jl(il) )
!         HL(IL2)=JL(IL2)+ CI*NL(IL2)
  hl(il2)=prefac * ( zval * hl(il-1)-(il+1)*hl(il) )
!         HL(IL2)=PREFAC * ( ZVAL * HL(IL-1)-(IL+1)*HL(IL) )
!         write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
END DO

END SUBROUTINE beshank_smallcomp


SUBROUTINE chebint(cslc1,csrc1,slc1sum,c1,n)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-19  Time: 14:23:20
 
!---------------------------------------------------------------------
! this subroutine calculates the matrices for the Chebyshev integration
! as defined on page 141 and 142 of the article:
! Integral Equation Method for the Continuous Spectrum Radial
! Schroedinger Equation by R. A. Gonzales et al
! in Journal of computational physics 134, 134-149 (1997)

! the matrix C is the discrete cosine transform matrix
! the matrix C1 is the inverse of C
! the matrix SL is the left spectral integration matrix
! the matrix SR is the right spectral integration matrix
! the matrix CSLC1 is the product of C, SL and C1
! the matrix CSRC1 is the product of C, SR and C1
!---------------------------------------------------------------------
!     .. Local Scalars ..
DOUBLE PRECISION :: pi
INTEGER :: j,k
!     ..
!     .. Local Arrays ..
DOUBLE PRECISION :: c(0:n,0:n),c1(0:n,0:n),s1(0:n,0:n),s2(0:n,0:n),  &
    sl(0:n,0:n),slc1(0:n,0:n),sr(0:n,0:n), src1(0:n,0:n)
!     ..
!     .. External Subroutines ..
EXTERNAL dgemm
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ATAN,COS
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION :: cslc1(0:n,0:n),csrc1(0:n,0:n),slc1sum(0:n)
!     ..
!     .. Scalar Arguments ..
INTEGER :: n
!     ..
pi = 4.d0*ATAN(1.d0)
!---------------------------------------------------------------------
! determine the discrete cosine transform matrix from the zeros of the
! Chebyshev polynomials
DO j = 0,n
  DO k = 0,n
    c(k,j) = COS(((2*k+1)*j*pi)/ (2* (n+1)))
  END DO
END DO
!---------------------------------------------------------------------
! determine the inverse of the discrete cosine transform matrix from
! the transpose of the discrete cosine transform matrix
DO j = 0,n
  DO k = 0,n
    c1(k,j) = c(j,k)*2.d0/ (n+1)
  END DO
  c1(0,j) = c1(0,j)*0.5D0
END DO
!---------------------------------------------------------------------
! next to statements can be used to check the products CT*C and C1*C
CALL dgemm('T','N',n+1,n+1,n+1,1.d0,c,n+1,c,n+1,0.d0,sr,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,c1,n+1,c,n+1,0.d0,sr,n+1)
!---------------------------------------------------------------------
! preparation of the left and right
! spectral integration matrices SL and SR
DO j = 0,n
  DO k = 0,n
    s1(k,j) = 0.0D0
    s2(k,j) = 0.0D0
  END DO
END DO
DO j = 0,n
  s1(0,j) = (-1.d0)** (j+1)
  s1(j,j) = 1.d0
END DO
DO j = 2,n - 1
  s2(j,j-1) = 0.5D0/j
  s2(j,j+1) = -0.5D0/j
END DO
s2(n,n-1) = 0.5D0/n
s2(1,0) = 1.d0
s2(1,2) = -0.5D0
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,s1,n+1,s2,n+1,0.d0,sl,n+1)
DO j = 0,n
  DO k = 0,n
    s1(k,j) = 0.0D0
  END DO
END DO
DO j = 0,n
  s1(j,j) = -1.d0
  s1(0,j) = 1.d0
END DO
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,s1,n+1,s2,n+1,0.d0,sr,n+1)
!---------------------------------------------------------------------
! determination of the products C*SL*C1 and C*SR*C1
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,sl,n+1,c1,n+1,0.d0,slc1,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,c,n+1,slc1,n+1,0.d0,cslc1,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,sr,n+1,c1,n+1,0.d0,src1,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,c,n+1,src1,n+1,0.d0,csrc1,n+1)
!---------------------------------------------------------------------
DO k = 0,n
  slc1sum(k) = 0.0D0
  DO j = 0,n
    slc1sum(k) = slc1sum(k) + slc1(j,k)
  END DO
END DO
RETURN
END SUBROUTINE

subroutine getLambda(Ncheb,Lambda)
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion 
implicit none
integer          :: Ncheb
double precision :: Lambda(0:Ncheb,0:Ncheb)
!local
integer icheb,icheb2
do icheb2=1,Ncheb,2
  Lambda(0,icheb2)=icheb2
end do
do icheb=1,Ncheb
  do icheb2=icheb+1,Ncheb,2
    Lambda(icheb,icheb2)=icheb2*2
  end do
end do
end subroutine

subroutine getCinvmatrix(Ncheb,Cinvmatrix)
! calculates the C**-1 matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
implicit none
integer, intent(in)           :: ncheb
double precision, intent(out) :: Cinvmatrix(0:Ncheb,0:Ncheb)
!local
double precision              :: pi
integer                       :: icheb1,icheb2
double precision              :: fac

pi=4d0*datan(1d0)
fac=1.0D0/(Ncheb+1)
do icheb1=0,ncheb
  do icheb2=0,ncheb
    Cinvmatrix(icheb1,icheb2)=fac*dcos(icheb1*pi*((Ncheb-icheb2)+0.5D0)/(Ncheb+1))
  end do
  fac=2.0D0/(Ncheb+1)
end do

end subroutine getCinvmatrix

subroutine getCmatrix(Ncheb,Cmatrix)
! calculates the C matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
implicit none
integer, intent(in)           :: ncheb
double precision, intent(out) :: Cmatrix(0:Ncheb,0:Ncheb)
double precision              :: pi
!local
integer                       :: icheb1,icheb2

pi=4d0*datan(1d0)
do icheb1=0,ncheb
  do icheb2=0,ncheb
    ! maybe incorrect
    Cmatrix(icheb2,icheb1)=dcos(icheb1*pi*((Ncheb-icheb2)+0.5D0)/(Ncheb+1))
  end do
end do
end subroutine getCmatrix

subroutine create_Umatrix(theta,phi,lmmax,Umat,Udeggamat)
implicit none
!***********************************************************************
! create the rotation matrix:
!     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
!  U= |                                                          |
!     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
!
!  Udegga = transpose(complex conjug ( U ) )
!***********************************************************************double
!precision :: phi
!interface
double precision,intent(in)     :: phi 
double precision,intent(in)     :: theta
integer,intent(in)              :: lmmax
double complex,intent(out)      :: Umat(2*lmmax,2*lmmax)
double complex,intent(out)      :: Udeggamat(2*lmmax,2*lmmax)
!local
double complex                  :: Umat11,Umat12,Umat21,Umat22
double complex                  :: Udeggamat11,Udeggamat12,Udeggamat21,Udeggamat22
integer                         :: ival
double complex,parameter        :: ci=(0.0D0,1.0D0)
character*25               :: spinmode

spinmode='kkr'
if (spinmode=='regular') then
  Umat11      =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
  Umat12      = -sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
  Umat21      =  sin(theta/2.0D0)*exp( ci/2.0D0*phi)
  Umat22      =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
else if (spinmode=='kkr') then
  Umat11      =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
  Umat12      =  sin(theta/2.0D0)*exp( ci/2.0D0*phi)
  Umat21      = -sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
  Umat22      =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
else 
  stop '[create_Umatrix] mode not known'
end if

Umat=(0.0D0,0.0D0)
do ival=1,lmmax
  Umat(      ival,      ival) = Umat11
  Umat(      ival,lmmax+ival) = Umat12
  Umat(lmmax+ival,ival)       = Umat21
  Umat(lmmax+ival,lmmax+ival) = Umat22
end do

if (spinmode=='regular') then
Udeggamat11 =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
Udeggamat12 =  sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
Udeggamat21 = -sin(theta/2.0D0)*exp( ci/2.0D0*phi)
Udeggamat22 =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
else if (spinmode=='kkr') then
Udeggamat11 =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
Udeggamat12 = -sin(theta/2.0D0)*exp( ci/2.0D0*phi)
Udeggamat21 =  sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
Udeggamat22 =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
else 
  stop '[create_Umatrix] mode not known'
end if



Udeggamat=(0.0D0,0.0D0)
do ival=1,lmmax
  Udeggamat(      ival,      ival) = Udeggamat11
  Udeggamat(      ival,lmmax+ival) = Udeggamat12
  Udeggamat(lmmax+ival,ival)       = Udeggamat21
  Udeggamat(lmmax+ival,lmmax+ival) = Udeggamat22
end do

end subroutine create_Umatrix

SUBROUTINE spin_orbit_one_l(lmax,l_s)

IMPLICIT NONE

INTEGER, INTENT(IN)                  :: lmax
DOUBLE COMPLEX, INTENT(OUT)    :: l_s((2*lmax+1)*2,(2*lmax+1)*2)
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics

!      schematically it has the form
!      (  -L_z    L_+  )
!      (  L_-     L_z  )




!  local variables
INTEGER                     ::    i1,i2,i1l
DOUBLE COMPLEX              ::    icompl
DOUBLE COMPLEX,allocatable  ::    l_min(:,:)
DOUBLE COMPLEX,allocatable  ::    l_up(:,:)
DOUBLE PRECISION            ::    lfac



icompl=(0D0,1D0)


allocate(l_min(-lmax:lmax,-lmax:lmax))
allocate(l_up(-lmax:lmax,-lmax:lmax))

!  initialize the matrix

DO i1=1,(2*lmax+1)*2
  DO i2=1,(2*lmax+1)*2
    l_s(i2,i1)=0D0
  END DO
END DO

DO i1=-lmax,lmax
  DO i2=-lmax,lmax
    l_min(i2,i1)=0D0
    l_up(i2,i1)=0D0
  END DO
END DO

!  fill the second and the forth quadrant with L_z
! (-L_z,respectively)


DO i1=1,2*lmax+1
  i1l=i1-lmax-1       ! the value of m (varies from -l to +l)
  i2=2*lmax+1-(i1-1)
  
!         L_S(i2,i1)=icompl*i1l
  l_s(i2,i1)=-icompl*i1l
  
END DO

DO i1=2*lmax+2,(2*lmax+1)*2
  i1l=i1-lmax-1-(2*lmax+1)       ! the value of m (varies from -l to +l)
  i2=(2*lmax+1)*2-(i1-(2*lmax+2))
  
!         L_S(i2,i1)=-icompl*i1l
  l_s(i2,i1)=icompl*i1l
  
END DO


!  implement now L_- in the third quadrant

IF (lmax>0) THEN
  
  lfac=SQRT(lmax*(lmax+1D0))/SQRT(2D0)
  l_min(0,-1)=-icompl*lfac
!         l_min(0,-1)=icompl*lfac
  l_min(0,1)=lfac
  l_min(-1,0)=icompl*lfac
  l_min(1,0)=-lfac
  
  IF (lmax > 1) THEN
    
    DO i1=2,lmax
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-i1*(i1-1D0))
      l_min(-i1,-i1+1)=-lfac
      l_min(-i1,i1-1)=icompl*lfac
      l_min(i1,-i1+1)=-icompl*lfac
      l_min(i1,i1-1)=-lfac
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-(i1-1)*(i1))
      l_min(-i1+1,-i1)=lfac
      l_min(-i1+1,i1)=icompl*lfac
      l_min(i1-1,-i1)=-icompl*lfac
      l_min(i1-1,i1)=lfac
      
    END DO
    
  END IF
END IF


DO i1=-lmax,lmax
  DO i2=-lmax,lmax
    l_s(i2+3*lmax+2,i1+lmax+1)=l_min(i1,i2)
  END DO
END DO


!  implement now L_+ in the   quadrant

IF (lmax>0) THEN
  
  lfac=SQRT(lmax*(lmax+1D0))/SQRT(2D0)
  l_up(0,-1)=-icompl*lfac
  l_up(0,1)=-lfac
  l_up(-1,0)=icompl*lfac
  l_up(1,0)=lfac
  
  IF (lmax > 1) THEN
    
    DO i1=2,lmax
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-i1*(i1-1D0))
      l_up(-i1,-i1+1)=lfac
      l_up(-i1,i1-1)=icompl*lfac
      l_up(i1,-i1+1)=-icompl*lfac
      l_up(i1,i1-1)=lfac
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-(i1-1)*(i1))
      l_up(-i1+1,-i1)=-lfac
      l_up(-i1+1,i1)=icompl*lfac
      l_up(i1-1,-i1)=-icompl*lfac
      l_up(i1-1,i1)=-lfac
      
    END DO
    
  END IF
END IF


DO i1=-lmax,lmax
  DO i2=-lmax,lmax
    l_s(i2+lmax+1,i1+3*lmax+2)=l_up(i1,i2)
  END DO
END DO



deallocate(l_min)
deallocate(l_up)


END SUBROUTINE spin_orbit_one_l

SUBROUTINE rhovalnew(ldorhoef,ielast,nsra,nspin,lmax,ez,wez,zat,  &
        socscale,cleb,icleb,iend,ifunm,lmsp,ncheb,  &
        npan_tot,npan_log,npan_eq,rmesh,irws,  &
        rpan_intervall,ipan_intervall,  &
        rnew,vinsnew,thetasnew,theta,phi,angle_fixed, &
        moment_x,moment_y,moment_z, &
        ipot,  &
        den_out,espv,rho2ns,r2nef,gmatn, muorb,  &
        lpotd,lmaxd,irmd,irmd_new,iemxd,soc) ! new parameters
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-21  Time: 11:39:57

 
#ifdef CPP_OMP
       use omp_lib
#endif

IMPLICIT NONE

LOGICAL, INTENT(IN)                      :: ldorhoef
INTEGER, INTENT(IN)                      :: ielast
INTEGER, INTENT(IN)                      :: nsra
INTEGER, INTENT(IN)                      :: nspin
INTEGER, INTENT(IN)                      :: lmax
DOUBLE COMPLEX, INTENT(IN)               :: ez(:)
DOUBLE COMPLEX, INTENT(IN)               :: wez(:)
DOUBLE PRECISION, INTENT(IN)             :: zat
DOUBLE PRECISION, INTENT(IN)             :: socscale
DOUBLE PRECISION, INTENT(IN)             :: cleb(:)
INTEGER, INTENT(IN)                      :: icleb(:,:)
INTEGER, INTENT(IN)                      :: iend
INTEGER, INTENT(IN)                      :: ifunm(:)
INTEGER, INTENT(IN)                      :: lmsp(:)
INTEGER, INTENT(IN)                      :: ncheb
INTEGER, INTENT(IN)                      :: npan_tot
INTEGER, INTENT(IN)                      :: npan_log
INTEGER, INTENT(IN)                      :: npan_eq
DOUBLE PRECISION, INTENT(IN)             :: rmesh(:)
INTEGER, INTENT(IN)                      :: irws
DOUBLE PRECISION, INTENT(IN)             :: rpan_intervall(0:)
INTEGER, INTENT(IN)                      :: ipan_intervall(0:)
DOUBLE PRECISION, INTENT(IN)             :: rnew(:)
DOUBLE PRECISION, INTENT(IN)             :: vinsnew(:,:,:)
DOUBLE PRECISION, INTENT(IN)             :: thetasnew(:,:)
DOUBLE PRECISION, INTENT(INOUT)          :: theta
DOUBLE PRECISION, INTENT(INOUT)          :: phi
INTEGER (kind=1), INTENT(IN)             :: angle_fixed
DOUBLE PRECISION, INTENT(OUT)            :: moment_x
DOUBLE PRECISION, INTENT(OUT)            :: moment_y
DOUBLE PRECISION, INTENT(OUT)            :: moment_z
!INTEGER, INTENT(IN)                      :: i1
INTEGER, INTENT(IN)                      :: ipot
DOUBLE COMPLEX, INTENT(OUT)              :: den_out(0:,:,:)
DOUBLE PRECISION, INTENT(OUT)            :: espv(0:,:)
DOUBLE PRECISION, INTENT(OUT)            :: rho2ns(:,:,:)
DOUBLE PRECISION, INTENT(OUT)            :: r2nef(:,:,:)
DOUBLE COMPLEX, INTENT(IN)               :: gmatn(:,:,:)
DOUBLE PRECISION, INTENT(OUT)            :: muorb(0:,:)
INTEGER, INTENT(IN)                      :: lpotd
INTEGER, INTENT(IN)                      :: lmaxd
INTEGER, INTENT(IN)                      :: irmd
INTEGER, INTENT(IN)                      :: irmd_new
INTEGER, INTENT(IN)                      :: iemxd
LOGICAL, INTENT(IN)                      :: soc
!INCLUDE 'inc.p'


!INTEGER, PARAMETER :: lmmaxd= (lmaxd+1)**2

!INTEGER, PARAMETER :: lmaxd1= lmaxd+1

!INTEGER, PARAMETER :: lmmaxso=2*lmmaxd
!INTEGER :: lmpotd
!DOUBLE PRECISION, INTENT, PARAMETER :: lmpotd= (lpotd+1)**2

!INTEGER, PARAMETER :: lmxspd= (2*lpotd+1)**2

DOUBLE PRECISION, PARAMETER :: cvlight=274.0720442D0
DOUBLE COMPLEX, PARAMETER :: czero=(0D0,0D0)
DOUBLE COMPLEX, PARAMETER :: cone=(1D0,0D0)

!INTEGER, PARAMETER :: nrmaxd=ntotd*(nchebd+1)

INTEGER :: lmmaxd, lmaxd1, lmmaxso, lmpotd, lmxspd, nrmaxd
DOUBLE COMPLEX eryd, ek,df

    
DOUBLE COMPLEX, allocatable :: tmatll(:,:),  &
    tmattemp(:,:)
DOUBLE COMPLEX, allocatable :: gmatll(:,:,:), gmat0(:,:)
INTEGER :: ir,use_sratrick,nvec,lm1,lm2,ie,irmdnew,imt1,  &
    jspin,idim,iorb
DOUBLE PRECISION :: pi,thetanew,phinew
DOUBLE COMPLEX gmatprefactor
DOUBLE PRECISION, allocatable :: vins(:,:,:)
DOUBLE COMPLEX,allocatable :: vnspll0(:,:,:),vnspll1(:,:,:,:), vnspll(:,:,:,:)
DOUBLE COMPLEX, allocatable :: hlk(:,:,:),jlk(:,:,:), hlk2(:,:,:),jlk2(:,:,:)
DOUBLE COMPLEX, allocatable :: rll(:,:,:,:),  &
    rllleft(:,:,:,:),sllleft(:,:,:,:)
DOUBLE COMPLEX, allocatable :: tmatsph(:,:)
DOUBLE COMPLEX, allocatable :: cden(:,:,:,:),  &
    cdenlm(:,:,:,:),cdenns(:,:,:),rho2nsc(:,:,:),r2nefc(:,:,:),  &
    rho2nsnew(:,:,:),r2nefnew(:,:,:),r2orbc(:,:,:,:),  &
    gflle_part(:,:,:),gflle(:,:,:,:),rho2nsc_loop(:,:,:,:), r2nefc_loop(:,:,:,:)

DOUBLE COMPLEX, allocatable:: den(:,:,:,:),denlm(:,:,:,:)
DOUBLE COMPLEX rho2(4),rho2int(4),temp1

DOUBLE COMPLEX rho2ns_temp(2,2),dentemp
DOUBLE PRECISION :: moment(3),totmoment,totxymoment
DOUBLE PRECISION :: denorbmom(3),denorbmomsp(2,4),  &
    denorbmomlm(0:lmaxd,3),denorbmomns(3)
DOUBLE COMPLEX, allocatable :: cdentemp(:,:), rhotemp(:,:),rhonewtemp(:,:)
INTEGER, allocatable :: jlk_index(:)

LOGICAL :: test,opt
EXTERNAL test,opt
!DOUBLE PRECISION :: qvec(:,:)       ! qdos ruess: q-vectors for qdos
!allocatable qvec                 ! qdos ruess
!DOUBLE COMPLEX dentot(2)         ! qdos ruess
!DOUBLE COMPLEX, allocatable :: dentmp(:,:) ! qdos ruess
INTEGER :: iq,nqdos ! qdos ruess: number of qdos points
!INTEGER :: m1,lmshift1(4),lmshift2(4) !, ix       ! qdos ruess
!INTEGER :: lrecgflle,ierr                           ! lmlm-dos
! OMP - number of threads, thread id
INTEGER :: nth,ith

lmmaxd = (lmaxd+1)**2
lmaxd1 = lmaxd+1
lmmaxso = 2*lmmaxd
lmpotd = (lpotd+1)**2
lmxspd = (2*lpotd+1)**2
nrmaxd=irmd_new

allocate(tmatll(lmmaxso,lmmaxso))
allocate(tmattemp(lmmaxso,lmmaxso))
allocate(gmatll(lmmaxso,lmmaxso,iemxd))
allocate(gmat0(lmmaxso,lmmaxso))
!allocate(dentmp(0:lmaxd1,2))
allocate(jlk_index(2*lmmaxso))


! determine if omp is used
       ith = 0
       nth = 1
#ifdef CPP_OMP
!$omp parallel shared(nth,ith)
!$omp single
       nth = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif

pi=4D0*DATAN(1D0)
irmdnew= npan_tot*(ncheb+1)
imt1=ipan_intervall(npan_log+npan_eq)+1
allocate(vins(irmdnew,lmpotd,nspin))
vins=0D0
DO lm1=1,lmpotd
  DO ir=1,irmdnew
    vins(ir,lm1,1)=vinsnew(ir,lm1,ipot)
    vins(ir,lm1,nspin)=vinsnew(ir,lm1,ipot+nspin-1)
  END DO
END DO

!c set up the non-spherical ll' matrix for potential VLL'
IF (NSRA.EQ.2) THEN
USE_SRATRICK=1
ELSE
USE_SRATRICK=0
ENDIF
allocate(vnspll0(lmmaxso,lmmaxso,irmdnew))
allocate(vnspll1(lmmaxso,lmmaxso,irmdnew,0:nth-1))
vnspll0=czero
CALL vllmat(1,irmdnew,lmmaxd,lmmaxso,vnspll0,vins,  &
    cleb,icleb,iend,nspin,zat,rnew,use_sratrick)

! initial allocate
IF (nsra == 2) THEN
  allocate(vnspll(2*lmmaxso,2*lmmaxso,irmdnew,0:nth-1))
ELSE
  allocate(vnspll(lmmaxso,lmmaxso,irmdnew,0:nth-1))
END IF

allocate(hlk(4*(lmax+1),irmdnew,0:nth-1))
allocate(jlk(4*(lmax+1),irmdnew,0:nth-1))
allocate(hlk2(4*(lmax+1),irmdnew,0:nth-1))
allocate(jlk2(4*(lmax+1),irmdnew,0:nth-1))
allocate(tmatsph(2*(lmax+1),0:nth-1))
allocate(rll(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1))
allocate(rllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1))
allocate(sllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1))
allocate(cden(irmdnew,0:lmaxd,4,0:nth-1))
allocate(cdenlm(irmdnew,lmmaxd,4,0:nth-1))
allocate(cdenns(irmdnew,4,0:nth-1))
allocate(rho2nsc(irmdnew,lmpotd,4))
allocate(rho2nsc_loop(irmdnew,lmpotd,4,ielast))
allocate(rho2nsnew(irmd,lmpotd,4))
allocate(r2nefc(irmdnew,lmpotd,4))
allocate(r2nefc_loop(irmdnew,lmpotd,4,0:nth-1))
allocate(r2nefnew(irmd,lmpotd,4))
allocate(r2orbc(irmdnew,lmpotd,4,0:nth-1))
allocate(cdentemp(irmdnew,0:nth-1))
allocate(gflle_part(lmmaxso,lmmaxso,0:nth-1))
allocate(gflle(lmmaxso,lmmaxso,ielast,1))
allocate(den(0:lmaxd1,iemxd,2,1),denlm(lmmaxd,iemxd,2,1))
rho2nsc=czero
rho2nsc_loop=czero
r2nefc=czero
r2nefc_loop=czero
r2orbc=czero
rho2ns=0.d0  ! fivos 19.7.2014, this was CZERO
r2nef=0.d0   ! fivos 19.7.2014, this was CZERO
rho2nsnew=czero
r2nefnew=czero
den=czero
espv=0D0
rho2int=czero
denorbmom=0D0
denorbmomsp=0D0
denorbmomlm=0D0
denorbmomns=0D0
thetanew=0D0
phinew=0D0
gflle_part=czero
gflle=czero
! LM shifts for correct density summation
!lmshift1(1)=0                                                   ! qdos ruess
!lmshift1(2)=lmmaxd                                              ! qdos ruess
!lmshift1(3)=0                                                   ! qdos ruess
!lmshift1(4)=lmmaxd                                              ! qdos ruess
!lmshift2(1)=0                                                   ! qdos ruess
!lmshift2(2)=lmmaxd                                              ! qdos ruess
!lmshift2(3)=lmmaxd                                              ! qdos ruess
!lmshift2(4)=0                                                   ! qdos ruess

GMAT0 = czero
gmatll = czero

DO ir=1,3
  DO lm1=0,lmaxd1+1
    muorb(lm1,ir)=0D0
  END DO
END DO

 nqdos = 1                                                         ! qdos ruess
!IF (opt('qdos    ')) THEN                                         ! qdos ruess
!        Read BZ path for qdos calculation:                             ! qdos ruess
!  OPEN(67,FILE='qvec.dat',STATUS='old',IOSTAT=ierr,ERR=3000)     ! qdos ruess
!  READ(67,*) nqdos                                               ! qdos ruess
!  allocate(qvec(3,nqdos))                                        ! qdos ruess
!  DO iq = 1,nqdos                                                ! qdos ruess
!    READ(67,*) (qvec(ix,iq),ix=1,3)                             ! qdos ruess
!  END DO                                                          ! qdos ruess
!  CLOSE(67)                                                      ! qdos ruess
!        Change allocation for GFLLE to be suitabel for qdos run        ! qdos ruess
!  deallocate(gflle,den,denlm)                                    ! qdos ruess
!  allocate(gflle(lmmaxso,lmmaxso,ielast,nqdos))                  ! qdos ruess
!  allocate(den(0:lmaxd1,iemxd,2,nqdos), denlm(lmmaxd,iemxd,2,nqdos))
!  3000  IF (ierr /= 0) STOP 'ERROR READING ''QVEC.DAT'''                  ! QDOS Ruess
!END IF  ! OPT('qdos    ')                                         ! qdos ruess

!IF ((opt('lmlm-dos')).AND.(i1 == 1)) THEN                         ! lmlm-dos ruess
!  lrecgflle = 4*lmmaxso*lmmaxso*ielast*nqdos                     ! lmlm-dos ruess
!  OPEN(91,ACCESS='direct',RECL=lrecgflle,FILE='gflle',           ! lmlm-dos ruess  &
!      FORM='unformatted',STATUS='replace',ERR=3001,IOSTAT=ierr)! lmlm-dos ruess
!  3001 IF (ierr /= 0) STOP 'ERROR CREATING ''GFLLE'''                    ! LMLM-DOs ruess
!END IF                                                             ! lmlm-dos ruess

! energy loop
!WRITE(6,*) 'atom: ',i1


#ifdef CPP_OMP
! omp: start parallel region here
!$omp parallel do default(none) ,&
!$omp& private(eryd,ie,ir,lm1,lm2,gmatprefactor,nvec) ,&
!$omp& private(jlk_index,tmatll,ith) ,&
!$omp& shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmax) ,&
!$omp& shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew) ,&
!$omp& shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop) ,&
!$omp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0) ,&
!$omp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,cdentemp) ,&
!$omp& shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft) ,&
!$omp& private(iq,df,ek,tmattemp,gmatll,gmat0,iorb,dentemp) ,&
!$omp& private(rho2ns_temp,rho2,temp1,jspin) ,&
!$omp& shared(ldorhoef,nqdos,wez,lmsp,imt1,ifunm) ,&
!$omp& shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop) ,&
!$omp& shared(lmaxd,lmaxd1,lmmaxd,lmpotd,nrmaxd,soc,lmmaxso,gmatn) ,&
!$omp& reduction(+:rho2int,espv) reduction(-:muorb) ,&
!$omp& reduction(-:denorbmom,denorbmomsp,denorbmomlm,denorbmomns)
#endif

DO ie=1,ielast

#ifdef CPP_OMP
    ith = omp_get_thread_num()
#else
        ith = 0
#endif
  
  eryd=ez(ie)
  ek=SQRT(eryd)
  df=wez(ie)/DBLE(nspin)
  IF (nsra == 2) ek = SQRT( eryd + eryd*eryd/(cvlight*cvlight) ) *  &
      ( 1D0 + eryd/(cvlight*cvlight) )
!!$noomp critical
!  WRITE(6,*) 'energy:',ie,'',eryd
!!$noomp end critical
!
!         IREC=IE+IELAST*(I1-1)
!         READ(69,REC=IREC) GMAT0
!
!  rotate gmat from global frame to local frame
        
!       GMAT0 = gmatn(:,:,ie)
!       CALL ROTATEMATRIX(GMAT0,THETA,PHI,LMMAXD,1)

!         DO LM1=1,LMMAXSO
!          DO LM2=1,LMMAXSO
!           GMATLL(LM1,LM2,IE)=GMAT0(LM1,LM2)
!          ENDDO
!         ENDDO
  
! recalculate wavefuntions, also include left solution
! contruct the spin-orbit coupling hamiltonian and add to potential
  CALL spinorbit_ham(lmax,lmmaxd,vins,rnew,  &
      eryd,zat,cvlight,socscale,nspin,lmpotd,  &
      theta,phi,ipan_intervall,rpan_intervall, npan_tot,ncheb,irmdnew,nrmaxd,  &
      vnspll0,vnspll1(:,:,:,ith),'1',soc)
  
!c extend matrix for the SRA treatment
  vnspll(:,:,:,ith)=czero
  IF (nsra == 2) THEN
    IF (use_sratrick == 0) THEN
      CALL vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,  &
          lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'Ref=0')
    ELSE IF (use_sratrick == 1) THEN
      CALL vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,  &
          lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'Ref=Vsph')
    END IF
  ELSE
    vnspll(:,:,:,ith)=vnspll1(:,:,:,ith)
  END IF
  
!c calculate the source terms in the Lippmann-Schwinger equation
!c these are spherical hankel and bessel functions
  hlk(:,:,ith)=czero
  jlk(:,:,ith)=czero
  hlk2(:,:,ith)=czero
  jlk2(:,:,ith)=czero
  gmatprefactor=czero
  jlk_index=0
  CALL rllsllsourceterms(nsra,nvec,eryd,rnew,irmdnew,nrmaxd,lmax,  &
      lmmaxso,1,jlk_index,hlk(:,:,ith),  &
      jlk(:,:,ith),hlk2(:,:,ith),jlk2(:,:,ith), gmatprefactor)
  
! using spherical potential as reference
  IF (use_sratrick == 1) THEN
    CALL calcsph(nsra,irmdnew,nrmaxd,lmax,nspin,zat,cvlight,eryd,  &
        rnew,vins,ncheb,npan_tot,rpan_intervall,  &
        jlk_index,hlk(:,:,ith),jlk(:,:,ith),hlk2(:,:,ith),  &
        jlk2(:,:,ith),gmatprefactor,tmatsph(:,ith), use_sratrick)
  END IF
  
!c calculate the tmat and wavefunctions
  rllleft(:,:,:,ith)=czero
  sllleft(:,:,:,ith)=czero
  
!c right solutions
  tmatll=czero
  CALL rll_only(rpan_intervall,rnew,vnspll(:,:,:,ith),  &
      rll(:,:,:,ith),tmatll,  &
      ncheb,npan_tot,lmmaxso,nvec*lmmaxso,4*(lmax+1),  &
      irmdnew,nsra,jlk_index,hlk(:,:,ith),jlk(:,:,ith),  &
      hlk2(:,:,ith),jlk2(:,:,ith), gmatprefactor,'1','1',use_sratrick)
  IF (nsra == 2) THEN
    rll(lmmaxso+1:nvec*lmmaxso,:,:,ith)=  &
        rll(lmmaxso+1:nvec*lmmaxso,:,:,ith)/cvlight
  END IF
  
! left solutions
! contruct the TRANSPOSE spin-orbit coupling hamiltonian and add to potential
  CALL spinorbit_ham(lmax,lmmaxd,vins,rnew,eryd,zat,  &
      cvlight,socscale,nspin,lmpotd,theta,phi,  &
      ipan_intervall,rpan_intervall,npan_tot,ncheb,  &
      irmdnew,nrmaxd,vnspll0,vnspll1(:,:,:,ith), 'transpose',soc)
  
!c extend matrix for the SRA treatment
  vnspll(:,:,:,ith)=czero
  IF (nsra == 2) THEN
    IF (use_sratrick == 0) THEN
      CALL vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,  &
          lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'Ref=0')
    ELSE IF (use_sratrick == 1) THEN
      CALL vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,  &
          lmmaxso,irmdnew,nrmaxd,eryd,cvlight,lmax,0,'Ref=Vsph')
    END IF
  ELSE
    vnspll(:,:,:,ith)=vnspll1(:,:,:,ith)
  END IF
  
!c calculate the source terms in the Lippmann-Schwinger equation
!c these are spherical hankel and bessel functions
  hlk(:,:,ith)=czero
  jlk(:,:,ith)=czero
  hlk2(:,:,ith)=czero
  jlk2(:,:,ith)=czero
  gmatprefactor=czero
  jlk_index=0
  CALL rllsllsourceterms(nsra,nvec,eryd,rnew,irmdnew,nrmaxd,lmax,  &
      lmmaxso,1,jlk_index,hlk(:,:,ith),  &
      jlk(:,:,ith),hlk2(:,:,ith),jlk2(:,:,ith), gmatprefactor)
  
!c using spherical potential as reference
! notice that exchange the order of left and right hankel/bessel functions
  IF (use_sratrick == 1) THEN
    CALL calcsph(nsra,irmdnew,nrmaxd,lmax,nspin,zat,cvlight,eryd,  &
        rnew,vins,ncheb,npan_tot,rpan_intervall,  &
        jlk_index,hlk2(:,:,ith),jlk2(:,:,ith),  &
        hlk(:,:,ith),jlk(:,:,ith),gmatprefactor, tmatsph(:,ith),use_sratrick)
  END IF
  
!c calculate the tmat and wavefunctions
  rllleft(:,:,:,ith)=czero
  sllleft(:,:,:,ith)=czero
  
!c left solutions
! notice that exchange the order of left and right hankel/bessel functions
  tmattemp=czero
  CALL rllsll(rpan_intervall,rnew,vnspll(:,:,:,ith),  &
      rllleft(:,:,:,ith),sllleft(:,:,:,ith),tmattemp,  &
      ncheb,npan_tot,lmmaxso,nvec*lmmaxso,4*(lmax+1),  &
      irmdnew,nsra,jlk_index,hlk2(:,:,ith),jlk2(:,:,ith),  &
      hlk(:,:,ith),jlk(:,:,ith), gmatprefactor,'1','1',use_sratrick)
  IF (nsra == 2) THEN
    rllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)=  &
        rllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)/cvlight
    sllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)=  &
        sllleft(lmmaxso+1:nvec*lmmaxso,:,:,ith)/cvlight
  END IF
  
  DO  iq = 1,nqdos                                       ! qdos
    den(:,ie,:,iq)=czero
! read in gf
!    irec = iq + nqdos * (ie-1) +  nqdos * ielast * (i1-1)     ! qdos
!!$noomp critical
!    READ(69,REC=irec) gmat0
!!$noomp end critical
   
     GMAT0 = gmatn(:,:,ie)
! rotate gmat from global frame to local frame
    CALL rotatematrix(gmat0,theta,phi,lmmaxd,1)
    
    DO lm1=1,lmmaxso
      DO lm2=1,lmmaxso
        gmatll(lm1,lm2,ie)=gmat0(lm1,lm2)
      END DO
    END DO
! calculate density
    CALL rhooutnew(nsra,lmmaxd,lmmaxso,lmax,gmatll(:,:,ie),ek,  &
        df,cleb,icleb,iend,  &
        irmdnew,thetasnew,ifunm,imt1,lmsp,  &
        rll(:,:,:,ith), rllleft(:,:,:,ith),sllleft(:,:,:,ith),  &
        cden(:,:,:,ith),cdenlm(:,:,:,ith),  &
        cdenns(:,:,ith),rho2nsc_loop(:,:,:,ie),0,  &
        lmaxd)
    
    DO jspin=1,4
      
      DO lm1 = 0,lmax
        cdentemp(:,ith)=czero
        dentemp=czero
        DO ir=1,irmdnew
          cdentemp(ir,ith)=cden(ir,lm1,jspin,ith)
        END DO
        CALL intcheb_cell(cdentemp(:,ith),dentemp,  &
            rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
        rho2(jspin)=dentemp
        rho2int(jspin)=rho2int(jspin)+rho2(jspin)*df
        IF (jspin <= 2) THEN
          den(lm1,ie,jspin,iq)=rho2(jspin)
        END IF
      END DO
      
      IF (jspin <= 2) THEN
        DO lm1 = 1,lmmaxd
          cdentemp(:,ith)=czero
          dentemp=czero
          DO ir=1,irmdnew
            cdentemp(ir,ith)=cdenlm(ir,lm1,jspin,ith)
          END DO
          CALL intcheb_cell(cdentemp(:,ith),dentemp,  &
              rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
          denlm(lm1,ie,jspin,iq)=dentemp
        END DO
        cdentemp(:,ith)=czero
        dentemp=czero
        DO ir=1,irmdnew
          cdentemp(ir,ith)=cdenns(ir,jspin,ith)
        END DO
        CALL intcheb_cell(cdentemp(:,ith),dentemp,  &
            rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
        den(lmaxd1,ie,jspin,iq)=dentemp
        rho2int(jspin)=rho2int(jspin)+den(lmaxd1,ie,jspin,iq)*df
      END IF
    END DO ! JSPIN
    
    DO jspin=1,4
      IF (jspin <= 2) THEN
        DO lm1=0,lmaxd1
          espv(lm1,jspin)=espv(lm1,jspin)+  &
              DIMAG( eryd * den(lm1,ie,jspin,iq) * df )
        END DO
      END IF
    END DO
  END DO   ! IQ = 1,NQDOS
!END DO

! get charge at the Fermi energy (IELAST)

IF (ie == ielast.AND.ldorhoef) THEN
  CALL rhooutnew(nsra,lmmaxd,lmmaxso,lmax,gmatll(:,:,ie),ek,  &
      cone,cleb,icleb,iend,  &
      irmdnew,thetasnew,ifunm,imt1,lmsp,  &
      rll(:,:,:,ith), rllleft(:,:,:,ith),sllleft(:,:,:,ith),  &
      cden(:,:,:,ith),cdenlm(:,:,:,ith),  &
      cdenns(:,:,ith),r2nefc_loop(:,:,:,ith),0,  &
      lmaxd)
END IF


! get orbital moment
DO iorb=1,3
  CALL rhooutnew(nsra,lmmaxd,lmmaxso,lmax,gmatll(:,:,ie),ek,  &
      cone,cleb,icleb,iend,  &
      irmdnew,thetasnew,ifunm,imt1,lmsp,  &
      rll(:,:,:,ith), rllleft(:,:,:,ith),sllleft(:,:,:,ith),  &
      cden(:,:,:,ith),cdenlm(:,:,:,ith),  &
      cdenns(:,:,ith),r2orbc(:,:,:,ith),iorb,  &
      lmaxd)
  DO jspin=1,4
    IF (jspin <= 2) THEN
      DO lm1=0,lmax
        cdentemp(:,ith)=czero
        dentemp=czero
        DO ir=1,irmdnew
          cdentemp(ir,ith)=cden(ir,lm1,jspin,ith)
        END DO
        CALL intcheb_cell(cdentemp(:,ith),dentemp,rpan_intervall,  &
            ipan_intervall,npan_tot,ncheb,irmdnew)
        rho2(jspin)=dentemp
        muorb(lm1,jspin)=muorb(lm1,jspin)-DIMAG(rho2(jspin)*df)
        denorbmom(iorb)=denorbmom(iorb)-DIMAG(rho2(jspin)*df)
        denorbmomsp(jspin,iorb)=denorbmomsp(jspin,iorb)- DIMAG(rho2(jspin)*df)
        denorbmomlm(lm1,iorb)=denorbmomlm(lm1,iorb)- DIMAG(rho2(jspin)*df)
        cdentemp(:,ith)=czero
        DO ir=1,irmdnew
          cdentemp(ir,ith)=cdenns(ir,jspin,ith)
        END DO
        CALL intcheb_cell(cdentemp(:,ith),temp1,  &
            rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
        denorbmomns(iorb)=denorbmomns(iorb)-DIMAG(temp1*df)
      END DO
    END IF
  END DO
END DO ! IORB
END DO ! IE loop

#ifdef CPP_OMP
!$omp end parallel do
#endif
! omp: move sum from rhooutnew here after parallel calculation
DO ir=1,irmdnew
  DO lm1=1,lmpotd
    DO jspin=1,4
      DO ie=1,ielast
        rho2nsc(ir,lm1,jspin) = rho2nsc(ir,lm1,jspin) +  &
            rho2nsc_loop(ir,lm1,jspin,ie)
      END DO
    END DO
  END DO
END DO
! omp: don't forget to do the same with density at fermi energy:
DO ith=0,nth-1
  r2nefc(:,:,:) = r2nefc(:,:,:) + r2nefc_loop(:,:,:,ith)
END DO

! omp: moved write-out of dos files out of parallel energy loop
! Write out qdos and lm-dos:                                            ! lm-dos
!DO ie=1,ielast                                                   ! lm-dos
!  DO iq=1,nqdos                                                    ! lm-dos
!    IF ((iq == 1).AND.(ie == 1)) THEN                                ! lm-dos
!      OPEN(29,                                                     ! lm-dos  &
!          FILE="lmdos."//CHAR(48+i1/10)//CHAR(48+MOD(i1,10))//"."//    ! lm-dos  &
!          CHAR(48+1)//".dat")                                          ! lm-dos
!     WRITE (29,*) ' '                                             ! lm-dos
!      WRITE (29,8600) '# ISPIN=',1,' I1=',i1                       ! lm-dos
!      OPEN(30,                                                     ! lm-dos  &
!          FILE="lmdos."//CHAR(48+i1/10)//CHAR(48+MOD(i1,10))//"."//    ! lm-dos  &
!          CHAR(48+2)//".dat")                                          ! lm-dos
!      WRITE (30,*) ' '                                             ! lm-dos
!      WRITE (30,8600) '# ISPIN=',2,' I1=',i1                       ! lm-dos
!    END IF                                                            ! lm-dos
    
!    IF (opt('qdos    ')) THEN                                        ! qdos ruess
!      IF ((iq == 1).AND.(ie == 1)) THEN                              ! qdos ruess
!        OPEN(31,                                                    ! qdos ruess  &
!            FILE="qdos."//CHAR(48+i1/10)//CHAR(48+MOD(i1,10))//"."//  ! qdos ruess  &
!            CHAR(48+1)//".dat")                                       ! qdos ruess
!        WRITE (31,*) ' '                                            ! qdos ruess
!        WRITE (31,8600) '# ISPIN=',1,' I1=',i1                      ! qdos ruess
!        WRITE(31,'(7(A,3X))') '#   Re(E)','Im(E)','k_x','k_y','k_z',! qdos  &
!            'DEN_tot','DEN_s,p,...'                     ! qdos
!        OPEN(32,                                                    ! qdos ruess  &
!            FILE="qdos."//CHAR(48+i1/10)//CHAR(48+MOD(i1,10))//"."//  ! qdos ruess  &
!            CHAR(48+2)//".dat")                                       ! qdos ruess
!        WRITE (32,*) ' '                                            ! qdos ruess
!        WRITE (32,8600) '# ISPIN=',2,' I1=',i1                      ! qdos ruess
!        WRITE(32,'(7A)') '#   Re(E)','Im(E)','k_x','k_y','k_z',     ! qdos  &
!            'DEN_tot','DEN_s,p,...'                     ! qdos
        
!        8600 FORMAT (a8,i3,a4,i5)                                              ! qdos ruess
!      END IF   ! IQ.EQ.1                                              ! qdos ruess
!      DO jspin =1,2                                               ! qdos ruess
!        dentot(jspin) = DCMPLX(0.d0,0.d0)                         ! qdos ruess
!        DO l1 = 0,lmaxd1                                          ! qdos ruess
!          dentot(jspin) = dentot(jspin) + den(l1,ie,1,iq)         ! qdos ruess
!        END DO                                                     ! qdos ruess
!      END DO                                                       ! qdos ruess
!    write qdos.nn.s.dat                                                ! qdos ruess
!    and lmdos.nn.s.dat                                                 ! qdos ruess
!      WRITE(29,9000) ez(ie),qvec(1,iq),qvec(2,iq),qvec(3,iq),     ! qdos ruess  &
!          (-DIMAG(denlm(l1,ie,1,iq))/pi,l1=1,lmmaxd)         ! qdos ruess
!      WRITE(30,9000) ez(ie),qvec(1,iq),qvec(2,iq),qvec(3,iq),     ! qdos ruess  &
!          (-DIMAG(denlm(l1,ie,2,iq))/pi,l1=1,lmmaxd)         ! qdos ruess
!      WRITE(31,9000) ez(ie),qvec(1,iq),qvec(2,iq),qvec(3,iq),     ! qdos ruess  &
!          -DIMAG(dentot(1))/pi,(-DIMAG(den(l1,ie,1,iq))/pi,l1=0,lmaxd1)! qdos ruess
!      WRITE(32,9000) ez(ie),qvec(1,iq),qvec(2,iq),qvec(3,iq),     ! qdos ruess  &
!          -DIMAG(dentot(2))/pi,(-DIMAG(den(l1,ie,2,iq))/pi,l1=0,lmaxd1)! qdos ruess
!    ELSE                                                             ! lm-dos
!      WRITE(29,9001) ez(ie),                                      ! lm-dos  &
!          (-DIMAG(denlm(l1,ie,1,iq))/pi,l1=1,lmmaxd)         ! lm-dos
!      WRITE(30,9001) ez(ie),                                      ! lm-dos  &
!          (-DIMAG(denlm(l1,ie,2,iq))/pi,l1=1,lmmaxd)         ! lm-dos
!      9001       FORMAT(30E12.4)                                             ! lm-dos
!    END IF      ! OPT('qdos    ')                                     ! qdos ruess
!    9000      FORMAT(5F10.6,40E16.8)                                       ! qdos ruess
!  END DO !IQ
!END DO !IE

! write
!IF (opt('lmlm-dos')) THEN                                         ! lmlm-dos
!          DO JSPIN = 1,2                                                  ! lmlm-dos
!           OPEN(90,                                                       ! lmlm-dos
!      &    FILE="lmlmdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//    ! lmlm-dos
!      &                                          char(48+JSPIN)//".dat")  ! lmlm-dos
!           DO IE = 1,IELAST                                               ! lmlm-dos
!              DO LM1 = 1,LMMAXD                                           ! lmlm-dos
!                 IF (.NOT.(OPT('qdos    '))) THEN                                   ! qdos
!                    WRITE(90,1000) EZ(IE),                                ! lmlm-dos
!      &                            (-DIMAG(GFLLE(LM1+LMSHIFT1(JSPIN),     ! lmlm-dos
!      &                     LM2+LMSHIFT2(JSPIN),IE,1))/PI,LM2 = 1,LMMAXD) ! lmlm-dos
!                 ELSE                                                               ! qdos
!                   DO IQ=1,NQDOS                                                    ! qdos
!                    WRITE(90,1000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),                    ! qdos
!      &                            QVEC(3,IQ),(-DIMAG(GFLLE(LM1+                    ! qdos
!      &                            LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),             ! qdos
!      &                            IE,IQ))/PI,LM2 = 1,LMMAXD)                       ! qdos
!                   ENDDO ! IQ=1,NQDOS                                               ! qdos
!                 ENDIF                                                              ! qdos
!              ENDDO                                                       ! lmlm-dos
!           ENDDO !IE                                                      ! lmlm-dos
!           CLOSE(90)                                                      ! lmlm-dos
!          ENDDO !JSPIN                                                    ! lmlm-dos
!  1000  FORMAT(5F10.6,I3,40E16.8)                                         ! lmlm-dos
! write gflle to file                                                    ! lmlm-dos
!  WRITE(91,REC=i1) gflle                                            ! lmlm-dos
!END IF                                                             ! lmlm-dos

allocate(rhotemp(irmdnew,lmpotd))
allocate(rhonewtemp(irws,lmpotd))
DO jspin=1,4
  rhotemp=czero
  rhonewtemp=czero
  DO lm1=1,lmpotd
    DO ir=1,irmdnew
      rhotemp(ir,lm1)=rho2nsc(ir,lm1,jspin)
    END DO
  END DO
  CALL cheb2oldgrid(irws,irmdnew,lmpotd,rmesh,ncheb,npan_tot,  &
      rpan_intervall,ipan_intervall, rhotemp,rhonewtemp,irmd)
  DO lm1=1,lmpotd
    DO ir=1,irws
      rho2nsnew(ir,lm1,jspin)=rhonewtemp(ir,lm1)
    END DO
  END DO
  
  rhotemp=czero
  rhonewtemp=czero
  DO lm1=1,lmpotd
    DO ir=1,irmdnew
      rhotemp(ir,lm1)=r2nefc(ir,lm1,jspin)
    END DO
  END DO
  CALL cheb2oldgrid(irws,irmdnew,lmpotd,rmesh,ncheb,npan_tot,  &
      rpan_intervall,ipan_intervall, rhotemp,rhonewtemp,irmd)
  DO lm1=1,lmpotd
    DO ir=1,irws
      r2nefnew(ir,lm1,jspin)=rhonewtemp(ir,lm1)
    END DO
  END DO
END DO
deallocate(rhotemp)
deallocate(rhonewtemp)
! calculate new THETA and PHI for non-colinear
!IF (.NOT.test('FIXMOM  ')) THEN
if (angle_fixed == 0) then ! angle not fixed
  rho2ns_temp(1,1)=rho2int(1)
  rho2ns_temp(2,2)=rho2int(2)
  rho2ns_temp(1,2)=rho2int(3)
  rho2ns_temp(2,1)=rho2int(4)
  
  CALL rotatematrix(rho2ns_temp,theta,phi,1,0)
  
  rho2int(1)=rho2ns_temp(1,1)
  rho2int(2)=rho2ns_temp(2,2)
  rho2int(3)=rho2ns_temp(1,2)
  rho2int(4)=rho2ns_temp(2,1)
  
  
  moment(1)=DIMAG(rho2int(3)+rho2int(4))
  moment(2)=-REAL(rho2int(3)-rho2int(4))
  moment(3)=DIMAG(-rho2int(1)+rho2int(2))

  moment_x=moment(1)
  moment_y=moment(2)
  moment_z=moment(3)
  
  totmoment=SQRT(moment(1)**2+moment(2)**2+moment(3)**2)
  totxymoment=SQRT(moment(1)**2+moment(2)**2)
  
  IF (ABS(totxymoment) > 1D-05) THEN
    IF (ABS(moment(3)) < 1D-05) THEN
      thetanew=pi/2D0
    ELSE
      thetanew=ACOS(moment(3)/totmoment)
    END IF
    IF (totxymoment < 1D-05) THEN
      phinew=0D0
    ELSE
      phinew=DATAN2(moment(2),moment(1))
    END IF
  END IF

  ! UPDATE ANGLES
!  phi   = phinew
!  theta = thetanew

  !          THETANEW=ACOS(MOMENT(3)/TOTMOMENT)
!          PHINEW=DATAN2(MOMENT(2),MOMENT(1))
!  WRITE(6,*) 'moment',moment(1),moment(2),moment(3)
!        WRITE(6,*) 'total moment',TOTMOMENT,TOTXYMOMENT
!  WRITE(6,*) 'angles', thetanew,phinew
!  WRITE(11,*) thetanew,phinew
!  WRITE(12,*) thetanew,phinew

! Use old angles for rotation
!if (angle_fixed == 1) then
!  thetanew = theta
!  phinew   = phi
!endif 

  CALL rotatevector(rho2nsnew,rho2ns,irws,lmpotd,thetanew,phinew,  &
      theta,phi,irmd)
  CALL rotatevector(r2nefnew,r2nef,irws,lmpotd,thetanew,phinew,  &
      theta,phi,irmd)

else ! angle fixed

  rho2ns_temp(1,1)=rho2int(1)
  rho2ns_temp(2,2)=rho2int(2)
  rho2ns_temp(1,2)=rho2int(3)
  rho2ns_temp(2,1)=rho2int(4)
  
  CALL rotatematrix(rho2ns_temp,theta,phi,1,0)
  
  rho2int(1)=rho2ns_temp(1,1)
  rho2int(2)=rho2ns_temp(2,2)
  rho2int(3)=rho2ns_temp(1,2)
  rho2int(4)=rho2ns_temp(2,1)

  moment(1)=DIMAG(rho2int(3)+rho2int(4))
  moment(2)=-REAL(rho2int(3)-rho2int(4))
  moment(3)=DIMAG(-rho2int(1)+rho2int(2))

  moment_x=moment(1)
  moment_y=moment(2)
  moment_z=moment(3)
  
  rho2ns(:,:,:)=DIMAG(rho2nsnew(:,:,:))
  r2nef(:,:,:)=DIMAG(r2nefnew(:,:,:))
endif

idim = irmd*lmpotd
CALL dscal(idim,2.d0,rho2ns(1,1,1),1)
CALL daxpy(idim,-0.5D0,rho2ns(1,1,1),1,rho2ns(1,1,2),1)
CALL daxpy(idim,1.0D0,rho2ns(1,1,2),1,rho2ns(1,1,1),1)

! --> do the same at the Fermi energy

CALL dscal(idim,2.d0,r2nef(1,1,1),1)
CALL daxpy(idim,-0.5D0,r2nef(1,1,1),1,r2nef(1,1,2),1)
CALL daxpy(idim,1.0D0,r2nef(1,1,2),1,r2nef(1,1,1),1)

DO lm1=0,lmaxd1
  DO ie=1,iemxd
    DO jspin=1,nspin
      den_out(lm1,ie,jspin) =  den(lm1,ie,jspin,1)
    END DO
  END DO
END DO

! UPDATE ANGLES
if (angle_fixed == 0) then
phi   = phinew
theta = thetanew        
endif

deallocate(vins)
deallocate(vnspll0)
deallocate(vnspll1)
deallocate(vnspll)
deallocate(hlk)
deallocate(jlk)
deallocate(hlk2)
deallocate(jlk2)
deallocate(tmatsph)
deallocate(rll)
deallocate(rllleft)
deallocate(sllleft)
deallocate(cden)
deallocate(cdenlm)
deallocate(cdenns)
deallocate(rho2nsc,rho2nsc_loop)
deallocate(rho2nsnew)
deallocate(r2nefc,r2nefc_loop)
deallocate(r2nefnew)
deallocate(r2orbc)
deallocate(cdentemp)
deallocate(gflle_part)
deallocate(gflle)
deallocate(den,denlm)
END SUBROUTINE rhovalnew

SUBROUTINE rhooutnew(nsra,lmmaxd,lmmaxso,lmax,gmatll,ek,  &
        df,cleb,icleb,iend,  &
        irmdnew,thetasnew,ifunm,imt1,  &
        lmsp,rll,rllleft,sllleft,  &
        cden,cdenlm,cdenns,rho2nsc,corbital,  &
        lmaxd)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-04-21  Time: 16:24:21

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nsra
INTEGER, INTENT(IN)                      :: lmmaxd
INTEGER, INTENT(IN)                      :: lmmaxso
INTEGER, INTENT(IN)                      :: lmax
DOUBLE COMPLEX, INTENT(IN)               :: gmatll(:,:)
DOUBLE COMPLEX, INTENT(IN)               :: ek
!INTEGER, INTENT(IN)                      :: lmpotd
DOUBLE COMPLEX, INTENT(IN)               :: df
!INTEGER, INTENT(IN)                      :: npan_tot
!INTEGER, INTENT(IN)                      :: ncheb
DOUBLE PRECISION, INTENT(IN)             :: cleb(:)
INTEGER, INTENT(IN)                      :: icleb(:,:)
INTEGER, INTENT(IN)                      :: iend
INTEGER, INTENT(IN)                      :: irmdnew
!INTEGER, INTENT(IN)                      :: nrmaxd
DOUBLE PRECISION, INTENT(IN)             :: thetasnew(:,:)
INTEGER, INTENT(IN)                      :: ifunm(:)
!DOUBLE PRECISION, INTENT(IN)             :: rnew(:)
INTEGER, INTENT(IN)                      :: imt1
INTEGER, INTENT(IN)                      :: lmsp(:)
DOUBLE COMPLEX, INTENT(IN)               :: rll(:,:,:)
!DOUBLE COMPLEX, INTENT(IN)               :: sll(:,:,:)
DOUBLE COMPLEX, INTENT(IN)               :: rllleft(:,:,:)
DOUBLE COMPLEX, INTENT(IN)               :: sllleft(:,:,:)
DOUBLE COMPLEX, INTENT(OUT)              :: cden(:,0:,:)
DOUBLE COMPLEX, INTENT(OUT)              :: cdenlm(:,:,:)
DOUBLE COMPLEX, INTENT(OUT)              :: cdenns(:,:)
DOUBLE COMPLEX, INTENT(OUT)              :: rho2nsc(:,:,:)
INTEGER, INTENT(IN)                      :: corbital
!DOUBLE COMPLEX, INTENT(OUT)              :: gflle_part(:,:)
!DOUBLE PRECISION, INTENT(IN)             :: rpan_intervall(:)
!INTEGER, INTENT(IN)                      :: ipan_intervall(:)
INTEGER, INTENT(IN)                      :: lmaxd  ! new parameter 

!INCLUDE 'inc.p'


DOUBLE COMPLEX, PARAMETER :: czero=(0D0,0D0)
DOUBLE COMPLEX, PARAMETER :: cone=(1D0,0D0)
DOUBLE COMPLEX cltdf

INTEGER :: ir,jspin,lm1,lm2,lm3,m1,l1,j,ifun
DOUBLE PRECISION :: c0ll

DOUBLE COMPLEX, allocatable :: wr(:,:,:),qnsi(:,:),pnsi(:,:),  &
    cwr(:)                                        ! lmlm-dos
INTEGER :: lmshift1(4),lmshift2(4)
DOUBLE COMPLEX, allocatable :: loperator(:,:,:)
LOGICAL :: test,opt
EXTERNAL test,opt
allocate(wr(lmmaxso,lmmaxso,irmdnew))
allocate(cwr(irmdnew))
allocate(qnsi(lmmaxso,lmmaxso))
allocate(pnsi(lmmaxso,lmmaxso))
allocate(loperator(lmmaxso,lmmaxso,3))

wr=czero
cwr=czero
qnsi=czero
pnsi=czero
! set LMSHIFT value which is need to construct CDEN
lmshift1(1)=0
lmshift1(2)=lmmaxd
lmshift1(3)=0
lmshift1(4)=lmmaxd
lmshift2(1)=0
lmshift2(2)=lmmaxd
lmshift2(3)=lmmaxd
lmshift2(4)=0

! for orbital moment
IF (corbital /= 0) THEN
  CALL calc_orbitalmoment(lmaxd,lmmaxso,loperator)
END IF

c0ll=1D0/SQRT(16D0*ATAN(1D0))
cden=czero
cdenlm=czero

DO ir = 1,irmdnew
  
  DO lm1 = 1,lmmaxso
    DO lm2 = 1,lmmaxso
      qnsi(lm1,lm2)=sllleft(lm1,lm2,ir)
!          PNSI(LM1,LM2)=RLL(LM1,LM2,IR)
      pnsi(lm1,lm2)=rllleft(lm1,lm2,ir)
    END DO
  END DO
!        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
!     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
  CALL zgemm('N','T',lmmaxso,lmmaxso,lmmaxso,cone,pnsi,  &
      lmmaxso,gmatll,lmmaxso,ek,qnsi,lmmaxso)
  DO lm1 = 1,lmmaxso
    DO lm2 = 1,lmmaxso
      pnsi(lm1,lm2)=rll(lm1,lm2,ir)
    END DO
  END DO
  CALL zgemm('N','T',lmmaxso,lmmaxso,lmmaxso,cone,pnsi,  &
      lmmaxso,qnsi,lmmaxso,czero,wr(1,1,ir),lmmaxso)
  
  IF (nsra == 2) THEN
    DO lm1 = 1,lmmaxso
      DO lm2 = 1,lmmaxso
!          QNSI(LM1,LM2)=SLLLEFT(LM1+LMMAXSO,LM2,IR)
        qnsi(lm1,lm2)=-sllleft(lm1+lmmaxso,lm2,ir)
!          PNSI(LM1,LM2)=RLLLEFT(LM1+LMMAXSO,LM2,IR)
        pnsi(lm1,lm2)=-rllleft(lm1+lmmaxso,lm2,ir)
      END DO
    END DO
!        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
!     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
    CALL zgemm('N','T',lmmaxso,lmmaxso,lmmaxso,cone,pnsi,  &
        lmmaxso,gmatll,lmmaxso,ek,qnsi,lmmaxso)
    DO lm1 = 1,lmmaxso
      DO lm2 = 1,lmmaxso
        pnsi(lm1,lm2)=rll(lm1+lmmaxso,lm2,ir)
      END DO
    END DO
    CALL zgemm('N','T',lmmaxso,lmmaxso,lmmaxso,cone,pnsi,  &
        lmmaxso,qnsi,lmmaxso,cone,wr(1,1,ir),lmmaxso)
  END IF
  
! for orbital moment
  IF (corbital /= 0) THEN
    CALL zgemm('N','N',lmmaxso,lmmaxso,lmmaxso,cone,  &
        loperator(1,1,corbital),lmmaxso,wr(1,1,ir), lmmaxso,czero,pnsi,lmmaxso)
    DO lm1=1,lmmaxso
      DO lm2=1,lmmaxso
        wr(lm1,lm2,ir)=pnsi(lm1,lm2)
      END DO
    END DO
  END IF
  
  DO jspin = 1,4
    DO lm1 = 1,lmmaxd
      DO lm2 = 1,lm1-1
        wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir)=  &
            wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir)+  &
            wr(lm2+lmshift1(jspin),lm1+lmshift2(jspin),ir)
      END DO
    END DO
  END DO ! JSPIN
  
END DO !IR




!IF (opt('lmlm-dos')) THEN                                                          ! lmlm-dos
! Integrate only up to muffin-tin radius.                                                 ! lmlm-dos
!  gflle_part = czero                                                                 ! lmlm-dos
!  DO lm2 = 1,lmmaxso                                                                 ! lmlm-dos
!    DO lm1 = 1,lmmaxso                                                               ! lmlm-dos
! For integration up to MT radius do this:                                                ! lmlm-dos
!              CWR(1:IMT1) = WR(LM1,LM2,1:IMT1)                                             ! lmlm-dos
!              CWR(IMT1+1:IRMDNEW) = CZERO                                                  ! lmlm-dos
!              CALL INTCHEB_CELL(CWR,GFLLE_PART(LM1,LM2),RPAN_INTERVALL,                    ! lmlm-dos
!      +                            IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)                  ! lmlm-dos
! For full cell integration replace loop content with this:                               ! lmlm-dos
!      cwr(1:irmdnew) = wr(lm1,lm2,1:irmdnew)                                        ! lmlm-dos
!      DO ir=imt1+1,irmdnew                                                       ! lmlm-dos
!        cwr(ir) = cwr(ir)*thetasnew(ir,1)*c0ll                                  ! lmlm-dos
!      END DO                                                                      ! lmlm-dos
!      CALL intcheb_cell(cwr,gflle_part(lm1,lm2),rpan_intervall, &                    ! lmlm-dos  &
!          ipan_intervall,npan_tot,ncheb,irmdnew)                  ! lmlm-dos
!    END DO                                                                           ! lmlm-dos
!  END DO                                                                              ! lmlm-dos
!END IF  ! OPT('lmlm-dos')


!      DO IR = 1,IRMDNEW
!       DO JSPIN = 1,4
!        DO LM1 = 1,LMMAXD
!         DO LM2 = 1,LM1-1
!          WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)=
!    +           WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)+
!    +           WR(LM2+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
!         ENDDO
!        ENDDO
!       ENDDO ! JSPIN
!      ENDDO !IR


! first calculate the spherical symmetric contribution

DO l1 = 0,lmax
  
  DO m1 = -l1,l1
    lm1 = l1*(l1+1)+m1+1
    DO ir = 1,irmdnew
      DO jspin=1,4
        cden(ir,l1,jspin) = cden(ir,l1,jspin)+  &
            wr(lm1+lmshift1(jspin),lm1+lmshift2(jspin),ir)
        cdenlm(ir,lm1,jspin) = wr(lm1+lmshift1(jspin),lm1+lmshift2(jspin),ir)
      END DO ! JPSIN
    END DO ! IR
  END DO ! M1
  
  DO jspin = 1,4
    DO ir = 1,irmdnew
      rho2nsc(ir,1,jspin) = rho2nsc(ir,1,jspin)+ c0ll*(cden(ir,l1,jspin)*df)
    END DO ! IR
    
    DO ir=imt1+1,irmdnew
      cden(ir,l1,jspin) = cden(ir,l1,jspin)*thetasnew(ir,1)*c0ll
      
      DO m1 = -l1,l1
        lm1 = l1*(l1+1)+m1+1
        cdenlm(ir,lm1,jspin) = cdenlm(ir,lm1,jspin) *thetasnew(ir,1)*c0ll
      END DO ! M1
    END DO ! IR
    
  END DO ! JSPIN
  
END DO ! L1

cdenns=czero

DO j = 1,iend
  lm1 = icleb(j,1)
  lm2 = icleb(j,2)
  lm3 = icleb(j,3)
  cltdf = df*cleb(j)
  
  DO jspin = 1,4
    DO ir = 1,irmdnew
      rho2nsc(ir,lm3,jspin) = rho2nsc(ir,lm3,jspin) +  &
          (cltdf*wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir))
    END DO
    
    IF (lmsp(lm3) > 0) THEN
      ifun = ifunm(lm3)
      DO ir=imt1+1,irmdnew
        cdenns(ir,jspin) = cdenns(ir,jspin)+  &
            cleb(j)*wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir)*  &
            thetasnew(ir,ifun)
      END DO
    END IF
  END DO ! JSPIN
END DO ! J


deallocate(wr)
deallocate(cwr)
deallocate(qnsi)
deallocate(pnsi)
END SUBROUTINE rhooutnew


subroutine intcheb_cell(cden,den,rpan_intervall,ipan_intervall, &
                        npan_tot,ncheb,irmdnew)
!***********************************************************************
! integrate the complex density of states for LM=1 
! gives the total complex charge which is then
! transformed to the xyz component of the magnetic 
! moment
!***********************************************************************
implicit none

integer           :: ncheb,npan_tot,irmdnew
integer           :: ipan_intervall(0:npan_tot)
double precision  :: rpan_intervall(0:npan_tot)
double complex    :: cden(irmdnew),den
integer           :: irstart,irstop,ipan
double precision  :: widthfac
double complex    :: int1

den=(0.0D0,0.0D0)

  do ipan=1,npan_tot
    irstart=ipan_intervall(ipan-1)+1
    irstop = ipan_intervall(ipan)
    widthfac = 0.5D0*(rpan_intervall(ipan)-rpan_intervall(ipan-1))
    call intcheb_complex(ncheb,cden(irstart:irstop),int1)
    den=den+int1*widthfac
    end do

end subroutine intcheb_cell 

subroutine intcheb_complex(ncheb,arr1,result1)
implicit none
integer, intent(in)         :: ncheb
double complex, intent(in)  :: arr1(0:ncheb)
double complex, intent(out) :: result1
double precision            :: pi
double precision  :: intweight(0:ncheb)
integer :: icheb1,icheb2

pi=4d0*datan(1d0)
intweight=1.0D0
  do icheb1=0,ncheb
    do icheb2=2,ncheb,2
      intweight(icheb1)=intweight(icheb1)+(-2.0D0/(icheb2**2-1.0D0))*dcos(icheb2*pi*(icheb1+0.5D0)/(Ncheb+1))
    end do
    intweight(icheb1)=intweight(icheb1)*2.0D0/(Ncheb+1)
  end do

result1=(0.0D0,0.0D0)
do icheb1=0,ncheb
  result1=result1+intweight(icheb1)*arr1(icheb1)
end do

end subroutine

subroutine cheb2oldgrid(nrmax,nrmaxnew,lmmaxpot,rmesh,ncheb, &
                        npan_tot,rpan_intervall,ipan_intervall, &
                        arrayin,arrayout,irmd)
implicit none
! Interpolate from Chebyshev mesh to the old radial mesh
! Programmed by David Bauer, 2011-2013
integer  :: ncheb,npan_tot,nrmax,nrmaxnew,lmmaxpot,irmd
double precision :: rmesh(irmd)
double precision :: rpan_intervall(0:npan_tot)
integer :: ipan_intervall(0:npan_tot)
double complex :: arrayin(nrmaxnew,lmmaxpot)
double complex :: arrayout(nrmax,lmmaxpot)
!local
integer :: in,ir,ir2,ir3,it,ilm
integer :: intsub(2,nrmax)
double precision :: Cinvmatrix(0:ncheb,0:ncheb)
double precision,allocatable :: CCmatrix(:,:)
double complex :: alphaparams(0:ncheb,lmmaxpot)
double precision  :: rmeshnorm(nrmax)
double precision halfsum,halfdiffinv,tol
parameter(tol=1.d-13)

! divide the mesh into subintervals
intsub(1,:)=0
intsub(2,:)=-1
in=1

ir=1
it=1
do while (ir.LE.nrmax .and. in.LE.npan_tot)
  if (abs(rmesh(ir)-rpan_intervall(in)).LT.tol) then
    intsub(it,in)=ir
    intsub(2,in)=ir
    it=1
    ir=ir+1
    in=in+1
  elseif (((rmesh(ir)-rpan_intervall(in).LT.-1d-10) ) .and. &
          ((rmesh(ir)-rpan_intervall(in-1)).GT.-1d-10) ) then
  intsub(it,in)=ir
  intsub(2,in)=ir
  it=2
  ir=ir+1
  elseif ((rmesh(ir)-rpan_intervall(in-1)).LT.-1d-10)  then
    intsub(1,in)=0
    intsub(2,in)=-1
    it=1
    ir=ir+1
  else
    it=1
    in=in+1
  end if
end do
if (abs(rmesh(nrmax)- rpan_intervall(npan_tot))>tol) then
   write(*,*) rmesh(nrmax),rpan_intervall(npan_tot),rmesh(nrmax)-rpan_intervall(npan_tot)
  stop 'error with the new and old mesh'
else
  
end if


call getCinvmatrix(Ncheb,Cinvmatrix)

in=0
do while (in<=npan_tot)
  in=in+1
  if (intsub(2,in)<intsub(1,in)) cycle
  alphaparams=0.0D0
  do ilm=1,lmmaxpot
  do ir2=0,ncheb
    do ir=0,ncheb
      alphaparams(ir2,ilm)=alphaparams(ir2,ilm)+ Cinvmatrix(ir2,ir)* arrayin(ipan_intervall(in-1)+1+ir,ilm)
    end do
  end do
  end do

  ! Transform to normalized coordinates between [-1,1]
  ! Shift upper & lower end to be centered around zero
  halfsum = 0.5d0 * (rpan_intervall(in) + rpan_intervall(in-1)) ! 0.5 * (upper + lower)
  halfdiffinv = 1.d0/(rpan_intervall(in) - halfsum)               ! 1. / (0.5*(upper - lower))  ! change Long 1

  do ir=intsub(1,in),intsub(2,in)
    ir2=ir+1-intsub(1,in)
    rmeshnorm(ir2)=(2*rmesh(ir)-(rpan_intervall(in)+rpan_intervall(in-1))) &
                                    /(rpan_intervall(in)-rpan_intervall(in-1))
    if (abs(rmeshnorm(ir2))>1.0D0 .and. abs(rmeshnorm(ir2))-1.0D0<10e-9) then
      rmeshnorm(ir2)=sign(1.0D0,rmeshnorm(ir2))
    end if
   enddo

  allocate(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

  call getCCmatrix(Ncheb,rmeshnorm(1:intsub(2,in)-intsub(1,in)+1),intsub(2,in)-intsub(1,in)+1,CCmatrix)
  do ilm=1,lmmaxpot
  do ir2=intsub(1,in),intsub(2,in)
    do ir=0,ncheb
      ir3=ir2-intsub(1,in)+1
      arrayout(ir2,ilm)=arrayout(ir2,ilm)+ CCmatrix(ir3,ir)* alphaparams(ir,ilm)
    end do
  end do
  end do !ilm=1,lmmaxpot
  deallocate(CCmatrix) !(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

end do !in

end subroutine cheb2oldgrid

subroutine rotatevector(rho2nsc,rho2ns,nrmax,lmpotd,theta,phi,theta_old,phi_old,nrmaxd)
implicit none
!***********************************************************************
!    does the rotation from the old local to the new local spin frame reference
!    for densities and charges 
!     rho_loc(ir,lm)= W1 * rho_glob(ir,lm) * W2 
!     where rho and W are matricies in spin space
!***********************************************************************
!interface
double complex   :: rho2nsc(nrmaxd,lmpotd,4)
double precision :: rho2ns(nrmaxd,lmpotd,4)
integer          :: nrmaxd,lmpotd,nrmax
double precision :: theta,phi
double precision :: theta_old,phi_old
!local
!double precision :: dcostheta2,dsintheta2
!double complex   :: im,cimphi,imphi
integer          :: ir,ilm
double complex   :: W1(2,2),W2(2,2)
double complex   :: W1_11W2_11, W1_11W2_22, W1_11W2_12, W1_11W2_21
double complex   :: W1_22W2_11, W1_22W2_22, W1_22W2_12, W1_22W2_21
double complex   :: W1_12W2_11, W1_12W2_22, W1_12W2_12, W1_12W2_21
double complex   :: W1_21W2_11, W1_21W2_22, W1_21W2_12, W1_21W2_21

call create_Wmatrix(theta,phi,theta_old,phi_old,1,W1,W2)

W1_11W2_11=W1(1,1)*W2(1,1)
W1_11W2_22=W1(1,1)*W2(2,2)
W1_11W2_12=W1(1,1)*W2(1,2)
W1_11W2_21=W1(1,1)*W2(2,1)
W1_22W2_11=W1(2,2)*W2(1,1)
W1_22W2_22=W1(2,2)*W2(2,2)
W1_22W2_12=W1(2,2)*W2(1,2)
W1_22W2_21=W1(2,2)*W2(2,1)
W1_12W2_11=W1(1,2)*W2(1,1)
W1_12W2_22=W1(1,2)*W2(2,2)
W1_12W2_12=W1(1,2)*W2(1,2)
W1_12W2_21=W1(1,2)*W2(2,1)
W1_21W2_11=W1(2,1)*W2(1,1)
W1_21W2_22=W1(2,1)*W2(2,2)
W1_21W2_12=W1(2,1)*W2(1,2)
W1_21W2_21=W1(2,1)*W2(2,1)



do ir=1,nrmax
  do ilm=1,lmpotd
    rho2ns(ir,ilm,1)= dimag(&
    +(rho2nsc(ir,ilm,1)*W1_11W2_11) &
    +(rho2nsc(ir,ilm,3)*W1_11W2_21) &
    +(rho2nsc(ir,ilm,4)*W1_12W2_11) &
    +(rho2nsc(ir,ilm,2)*W1_12W2_21)) 

    rho2ns(ir,ilm,2)= dimag (&
    +(rho2nsc(ir,ilm,1)*W1_21W2_12)  &
    -(rho2nsc(ir,ilm,3)*W1_21W2_22) &
    -(rho2nsc(ir,ilm,4)*W1_22W2_12) &
    +(rho2nsc(ir,ilm,2)*W1_22W2_22))
  end do
enddo 

end subroutine


subroutine calc_orbitalmoment(lmax,lmsize,Loperator)
implicit none
integer                         :: lmax,lmsize
double complex                  :: Loperator(lmsize,lmsize,3)
integer,save                    :: first=1
integer                         :: lval
double complex,allocatable      :: lorbit_onel(:,:,:)
integer                         :: lmmax,lstart,lstop
lmmax=(lmax+1)**2
Loperator=(0.0D0,0.0D0)

Loperator(1,1,1)=(0.0D0,0.0D0)
Loperator(1,1,2)=(0.0D0,0.0D0)
Loperator(1,1,3)=(0.0D0,0.0D0)

do lval=1,lmax
  allocate(lorbit_onel(  2*lval+1,2*lval+1,3 )  )  
  lorbit_onel=(0.0D0,0.0D0)
  call calc_orbit_onel(lval,Lorbit_onel)
  lstart=((lval-1)+1)**2+1
  lstop=((lval  )+1)**2
  Loperator(lstart:lstop,lstart:lstop,1) = Lorbit_onel(:,:,1)
  Loperator(lstart:lstop,lstart:lstop,2) = Lorbit_onel(:,:,2)
  Loperator(lstart:lstop,lstart:lstop,3) = Lorbit_onel(:,:,3)
  deallocate(lorbit_onel)
end do
if (lmsize/=lmmax) then
  Loperator(lmmax+1:lmsize,lmmax+1:lmsize,:) = Loperator(:lmmax,:lmmax,:)
end if

if (first==1) then
!  open(unit=423492157,file='out_Lx')
!  open(unit=423492158,file='out_Ly')
!  open(unit=423492159,file='out_Lz')
!  do ilm=1,lmsize
!    write(423492157,'(5000F)'),Loperator(ilm,:,1)
!    write(423492158,'(5000F)'),Loperator(ilm,:,2)
!    write(423492159,'(5000F)'),Loperator(ilm,:,3)
!  end do
!  close(423492157);close(423492158);close(423492159)
end if

first=0
end subroutine calc_orbitalmoment





subroutine calc_orbit_onel(lval,Lorbit_onel)
implicit none
integer        :: lval
double complex Lorbit_onel(2*lval+1,2*lval+1,3)
double complex L_z (-lval:lval,-lval:lval),L_x (-lval:lval,-lval:lval),L_y(-lval:lval,-lval:lval)
double complex L_dn(-lval:lval,-lval:lval),L_up(-lval:lval,-lval:lval)
integer :: i1
double precision :: lfac
double complex,parameter :: icompl=(0.0D0,1.0D0)


L_z=(0.0D0,0.0D0)
L_x=(0.0D0,0.0D0)
L_y=(0.0D0,0.0D0)
L_dn=(0.0D0,0.0D0)
L_up=(0.0D0,0.0D0)

Lorbit_onel=(0.0D0,0.0D0)
       do i1=-lval,lval
         L_z(-i1,i1)=-icompl*i1 
       end do 





       IF (lval>0) then

         lfac=sqrt(lval*(lval+1d0))/sqrt(2d0)
         l_dn(0,-1)=-icompl*lfac
         l_dn(0,1)=lfac
         l_dn(-1,0)=icompl*lfac
         l_dn(1,0)=-lfac
       
         IF (lval > 1) then

            do i1=2,lval

              lfac=0.5d0*SQRT(lval*(lval+1d0)-i1*(i1-1d0))
              l_dn(-i1,-i1+1)=-lfac
              l_dn(-i1,i1-1)=icompl*lfac
              l_dn(i1,-i1+1)=-icompl*lfac
              l_dn(i1,i1-1)=-lfac

              lfac=0.5d0*SQRT(lval*(lval+1d0)-(i1-1d0)*i1)
              l_dn(-i1+1,-i1)=lfac
              l_dn(-i1+1,i1)=icompl*lfac
              l_dn(i1-1,-i1)=-icompl*lfac
              l_dn(i1-1,i1)=lfac

            end do

         END IF
       END IF




       IF (lval>0) then

         lfac=sqrt(lval*(lval+1d0))/sqrt(2d0)
         l_up(0,-1)=-icompl*lfac
         l_up(0,1)=-lfac
         l_up(-1,0)=icompl*lfac
         l_up(1,0)=lfac
       
         IF (lval > 1) then

            do i1=2,lval

              lfac=0.5d0*SQRT(lval*(lval+1d0)-i1*(i1-1d0))
              l_up(-i1,-i1+1)=lfac
              l_up(-i1,i1-1)=icompl*lfac
              l_up(i1,-i1+1)=-icompl*lfac
              l_up(i1,i1-1)=lfac

              lfac=0.5d0*SQRT(lval*(lval+1d0)-(i1-1d0)*i1)
              l_up(-i1+1,-i1)=-lfac
              l_up(-i1+1,i1)=icompl*lfac
              l_up(i1-1,-i1)=-icompl*lfac
              l_up(i1-1,i1)=-lfac

            end do

         END IF
       END IF


L_x =  0.5D0          * (L_up+L_dn)
L_y = -0.5D0 * icompl * (L_up-L_dn)


Lorbit_onel(:,:,1)=L_x
Lorbit_onel(:,:,2)=L_y
Lorbit_onel(:,:,3)=L_z

end subroutine calc_orbit_onel


subroutine getCCmatrix(Ncheb,rmesh,nrmesh,Cmatrix)
! calculates the C matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
implicit none
integer  :: ncheb
double precision :: rmesh(nrmesh)
double precision :: Cmatrix(1:nrmesh,0:Ncheb)
integer  :: icheb,nrmesh,ir

do ir=1,nrmesh
  do icheb=0,ncheb
    Cmatrix(ir,icheb)=cos(icheb*acos(rmesh(ir)))
  end do
end do
end subroutine getCCmatrix


subroutine create_Wmatrix(theta,phi,theta_old,phi_old,lmmax,Wmat1,Wmat2)
implicit none
!***********************************************************************
! create the rotation matrix W:
!
!  W1= U_degga_locnew * U_locold
!  W2= U_degga_locold * U_locnew
!
!  the rotation matrix is created such that it rotates an operator
!  which is in a local frame (locold) to another local frame (locnew)
!  This is done by first transforming the old local frame to the
!  global frame using the U matrix and then transforming the global
!  frame to the new local frame
!
!  A_locnew = W1 * A_locold * W2

!  Udegga = transpose(complex conjug ( U ) )
!***********************************************************************double precision :: phi
!interface
double precision,intent(in)     :: phi
double precision,intent(in)     :: theta
double precision,intent(in)     :: phi_old
double precision,intent(in)     :: theta_old
integer,intent(in)              :: lmmax
double complex,intent(out)      :: Wmat1(2*lmmax,2*lmmax)
double complex,intent(out)      :: Wmat2(2*lmmax,2*lmmax)
!local
double complex                  :: Umat1     (2*lmmax,2*lmmax)
double complex                  :: Udeggamat1(2*lmmax,2*lmmax)
double complex                  :: Umat2     (2*lmmax,2*lmmax)
double complex                  :: Udeggamat2(2*lmmax,2*lmmax)
!double precision                :: matmat_zmzm

call create_Umatrix(theta_old,phi_old,lmmax,Umat1,Udeggamat1)

call create_Umatrix(theta,phi,lmmax,Umat2,Udeggamat2)

!Wmat1 = matmat_zmzm(Udeggamat2,Umat1)
!Wmat2 = matmat_zmzm(Udeggamat1,Umat2)
  call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Udeggamat2,2*lmmax,Umat1,2*lmmax,(0d0,0d0),Wmat1,2*lmmax)
  call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Udeggamat1,2*lmmax,Umat2,2*lmmax,(0d0,0d0),Wmat2,2*lmmax)


end subroutine create_Wmatrix


endmodule !NonCollinearMagnetism_mod
