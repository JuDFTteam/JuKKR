
module JelliumPotentials_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: jellstart12
  
  contains
  

!---------- Routines for creation of Jellium potentials -----------

SUBROUTINE jellstart12(nspin,ins,natoms,z,idshape,  &
        rwscl,rmtcl,meshn,xrn,drn,  &
        irws,irns,  &
        alatnew,qbound,dims,atom_index)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-21  Time: 16:35:08
 
! ******************************************************
! * This subroutine reads a jellium potential from the database
! * file. and interpolates to the new mesh
! ******************************************************

use DimParams_mod, only: DimParams

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nspin
INTEGER, INTENT(IN)                      :: ins
INTEGER, INTENT(IN)                      :: natoms
REAL*8, INTENT(IN)                       :: z(:)
INTEGER, INTENT(IN)                      :: idshape(*)
REAL*8, INTENT(IN)                       :: rwscl(*)
REAL*8, INTENT(IN)                       :: rmtcl(*)
INTEGER, INTENT(IN)                      :: meshn(:)
REAL*8, INTENT(IN)                       :: xrn(:,:)
REAL*8, INTENT(IN)                       :: drn(:,:)
INTEGER, INTENT(IN)                      :: irws(:)
INTEGER, INTENT(IN)                      :: irns(:)
REAL*8, INTENT(IN)                       :: alatnew
REAL*8, INTENT(IN OUT)                   :: qbound
type(DimParams), intent(in)              :: dims
INTEGER, INTENT(IN)                      :: atom_index


! Parameters that are no longer taken from 'inc.geometry' but depend on dims
INTEGER :: npotd,lmpotd,irmind,inslpd,lmxspd,irmdjj

!     ..
!     .. Scalar Arguments ..
REAL*8           efermi
INTEGER :: kshape
!     ..
!     .. Array Arguments ..

INTEGER :: lcore(20),ncore
!     ..
!     .. Local Scalars ..
REAL*8           ea,s1,z1,  &
    vbc(2),ecore1(20),maxa,aout,bout,rmtout,  &
    parsum,parsumderiv,r0,rmaxout,rmtnew
REAL*8           rws0,br,za,zvali,einf,ar,amsh
INTEGER :: i,ir,iri,  &
    ispin,  &
    ipot,id,lm,lm1,irnsout,irmtout,irwsout,  &
    nr,iat,ncore1,lcore1(20),irc,nz,  &
    irs1,nsec,nzvali,nc,ii,i1,i2
LOGICAL, allocatable :: potlm(:)
!     ..
!     .. Local Arrays ..
REAL*8, allocatable  ::    u(:),drdi(:),ecore(:),  &
                           rmesh(:),vins(:,:),  &
                           vm2z(:),vinsout(:,:),  &
                           vm2zout(:),vm2zb(:),rout(:), vinsb(:,:),drdiout(:),  &
                           work(:,:),ra(:)
CHARACTER (LEN=40) :: baner
CHARACTER (LEN=4) :: aaaa,tran

CHARACTER (LEN=4) :: elem_file(0:113)
CHARACTER (LEN=26) :: atompot
CHARACTER (LEN=2) :: txtc(20)
character(len=17) :: filename

DATA elem_file/'Vac0',  &
    'H_01','He02','Li03','Be04','B_05','C_06','N_07','O_08', 'F_09','Ne10',  &
    'Na11','Mg12','Al13','Si14','P_15','S_16','Cl17','Ar18',  &
    'K_19','Ca20','Sc21','Ti22', 'V_23','Cr24','Mn25','Fe26','Co27','Ni28',  &
    'Cu29','Zn30', 'Ga31','Ge32','As33','Se34','Br35','Kr36','Rb37','Sr38',  &
    'Y_39','Zr40', 'Nb41','Mo42','Tc43','Ru44','Rh45','Pd46','Ag47','Cd48',  &
    'In49','Sn50', 'Sb51','Te52','I_53','Xe54','Cs55','Ba56','La57','Ce58',  &
    'Pr59','Nd60', 'Pm61','Sm62','Eu63','Gd64','Tb65','Dy66','Ho67','Er68',  &
    'Tm69','Yb70', 'Lu71','Hf72','Ta73','W_74','Re75','Os76','Ir77','Pt78',  &
    'Au79','Hg80', 'Tl81','Pb82','Bi83','Po84','At85','Rn68','Fr87','Ra88',  &
    'Ac89','Th90', 'Pa91','U_92','Np93','Pu94','Am95','Cm96','Bk97','Cf98',  &
    'Es99','Fm__', 'Md__','No__','Lr__','Rf__','Db__','Sg__','Bh__','Hs__',  &
    'Mt__','Uun_', 'Uuu_','Uub_','NoE_'/
!     --------------------------------------------------------------

! set parameters depending on 'dims' and 'num_local_atoms'
npotd=dims%nspind*natoms
lmpotd= (dims%lpot+1)**2
irmind=dims%irmd-dims%irnsd
inslpd= (dims%irnsd+1)*lmpotd
lmxspd= (2*dims%lpot+1)**2
irmdjj=1501
kshape=2  ! always full-pot calculations
WRITE(*,*) 'dims%irmd', dims%irmd

! allocate arrays
allocate(u(dims%irmd))
allocate(drdi(irmdjj))
allocate(ecore(20))
allocate(rmesh(irmdjj))
allocate(vins(irmind:dims%irmd,lmpotd))
allocate(vm2z(irmdjj))
allocate(vinsout(irmind:dims%irmd,lmpotd))
allocate(vm2zout(dims%irmd))
allocate(vm2zb(irmdjj))
allocate(rout(dims%irmd))
allocate(vinsb(dims%irmd,lmpotd))
allocate(drdiout(dims%irmd))
allocate(work(dims%irmd,lmpotd))
allocate(ra(dims%irmd))
allocate(potlm(lmpotd))


WRITE(6,*) ' ****  READING  POTENTIAL  **** '

!OPEN(19,STATUS='UNKNOWN',FILE='output.pot')
write(filename, fmt="(a,i7.7)") "potential.",atom_index
open(19, file=filename, form="formatted", action='write')

DO i2=1,lmpotd
  DO i1=irmind,dims%irmd
    vins(i1,i2) = 0.d0
  END DO
END DO
DO i1=1,irmdjj
  vm2z(i1) = 0.d0
END DO

DO iat = 1,natoms
  DO ispin=1,nspin
    DO lm=1,lmpotd
      potlm(lm) =.false.
    END DO
    ipot =  nspin* (iat-1) + ispin
    
! Find out what atom is needed
    
    nz = z(iat)
    IF (((nz >= 24.AND.nz <= 28).OR.(nz >= 57.AND.nz <= 70))  &
          .AND.ispin == 2) THEN
      atompot = 'ElementDataBase/'//elem_file(nz)//'.pots2'
    ELSE
      atompot = 'ElementDataBase/'//elem_file(nz)//'.pot  '
    END IF
    WRITE(6,*) 'Using database ....: ',atompot
    OPEN(21,STATUS='OLD',FILE=atompot,ERR=1010)
!           IRWS1 =  NR
    
! --------------------------------------------------------------------
    
    efermi = .409241D+00
    vbc(1)    = .500D0
    vbc(2)    = .500D0
!           read potential from jellium
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    READ(21,141) baner,aaaa,aaaa
    READ(21,142) rws0,s1,irs1,br
    READ(21,142) za,zvali,nsec,einf
!     Calculate number of core states
    nz = za             ! make integer
    nzvali = zvali      ! make integer
    nc = nz - nzvali
    ncore = 0
    IF (nc == 2 ) ncore = 1 ! 1s
    IF (nc == 4 ) ncore = 2 ! 1s2s
    IF (nc == 10) ncore = 3 ! 1s2s2p
    IF (nc == 12) ncore = 4 ! 1s2s2p3s
    IF (nc == 18) ncore = 5 ! 1s2s2p3s3p
    IF (nc == 28) ncore = 6 ! 1s2s2p3s3p3d
    IF (nc == 30) ncore = 7 ! 1s2s2p3s3p3d4s
    IF (nc == 36) ncore = 8 ! 1s2s2p3s3p3d4s4p
    IF (nc == 46) ncore = 9 ! 1s2s2p3s3p3d4s4p4d
    IF (nc == 48) ncore = 10 ! 1s2s2p3s3p3d4s4p4d4s
    IF (nc == 54) ncore = 11 ! 1s2s2p3s3p3d4s4p4d4s4p
    IF (nc == 68) ncore = 12 ! 1s2s2p3s3p3d4s4p4d4s4p4f
    IF (nc == 78) ncore = 13 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d
    IF (nc == 80) ncore = 14 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s
    IF (nc == 86) ncore = 15 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s4p
    WRITE(6,*) '*************************************'
    WRITE(6,*) '   Potential Interpolation Program   '
    WRITE(6,*) '   Using the Jellium Database v1.0   '
    WRITE(6,*) '*************************************'
    WRITE(6,163) efermi
    WRITE(6,161) za
    WRITE (6,162) ncore
    READ(21,133)(lcore(nc),txtc(nc),ecore(nc),nc=1,ncore)
    WRITE(6,*) ' ** Position of the Core States ** '
    DO i=1,ncore
      WRITE(6,135) lcore(i),txtc(i),ecore(i)
      IF (txtc(i) == 's ') lcore(i) = 0
      IF (txtc(i) == 'p ') lcore(i) = 1
      IF (txtc(i) == 'd ') lcore(i) = 2
      IF (txtc(i) == 'f ') lcore(i) = 3
    END DO
    ncore1 = ncore
    DO i=1,ncore1
      lcore1(i) = lcore(i)
      ecore1(i) = ecore(i)
    END DO
    WRITE(6,134) einf
    WRITE(6,*) '**********************************************'
    
    READ(21,131)(vm2z(ii),ii=1,irs1)
    READ(21,132)tran
    CLOSE (21)
!     make mesh r0
    ar = LOG(s1/br+1.d0)/FLOAT(irs1-1)
    ea=EXP(ar)
    amsh=1.d0
    rmesh(1)=0.d0
    drdi(1)=br*ar
    DO  i=2,irs1
      amsh=amsh*ea
      rmesh(i)=br*(amsh-1.d0)
      drdi(i)=drdi(1)*amsh
    END DO
    WRITE(6,*) 'Jellium Potential Read In ',irs1,' points'
    nr = irs1
    z1 = za
    
    131        FORMAT(4D15.8)
    132        FORMAT(a4)
    133        FORMAT(i3,a2,d15.8)
    135        FORMAT(i3,a2,f15.6,' Ry')
    134        FORMAT('All other states are above :',f8.4,' Ry in Energy')
    141        FORMAT(3X,a40,3X,a4,3X,a4)
    142        FORMAT(7X,f8.4,7X,f8.4,7X,i5,7X,f8.4)
    161        FORMAT('Potential Atomic Number :',f7.2)
    162        FORMAT('Number of Core   States :',i4)
    163        FORMAT('Jellium Fermi Energy :',f10.5,' Ry')
! --------------------------------------------------------------------
    
!     The input mesh has been constructed. Now construct the output mesh.
    
    id = idshape(iat)
    rout(1) = 0.d0
    aout = 0.025d0
    rmaxout = rwscl(id)
    rmtout  = rmtcl(id)
    irwsout = irws(iat)
    irmtout = irws(iat) - meshn(id)
    irnsout = irns(iat)  ! 22.1.12 Changed from IRNS(ID) to IRNS(IAT)
    
    
    IF (ins == 0) THEN
      bout = rmaxout / (EXP(aout*REAL(irwsout-1))-1.0D0)
      DO ir=2,irwsout
        ea = EXP(aout*REAL(ir-1))
        rout(ir) = bout* (ea-1.0D0)
        drdiout(ir) = aout*bout*ea
      END DO
      DO i=1,irwsout
        IF (rout(i) < rmtout) irmtout = i
      END DO
      IF (MOD(irmtout,2) == 0) irmtout = irmtout+1
      rmtnew = rout(irmtout)
      rmaxout = rout(irwsout)
    ELSE
      bout = rmtout /  (EXP(aout*REAL(irmtout-1))-1.0D0)
      DO ir=2,irmtout
        ea = EXP(aout*REAL(ir-1))
        rout(ir) = bout* (ea-1.0D0)
        drdiout(ir) = aout*bout*ea
      END DO
      DO iri=1,meshn(id)
        ir = iri + irmtout
        rout(ir) = alatnew*xrn(iri,id)   ! scaling is always 1.0d0
        drdiout(ir) = alatnew*drn(iri,id)
      END DO
      rmtnew = rout(irmtout)
      rmaxout = rout(irwsout)
    END IF
    
!  Ok now interpolate
    
    
    maxa = 1.d35
    CALL spline(irmdjj,rmesh,vm2z,nr,maxa,maxa,vm2zb)
    
! OK with spline
    
    vm2zout(1) = vm2z(1)
    DO ir = 2,irwsout
      r0 = rout(ir)
      CALL splint(rmesh,vm2z,vm2zb,nr,r0,parsum,parsumderiv)
      vm2zout(ir) = parsum
    END DO
    
    
    
    
    IF (ins > 0) THEN
      irc = irwsout - irnsout
      DO lm1=1,lmpotd
        DO ir = irc,irwsout
          vinsout(ir,lm1) = 0.d0
        END DO
      END DO
    END IF
    
    CALL ritesone12(19,ispin,z1,alatnew,rmtout,rmtnew,rmaxout,  &
        rout,drdiout,vm2zout,irwsout,aout,bout,ins,irnsout,  &
        vinsout,qbound,irwsout,kshape,efermi,vbc,  &
        ecore1,lcore1,ncore1,elem_file(nz),nspin,dims)
    
! Next atom or next spin
    
  END DO
END DO

RETURN


!1000  WRITE(6,*) 'Error read file......... ',i13
WRITE(6,*) 'Error occured on atom... ',iat
STOP
1010  WRITE(6,*) ' Error in JELLSTART '
WRITE(6,*) ' Potential.............',elem_file(nz)
WRITE(6,*) ' Does not exist in the database'
STOP
8000 FORMAT (a40)
8010 FORMAT (3X,a4,26X,f8.3)
8011 FORMAT ('#  ',a4,'POTENTIAL             Z = ',f8.3)
8012 FORMAT ('#  ',a4,'POTENTIAL SPIN UP     Z=  ',f8.3)
8013 FORMAT ('#  ',a4,'POTENTIAL SPIN DOWN   Z=  ',f8.3)
8020  FORMAT(a40)
8030 FORMAT (4F12.8)
8040 FORMAT (1X,3I6)
8050 FORMAT (2D15.8)
8060 FORMAT (3F12.8)
8070 FORMAT (2I5)
9051 FORMAT (1P,4D20.12)
9000 FORMAT (7A4,6X,'  exc:',a24,3X,a10)
9010 FORMAT (3F12.8)
9020 FORMAT (f10.5,/,f10.5,2F15.10)
9030 FORMAT (i3,/,2D15.8,/,2I2)
9140 FORMAT (i5,1P,d20.11)
9040 FORMAT (f10.5,/,f10.5,2F15.10)
9050 FORMAT (i3,/,2D15.8,/,2I2)
9060 FORMAT (1P,2D15.6,1P,d15.8)
9160 FORMAT (10I5)
9061 FORMAT (1P,5D15.8)
9070 FORMAT (i5,1P,d20.11)
! 9080 FORMAT (10x,20a4)
9080 FORMAT (' < ',20A4)
9081 FORMAT (' <#',20A4)
9090 FORMAT (10I5)
9100 FORMAT (1P,4D20.13)
END SUBROUTINE jellstart12



SUBROUTINE ritesone12(ifile,is,z,alat,rmt,rmtnew,rws,  &
        r,drdi,vm2z,irws,a,b,ins,irns,  &
        vins,qbound,irc,kshape,efermi,vbc,ecore,  &
        lcore,ncore,elem_name,nspin,dims)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-21  Time: 16:59:42
 
! ************************************************************************
!      this subroutine stores in 'ifile' the necessary results
!      (potentials e.t.c.) to start self-consistency iterations

!      modified for the full potential case - if ins .gt. 0 there
!       is written a different potential card
!       if the sum of absolute values of an lm component of vins (non
!       spher. potential) is less than the given rms error qbound this
!       component will not be stored .

!        (see to subroutine start , where most of the arrays are
!         described)

!                            modified by b. drittler  aug. 1988
!-----------------------------------------------------------------------
!     .. Parameters ..


use DimParams_mod, only: DimParams

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ifile
INTEGER, INTENT(IN)                      :: is
REAL*8, INTENT(IN)                       :: z
REAL*8, INTENT(IN)                       :: alat
REAL*8, INTENT(IN)                       :: rmt
REAL*8, INTENT(IN)                       :: rmtnew
REAL*8, INTENT(IN)                       :: rws
REAL*8, INTENT(IN)                       :: r(:)
REAL*8, INTENT(IN)                       :: drdi(:)
REAL*8, INTENT(IN)                       :: vm2z(:)
INTEGER, INTENT(IN)                      :: irws
REAL*8, INTENT(IN)                       :: a
REAL*8, INTENT(IN)                       :: b
INTEGER, INTENT(IN)                      :: ins
INTEGER, INTENT(IN)                      :: irns
REAL*8, INTENT(IN)                       :: qbound
INTEGER, INTENT(IN)                      :: irc
INTEGER, INTENT(IN)                      :: kshape
REAL*8, INTENT(IN)                       :: efermi
REAL*8, INTENT(IN)                       :: vbc(2)
REAL*8, INTENT(IN)                       :: ecore(20)
INTEGER, INTENT(IN)                      :: lcore(20)
INTEGER, INTENT(IN)                      :: ncore
CHARACTER (LEN=4), INTENT(IN)            :: elem_name
INTEGER, INTENT(IN)                      :: nspin
type(DimParams), intent(in)              :: dims
REAL*8, INTENT(IN)                       :: vins((dims%irmd-dims%irnsd):dims%irmd,(dims%lpot+1)**2)

!     ..
!     .. Local Scalars ..
REAL*8           a1,b1,rmax,rmt1,rmtnw1,rv,sum,z1
INTEGER :: icore,inew,ir,irmin,irns1,isave,j,lm,lmnr,lmpot, &
ncore1,nr,lpot
!     ..
!     .. Local Arrays ..
REAL*8, allocatable :: dradi(:),ecore1(:),ra(:),vm2za(:)
INTEGER :: lcore1(20)
!     ..
!     .. Intrinsic Functions ..
INTRINSIC SQRT
!     ..
! -------------------------------------------------------------------

allocate(dradi(dims%irmd))
allocate(ecore1(20))
allocate(ra(dims%irmd))
allocate(vm2za(dims%irmd))

lpot=dims%lpot
isave = 1
inew  = 1



lmpot = (lpot+1)* (lpot+1)

rmt1 = rmt
rmtnw1 = rmtnew
z1 = z
rmax = rws
IF (kshape == 0) THEN
  nr = irws
  
ELSE
  nr = irc
END IF

irns1 = irns
irmin = nr - irns1
a1 = a
b1 = b
ncore1 = ncore

DO  j = 1,nr
  ra(j) = r(j)
  dradi(j) = drdi(j)
  
!--->       store only lm=1 component of the potential
  
  vm2za(j) = vm2z(j)
END DO

IF (ncore1 >= 1) THEN
  
  DO  j = 1,ncore1
    lcore1(j) = lcore(j)
    ecore1(j) = ecore(j)
  END DO
END IF


IF (nspin == 1) THEN
  WRITE (ifile,FMT=8999) elem_name,''
ELSE
  IF (is == 1) THEN
    WRITE (ifile,FMT=8995) elem_name,''
  ELSE
    WRITE (ifile,FMT=8996) elem_name,''
  END IF
END IF
!     WRITE (IFILE,FMT=9000) (ITITLE(I),I=1,7),TXC(KXC+1)
WRITE (ifile,FMT=9010) rmt1,alat,rmtnw1
WRITE (ifile,FMT=9020) z1,rmax,efermi,vbc(is)
IF (nr <= 999) THEN
  WRITE (ifile,FMT=9030) nr,a1,b1,ncore1,inew
ELSE
  WRITE (ifile,FMT=9031) nr,a1,b1,ncore1,inew
END IF
IF (ncore1 >= 1) WRITE (ifile,FMT=9040) (lcore1(icore),  &
    ecore1(icore),icore=1,ncore1)


IF (ins == 0 ) THEN
  
!--->       store only the spherically averaged potential
!           (in mt or as - case)
!           this is done always for the host
  
  IF (inew == 0) THEN
    WRITE (ifile,FMT=9050) (ra(ir),dradi(ir),vm2za(ir),ir=1,nr)
  ELSE
    WRITE (ifile,FMT=9051) (vm2za(ir),ir=1,nr)
  END IF
  
ELSE
  
!--->     store the full potential , but the non spherical contribution
!         only from irns1 up to irws1 ;
!         remember that the lm = 1 contribution is multiplied
!         by a factor 1/sqrt(4 pi)
  
  WRITE (ifile,FMT=9060) nr,irns1,lmpot,isave
  WRITE (ifile,FMT=9070) (vm2za(ir),ir=1,nr)
  IF (lpot > 0) THEN
    lmnr = 1
    DO  lm = 2,lmpot
      sum = 0.0D0
      DO  ir = irmin,nr
        rv = vins(ir,lm)*ra(ir)
        sum = sum + rv*rv*dradi(ir)
      END DO
      
      IF (SQRT(sum) > qbound) THEN
        lmnr = lmnr + 1
        WRITE (ifile,FMT=9060) lm
        WRITE (ifile,FMT=9070) (vins(ir,lm),ir=irmin,nr)
      END IF
      
    END DO
    
!--->         write a one to mark the end
    
    IF (lmnr < lmpot) WRITE (ifile,FMT=9060) isave
  END IF
  
END IF
!write(6,*) ' Potential finished go on'
50   CONTINUE
60 CONTINUE
8995 FORMAT (a4,' POTENTIAL SPIN DOWN',10X,'  exc:',a24)
8996 FORMAT (a4,' POTENTIAL SPIN UP  ',10X,'  exc:',a24)
8999 FORMAT (a4,' POTENTIAL ',19X,'  exc:',a24)
9000 FORMAT (7A4,6X,'  exc:',a24,3X,a10)
9010 FORMAT (3F12.8)
9020 FORMAT (f10.5,/,f10.5,2F15.10)
9030 FORMAT (i3,/,2D15.8,/,2I2)
9031 FORMAT (i4,/,2D15.8,/,2I2)
9040 FORMAT (i5,1P,d20.11)
9050 FORMAT (1P,2D15.6,1P,d15.8)
9051 FORMAT (1P,4D20.12)
9060 FORMAT (10I5)
9070 FORMAT (1P,4D20.13)
END SUBROUTINE ritesone12


!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-12  Time: 14:48:44

SUBROUTINE spline(nmax,x,y,n,yp1,ypn,y2)

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nmax
REAL*8, INTENT(IN)                       :: x(nmax)
REAL*8, INTENT(IN)                       :: y(nmax)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN)                       :: yp1
REAL*8, INTENT(IN)                       :: ypn
REAL*8, INTENT(OUT)                      :: y2(nmax)


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

IF (n > nmax) STOP 'SPLINE: n > NMAX.'
IF (yp1 > 0.99D30) THEN
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

IF (ypn > 0.99D30) THEN
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
END SUBROUTINE spline


!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-12  Time: 14:48:43

SUBROUTINE splint(xa,ya,y2a,n,x,y,yderiv)

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
IF (h == 0.d0) PAUSE 'bad xa input in splint'
! Cubic spline polynomial is now evaluated.
a = (xa(khi)-x)/h
b = (x-xa(klo))/h
y = a*ya(klo) + b*ya(khi) +  &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (h**2)/6.d0
yderiv = (ya(khi)-ya(klo))/h -  &
    ((3.d0*a*a-1.d0)*y2a(klo) - (3.d0*b*b-1.d0)*y2a(khi))*h/6.d0

RETURN
END SUBROUTINE splint


endmodule ! JelliumPotentials_mod
