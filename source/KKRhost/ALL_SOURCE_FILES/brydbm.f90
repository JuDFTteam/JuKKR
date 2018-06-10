!*********************************************************************
SUBROUTINE brydbm(visp,v,vins,vspsme,vspsmo,ins,lmpot,  &
    r,drdi,alpha,atwght,irc,irmin,nspin,  &
    natps,natyp,itdept,imix,iobroy,ipf,lsmear)
!*********************************************************************
!     imix :
!       3      broyden's            f i r s t  m e t h o d
!       4      broyden's          s e c o n d  m e t h o d
!       5      anderson's     g e n e r a l i z e d   m e t h o d

!     implemented here according to notes of s.b.
!     broyden's iteration scheme following the papers of :
!     srivastava, j. phys. , 17 (1984) , pp l317
!     c.g. broyden in math.comput., 19 , pp 577, 1965
!     c.g. broyden in ibid, 21 ,pp 368 ,1967
!     the method has been generalized to include a metric.the
!     definition of the necessary inner products are similar to the
!     discription given in the notes of m.weinert.the algorithm
!     discribed in the paper srivastava  has been simplified
!     ( see notes of s.b.)
!     the files ui,vi are stored on high speed ssd memory.
!     broyden's update treats charge and spin on the same footing
!                  s. bluegel , kfa , may 1987
!     the anderson method (d.g. anderson, j. acm 12, 547 (1964)) has
!     been generalized and reformulated as an improvement of broyden's
!     second method. successive linesearch is replaced by successive
!     search on hyperplanes. ( see notes of s.b. )
!                  s. bluegel , issp , july 1989

!     modified for non spherical potential
!                  b. drittler , aug. 1988
!*********************************************************************

use mod_types

!     .. Parameters ..
INCLUDE 'inc.p'

INTEGER LMPOTD
PARAMETER (LMPOTD= (LPOTD+1)**2)
INTEGER IRMIND
PARAMETER (IRMIND=IRMD-IRNSD)
INTEGER NSPINDD
PARAMETER (NSPINDD=2*KREL + (1-KREL)*NSPIND)
INTEGER NTIRD
PARAMETER (NTIRD= (IRMD*NTPERD+ (IRNSD+1)* (LMPOTD-1)*NATYPD)* &
          NSPINDD)
INTEGER ITDTHD
PARAMETER (ITDTHD=40)
!..
!.. Array Arguments ..
DOUBLE PRECISION ATWGHT(*),DRDI(IRMD,*),R(IRMD,*), &
                 V(IRMD,LMPOTD,*),VINS(IRMIND:IRMD,LMPOTD,*), &
                 VISP(IRMD,*),VSPSMO(IRMD,*),VSPSME(IRMD,*)
INTEGER IRC(*),IRMIN(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION CMM,ONE,RMIXIV,SMNORM,VMDENO,VMNORM,VOLINV,ZERO
INTEGER IA,IJ,IMAP,IR,IRC1,IRMIN1,ISP,IT,LM,MIT,LSMEAR
!..
!.. External Functions ..
DOUBLE PRECISION DDOT
EXTERNAL DDOT
!..
!.. External Subroutines ..
EXTERNAL BRYSH1,BRYSH2,BRYSH3,DAXPY,DSCAL,RCSTOP
!..
!.. Intrinsic Functions ..
INTRINSIC ABS
!..
!.. Save statement ..
SAVE MIT,ZERO,ONE,WIT
!..
!.. Local Arrays ..
DOUBLE PRECISION, allocatable :: AM(:),BM(:),FM(:), &
                 FM1(:),G(:),SM(:),SM1(:), &
                 VI3(:),WIT(:),UI2(:),UI3(:),VI2(:)
!..
!.. Scalar Arguments ..
DOUBLE PRECISION ALPHA
INTEGER IMIX,INS,IOBROY,IPF,ITDEPT,LMPOT,NATPS,NATYP,NSPIN
!..
!.. Data statements ..
DATA MIT/1/,ZERO,ONE/0.0D0,1.0D0/


allocate(am(2:itdthd-1),bm(2:itdthd-1),fm(ntird),  &
    fm1(ntird),g(ntird),sm(ntird),sm1(ntird),  &
    vi3(ntird),wit(2:200),ui2(ntird),ui3(ntird),vi2(ntird))

mit = t_inc%mit_bry

IF (itdept > itdthd .OR. itdthd > 200) CALL rcstop('itdbry  ')

IF (imix <= 2 .OR. imix > 5) CALL rcstop('IMIXD   ')

IF (mit > itdept) mit = 1
IF (imix == 3) WRITE (ipf,FMT='('' BROYDEN"S 1ST METHOD USED '')')
IF (imix == 4) WRITE (ipf,FMT='('' BROYDEN"S 2ND METHOD USED '')')
IF (imix == 5) WRITE (ipf,FMT= '('' GENERALIZED ANDERSON METHOD USED '')')
WRITE(ipf,'(A,i4)') ' Iteration index (read in):',mit

rmixiv = one/alpha

!---->  the following block is activated only one iteration before
!        broyden iteration scheme is used
!---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
!                   metric  g := r*r*drdi
!---->  map data of all muffin-tin spheres into one single vector


CALL brysh3(sm1,visp,vins,vspsme,ins,irmin,irc,natps,  &
    natyp,nspin,imap,lmpot,lsmear)
CALL brysh1(fm1,v,vspsmo,ins,irmin,irc,natps, natyp,nspin,imap,lmpot,lsmear)

IF (imap > ntird) CALL rcstop('NIRDBRY ')

DO  ij = 1,imap
  fm1(ij) = rmixiv* (fm1(ij)-sm1(ij))
END DO

ij = 0
DO  isp = 1,nspin
  DO  ia = natps,natyp
    irc1 = irc(ia)
    volinv = 3.0D0/ (r(irc1,ia)**3)
    DO  ir = 1,irc1
      ij = ij + 1
      g(ij) = atwght(ia)*volinv*r(ir,ia)*r(ir,ia)*drdi(ir,ia)
    END DO
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    
!     Next for SMEARED spherical potential
    
    IF ( lsmear > 0 ) THEN
      DO  ir = 1, irc1
        ij = ij + 1
        g(ij) = atwght(ia)*volinv*r(ir,ia)*r(ir,ia)*drdi(ir,ia)
      END DO
    END IF
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    IF (ins >= 1 .AND. lmpot > 1) THEN
      irmin1 = irmin(ia)
      DO  lm = 2,lmpot
        DO  ir = irmin1,irc1
          ij = ij + 1
          g(ij) = atwght(ia)*volinv*r(ir,ia)*r(ir,ia)*drdi(ir,ia)
        END DO
      END DO
    END IF
    
  END DO
  
END DO


IF (mit > 1) THEN
  REWIND iobroy + 2
  READ (iobroy+2) (sm1(ij),ij=1,imap), (fm1(ij),ij=1,imap)
  
  
!----> map rho(m) of all mt-spheres into one single vector
  
  CALL brysh3(sm,visp,vins,vspsme,ins,irmin,irc,natps,  &
      natyp,nspin,imap,lmpot,lsmear)
  
!----> map f[m] = f(m) - rho(m) = f(rho(m)) - rho(m) of all mt-spheres
!      into one single vector
  
  CALL brysh1(fm,v,vspsmo,ins,irmin,irc,natps, natyp,nspin,imap,lmpot,lsmear)
  DO  ij = 1,imap
    fm(ij) = rmixiv* (fm(ij)-sm(ij))
  END DO
  
!----> calculate  sm = rho(m) - rho(m-1)
!----> calculate dfm = f[m] - f[m-1]
  
  DO  ij = 1,imap
    sm1(ij) = sm(ij) - sm1(ij)
    fm1(ij) = fm(ij) - fm1(ij)
  END DO
  
!----> loop to generate u[m] = u(ij,mit)
  
  DO  ij = 1,imap
    ui3(ij) = alpha*fm1(ij) + sm1(ij)
  END DO
  REWIND iobroy
  DO  it = 2,mit - 1
    READ (iobroy) (ui2(ij),ij=1,imap), (vi2(ij),ij=1,imap),wit
    
    am(it) = ddot(imap,fm1,1,vi2,1)
    CALL daxpy(imap,-am(it),ui2,1,ui3,1)
  END DO
  
!----> print amj = the importance of the history of ui
  
  WRITE (ipf,FMT='(5x,'' AMJ , ---> J=2,'',I3,/,(9X,1P,7D10.2))')  &
      mit - 1, (am(it),it=2,mit-1)
  
  
  IF (imix == 3) THEN
!-------->     b r o y d e n ' s   f i r s t   m e t h o d
    
    
!----> calculate dsmnorm
    
    smnorm = zero
    DO  ij = 1,imap
      smnorm = smnorm + sm1(ij)*g(ij)*sm1(ij)
    END DO
    
!----> convolute dsm with the metric g
    
    DO  ij = 1,imap
      sm1(ij) = g(ij)*sm1(ij)
    END DO
    
!----> loop to generate v[m] = v(ij,mit)
    
    DO  ij = 1,imap
      vi3(ij) = alpha*sm1(ij)
    END DO
    REWIND iobroy
    DO  it = 2,mit - 1
      READ (iobroy) (ui2(ij),ij=1,imap), (vi2(ij),ij=1,imap),wit
      
      bm(it) = ddot(imap,sm1,1,ui2,1)
      CALL daxpy(imap,-bm(it),vi2,1,vi3,1)
    END DO
    
!----> complete the evaluation of v[m]
    
    vmdeno = ddot(imap,sm1,1,ui3,1) - smnorm
    
    IF (ABS(vmdeno) < 1D-70) CALL rcstop('BRY0SN  ')
    
    CALL dscal(imap,one/vmdeno,vi3,1)
    
!----> print bmj = the importance of the history of vi
    
    WRITE (ipf,FMT='(5x,'' BMJ , ---> J=2,'',I3,/,(9X,1P,7D10.2))'  &
        ) mit - 1, (bm(it),it=2,mit-1)
    
  ELSE IF (imix == 4) THEN
!-------->     b r o y d e n ' s   s e c o n d    m e t h o d
    
!----> calculate v[m] ; convoluted with the metric g
    
    DO  ij = 1,imap
      vi3(ij) = g(ij)*fm1(ij)
    END DO
    
!----> calculate #vm# and normalize v[m]
    
    vmnorm = ddot(imap,vi3,1,fm1,1)
    CALL dscal(imap,one/vmnorm,vi3,1)
    
  ELSE IF (imix == 5) THEN
!-------->     g e n e r a l i z e d   a n d e r s o n   m e t h o d
    
!----> calculate v[m] ; convoluted with the metric g
    
    DO  ij = 1,imap
      vi3(ij) = g(ij)*fm1(ij)
    END DO
    REWIND iobroy
    DO  it = 2,mit - 1
      READ (iobroy) (ui2(ij),ij=1,imap), (vi2(ij),ij=1,imap),wit
      
      CALL daxpy(imap,-am(it)*wit(it),vi2,1,vi3,1)
    END DO
    
!----> complete the evaluation of v[m]
    
    vmdeno = ddot(imap,fm1,1,vi3,1)
    
    IF (ABS(vmdeno) < 1D-70) CALL rcstop('BRY1SN  ')
    
    CALL dscal(imap,one/vmdeno,vi3,1)
    
!----> save wit(mit) for next iteration
    
    wit(mit) = vmdeno
    
  END IF
  
!----> write u3(ij) and v3(ij) on disk
  
  WRITE (iobroy) (ui3(ij),ij=1,imap), (vi3(ij),ij=1,imap),wit
  
!----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
  
  DO  ij = 1,imap
    fm1(ij) = fm(ij)
    sm1(ij) = sm(ij)
  END DO
  
!----> calculate cmm
  
  cmm = ddot(imap,fm,1,vi3,1)
!           WRITE (IPF,FMT='(5X,'' CMM = '',1P,D12.4)') CMM
  
!----> update rho(m+1)
  
  CALL daxpy(imap,one-cmm,ui3,1,sm,1)
  
!----> map solution back into each mt-sphere
  
  CALL brysh2(sm,v,vspsmo,ins,irmin,irc,natps, natyp,nspin,imap,lmpot,lsmear)
  
END IF
mit = mit + 1
t_inc%mit_bry = mit

REWIND iobroy + 2
WRITE (iobroy+2) (sm1(ij),ij=1,imap), (fm1(ij),ij=1,imap)

deallocate(am,bm,fm,fm1,g,sm,sm1,vi3,wit,ui2,ui3,vi2)

RETURN


END SUBROUTINE brydbm
