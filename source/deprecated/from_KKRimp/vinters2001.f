
      SUBROUTINE VINTERS2010(NATOM,NKILLATOM,CMOM,CMOMREF


,IREF,LMAX,NATPER,NATREF,
     $     NSPIN,NSTART,NEND,V,Z,R,IRWS,IRCUT,IPAN,KSHAPE,
     +     CMINST,irm,ntim,alat,ioper,ND,NSHELL,RM,W,YR,YRA,WTYR,
     &     RIJ,IJEND)
!       SUBROUTINE VINTERS2001(CMOM,IREF,LMAX,NATPER,NATREF,
!      $     NSPIN,NSTART,NEND,V,Z,R,IRWS,IRCUT,IPAN,KSHAPE,
!      +     CMINST,irm,ntim,alat,ioper,ND,NSHELL,RM,W,YR,YRA,WTYR,
!      &     RIJ,IJEND)
c-----------------------------------------------------------------------
c     calculate the intercell-potentials and add these to the poten-
c     tial v  (in the spin-polarized case for each spin-direction
c     the intercell-potential is the same . )
c     it uses the structure dependent matrices amat and bmat which
c     are calculate once in the subroutine amn .
c     the charge-moments are calculated in the subroutine vintra2 ,
c     therefore vintra2 has to be called first .
c     the intercell-potential is expanded into spherical harmonics .
c     the lm-term of the intercell-potential v of the representive
c     atom i is given by
c     
c     v(r,lm,i) =  (-r)**l * {amat(i1,i2,lm,l'm')*cmom(i2,l'm')
c     +bmat(i1,i2,lm)*z(i2)}
c     
c     summed over i2 (all shells) and l'm' .    (i1=i-natref)
c     (see notes by b.drittler)
c     
c     in case of shape correction the madelung potential of the host
c     is taken into account . in all other case the madelung poten-
c     tial of the host is set to be zero .
c     as actual values for z and cmom the differences between the
c     values of the  representive atoms and those of the references
c     are used .
c     
c     attention : the first index of cmom (moment of the charge
c     density - calculated in vintr2) and of z (nuclear
c     charge of the atoms in the shell) is in the program
c     different defined , there one has to use :
c     cmom(natref+i2,l'm')  and  z(natref+i2)
c     
c     b.drittler   june 1987
c--------------------------------------------------------------------
c In the case of impurities on surfaces, the intercell potential of the
c reference system is read in from the fxdr file 'intercell_ref'; which
c is calculated in the surface program. 
c                  july 1998
c-----------------------------------------------------------------------
C     .. Parameters ..
!       include 'gaunt.param'
      INTEGER IJD,ICJD,JCOEFF
      PARAMETER(IJD=434, ICJD=26670, JCOEFF=17261)
      include 'parameters.file'
c      INTEGER NATYPD,NTREFD,NTPERD
c      PARAMETER (NATYPD=19,NTREFD=4,NTPERD=NATYPD-NTREFD)
c      INTEGER LPOTD
c      PARAMETER (LPOTD=6)
c      INTEGER IPAND
c      PARAMETER (IPAND=5)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER N,LM3D
      PARAMETER (N=4*LMAXD,LM3D=(2*LPOTD+1)**2)
      INTEGER NCLEB
      PARAMETER (NCLEB=LM3D*LMPOTD)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,LMAX,NATPER,NATREF,NEND,NSPIN,NSTART,
     $     irm, ntim,IJEND
      REAL*8 alat 
C     ..
C     .. Array Arguments ..
      REAL*8 CMINST(LMPOTD,*),CMOM(LMPOTD,*),R(IRM,*),
     +     V(IRM,LMPOTD,*),Z(*)
      REAL*8 RIJ(IJD,3),RM(3,*),WTYR(IJD,LMPOTD),YRA(IJD,LMPOTD) ! www
      REAL*8 W(N),YR(N,0:N,0:N)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IREF(*),IRWS(*),IOPER(*),
     &        ND(48,3,3),NSHELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 AC,PI,SUM
      INTEGER I,I1,I2,IATYP,IATYP2,IPOT,IRS1,ISPIN,L,LM,LM2,LMMAX,M,
     $     NREF,ihandle,lmmax_file,natref_file,IEND
      logical first_call(2)
      
C     ..
C     .. Local Arrays ..
      REAL*8 ACH(LMPOTD,NTREFD,2),CLEB(NCLEB,2)
      INTEGER ICLEB(NCLEB,4),LOFLM(LM3D)
      REAL*8 AMAT(NTPERD,LMPOTD,LMPOTD),BMAT(NTPERD,LMPOTD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DSQRT
C     ..
C     .. Save statement ..
      SAVE ACH,first_call,CLEB,ICLEB,IEND
      
      DATA first_call /.true.,.true./
C     ..
      PI = 4.D0*DATAN(1.D0)
      LMMAX = (LMAX+1)* (LMAX+1)
c     read intercellpotential of reference system in first iteration
c     and calculate intercell contribution to reference potential 
c     (needed for total energy calculations)
!       if (first_call(ntim)) then
! c         call fxdropn('intercell_ref','DECODE',ihandle)
! 
!          CALL AMNGAUNT(LMAX,CLEB,ICLEB,IEND,W,YR)
! 
!          open(1,FILE='intercell_ref',STATUS='old',FORM='formatted')
!          read(1,*) natref_file,lmmax_file
! 
!          if (natref_file.ne.natref) then
!              STOP 'intercell_ref not correct natref'
!          end if
!          if (lmmax_file.gt.lmpotd) then
!              STOP 'intercell_ref not correct lmmax'
!          end if
!          do i1 =1,natref
! c---  > determine shell index             ! Corrected 1.2.2000
!          IF (KSHAPE.NE.0) THEN
!             IRS1 = IRCUT(IPAN(I1),I1)
!          ELSE
!             IRS1 = IRWS(I1)
!          END IF
! 
!             read(1,*) (ach(lm,i1,ntim),lm=1,lmmax_file)
! 
!             WRITE (6,FMT=9000) I1, (ACH(1,i1,ntim)/DSQRT(4.D0*PI))
!             do ispin=1,nspin
!                ipot= (i1-1)*nspin +ispin
!                do l=0,lmax
!                   do m=-l,l
!                      LM = L*L + L + M + 1
!                      IF (L.EQ.0)
!      $                    V(1,1,IPOT) = V(1,1,IPOT) + ach(1,i1,ntim)
!                      DO I = 2,IRS1
!                         V(I,LM,IPOT) = V(I,LM,IPOT) +
!      $                       (-R(I,I1))**L*ach(lm,i1,ntim)
!                      end do
!                   end do
!                end do
!             end do
!          end do
! c         call fxdrcls(ihandle)
!          close (1)
!          first_call(ntim) = .false.
!       end if
! c
c            
      DO 120 I1=1,NATOM
         iatyp = i1
c Calculate A(I1,I2,lm1,lm2), b(i1,i2,lm1) for fixed I1 !
         CALL AMN2001(I1,AMAT,ALAT,BMAT,IOPER,LMAX,NATPER,ND,NSHELL,
     &     RM,YRA,WTYR,RIJ,IJEND,CLEB,ICLEB,IEND)
c     

         DO 110 L = 0,LMAX
            DO 100 M = -L,L
               LM = L*L + L + M + 1
               AC = 0.0D0
               DO 50 I2 = 1,NATOM
c---  > determine the reference and representive atom index
                  NREF = IREF(I2)
                  IATYP2 = I2
                  DO 40 LM2 = 1,LMMAX
c---  > take moments of mt sphere and interstial
                     SUM = CMOM(LM2,IATYP2) - CMOM(LM2,NREF) 
                     if (kshape .gt. 0 ) then
                        sum = sum + 
     +                       CMINST(LM2,IATYP2) - CMINST(LM2,NREF)
                     end if
                     AC = AC + AMAT(I2,LM,LM2)*SUM
 40               CONTINUE
                  AC = AC + BMAT(I2,LM)* (Z(IATYP2)-Z(NREF))
 50            CONTINUE
               nref = iref(i1)
               AC = AC + ACH(LM,NREF,ntim)


               IF (LM.EQ.1) WRITE (6,FMT=9000) I1, (AC/DSQRT(4.D0*PI))
c     
c---  > add to v the intercell-potential
c     
               DO 90 ISPIN = 1,NSPIN
c     
c---  > determine the right potential number
c     
                  IPOT = NSPIN* (IATYP-1) + ISPIN
c     
c---  > in the case of l=0 : r(1)**l is not defined
c     
                  IF (L.EQ.0) V(1,1,IPOT) = V(1,1,IPOT) + AC
                  DO 80 I = 2,IRS1
                     V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IATYP))**L*AC
 80               CONTINUE
 90            CONTINUE

 100        CONTINUE  ! M
 110     CONTINUE     ! L
 120  CONTINUE        ! IATYP



 9000 FORMAT (1x,'spherically averaged intercell-potential for shell',
     +     i2,' :',1p,d14.6)
      END
