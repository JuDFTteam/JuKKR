!------------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> 
!------------------------------------------------------------------------------------
      MODULE MOD_MIXBROYDEN
      INTEGER ::  MIT=1
       CONTAINS
! c*********************************************************************
   !-------------------------------------------------------------------------------
   !> Summary: Broyden mixing wrapper subroutine
   !> Author: 
   !> Category: KKRimp, potential
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> 
   !-------------------------------------------------------------------------------
      SUBROUTINE MIXBROYDEN(VPOT,VPOT_OUT,INS,MIXFAC,NSPIN,CELL,LMAXATOM, &
                        NATOM,ITDBRY,IMIX,LMPOTD,IRMD)
! c*********************************************************************
      USE TYPE_CELL
      IMPLICIT NONE
     REAL(8)            :: VPOT(IRMD,LMPOTD,NSPIN,NATOM)
     REAL(8)            :: VPOT_OUT(IRMD,LMPOTD,NSPIN,NATOM)
     INTEGER            :: INS
     REAL(8)            :: MIXFAC
     INTEGER            :: NSPIN
     TYPE(CELL_TYPE)    :: CELL(NATOM)
     INTEGER            :: LMAXATOM(NATOM)
     INTEGER            :: NATOM
     INTEGER            :: ITDBRY
     INTEGER            :: IMIX
     INTEGER            :: LMPOTD
     INTEGER            :: IRMD


     INTEGER               :: NTIRD
     REAL(8)   :: DUMMY(100)
     integer,parameter     :: iobroy=23420
     integer,save         :: first=1
! write(*,*) '1'
     if (first==1) then
       first=0
       OPEN (IOBROY,FORM='unformatted',STATUS='unknown',file='storage_broy1')
       OPEN (IOBROY+2,FORM='unformatted',STATUS='unknown',file='storage_broy2')
       mit=1
     end if
!     writE(*,*) 'first',first
!     writE(*,*) 'mit',mit
! stop
! ntird=100000000
! write(*,*) '2'
     CALL MIXBROYDEN_SHIFT(DUMMY,'.',VPOT_OUT,INS,CELL,NATOM,NSPIN,NTIRD,LMAXATOM,LMPOTD,IRMD,1)
! write(*,*) 'ntird= ',ntird
! write(*,*) '3'
! write(*,*) 'ntird',ntird

! write(*,*) '4',ntird
! write(*,*) INS,MIXFAC,NSPIN,LMAXATOM, &
!                         NATOM,ITDBRY,IMIX,NTIRD,LMPOTD,IRMD,iobroy
     CALL MIXBROYDEN_START(VPOT,VPOT_OUT,INS,MIXFAC,NSPIN,CELL,LMAXATOM, &
                        NATOM,ITDBRY,IMIX,NTIRD,LMPOTD,IRMD,iobroy)
! write(*,*) '5000',ntird
     END SUBROUTINE MIXBROYDEN
! c*********************************************************************



  
! c*********************************************************************
   !-------------------------------------------------------------------------------
   !> Summary: Broyden mixing of the potential
   !> Author: S. Bluegel May 1987, Jul 1989; B. Drittler Aug 1988
   !> Category: KKRimp, potential
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> 
   !>   imix :
   !>     3      broyden's            f i r s t  m e t h o d
   !>     4      broyden's          s e c o n d  m e t h o d
   !>     5      anderson's     g e n e r a l i z e d   m e t h o d
   !>
   !>   implemented here according to notes of s.b.
   !>   broyden's iteration scheme following the papers of :
   !>   srivastava, j. phys. , 17 (1984) , pp l317
   !>   c.g. broyden in math.comput., 19 , pp 577, 1965
   !>   c.g. broyden in ibid, 21 ,pp 368 ,1967
   !>   the method has been generalized to include a metric.the
   !>   definition of the necessary inner products are similar to the
   !>   discription given in the notes of m.weinert.the algorithm
   !>   discribed in the paper srivastava  has been simplified
   !>   ( see notes of s.b.)
   !>   the files ui,vi are stored on high speed ssd memory.
   !>   broyden's update treats charge and spin on the same footing
   !>                s. bluegel , kfa , may 1987
   !>   the anderson method (d.g. anderson, j. acm 12, 547 (1964)) has
   !>   been generalized and reformulated as an improvement of broyden's
   !>   second method. successive linesearch is replaced by successive
   !>   search on hyperplanes. ( see notes of s.b. )
   !>                s. bluegel , issp , july 1989
   !>
   !>   modified for non spherical potential
   !>                  b. drittler , aug. 1988
   !-------------------------------------------------------------------------------
      SUBROUTINE MIXBROYDEN_START(VPOT,VPOT_OUT,INS,MIXFAC,NSPIN,CELL,LMAXATOM, &
                        NATOM,ITDBRY,IMIX,NTIRD,LMPOTD,IRMD,iobroy)
! C     .. Parameters ..
      USE TYPE_CELL
      use mod_types, only: t_inc
      IMPLICIT NONE
! !       INCLUDE 'inc.p'
     REAL(8)            :: VPOT(IRMD,LMPOTD,NSPIN,NATOM)
     REAL(8)            :: VPOT_OUT(IRMD,LMPOTD,NSPIN,NATOM)
     INTEGER            :: INS
     REAL(8)            :: MIXFAC
     INTEGER            :: NSPIN
     TYPE(CELL_TYPE)    :: CELL(NATOM)
     INTEGER            :: LMAXATOM(NATOM)
     INTEGER            :: NATOM
     INTEGER            :: ITDBRY
     INTEGER            :: IMIX
     INTEGER            :: NTIRD
     INTEGER            :: LMPOTD
     INTEGER            :: IRMD
     INTEGER            :: IOBROY

      INTEGER ITDTHD
      PARAMETER (ITDTHD=40)
     INTEGER,parameter :: IPF=1337
real(8),parameter      :: zero=0.0D0
real(8),parameter      :: one=1.0D0
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
!       INTEGER IRMIND
!       PARAMETER (IRMIND=IRMD-IRNSD)
!       INTEGER NSPINDD
!       PARAMETER (NSPINDD=2*KREL + (1-KREL)*NSPIND)
!       INTEGER NTIRD
!       PARAMETER (NTIRD= (IRMD*NTPERD+ (IRNSD+1)* (LMPOTD-1)*NATYPD)*
!      +          NSPINDD)

! C     ..
! C     .. Array Arguments ..
!       DOUBLE PRECISION ATWGHT(*),DRDI(IRMD,*),R(IRMD,*),
!      +                 V(IRMD,LMPOTD,*),VINS(IRMIND:IRMD,LMPOTD,*),
!      +                 VISP(IRMD,*),VSPSMO(IRMD,*),VSPSME(IRMD,*)
!       INTEGER IRC(*),IRMIN(*)
! C     ..
! C     .. Local Scalars ..
      DOUBLE PRECISION CMM,RMIXIV,SMNORM,VMDENO,VMNORM,VOLINV
      INTEGER IA,IJ,IMAP,IR,IRC1,IRMIN1,ISP,IT,LM
! C     ..
! C     .. External Functions ..
      DOUBLE PRECISION DDOT
!       EXTERNAL DDOT
! C     ..
! C     .. External Subroutines ..
!       EXTERNAL BRYSH1,BRYSH2,BRYSH3,DAXPY,DSCAL,RCSTOP
! C     ..
! C     .. Intrinsic Functions ..
!       INTRINSIC ABS
! C     ..
! C     .. Save statement ..
!       SAVE MIT,ZERO,ONE,WIT
! C     ..
! C     .. Local Arrays ..
      DOUBLE PRECISION,allocatable :: AM(:),BM(:),FM(:),WIT(:) 
!       DOUBLE PRECISION AM(2:ITDTHD-1),BM(2:ITDTHD-1),FM(NTIRD),WIT(2:200) 


 DOUBLE PRECISION,allocatable :: FM1(:),G(:),SM(:),SM1(:), &
                       VI3(:),UI2(:),UI3(:),VI2(:)

allocate(AM(2:ITDTHD-1),BM(2:ITDTHD-1),FM(NTIRD),WIT(2:200))


allocate(FM1(NTIRD),G(NTIRD),SM(NTIRD),SM1(NTIRD),VI3(NTIRD),UI2(NTIRD),UI3(NTIRD),VI2(NTIRD))

! write(*,*) '4.1'! C     ..
! C     .. Scalar Arguments ..
!       DOUBLE PRECISION MIXFAC
!       INTEGER IMIX,INS,IOBROY,IPF,ITDBRY,LMPOT,NATPS,NATOM,NSPIN
! C     ..
! C     .. Data statements ..
!       DATA MIT/1/,ZERO,ONE/0.0D0,1.0D0/
! C     ..
!       READ(28,FMT='(I5)') MIT
!       REWIND 28
      IF (ITDBRY.GT.ITDTHD .OR. ITDTHD.GT.200) STOP 'ITDBRY  '

      IF (IMIX.LE.2 .OR. IMIX.GT.5) STOP 'IMIXD   '

      IF (MIT.GT.ITDBRY) MIT = 1
      IF (IMIX.EQ.3 .and. t_inc%i_write>0)  WRITE (IPF,FMT='('' broyden"s 1st method used '')')
      IF (IMIX.EQ.4 .and. t_inc%i_write>0)  WRITE (IPF,FMT='('' broyden"s 2nd method used '')')
      IF (IMIX.EQ.5 .and. t_inc%i_write>0)  WRITE (IPF,FMT= &
          '('' generalized anderson method used '')')
      if (t_inc%i_write>0) WRITE(IPF,'(A,i4)') ' Iteration index (read in):',MIT
! c
! write(*,*) '4.2'
      RMIXIV = ONE/MIXFAC
! c
! c---->  the following block is activated only one iteration before
! c        broyden iteration scheme is used
! c---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
! c                   metric  g := r*r*drdi
! c---->  map data of all muffin-tin spheres into one single vector
! c
! c
! write(*,*) '5',NTIRD
      CALL MIXBROYDEN_SHIFT(SM1,'<',VPOT,INS,CELL,NATOM,NSPIN,IMAP,LMAXATOM,LMPOTD,IRMD,NTIRD)
!       CALL BRYSH3(SM1,VISP,VINS,VSPSME,INS,IRMIN,IRC,NATPS, &
!                   NATOM,NSPIN,IMAP,LMPOT,LSMEAR)

! write(*,*) '6',IMAP

      CALL MIXBROYDEN_SHIFT(FM1,'<',VPOT_OUT,INS,CELL,NATOM,NSPIN,IMAP,LMAXATOM,LMPOTD,IRMD,NTIRD)
!       CALL BRYSH1(FM1,V,VSPSMO,INS,IRMIN,IRC,NATPS, &
!                   NATOM,NSPIN,IMAP,LMPOT,LSMEAR)

! write(*,*) '7',IMAP
      IF (IMAP.GT.NTIRD) STOP 'NIRDBRY '

      DO IJ = 1,IMAP
        FM1(IJ) = RMIXIV* (FM1(IJ)-SM1(IJ))
      END DO !IJ
! c
! write(*,*) '8'
      IJ = 0
      DO ISP = 1,NSPIN
        DO IA = 1,NATOM
          IRC1 = CELL(IA)%NRMAX !IRC(IA)
          VOLINV = 3.0D0/ (CELL(IA)%RMESH(IRC1)**3) !R(IRC1,IA)**3)
          DO IR = 1,IRC1
            IJ = IJ + 1
            G(IJ) = VOLINV*(CELL(IA)%RMESH(IR)**2)*CELL(IA)%DRMESHDI(IR)
!             G(IJ) = ATWGHT(IA)*VOLINV*R(IR,IA)*R(IR,IA)*DRDI(IR,IA)
          END DO !IR
! C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
          IF (INS.GE.1 .AND. LMAXATOM(IA).GT.0) THEN
            IRMIN1 = CELL(IA)%NRMIN_NS !IRMIN(IA)
            DO LM = 2,(2*LMAXATOM(IA)+1)**2
              DO IR = IRMIN1,IRC1
                IJ = IJ + 1
                  G(IJ) = VOLINV*(CELL(IA)%RMESH(IR)**2)*CELL(IA)%DRMESHDI(IR)
!                 G(IJ) = ATWGHT(IA)*VOLINV*R(IR,IA)*R(IR,IA)*DRDI(IR,IA)
              END DO !IR
            END DO !LM
          END IF
        END DO !IA
      END DO !ISP
! c
! write(*,*) '9'

      IF (MIT.GT.1) THEN
        REWIND IOBROY + 2
        READ (IOBROY+2) (SM1(IJ),IJ=1,IMAP), (FM1(IJ),IJ=1,IMAP)

! c
! c----> map rho(m) of all mt-spheres into one single vector
! c
!  write(*,*) '10'
       CALL MIXBROYDEN_SHIFT(SM,'<',VPOT,INS,CELL,NATOM,NSPIN,IMAP,LMAXATOM,LMPOTD,IRMD,NTIRD)
!         CALL BRYSH3(SM,VISP,VINS,VSPSME,INS,IRMIN,IRC,NATPS, &
!                     NATOM,NSPIN,IMAP,LMPOT,LSMEAR)
! c
! c----> map f[m] = f(m) - rho(m) = f(rho(m)) - rho(m) of all mt-spheres
! c      into one single vector
! c
! write(*,*) '11'
        CALL MIXBROYDEN_SHIFT(FM,'<',VPOT_OUT,INS,CELL,NATOM,NSPIN,IMAP,LMAXATOM,LMPOTD,IRMD,NTIRD)
!         CALL BRYSH1(FM,V,VSPSMO,INS,IRMIN,IRC,NATPS,&
!                     NATOM,NSPIN,IMAP,LMPOT,LSMEAR)
        DO 70 IJ = 1,IMAP
          FM(IJ) = RMIXIV* (FM(IJ)-SM(IJ))
   70   CONTINUE
! c
! c----> calculate  sm = rho(m) - rho(m-1)
! c----> calculate dfm = f[m] - f[m-1]
! c
        DO 80 IJ = 1,IMAP
          SM1(IJ) = SM(IJ) - SM1(IJ)
          FM1(IJ) = FM(IJ) - FM1(IJ)
   80   CONTINUE
! c
! c----> loop to generate u[m] = u(ij,mit)
! c
        DO 90 IJ = 1,IMAP
          UI3(IJ) = MIXFAC*FM1(IJ) + SM1(IJ)
   90   CONTINUE
        REWIND IOBROY
        DO 100 IT = 2,MIT - 1
            READ (IOBROY) (UI2(IJ),IJ=1,IMAP), (VI2(IJ),IJ=1,IMAP),WIT

          AM(IT) = DDOT(IMAP,FM1,1,VI2,1)
          CALL DAXPY(IMAP,-AM(IT),UI2,1,UI3,1)
  100   CONTINUE
! c
! c----> print amj = the importance of the history of ui
! c
        if (t_inc%i_write>0) WRITE (IPF,FMT='(5x,'' amj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')&
          MIT - 1, (AM(IT),IT=2,MIT-1)
! c
! c
        IF (IMIX.EQ.3) THEN
! c-------->     b r o y d e n ' s   f i r s t   m e t h o d
! c
! c
! c----> calculate dsmnorm
! c
          SMNORM = ZERO
          DO 110 IJ = 1,IMAP
            SMNORM = SMNORM + SM1(IJ)*G(IJ)*SM1(IJ)
  110     CONTINUE
! c
! c----> convolute dsm with the metric g
! c
          DO 120 IJ = 1,IMAP
            SM1(IJ) = G(IJ)*SM1(IJ)
  120     CONTINUE
! c
! c----> loop to generate v[m] = v(ij,mit)
! c
          DO 130 IJ = 1,IMAP
            VI3(IJ) = MIXFAC*SM1(IJ)
  130     CONTINUE
          REWIND IOBROY
          DO 140 IT = 2,MIT - 1
              READ (IOBROY) (UI2(IJ),IJ=1,IMAP), (VI2(IJ),IJ=1,IMAP),WIT

            BM(IT) = DDOT(IMAP,SM1,1,UI2,1)
            CALL DAXPY(IMAP,-BM(IT),VI2,1,VI3,1)
  140     CONTINUE
! c
! c----> complete the evaluation of v[m]
! c
          VMDENO = DDOT(IMAP,SM1,1,UI3,1) - SMNORM

          IF (ABS(VMDENO).LT.1D-70) STOP 'BRY0SN  '

          CALL DSCAL(IMAP,ONE/VMDENO,VI3,1)
! c
! c----> print bmj = the importance of the history of vi
! c
          if (t_inc%i_write>0) WRITE (IPF,FMT='(5x,'' bmj , ---> j=2,'',i3,/,(9x,1p,7d10.2))' &
            ) MIT - 1, (BM(IT),IT=2,MIT-1)
! c
        ELSE IF (IMIX.EQ.4) THEN
! c-------->     b r o y d e n ' s   s e c o n d    m e t h o d
! c
! c----> calculate v[m] ; convoluted with the metric g
! c
          DO 150 IJ = 1,IMAP
            VI3(IJ) = G(IJ)*FM1(IJ)
  150     CONTINUE
! c
! c----> calculate #vm# and normalize v[m]
! c
          VMNORM = DDOT(IMAP,VI3,1,FM1,1)
          CALL DSCAL(IMAP,ONE/VMNORM,VI3,1)
! c
        ELSE IF (IMIX.EQ.5) THEN
! c-------->     g e n e r a l i z e d   a n d e r s o n   m e t h o d
! c
! c----> calculate v[m] ; convoluted with the metric g
! c
          DO 160 IJ = 1,IMAP
            VI3(IJ) = G(IJ)*FM1(IJ)
  160     CONTINUE
          REWIND IOBROY
          DO 170 IT = 2,MIT - 1
              READ (IOBROY) (UI2(IJ),IJ=1,IMAP), (VI2(IJ),IJ=1,IMAP),WIT

            CALL DAXPY(IMAP,-AM(IT)*WIT(IT),VI2,1,VI3,1)
  170     CONTINUE
! c
! c----> complete the evaluation of v[m]
! c
          VMDENO = DDOT(IMAP,FM1,1,VI3,1)

          IF (ABS(VMDENO).LT.1D-70) STOP 'BRY1SN  '

          CALL DSCAL(IMAP,ONE/VMDENO,VI3,1)
! c
! c----> save wit(mit) for next iteration
! c
          WIT(MIT) = VMDENO
! c
        END IF
! c
! c----> write u3(ij) and v3(ij) on disk
! c
          WRITE (IOBROY) (UI3(IJ),IJ=1,IMAP), (VI3(IJ),IJ=1,IMAP),WIT
! c
! c----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
! c
        DO 180 IJ = 1,IMAP
          FM1(IJ) = FM(IJ)
          SM1(IJ) = SM(IJ)
  180   CONTINUE
! c
! c----> calculate cmm
! c
        CMM = DDOT(IMAP,FM,1,VI3,1)
! C           WRITE (IPF,FMT='(5X,'' CMM = '',1P,D12.4)') CMM
! c
! c----> update rho(m+1)
! c
        CALL DAXPY(IMAP,ONE-CMM,UI3,1,SM,1)
! c
! c----> map solution back into each mt-sphere
! c
        CALL MIXBROYDEN_SHIFT(SM,'>',VPOT_OUT,INS,CELL,NATOM,NSPIN,IMAP,LMAXATOM,LMPOTD,IRMD,NTIRD)
!         CALL BRYSH2(SM,V,VSPSMO,INS,IRMIN,IRC,NATPS, &
!                     NATOM,NSPIN,IMAP,LMPOT,LSMEAR)
! c
      END IF
      MIT = MIT + 1
!       WRITE(28,FMT='(I5)') MIT
      REWIND IOBROY + 2
      WRITE (IOBROY+2) (SM1(IJ),IJ=1,IMAP), (FM1(IJ),IJ=1,IMAP)

      RETURN

! c  190   CALL RCSTOP('broy10  ')
! c  200   CALL RCSTOP('broy11  ')
! c  210   CALL RCSTOP('broy12  ')
! c  220   CALL RCSTOP('broy13  ')

      END SUBROUTINE MIXBROYDEN_START










   !-------------------------------------------------------------------------------
   !> Summary: Broyden mixing of the potential
   !> Author: S. Bluegel 1987; D. Bauer
   !> Category: KKRimp, potential
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> 
   !> @Bauer
   !> rewriten for KKRFLEX impurity calculatuions
   !> this is a combination of brysh1 and brysh2 routines using variable
   !> mode. If mode is '>' Y will be copied in X
   !> mode. If mode is '<' X will be copied in Y
   !> c*********************************************************************
   !> old brysh1:
   !> c*********************************************************************
   !> c     shifts the density or potential of all mt-cell into one single
   !> c     vector and projects out the coulomb part only.
   !> c                                    s. bluegel , kfa , 1987
   !> c
   !> c ------------------------------------------------------------------------
   !> old brysh2:
   !> c*********************************************************************
   !> c     maps the density or potential back from one single vector into
   !> c     the proper bins of each single mt-cell . the magnetization
   !> c     density is also added in.
   !> c                                    s. bluegel , kfa , 1987
   !> c
   !-------------------------------------------------------------------------------
      SUBROUTINE MIXBROYDEN_SHIFT(Y,mode,X,INS,CELL, &
                        NATYP,NSPIN,IMAP,LMAXATOM,LMPOTD,IRMD,NTIRD)

! c*********************************************************************
      USE TYPE_CELL
      IMPLICIT NONE
      DOUBLE PRECISION Y(NTIRD)
      CHARACTER MODE
      DOUBLE PRECISION X(IRMD,LMPOTD,NSPIN,NATYP)
      TYPE(CELL_TYPE) :: CELL(NATYP)
      INTEGER INS,NATYP,NSPIN,IMAP,LMAXATOM(NATYP),LMPOTD,IRMD
! C     .. Local Scalars ..
      INTEGER IA,IR,IRC1,IRMIN1,IS,LM,NTIRD
! C     ..
      IMAP = 0

      DO IS = 1,NSPIN
        DO IA = 1,NATYP
!           IP = NSPIN* (IA-1) + IS
          IRC1 = CELL(IA)%NRMAX !IRC(IA)
          DO IR = 1,IRC1
            IMAP = IMAP + 1
            if (mode=='<') then
              Y(IMAP) = X(IR,1,IS,IA)  !brysh1
            else if (mode=='>') then
              X(IR,1,IS,IA) = Y(IMAP)  !brysh1
            else if (mode=='.') then
              ! do nothing just determine IMAP
            else 
              stop 'mode error'
            end if
          END DO !IR
! c
          IF (INS.GT.0 .AND. LMAXATOM(IA).GT.1) THEN
            IRMIN1 = CELL(IA)%NRMIN_NS !IRMIN(IA)
            DO LM = 2,(2*LMAXATOM(IA)+1)**2 !LMPOT
              DO IR = IRMIN1,IRC1
                IMAP = IMAP + 1
                if (mode=='<') then
                    Y(IMAP) = X(IR,LM,IS,IA) !brysh1
                else if (mode=='>') then
                    X(IR,LM,IS,IA) = Y(IMAP) !brysh2
                else if (mode=='.') then
                  ! do nothing just determine IMAP
                end if
              END DO !IR
            END DO !LM
          END IF
! c
        END DO !IA
      END DO !IS
! c
  if (mode/='.') then
    if (imap>NTIRD) stop '[MIXBROYDEN_SHIFT] IMAP > NTIRD'
  end if
      END SUBROUTINE MIXBROYDEN_SHIFT

END MODULE MOD_MIXBROYDEN



