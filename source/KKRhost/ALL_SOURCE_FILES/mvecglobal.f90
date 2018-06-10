SUBROUTINE mvecglobal(it,iq,natyp,qmphi,qmtet,  &
        mvevi,mvevil,mvevief,  &
        natypd,lmaxd,nmvecmax)
!   ********************************************************************
!   *                                                                  *
!   *  this routine has been build up from the last part of the        *
!   *  original Munich CALCMVEC routine.                               *
!   *  on exit, MVEVI,MVEVIL,MVEVIEF are in the CARTESIAN GLOBAL       *
!   *                                           frame of reference     *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

!Parameter definitions
INTEGER LMAXDLOC
PARAMETER (LMAXDLOC=8)
COMPLEX*16 CI,CZERO
PARAMETER (CI = (0.0D0,1.0D0), CZERO = (0.0D0,0.0D0))

!Scalar Arguments
INTEGER IT,IQ,NATYP,NATYPD,LMAXD,NMVECMAX
DOUBLE PRECISION QMPHI,QMTET

!Array Arguments
DOUBLE COMPLEX MVEVI(NATYPD,3,NMVECMAX), &
           MVEVIL(0:LMAXD,NATYPD,3,NMVECMAX) 
DOUBLE COMPLEX MVEVIEF(NATYPD,3,NMVECMAX)

!Local Scalars
INTEGER ICALL,I,J,K,L,IMV,NMVEC
DOUBLE COMPLEX CS
DOUBLE COMPLEX AMIN,APLS
DOUBLE PRECISION PI,MV,MVX,MVXY,MVY,MVZ,WSQ2

!Local Arrays
DOUBLE COMPLEX USC(3,3),DROT4(4,4),W3X3(3,3)
DOUBLE COMPLEX MVG(3,NMVECMAX),MVGEF(3,NMVECMAX)
DOUBLE COMPLEX MVGL(0:LMAXD,3,NMVECMAX)
DOUBLE PRECISION MROT(3,3),FACT(0:100)
DOUBLE PRECISION MVGLO(3,NMVECMAX),MVGLOL(0:LMAXD,3,NMVECMAX)
DOUBLE PRECISION MVPHI(NMVECMAX),MVTET(NMVECMAX)
CHARACTER*1  TXTL(0:LMAXDLOC)

!Intrinsic Functions
INTRINSIC ABS,ACOS,ATAN,DREAL,DIMAG,DCONJG

!External subroutines
EXTERNAL CALCROTMAT      

!Data Statements
DATA ICALL /0/

!Save Statements
SAVE ICALL,NMVEC,USC,FACT,TXTL,PI

icall = icall + 1
!=======================================================================
IF (icall == 1) THEN
  
  IF (lmaxd > lmaxdloc) THEN
    WRITE (6,*)
    WRITE (6,*) ' Please increase parameter LMAXDLOC to ',lmaxd
    WRITE (6,*) ' in the < MVECGLOBAL > routine.'
    STOP ' < TBKKR2 > '
  END IF
  
  txtl(0) = 's'
  txtl(1) = 'p'
  txtl(2) = 'd'
  IF( lmaxd >= 3 ) THEN
    DO l=3,lmaxd
      txtl(l) = CHAR( ICHAR('f') + l - 3 )
    END DO
  END IF
  
  WRITE (1337,'(78(1H#))')
  WRITE (1337,99001)
  WRITE (1337,'(78(1H#))')
  WRITE (1337,*)
  WRITE (1337,99002)
  
  nmvec = 2
  pi = 4.d0 * ATAN(1.d0)
  
  fact(0) = 1.0D0
  DO i=1,100
    fact(i) = fact(i-1) * DBLE(i)
  END DO
!-----------------------------------------------------------------------
!  create transformation matrix   U  cartesian/sperical ccordinates
!-----------------------------------------------------------------------
!  RC,RCP  vectors in cartesian coordinates
!  RS,RSP  vectors in spherical coordinates
!         RS  = USC * R!..                          (4.40)
!         RSP = MS  * RS                                 (4.37)
!     MS(i,j) = D(j,i)                                   (4.42)
!     D  rotation matrix for complex spherical harmonics
  
! ordering of: m=-1,0,+1 >>> row 1 and 3 interchanged compared to (4.44)
  
  wsq2 = 1.0D0/SQRT(2.0D0)
  
  usc(1,1) = wsq2
  usc(1,2) = -ci*wsq2
  usc(1,3) = 0.0D0
  usc(2,1) = 0.0D0
  usc(2,2) = 0.0D0
  usc(2,3) = 1.0D0
  usc(3,1) = -wsq2
  usc(3,2) = -ci*wsq2
  usc(3,3) = 0.0D0
!-----------------------------------------------------------------------
END IF
!=======================================================================

!-----------------------------------------------------------------------
!   create the rotation matrices  DROT4 for complex spherical harmonics
!-----------------------------------------------------------------------

CALL calcrotmat(2,1,qmphi,qmtet,0.0D0,drot4,fact,4)

!-----------------------------------------------------------------------
! create the rotation matrix  MROT for vectors in cartesian coordinates
! NOTE:  U^+ D^T U gives the inverse of the real matrix  M
!        for that reason  the transposed matrix is stored as  MROT(J,I)
!-----------------------------------------------------------------------

DO i = 1,3
  DO j = 1,3
    cs = 0.0D0
    DO k = 1,3
      cs = cs + drot4(k+1,i+1)*usc(k,j)
    END DO
    w3x3(i,j) = cs
  END DO
END DO

DO i = 1,3
  DO j = 1,3
    cs = 0.0D0
    DO k = 1,3
      cs = cs + DCONJG(usc(k,i))*w3x3(k,j)
    END DO
    IF ( DIMAG(cs) > 1D-8 ) WRITE (*,*) ' MROT',i,j,cs, ' ???????????'
!     see above >> MROT(I,J) = DREAL(CS)
    mrot(j,i) = dreal(cs)
  END DO
END DO
!-----------------------------------------------------------------------

! **********************************************************************
DO imv = 1,nmvec
!-----------------------------------------------------------------------
!     transform from (+,-,z) to cartesian coordinates  (x,y,z)
!     note the convention
!-----------------------------------------------------------------------
  apls = mvevi(it,1,imv)
  amin = mvevi(it,2,imv)
  mvevi(it,1,imv) = (amin+apls)*0.5D0
  mvevi(it,2,imv) = (amin-apls)*0.5D0*ci
  
  apls = mvevief(it,1,imv)
  amin = mvevief(it,2,imv)
  mvevief(it,1,imv) = (amin+apls)*0.5D0
  mvevief(it,2,imv) = (amin-apls)*0.5D0*ci
  
  DO l = 0,lmaxd
    apls = mvevil(l,it,1,imv)
    amin = mvevil(l,it,2,imv)
    mvevil(l,it,1,imv) = (amin+apls)*0.5D0
    mvevil(l,it,2,imv) = (amin-apls)*0.5D0*ci
  END DO
!-----------------------------------------------------------------------
!     transform from LOCAL cartesian coordinates (x,y,z)
!               to  GLOBAL cartesian coordinates
!-----------------------------------------------------------------------
  DO i = 1,3
    mvg(i,imv) = czero
    mvgef(i,imv) = czero
    DO j = 1,3
      mvg  (i,imv) = mvg  (i,imv) + mrot(i,j)*mvevi(it,j,imv)
      mvgef(i,imv) = mvgef(i,imv) + mrot(i,j)*mvevief(it,j,imv)
    END DO
    mvglo(i,imv) = DIMAG(mvg(i,imv))
    
    DO l = 0,lmaxd
      mvgl(l,i,imv) = czero
      DO j = 1,3
        mvgl(l,i,imv) = mvgl(l,i,imv) + mrot(i,j)*mvevil(l,it,j,imv)
      END DO
      mvglol(l,i,imv) = DIMAG(mvgl(l,i,imv))
    END DO
    
  END DO
! ......................................................................
  DO i = 1,3
    mvevi(it,i,imv)   = mvg  (i,imv)
    mvevief(it,i,imv) = mvgef(i,imv)
    
    DO l = 0,lmaxd
      mvevil(l,it,i,imv)   = mvgl  (l,i,imv)
    END DO
  END DO
!-----------------------------------------------------------------------
!        calculate the angles
!-----------------------------------------------------------------------
  mvx = mvglo(1,imv)
  mvy = mvglo(2,imv)
  mvz = mvglo(3,imv)
  
  mv = SQRT(mvx**2+mvy**2+mvz**2)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF ( mv < 1D-8 ) THEN
    mvphi(imv) = 0D0
    mvtet(imv) = 0D0
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ELSE
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    mvxy = SQRT(mvx**2+mvy**2)
! ======================================================================
    IF ( ABS(mvxy) < 1D-8 ) THEN
      mvphi(imv) = 0D0
! ======================================================================
    ELSE
! ======================================================================
      IF ( mvy >= 0D0 ) THEN
        mvphi(imv) = ACOS(mvx/mvxy)
      ELSE IF ( mvx < 0D0 ) THEN
        mvphi(imv) = pi + ACOS(-mvx/mvxy)
      ELSE
        mvphi(imv) = 2*pi - ACOS(mvx/mvxy)
      END IF
      mvphi(imv) = mvphi(imv)*180D0/pi
      IF ( ABS(mvphi(imv)-360.0D0) < 1D-8 ) mvphi(imv) = 0D0
    END IF
! ======================================================================
    IF (mvphi(imv) >= 345.d0) mvphi(imv) = 360.d0 - mvphi(imv)
    mvtet(imv) = ACOS(mvz/mv)*180D0/pi
  END IF
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
END DO
! **********************************************************************
! output vector components, in and out angles
! ----------------------------------------------------------------------
l = 0
WRITE (1337,99003) it,iq,txtl(l),((mvglol(l,i,imv),i=1,3),imv=1,2)
WRITE (1337,99004) (txtl(l),((mvglol(l,i,imv),i=1,3),imv=1,2), l=1,lmaxd)
WRITE (1337,99005) ((mvglo(i,imv),i=1,3),imv=1,2)
WRITE (1337,99006) qmphi,qmtet,(mvphi(imv),mvtet(imv),imv=1,2)
! ----------------------------------------------------------------------
IF ( it < natyp ) THEN
  WRITE (1337,'(3X,75(1H=))')
ELSE
  WRITE (1337,*)
  WRITE (1337,'(78(1H#))')
END IF
! ----------------------------------------------------------------------

99001 FORMAT (15X,'vectorial magnetic properties given with respect',/,  &
    15X,'   to the GLOBAL (crystal) frame of reference')
99002 FORMAT (29X,'m_spin',27X,'m_orb',/, 3X,'ATOM/SITE     ',  &
    '    x         y         z',8X,'    x         y         z',/, 3X,75(1H=))
99003 FORMAT (3X,i3,'/',i3,2X,a1,' =',3F10.5,3X,3F10.5)
99004 FORMAT (12X,a1,' =',3F10.5,3X,3F10.5)
99005 FORMAT (12X,66(1H-),/,12X,'sum',3F10.5,3X,3F10.5,/)
99006 FORMAT (3X,'angles (IN)   TET =',f9.4,' PHI =',f9.4,/,  &
    3X,'angles (calc) TET =',f9.4,' PHI =',f9.4, '   TET =',f9.4,' PHI =',f9.4)
END SUBROUTINE mvecglobal
