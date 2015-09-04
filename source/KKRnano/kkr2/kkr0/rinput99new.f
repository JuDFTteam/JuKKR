!>    Reads atominfo and rbasis files.
      SUBROUTINE RINPUTNEW99(RBASIS, NAEZ)
      IMPLICIT NONE

!     .. Array Arguments ..
!      INTEGER IRNS_dummy,KFGdummy(4),LMXCdummy,
!     &        NTCELL(naez),CLS, refpot

!      DOUBLE PRECISION Z(*),MTFACdummy,RBASIS(3,*),radius_muffin_tin(*)
      DOUBLE PRECISION, intent(out) :: RBASIS(3,*)
      INTEGER, intent(in) :: NAEZ
      
!      double precision temp
!       INTEGER :: NREF, NCLS
                                 ! atom types located at a given site
!-----------------------------------------------------------------------
!     .. Local Scalars ..
      INTEGER I!,J

!------------ array set up and definition of input parameter -----------

      WRITE (6,2004) 'Jun 2013'

!      The file 'atominfo' is no longer needed in KKRnano (Z and RMT are read from 'potential'
!      and NTCELL is read from 'shapefun')

!      OPEN(77,FILE='atominfo',FORM='formatted')
!
!      DO I=1,NAEZ
!
!                           READ (UNIT=77,FMT=*)    Z(I),
!     +                        LMXCdummy,
!     +                       (KFGdummy(J),J=1,4),
!     +                        CLS,
!     +                        REFPOT,
!     +                        NTCELL(I),
!     +                        MTFACdummy,
!     +                        IRNS_dummy,
!     +                        temp
!
!      radius_muffin_tin(I) = temp
!
!      END DO
!      CLOSE (77)

!----------------------------------------------------------------------
!
      OPEN(77,FILE='rbasis',FORM='formatted')
      DO I=1,NAEZ
             READ (UNIT=77,FMT=*) RBASIS(1:3,I)
      ENDDO                         
      CLOSE (77)

!     WRITE(6,2110)
!     WRITE(6,2103)
      WRITE(6,fmt="( 3(7(1H-),1H+) ,55(1H-))")
      WRITE(6,fmt="(10(3(1H-),1H+) ,39(1H-))")

! *********************************************Input-End ********
!
!     No cleanup of Format statements due to nostalgic reasons.
!
 2004 FORMAT( /80(1H=)/
     & '|',78X,'|'/
     & '|',1X,'KKRnano ',69X,'|'/
     & '|',1X,'Massively Parallel Screened Korringa-Kohn-Rostoker ',
     &         'Electronic Structure Code',1X,'|'/
     & '|',1X,'for Bulk',69X,'|'/
     & '|',78X,'|'/
     & '|',1X,'established Juelich 2008',26X,
     & 'Version : ',A8,9X,'|'/
     & '|',78X,'|'/80(1H=))
 2010 FORMAT(' NSPIN '/I4)
 2011 FORMAT(' NSTEPS'/I4)
 2014 FORMAT('          ALAT = ',F15.8)
 2015 FORMAT('   INTERVX   INTERVY   INTERVZ'/3I10)
 2016 FORMAT('    NCLS    NREF   '/,2I8)
 2018 FORMAT(' RBASIS'/,
     &     'SITE                BASIS VECTORS                 ')
 2019 FORMAT('         ABASIS         BBASIS         CBASIS'/3F15.8)
 2025 FORMAT((i4,3F15.8))
 2028 FORMAT(' NAEZ ',/,I8)
! ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2101 format(   3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
 2102 format( 3(9(1H-),1H+) ,49(1H-))
 2103 FORMAT(10(3(1H-),1H+) ,39(1H-))
 2104 format(   3(1H-),1H+  ,75(1H-))
 2107 format( 3(14(1H-),1H+),34(1H-))
 2110 format( 3(7(1H-),1H+) ,55(1H-))
 9000 FORMAT (I2,3X,4I5)
 9010 FORMAT (1X,I1,1X,4I1)
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,
     +       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,
     +       'NAEZ  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,
     +       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',
     +       f8.5)
 9040 FORMAT (3f12.7)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9060 FORMAT (8i4)
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/,
     +        ' convergence quality required :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3,
     +       ' is used up to iteration-      ',/,20x,'depth :',i3,
     +       '  then jacobian is fixed and potential      ',/,20x,
     +       'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts ',
     +       'the parameter lmaxd has to be set equal lmax ')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation ',
     +       '- cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43,/,21x,' using',i3,'-th. born approximation ')
 9170 FORMAT (21x,a43)
 9180 FORMAT (2i5)
 9190 FORMAT (3f12.7,/,4i4)
 9200 FORMAT (3i4,1f12.7)
 9210 FORMAT (' lmax'/,i4)
 9220 FORMAT ('          E1          E2          TK'/,3f12.6)
 9230 FORMAT ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
 9250 FORMAT ('  IFILE    IPE ISHIFT'/,3i7)
 9260 FORMAT (' KSHAPE    IRM    ICST'/,3i7)
 9270 FORMAT ('   KCOR  KVREL    KWS   KHFELD    KXC'/,5i7)
 9280 FORMAT (' external magnetic hfield     :',f15.4/,
     +        ' VCONST                       :',f15.6)
 9290 FORMAT ('   IMIX   '/,i7)
 9300 FORMAT (' ITDBRY'/,i7)
 9310 FORMAT ('      STRMIX        FCM       QBOUND'/,3f12.6)
 9320 FORMAT ('      BRYMIX'/,f12.6)
 9330 FORMAT ('    KTE   KPRE   KVMAD'/,3i7)
 9301 format(   3(1H-),1H+  ,75(1H-))
 9302 format( 3(11(1H-),1H+),43(1H-))
 9303 format(3(6(1H-),1H+) ,58(1H-))
 9304 format(4(6(1H-),1H+) ,51(1H-))
 9305 format(3(6(1H-),1H+),11(1H-),1H+ ,46(1H-))
 9306 format(6(6(1H-),1H+) ,37(1H-))
 9307 format(6(1H-),1H+,72(1H-))
 9308 format(11(1H-),1H+,67(1H-))
 9309 format(5(6(1H-),1H+) ,44(1H-))
!
      END
