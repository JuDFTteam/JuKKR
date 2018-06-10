! ************************************************************************
SUBROUTINE brysh3(y,x,z,xsme,ins,irmin,irc,natps,  &
    natyp,nspin,imap,lmpot,lsmear)
!*********************************************************************
!     shifts the density or potential of all mt-cell into one single
!     vector and projects out the coulomb part only.

!                                    s. bluegel , kfa , 1987

! ------------------------------------------------------------------------
!.. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
!..
!.. Scalar Arguments ..
      INTEGER IMAP,INS,LMPOT,NATPS,NATYP,NSPIN,LSMEAR
!..
!.. Array Arguments ..
      DOUBLE PRECISION X(IRMD,*),Y(*),Z(IRMIND:IRMD,LMPOTD,*), &
                       XSME(IRMD,*)
      INTEGER IRC(*),IRMIN(*)
!..
!.. Local Scalars ..
      INTEGER IA,IP,IR,IRC1,IRMIN1,IS,LM
!     ..
imap = 0
DO  is = 1,nspin
  DO  ia = natps,natyp
    ip = nspin* (ia-1) + is
    irc1 = irc(ia)
    DO  ir = 1,irc1
      imap = imap + 1
      y(imap) = x(ir,ip)
    END DO
    
!     SMEARed spherical potential
    IF ( lsmear > 0 ) THEN
      DO ir = 1,irc1
        imap = imap + 1
        y(imap) = xsme(ir,ip)
      END DO
    END IF
    
    IF (ins > 0 .AND. lmpot > 1) THEN
      irmin1 = irmin(ia)
      DO  lm = 2,lmpot
        DO  ir = irmin1,irc1
          imap = imap + 1
          y(imap) = z(ir,lm,ip)
        END DO
      END DO
    END IF
    
  END DO
END DO

END SUBROUTINE brysh3
