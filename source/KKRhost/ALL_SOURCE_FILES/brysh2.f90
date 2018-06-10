! ************************************************************************
SUBROUTINE brysh2(y,x,xsme,ins,irmin,irc,natps, natyp,nspin,imap,lmpot,lsmear)
!*********************************************************************
!     maps the density or potential back from one single vector into
!     the proper bins of each single mt-cell . the magnetization
!     density is also added in.
!                                    s. bluegel , kfa , 1987

! ------------------------------------------------------------------------
!.. Parameters ..
      include 'inc.p'
!      INTEGER IRMD,LPOTD
!      PARAMETER (IRMD=424,LPOTD=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
!..
!.. Scalar Arguments ..
      INTEGER IMAP,INS,LMPOT,NATPS,NATYP,NSPIN,LSMEAR
!..
!.. Array Arguments ..
      DOUBLE PRECISION X(IRMD,LMPOTD,*),Y(*),XSME(IRMD,*)
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
      x(ir,1,ip) = y(imap)
    END DO
    
!     Next for SMEARed spherical potential
    IF ( lsmear > 0 ) THEN
      DO ir = 1, irc1
        imap = imap + 1
        xsme(ir,ip) = y(imap)
      END DO
    END IF
    
    IF (ins > 0 .AND. lmpot > 1) THEN
      irmin1 = irmin(ia)
      DO  lm = 2,lmpot
        DO  ir = irmin1,irc1
          imap = imap + 1
          x(ir,lm,ip) = y(imap)
        END DO
      END DO
    END IF
    
  END DO
END DO

END SUBROUTINE brysh2
