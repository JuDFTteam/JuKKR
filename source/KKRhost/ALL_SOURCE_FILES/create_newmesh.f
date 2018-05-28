       module mod_create_newmesh

       contains

       SUBROUTINE CREATE_NEWMESH(NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
     +                     IRID,NTOTD,NFUND,NCHEBD,IRMDNEW,
     +                     NSPIN,R,IRMIN,IPAN,IRCUT,
     +                     R_LOG,NPAN_LOG,NPAN_EQ,NCHEB,
     +                     NPAN_LOGNEW,NPAN_EQNEW,
     +                     NPAN_TOT,RNEW,RPAN_INTERVALL,IPAN_INTERVALL,
     +                     NCELLD,NTCELL,THETAS,THETASNEW) !< optional arguments
       IMPLICIT NONE
       !INCLUDE 'inc.p'
       !changed interface to get rid of inc.p and to be able to use i
       !create_newmesh in tmatimp routine for GREENIMP option
       !this is the list of  array dimensions previously importted from inc.p
       INTEGER, intent(in) :: NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND,IRID,
     &                        NTOTD,NFUND,NCHEBD
       INTEGER NSPIN,IRMIN(NATYPD),IPAN(NATYPD)
       INTEGER LMMAXD
       !PARAMETER (LMMAXD= (LMAXD+1)**2)
       INTEGER LMPOTD
       !PARAMETER (LMPOTD= (LPOTD+1)**2)
       INTEGER IRMIND
       !PARAMETER (IRMIND= IRMD-IRNSD)
       INTEGER NPAN_LOG,NPAN_EQ,NCHEB,NPAN_INST,NPAN_TOT(NATYPD)
       INTEGER NPAN_LOGNEW(NATYPD),NPAN_EQNEW(NATYPD)
       INTEGER IRCUT(0:IPAND,NATYPD)
       DOUBLE PRECISION R(IRMD,NATYPD)
       DOUBLE PRECISION FAC
       PARAMETER (FAC=2d0)
       INTEGER I1,IPOT,IPOTM,IR2,IP,ICELL,
     +         ISHIFT,ILOGPANSHIFT,ILINPANSHIFT,NPAN_LOGTEMP,
     +         IMIN,IMAX,IMINNEW,IMAXNEW,LM1,IRMDNEW
       DOUBLE PRECISION R_LOG,RMIN,RMAX,RVAL
       DOUBLE PRECISION RNEW(IRMDNEW,NATYPD),
     +                  RPAN_INTERVALL(0:NTOTD,NATYPD)
       INTEGER IPAN_INTERVALL(0:NTOTD,NATYPD)
       ! optional arguments, do interpolation when given
       INTEGER, intent(in):: NCELLD ! needs to be given in any case for interface to work (choose NCELLD=1 then)
       INTEGER, intent(in), optional :: NTCELL(NATYPD)
       DOUBLE PRECISION, intent(in), optional :: 
     +       THETAS(IRID,NFUND,NCELLD)
       DOUBLE PRECISION, intent(inout), optional ::
     +       THETASNEW(NTOTD*(NCHEBD+1),NFUND,NCELLD)
       ! allocatable arrays
       double precision, allocatable :: THETASIN(:,:,:)

       LMMAXD= (LMAXD+1)**2
       LMPOTD= (LPOTD+1)**2
       IRMIND= IRMD-IRNSD

       ! checks for optional arguments
       if( present(NTCELL) .and. 
     +     (.not.present(THETAS) .or. .not.present(THETASNEW)) ) then
           write(*,*) 'Error in create_newmesh:'
           write(*,*) 'List of optional arguments not complete'
           stop 
       end if

! allocations
       if( present(NTCELL) ) allocate(THETASIN(IRID,NFUND,NCELLD))
       if( present(NTCELL) ) THETASNEW=0d0

       IPOTM=0

       DO I1 = 1,NATYPD

        IPOT=NSPIN*(I1-1)+1
        NPAN_INST= IPAN(I1)-1
        NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST


c log panel
       RMIN=R(2,I1)
       RMAX=R_LOG
       RVAL=0d0
       ISHIFT=0
       IF (R_LOG.GT.R(IRMIN(I1),I1)) THEN
         ILOGPANSHIFT=1
         ILINPANSHIFT=0
       ELSE
         ILOGPANSHIFT=0
         ILINPANSHIFT=1
       ENDIF 
  
       IF (ILINPANSHIFT.EQ.1) THEN
         write(*,*) 'ERORR: non-spherical part of the potential needs'
         write(*,*) 'to be inside the log panel'
         write(*,*) 'atom (I1):', I1
         write(*,*) 'R_LOG', R_LOG
         write(*,*) 'R(IRMIN(I1), I1)', R(IRMIN(I1),I1)
         write(*,*) 'IRMIN(I1)', IRMIN(I1)
         STOP 'Error creating newmesh' 
       ENDIF
      
       DO IP=0,NPAN_LOG-ILOGPANSHIFT
        RVAL=(FAC**IP-1d0)/(FAC**(NPAN_LOG-ILOGPANSHIFT)-1d0)       
        RPAN_INTERVALL(IP+ISHIFT,I1)= RMIN+RVAL*(RMAX-RMIN)
        IPAN_INTERVALL(IP+ISHIFT,I1)= (IP+ISHIFT)*(NCHEB+1)
        IF (ISHIFT.EQ.0.AND.
     +      RPAN_INTERVALL(IP,I1).GT.R(IRMIN(I1),I1)) THEN
         ISHIFT=1
         NPAN_LOGTEMP=IP
         RPAN_INTERVALL(IP+1,I1)=RPAN_INTERVALL(IP,I1)
         IPAN_INTERVALL(IP+1,I1)=(IP+ISHIFT)*(NCHEB+1)
         RPAN_INTERVALL(IP,I1)=R(IRMIN(I1),I1)
         IPAN_INTERVALL(IP,I1)=IP*(NCHEB+1)
        ENDIF
       ENDDO ! NPAN_LOG

c equivalent panel
       ISHIFT=0
       RMIN=R_LOG
       RMAX=R(IRCUT(1,I1),I1)
       DO IP=0,NPAN_EQ-ILINPANSHIFT
        RPAN_INTERVALL(IP+ISHIFT+NPAN_LOG,I1)=RMIN+IP*(RMAX-RMIN)/
     +          (NPAN_EQ-ILINPANSHIFT)
        IPAN_INTERVALL(IP+ISHIFT+NPAN_LOG,I1)=(NPAN_LOG+IP+ISHIFT)*
     +          (NCHEB+1)
       ENDDO ! NPAN_EQ

c intersection zone
       DO IP=1,NPAN_INST
         RPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP,I1)=R(IRCUT(IP+1,I1),I1)
         IPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP,I1)=(NPAN_LOG+NPAN_EQ+IP)*
     +           (NCHEB+1)
       ENDDO ! NPAN_INST

       NPAN_EQNEW(I1)=NPAN_EQ+NPAN_LOG-NPAN_LOGTEMP
       NPAN_LOGNEW(I1)=NPAN_LOGTEMP

       CALL CHEBMESH(NPAN_TOT(I1),NCHEB,RPAN_INTERVALL(0:,I1),
     +               RNEW(1,I1))

       ! do interpolation only when optional arguments are given
       if( present(NTCELL) ) then
c interpolate shape function THETAS to new shape function THETASNEW
c save THETAS to THETASIN
        ICELL = NTCELL(I1)
        DO LM1=1,NFUND
         THETASIN(:,LM1,ICELL)=THETAS(:,LM1,ICELL)
         IR2=0
         DO IP=NPAN_LOGNEW(I1)+NPAN_EQNEW(I1)+1,NPAN_TOT(I1)
          IR2=IR2+1
          IMIN=IRCUT(IR2,I1)+1
          IMAX=IRCUT(IR2+1,I1)
          IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
          IMAXNEW=IPAN_INTERVALL(IP,I1)
          CALL INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),
     +           THETASIN(IMIN-IRCUT(1,I1):IMAX-IRCUT(1,I1),LM1,ICELL),
     +           THETASNEW(IMINNEW:IMAXNEW,LM1,ICELL),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
         ENDDO
        ENDDO
       end if ! present(NTCELL) 



       ENDDO ! I1
       
       if( present(NTCELL) ) deallocate(thetasin)
       
       END SUBROUTINE CREATE_NEWMESH

c-----------------------------------------------
       SUBROUTINE CHEBMESH(NPAN,NCHEB,RI,RO)
       IMPLICIT NONE
       INTEGER NPAN,NCHEB
       DOUBLE PRECISION RI(0:NPAN)
       DOUBLE PRECISION RO(NPAN*(NCHEB+1))
       INTEGER I,K,IK
       DOUBLE PRECISION TAU,PI
       PI=4d0*DATAN(1d0)
       DO I=1,NPAN
        DO K=0,NCHEB
         IK=I*NCHEB+I-K
         TAU=DCOS(((2*K+1)*PI)/(2*(NCHEB+1)))
         TAU=0.5d0*((RI(I)-RI(I-1))*TAU+RI(I)+RI(I-1))
         RO(IK)=TAU
        ENDDO
       ENDDO
       END SUBROUTINE CHEBMESH
       
       end module mod_create_newmesh
