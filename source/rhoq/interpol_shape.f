       SUBROUTINE INTERPOL_SHAPE(R,IRMIN,IRWS,IPAN,IRCUT,
     +                     R_LOG,NPAN_LOG,NPAN_EQ,NCHEB,
     +                     NPAN_TOT,RNEW,RPAN_INTERVALL,IPAN_INTERVALL,
     +                     THETAS,THETASNEW,NFUN,NPAN_LOGNEW,NPAN_EQNEW)
       use mod_interpolspline, only: interpolspline
       IMPLICIT NONE
       INTEGER IRMIN,IPAN,IRWS,NFUN
       INTEGER NPAN_LOG,NPAN_EQ,NCHEB,NPAN_INST,NPAN_TOT
       INTEGER NPAN_LOGNEW,NPAN_EQNEW
       INTEGER IRCUT(1:IPAN)
       DOUBLE PRECISION R(IRWS)
       DOUBLE PRECISION FAC
       PARAMETER (FAC=2d0)
       DOUBLE PRECISION THETAS(IRWS-IRMIN,NFUN),
     +   THETASNEW((IPAN+30)*(NCHEB+1),NFUN)
       double precision, allocatable :: THETASIN(:,:)
       INTEGER IPOT,IPOTM,IR2,IP,
     +         ISHIFT,ILOGPANSHIFT,ILINPANSHIFT,NPAN_LOGTEMP,
     +         IMIN,IMAX,IMINNEW,IMAXNEW,LM1
       DOUBLE PRECISION R_LOG,RMIN,RMAX,RVAL
       DOUBLE PRECISION RNEW((IPAN+30)*(NCHEB+1)),
     +                  RPAN_INTERVALL(0:(IPAN+30))
       INTEGER IPAN_INTERVALL(0:(IPAN+30))


! allocations
       allocate(THETASIN(IRWS-IRMIN,NFUN))

       THETASNEW=0d0
       IPOTM=0
       
       
        IPOT=1
        NPAN_INST= IPAN-1
        NPAN_TOT= NPAN_LOG+NPAN_EQ+NPAN_INST
!         write(*,*) 'npan_inst',npan_inst, npan_tot


c log panel
       RMIN=R(2)
       RMAX=R_LOG
       RVAL=0d0
       ISHIFT=0
       IF (R_LOG.GT.R(IRMIN)) THEN
         ILOGPANSHIFT=1
         ILINPANSHIFT=0
       ELSE
         ILOGPANSHIFT=0
         ILINPANSHIFT=1
       ENDIF 
  
       IF (ILINPANSHIFT.EQ.1) THEN
        STOP 'non-spherical part of the potential needs to be inside 
     +        the log panel'
       ENDIF
      
       DO IP=0,NPAN_LOG-ILOGPANSHIFT
        RVAL=(FAC**IP-1d0)/(FAC**(NPAN_LOG-ILOGPANSHIFT)-1d0)       
        RPAN_INTERVALL(IP+ISHIFT)= RMIN+RVAL*(RMAX-RMIN)
        IPAN_INTERVALL(IP+ISHIFT)= (IP+ISHIFT)*(NCHEB+1)
        IF (ISHIFT.EQ.0.AND.
     +      RPAN_INTERVALL(IP).GT.R(IRMIN)) THEN
         ISHIFT=1
         NPAN_LOGTEMP=IP
         RPAN_INTERVALL(IP+1)=RPAN_INTERVALL(IP)
         IPAN_INTERVALL(IP+1)=(IP+ISHIFT)*(NCHEB+1)
         RPAN_INTERVALL(IP)=R(IRMIN)
         IPAN_INTERVALL(IP)=IP*(NCHEB+1)
        ENDIF
       ENDDO ! NPAN_LOG

c equivalent panel
       ISHIFT=0
       RMIN=R_LOG
       RMAX=R(IRCUT(1))
       DO IP=0,NPAN_EQ-ILINPANSHIFT
        RPAN_INTERVALL(IP+ISHIFT+NPAN_LOG)=RMIN+IP*(RMAX-RMIN)/
     +          (NPAN_EQ-ILINPANSHIFT)
        IPAN_INTERVALL(IP+ISHIFT+NPAN_LOG)=(NPAN_LOG+IP+ISHIFT)*
     +          (NCHEB+1)
       ENDDO ! NPAN_EQ

c intersection zone
       DO IP=1,NPAN_INST
!          write(*,*) ip,ircut(ip),npan_log+npan_eq+ip,ircut(ip+1)
         RPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP)=R(IRCUT(IP+1))
         IPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP)=(NPAN_LOG+NPAN_EQ+IP)*
     +           (NCHEB+1)
       ENDDO ! NPAN_INST

       NPAN_EQNEW=NPAN_EQ+NPAN_LOG-NPAN_LOGTEMP
       NPAN_LOGNEW=NPAN_LOGTEMP

       CALL CHEBMESH2(NPAN_TOT,NCHEB,RPAN_INTERVALL(0:),
     +               RNEW(1))

c interpolate shape function THETAS to new shape function THETASNEW
c save THETAS to THETASIN
        DO LM1=1,NFUN
         THETASIN(:,LM1)=THETAS(:,LM1)
         IR2=0
         DO IP=NPAN_LOGNEW+NPAN_EQNEW+1,NPAN_TOT
          IR2=IR2+1
          IMIN=IRCUT(IR2)+1
          IMAX=IRCUT(IR2+1)
          IMINNEW=IPAN_INTERVALL(IP-1)+1
          IMAXNEW=IPAN_INTERVALL(IP)
!           write(*,*) IMIN,IMAX,IMINNEW,IMAXNEW
          CALL INTERPOLSPLINE(R(IMIN:IMAX),RNEW(IMINNEW:IMAXNEW),
     +           THETASIN(IMIN-IRCUT(1):IMAX-IRCUT(1),LM1),
     +           THETASNEW(IMINNEW:IMAXNEW,LM1),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
         ENDDO
        ENDDO
       
       deallocate(thetasin)
       
       END 

c-----------------------------------------------
       SUBROUTINE CHEBMESH2(NPAN,NCHEB,RI,RO)
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
       END SUBROUTINE CHEBMESH2
       
