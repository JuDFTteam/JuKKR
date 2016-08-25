       SUBROUTINE CREATE_NEWMESH(INS,NSPIN,R,LMMAXD,IRMIND,IRWS,IRMD,
     +                     LMPOTD,IPAN,IRCUT,IRMIN,VINS,VM2Z,R_LOG,
     +                     NPAN_LOG,NPAN_EQ,NPAN_INST,NCHEB,NPAN_TOT,
     +                     RNEW,RPAN_INTERVALL,IPAN_INTERVALL,VINSNEW)
       IMPLICIT NONE
       INTEGER NSPIN,INS,IRMIN,IPAN,IRWS,IRMD
       INTEGER LMMAXD,LMPOTD,IRMIND
       INTEGER NPAN_LOG,NPAN_EQ,NPAN_INST,NCHEB,NPAN_TOT
       INTEGER IRCUT(0:IPAN)
       DOUBLE PRECISION R(IRWS)
       DOUBLE PRECISION FAC
       PARAMETER (FAC=2d0)
       DOUBLE PRECISION 
     +   VINS(IRMIND:IRMD,LMPOTD,NSPIN),
     +   VM2Z(IRMD,NSPIN),
     +   VINSIN(IRWS,LMPOTD,NSPIN)
       INTEGER I1,IR,ISPIN,IR2,IP,
     +         ISHIFT,ILOGPANSHIFT,ILINPANSHIFT,NPAN_LOGTEMP,
     +         IMIN,IMAX,IMINNEW,IMAXNEW,LM1,
     +         NPAN_LOGNEW,NPAN_EQNEW
       DOUBLE PRECISION R_LOG,RMIN,RMAX,RVAL
       DOUBLE PRECISION RNEW(NPAN_TOT*(NCHEB+1)),
     +                  RPAN_INTERVALL(0:NPAN_TOT)
       INTEGER IPAN_INTERVALL(0:NPAN_TOT)
       DOUBLE PRECISION VINSNEW(NPAN_TOT*(NCHEB+1),LMPOTD,NSPIN)
       NPAN_LOGTEMP=0
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
       DO IP=0,NPAN_LOG-ILOGPANSHIFT
        RVAL=(FAC**IP-1d0)/(FAC**(NPAN_LOG-ILOGPANSHIFT)-1d0)       
        RPAN_INTERVALL(IP+ISHIFT)= RMIN+RVAL*(RMAX-RMIN)
        IPAN_INTERVALL(IP+ISHIFT)= (IP+ISHIFT)*(NCHEB+1)
        IF (ISHIFT.EQ.0.AND.RPAN_INTERVALL(IP).GT.R(IRMIN)) THEN
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
         RPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP)=R(IRCUT(IP+1))
         IPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP)=(NPAN_LOG+NPAN_EQ+IP)*
     +           (NCHEB+1)
       ENDDO ! NPAN_INST

       NPAN_EQNEW=NPAN_EQ+NPAN_LOG-NPAN_LOGTEMP
       NPAN_LOGNEW=NPAN_LOGTEMP
       CALL CHEBMESH(NPAN_TOT,NCHEB,RPAN_INTERVALL,RNEW)

c interpolate potential to new mesh
c save input potential to VINSIN
       VINSIN=0d0
       DO ISPIN=1,NSPIN
       DO IR=1,IRWS
        IF (IR.LT.IRMIN) THEN
         VINSIN(IR,1,ISPIN)=VM2Z(IR,ISPIN)
        ELSE  
         DO LM1=1,LMPOTD
          IF (LM1.EQ.1) THEN 
           VINSIN(IR,LM1,ISPIN)=VM2Z(IR,ISPIN)
          ELSE
           VINSIN(IR,LM1,ISPIN)=VINS(IR,LM1,ISPIN)
          ENDIF
         ENDDO
        ENDIF
       ENDDO
       ENDDO
       VINSNEW=0d0
       DO ISPIN=1,NSPIN
        DO LM1=1,LMPOTD
         IMIN=1
         IMAX=IRMIN
         DO IP=1,NPAN_LOGNEW
          IMINNEW=IPAN_INTERVALL(IP-1)+1
          IMAXNEW=IPAN_INTERVALL(IP)
c           IF (LM1.EQ.1) THEN
c            CALL INTERPOLSPLINE(R(IMIN:IMAX),RNEW(IMINNEW:IMAXNEW),
c     +        VM2Z(IMIN:IMAX,ISPIN),VINSNEW(IMINNEW:IMAXNEW,LM1,ISPIN),
c     +        IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
c           ELSE
             CALL INTERPOLSPLINE(R(IMIN:IMAX),RNEW(IMINNEW:IMAXNEW),
     +           VINSIN(IMIN:IMAX,LM1,ISPIN),
     +           VINSNEW(IMINNEW:IMAXNEW,LM1,ISPIN),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
c           ENDIF
         ENDDO

         IMIN=IRMIN
         IMAX=IRCUT(1)
          DO IP=NPAN_LOGNEW+1,NPAN_LOGNEW+NPAN_EQNEW
           IMINNEW=IPAN_INTERVALL(IP-1)+1
           IMAXNEW=IPAN_INTERVALL(IP)
c            IF (LM1.EQ.1) THEN
c             CALL INTERPOLSPLINE(R(IMIN:IMAX),RNEW(IMINNEW:IMAXNEW),
c     +           VM2Z(IMIN:IMAX,ISPIN),
c     +           VINSNEW(IMINNEW:IMAXNEW,LM1,ISPIN),
c     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
c            ELSE
             CALL INTERPOLSPLINE(R(IMIN:IMAX),RNEW(IMINNEW:IMAXNEW),
     +           VINSIN(IMIN:IMAX,LM1,ISPIN),
     +           VINSNEW(IMINNEW:IMAXNEW,LM1,ISPIN),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
c            ENDIF
          ENDDO
         IR2=0
         DO IP=NPAN_LOGNEW+NPAN_EQNEW+1,NPAN_TOT
          IR2=IR2+1
          IMIN=IRCUT(IR2)+1
          IMAX=IRCUT(IR2+1)
          IMINNEW=IPAN_INTERVALL(IP-1)+1
          IMAXNEW=IPAN_INTERVALL(IP)
c          IF (LM1.EQ.1) THEN
c           CALL INTERPOLSPLINE(R(IMIN:IMAX),RNEW(IMINNEW:IMAXNEW),
c     +           VM2Z(IMIN:IMAX,ISPIN),
c     +           VINSNEW(IMINNEW:IMAXNEW,LM1,ISPIN),
c     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
c          ELSE
           CALL INTERPOLSPLINE(R(IMIN:IMAX),RNEW(IMINNEW:IMAXNEW),
     +           VINSIN(IMIN:IMAX,LM1,ISPIN),
     +           VINSNEW(IMINNEW:IMAXNEW,LM1,ISPIN),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
c          ENDIF
         ENDDO
        ENDDO ! LM1
       ENDDO ! ISPIN
       END 

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
       
