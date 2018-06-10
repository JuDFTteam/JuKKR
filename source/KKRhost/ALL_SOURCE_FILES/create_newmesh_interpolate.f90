SUBROUTINE create_newmesh(nspin,r,irmin,irws,ipan,ircut,  &
        vins,visp,r_log,npan_log,npan_eq,ncheb,  &
        npan_tot,rnew,rpan_intervall,ipan_intervall,  &
        vinsnew,ntcell,thetas,thetasnew)

IMPLICIT NONE
INCLUDE 'inc.p'
INTEGER NSPIN,IRMIN(NATYPD),IPAN(NATYPD),IRWS(NATYPD)
INTEGER LMMAXD
PARAMETER (LMMAXD= (LMAXD+1)**2)
INTEGER LMPOTD
PARAMETER (LMPOTD= (LPOTD+1)**2)
INTEGER IRMIND
PARAMETER (IRMIND= IRMD-IRNSD)
INTEGER NPAN_LOG,NPAN_EQ,NCHEB,NPAN_INST,NPAN_TOT(NATYPD)
INTEGER IRCUT(0:IPAND,NATYPD)
DOUBLE PRECISION R(IRMD,NATYPD)
DOUBLE PRECISION FAC
PARAMETER (FAC=2d0)
DOUBLE PRECISION &
  VINS(IRMIND:IRMD,LMPOTD,NSPOTD), &
  VISP(IRMD,NSPOTD), &
  VINSIN(IRMD,LMPOTD,NSPIN)
DOUBLE PRECISION THETAS(IRID,NFUND,NCELLD), &
  THETASIN(IRID,NFUND,NCELLD), &
  THETASNEW(NTOTD*(NCHEBD+1),NFUND,NCELLD)
INTEGER NTCELL(NATYPD)
INTEGER I1,IPOT,IPOTM,IR,ISPIN,IR2,IP,ICELL, &
        ISHIFT,ILOGPANSHIFT,ILINPANSHIFT,NPAN_LOGTEMP, &
        IMIN,IMAX,IMINNEW,IMAXNEW,LM1
DOUBLE PRECISION R_LOG,RMIN,RMAX,RVAL
DOUBLE PRECISION RNEW(NTOTD*(NCHEBD+1),NATYPD), &
                 RPAN_INTERVALL(0:NTOTD,NATYPD)
INTEGER IPAN_INTERVALL(0:NTOTD,NATYPD)
DOUBLE PRECISION VINSNEW(NTOTD*(NCHEBD+1),LMPOTD,NSPOTD)

vinsnew=0D0
thetasnew=0D0
ipotm=0

DO i1 = 1,natypd
  
  ipot=nspin*(i1-1)+1
  npan_inst= ipan(i1)-1
  npan_tot(i1)= npan_log+npan_eq+npan_inst
  
  
! log panel
  rmin=r(2,i1)
  rmax=r_log
  rval=0D0
  ishift=0
  IF (r_log > r(irmin(i1),i1)) THEN
    ilogpanshift=1
    ilinpanshift=0
  ELSE
    ilogpanshift=0
    ilinpanshift=1
  END IF
  
  IF (ilinpanshift == 1) THEN
    WRITE(*,*) 'ERORR: non-spherical part of the potential needs'
    WRITE(*,*) 'to be inside the log panel'
    WRITE(*,*) 'atom (I1):', i1
    WRITE(*,*) 'R_LOG', r_log
    WRITE(*,*) 'R(IRMIN(I1), I1)', r(irmin(i1),i1)
    WRITE(*,*) 'IRMIN(I1)', irmin(i1)
    STOP 'Error creating newmesh'
  END IF
  
  DO ip=0,npan_log-ilogpanshift
    rval=(fac**ip-1D0)/(fac**(npan_log-ilogpanshift)-1D0)
    rpan_intervall(ip+ishift,i1)= rmin+rval*(rmax-rmin)
    ipan_intervall(ip+ishift,i1)= (ip+ishift)*(ncheb+1)
    IF (ishift == 0.AND. rpan_intervall(ip,i1) > r(irmin(i1),i1)) THEN
      ishift=1
      npan_logtemp=ip
      rpan_intervall(ip+1,i1)=rpan_intervall(ip,i1)
      ipan_intervall(ip+1,i1)=(ip+ishift)*(ncheb+1)
      rpan_intervall(ip,i1)=r(irmin(i1),i1)
      ipan_intervall(ip,i1)=ip*(ncheb+1)
    END IF
  END DO ! NPAN_LOG
  
! equivalent panel
  ishift=0
  rmin=r_log
  rmax=r(ircut(1,i1),i1)
  DO ip=0,npan_eq-ilinpanshift
    rpan_intervall(ip+ishift+npan_log,i1)=rmin+ip*(rmax-rmin)/  &
        (npan_eq-ilinpanshift)
    ipan_intervall(ip+ishift+npan_log,i1)=(npan_log+ip+ishift)* (ncheb+1)
  END DO ! NPAN_EQ
  
! intersection zone
  DO ip=1,npan_inst
    rpan_intervall(npan_log+npan_eq+ip,i1)=r(ircut(ip+1,i1),i1)
    ipan_intervall(npan_log+npan_eq+ip,i1)=(npan_log+npan_eq+ip)* (ncheb+1)
  END DO ! NPAN_INST
  
  npan_eq=npan_eq+npan_log-npan_logtemp
  npan_log=npan_logtemp
  
  CALL chebmesh(npan_tot,ncheb,rpan_intervall(0:,i1),rnew(1,i1))
  
! interpolate potential to new mesh
! save input potential to VINSIN
  
  vinsin=0D0
  DO ir=1,irws(i1)
    IF (ir < irmin(i1)) THEN
      vinsin(ir,1,1)=visp(ir,ipot)
      vinsin(ir,1,nspin)=visp(ir,ipot+nspin-1)
    ELSE
      DO lm1=1,lmpotd
        IF (lm1 == 1) THEN
          vinsin(ir,lm1,1)=visp(ir,ipot)
          vinsin(ir,lm1,nspin)=visp(ir,ipot+nspin-1)
        ELSE
          vinsin(ir,lm1,1)=vins(ir,lm1,ipot)
          vinsin(ir,lm1,nspin)=vins(ir,lm1,ipot+nspin-1)
        END IF
      END DO
    END IF
  END DO
  
  DO ispin=1,nspin
    
    ipotm=ipotm+1
    
    DO lm1=1,lmpotd
      imin=1
      imax=irmin(i1)
      DO ip=1,npan_log
        iminnew=ipan_intervall(ip-1,i1)+1
        imaxnew=ipan_intervall(ip,i1)
        CALL interpolspline(r(imin:imax,i1),rnew(iminnew:imaxnew,i1),  &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm),  &
            imax-imin+1,imaxnew-iminnew+1)
      END DO
      
      imin=irmin(i1)
      imax=ircut(1,i1)
      DO ip=npan_log+1,npan_log+npan_eq
        iminnew=ipan_intervall(ip-1,i1)+1
        imaxnew=ipan_intervall(ip,i1)
        CALL interpolspline(r(imin:imax,i1),rnew(iminnew:imaxnew,i1),  &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm),  &
            imax-imin+1,imaxnew-iminnew+1)
      END DO
      ir2=0
      DO ip=npan_log+npan_eq+1,npan_tot(i1)
        ir2=ir2+1
        imin=ircut(ir2,i1)+1
        imax=ircut(ir2+1,i1)
        iminnew=ipan_intervall(ip-1,i1)+1
        imaxnew=ipan_intervall(ip,i1)
        CALL interpolspline(r(imin:imax,i1),rnew(iminnew:imaxnew,i1),  &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm),  &
            imax-imin+1,imaxnew-iminnew+1)
      END DO
    END DO ! LM1
    
  END DO ! ISPIN
! interpolate shape function THETAS to new shape function THETASNEW
! save THETAS to THETASIN
  icell = ntcell(i1)
  DO lm1=1,nfund
    thetasin(:,lm1,icell)=thetas(:,lm1,icell)
    ir2=0
    DO ip=npan_log+npan_eq+1,npan_tot(i1)
      ir2=ir2+1
      imin=ircut(ir2,i1)+1
      imax=ircut(ir2+1,i1)
      iminnew=ipan_intervall(ip-1,i1)+1
      imaxnew=ipan_intervall(ip,i1)
      CALL interpolspline(r(imin:imax,i1),rnew(iminnew:imaxnew,i1),  &
          thetasin(imin-ircut(1,i1):imax-ircut(1,i1),lm1,icell),  &
          thetasnew(iminnew:imaxnew,lm1,icell), imax-imin+1,imaxnew-iminnew+1)
    END DO
  END DO
END DO ! I1
END SUBROUTINE create_newmesh

!-----------------------------------------------

SUBROUTINE chebmesh(npan,ncheb,ri,ro)

INTEGER, INTENT(IN)                      :: npan
INTEGER, INTENT(IN)                      :: ncheb
DOUBLE PRECISION, INTENT(IN)             :: ri(0:npan)
DOUBLE PRECISION, INTENT(OUT)            :: ro(npan*(ncheb+1))
IMPLICIT NONE



INTEGER :: i,k,ik
DOUBLE PRECISION :: tau,pi

pi=4D0*DATAN(1D0)
DO i=1,npan
  DO k=0,ncheb
    ik=i*ncheb+i-k
    tau=DCOS(((2*k+1)*pi)/(2*(ncheb+1)))
    tau=0.5D0*((ri(i)-ri(i-1))*tau+ri(i)+ri(i-1))
    ro(ik)=tau
  END DO
END DO
END SUBROUTINE chebmesh

