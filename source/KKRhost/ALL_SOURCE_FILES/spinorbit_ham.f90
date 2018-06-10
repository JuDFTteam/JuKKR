SUBROUTINE spinorbit_ham(lmax,lmmaxd,vins,rnew,e,z,c,socscale,  &
        nspin,lmpotd,theta,phi,  &
        ipan_intervall,rpan_intervall,  &
        npan_tot,ncheb,irmdnew,nrmaxd,vnspll,vnspll1,  &
        mode)
IMPLICIT NONE 

INTEGER LMAX,LMMAXD,NSPIN,NPAN_TOT,NCHEB,IRMDNEW,NRMAXD
INTEGER LMPOTD
DOUBLE PRECISION C,Z
DOUBLE COMPLEX   E
DOUBLE PRECISION SOCSCALE 
DOUBLE PRECISION VINS(IRMDNEW,LMPOTD,NSPIN), &
                 RNEW(NRMAXD), &
                 RPAN_INTERVALL(0:NPAN_TOT)
DOUBLE COMPLEX   VNSPLL(2*LMMAXD,2*LMMAXD,IRMDNEW)
DOUBLE COMPLEX   VNSPLL1(2*LMMAXD,2*LMMAXD,IRMDNEW)
INTEGER IPAN_INTERVALL(0:NPAN_TOT)
DOUBLE PRECISION VR(IRMDNEW),DVDR(IRMDNEW)
DOUBLE PRECISION RMASS(IRMDNEW),HSOFAC(IRMDNEW)
DOUBLE PRECISION RNUCL,ATN,WIDTHFAC,PHI,THETA
INTEGER  IR,IP,LM1,LM2,ISPIN,IRMIN,IRMAX,NCOLL
DOUBLE COMPLEX LSMH(2*LMMAXD,2*LMMAXD),TEMP
DOUBLE PRECISION CLAMBDACINV(0:NCHEB,0:NCHEB)
CHARACTER(LEN=*) :: MODE
LOGICAL TEST,OPT
EXTERNAL TEST,OPT

vnspll1=(0D0,0D0)
vr=0D0
DO ispin=1,nspin
  DO ir=1,ipan_intervall(npan_tot)
    vr(ir)=vr(ir)+vins(ir,1,ispin)/nspin
  END DO
END DO
! derivative of potential
dvdr=0D0
CALL getclambdacinv(ncheb,clambdacinv)
DO ip=1,npan_tot
  irmin=ipan_intervall(ip-1)+1
  irmax=ipan_intervall(ip)
  widthfac= 2D0/(rpan_intervall(ip)-rpan_intervall(ip-1))
  CALL dgemv('N',ncheb+1,ncheb+1,1D0,clambdacinv,ncheb+1,  &
      vr(irmin:irmax),1,0D0,dvdr(irmin:irmax),1)
  dvdr(irmin:irmax)= dvdr(irmin:irmax)*widthfac
END DO
! core potential
IF (z > 24D0) THEN
  atn=-16.1532921+2.70335346*z
ELSE
  atn=0.03467714+2.04820786*z
END IF
rnucl=1.2D0/0.529177D0*atn**(1./3D0)*1.d-5

DO ir=1,ipan_intervall(npan_tot)
  IF (rnew(ir) <= rnucl) THEN
!        DVDR(IR)=DVDR(IR)+2d0*Z*RNEW(IR)/RNUCL**3d0
  ELSE
!        DVDR(IR)=DVDR(IR)+2d0*Z/RNEW(IR)**2d0
  END IF
  dvdr(ir)=dvdr(ir)+2D0*z/rnew(ir)**2D0
END DO
! contruct LS matrix

CALL spin_orbit_compl(lmax,lmmaxd,lsmh)

! roate LS matrix
ncoll=1
IF (ncoll == 1) THEN
  CALL rotatematrix(lsmh,theta,phi,lmmaxd,1)
END IF

IF (mode == 'transpose') THEN
  DO lm1=1,2*lmmaxd
    DO lm2=1,lm1-1
      temp=lsmh(lm2,lm1)
      lsmh(lm2,lm1)=lsmh(lm1,lm2)
      lsmh(lm1,lm2)=temp
    END DO
  END DO
ELSE IF (mode == '1') THEN
END IF
! contruct prefactor of spin-orbit hamiltonian

hsofac=0D0
IF (test('NOSOC   ').OR.z < 1D-6) THEN
  DO ir=1,irmdnew
    DO lm1=1,2*lmmaxd
      DO lm2=1,2*lmmaxd
        vnspll1(lm1,lm2,ir)=vnspll(lm1,lm2,ir)
      END DO
    END DO
  END DO
ELSE
  DO ir=1,irmdnew
    rmass(ir)=0.5D0-0.5D0/c**2*((vr(ir)-REAL(e))-2D0*z/rnew(ir))
    hsofac(ir)=socscale/(2D0*rmass(ir)**2*c**2*rnew(ir))*dvdr(ir)
    
! add to potential
    DO lm1=1,2*lmmaxd
      DO lm2=1,2*lmmaxd
        vnspll1(lm1,lm2,ir)=vnspll(lm1,lm2,ir)+hsofac(ir)*lsmh(lm1,lm2)
      END DO
    END DO
  END DO
END IF
END SUBROUTINE spinorbit_ham


