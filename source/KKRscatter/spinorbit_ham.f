      SUBROUTINE SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,E,Z,C,NSRA,NSPIN,
     +                  LMPOTD,IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,
     +                  NCHEB,IRMDNEW,VNSPLL,THETA,PHI,MODE)

      IMPLICIT NONE 

      INTEGER LMAX,LMMAXD,NSPIN,NSRA,NPAN_TOT,NCHEB,IRMDNEW,LMPOTD
      DOUBLE PRECISION C,Z
      DOUBLE COMPLEX   E
      DOUBLE PRECISION VINS(IRMDNEW,LMPOTD,NSPIN),
     +                 RNEW(IRMDNEW),
     +                 RPAN_INTERVALL(0:NPAN_TOT)
      DOUBLE COMPLEX   VNSPLL(2*LMMAXD,2*LMMAXD,IRMDNEW)
      INTEGER IPAN_INTERVALL(0:NPAN_TOT)
      DOUBLE PRECISION VR(IRMDNEW),DVDR(IRMDNEW),RMASS(IRMDNEw)
      DOUBLE PRECISION RNUCL,ATN,WIDTHFAC,PHI,THETA
      INTEGER  IR,IP,LM1,LM2,ISPIN,IRMIN,IRMAX,NCOLL,I1
      DOUBLE COMPLEX LSMH(2*LMMAXD,2*LMMAXD),TEMP
c      DOUBLE COMPLEX HSOFAC(IRMDNEW)
      DOUBLE PRECISION HSOFAC(IRMDNEW)
      DOUBLE PRECISION CLAMBDACINV(0:NCHEB,0:NCHEB)
      DOUBLE PRECISION MATVEC_DMDM
      CHARACTER(LEN=*) :: MODE
      LOGICAL TEST,OPT
      EXTERNAL TEST,OPT

      VR=0d0
      DO ISPIN=1,NSPIN
        DO IR=1,IPAN_INTERVALL(NPAN_TOT)
          VR(IR)=VR(IR)+VINS(IR,1,ISPIN)/NSPIN
        END DO
      END DO
     
c derivative of potential
      DVDR=0d0
      CALL GETCLAMBDACINV(NCHEB,CLAMBDACINV)
      DO IP=1,NPAN_TOT
       IRMIN=IPAN_INTERVALL(IP-1)+1
       IRMAX=IPAN_INTERVALL(IP)
       WIDTHFAC= 2d0/(RPAN_INTERVALL(IP)-RPAN_INTERVALL(IP-1))
       CALL DGEMV('N',NCHEB+1,NCHEB+1,1d0,CLAMBDACINV,NCHEB+1,
     +   VR(IRMIN:IRMAX),1,0d0,DVDR(IRMIN:IRMAX),1)
c       DVDR(IRMIN:IRMAX)= MATVEC_DMDM(CLAMBDACINV,VR(IRMIN:IRMAX))
       DVDR(IRMIN:IRMAX)= DVDR(IRMIN:IRMAX)*WIDTHFAC
      ENDDO
c core potential
      IF (Z.GT.24D0) THEN
        ATN=-16.1532921+2.70335346*Z
      ELSE
        ATN=0.03467714+2.04820786*Z
      ENDIF       
      RNUCL=1.2d0/0.529177d0*ATN**(1./3d0)*1.D-5  

      DO IR=1,IPAN_INTERVALL(NPAN_TOT)
       IF (RNEW(IR) .LE. RNUCL) THEN
c        DVDR(IR)=DVDR(IR)+2d0*Z*RNEW(IR)/RNUCL**3d0
       ELSE
c        DVDR(IR)=DVDR(IR)+2d0*Z/RNEW(IR)**2d0
       END IF
        DVDR(IR)=DVDR(IR)+2d0*Z/RNEW(IR)**2d0
      END DO
c contruct LS matrix

      CALL SPIN_ORBIT_COMPL(LMAX,LMMAXD,LSMH)

c rotate LS matrix 
      CALL ROTATEMATRIX(LSMH,THETA,PHI,LMMAXD,1)

      IF (MODE.EQ.'transpose') THEN
       DO LM1=1,2*LMMAXD
        DO LM2=1,LM1-1
         TEMP=LSMH(LM2,LM1)
         LSMH(LM2,LM1)=LSMH(LM1,LM2)
         LSMH(LM1,LM2)=TEMP
        ENDDO
       ENDDO
      ELSEIF (MODE.EQ.'1') THEN
      ENDIF

c contruct prefactor of spin-orbit hamiltonian

      HSOFAC=0d0
      DO IR=1,IPAN_INTERVALL(NPAN_TOT)
       RMASS(IR)=0.5d0-0.5d0/C**2*((VR(IR)-REAL(E))-2d0*Z/RNEW(IR))
       IF (TEST('NOSOC   ').OR.Z.LT.1D-6) THEN 
        HSOFAC(IR)=0d0
       ELSE
        HSOFAC(IR)=1d0/(2d0*RMASS(IR)**2*C**2*RNEW(IR))*DVDR(IR)
       ENDIF
      
c add to potential
       
       DO LM1=1,2*LMMAXD
        DO LM2=1,2*LMMAXD
         VNSPLL(LM1,LM2,IR)=VNSPLL(LM1,LM2,IR)+HSOFAC(IR)*LSMH(LM1,LM2)
        ENDDO
       ENDDO
      END DO

      END SUBROUTINE SPINORBIT_HAM


