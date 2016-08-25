      SUBROUTINE HAM_SO(LMAX,LMMAXD,IRMD,ISPOEND,VIN,RM,IE,E,Z,C,
     +                  NSRA, LSM,HSOFAC,IRMINSO,IPAN,IRCUT,NSPIN,RDVDR)

      implicit none 

c  input parameters: 

      integer          ::   IRMD,LMAX,LMMAXD,IRMINSO,IE,NSPIN,
     +                      ISPOEND,        !maximal r for the SP-O
                                            !  interaction (chosen as MT-Radius)
     +                      IPAN,IRCUT(0:IPAN),NSRA
      real             ::   VIN(IRMD,NSPIN),
     +                      RM(IRMD)
      DOUBLE COMPLEX   ::   E
      DOUBLE PRECISION ::   C,Z

c  output parameters: 
      double precision,intent(out) ::  HSOFAC(IRMINSO:IRMD)
      double complex,intent(out)   ::  LSM(2*LMMAXD,2*LMMAXD)
      double complex               ::  LSMH(2*LMMAXD,2*LMMAXD)
      double precision             ::  RDVDR(IRMINSO:IRMD)

c local variables 
      integer     ::     nr,LM1,LM2,NRNUCL,ISP,FACTOR
      real        ::     RMASS(IRMD),RNUCL,ATN,VR(IRMD)
      real        ::     DVDR(IRMD) ! radial derivative of the spherical input potential

      VR=0d0
      LSM=0d0
      LSMH=0d0

      DO ISP=1,NSPIN
        DO NR=1,IRMD
          VR(NR)=VR(NR)+VIN(NR,ISP)/NSPIN
c          WRITE(6,"(I5,e17.9)") NR,VR(NR)
        END DO
      END DO

      IF (NSPIN==2) THEN

        WRITE(6,*) "calculation of the mean of the spin up and 
     +                                  the spin down potential:" 
        WRITE(6,*) "required for the spin-orbit Hamiltonian"

      END IF

      CALL DER_POT(DVDR,VR,RM,IRMD,IPAN,IRCUT)

      IF (NSRA == 2 ) THEN 
        FACTOR=10d0

c        WRITE(6,*) "z",Z
        IF (Z .EQ. 27.0d0) THEN
          ATN=58.933d0   ! for Co
        ELSE IF (Z .EQ. 3.0d0) THEN
          ATN=6.94d0   ! for Li
          FACTOR=40d0
        ELSE IF (Z .EQ. 11.0d0) THEN
          ATN=22.99d0   ! for Na
          FACTOR=30d0
        ELSE IF (Z .EQ. 19.0d0) THEN
          ATN=39.10d0   ! for K
          FACTOR=20d0
        ELSE IF (Z .EQ. 37.0d0) THEN
          ATN=85.47d0   ! for Rb
        ELSE IF (Z .EQ. 55.0d0) THEN
          ATN=132.91d0   ! for Cs
        ELSE IF (Z .EQ. 26.0d0) THEN
          ATN=55.847d0   ! for Fe
        ELSE IF (Z .EQ. 29.0d0) THEN
          ATN=63.546d0   ! for Cu
        ELSE IF (Z .EQ. 30.0d0) THEN
          ATN=65.380d0   ! for Zn
        ELSE IF (Z .EQ. 42.0d0) THEN
          ATN=95.960d0   ! for Mo
        ELSE IF (Z .EQ. 46.0d0 ) THEN        
          ATN=106.42d0   ! for Pd
        ELSE IF (Z .EQ. 47.0d0) THEN
          ATN=107.87d0   ! for Ag
        ELSE IF (Z .EQ. 48.0d0) THEN
          ATN=112.41d0   ! for Cd
        ELSE IF (Z .EQ. 49.0d0) THEN
          ATN=114.82d0   ! for In
        ELSE IF (Z .EQ. 50.0d0) THEN
          ATN=118.71d0   ! for Sn
        ELSE IF (Z .EQ. 51.0d0 ) THEN        
          ATN=121.76d0   ! for Sb
        ELSE IF (Z .EQ. 52.0d0 ) THEN        
          ATN=127.60d0   ! for Te
        ELSE IF (Z .EQ. 56.0d0 ) THEN        
          ATN=138.91d0   ! for La
        ELSE IF (Z .EQ. 72.0d0 ) THEN        
          ATN=178.49d0   ! for Hf
        ELSE IF (Z .EQ. 73.0d0 ) THEN        
          ATN=180.95d0   ! for Ta
        ELSE IF (Z .EQ. 74.0d0 ) THEN        
          ATN=183.85d0   ! for W
        ELSE IF (Z .EQ. 75.0d0 ) THEN        
          ATN=186.21d0   ! for Re
        ELSE IF (Z .EQ. 76.0d0 ) THEN        
          ATN=190.20d0   ! for Os
        ELSE IF (Z .EQ. 77.0d0 ) THEN        
          ATN=192.22d0   ! for Ir
        ELSE IF (Z .EQ. 78.0d0) THEN
          ATN=195.08d0   ! for Pt
        ELSE IF (Z .EQ. 79.0d0) THEN
          ATN=196.97d0   ! for Au
        ELSE IF (Z .EQ. 82.0d0) THEN
          ATN=207.20d0   ! for Pb
        ELSE IF (Z .EQ. 83.0d0) THEN
          ATN=208.98d0   ! for Bi
        ELSE IF (Z .EQ. 84.0d0) THEN
          ATN=209.00d0   ! for Po
        ELSE
          STOP "CHANGE ATN==atomic number"
        END IF

c ---> estimate the nucleus of the core        
        RNUCL=1.2d0/0.529177d0*ATN**(1./3d0)*1.D-5  
        RNUCL=FACTOR*RNUCL
c        RNUCL=10.0d0*RNUCL
c      RNUCL=3d0*RNUCL
      ELSE 

        RNUCL=0d0

      END IF

      do nr=1,IRMD
        IF (RM(NR) .LE. RNUCL) THEN
          DVDR(NR)=DVDR(NR)+2d0*Z*RM(NR)/RNUCL**3d0
          NRNUCL=NR
        ELSE
c          WRITE(6,"(I5,3e17.9)") NR,RM(NR),DVDR(NR),2d0*Z/RM(NR)**2d0
          DVDR(NR)=DVDR(NR)+2d0*Z/RM(NR)**2d0
c          WRITE(6,"(I5,3e17.9)") NR,RM(NR),DVDR(NR)
        END IF
      end do

      CALL REL_MASS(RMASS,RM,VR,Z,E,C,IRMD)
      
      RMASS(1)=RMASS(2)
      DVDR(1)=DVDR(2)
      DVDR(IRMD)=DVDR(IRMD-1)

      CALL SPIN_ORBIT_COMPL(LMAX,LMMAXD,LSMH)

      DO LM2=1,2*LMMAXD
        DO LM1=1,2*LMMAXD
          LSM(LM2,LM1)=LSMH(LM2,LM1)
c          LSM(LM2,LM1)=LSMH(LM1,LM2)
        END DO
      END DO

      HSOFAC=0d0
c      WRITE(16,*) "HSOFAC"
      DO NR=IRMINSO,ISPOEND
c      DO NR=IRMINSO,ISPOEND-1
        HSOFAC(NR)=1d0/(2d0*RMASS(NR)**2*C**2*RM(NR))*DVDR(NR)
        RDVDR(NR)=1d0/(2d0*C**2*RM(NR))*DVDR(NR)
c        WRITE(16,"((I5),(8e17.9))") NR,VR(NR),HSOFAC(NR),
c     +                              RMASS(NR),DVDR(NR),RM(NR)
      END DO

      IF (IRMINSO .EQ. 1)  HSOFAC(1)=HSOFAC(2)
      IF (IRMINSO .EQ. 1)  RDVDR(1)=RDVDR(2)


      END SUBROUTINE HAM_SO


c************************************************************************************
cc the subroutine REL_MASS calculates the relativistic mass,
cc depending on the radial position.   
c************************************************************************************

      SUBROUTINE REL_MASS(RMASS,RM,VM,Z,E,C,IRMD)

      implicit none 

c  input parameters: 

      integer          ::   IRMD
      real             ::   RM(IRMD),VM(IRMD)
      DOUBLE COMPLEX   ::   E
      DOUBLE PRECISION ::   C,Z

c  output parameters: 
      real        ::     RMASS(IRMD), ! radial derivative of the spherical input potential
     +                   DIFFV
c local variables 
      integer     ::     nr

      do nr=1,IRMD

        DIFFV = (VM(NR)-E) - 2.0D0*Z/RM(NR)
        RMASS(nr)=0.5d0 - 0.5d0/C**2*DIFFV

      end do

      END SUBROUTINE REL_MASS

c************************************************************************************
cc the subroutine DER_POT calculates the radial derivative of the spherical potential 
c************************************************************************************

      SUBROUTINE DER_POT(DVDR,VR,RM,IRMD,IPAN,IRCUT)

      implicit none 

c  input parameters: 

      integer     ::     IRMD,IPAN,IRCUT(0:IPAN)
      real        ::     VR(IRMD),
     +                   RM(IRMD)

c  output parameters: 
      real        ::     DVDR(IRMD) ! radial derivative of the spherical input potential

c local variables 
      integer     ::     nr,ip

c      WRITE(6,*) "IP,IRCUT",0,IRCUT(0)
      DO IP=1,IPAN

c        WRITE(6,*) "IP,IRCUT",IP,IRCUT(IP)
        IF (IP==1) THEN

          dVdr(1)=(VR(2)-VR(1))/(RM(2)-RM(1))
          do nr=2,IRCUT(1)-1
            dVdr(nr)=0.5d0*( 
     +            (VR(nr)-VR(nr-1))/(RM(nr)-RM(nr-1)) +
     +            (VR(nr+1)-VR(nr))/(RM(nr+1)-RM(nr)) )
          end do
          dVdr(IRCUT(1))=0.5d0*(
     +           (VR(IRCUT(1))-VR(IRCUT(1)-1))/
     +               (RM(IRCUT(1))-RM(IRCUT(1)-1))  +
     +           (VR(IRCUT(1)+2)-VR(IRCUT(1)))/
     +               (RM(IRCUT(1)+2)-RM(IRCUT(1)))  )

        ELSE 

          dVdr(IRCUT(IP-1)+1)=dVdr(IRCUT(IP-1)) 
          do nr=IRCUT(IP-1)+2,IRCUT(IP)-1
            dVdr(nr)=0.5d0*( 
     +            (VR(nr)-VR(nr-1))/(RM(nr)-RM(nr-1)) +
     +            (VR(nr+1)-VR(nr))/(RM(nr+1)-RM(nr)) )
          end do
          nr=IRCUT(IP)
c          WRITE(6,*) "NR",NR
          IF (nr==IRMD) THEN
            dVdr(nr)=
     +            (VR(nr)-VR(nr-1))/(RM(nr)-RM(nr-1))      
c            WRITE(116,"(I5,2e17.9)") NR,DVDR(NR)
          ELSE
            dVdr(nr)=0.5d0*( 
     +            (VR(nr)-VR(nr-1))/(RM(nr)-RM(nr-1)) +
     +            (VR(nr+2)-VR(nr))/(RM(nr+2)-RM(nr)) )
c            WRITE(116,"(I5,2e17.9)") NR,DVDR(NR)
          END IF
        END IF

      END DO

c     WRITE(6,*) "IRMD",IRMD
c     WRITE(6,*) "VR(IRMD),VR(IRMD-1)"
c     WRITE(6,"(2e17.9)") VR(IRMD),VR(IRMD-1)
c     WRITE(6,*) "RM(IRMD),RM(IRMD-1)"
c     WRITE(6,"(2e17.9)") RM(IRMD),RM(IRMD-1)
      DVDR(IRMD)=(VR(IRMD)-VR(IRMD-1))/(RM(IRMD)-RM(IRMD-1))      
c      WRITE(6,*) "DVDR(IRMD)"
c      WRITE(6,"(2e17.9)") DVDR(IRMD)
c      DVDR(IRMD)=DVDR(IRMD-1)
c      DVDR(IRMD)=(DVDR(IRMD)-DVDR(IRMD-1))/(RM(IRMD)-RM(IRMD-1))      
c      WRITE(6,*) "DVDR(IRMD)"
c      WRITE(6,"(2e17.9)") DVDR(IRMD)

c     DO NR=1,IRMD
c       WRITE(116,"(I5,2e17.9)") NR,DVDR(NR)
c     END DO

      END SUBROUTINE DER_POT
