      SUBROUTINE LIFETIME_2D_SO(NKPOID,LMAX,VOLBZ,SCAT_MAT,
     +         KF_IRR,NSPO,NSPOH,NSPIN,ISPIN,NCL,IMPLAYER,
     +         NAEZ,WEIGHT,TAU_K_0,RN,ATOMIMP,DOSW,FERMI_VL)
      

      IMPLICIT NONE

      INTEGER,INTENT(IN)          ::  NKPOID,NAEZ,   
     +                                LMAX,NSPO,NSPIN,
     +                                NCL,IMPLAYER,
     +                                ISPIN,NSPOH

      DOUBLE PRECISION,INTENT(IN) ::  VOLBZ,
     +                                KF_IRR(3,NKPOID),
     +                                TAU_K_0(NKPOID,2,2),
     +                                WEIGHT(NKPOID),
     +                                FERMI_VL(NKPOID)


      DOUBLE PRECISION,INTENT(OUT)::  SCAT_MAT(NKPOID,2,NKPOID,2)
      DOUBLE PRECISION            ::  BZKP(3,NKPOID),RN(3,NCL),DOSW
      INTEGER                     ::  ATOMIMP(NCL)

c  local variables

      INTEGER                     ::  NF,LMMAX,NK,ENTG(NKPOID),LMMAXSO,
     +                                LMCLSO,LM1,LM2,NC,J,
     +                                INDEX_LAYER(NCL)


      COMPLEX,ALLOCATABLE         ::  GLL(:,:),COEFF_0(:,:,:,:),
     +                                DTMATLM(:,:),
     +                                DELTAMAT(:,:)
      INTEGER                     ::  NTHETA,NPHI,ITHETA,IPHI

      LOGICAL                     ::  OPT,TEST
      EXTERNAL                    ::  OPT,TEST

      WRITE(6,*) "In lifetime_2D_SO"

      LMMAX=(LMAX+1)**2
      WRITE(6,*) "NKPOID",NKPOID
      LMMAXSO=NSPOH*LMMAX
c      LMCLSO=NSPOH*LMMAX*NCL
      LMCLSO=NSPO*LMMAX*NCL

      write(6,*) "read in the Green function of the perturbed" 
      write(6,*) "system at the impurity side"
      ALLOCATE(GLL(LMCLSO,LMCLSO))
      CALL READ_GREEN_LL(GLL,LMCLSO)

      write(6,*) "read in the delta t matrix"
      ALLOCATE(DTMATLM(LMCLSO,LMCLSO))
      ALLOCATE(DELTAMAT(LMCLSO,LMCLSO))
      CALL READ_DELTA_TMAT(DTMATLM,DELTAMAT,LMCLSO,NCL)

      write(6,*) "read in the wave function coefficients of the host"
      ALLOCATE(COEFF_0(LMMAXSO,NAEZ,4,NKPOID))
      CALL READ_COEFF_2D_SO(COEFF_0,BZKP,LMMAXSO,NKPOID,
     +                      NAEZ,ENTG,NSPOH)


      write(6,*) "before calculation of the scattering matrix"
        IF (NSPOH.EQ.2.OR.NSPO.EQ.2) THEN
        CALL SCATTERING_MATRIX_2D(SCAT_MAT,NKPOID,
     +      NAEZ,LMAX,LMMAX,DTMATLM,DELTAMAT,GLL,
     +      COEFF_0,BZKP,NCL,IMPLAYER,INDEX_LAYER,ISPIN,NSPO,
     +      NSPOH,ENTG,VOLBZ,WEIGHT,TAU_K_0,RN,ATOMIMP,DOSW,FERMI_VL)
      IF (TEST('COEFF   ')) THEN   
      OPEN(unit=81,file="Aniso_spinconserve",FORM="FORMATTED")
      OPEN(unit=82,file="Aniso_spinrelax",FORM="FORMATTED")
      OPEN(unit=83,file="Aniso_momentum",FORM="FORMATTED")
      NPHI=10
      NTHETA=10
      DO ITHETA=1,NTHETA
       DO IPHI=1,NPHI
      CALL READ_COEFF_2D_SO_ANI(COEFF_0,BZKP,LMMAXSO,NKPOID,
     +                      NAEZ,ENTG,NSPOH,ITHETA,IPHI)


      write(6,*) "before calculation of the scattering matrix"

        CALL SCATTERING_MATRIX_2D_ANI(SCAT_MAT,NKPOID,
     +      NAEZ,LMAX,LMMAX,DTMATLM,DELTAMAT,GLL,
     +      COEFF_0,BZKP,NCL,IMPLAYER,INDEX_LAYER,ISPIN,NSPO,
     +      NSPOH,ENTG,VOLBZ,WEIGHT,TAU_K_0,RN,ATOMIMP,
     +      ITHETA,IPHI,NTHETA,NPHI,DOSW,FERMI_VL)
       ENDDO
      ENDDO
      CLOSE(81)
      CLOSE(82)
      ENDIF
      ELSE
      WRITE(6,*) 'relaxation in NOSOC'
        CALL SCATTERING_MATRIX_2D_NOSOC(SCAT_MAT,NKPOID,
     +      NAEZ,LMAX,LMMAX,DTMATLM,DELTAMAT,GLL,
     +      COEFF_0,BZKP,NCL,IMPLAYER,INDEX_LAYER,ISPIN,NSPO,
     +      NSPOH,ENTG,VOLBZ,WEIGHT,TAU_K_0,RN,ATOMIMP,FERMI_VL)
      ENDIF

      write(6,*) "after calculation of the scattering matrix"


 5678 FORMAT(9E17.9)

      DEALLOCATE(GLL)
      DEALLOCATE(DTMATLM)
      DEALLOCATE(DELTAMAT)
      DEALLOCATE(COEFF_0)

      WRITE(6,*) "END OF LIFETIME 2D"


      END SUBROUTINE LIFETIME_2D_SO
