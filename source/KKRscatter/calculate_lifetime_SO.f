      SUBROUTINE CALCULATE_LIFETIME(LMAX,NSPO,NSPOH,NSPIN,
     +                     NCL_IMP,ALATC,INTERF,NAEZ,IMPLAYER,
     +                     RR,RBASIS,NR,LATT,NBX,NBY,NBZ,BBYA,CBYA)
      
      IMPLICIT NONE
      include 'inc.p'
      include 'inc.cls'

      INTEGER NKPOID
      PARAMETER (NKPOID=100000)
      LOGICAL,INTENT(IN)           ::  INTERF
      INTEGER,INTENT(IN)           ::  LMAX,NSPO,NSPIN,NSPOH,
     +                                 NCL_IMP,NAEZ,IMPLAYER,NR
      DOUBLE PRECISION             ::  VOLBZ,ALATC,
     +                                 RR(3,0:NRD),
     +                                 RBASIS(3,*) 

c local variables

      DOUBLE PRECISION,ALLOCATABLE ::  KF_IRR(:,:),    
     +                                 TAU_K_0(:,:,:),   
     +                                 SCAT_MAT(:,:,:,:), ! scattering matrix
     +                                 WEIGHTS_2D(:)
      DOUBLE PRECISION              :: RCLSIMP(3,NATOMIMPD),DIFF,
     +                                 RN(3,NATOMIMPD),TMP(3),
     +                                 V_FERMI(3,NKPOID),
     +                                 FERMI_VL(NKPOID),
     +                                 DOSW
      INTEGER                       :: ATOMIMP(NATOMIMPD)
      INTEGER                       :: STATUS_OPEN,STATUS_READ,
     +                                 NKFS,ISPIN,J,ISH,ISIMP,I,
     +                                 NATOMIMP,N,NA,JATOM,
     +                                 LATT,NBX,NBY,NBZ,NOFKS,NK
      DOUBLE PRECISION              :: BBYA,CBYA,BSC,CSC
      LOGICAL                       :: OPT
      EXTERNAL                      :: OPT      
      CHARACTER*40                  :: NAME,NEW

c read in BZ volume
      IF (NBY.EQ.0) NBY = NBX/BBYA
      IF (NBZ.EQ.0) NBZ = NBX/CBYA
  
        NAME='mesh/kp'

        CALL SNAME(NAME,NEW,LATT)
        CALL SNAME(NEW,NAME,NBX)
        CALL SNAME(NAME,NEW,NBY)
        CALL SNAME(NEW,NAME,NBZ)


        OPEN(52,FILE=NAME,status='old')
        READ (52,*) NOFKS,VOLBZ,BSC,CSC
        WRITE(6,*) 'VOLBZ',VOLBZ
     

      IF (INTERF .EQ. .FALSE.) THEN

      ELSE     ! (INTERF .EQ. TRUE)

       ALLOCATE(KF_IRR(3,NKPOID))
       ALLOCATE(WEIGHTS_2D(NKPOID))

       OPEN(unit=76,file="kpoints_FS",status="old",
     +              action="read",iostat=STATUS_OPEN)
       IF (STATUS_OPEN /= 0) STOP "Problems to open kpoints_FS"
       NKFS=0
       ein: DO
        NKFS=NKFS+1
        READ(76,"(4(E17.9))",iostat=STATUS_READ)
     +    (KF_IRR(J,NKFS),J=1,3),WEIGHTS_2D(NKFS) 
        IF (STATUS_READ /= 0 ) EXIT
       IF (NKFS.GE.NKPOID) STOP 'change NKPOID to lager value'
       END DO ein
       NKFS=NKFS-1
       CLOSE(76)
c       OPEN(unit=77,file="FERMI_VELOCITY",status="old",
c     +              action="read")
       OPEN(unit=78,file="vectorFERMIVL",status="old",
     +              action="read")
       DO NK=1,NKFS
c        READ(77,"(4(E17.9))")
c     +    (BZKP(J,NK),J=1,3),FERMI_VL(NK) 
        READ(78,"(3(E17.9))")
     +    (V_FERMI(J,NK),J=1,3) 
       FERMI_VL(NK)=sqrt(V_FERMI(1,NK)**2+V_FERMI(2,NK)**2+
     +                   V_FERMI(3,NK)**2)
       END DO 
c       CLOSE(77)
       CLOSE(78)
c      calculate density of states
       DOSW=0d0
       DO NK=1,NKFS
         DOSW=DOSW+WEIGHTS_2D(NK)/FERMI_VL(NK)/VOLBZ
       ENDDO
       WRITE(6,*) 'Density of states',DOSW
c      read in cluster from scoef for calculate impurity RR
       REWIND(25)
       READ(25,FMT=*) NATOMIMP
       DO I=1,NATOMIMP
        READ(25,FMT=*) (RCLSIMP(J,I),J=1,3),ATOMIMP(I)
        write(6,FMT=*) (RCLSIMP(J,I),J=1,3),ATOMIMP(I)
       ENDDO
       DO JATOM=1,NAEZ
        IF (JATOM.EQ.ATOMIMP(1)) THEN
        DO NA=1,NAEZ
         DO N=0,NR
          DO J=1,3
            TMP(J)=RR(J,N)+RBASIS(J,NA)-RBASIS(J,JATOM)
          ENDDO
          DO I=1,NATOMIMP
           DIFF=0d0
           DO J=1,3 
            DIFF=DIFF+(RCLSIMP(J,I)-TMP(J))**2
           ENDDO
           IF (DIFF.LT.1D-04) THEN
            DO J=1,3
             RN(J,I)=RR(J,N)
            ENDDO
           ENDIF
          ENDDO ! I loop
          ENDDO ! N loop
         ENDDO ! NA loop
         ENDIF
       ENDDO ! JATOM loop
       
     
       ALLOCATE(TAU_K_0(NKFS,2,2))
       ALLOCATE(SCAT_MAT(NKFS,2,NKFS,2))
       TAU_K_0=0d0

        IF (NSPIN==2) THEN
          WRITE (6,*) "magnetic impurities"
          WRITE (6,*) "Spin loop starts..."
        ELSE
          WRITE (6,*) "non-(or para-) magnetic impurities"
        END IF
  
        DO ISPIN=1,NSPIN
  
          IF (OPT('LIFETIME')) THEN
 
            WRITE(6,*) "before Lifetime_2D_SO"
            CALL LIFETIME_2D_SO(NKFS,LMAX,VOLBZ,SCAT_MAT,
     +           KF_IRR,NSPO,NSPOH,NSPIN,ISPIN,NCL_IMP,IMPLAYER,
     +           NAEZ,WEIGHTS_2D,TAU_K_0,RN,ATOMIMP,DOSW,FERMI_VL)
        
                write(6,*) "after lifetime  2F " 

          END IF

        END DO


        DEALLOCATE(KF_IRR)
        DEALLOCATE(WEIGHTS_2D)
        DEALLOCATE(TAU_K_0)
        DEALLOCATE(SCAT_MAT)


      END IF
      STOP 'after lifetime'

      END SUBROUTINE CALCULATE_LIFETIME     
