      SUBROUTINE FERMIVL(LSTART,NBX,TMATLL,
     +                 EZ,DF,E2,ALAT,NSPIN,NSPOH,MAXMESH,NMESH,LMMAXSO,
     +                 IE,IELAST,IGF,KSCOEF,NSHELL,NBY,NBZ,
     +                 BBYA,CBYA,NAEZ,NATYP,CLS,EQINV,NACLS,RR,
     +                 RBASIS,EZOA,ATOM,RCLS,KAOEZ,LATT,ICC,GINP,
     +                 BRAVAIS,RECBV,LPOT,YR,WTYR,RIJ,IJEND,
     &                 ATOMIMP,IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT,
     &                 VOLBZ,INTERF,INS,IMPLAYER,
     &                 DELTALAMDAK,TMATLLE,DELTALAMDADELE,TSTAR,
     &                 PMUP,PMDOWN,SPERTURB)
      implicit none
c ************************************************************************
c calculate fermi velocity,modified by N. H. Long, Juelich, 24.02.11
c-----------------------------------------------------------------------
c     .. parameters ..
      include 'inc.p'
      include 'inc.cls'
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      DOUBLE PRECISION PI,RFCTOR
      PARAMETER (PI= 3.14159265358979312D0)
      INTEGER NKPOID
      PARAMETER (NKPOID=9000)
      DOUBLE COMPLEX CONEM,FAC
      PARAMETER (CONEM= (-1d0,0d0))
      INTEGER ALM,ALMSO
      PARAMETER (ALM=LMAXSQ*NAEZD,ALMSO=NAEZD*NSPD*LMAXSQ)
      INTEGER NSPBLOCK,NDIM
      PARAMETER (NSPBLOCK= 1, NDIM= NSPBLOCK*LMAXSQ)
      INTEGER NAUX
      PARAMETER (NAUX= 2*ALMSO**2+5*ALMSO)
c     ..
c     .. scalar arguments ..
      DOUBLE COMPLEX DF,EZ(IEMXD)          ! E-WEIGHT AND E-POINT
      DOUBLE PRECISION 
     +       BBYA,CBYA,             ! b/a, c/a
     +       E2,                    ! Fermi energy
     +       STIME0
      INTEGER ICC,IDGRP,IE,IELAST,IGF,IORIGIN,
     +        IJEND,INS,I1,
     +        KSCOEF,
     +        LATT,LPOT,
     +        MAXMESH,LMMAXSO,NSPOH,
     +        NAEZ,NATYP,NBY,NBZ,NSPIN,
     +        NMESH,M,
     +        NZ,                    ! number of atoms at centers of inversion
     +        NSYMAT,
     +        ILMC,NSYM,INFO,KSUM,LMSP1,IMPLAYER,
     +        IM,JN,IK,IS1,IS2,IL1,IL2,IL1T,IL2T,IENT,LM2TAKE,ENT0
      INTEGER  IATCONDL(*),IATCONDR(*),NCONDPAIR,IEGFOUT
      LOGICAL LSTART,LIRR,INTERF
c     ..
c     .. array arguments ..
      DOUBLE COMPLEX GINP(LMAXSQ*NACLSD,LMAXSQ,NCLSD,3),
     +               TMATLLE(LMMAXSO,LMMAXSO,NAEZD,3),
     +               TMATLL(LMMAXSO,LMMAXSO,*),RES,
     +               TSTAR(LMMAXSO,LMMAXSO),PMUP(NKPOID),PMDOWN(NKPOID)
      INTEGER ATOM(NACLSD,*),
     +        CLS(*),
     +        EQINV(*),EZOA(NACLSD,*),
     +        KAOEZ(*),
     +        NACLS(*),
     +        NSHELL(0:NSHELD),ATOMIMP(NATOMIMPD),
     &        HOSTIMP(0:natypd),
     +        IPVT(LMMAXSO),ENT(ALMSO),ENTW(ALMSO,ALMSO)
       DOUBLE PRECISION RBASIS(3,*),RR(3,*),RCLS(3,NACLSD,*),
     +                  BRAVAIS(3,3),RECBV(3,3)
       DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*),
     +                  DAUX(2*ALMSO)
C     ..
c     .. local scalars ..
      DOUBLE COMPLEX CCONST,DFZ,EK,TAUVBZ,NORM,NORMT,PROD
      DOUBLE PRECISION ALAT,BSC,CSC,R,STIME,VOLBZ,FACTINV,NORMABS,
     +                 MINLM1,MINLM2
      INTEGER I,II,IC,ID,IN,IU,IFILE1,ISHELL,
     +        J,IDK,DIMENK,IDELK,
     +        KS,L,LM,LM1,LM2,LMAXIN,M
     +        N,NBX,NDUM,NMESHOLD,NOFKS,NS,NR,NA,I1SP1,I2SP1,ISIGMA,
     +        I1SP2,I2SP2,IWRITE,
     +        STATUS_READ,STATUS_OPEN,LM1TAKE,LM3,LM4,NK1,NKNEAR,KSUM1
      INTEGER NATOMIMP,N1,N2,N3,K,NKF,NK,NKF_EV,ENTSO
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD),PHI
      CHARACTER*5 STRUCT
      CHARACTER*40 NAME,NEW
      CHARACTER*10 ROTNAME(64)
C     ..
C     .. LOCAL ARRAYS ..
      DOUBLE COMPLEX DSYMLL(LMMAXSO,LMMAXSO,NSYMD),
     +               AUX(NAUX)
      DOUBLE PRECISION BZKP1(3,KPOIBZ),BZKP(3,NKPOID),KP(3,NKPOID),
     +               RSYMAT(64,3,3),
     +               RATOM(3,NSHELD),
     +               RATOMS(3,NSHELD),
     +               RROT(48,3,NSHELD),
     +               VOLCUB(KPOIBZ),DKABS,DK(3,NKPOID),
     +               DKSOF1(3,NKPOID),DKSOF2(3,NKPOID),DV
      DOUBLE COMPLEX, ALLOCATABLE   ::  
c     +                                  GLLKE1(:,:),
c     +                                  GLLKE(:,:),
c     +                                  GLLKEDK(:,:),
c     +                                  GLLKEDE(:,:),
c     +                                  HILF_1(:,:),
c     +                                  HILF_2(:,:),
c     +                                  TMAT_1(:,:),
c     +                                  TMAT_2(:,:),
c     +                                  LV(:,:),RV(:,:),W(:),
     +                                  DENS(:,:,:,:,:,:,:),
     +                                  CKN(:,:),CKGES(:,:,:,:),
     +                                  RHOD(:,:,:,:),CKHELP(:),
     +                                  DENS_GESAMT_I1(:,:,:,:),
     +                                  DENS_GESAMT(:,:,:),S_12(:,:),
     +                                  CKGES1(:,:,:,:)
      DOUBLE COMPLEX GLLKE1(NAEZD*LMAXSQ,LMAXSQ),GLLKE(ALMSO,ALMSO),
     +               GLLKEDK(ALMSO,ALMSO),GLLKEDE(ALMSO,ALMSO),
     +               HILF_1(ALMSO,ALMSO),HILF_2(ALMSO,ALMSO),
     +               TMAT_1(ALMSO,ALMSO),TMAT_2(ALMSO,ALMSO),
     +               LV(ALMSO,ALMSO),RV(ALMSO,ALMSO),W(ALMSO)
      DOUBLE PRECISION V_FERMI(NKPOID,3)
      DOUBLE COMPLEX DELTALAMDAK(NKPOID,3),
     +               DELTALAMDAE(NKPOID),
     +               DELTALAMDADELK(NKPOID,3,2),
     +               DELTALAMDADELE(NKPOID,3),
     +               SPERTURB(2,4,NKPOID),CKHELP1
      INTEGER LF(LMAXSQ),
     +               NSH1(NSHELD),NSH2(NSHELD),
     +               ISYMINDEX(NSYMAXD),nxyz(3)
      logical SYMUNITARY(NSYMAXD)

C     ..
C     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION DCLOCK
      LOGICAL OPT,TEST,MAX
      EXTERNAL DCLOCK,DSORT,SHELLGEN2k,OPT,TEST,CINIT
      INTRINSIC DATAN,IDINT
C     ..
      EXTERNAL BZIRR3D,KKRMAT2001,RCSTOP,ROTBRILL,POINTGRP,FINDGROUP,
     +     TAUTOG1,ZAXPY,ZCOPY,ZGEMM
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DSQRT,SQRT
C     ..
C     .. SAVE STATEMENT ..
      SAVE
C     ..
C     ..
      IF (LSTART) THEN
        STIME = DCLOCK()
        RFCTOR = ALAT/(8.D0*DATAN(1.0D0))

        CALL POINTGRP(rsymat,rotname)

        IF (TEST('SETGROUP').OR.TEST('fullBZ  ')) THEN 

          CALL SETGROUP(bravais,recbv,rbasis,rfctor,naez,
     &                 rsymat,rotname,isymindex,nsymat)
          write(6,*) "Option setgroup is used"
          write(6,*) "nsymat:", nsymat
          write(6,*) "isymindex(1):", isymindex(1)
        ELSE       
          CALL FINDGROUP(bravais,recbv,rbasis,rfctor,naez,
     &                 rsymat,rotname,isymindex,nsymat)
          do i=1,NSYMAXD
             symunitary(i) = .false.
          end do
        END IF


c     find if the inversion is a symmetry op of the real lattice
c     This is going in sub kkrmat99 so that the transpose is
c     added or not !
c
        FACTINV = 1.d0
        do i=1,nsymat
          if (ISYMINDEX(i).eq.25) FACTINV=0.d0
c And in case of 2d              ! Changed on 20.01.2000 
          if (bravais(1,3).eq.0.d0.and.bravais(2,3).eq.0.d0.and.
     &        bravais(3,3).eq.0.d0) THEN
            IF (ISYMINDEX(i).eq.12) FACTINV= 0.d0
          END IF
        end do

        if (TEST('noinvers').OR.TEST('fullBZ  ')) FACTINV=0.d0
c        if (TEST('noinvers')) write(6,*) '!!!>>> Option NOINVERS Used'
        if (NSYMAT.GT.NSYMD) THEN
            write(6,*) 'Increase NSYMD in inc.p to : ',NSYMAT
        end if


        DSYMLL=0d0

        CALL ROTBRILL(DSYMLL(1:LMAXSQ,1:LMAXSQ,:),LMAX,NSYMAT,RSYMAT,
     &                ISYMINDEX,LPOT,YR,WTYR,RIJ,IJEND)
c
c Now DSYMLL hold NSYMAT symmetrization matrices
c
        CALL GFSHELLS(ICC,NATOMIMP,NSH1,NSH2,NSHELL,NAEZ,NATYP,
     &                    RBASIS,BRAVAIS,RATOM,RATOMS,RCLSIMP,
     &                    NSYMAT,ISYMINDEX,RSYMAT,RCLS,CLS,
     &                    EQINV,KAOEZ,NACLS,ATOM,ATOMIMP,RFCTOR,
     &                    IATCONDL,IATCONDR,NCONDPAIR,ROTNAME,
     &                    HOSTIMP)

c
c --->  creates difference vectors RROT for BZ integration in KKRMAT1
c
        DO 107 I=1,NSHELL(0)
          DO 108 IN=1,3
            RATOMS(IN,I) =RATOM(IN,I) 
 108      END DO
 107    END DO
 888    format(3F12.6)
c
        CALL CRTSTAR(RATOMS,NSHELL(0),RSYMAT,NSYMAT,ISYMINDEX,RROT)
c ------------------------------------------------------------------------
        IF (TEST('RROT    ')) THEN
          DO ISHELL=1,NSHELL(0)
            WRITE(6,FMT='(I4)') ISHELL
            WRITE(6,FMT='((I4,3F10.1))') 
     +           (IU,(RROT(IU,I,ISHELL),I=1,3),
     +           IU=1,NSYMAT)   ! 2*IUMAX*IVMAX)
          END DO
        END IF
c ------------------------------------------------------------------------
        LM = 0
        DO 100 L = 1,LMAX + 1
          DO 90 M = 1,L + L - 1
            LM = LM + 1
            LF(LM) = L
   90     CONTINUE
  100   CONTINUE

        IF (NBY.EQ.0) NBY = NBX/BBYA
        IF (NBZ.EQ.0) NBZ = NBX/CBYA
        IF (OPT('NBZ=0   ') .OR. OPT('SLAB    ')) NBZ = 0
        IF (OPT('WIRE    ')) THEN
          NBX = 0
          NBY = 0
        END IF
c
c k-point calculation or recover from file!
c
        NAME='mesh/kp'

        CALL SNAME(NAME,NEW,LATT)
        CALL SNAME(NEW,NAME,NBX)
        CALL SNAME(NAME,NEW,NBY)
        CALL SNAME(NEW,NAME,NBZ)

        IF (TEST('SETGROUP').OR.TEST('fullBZ  ')) GOTO 1000

        OPEN(52,FILE=NAME,status='old',ERR=1000)

        write(6,*) 'k-mesh file exist : ',name
        DO 242 L=1,MAXMESH
          READ (52,FMT='(I8,3f17.9)') NOFKS,VOLBZ,BSC,CSC
c ------------------------------------------------------------------------
          IF (NOFKS.GT.KPOIBZ) THEN
            write(6,*) 'Error in GLL: (NOFKS.GT.KPOIBZ)'
            write(6,*) 'Please increase the parameter KPOIBZ in ',
     +           'file ''inc.p''.'
            STOP 'GLL'
          END IF
c  ------------------------------------------------------------------------
          READ (52,*) ((BZKP1(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
          write(6,FMT=9020) L,NOFKS,VOLBZ
          IF (BSC .NE. 0.D0) THEN
            DO KS = 1,NOFKS
              BZKP1(2,KS)=BZKP1(2,KS) ! *BSC/BBYA        30.1.2002
              BZKP1(3,KS)=BZKP1(3,KS) ! *CSC/CBYA
            END DO
          END IF
          IF (TEST('k-net   ')) THEN
            DO KS = 1,NOFKS
              WRITE(6,9000) (BZKP1(I,KS),I=1,3),VOLCUB(KS)
            END DO
          END IF
 242    END DO
        GOTO 1010
       
 1000   write(6,*) 'Create k-mesh, write to file ',name

        DO 230 L=1,MAXMESH
          IF (L.GT.1) THEN
            NBX = (NBX)/1.4
            NBY = (NBY)/1.4
            NBZ = (NBZ)/1.4
            IF (NBX.EQ.0) NBX =1
            IF (NBY.EQ.0) NBY =1
            IF (NBZ.EQ.0) NBZ =1
          END IF

          LIRR=.TRUE.
          nxyz(1) = nbx
          nxyz(2) = nby
          nxyz(3) = nbz

          call BZIRR3D(NOFKS,nxyz,KPOIBZ,BZKP1,RECBV,BRAVAIS,RFCTOR,
     &                 VOLCUB,volbz,rsymat,nsymat,ISYMINDEX,LIRR)
c
           write(6,*) 'volbz',volbz 
c
c
           IF (L.EQ.1) OPEN(52,FILE=NAME)

           WRITE(52,FMT='(I8,3f17.9,/,(3f12.8,d20.10))') 
     +         NOFKS,VOLBZ,BBYA,CBYA,
     +         ((BZKP1(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
c ------------------------------------------------------------------------
c
c --->      output of k-mesh created
c
           IF (TEST('k-net   ')) THEN
             DO 200 KS = 1,NOFKS
               WRITE(6,9000) (BZKP1(I,KS),I=1,3),VOLCUB(KS) ! /VOLCUB(1)
 200         CONTINUE
           END IF
c ------------------------------------------------------------------------
 230     END DO                      !  L=1,MAXMESH
c
 1010    NMESHOLD = 0
c
C
C --->  ALAT=2*PI*RFCTOR
C
         ALAT = RFCTOR*8.D0*DATAN(1.0D0)       
C
C
         WRITE (6,*) 'IGF    = ',IGF
         IF (IGF.EQ.3) READ (81) NDUM

         WRITE (6,FMT=9070) DCLOCK()-STIME
 9070    format (' TIME IN BZMESH            : ',f9.2)



       END IF                        ! (LSTART)

       LSTART = .FALSE.

       IF (NMESH.NE.NMESHOLD) THEN
C         WRITE(6,*) '>>> READ K-MESH FROM FILE ''mesh.kp'' :'
         REWIND(52)
         DO 240 L=1,NMESH
          ! write(6,*) 'l, nmesh ',l,nmesh
           READ (52,FMT='(I8,3f17.9)') NOFKS,VOLBZ,BSC,CSC
           READ (52,*) ((BZKP1(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
           IF (BSC.NE.0.D0 .AND. L.EQ.NMESH) THEN
             DO KS = 1,NOFKS
               BZKP1(2,KS)=BZKP1(2,KS) ! *BSC/BBYA    30.1.2002
               BZKP1(3,KS)=BZKP1(3,KS) ! *CSC/CBYA
             END DO
           END IF                    ! BSC.NE.0.D0 .AND. L.EQ.NMESH
 240     END DO

         IF(TEST('flow    '))  
     +     write(6,*) "before test k mesh"

         IF (TEST('k-mesh  ')) THEN
           write(6,*) 'NMESH : ',NMESH
           WRITE(6,*) 'NOFKS : ',NOFKS
         END IF
       END IF                        ! (NMESH.NE.NMESHOLD)

       TAUVBZ = 1.D0/VOLBZ
C ------------------------------------------------------------------------
C ---> inverted t-matrix, save TMALL=TMALLE(E) 
C
       DO 140 I=1,NAEZ
                 
         DO LM1 = 1,LMMAXSO
         DO LM2 = 1,LMMAXSO
          TMATLLE(LM1,LM2,I,IE)=TMATLL(LM1,LM2,I) 
         ENDDO
         ENDDO

 140  CONTINUE                      ! 140 I=1,NAEZ

      DFZ = DF*RFCTOR**2
      STIME0 = DCLOCK()

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE(6,*) "before BAND-STR etc"

c      ALLOCATE (GLLKE1(NAEZD*LMAXSQ,LMAXSQ))
c      ALLOCATE (GLLKE(ALMSO,ALMSO))
c      ALLOCATE (GLLKEDK(ALMSO,ALMSO))
c      ALLOCATE (GLLKEDE(ALMSO,ALMSO))
c      ALLOCATE (HILF_1(ALMSO,ALMSO))
c      ALLOCATE (HILF_2(ALMSO,ALMSO))
c      ALLOCATE (TMAT_1(ALMSO,ALMSO))
c      ALLOCATE (TMAT_2(ALMSO,ALMSO))
c      ALLOCATE (RV(ALMSO,ALMSO))
c      ALLOCATE (LV(ALMSO,ALMSO))
c      ALLOCATE (W(ALMSO))
      IF (OPT('FERMI-VL')) THEN
c      open(unit=18, file="eigenvector", form='formatted')
      open(unit=29, file="vectorFERMIVL", form='formatted')
       IF (OPT('SOFIELD ')) THEN
        open(unit=25, file="SO_FIELD", form='formatted')
        open(unit=26, file="kshiftedSOF", form='formatted')
        open(unit=27, file="S_perturbed_up", form='formatted')
        open(unit=28, file="S_perturbed_dn", form='formatted')
        DV=1d-1
        WRITE(25,'(1e17.9)') DV
       ENDIF
c read the k_F from kpoint_FS 
      BZKP=0d0
      OPEN(unit=20, file="kpoints_FS",status="old",
     +            action="read",iostat=status_open)
      
      if ( status_open /= 0 ) stop "problems to open kpoints_FS"
      
      NKF=0      
      einlese_schleife: do
        NKF=NKF+1
        read (20,"(3(E17.9))",iostat=status_read)
     +       (BZKP(J,NKF),J=1,3)
      
        if ( status_read /= 0 ) exit   
      end do einlese_schleife

      NKF=NKF-1
      If(status_read < 0) then
        write (6,*) "Number of kpoints:", NKF
      end if

      CLOSE(20)
      KSUM=NKF
        write (6,*) "after read in kpoints fermi surface" 
        

c     calculate Green's function from GINP (summation in k space)
c     G(E,k)=sum_n{(exp(ikR)G_on(E)}
c
               
        IF (IE.EQ.1) THEN 
c        open(unit=19, file="coefficients_0", form='formatted')
        open(unit=19, file="coefficients_0", form='unformatted')
        open(unit=21, file="S_12", form='formatted')
        open(unit=22, file="S_xyz", form='formatted')
c        open(unit=23, file="rhod", form='formatted')
        open(unit=23, file="rhod", form='unformatted')
        open(unit=24, file="DENS", form='unformatted')
        IF (OPT('SOFIELD ').AND.TEST('wrtcoeff')) 
     +  open(unit=30, file="coefficients_so", form='unformatted')
        ALLOCATE (DENS(LMAXSQ,LMAXSQ,NSPOD,NSPOD,NSPOD,NSPOD,NATYPD))
        ALLOCATE (CKN(ALMSO,5))
        ALLOCATE (CKHELP(LMMAXSO))
        ALLOCATE (RHOD(LMMAXSO,LMMAXSO,NATYPD,4))
        ALLOCATE (CKGES(LMMAXSO,NATYPD,4,KSUM))
        ALLOCATE (CKGES1(LMMAXSO,NATYPD,2,KSUM))        
        ALLOCATE (DENS_GESAMT(4,KSUM,4))
        ALLOCATE (DENS_GESAMT_I1(NATYPD,4,KSUM,4))
        ALLOCATE (S_12(KSUM,4))
        DENS=CZERO
        CKN=CZERO
        CKGES=CZERO
        CKGES1=CZERO
        RHOD=CZERO
        DENS_GESAMT=CZERO
        DENS_GESAMT_I1=CZERO
        S_12=CZERO
c read in DENS
        READ(24) IWRITE
        DO I1=1,NATYPD
         DO I1SP1=1,NSPOD
          DO I2SP1=1,NSPOD
           DO I1SP2=1,NSPOD
            DO I2SP2=1,NSPOD
             DO LM1=1,LMAXSQ
              DO LM2=1,LMAXSQ
               READ(24) DENS(LM1,LM2,I1SP1,I2SP1,I1SP2,I2SP2,I1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CLOSE(24) 

        IF (NSPOH.NE.1) THEN
          DO ISIGMA=1,4
           DO I1=1,NATYPD
            IF (ISIGMA.EQ.1) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMAXSQ
                DO LM2=1,LMAXSQ
          RHOD((I1SP2-1)*LMAXSQ+LM2,(I1SP1-1)*LMAXSQ+LM1,I1,ISIGMA)=
     +           DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+ 
     +           DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ELSEIF (ISIGMA.EQ.2) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMAXSQ
                DO LM2=1,LMAXSQ
          RHOD((I1SP2-1)*LMAXSQ+LM2,(I1SP1-1)*LMAXSQ+LM1,I1,ISIGMA)=
     +             DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+ 
     +             DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ELSEIF (ISIGMA.EQ.3) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMAXSQ
                DO LM2=1,LMAXSQ
          RHOD((I1SP2-1)*LMAXSQ+LM2,(I1SP1-1)*LMAXSQ+LM1,I1,ISIGMA)=
c     +          (0D0,1D0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)- 
c     +                     DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)) 
     +          -(0D0,1D0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)-
     +                     DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)) 
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ELSEIF (ISIGMA.EQ.4) THEN
              DO I1SP1=1,NSPOD
               DO I1SP2=1,NSPOD
                DO LM1=1,LMAXSQ
                 DO LM2=1,LMAXSQ
          RHOD((I1SP2-1)*LMAXSQ+LM2,(I1SP1-1)*LMAXSQ+LM1,I1,ISIGMA)=
c     +            DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)- 
c     +            DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
     +           -DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+
     +            DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDIF
            ENDDO
           ENDDO
          DO ISIGMA=1,4
           DO I1=1,NATYPD
            DO LM1=1,LMMAXSO
             DO LM2=1,LMMAXSO
c              WRITE(23,'(2e17.9)') RHOD(LM1,LM2,I1,ISIGMA)
              WRITE(23) RHOD(LM1,LM2,I1,ISIGMA)
             ENDDO
            ENDDO
           ENDDO
          ENDDO        
         ENDIF
         CLOSE(23)
c        WRITE(19,'(I9)') KSUM
        WRITE(19) KSUM
        DKSOF1=0d0
        DKSOF2=0d0
        DO K=1,KSUM
c        DO K=1,1
        write(*,*) 'IE',K,IE
        GLLKE=CZERO

         DO I=1,NAEZ
          IC=CLS(KAOEZ(I))
          FAC=CONE
           IF (I.NE.EQINV(I)) FAC = CONEM      
             CALL CINIT(NAEZD*LMAXSQ*LMAXSQ,GLLKE1)
             CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,EZOA(1,I),ATOM(1,I),
     +             BZKP(1,K),IE,IC,FAC,GINP(1,1,IC,1),RCLS(1,1,IC))   
             DO M=1,LMAXSQ
               IM=(I-1)*LMAXSQ+M
                DO JN=1,LMAXSQ*NAEZ
                 GLLKE(JN,IM) = GLLKE(JN,IM)+ GLLKE1(JN,M)
                ENDDO
             ENDDO
          ENDDO ! atom loop
c  for spin-orbit coupling, G_LL'(k) 
        IF (NSPD == 2) THEN
          DO LM1=1,ALM
            DO LM2=1,ALM
              GLLKE(ALM+LM2,ALM+LM1) = GLLKE(LM2,LM1) 
            END DO
          END DO
        END IF
c       the matrix M=[1 - G^r*delta t] is constructed
        IS1=1
        IS2=1
        TMAT_1=CZERO
        TMAT_2=CZERO
        DO I1=1,NAEZ
          DO IS1=1,NSPD
            DO LM1=1,LMAXSQ
              DO IS2=1,NSPD
                DO LM2=1,LMAXSQ
                  IL1=(IS1-1)*ALM+LMAXSQ*(I1-1)+LM1
                  IL2=(IS2-1)*ALM+LMAXSQ*(I1-1)+LM2
                  IL1T=(IS1-1)*LMAXSQ+LM1
                  IL2T=(IS2-1)*LMAXSQ+LM2
                  TMAT_1(IL2,IL1)=TMATLLE(IL2T,IL1T,I1,IE)/RFCTOR
                  IF (I1.EQ.IMPLAYER) THEN
                   TMAT_2(IL2,IL1)=TSTAR(IL2T,IL1T)/RFCTOR
                  ELSE
                   TMAT_2(IL2,IL1)=CZERO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        LV=CZERO
        RV=CZERO
        W=CZERO
        HILF_1=CZERO
        
        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,GLLKE,ALMSO,
     +                           TMAT_1,ALMSO,CZERO,HILF_1,ALMSO)
        DO LM1=1,ALMSO
          DO LM2=1,ALMSO
            IF (LM1.EQ.LM2) THEN
              HILF_1(LM1,LM2)=CONE-HILF_1(LM1,LM2)
            ELSE
              HILF_1(LM1,LM2)=-HILF_1(LM1,LM2)
            ENDIF
          ENDDO
        ENDDO
c       eigenvalues determination
c       use lapack routine to calculate the eigenvalues of the
c       constructred matrix
c LV and RV are left or right eigen values; A*RV=lamda*RV;LV**H*A=lamda*LV**H 

        CALL ZGEEV('V','V',ALMSO,HILF_1,ALMSO,W,LV,ALMSO,RV,ALMSO,
     +              AUX,NAUX,DAUX,INFO)
        IENT=0
          MINLM1=1D6
        DO LM1=1,ALMSO
          IF (ABS(W(LM1)).LE.MINLM1) THEN
          MINLM1=ABS(W(LM1))
          LM1TAKE=LM1
          ENDIF
        ENDDO

        ENT(:)=1
        ENTW(:,:)=1
 
        DO LM2=1,ALMSO
          ENTW(1,LM2)=LM2
          DO LM1=LM2+1,ALMSO
            IF (ABS(W(LM1)-W(LM2)) .LT. 1.D-4 ) THEN 
              ENT(LM2)=ENT(LM2)+1
              ENT(LM1)=ENT(LM1)+1
              ENTW(ENT(LM2),LM2)=LM1
              ENTW(ENT(LM1),LM1)=LM2
            END IF
          END DO
        END DO
           
c calculate spin-orbit strength by purturbed t*
        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,GLLKE,ALMSO,
     +                           TMAT_2,ALMSO,CZERO,HILF_2,ALMSO)
        DO LM1=1,ALMSO
         DO LM2=1,ALMSO
          HILF_2(LM1,LM2)=-HILF_2(LM1,LM2)
         ENDDO
        ENDDO
       IF (OPT('SOFIELD ')) THEN
         CALL CALC_WAPR_SOF(W,ALMSO,LV,RV,HILF_2,
     +      ENT,ENTW,MAXVAL(ENT),PMUP(K),PMDOWN(K),RHOD,NATYPD,LMMAXSO,
     +      LMAXSQ,SPERTURB(1,1,K),CKGES(1,1,1,K))
        IF (TEST('wrtcoeff')) THEN
         DO IENT=1,2
          DO I1=1,NATYPD
           DO LM1=1,LMMAXSO
            WRITE(30) CKGES(LM1,I1,IENT,K)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDIF
       CKGES(:,:,:,K)=CZERO 
c normalized coefficients and write down to the files
        IENT=0
        MINLM1=1D6
        DO LM1=1,ALMSO
          IF (ABS(W(LM1)).LE.MINLM1) THEN
          MINLM1=ABS(W(LM1))
          LM1TAKE=LM1
          ENDIF
        ENDDO
c        WRITE(19,'((I9),(3e17.9))') K,(BZKP(J,K),J=1,3)
c        WRITE(19,'(I5)') ENT(LM1TAKE)
        WRITE(19) K,(BZKP(J,K),J=1,3)
        WRITE(19) ENT(LM1TAKE)
        DO IENT=1,ENT(LM1TAKE)
          IF (ABS(AIMAG(W(LM1TAKE))). LE .1d-02 
     +        .AND. ABS(REAL(W(LM1TAKE))). LE. 1D-02) THEN
           IF (NSPOH.EQ.1) THEN
            LM2=0
            DO I1=1,NATYPD
             DO LM1=1,LMAXSQ
              LM2=LM2+1
              CKGES(LM1,I1,IENT,K)=RV(LM2,ENTW(IENT,LM1TAKE))
             ENDDO
            ENDDO
           ELSE
            LM2=0
            DO I1=1,NATYPD
             DO LM1=1,LMAXSQ
              LM2=LM2+1
              CKGES(LM1,I1,IENT,K)=RV(LM2,ENTW(IENT,LM1TAKE))
             ENDDO
            ENDDO
            DO I1=1,NATYPD
             DO LM1=1,LMAXSQ
              LM2=LM2+1
              CKGES(LM1+LMAXSQ,I1,IENT,K)=RV(LM2,ENTW(IENT,LM1TAKE))
             ENDDO 
            ENDDO   
           ENDIF
           DO I1=1,NATYPD
            DO LM1=1,LMMAXSO
             LMSP1=(I1-1)*LMMAXSO+LM1
             CKN(LMSP1,IENT)=CKGES(LM1,I1,IENT,K)
            ENDDO
           ENDDO
          ENDIF
        ENDDO
          ENT0=2
           DO 101 IENT=2,ENT(LM1TAKE)          
            PROD=DOT_PRODUCT(CKN(:,1),CKN(:,J))
            IF (ABS(PROD).GT.1D-07) THEN
             ENT0=IENT
             GO TO 101
            ENDIF
101      CONTINUE 
           IF (ENT0.NE.2) THEN
            CKN(:,5)=CKN(:,2)
            CKN(:,2)=CKN(:,ENT0)
            CKN(:,ENT0)=CKN(:,5)
           ENDIF
c orthogonalize the coefficients if they are degenerated
           IF (ENT(LM1TAKE).NE.1) THEN
            PROD=DOT_PRODUCT(CKN(:,1),CKN(:,2)) 
             DO LM1=1,ALMSO
              CKN(LM1,2)=CKN(LM1,2)-PROD*CKN(LM1,1)
             ENDDO
            IF (ENT(LM1TAKE).GT.2) THEN
             PROD=DOT_PRODUCT(CKN(:,3),CKN(:,4)) 
              DO LM1=1,ALMSO
               CKN(LM1,4)=CKN(LM1,4)-PROD*CKN(LM1,3)
              ENDDO
             ENDIF
           ENDIF
         DO IENT=1,ENT(LM1TAKE)
           DO I1=1,NATYPD
            DO LM1=1,LMMAXSO
             LMSP1=(I1-1)*LMMAXSO+LM1
             CKGES(LM1,I1,IENT,K)=CKN(LMSP1,IENT)
            ENDDO
           ENDDO
c write down to the files
c           WRITE(19,'(I5)') IENT
           WRITE(19) IENT
            IF (NSPOH.EQ.1) THEN
             DO I1=1,NATYPD      
              DO LM1=1,LMAXSQ
               DO LM2=1,LMAXSQ
                DENS_GESAMT_I1(I1,IENT,K,1)=DENS_GESAMT_I1(I1,IENT,K,1)+
     +            CONJG(CKGES(LM1,I1,IENT,K))*CKGES(LM2,I1,IENT,K)*
     +                 DENS(LM1,LM2,1,1,1,1,I1) 
                 DENS_GESAMT(IENT,K,1)=DENS_GESAMT(IENT,K,1)+
     +             CONJG(CKGES(LM1,I1,IENT,K))*CKGES(LM2,I1,IENT,K)*
     +                 DENS(LM1,LM2,1,1,1,1,I1) 
               ENDDO
              ENDDO
             ENDDO
            ELSE 
             DO I1=1,NATYPD
              CKHELP=CZERO
               CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,CONE,RHOD(:,:,I1,1),
     +          LMMAXSO,CKGES(:,I1,IENT,K),LMMAXSO,CZERO,CKHELP,LMMAXSO)

               DENS_GESAMT_I1(I1,IENT,K,1)=
     +                        DOT_PRODUCT(CKGES(:,I1,IENT,K),CKHELP)
               DENS_GESAMT(IENT,K,1)=DENS_GESAMT(IENT,K,1)+
     +                                 DENS_GESAMT_I1(I1,IENT,K,1)
             ENDDO
            ENDIF
            DO I1=1,NATYPD
             DO LM1=1,LMMAXSO
               CKGES(LM1,I1,IENT,K)=CKGES(LM1,I1,IENT,K)/
     +                               CDSQRT(DENS_GESAMT(IENT,K,1))
c               WRITE(19,'(4(E20.12,2x))') CKGES(LM1,I1,IENT,K)
               WRITE(19) CKGES(LM1,I1,IENT,K)
             ENDDO
            ENDDO
          ENDDO
c calculate S_12 and S_xyz
         IF (NSPOH.NE.1) THEN
         DENS_GESAMT(:,K,:)=CZERO
         DENS_GESAMT_I1(:,:,K,:)=CZERO
          DO IENT=1,ENT(LM1TAKE)
           DO ISIGMA=1,4
             DO I1=1,NATYPD
              CKHELP=CZERO
          CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,CONE,RHOD(:,:,I1,ISIGMA),
     +          LMMAXSO,CKGES(:,I1,IENT,K),LMMAXSO,CZERO,CKHELP,LMMAXSO)

               DENS_GESAMT_I1(I1,IENT,K,ISIGMA)=
     +                        DOT_PRODUCT(CKGES(:,I1,IENT,K),CKHELP)
               DENS_GESAMT(IENT,K,ISIGMA)=DENS_GESAMT(IENT,K,ISIGMA)+
     +                                 DENS_GESAMT_I1(I1,IENT,K,ISIGMA)
             IF (IENT.EQ.2) THEN
               CKHELP=CZERO
          CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,CONE,RHOD(:,:,I1,ISIGMA),
     +          LMMAXSO,CKGES(:,I1,2,K),LMMAXSO,CZERO,CKHELP,LMMAXSO)
                 S_12(K,ISIGMA)=S_12(K,ISIGMA)+
     +                  DOT_PRODUCT(CKGES(:,I1,1,K),CKHELP)
             ENDIF
            ENDDO
           ENDDO
          ENDDO
          WRITE(21,'((I5),(7e17.9))') K,(S_12(K,ISIGMA),ISIGMA=2,4)
          WRITE(22,'((I9),(7e17.9))') K,(BZKP(J,K),J=1,2),
     +                      (REAL(DENS_GESAMT(1,K,ISIGMA)),ISIGMA=2,4)
         ENDIF
c normalized coeff

c  ---> check whether LV and RV are normed correctly
        
          DO LM2=1,ALMSO
            DO LM3=1,ALMSO
              NORM=0d0
              CALL ZGEMM('C','N',1,1,ALMSO,CONE,LV(:,LM3),ALMSO,
     +                     RV(:,LM2),ALMSO,CZERO,NORM,1)
              IF (LM2 .EQ. LM3 ) THEN
                PHI=atan(aimag(NORM)/real(NORM))
                NORMT=NORM*exp(-(0,1.d0)*PHI)
                NORMABS=SQRT(REAL(NORMT)**2+AIMAG(NORMT)**2)
                IF(REAL(NORMT)<0.D0) THEN
                  DO LM4=1,ALMSO
                    RV(LM4,LM2)=-EXP(-(0,1.D0)*PHI)/NORMABS*RV(LM4,LM2) 
                  END DO
                ELSE
                  DO LM4=1,ALMSO
                    RV(LM4,LM2)=EXP(-(0,1.D0)*PHI)/NORMABS*RV(LM4,LM2) 
                  END DO
                END IF
              END IF
            END DO  ! LM3=1,ALMSO
          END DO    ! LM2=1,ALMSO

c --->  calculate the gradient of the Green function at k times dk
c       GLLKEDK for extrapolating the eigenvalue at k+dk
        IF (INTERF .EQ. .TRUE.) THEN
        DIMENK=2
        ELSE
        DIMENK=3
        ENDIF
        IDK=0
        DO IDELK=1,2
         IF (IDELK.EQ.1) THEN
         DKABS=0.001
         ELSEIF (IDELK.EQ.2) THEN
         DKABS=-0.001
         ENDIF     

        DO IDK=1,DIMENK

        GLLKEDK=CZERO
        DELTALAMDAK(K,IDK)=CZERO
        KP=0d0
        DO J=1,3
        DK(J,K)=0d0
        ENDDO
         IF (IDK.EQ.1) THEN
          DK(IDK,K)=DKABS
         ELSE IF (IDK.EQ.2) THEN
          DK(IDK,K)=DKABS
         ELSE IF (IDK.EQ.3) THEN
          DK(IDK,K)=DKABS
         ENDIF
         DO IK=1,3
c          KP(IK,K) = BZKP(IK,K)+DKABS
          KP(IK,K) = BZKP(IK,K)+DK(IK,K)
         ENDDO
         DO I=1,NAEZ
          IC=CLS(KAOEZ(I))
          FAC=CONE
           IF (I.NE.EQINV(I)) FAC = CONEM      
           CALL CINIT(NAEZD*LMAXSQ*LMAXSQ,GLLKE1)
           CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,EZOA(1,I),ATOM(1,I),
     +                KP(1,K),IE,IC,FAC,GINP(1,1,IC,1),RCLS(1,1,IC))   
             DO M=1,LMAXSQ
               IM=(I-1)*LMAXSQ+M
                DO JN=1,LMAXSQ*NAEZ
                 GLLKEDK(JN,IM) = GLLKEDK(JN,IM)+ GLLKE1(JN,M)
                ENDDO
             ENDDO
         ENDDO ! atom loop
          DO LM2=1,ALM
            DO LM1=1,ALM
              GLLKEDK(LM1,LM2) = GLLKEDK(LM1,LM2)-GLLKE(LM1,LM2)
            END DO
          END DO
       
c  for spin-orbit coupling, dG_LL'(k) 
        IF (NSPD == 2) THEN
          DO LM1=1,ALM
            DO LM2=1,ALM
              GLLKEDK(ALM+LM2,ALM+LM1) = GLLKEDK(LM2,LM1) 
            END DO
          END DO
        END IF

c       the matrix dM=[ - dG*delta t] is constructed
        IS1=1
        IS2=1
        TMAT_1=CZERO
        HILF_1=CZERO
        DO I1=1,NAEZ
          DO IS1=1,NSPD
            DO LM1=1,LMAXSQ
              DO IS2=1,NSPD
                DO LM2=1,LMAXSQ
                  IL1=(IS1-1)*ALM+LMAXSQ*(I1-1)+LM1
                  IL2=(IS2-1)*ALM+LMAXSQ*(I1-1)+LM2
                  IL1T=(IS1-1)*LMAXSQ+LM1
                  IL2T=(IS2-1)*LMAXSQ+LM2
                  TMAT_1(IL2,IL1)=TMATLLE(IL2T,IL1T,I1,IE)/RFCTOR
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO


c calculate delta(lamda)/dk

        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,GLLKEDK,ALMSO,
     +                           TMAT_1,ALMSO,CZERO,HILF_1,ALMSO)
        DO LM1=1,ALMSO
          DO LM2=1,ALMSO
              HILF_1(LM1,LM2)=-HILF_1(LM1,LM2)
          ENDDO
        ENDDO


c pertubation method for dM --> dlamda/dk
         CALL CALC_WAPR_FV(W,ALMSO,LV,RV,HILF_1,
     +       ENT,ENTW,MAXVAL(ENT),EZ,DELTALAMDADELK(K,IDK,IDELK))       

       ENDDO ! dkx,dky,dkz loop
       ENDDO !IDELK loop
         DO J=1,3
          DELTALAMDAK(K,J)=(DELTALAMDADELK(K,J,1)-DELTALAMDADELK(K,J,2))
     +                        /(2d0*0.001)
         ENDDO
       ENDDO ! k loop 
      CLOSE(19)
      CLOSE(21)
      CLOSE(22)
      CLOSE(30)
      DEALLOCATE (DENS)
      DEALLOCATE (CKN)
      DEALLOCATE (CKHELP)
      DEALLOCATE (RHOD)
      DEALLOCATE (CKGES)
      DEALLOCATE (DENS_GESAMT)
      DEALLOCATE (DENS_GESAMT_I1)
      DEALLOCATE (S_12)

      ELSE IF (IE.GE.2) THEN  !energy E+dE
       IF (IE.EQ.3) THEN
        IF (OPT('SOFIELD ').AND.TEST('wrtcoeff')) THEN
        open(unit=30, file="coefficients_so", form='unformatted')
        open(unit=31, file="coefficients_so_rot", form='unformatted')
         WRITE(31) KSUM
        ENDIF
       ENDIF
        DO K=1,KSUM
c        DO K=1,1
        write(*,*) 'IE',K,IE

        GLLKE=CZERO
        GLLKEDE=CZERO

        DO I=1,NAEZ

          IC=CLS(KAOEZ(I))
          FAC=CONE
           IF (I.NE.EQINV(I)) FAC = CONEM     
           CALL CINIT(NAEZD*LMAXSQ*LMAXSQ,GLLKE1)
           CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,EZOA(1,I),ATOM(1,I),
     +               BZKP(1,K),IE,IC,FAC,GINP(1,1,IC,1),RCLS(1,1,IC))   
           DO M=1,LMAXSQ
             IM=(I-1)*LMAXSQ+M
              DO JN=1,LMAXSQ*NAEZ
               GLLKE(JN,IM) = GLLKE(JN,IM)+ GLLKE1(JN,M)
              ENDDO
           ENDDO
           CALL CINIT(NAEZD*LMAXSQ*LMAXSQ,GLLKE1)
           CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,EZOA(1,I),ATOM(1,I),
     +                BZKP(1,K),IE,IC,FAC,GINP(1,1,IC,IE),RCLS(1,1,IC)) 
             DO M=1,LMAXSQ
               IM=(I-1)*LMAXSQ+M
                DO JN=1,LMAXSQ*NAEZ
                 GLLKEDE(JN,IM) = GLLKEDE(JN,IM)+ GLLKE1(JN,M)
                ENDDO
             ENDDO
         ENDDO ! atom loop
c  for spin-orbit coupling, G_LL'(E),G_LL'(E+DE) 
        IF (NSPD == 2) THEN
          DO LM1=1,ALM
            DO LM2=1,ALM
              GLLKE(ALM+LM2,ALM+LM1) = GLLKE(LM2,LM1) 
              GLLKEDE(ALM+LM2,ALM+LM1) = GLLKEDE(LM2,LM1) 
            END DO
          END DO
        END IF
c       the matrix M=[1 - G^r*delta t] is constructed
        IS1=1
        IS2=1
        TMAT_1=CZERO
        TMAT_2=CZERO

        DO I1=1,NAEZ
          DO IS1=1,NSPD
            DO LM1=1,LMAXSQ
              DO IS2=1,NSPD
                DO LM2=1,LMAXSQ
                  IL1=(IS1-1)*ALM+LMAXSQ*(I1-1)+LM1
                  IL2=(IS2-1)*ALM+LMAXSQ*(I1-1)+LM2
                  IL1T=(IS1-1)*LMAXSQ+LM1
                  IL2T=(IS2-1)*LMAXSQ+LM2
                  TMAT_1(IL2,IL1)=TMATLLE(IL2T,IL1T,I1,1)/RFCTOR
                  TMAT_2(IL2,IL1)=TMATLLE(IL2T,IL1T,I1,IE)/RFCTOR
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
       DO LM1=1,ALMSO
        DO LM2=1,ALMSO
          GLLKEDE(LM1,LM2)=GLLKEDE(LM1,LM2)-GLLKE(LM1,LM2)
          TMAT_2(LM1,LM2)=TMAT_2(LM1,LM2)-TMAT_1(LM1,LM2)
        ENDDO
       ENDDO

c LV and RV are left or right eigen values; A*RV=lamda*RV;LV**H*A=lamda*LV**H 

        LV=CZERO
        RV=CZERO
        W=CZERO
        HILF_1=CZERO
        
        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,GLLKE,ALMSO,
     +                           TMAT_1,ALMSO,CZERO,HILF_1,ALMSO)
        DO LM1=1,ALMSO
          DO LM2=1,ALMSO
            IF (LM1.EQ.LM2) THEN
              HILF_1(LM1,LM2)=1d0-HILF_1(LM1,LM2)
            ELSE
              HILF_1(LM1,LM2)=-HILF_1(LM1,LM2)
            ENDIF
          ENDDO
        ENDDO
        CALL ZGEEV('V','V',ALMSO,HILF_1,ALMSO,W,LV,ALMSO,RV,ALMSO,
     +              AUX,NAUX,DAUX,INFO)

c  ---> check whether LV and RV are normed correctly
        
          DO LM2=1,ALMSO
            DO LM3=1,ALMSO
              NORM=0d0
              CALL ZGEMM('C','N',1,1,ALMSO,CONE,LV(:,LM3),ALMSO,
     +                     RV(:,LM2),ALMSO,CZERO,NORM,1)
              IF (LM2 .EQ. LM3 ) THEN
                PHI=atan(aimag(NORM)/real(NORM))
                NORMT=NORM*exp(-(0,1.d0)*PHI)
                NORMABS=SQRT(REAL(NORMT)**2+AIMAG(NORMT)**2)
                IF(REAL(NORMT)<0.D0) THEN
                  DO LM4=1,ALMSO
                    RV(LM4,LM2)=-EXP(-(0,1.D0)*PHI)/NORMABS*RV(LM4,LM2) 
                  END DO
                ELSE
                  DO LM4=1,ALMSO
                    RV(LM4,LM2)=EXP(-(0,1.D0)*PHI)/NORMABS*RV(LM4,LM2) 
                  END DO
                END IF
              END IF
            END DO  ! LM3=1,ALMSO
          END DO    ! LM2=1,ALMSO

        ENT(:)=1
        ENTW(:,:)=1
 
        DO LM2=1,ALMSO
          ENTW(1,LM2)=LM2
          DO LM1=LM2+1,ALMSO
            IF (ABS(W(LM1)-W(LM2)) .LT. 1.D-4 ) THEN 
              ENT(LM2)=ENT(LM2)+1
              ENT(LM1)=ENT(LM1)+1
              ENTW(ENT(LM2),LM2)=LM1
              ENTW(ENT(LM1),LM1)=LM2
            END IF
          END DO
        END DO

c       the matrix dM=[ - dG*delta t - G*dt] is constructed

        HILF_1=CZERO
        HILF_2=CZERO

        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,GLLKEDE,ALMSO,
     +                           TMAT_1,ALMSO,CZERO,HILF_1,ALMSO)
        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,GLLKE,ALMSO,
     +                           TMAT_2,ALMSO,CZERO,HILF_2,ALMSO)
        DO LM1=1,ALMSO
          DO LM2=1,ALMSO
              HILF_2(LM1,LM2)=-HILF_2(LM1,LM2)-HILF_1(LM1,LM2)
          ENDDO
        ENDDO

c pertubation method for dM --> dlamda/dE

         CALL CALC_WAPR_FV(W,ALMSO,LV,RV,HILF_2,
     +        ENT,ENTW,MAXVAL(ENT),EZ,DELTALAMDADELE(K,IE))

      IF (IE.EQ.3) THEN
       DELTALAMDAE(K)=(DELTALAMDADELE(K,2)-DELTALAMDADELE(K,3))
     +                /(2d0*0.001)
c     DELTALAMDAE(K)=DELTALAMDAE(K)/(EZ(2)-EZ(1))i
       DO J=1,3
        V_FERMI(K,J)=-DBLE(DELTALAMDAK(K,J)/DELTALAMDAE(K))
       ENDDO
        WRITE(29,"(3e17.9)") (V_FERMI(K,J),J=1,3)
      IF (OPT('SOFIELD ')) THEN
      IF (TEST('wrtcoeff')) THEN
      DO IENT=1,2
       DO I1=1,NATYPD
        DO LM1=1,LMMAXSO
         READ(30) CKGES1(LM1,I1,IENT,K)
        ENDDO
       ENDDO
      ENDDO
      ENDIF
      IF (-DBLE(PMDOWN(K)/DELTALAMDAE(K))
     +     .GT.-DBLE(PMUP(K)/DELTALAMDAE(K))) THEN
       WRITE(27,'((I9),(5e17.9))') K,(BZKP(J,K),J=1,2),
     +                      (REAL(SPERTURB(2,ISIGMA,K)),ISIGMA=2,4)
       WRITE(28,'((I9),(5e17.9))') K,(BZKP(J,K),J=1,2),
     +                      (REAL(SPERTURB(1,ISIGMA,K)),ISIGMA=2,4)
       DO I1=1,NATYPD
        DO LM1=1,LMMAXSO
         CKHELP1=CZERO
         CKHELP1=CKGES1(LM1,I1,1,K)
         CKGES1(LM1,I1,1,K)=CKGES1(LM1,I1,2,K)
         CKGES1(LM1,I1,2,K)=CKHELP1
        ENDDO
       ENDDO       
      ELSE
       WRITE(27,'((I9),(5e17.9))') K,(BZKP(J,K),J=1,2),
     +                      (REAL(SPERTURB(1,ISIGMA,K)),ISIGMA=2,4)
       WRITE(28,'((I9),(5e17.9))') K,(BZKP(J,K),J=1,2),
     +                      (REAL(SPERTURB(2,ISIGMA,K)),ISIGMA=2,4)
      ENDIF  

c write down coefficients SOFIELD 
       IF (TEST('wrtcoeff')) THEN
        ENTSO=2
        WRITE(31) K,(BZKP(J,K),J=1,3)
        WRITE(31) ENTSO
        DO IENT=1,ENTSO
         WRITE(31) IENT
         DO I1=1,NATYPD
          DO LM1=1,LMMAXSO
           WRITE(31) CKGES1(LM1,I1,IENT,K)
          ENDDO
         ENDDO
        ENDDO
        ENDIF
      RES=CZERO
      IF (-DBLE(PMDOWN(K)/DELTALAMDAE(K))
     +     .GT.-DBLE(PMUP(K)/DELTALAMDAE(K))) THEN
       RES=PMUP(K)
       PMUP(K)=PMDOWN(K)
       PMDOWN(K)=RES
      ENDIF  
      WRITE(25,'((I5),(2e17.9))') K,
     +         -DBLE(PMUP(K)/DELTALAMDAE(K)),
     +         -DBLE(PMDOWN(K)/DELTALAMDAE(K))
c       DV=1D-1
         DO J=1,2
          DKSOF1(J,K)=BZKP(J,K)+DBLE(PMUP(K)/DELTALAMDAE(K))*
     +                  DV*V_FERMI(K,J)/
     +                  (V_FERMI(K,1)**2+V_FERMI(K,2)**2)
          DKSOF2(J,K)=BZKP(J,K)+DBLE(PMDOWN(K)/DELTALAMDAE(K))*
     +                  DV*V_FERMI(K,J)/
     +                  (V_FERMI(K,1)**2+V_FERMI(K,2)**2)
         ENDDO
         WRITE(26,'(3e17.9)') (DKSOF1(J,K),J=1,3)
         WRITE(26,'(3e17.9)') (DKSOF2(J,K),J=1,3)
      ENDIF ! write down spin-orbit field
       ENDIF
      ENDDO   ! k loop
  
            
      CLOSE(30)
      CLOSE(31)
      ENDIF ! energy loop
c      CLOSE(18)
      CLOSE(25)
      CLOSE(26)
      CLOSE(27)
      CLOSE(28)
      CLOSE(29)
      ENDIF
c      DEALLOCATE (GLLKE)
c      DEALLOCATE (GLLKE1)
c      DEALLOCATE (GLLKEDK)
c      DEALLOCATE (GLLKEDE)
c      DEALLOCATE (HILF_1)
c      DEALLOCATE (HILF_2)
c      DEALLOCATE (TMAT_1)
c      DEALLOCATE (TMAT_2)
c      DEALLOCATE (RV)
c      DEALLOCATE (LV)
c      DEALLOCATE (W)


 9000 FORMAT(3F12.5,F15.8)
 9010 FORMAT(2F12.6)
 9020 FORMAT(' NMESH : ',I4,'  NOFKS : ',I7,'  VOLBZ :',f14.8)
 9030 FORMAT(I3,I7,I5,F14.3,2F11.3,I8,F8.3)
c      STOP

      RETURN

      END






