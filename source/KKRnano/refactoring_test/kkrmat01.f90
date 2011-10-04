    subroutine KKRMAT01(BZKP,NOFKS,GS,VOLCUB,VOLBZ, &
    TMATLL,MSSQ, &
    IE,IELAST,ITER, &
    ALAT,NSYMAT,NAEZ,CLS,NACLS,RR,EZOA,ATOM, &
    GINP,DGINP, &
    NUMN0,INDN0,IAT, &
    NATRC,ATTRC,EZTRC,NUTRC,INTRC, &
    SPRS,PRSC,EKM,NOITER, &
    EZ,QMRBOUND,IGUESS,BCP,CNVFAC, &
    DTDE_LOCAL, &
    GSXIJ, &
    NXIJ,XCCPL,IXCP,ZKRXIJ, &
    LLY_GRDT,TR_ALPH, &
    LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE, &
    LSMYRANK,LSRANK,LSMPIB,LSMPIC)

    use inc_p_wrapper_module
    use mpi

    implicit none
    ! ************************************************************************
    !   performs k-space integration,
    !   determines scattering path operator (g(k,e)-t**-1)**-1 and
    !   Greens function of the real system -> GS(*,*,*,*),

    !   NEW VERSION 10.99
    !   up -> left , down -> right, for decimation
    ! ------------------------------------------------------------------------

    !     .. parameters ..

    integer::        LMAX
    integer::NSYMAXD
    parameter      (LMAX=LMAXD,NSYMAXD=48)
    integer::        LMGF0D
    parameter      (LMGF0D= (LMAXD+1)**2)
    integer::        LMMAXD
    parameter      (LMMAXD= (LMAX+1)**2)
    integer::        ALM
    parameter      (ALM = NAEZD*LMMAXD)
    integer::        ALMGF0
    parameter      (ALMGF0 = NAEZD*LMGF0D)
    integer::        NGTBD
    parameter      (NGTBD = NACLSD*LMMAXD)
    integer::        NBLCKD
    parameter     (NBLCKD = XDIM*YDIM*ZDIM)
    integer::        LLYALM
    parameter      (LLYALM=LLY*(NAEZD*LMMAXD-1)+1)
    double complex :: CI
    parameter      (CI=(0.D0,1.D0))
    double complex :: CONE
    parameter      (CONE  = ( 1.0D0,0.0D0))
    double complex :: CIONE
    parameter      (CIONE  = ( 0.0D0,-1.0D0))
    double complex :: CZERO
    parameter      (CZERO=(0.0D0,0.0D0))
    !     ..
    !     .. GLOBAL SCALER ARGUMENTS ..
    double precision:: ALAT
    double precision::VOLBZ
    integer::NAEZ
    integer::NOFKS
    integer::NSYMAT
    integer::IGUESS
    integer::BCP
    integer::IE
    integer::IELAST
    integer::ITER
    integer::NXIJ
    integer::IXCP(NXIJD)
    integer::EKM
    integer::NOITER
    integer::IAT
    !     ..
    !     .. GLOBAL ARRAY ARGUMENTS ..
    double complex :: DGINP(LMGF0D,LMGF0D,NACLSD,NCLSD)
    double complex :: GINP(LMGF0D,LMGF0D,NACLSD,NCLSD)
    double complex :: GS(LMMAXD,LMMAXD,NSYMAXD)
    double complex :: GSXIJ(LMMAXD,LMMAXD,NSYMAXD,NXIJD)
    double complex :: EIKRP(NACLSD)
    double complex :: EIKRM(NACLSD)
    ! .. Lloyd
    double complex :: DTDE_LOCAL(LMMAXD,LMMAXD)
    double complex :: DGDE(LLYALM,LMMAXD)
    double complex :: GLLKE_X(LLYALM,LMMAXD)
    double complex :: MSSQ(LMMAXD,LMMAXD)
    double complex :: TR_ALPH
    double complex :: DPDE_LOCAL(LLYALM,LMMAXD)
    double complex :: LLY_GRDT
    ! .. precond
    double complex :: EZ(IEMXD)
    complex        :: PRSC(NGUESSD*LMMAXD,EKMD)
    double precision::BZKP(3,KPOIBZ)
    double precision::VOLCUB(KPOIBZ)
    double precision::RR(3,0:NRD)
    double precision::ZKRXIJ(48,3,NXIJD)
    double precision::CNVFAC(EKMD)
    integer:: NUMN0(NAEZD)
    integer:: INDN0(NAEZD,NACLSD)
    integer:: ATOM(NACLSD,*)
    integer:: CLS(*)
    integer:: EZOA(NACLSD,*)
    integer:: NACLS(*)
    integer:: SPRS(NGUESSD*LMMAXD+1,EKMD+1)
    integer:: NUTRC          !number of inequivalent atoms in the trunation cluster
    integer:: INTRC(NATRCD)  !pointer to atoms in the unit cell
    integer:: NATRC          !number of atoms in truncation cluster
    integer:: ATTRC(NATRCD)  !index to atom in elem/cell at site in truncation cluster
    integer:: EZTRC(NATRCD)  !index to bravais lattice  at site in truncation cluster
    !     ..
    !     .. LOCAL SCALARS ..
    double complex::   TRACE
    double complex::TRACEK
    double complex::GTDPDE
    double complex::BZTR2
    double complex::CFCTORINV
    double precision::TPI
    double precision::QMRBOUND
    double precision::DZNRM2
    integer::IC
    integer::I1
    integer::I2
    integer::ILM
    integer::ISYM
    integer::IU
    integer::IL1
    integer::IL1B
    integer::IL2
    integer::IL2B
    integer::KPT
    integer::LM
    integer::LM1
    integer::LM2
    integer::LM3
    integer::XIJ
    integer::ITCOUNT
    logical::XCCPL
    !     ..
    !     .. LOCAL ARRAYS ..
    double complex::G(LMMAXD,LMMAXD)
    double complex::TGH(LMMAXD)
    double complex::GLLKE1(ALM,LMMAXD)
    double complex::DUMMY(ALM,LMMAXD)
    double complex::GLLH(LMMAXD,NGTBD,NAEZD)
    double complex::TMATLL(LMMAXD,LMMAXD,NAEZD)
    double complex::GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD)
    double complex::X0(ALM,LMMAXD)
    double precision:: N2B(LMMAXD)
!     ..
!     .. EXTERNAL SUBROUTINES ..
    external CINIT,DLKE0,OUTTIME,DZNRM2
!     ..
!     .. INTRINSIC FUNCTIONS ..
    intrinsic ATAN,EXP
    !     ..
    !     .. MPI ..
    
    !     .. L-MPI
    integer:: MYLRANK(LMPID*SMPID*EMPID)
    integer:: &
        LCOMM(LMPID*SMPID*EMPID)
    integer:: &
        LGROUP(LMPID*SMPID*EMPID)
    integer:: &
        LSIZE(LMPID*SMPID*EMPID)
    integer::LMPIC
    !     .. LS-MPI
    integer:: LSMYRANK(LMPID,NAEZD*SMPID*EMPID)
    integer:: &
        LSRANK(LMPID,NAEZD*SMPID*EMPID)
    integer::LSMPIB
    integer::LSMPIC
    !     .. N-MPI
    integer:: IERR

!     ..
!-----------------------------------------------------------------------

    TPI = 8.D0*ATAN(1.D0)    ! = 2*PI
    CFCTORINV = (CONE*TPI)/ALAT

    BZTR2 = CZERO

    do 20 IU = 1,NSYMAXD
        call CINIT(LMMAXD**2,GS(1,1,IU))
    20 end do

    if (XCCPL) then
        do XIJ = 1, NXIJ
            do ISYM = 1,NSYMAT
                call CINIT(LMMAXD*LMMAXD,GSXIJ(1,1,ISYM,XIJ))
            enddo             ! ISYM = 1,NSYMAT
        enddo
    endif


! ---> use sit
!      G(n,n',L,L')(-k) = G(n',n,L',L)(k)


!=======================================================================
    do 300 KPT = 1, NOFKS                               ! K-POINT-LOOP
    !=======================================================================
    
    ! ---> fourier transformation
    
    !     added by h.hoehler 3.7.2002
    
    !                                                     n   0          n
    !     define fourier transform as g mu mu'= ( sum_n g mu mu' exp(-iKR )
    !                                   L  L'             L   L'
    
    !                                             n   0           n
    !                                 +   sum_n g mu'mu exp(-iK(-R ))) *0.5
    !                                             L'  L
    
    !     this operation has to be done to satisfy e.g. the point symmetry!
    !     application of fourier transformation is just an approximation
    !     for the tb system, since the transl. invariance is not satisfied.
    

    !================= Lloyd's Formula ====================================

        if (LLY == 1) then
            call CINIT(ALM*NGTBD,GLLH)

            do I1 = 1,NAEZ

                IC = CLS(I1)
                call DLKE1(ALAT,NACLS,RR,EZOA(1,I1), &
                BZKP(1,KPT),IC,EIKRM,EIKRP)
                call DLKE0(I1,GLLH,EIKRP,EIKRM, &
                IC,NACLS,ATOM(1,I1),NUMN0,INDN0,DGINP(1,1,1,IC))
            end do

            do I1=1,NAEZ
                do I2=1,NUMN0(I1)
                    do LM2=1,LMMAXD
                        IL2=LMMAXD*(I2-1)+LM2
                        if (INDN0(I1,I2) == IAT) then
                            do LM1=1,LMMAXD
                                IL1=LMMAXD*(I1-1)+LM1
                                DGDE(IL1,LM2)= GLLH(LM1,IL2,I1)
                            enddo
                        endif
                    enddo
                enddo
            enddo
        endif

    !============ END Lloyd's Formula =====================================
    

    !       Fourier transform of reference clusters' Green's function
    !       (from real space to k-space GINP -> GLLH)

        call CINIT(ALM*NGTBD,GLLH)

        do I1 = 1,NAEZ
            IC = CLS(I1)
            call DLKE1(ALAT,NACLS,RR,EZOA(1,I1), &
            BZKP(1,KPT),IC,EIKRM,EIKRP)
            call DLKE0(I1,GLLH,EIKRP,EIKRM, &
            IC,NACLS,ATOM(1,I1),NUMN0,INDN0, &
            GINP(1,1,1,IC))
        end do

    !=========== Lloyd's Formula ==========================================

        if (LLY == 1) then
            do I1=1,NAEZ
                do I2=1,NUMN0(I1)
                    do LM2=1,LMMAXD
                        IL2=LMMAXD*(I2-1)+LM2
                        if (INDN0(I1,I2) == IAT) then
                            do LM1=1,LMMAXD
                                IL1=LMMAXD*(I1-1)+LM1
                                GLLKE_X(IL1,LM2)= GLLH(LM1,IL2,I1)
                            enddo
                        endif
                    enddo
                enddo
            enddo
        end if

    !===========END Lloyd's Formula=========================================


    ! (!) The same calculation as some lines above is done all over again(!)
    ! ???

        call CINIT(ALM*NGTBD,GLLH)

        do I1 = 1,NAEZ
            IC = CLS(I1)
            call DLKE1(ALAT,NACLS,RR,EZOA(1,I1), &
            BZKP(1,KPT),IC,EIKRM,EIKRP)
            call DLKE0(I1,GLLH,EIKRM,EIKRP, &
            IC,NACLS,ATOM(1,I1),NUMN0,INDN0, &
            GINP(1,1,1,IC))
        end do

    ! -------------- Calculation of (Delta_t * G_ref - 1) ---------------
    ! INDN0 stores the index of the atom in the basis corresponding to
    ! the reference cluster atom

        do I1=1,NAEZ
            IL1B=LMMAXD*(I1-1)
            do I2=1,NUMN0(I1)
                do LM2=1,LMMAXD
                    IL2=LMMAXD*(I2-1)+LM2
                    IL2B=LMMAXD*(INDN0(I1,I2)-1)+LM2
                    do LM1=1,LMMAXD
                        TGH(LM1) = CZERO
                        do LM3=1,LMMAXD
                            TGH(LM1)=TGH(LM1)+TMATLL(LM1,LM3,I1)*GLLH(LM3,IL2,I1)
                        enddo
                    enddo
                    do LM1=1,LMMAXD
                        IL1=IL1B+LM1
                        GLLH(LM1,IL2,I1) = TGH(LM1)
                        if (IL1 == IL2B) then
                        !                    substract 1 only at the 'diagonal'
                            GLLH(LM1,IL2,I1) = GLLH(LM1,IL2,I1) - CONE
                        endif
                    enddo
                enddo
            enddo
        enddo


    ! ####################################################################
    ! Start Lloyd's formula here
    
    ! ####################################################################
    
    ! dP(E,k)   dG(E,k)                   dT(E)
    ! ------- = ------- * T(E) + G(E,k) * -----
    !   dE        dE                       dE

        if (LLY == 1) then

            call CINIT(LLYALM*LMMAXD,DPDE_LOCAL)

            call ZGEMM('N','N',ALM,LMMAXD,LMMAXD,CONE, &
            DGDE,ALM, &
            TMATLL(1,1,IAT),LMMAXD,CZERO, &
            DPDE_LOCAL,ALM)
            call ZGEMM('N','N',ALM,LMMAXD,LMMAXD,CFCTORINV, &
            GLLKE_X,ALM, &
            DTDE_LOCAL,LMMAXD,CONE,DPDE_LOCAL,ALM)

        endif
    !#######################################################################
    ! LLOYD
    !#######################################################################


    ! Now solve the linear matrix equation A*X = b (b is also a matrix),
    ! where A = (Delta_t*G_ref - 1) (inverse of scattering path operator)
    ! and b = (-1) * Delta_t

    ! If the initial guess optimisation is used, X and b are modified, but
    ! the form of the matrix equation stays the same

    !===================================================================
    !===================================================================
    
    ! solve linear matrix equation in 5 subsequent steps
    
    !===================================================================
    
    ! 1) find true residual tolerance by calculation of |b|
    !    meaning of true residual tolerance:
    !    the norm of |b| is used as convergence criterion
    !    note that if the initial guess optimisation is used, the matrix
    !    equation is solved using a modified b' - nevertheless the norm
    !    |b| of the original right hand side of the equation is used for
    !    convergence
    
        call CINIT(ALM*LMMAXD,DUMMY)
    
        do LM1=1,LMMAXD
            IL1=LMMAXD*(IAT-1)+LM1
            do LM2=1,LMMAXD
                DUMMY(IL1,LM2)=-TMATLL(LM1,LM2,IAT)
            enddo
        enddo

    !     store the norms of the columns of b in an array

        do LM2=1,LMMAXD
            N2B(LM2) = DZNRM2(NAEZD*LMMAXD,DUMMY(1,LM2),1)
        enddo
    
    ! ..
    !===================================================================
    
    !===================================================================
    ! 2) if IGUESS is activated perform intitial step 'I' of
    !    intial guess - store original b in DUMMY and set up modified b'
    
        if (IGUESS == 1) then
        
            call CINIT(ALM*LMMAXD,X0)
            call CINIT(ALM*LMMAXD,DUMMY)
        
            do I1=1,NAEZD
                do LM1=1,LMMAXD
                    IL1=LMMAXD*(I1-1)+LM1
                    do LM2=1,LMMAXD
                        DUMMY(IL1,LM2)=TMATLL(LM1,LM2,I1)
                    enddo
                enddo
            enddo
        
            if (ITER > 1) then
                call PRINVIT( &
                'I',IAT, &
                NUMN0,INDN0, &
                TMATLL,GLLH,X0, &
                PRSC(1,EKM+KPT),SPRS(1,EKM+KPT), &
                GLLKE1)
            endif
        endif

    ! ..
    !===================================================================
    
    !===================================================================
    ! 3) if BCP is activated determine preconditioning matrix
    !    GLLHBLCK ..
    
        call CINIT(LMMAXD*NATBLD*LMMAXD*NATBLD*NBLCKD,GLLHBLCK)
    
        if (BCP == 1) &
        call BCPWUPPER(GLLH,GLLHBLCK,NAEZ,NUMN0,INDN0)
    ! ..
    !===================================================================

    !===================================================================
    ! 4) solve linear set of equations by iterative TFQMR scheme
    
        call CINIT(ALM*LMMAXD,GLLKE1)
    
        call MMINVMOD( &
        GLLH,GLLKE1,TMATLL,NUMN0,INDN0,N2B, &
        IAT,ITER,ITCOUNT, &
        GLLHBLCK,BCP,IGUESS,CNVFAC(EKM+KPT), &
        QMRBOUND)
    
        NOITER = NOITER + ITCOUNT
    
    ! ..
    !===================================================================

    !===================================================================
    ! 5) if IGUESS is activated perform final step 'F' of
    !    intial guess           ~
    !                  X = X  + X
    !                       0          ..
    
        if (IGUESS == 1) then
            call PRINVIT( &
            'F',IAT, &
            NUMN0,INDN0, &
            TMATLL,GLLH,X0, &
            PRSC(1,EKM+KPT),SPRS(1,EKM+KPT), &
            GLLKE1)
        
            do I1=1,NAEZD
                do LM1=1,LMMAXD
                    IL1=LMMAXD*(I1-1)+LM1
                    do LM2=1,LMMAXD
                        TMATLL(LM1,LM2,I1)=DUMMY(IL1,LM2)
                    enddo
                enddo
            enddo
        
        endif
    ! ..
    !===================================================================
    
    ! solved. Result in GLLKE1
    
    !===================================================================
    !===================================================================



    !#######################################################################
    ! LLOYD calculate Trace of matrix ...
    !#######################################################################
        if (LLY == 1) then
        !===================================================================
        !                /  -1    dM  \
        ! calculate  Tr  | M   * ---- |
        !                \        dE  /
        !===================================================================
        

            TRACEK=CZERO

            do LM1=1,LMMAXD
                do LM2=1,LMMAXD
                    GTDPDE = CZERO
                    do IL1 = 1,LLYALM
                        GTDPDE = GTDPDE + GLLKE1(IL1,LM2)*DPDE_LOCAL(IL1,LM1)
                    enddo
                    TRACEK = TRACEK + MSSQ(LM1,LM2)*GTDPDE
                enddo
            enddo
        
            BZTR2 = BZTR2 + TRACEK*VOLCUB(KPT)
        
        endif
    !#######################################################################
    ! LLOYD .
    !#######################################################################
    
    !===================================================================
    
    
    

        ILM = LMMAXD*(IAT-1) + 1

    !                                      nn
    !         Copy the diagonal elements G_LL' of the Green's-function,
    !         dependent on (k,E) into matrix G
    !         (n = n' = IAT)

        do 140 LM = 1,LMMAXD
            call ZCOPY(LMMAXD,GLLKE1(ILM,LM),1,G(1,LM),1)
        140 end do

    !         Perform the k-space integration for diagonal element of
    !         Green's function of atom IAT

        do 110 ISYM = 1,NSYMAT
            do LM1=1,LMMAXD
                do LM2=1,LMMAXD
                    GS(LM1,LM2,ISYM) = GS(LM1,LM2,ISYM) + &
                    VOLCUB(KPT) * G(LM1,LM2)
                end do
            end do
            continue             ! ISYM = 1,NSYMAT
        110 end do
    
    
    
    ! ================================================================
        if (XCCPL) then
        ! ================================================================
        !       XCCPL communicate off-diagonal elements and multiply with
        !       exp-factor
        ! ================================================================
        
            call KKRJIJ( &
            BZKP,VOLCUB,KPT, &
            NSYMAT,NAEZ,IAT, &
            NXIJ,IXCP,ZKRXIJ, &
            GLLKE1, &
            GSXIJ, &
            LMPIC,MYLRANK, &
            LGROUP,LCOMM,LSIZE)
        
        ! ================================================================
        endif
    ! ================================================================
    
    
    
    !===================================================================
    
    
    !=======================================================================
300 end do ! KPT = 1,NOFKS
!=======================================================================



!=========== Lloyd's Formula =====================================

    if(LLY == 1)  then
        BZTR2 = BZTR2*NSYMAT/VOLBZ + TR_ALPH
        TRACE=CZERO
        call MPI_ALLREDUCE(BZTR2,TRACE,1, &
        MPI_DOUBLE_COMPLEX,MPI_SUM, &
        LCOMM(LMPIC),IERR)
        LLY_GRDT = TRACE
    endif
!========== END Lloyd's Formula ==================================

    return

    9100 format (2e24.5)

    end subroutine KKRMAT01
