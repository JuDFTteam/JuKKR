module main2_aux_mod
  implicit none

  contains

  ! from dimension parameters, calculate some derived parameters
  !----------------------------------------------------------------------------
  subroutine getDerivedParameters(IGUESSD, IRMD, IRMIND, IRNSD, LASSLD, LM2D, &
                                  LMAXD, LMAXD1, LMMAXD, LMPOTD, LMXSPD, LPOTD, &
                                  LRECPOT, LRECRES2, MMAXD, NAEZD, NCLEB, &
                                  NGUESSD, NPOTD, NSPIND, NTIRD)
    implicit none
    integer :: IGUESSD
    integer :: IRMD
    integer :: IRMIND
    integer :: IRNSD
    integer :: LASSLD
    integer :: LM2D
    integer :: LMAXD
    integer :: LMAXD1
    integer :: LMMAXD
    integer :: LMPOTD
    integer :: LMXSPD
    integer :: LPOTD
    integer :: LRECPOT
    integer :: LRECRES2
    integer :: MMAXD
    integer :: NAEZD
    integer :: NCLEB
    integer :: NGUESSD
    integer :: NPOTD
    integer :: NSPIND
    integer :: NTIRD

    ! derived dimension parameters
    LPOTD = 2*LMAXD
    LMMAXD= (LMAXD+1)**2
    NPOTD=NSPIND*NAEZD
    LMAXD1=LMAXD+1
    MMAXD  = 2*LMAXD + 1
    LM2D= (2*LMAXD+1)**2
    LMXSPD= (2*LPOTD+1)**2
    LASSLD=4*LMAXD
    LMPOTD= (LPOTD+1)**2
    NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND
    IRMIND=IRMD-IRNSD

    NGUESSD = 1 + IGUESSD * ( NAEZD * (LMAXD+1)**2 - 1 )
    NCLEB = (LMAXD*2+1)**2 * (LMAXD+1)**2

    ! Record lengths
    LRECPOT=8*(LMPOTD*(IRNSD+1)+IRMD+20)
    LRECRES2=4+8*(NSPIND*(LMAXD+7)+2*LPOTD+4+2)
  end subroutine

  ! Consistency checks
  !----------------------------------------------------------------------------
  subroutine consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID)
    implicit none
    integer :: IEMXD
    integer :: LMAXD
    integer :: NSPIND
    integer :: SMPID

    ! -------------------------------------------------------------------------
    ! consistency checks
    if (IEMXD < 1) then
      write (*,*) "main2: IEMXD must be >= 1"
      stop
    end if

    if (LMAXD < 0) then
      write (*,*) "main2: LMAXD must be >= 0"
      stop
    end if

    if (SMPID /= 1 .and. SMPID /=2) then
      write (*,*) "main2: SMPID must be 1 or 2"
      stop
    end if

    if (NSPIND /= 1 .and. NSPIND /=2) then
      write (*,*) "main2: NSPIND must be 1 or 2"
      stop
    end if
  end subroutine

  !----------------------------------------------------------------------------
  subroutine consistencyCheck02(IELAST, IEMXD, IGUESS, IGUESSD, LMAX, LMAXD, &
                                NAEZ, NAEZD, NPNT1, NPNT2, NPNT3, NPOL, &
                                NR, NRD, NSPIN, NSPIND)
    implicit none
    integer :: IELAST
    integer :: IEMXD
    integer :: IGUESS
    integer :: IGUESSD
    integer :: LMAX
    integer :: LMAXD
    integer :: NAEZ
    integer :: NAEZD
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    integer :: NR
    integer :: NRD
    integer :: NSPIN
    integer :: NSPIND

    ! -------------- Consistency checks -------------------------------
    if (IGUESS /= IGUESSD .or. IGUESS < 0 .or. IGUESS > 1) then
      write (*,*) "main2: IGUESSD IGUESS inconsistent ", IGUESSD, IGUESS
      stop
    end if

    if (NSPIN /= NSPIND) then
      write (*,*) "main2: NSPIN /= NSPIND"
      stop
    end if

    if (NAEZ /= NAEZD) then
      write (*,*) "main2: NAEZ /= NAEZD"
      stop
    end if

    if (IEMXD /= IELAST) then
      write (*,*) "main2: IEMXD /= IELAST"
      stop
    end if

    if (NPOL /= 0) then
      if (NPNT1 + NPNT2 + NPNT3 + NPOL /= IELAST) then
        write(*,*) "main2: Energy point numbers inconsistent."
      end if
    end if

    if (LMAX /= LMAXD) then
      write (*,*) "main2: LMAX /= LMAXD"
      stop
    end if

    if (NSPIN /= NSPIND) then
      write (*,*) "main2: NSPIN /= NSPIND"
      stop
    end if

    if (NR > NRD .or. NR < 1) then
      write (*,*) "main2: NR inconsistent NR NRD ", NR, NRD
      stop
    end if
  end subroutine

  !----------------------------------------------------------------------------
  subroutine consistencyCheck03(ATOM, CLS, EZOA, INDN0, NACLS, NACLSD, NAEZ, NCLSD, NR, NUMN0)
    implicit none
    integer :: ATOM(:,:)
    integer :: CLS(:)
    integer :: EZOA(:,:)
    integer :: INDN0(:,:)
    integer :: NACLS(:)
    integer :: NACLSD
    integer :: NAEZ
    integer :: NCLSD
    integer :: NR
    integer :: NUMN0(:)

    integer :: I1
    integer :: IE

    do I1 = 1, NAEZ
      if (CLS(I1) < 1 .or. CLS(I1) > NCLSD) then
        write (*,*) "main2: CLS defect, site ", I1, " value ", CLS(I1)
        stop
      end if
    end do

    do I1 = 1, NCLSD
      if (NACLS(I1) < 1 .or. NACLS(I1) > NACLSD) then
        write (*,*) "main2: NACLS defect, cluster ", I1, " value ", NACLS(I1)
        stop
      end if
    end do

    do I1 = 1, NAEZ
      do IE = 1, NACLS(CLS(I1))
        if (ATOM(IE, I1) < 1 .or. ATOM(IE, I1) > NAEZ) then
          write (*,*) "main2: ATOM defect value ", IE, I1, ATOM(IE, I1)
          stop
        end if

        ! EZOA - index array to point to lattice vector (0 is allowed value)
        if (EZOA(IE, I1) < 0 .or. EZOA(IE, I1) > NR) then
          write (*,*) "main2: EZOA defect value ", IE, I1, EZOA(IE, I1)
          stop
        end if
      end do
    end do

    do I1 = 1, NAEZ
      if (NUMN0(I1) < 1 .or. NUMN0(I1) > NAEZ) then
        write (*,*) "main2: NUMN0 inconsistent site ", I1, " value ", NUMN0(I1)
        stop
      end if

      do IE = 1, NUMN0(I1)
        if (INDN0(I1, IE) < 1 .or. INDN0(I1, IE) > NAEZ) then
          write (*,*) "main2: INDN0 inconsistent site ", I1, " value ", INDN0(I1, IE)
          stop
        end if
      end do
    end do
  end subroutine


  !--------------------------------------------------------------------------------------
  ! read input generated by kkr0 from file inp.unf
  subroutine readKKR0Input(NSYMAXD, A, ALAT, ATOM, B, BCP, BRAVAIS, CLEB1C, &
                           CLS, DRDI, DSYMLL, EZOA, FCM, GMAX, GSH, ICLEB1C, ICST, &
                           IEND1, IFUNM, IGUESS, ILM, IMAXSH, IMIX, IMT, INDN0, IPAN, &
                           IRC, IRCUT, IRMIN, IRNS, IRWS, ISHIFT, ISYMINDEX, ITITLE, &
                           JEND, JIJ, KFORCE, KMESH, KPRE, KTE, KVMAD, KXC, LCORE, &
                           LDAU, LLMSP, LMAX, LMPOT, LMSP, LOFLM1C, LPOT, MAXMESH, &
                           MIXING, NACLS, NAEZ, NCLS, NCORE, NFU, NR, NREF, NSPIN, &
                           NSRA, NSYMAT, NTCELL, NUMN0, OPTC, QBOUND, QMRBOUND, R, &
                           RBASIS, RCLS, RCUTJIJ, RECBV, REFPOT, RMAX, RMT, RMTREF, &
                           RR, RWS, SCFSTEPS, TESTC, THETAS, VOLUME0, VREF, ZAT)

    implicit none
    integer :: NSYMAXD
    double precision, allocatable :: A(:)
    double precision :: ALAT
    integer, allocatable :: ATOM(:,:)
    double precision, allocatable :: B(:)
    integer :: BCP
    double precision :: BRAVAIS(3,3)
    double precision, allocatable :: CLEB1C(:,:)
    integer, allocatable :: CLS(:)
    double precision, allocatable :: DRDI(:,:)
    doublecomplex, allocatable :: DSYMLL(:,:,:)
    integer, allocatable :: EZOA(:,:)
    double precision :: FCM
    double precision :: GMAX
    double precision, allocatable :: GSH(:)
    integer, allocatable :: ICLEB1C(:,:)
    integer :: ICST
    integer :: IEND1
    integer, allocatable :: IFUNM(:,:)
    integer :: IGUESS
    integer, allocatable :: ILM(:,:)
    integer, allocatable :: IMAXSH(:)
    integer :: IMIX
    integer, allocatable :: IMT(:)
    integer, allocatable :: INDN0(:,:)
    integer, allocatable :: IPAN(:)
    integer, allocatable :: IRC(:)
    integer, allocatable :: IRCUT(:,:)
    integer, allocatable :: IRMIN(:)
    integer, allocatable :: IRNS(:)
    integer, allocatable :: IRWS(:)
    integer :: ISHIFT
    integer :: ISYMINDEX(NSYMAXD)
    integer, allocatable :: ITITLE(:,:)
    integer, allocatable :: JEND(:,:,:)
    logical :: JIJ
    integer :: KFORCE
    integer, allocatable :: KMESH(:)
    integer :: KPRE
    integer :: KTE
    integer :: KVMAD
    integer :: KXC
    integer, allocatable :: LCORE(:,:)
    logical :: LDAU
    integer, allocatable :: LLMSP(:,:)
    integer :: LMAX
    integer :: LMPOT
    integer, allocatable :: LMSP(:,:)
    integer, allocatable :: LOFLM1C(:)
    integer :: LPOT
    integer :: MAXMESH
    double precision :: MIXING
    integer, allocatable :: NACLS(:)
    integer :: NAEZ
    integer :: NCLS
    integer, allocatable :: NCORE(:)
    integer, allocatable :: NFU(:)
    integer :: NR
    integer :: NREF
    integer :: NSPIN
    integer :: NSRA
    integer :: NSYMAT
    integer, allocatable :: NTCELL(:)
    integer, allocatable :: NUMN0(:)
    character(len=8) :: OPTC(8)
    double precision :: QBOUND
    double precision :: QMRBOUND
    double precision, allocatable :: R(:,:)
    double precision, allocatable :: RBASIS(:,:)
    double precision, allocatable :: RCLS(:,:,:)
    double precision :: RCUTJIJ
    double precision :: RECBV(3,3)
    integer, allocatable :: REFPOT(:)
    double precision :: RMAX
    double precision, allocatable :: RMT(:)
    double precision, allocatable :: RMTREF(:)
    double precision, allocatable :: RR(:,:)
    double precision, allocatable :: RWS(:)
    integer :: SCFSTEPS
    character(len=8) :: TESTC(16)
    double precision, allocatable :: THETAS(:,:,:)
    double precision :: VOLUME0
    double precision, allocatable :: VREF(:)
    double precision, allocatable :: ZAT(:)

    ! ======================================================================
    ! =             read in variables from unformatted files               =
    ! ======================================================================
    open (67,file='inp.unf',form='unformatted')
    read (67) KMESH,MAXMESH,RR,EZOA,NUMN0,INDN0,NSYMAT,DSYMLL
    read (67) NAEZ,NSPIN,IPAN,IRNS,IRCUT,LCORE,NCORE,NTCELL, &
    LMAX,LPOT,LMPOT
    read (67) IMIX,MIXING,QBOUND,FCM,KPRE,KTE, &
    KVMAD,KXC,ISHIFT,KFORCE,IGUESS,BCP,QMRBOUND
    read (67) A,B,DRDI,R,THETAS,ZAT,IMT,IRC, &
    IRMIN,IRWS,RWS,RMT,ITITLE,LLMSP,NFU
    read (67) ALAT,GSH,ILM,IMAXSH,TESTC,OPTC
    read (67) BRAVAIS,RBASIS,RECBV,VOLUME0,RMAX,GMAX
    read (67) IEND1,CLEB1C,ICLEB1C,LOFLM1C,JEND,IFUNM,LMSP,NSRA,ICST
    read (67) NCLS,NREF,RMTREF,VREF,RCLS, &
    ATOM,CLS,NACLS,REFPOT
    read (67) NR,RCUTJIJ,JIJ,LDAU
    read (67) ISYMINDEX,SCFSTEPS
    close (67)
  end subroutine

  !----------------------------------------------------------------------------
  ! Read k-mesh file
  subroutine readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)
    implicit none
    double precision :: BZKP(:,:,:)
    integer :: MAXMESH
    integer :: NOFKS(:)
    double precision :: VOLBZ(:)
    double precision :: VOLCUB(:,:)

    ! -----------------------------
    integer :: I
    integer :: ID
    integer :: L

    open (52,file='kpoints',form='formatted')
    rewind (52)

    do L = 1,MAXMESH
      read (52,fmt='(I8,f15.10)') NOFKS(L),VOLBZ(L)
      read (52,fmt=*) (BZKP(ID,1,L),ID=1,3),VOLCUB(1,L)
      do I=2,NOFKS(L)
        read (52,fmt=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
      end do
    end do

    close (52)
  end subroutine

  !----------------------------------------------------------------------------
  !> Print Fermi-Energy information to screen.
  subroutine printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
    implicit none
    double precision :: DENEF
    double precision :: E2
    double precision :: E2SHIFT
    double precision :: EFOLD
    integer :: NAEZ

    write (6,fmt=9020) EFOLD,E2SHIFT

    ! --> divided by NAEZ because the weight of each atom has been already
    !     taken into account in 1c

    write (6,fmt=9030) E2,DENEF/DBLE(NAEZ)
    write(6,'(79(1H+),/)')
9020 format ('                old', &
    ' E FERMI ',F12.6,' Delta E_F = ',f12.6)
9030 format ('                new', &
    ' E FERMI ',F12.6,'  DOS(E_F) = ',f12.6)
  end subroutine

end module main2_aux_mod
