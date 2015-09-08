module EnergyMeshHelpers_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: emesht, epathtb
  
  public :: readEnergyMeshImpl, writeEnergyMeshImpl, updateEnergyMeshImpl, broadcastEnergyMeshImpl_com
  public :: readEnergyMeshImplSemi, writeEnergyMeshImplSemi, updateEnergyMeshImplSemi

  contains

  ! VALENCE CONTOUR ONLY!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh.0'
  subroutine readEnergyMeshImpl(e1, e2, efermi, ez, ielast, npnt1, npnt2, npnt3, npol, tk, wez)
    double precision, intent(out) :: e1
    double precision, intent(out) :: e2
    double precision, intent(out) :: efermi
    double complex, intent(out) :: ez(:)
    integer, intent(out) :: ielast
    integer, intent(out) :: npnt1
    integer, intent(out) :: npnt2
    integer, intent(out) :: npnt3
    integer, intent(out) :: npol
    double precision, intent(out) :: tk
    double complex, intent(out) :: wez(:)

    open (67, file='energy_mesh.0', form='unformatted', action='read', status='old')
    read (67) ielast,ez,wez,e1,e2
    read (67) npol,tk,npnt1,npnt2,npnt3
    ! if (npol == 0) read(67) efermi
    read (67) efermi
    close(67)
  endsubroutine ! read

  ! VALENCE CONTOUR ONLY!
  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMeshImpl(e1, e2, efermi, ez, ielast, npnt1, npnt2, npnt3, npol, tk, wez)
    double precision, intent(in) :: e1
    double precision, intent(in) :: e2
    double precision, intent(in) :: efermi
    double complex, intent(in) :: ez(:)
    integer, intent(in) :: ielast
    integer, intent(in) :: npnt1
    integer, intent(in) :: npnt2
    integer, intent(in) :: npnt3
    integer, intent(in) :: npol
    double precision, intent(in) :: tk
    double complex, intent(in) :: wez(:)

    open  (67, file='energy_mesh', form='unformatted', action='write')
    write(67) ielast,ez,wez,e1,e2
    write(67) npol,tk,npnt1,npnt2,npnt3
    write(67) efermi
    close (67)
    
  endsubroutine ! write

  ! VALENCE CONTOUR ONLY!
  !------------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EMESHT
  subroutine updateEnergyMeshImpl(ez,wez,ielast,e1,e2,tk,npol,npnt1,npnt2,npnt3)
  
    double complex, intent(out) :: ez(:)
    double complex, intent(out) :: wez(:)
    integer, intent(inout) :: ielast
    double precision, intent(in) :: e1
    double precision, intent(in) :: e2
    double precision, intent(in) :: tk
    integer, intent(in) :: npol
    integer, intent(in) :: npnt1
    integer, intent(in) :: npnt2
    integer, intent(in) :: npnt3
    
    integer :: ie, iemxd
    double precision :: pi
    
    pi = 4.0d0*atan(1.0d0)
    iemxd = ielast

    ! --> update energy contour
    call emesht(ez, wez, ielast, e1, e2, e2, tk, npol, npnt1, npnt2, npnt3, iemxd)
    ! ielast will be overwritten

    do ie = 1, ielast
      wez(ie) = -2.d0/pi*wez(ie)
    enddo ! ie
    
  endsubroutine ! update


  ! VALENCE CONTOUR ONLY!
  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from rank 'BCRANK' to all other ranks
  subroutine broadcastEnergyMeshImpl_com(actvcomm, bcrank, e1, e2, ez, iemxd, wez)
    integer, intent(in) :: actvcomm
    integer, intent(in) :: bcrank
    double precision, intent(inout) :: e1
    double precision, intent(inout) :: e2
    double complex, intent(inout) :: ez(:)
    integer, intent(in) :: iemxd
    double complex, intent(inout) :: wez(:)

    include 'mpif.h'
    integer :: ierr

    call MPI_Bcast(ez,iemxd,mpi_double_complex, bcrank,actvcomm,ierr)
    call MPI_Bcast(wez,iemxd,mpi_double_complex, bcrank,actvcomm,ierr)
    call MPI_Bcast(e1,1,mpi_double_precision, bcrank,actvcomm,ierr)
    call MPI_Bcast(e2,1,mpi_double_precision, bcrank,actvcomm,ierr)
  endsubroutine ! broadcast


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh.0'
  subroutine readEnergyMeshImplSemi(e1, e2, efermi, ez, ielast, npnt1, npnt2, npnt3, npol, &
                   tk, wez, ebotsemi, emusemi, fsemicore, iesemicore, n1semi, n2semi, n3semi)
    ! valence contour parameters
    double precision, intent(out) :: e1
    double precision, intent(out) :: e2
    double precision, intent(out) :: efermi
    double complex, intent(out) :: ez(:)
    integer, intent(out) :: ielast
    integer, intent(out) :: npnt1
    integer, intent(out) :: npnt2
    integer, intent(out) :: npnt3
    integer, intent(out) :: npol
    double precision, intent(out) :: tk
    double complex, intent(out) :: wez(:)

    ! semicore contour parameters
    double precision, intent(out) :: ebotsemi
    double precision, intent(out) :: emusemi
    double precision, intent(out) :: fsemicore
    integer, intent(out) :: iesemicore
    integer, intent(out) :: n1semi
    integer, intent(out) :: n2semi
    integer, intent(out) :: n3semi

    open (67,file='energy_mesh.0',form='unformatted', action='read', status='old')
    read (67) ielast,ez,wez,e1,e2
    read (67) npol,tk,npnt1,npnt2,npnt3
    read (67) efermi
    read (67) iesemicore,fsemicore,ebotsemi
    read (67) emusemi
    read (67) n1semi,n2semi,n3semi
    close(67)
    
  endsubroutine ! read

  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMeshImplSemi(e1, e2, efermi, ez, ielast, npnt1, npnt2, npnt3, npol, &
                    tk, wez, ebotsemi, emusemi, fsemicore, iesemicore, n1semi, n2semi, n3semi)
    ! valence contour parameters
    double precision, intent(in) :: e1
    double precision, intent(in) :: e2
    double precision, intent(in) :: efermi
    double complex, intent(in) :: ez(:)
    integer, intent(in) :: ielast
    integer, intent(in) :: npnt1
    integer, intent(in) :: npnt2
    integer, intent(in) :: npnt3
    integer, intent(in) :: npol
    double precision, intent(in) :: tk
    double complex, intent(in) :: wez(:)

    ! semicore contour parameters
    double precision, intent(in) :: ebotsemi
    double precision, intent(in) :: emusemi
    double precision, intent(in) :: fsemicore
    integer, intent(in) :: iesemicore
    integer, intent(in) :: n1semi
    integer, intent(in) :: n2semi
    integer, intent(in) :: n3semi

    open  (67,file='energy_mesh',form='unformatted', action='write')
    write(67) ielast,ez,wez,e1,e2
    write(67) npol,tk,npnt1,npnt2,npnt3
    write(67) efermi
    write(67) iesemicore,fsemicore,ebotsemi
    write(67) emusemi
    write(67) n1semi,n2semi,n3semi
    close (67)
    
  endsubroutine ! write

  ! VALENCE AND SEMICORE CONTOUR!
  !------------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EPATHTB
  subroutine updateEnergyMeshImplSemi(ez,wez,ielast,e1,e2,tk,npol,npnt1,npnt2,npnt3, &
                                      ebotsemi,emusemi,iesemicore,fsemicore,n1semi,n2semi,n3semi)
    ! valence contour parameters
    double complex, intent(out) :: ez(:)
    double complex, intent(out) :: wez(:)
    integer, intent(inout) :: ielast
    double precision, intent(in) :: e1
    double precision, intent(in) :: e2
    double precision, intent(in) :: tk
    integer, intent(in) :: npol
    integer, intent(in) :: npnt1
    integer, intent(in) :: npnt2
    integer, intent(in) :: npnt3
    
    ! semicore contour parameters
    double precision, intent(in) :: ebotsemi
    double precision, intent(in) :: emusemi
    double precision, intent(in) :: fsemicore
    integer, intent(out) :: iesemicore
    integer, intent(in) :: n1semi
    integer, intent(in) :: n2semi
    integer, intent(in) :: n3semi

    double precision :: pi
    integer :: iemxd
    integer :: ie

    pi = 4.0d0*atan(1.0d0)
    iemxd = ielast

    ! --> update energy contour
    call epathtb(ez, wez, e2, ielast, iesemicore, 1, e1, e2, tk, npol, npnt1, npnt2, npnt3, &
                                      ebotsemi, emusemi, tk, npol, n1semi, n2semi, n3semi, iemxd)

    do ie = 1, ielast
      wez(ie) = -2.d0/pi*wez(ie)
      if (ie <= iesemicore) wez(ie) = wez(ie)*fsemicore
    enddo ! ie

  endsubroutine ! update


  ! VALENCE AND SEMICORE CONTOUR!
  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from rank 'BCRANK' to all other ranks
  subroutine broadcastEnergyMeshImplSemi_com(actvcomm, bcrank, e1, e2, ez, iemxd, wez, ebotsemi, emusemi)
    ! valence contour parameters
    integer, intent(in) :: actvcomm
    integer, intent(in) :: bcrank
    double precision, intent(inout) :: e1
    double precision, intent(inout) :: e2
    double complex, intent(inout) :: ez(:)
    integer, intent(in) :: iemxd
    double complex, intent(inout) :: wez(:)

    ! valence contour parameters
    double precision, intent(inout) :: ebotsemi
    double precision, intent(inout) :: emusemi
    
    include 'mpif.h'
    integer :: ierr

    call MPI_Bcast(ez,iemxd,mpi_double_complex, bcrank,actvcomm,ierr)
    call MPI_Bcast(wez,iemxd,mpi_double_complex, bcrank,actvcomm,ierr)
    call MPI_Bcast(e1,1,mpi_double_precision, bcrank,actvcomm,ierr)
    call MPI_Bcast(e2,1,mpi_double_precision, bcrank,actvcomm,ierr)
    call MPI_Bcast(ebotsemi,1,mpi_double_precision, bcrank,actvcomm,ierr)
    call MPI_Bcast(emusemi,1,mpi_double_precision, bcrank,actvcomm,ierr)
    
  endsubroutine ! broadcast

  
! 06.10.09 *************************************************************
!     SUBROUTINE EMESHT(EZ,DF,NPNT,EBOT,EMU,EFERMI,TK,
!    &                  NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
! **********************************************************************
! *                                                                    *
! * This subroutine provides the energy mesh in array EZ and the       *
! * appropriate integration weights in array DF for contour            *
! * integration in the complex energy plane.                           *
! *                                                                    *
! * The contour consists of two straight lines determined by two       *
! * real valued input arguments: TK and EBOT and by a parameter NPOL   *
! * which is determined in this subroutine (see below).                *
! *                                                                    *
! *         TK   = temperature in K                                    *
! *         EBOT = bottom of contour in Ry                             *
! *         NPOL = number of used poles of the Fermi-Dirac function    *
! *         EMU =  chemical potential in Ry                            *
! *                                                                    *
! * One line begins at EBOT on the real energy axis and ends at        *
! * EBOT + 2*NPOL*i*pi*k*TK. On this line NPNT1 mesh points are used   *
! * which are determined according to a Gauss-Legendre rule.           *
! *                                                                    *
! * The second line begins at EBOT+2*NPOL*i*pi*k*TK and goes parallel  *
! * to the real axis to infinity. This line is divided into            *
! * NLEG + 2 parts (NLEG can be specified below, usually NLEG = 1 is   *
! * enough, but NLEG can be increased for large values of EMU-EBOT.    *
! * On the first NLEG parts Gauss-Legendre integration with NPNT2      *
! * mesh points each are used and the Fermi-Dirac function is replaced *
! * by one since it differs from one by less than 10**(-13).           *
! * On the remaining two parts Gauss-Fermi Dirac integrations rules    *
! * with NPNT3 and NPNT4 mesh points are used. These parts are         *
! * separated by EMU + 2*NPOL*i*pi*k*TK or by the point where the      *
! * weight function (the difference of two Fermi-Dirac functions)      *
! * changes sign. (see value of XZERO in GAUFDT2 and GAUFDT3)          *
! *                                                                    *
! * useful values for NPNT1, NPNT2 and NPNT3 depend on EBOT, TK and on *
! * the option EMESHT1,  EMESHT2 or EMESHT3. For Ni, Cu and Pd useful  *
! * values are NPNT1 = 7 and specifically                              *
! * for EMESHT1                                                        *
!        TK = 3200 K    EBOT = -0.56   NPNT2 =  8   NPNT3 = 12         *
!        TK = 2400 K    EBOT = -0.56   NPNT2 = 10   NPNT3 =  8         *
!        TK = 1600 K    EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
!        TK = 1200 K    EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
!        TK = 800 K     EBOT = -0.56   NPNT2 = 16   NPNT3 =  8         *
!        TK = 400 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
!        TK = 200 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
! * for EMESHT2                                                        *
!        TK = 3200 K    EBOT = -1.36   NPNT2 =  6   NPNT3 = 20         *
!        TK = 2400 K    EBOT = -0.96   NPNT2 =  6   NPNT3 = 16         *
!        TK = 1600 K    EBOT = -0.56   NPNT2 =  8   NPNT3 = 12         *
!        TK = 1200 K    EBOT = -0.56   NPNT2 = 10   NPNT3 =  8         *
!        TK = 800 K     EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
!        TK = 400 K     EBOT = -0.56   NPNT2 = 16   NPNT3 =  8         *
!        TK = 200 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
! * for EMESHT3                                                        *
!        TK = 3200 K    EBOT = -1.36   NPNT2 =  6   NPNT3 = 20         *
!        TK = 2400 K    EBOT = -0.96   NPNT2 =  6   NPNT3 = 20         *
!        TK = 1600 K    EBOT = -0.56   NPNT2 =  8   NPNT3 = 12         *
!        TK = 1200 K    EBOT = -0.56   NPNT2 = 10   NPNT3 =  8         *
!        TK = 800 K     EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
!        TK = 400 K     EBOT = -0.56   NPNT2 = 16   NPNT3 =  8         *
!        TK = 200 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
! * 
! *  There are two special cases determined by NPOL = 0 and NPOL < 0.  *
! *                                                                    *
! *  a) NPOL = 0 leads to density-of-states calculations               *
! *  with constant integration weights and equally distributed points  *
! *  between EBOT - pi*i*k*TK and EMU - pi*i*k*TK.                     *
! *                                                                    *
! *  The total number of integration points is given by:               *
! *              NPNT=NPNT2                                            *
! *                                                                    *
! *  b) NPOL < 0 is meant for calculations where the Fermi-Dirac       *
! *  function is replaced by a step function with step at EMU. When    *
! *  this option is used no poles of the Fermi-Dirac function are used *
! *  and the contour consists of the three straight lines:             *
! *                                                                    *
! *  1. the line from EBOT to EBOT-2*NPOL*pi*i*k*TK                    *
! *              with NPNT1 integration points (Gauss-Legendre rule)   *
! *                                                                    *
! *  2. the line from EBOT-2*NPOL*pi*i*k*TK to EMU-2*NPOL*pi*i*k*TK    *
! *              with NPNT2 integration points (Gauss-Legendre rule)   *
! *                                                                    *
! *  3. the line from EMU-2*NPOL*pi*i*k*TK to EMU                      *
! *              with NPNT3 integration points (Gauss-Legendre rule)   *
! *                                                                    *
! *                                                                    *
! **********************************************************************
  subroutine emesht(ez, df, npnt, ebot, emu, efermi, tk, npol, npnt1, npnt2, npnt3, iemxd)
    use GaussWeights_mod, only: gauleg => gauss_legendre_weights, gaufd => gauss_fermi_dirac_weights
    use Constants_mod, only: pi
    
    double precision, intent(in) :: ebot, emu, tk, efermi
    integer, intent(in) :: npnt1, npnt2, npnt3, npol, iemxd
    double complex, intent(out) :: df(*), ez(*)
    integer, intent(out) :: npnt

    double complex :: de
    double precision :: er, etk, febot, femu
    integer :: i, nleg, npnt4
    double precision, parameter :: kb=0.6333659d-5, ryd=13.6058d0!, pi=3.14159265358979312d0
    double precision :: wi(228), xi(228)
    
    integer :: ist
    logical, external :: opt
    
    character(len=*), parameter :: F9030="(13X,'points used with Gauss-Legendre rule: NPT2 =',I3,/,13X,'points near Fermi level: NPT3 =',I3,' and NPT4 =',I3)", &
    F9040="(13X,'points used with Gauss-Legendre rule: ',I1,'*NPT2 =',I2,/,13X,'points near Fermi level: NPT3 =',I3,' and NPT4 =',I3)", &
    F9000="(5X,'Density-of-States calculation',/,5X,'Number of energy points :',I4,4X,'broadening =',3P,F9.3,' ( mRy )',/,48X,' =',3P,F9.3,' ( meV )')", &
    F9090="(5X,'GF integration rectangular contour ( ImE ',a1,' 0 )',/,5X,'Number of energy points :',I4,13X,'poles =',I2,/,23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)", &
    F9010="(/,5X,'GF integration rectangular contour ( ImE < 0 )',/,5X,'Number of energy points :',I4,/,13X,'poles: NPOL  =',I3,/,13X,'points on upwards contour: NPT1 =',I3)", &
    F9020="(/,5X,'GF integration rectangular contour ( ImE < 0 )',/,5X,'Number of energy points :',I4,/,13X,'poles: NPOL  =',I3,3X,'NPOL2  =',I3,/,13X,'points on upwards contour: NPT1 =',I3)"
    
    nleg = 1
    npnt4 = 3
    if (tk < 2000.d0) npnt4 = 2

    write(6,'(5x,a,f10.6," (Ry)",8x,a,f10.6," (Ry)")')              'E min = ',ebot,'Fermi energy = ',efermi
    write(6,'(5x,a,f10.6," (Ry)",8x,a,f12.6," (K )",/,5x,62(1h-))') 'E max = ',emu, 'Temperature  = ',tk

    etk = pi*kb*tk
    selectcase(npol)
    case (0) ! density of states calculation
    
      open(87, file='FGRID', form='formatted', status='old', action='read', iostat=ist)
      if (ist == 0) then
        ! file FGRID for fine energy grid exists
        read(87, *) febot, femu
        close(87)

        de = femu - febot
        if (npnt2 > 1) then
          de = de/(npnt2-3)
        else
          de = dcmplx(1.d0, 0.d0)
        endif
        ez(1) = dcmplx(ebot, etk)
        df(1) = de
        npnt = 1
        do i = 2, npnt2-1
          npnt = npnt + 1
          assert(npnt <= iemxd)
          er = febot + (i-1)*de
          ez(npnt) = dcmplx(er, etk)
          df(npnt) = de
        enddo ! i
        ez(npnt+1) = dcmplx(emu, etk)
        df(npnt+1) = de
        npnt       = npnt + 1

      else ! ist == 0
      
        de = (emu - ebot)
        if (npnt2 > 1) then
          de = de/(npnt2-1)
        else
          de = dcmplx(1.d0, 0.d0)
        endif
        npnt = 0
        do i = 1, npnt2
          npnt = npnt + 1
          assert(npnt <= iemxd)
          er = ebot + (i-1)*de
          ez(npnt) = dcmplx(er, etk)
          df(npnt) = de
        enddo ! i
        
      endif ! ist == 0
      write(6, fmt=F9000) npnt, etk, etk*ryd
      
    case (1:) ! npol > 0)
! The following statements for NPOL > 0 are kept for backwards compatibility.
! They were used before the year 2009 in older versions of the 
! KKR-GF programs. (remark E.R.: kkrnano - this code is used if NPOL>0)
! *                                                                    *
! * The three lines in the backwards compatibility code are defined by:*
! *                                                                    *
! *  1. the line from EBOT to EBOT+2*NPOL*pi*i*k*TK                    *
! *              with NPNT1 integration points (Gauss-Legendre rule)   *
! *                                                                    *
! *  2. the line from EBOT+2*NPOL*pi*i*k*TK to                         *
! *                   EMU+(2*NPOL*pi*i-30)*k*TK                        *
! *              with NPNT2 integration points (Gauss-Legendre rule)   *
! *                                                                    *
! *  3. the line from EMU+(2*NPOL*pi*i-30)*k*TK to infinity            *
! *              with NPNT3 integration points (Gauss-Fermi-Dirac rule)*
! *                                                                    *
! *  The total number of integration points is given by:               *
! *              NPNT=NPNT1+NPNT2+NPNT3+NPOL                           *
! *                                                                    *
! *  The integration points and weights on three lines are chosen      *
! *  according to Gauss integration rules. Only in third interval      *
! *  the Fermi function matters since exp(x) < 10**(-10) for x < -25.  *
! *                                                                    *
      ist = gauleg(npnt1, xi, wi)
      de = npol*dcmplx(0.d0, etk)
      npnt = 0
      do i = 1, npnt1
        npnt = npnt + 1
        assert(npnt <= iemxd)
        ez(npnt) = xi(i)*de + de + ebot
        df(npnt) = wi(i)*de
      end do ! i
      ist = gauleg(npnt2, xi, wi)
      de = (emu-30*kb*tk-ebot)*0.5d0
      do i = 1, npnt2
        npnt = npnt + 1
        assert(npnt <= iemxd)
        ez(npnt) = xi(i)*de + de + ebot + 2*npol*dcmplx(0.d0, etk)
        df(npnt) = wi(i)*de
      enddo ! i
      ist = gaufd(npnt3, xi, wi)
      de = 30*kb*tk
      do i = 1, npnt3
        npnt = npnt + 1
        assert(npnt <= iemxd)
        ez(npnt) = xi(i)*de + emu + 2*npol*dcmplx(0.d0, etk)
        df(npnt) = wi(i)*de
      enddo ! i
      do i = npol, 1, -1
        npnt = npnt + 1
        assert(npnt <= iemxd)
        ez(npnt) = emu + (2*i-1)*dcmplx(0.d0, etk)
        df(npnt) = -2*dcmplx(0.d0, etk)
      enddo ! i
      write(6, fmt=F9090) '>', npnt, npol, npnt1, npnt2, npnt3
      
    case (:-1) ! npol < 0
    
      if (npnt1 > 0) ist = gauleg(npnt1, xi, wi)
      de = -npol*dcmplx(0.d0, etk)
      npnt = 0
      do i = 1, npnt1
        assert(npnt <= iemxd)
        npnt = npnt + 1
        ez(npnt) = xi(i)*de + de + ebot
        df(npnt) = wi(i)*de
      enddo ! i
      ist = gauleg(npnt2, xi, wi)
      de = (emu - ebot)*0.5d0
      do i = 1, npnt2
        npnt = npnt + 1
        assert(npnt <= iemxd)
        ez(npnt) = xi(i)*de + de + ebot - 2*npol*dcmplx(0.d0, etk)
!         if (opt('gf-ef   ')) ez(npnt) = emu + npol*dcmplx(0.d0,etk)
        df(npnt) = wi(i)*de
      enddo ! i
      if (npnt3 > 0) ist = gauleg(npnt3, xi, wi)
      de = -npol*dcmplx(0.d0, etk)
      do i = npnt3, 1, -1
        npnt = npnt + 1
        assert(npnt <= iemxd)
        ez(npnt) = xi(i)*de + de + emu
        df(npnt) = -wi(i)*de
      enddo ! i
      write(6, fmt=F9090) '<', npnt, -npol, npnt1, npnt2, npnt3
      
    endselect ! npol
    write(6, *) ! empty line
  endsubroutine ! emesht
  
  
  
  subroutine epathtb(ez, df, efermi, npnt, iesemicore, idosemicore, &
                    ebotval, emuval, tkval, npolval, n1val, n2val, n3val, &
                    ebotsem, emusem, tksem, npolsem, n1sem, n2sem, n3sem, &
                    iemxd)
! **********************************************************************
! * Generating the energy mesh.                                        *
! *                                                                    *
! * Calls the routine EMESHT once for the valence contour and once for *
! * the semicore contour.                                              *
! * In the semicore range, -NPOLSEM is used to create a rectangular    *
! * contour.                                                           *
! *              ph. mavropoulos, v.popescu Juelich/Munich 2004        *
! **********************************************************************
    integer, intent(in) :: iemxd
    double complex, intent(out) :: ez(*), df(*)
    double precision, intent(in) :: ebotsem,emusem,tksem,ebotval,emuval,tkval,efermi
    integer, intent(in) :: npolsem,n1sem,n2sem,n3sem
    integer, intent(in) :: npolval, n1val, n2val, n3val
    integer, intent(in) :: idosemicore
    integer, intent(out) :: iesemicore
    integer, intent(out) :: npnt
    
    double complex :: ezsemi(iemxd), dfsemi(iemxd), ezval(iemxd), dfval(iemxd)
    integer :: npntsemi,npntval,ie,je
    character(len=*), parameter :: F99001="(7x,'* ',a,/,7x,20(1h-),/)"

    write(6,'(/,79(1h=))')
    write(6,'(20x,a)') 'EPATHTB: generates a complex E contour'
    write(6,'(79(1h=),/)')

    iesemicore = 0 
    if (idosemicore == 1) then
      write(6, fmt=F99001) 'semi-core contour'
      call emesht(ezsemi, dfsemi, npntsemi, ebotsem, emusem, efermi, tksem, -npolsem, n1sem, n2sem, n3sem, iemxd)
      iesemicore = npntsemi
      write(6, fmt=F99001) 'valence contour'
    endif ! semi
    call emesht(ezval, dfval, npntval, ebotval, emuval, efermi, tkval, npolval, n1val, n2val, n3val, iemxd)

    npnt = iesemicore + npntval

    do ie = 1, iesemicore
      ez(ie) = ezsemi(ie)
      df(ie) = dfsemi(ie)
    enddo ! ie

    do ie = iesemicore+1, npnt
      je = ie - iesemicore
      ez(ie) = ezval(je)
      df(ie) = dfval(je)
    enddo ! ie

  endsubroutine ! epathtb
  
endmodule EnergyMeshHelpers_mod
