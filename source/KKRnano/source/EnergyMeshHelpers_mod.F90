!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module EnergyMeshHelpers_mod
!-------------------------------------------------------------------------------
!> Summary: Point mesh helpers for the energy contour integration
!> Author: Phivos Mavropoulos, Voicu Popescu, Elias Rabel, Marcel Bornemann, Paul F Baumeister
!> Category: KKRnano, input-output, initialization
!-------------------------------------------------------------------------------
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: load, store, update, broadcast
  public :: epathtb

  interface load
    module procedure readEnergyMesh
  endinterface

  interface store
    module procedure writeEnergyMesh
  endinterface

  interface update
    module procedure updateEnergyMesh
  endinterface

  interface broadcast
    module procedure broadcastEnergyMesh_com
  endinterface

  contains

  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file
  subroutine readEnergyMesh(e1, e2, efermi, ez, ielast, npnt123, npol, &
               tk, wez, ebotsemi, emusemi, fsemicore, iesemicore, npntsemi, kmesh, filename)
    ! valence contour parameters
    double precision, intent(out) :: e1, e2, efermi
    double complex, intent(out) :: ez(:)
    integer, intent(out) :: ielast
    integer, intent(out) :: npnt123(3)
    integer, intent(out) :: npol
    double precision, intent(out) :: tk
    double complex, intent(out) :: wez(:)

    ! semicore contour parameters
    double precision, intent(out) :: ebotsemi, emusemi, fsemicore
    integer, intent(out) :: iesemicore
    integer, intent(out) :: npntsemi(3)
    integer, intent(out) :: kmesh(:) !< dim(iemxd)
    character(len=*), intent(in) :: filename

    open(67, file=filename, form='unformatted', action='read', status='old')
    read(67) ielast,ez,wez,e1,e2
    read(67) npol,tk,npnt123
    read(67) efermi
    read(67) iesemicore,fsemicore,ebotsemi
    read(67) emusemi
    read(67) npntsemi
    read(67) kmesh
    close(67)

  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMesh(e1, e2, efermi, ez, ielast, npnt123, npol, &
         tk, wez, ebotsemi, emusemi, fsemicore, iesemicore, npntsemi, kmesh, filename)
    ! valence contour parameters
    double precision, intent(in) :: e1, e2, efermi
    double complex, intent(in) :: ez(:)
    integer, intent(in) :: ielast
    integer, intent(in) :: npnt123(3)
    integer, intent(in) :: npol
    double precision, intent(in) :: tk
    double complex, intent(in) :: wez(:)

    ! semicore contour parameters
    double precision, intent(in) :: ebotsemi, emusemi, fsemicore
    integer, intent(in) :: iesemicore
    integer, intent(in) :: npntsemi(3)
    integer, intent(in) :: kmesh(:) !< dim(iemxd)
    character(len=*), intent(in) :: filename

    open(67, file=filename, form='unformatted', action='write')
    write(67) ielast,ez,wez,e1,e2
    write(67) npol,tk,npnt123
    write(67) efermi
    write(67) iesemicore,fsemicore,ebotsemi
    write(67) emusemi
    write(67) npntsemi
    write(67) kmesh
    close(67)

  endsubroutine ! write

  !------------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EPATHTB
  subroutine updateEnergyMesh(ez, wez, ielast, e1, e2, tk, npol, npnt123v, &
                        ebotsemi, emusemi, iesemicore, fsemicore, npntsemi)
    use Constants_mod, only: pi
    double complex, intent(out) :: ez(:), wez(:)
    integer, intent(inout) :: ielast
    double precision, intent(in) :: e1, e2, tk
    integer, intent(in) :: npol, npnt123v(3) !< valence contour parameters
    
    ! semicore contour parameters
    double precision, intent(in) :: ebotsemi, emusemi, fsemicore
    integer, intent(in) :: npntsemi(3)
    integer, intent(out) :: iesemicore !< points #1 through #iesemicore belong to the semicore contour

    integer :: iemxd

    iemxd = ielast
    ! --> update energy contour for both, semicore contour (if present) and valence
    call epathtb(ez, wez, e2, ielast, iesemicore, e1, e2, tk, npol, npnt123v(1), npnt123v(2), npnt123v(3), &
                                       ebotsemi, emusemi, tk, npol, npntsemi(1), npntsemi(2), npntsemi(3), iemxd)

    wez(1:ielast) = -2.d0/pi*wez(1:ielast)
    wez(1:iesemicore) = wez(1:iesemicore)*fsemicore ! scale the weights for the semicore contour only

  endsubroutine ! update

  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from rank 'BCRANK' to all other ranks
  subroutine broadcastEnergyMesh_com(comm, e1, e2, e3, e4, iemxd, ez, wez)
    integer, intent(in) :: comm
    double precision, intent(inout) :: e1, e2 !< valence contour parameters
    double precision, intent(inout) :: e3, e4 !< semicore contour parameters
    integer, intent(in) :: iemxd
    double complex, intent(inout) :: ez(:), wez(:)
#ifndef NoMPI
    include 'mpif.h'
    integer :: ierr
    double precision :: es(4)
    es = [e1, e2, e3, e4]
    call MPI_Bcast(ez,  iemxd, MPI_DOUBLE_COMPLEX, 0, comm, ierr)
    call MPI_Bcast(wez, iemxd, MPI_DOUBLE_COMPLEX, 0, comm, ierr)
    call MPI_Bcast(es,    4, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    e1 = es(1); e2 = es(2); e3 = es(3); e4 = es(4)
#endif
  endsubroutine ! broadcast

  
  subroutine epathtb(ez, df, efermi, ieboth, iesemicore, &
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
    double complex, intent(out) :: ez(iemxd), df(iemxd)
    double precision, intent(in) :: efermi
    double precision, intent(in) :: ebotsem, emusem, tksem
    integer, intent(in) :: npolsem, n1sem, n2sem, n3sem
    double precision, intent(in) :: ebotval, emuval, tkval
    integer, intent(in) :: npolval, n1val, n2val, n3val
    integer, intent(out) :: iesemicore, ieboth
    
    integer :: ievalence
    character(len=*), parameter :: F99001="(7x,'* ',a,/,7x,20(1h-),/)"

!    assert(max(0, npolval) + n1val + n2val + n3val + n1sem + n2sem + n3sem <= iemxd)
    
    if (.false.) then
      write(6,'(/,79(1h=))')
      write(6,'(20x,a)') 'EPATHTB: generates a complex E contour'
      write(6,'(79(1h=),/)')
    endif

    iesemicore = 0 
    if (n1sem + n2sem + n3sem > 0) then
      write(6, fmt=F99001) 'semi-core contour'
      call emesht(ez, df, iesemicore, ebotsem, emusem, efermi, tksem, -npolsem, n1sem, n2sem, n3sem, iemxd)
      write(6, fmt=F99001) 'valence contour'
    endif ! semi
    call emesht(ez(iesemicore+1:), df(iesemicore+1:), ievalence, ebotval, emuval, efermi, tkval, npolval, n1val, n2val, n3val, iemxd-iesemicore)

    ieboth = iesemicore + ievalence

  endsubroutine ! epathtb

  
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
    npnt4 = 3; if (tk < 2000.d0) npnt4 = 2

    write(6,'(5x,a,f10.6," (Ry)",8x,a,f10.6," (Ry)")')              'E min = ',ebot,'Fermi energy = ',efermi
    write(6,'(5x,a,f10.6," (Ry)",8x,a,f12.6," (K )",/,5x,62(1h-))') 'E max = ',emu, 'Temperature  = ',tk

    etk = pi*kb*tk
    selectcase(npol)
    case (0) ! density of states calculation
    
      assert(iemxd >= npnt2)
    
      open(87, file='FGRID', form='formatted', status='old', action='read', iostat=ist)
      if (ist == 0) then
        ! file FGRID for fine energy grid exists
        read(87, *) febot, femu
        close(87)

        de = dcmplx(1.d0, 0.d0); if (npnt2 > 1) de = (femu - febot)/(npnt2 - 3.)
        
        ez(1) = dcmplx(ebot, etk)
        df(1) = de
        npnt = 1
        do i = 2, npnt2-1
          npnt = npnt + 1
          er = febot + (i-1)*dreal(de)
          ez(npnt) = dcmplx(er, etk)
          df(npnt) = de
        enddo ! i
        npnt = npnt + 1
        ez(npnt) = dcmplx(emu, etk)
        df(npnt) = de

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
          er = ebot + (i-1)*dreal(de)
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

      assert(iemxd >= npnt1 + npnt2 + npnt3 + npol)

      ist = gauleg(npnt1, xi, wi)
      de = npol*dcmplx(0.d0, etk)
      npnt = 0
      do i = 1, npnt1
        npnt = npnt + 1
        ez(npnt) = xi(i)*de + de + ebot
        df(npnt) = wi(i)*de
      enddo ! i
      ist = gauleg(npnt2, xi, wi)
      de = (emu-30*kb*tk-ebot)*0.5d0
      do i = 1, npnt2
        npnt = npnt + 1
        ez(npnt) = xi(i)*de + de + ebot + 2*npol*dcmplx(0.d0, etk)
        df(npnt) = wi(i)*de
      enddo ! i
      ist = gaufd(npnt3, xi, wi)
      de = 30*kb*tk
      do i = 1, npnt3
        npnt = npnt + 1
        ez(npnt) = xi(i)*de + emu + 2*npol*dcmplx(0.d0, etk)
        df(npnt) = wi(i)*de
      enddo ! i
      do i = npol, 1, -1
        npnt = npnt + 1
        ez(npnt) = emu + (2*i-1)*dcmplx(0.d0, etk)
        df(npnt) = -2*dcmplx(0.d0, etk)
      enddo ! i
      write(6, fmt=F9090) '>', npnt, npol, npnt1, npnt2, npnt3
      
    case (:-1) ! npol < 0
    
      assert(iemxd >= npnt1 + npnt2 + npnt3)
      
      if (npnt1 > 0) ist = gauleg(npnt1, xi, wi)
      de = -npol*dcmplx(0.d0, etk)
      npnt = 0
      do i = 1, npnt1
        npnt = npnt + 1
        ez(npnt) = xi(i)*de + de + ebot
        df(npnt) = wi(i)*de
      enddo ! i
      ist = gauleg(npnt2, xi, wi)
      de = (emu - ebot)*0.5d0
      do i = 1, npnt2
        npnt = npnt + 1
        ez(npnt) = xi(i)*de + de + ebot - 2*npol*dcmplx(0.d0, etk)
!         if (opt('gf-ef   ')) ez(npnt) = emu + npol*dcmplx(0.d0,etk)
        df(npnt) = wi(i)*de
      enddo ! i
      if (npnt3 > 0) ist = gauleg(npnt3, xi, wi)
      de = -npol*dcmplx(0.d0, etk)
      do i = npnt3, 1, -1
        npnt = npnt + 1
        ez(npnt) =  xi(i)*de + de + emu
        df(npnt) = -wi(i)*de
      enddo ! i
      write(6, fmt=F9090) '<', npnt, -npol, npnt1, npnt2, npnt3
      
    endselect ! npol
    write(6, *) ! empty line
  endsubroutine ! emesht
  
endmodule ! EnergyMeshHelpers_mod
