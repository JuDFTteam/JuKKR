!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_deciopt

  private
  public :: deciopt
  
contains

  !-------------------------------------------------------------------------------
  !> Summary: Read-in driver for decifiles
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This routine treats the DECIMATION case setting up the single-site
  !> (Delta t)^(-1) matrices and the charge moments of the host(s).    
  !>                                                                   
  !> This is realised in two ways:                                     
  !>      - either reading in the matrices (and moments - if SCFSTEPS  
  !>        is greater than 1) as written out in a previous (bulk) run 
  !>        -- DECIFILES token points to the files containing the nece-
  !>        ssary information                                          
  !>      - or reading in the self-consistent potential for each host  
  !>        and effectively calculating the matrices; the potential    
  !>        must have the specific format set in < OUTPOTHOST > routine
  !>        and the DECIPOTS token should point to the corresponding   
  !>        potential file(s)                                          
  !>
  !> Notes:                                                            
  !>        - DECIFILES token is considered by default and sought first
  !>                          is requiring the same energy mesh for the
  !>                          system as for the host                   
  !>        - DECIPOTS token is not restrictive in this sense          
  !>                         however, it does not calculate charge mo- 
  !>                         ments -- does not work in SCF mode        
  !>                         is not dealing with CPA host so far       
  !>                         is not dealing with NON-SPHERICAL poten-  
  !>                         tials so far                              
  !>
  !>                                     V. Popescu - Munich, Dec 04
  !-------------------------------------------------------------------------------
  subroutine deciopt(alat, ins, krel, kvrel, kmrot, nspin, naez, lmmax, bravais, tk, npol, npnt1, npnt2, npnt3, ez, ielast, kaoez, lefttinvll, righttinvll, vacflag, nlbasis, &
    nrbasis, cmomhost, vref, rmtref, nref, refpot, lmaxd, lmgf0d, lmmaxd, lm2d, nembd1, iemxd, nspind, lmpotd, natypd, irmd, ipand)

    use :: mod_datatypes, only: dp
    use :: mod_decitset, only: decitset
    use :: mod_decimaread, only: decimaread
    use :: mod_cmomsread, only: cmomsread
    use :: mod_ioinput, only: ioinput
    implicit none
    ! ..
    ! .. Scalar arguments
    integer :: lmmaxd, nembd1, iemxd, nspind, lmpotd, natypd, ipand, irmd, lmaxd
    integer :: lm2d, nref, lmgf0d
    integer :: ins, krel, kmrot, nspin, naez, lmmax, npol, npnt1, npnt2, npnt3
    integer :: ielast, nlbasis, nrbasis, kvrel
    real (kind=dp) :: alat, tk
    ! ..
    ! .. Array arguments
    integer :: kaoez(natypd, *), refpot(nembd1)
    real (kind=dp) :: bravais(3, 3), cmomhost(lmpotd, *)
    real (kind=dp) :: vref(*), rmtref(*)
    complex (kind=dp) :: lefttinvll(lmmaxd, lmmaxd, nembd1, nspind, iemxd), righttinvll(lmmaxd, lmmaxd, nembd1, nspind, iemxd)
    complex (kind=dp) :: ez(iemxd)
    logical :: vacflag(2)
    ! ..
    ! .. Local scalars
    integer :: ierror, il, ie, ispin, nspinso ! ruess: for tmat newsolver
    complex (kind=dp) :: cfctor
    character (len=40) :: fileleft, fileright
    character (len=256) :: uio                             ! NCOLIO=256

    ! ..                                  ! ruess: for NEWSOSOL running option
    ! .. External Functions ..
    logical, external :: opt

    ! ======================================================================
    write (1337, '(79("="))')
    write (1337, '(15X,A,/,79("="),/)') 'DECIOPT: reading left/right host decimation files'
    il = 1
    ierror = 0
    call ioinput('DECIFILES       ', uio, il, 7, ierror)
    ! :::::::::::::::::::::::::::::::::::::::::::::::: decifiles (tmatrices)
    if (ierror==0) then
      read (unit=uio, fmt='(A40)') fileleft
      call ioinput('DECIFILES       ', uio, il+1, 7, ierror)
      read (unit=uio, fmt='(A40)') fileright
      ! ----------------------------------------------------------------------

      ! --> first call to read the header ( IE = 0 )

      ie = 0
      call decimaread(ez, tk, npnt1, npnt2, npnt3, npol, nspin, lefttinvll(1,1,1,1,1), righttinvll(1,1,1,1,1), vacflag, ie, nlbasis, nrbasis, naez, kaoez, kmrot, ins, nspin, lmmax, &
        ielast, fileleft, fileright, krel, natypd, lmmaxd, nembd1)

      ! --> get the left and right host Delta_t matrices

      cfctor = alat/(8.e0_dp*atan(1.0e0_dp)) ! = ALAT/(2*PI)
      nspinso = nspin
      if (opt('NEWSOSOL')) nspinso = 1 ! ruess: only combined l-s index for
      ! newsolver
      do ispin = 1, nspinso
        do ie = 1, ielast
          call decimaread(ez, tk, npnt1, npnt2, npnt3, npol, ispin, lefttinvll(1,1,1,ispin,ie), righttinvll(1,1,1,ispin,ie), vacflag, ie, nlbasis, nrbasis, naez, kaoez, kmrot, ins, &
            nspin, lmmax, ielast, fileleft, fileright, krel, natypd, lmmaxd, nembd1)

          ! --> host matrices have been written out in true units
          ! they are used in p.u. units (see kloopz) --> convert them here

          call zscal(lmmaxd*lmmaxd*nembd1, cfctor, lefttinvll(1,1,1,ispin,ie), 1)
          call zscal(lmmaxd*lmmaxd*nembd1, cfctor, righttinvll(1,1,1,ispin,ie), 1)
        end do
      end do

      ! --> get the left and right host charge moments
      ! ( not needed in single-iteration mode calculations )

      ! fivos        IF ( SCFSTEPS.GT.1 )
      call cmomsread(nlbasis, nrbasis, naez, cmomhost, vacflag, kaoez, natypd, nembd1, lmpotd)
      close (37)
      close (38)
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    else
      ! :::::::::::::::::::::::::::::::::::::::::::::::: decipots (calc tmats)
      ierror = 0
      call ioinput('DECIPOTS        ', uio, il, 7, ierror)
      if (ierror/=0) then
        write (6, 100)
        stop
      end if
      read (unit=uio, fmt='(A40)') fileleft
      call ioinput('DECIPOTS        ', uio, il+1, 7, ierror)
      read (unit=uio, fmt='(A40)') fileright
      call decitset(alat, bravais, ez, ielast, nlbasis, nrbasis, fileleft, fileright, ins, kvrel, krel, nspin, kmrot, vref, rmtref, nref, refpot, lefttinvll, righttinvll, vacflag, &
        nembd1, iemxd, irmd, ipand, lmaxd, lmgf0d, lmmaxd, lm2d, nspind)
    end if
    ! ======================================================================

100 format (/, 6x, 'ERROR : Missing decimation files (t-mat or pot)', /, 14x, 'Please use one of the tokens DECIFILES/DECIPOTS', /)
  end subroutine deciopt

end module mod_deciopt
