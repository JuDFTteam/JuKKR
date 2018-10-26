!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_cmomsread
  
  private
  public :: cmomsread

contains

  !-------------------------------------------------------------------------------
  !> Summary: Read CMOMS for host for decimation
  !> Author: 
  !> Date: 29.10.99
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine reads in the CMOMHOST from the decimation files  
  !> Note that they are needed only in case of an SCF decimation 
  !> calculation (SCFSTEPS > 1 )                                         
  !>                                                                  
  !> The t-matrices are writen out in kloopz1 (option 'deci-out')    
  !>                                                                  
  !> This subroutine must be called after the t-matrices for all the  
  !> energies are read in (see < decimaread > )                       
  !> It returns the CMOMHOST array. First the left bulk (unit 37) then
  !> the right bulk (unit 38) are indexed.                            
  !> CMOMHOST(*,NEMBD1) =                                             
  !>                CMOMHOST(*,1..NLBASIS,NLBASIS+1..NLBASIS+NRBASIS) 
  !> Condider this mapping for further use.                           
  !>                                                                  
  !>                                                  29.10.99        
  !>                                                  05.06.04
  !-------------------------------------------------------------------------------
  subroutine cmomsread(nlbasis, nrbasis, naez, cmomhost, vacflag, kaoez, natypd, nembd1, lmpotd)
    use :: mod_datatypes, only: dp
    implicit none
    ! ..
    ! .. Arguments
    integer :: lmpotd, natypd, nembd1
    integer :: naez, nlbasis, nrbasis
    real (kind=dp) :: cmomhost(lmpotd, nembd1)
    integer :: kaoez(natypd, *)
    logical :: vacflag(2)
    ! ..
    ! .. Local variables ..
    real (kind=dp) :: c00(lmpotd)
    character (len=5) :: chhost(2)
    integer :: ih, ih1, ihl, ihost, lm, lmpotl, naezl, nathost
    ! ..
    ! .. Data statements
    data chhost/'LEFT ', 'RIGHT'/
    ! ..
    write (1337, '(5X,A,/,8X,30("-"),/,8X,3A6,A10,/,8X,30("-"))') 'Reading in host charge moments ( SCFSTEPS > 1 )', ' HOST ', '  IBAS', '  ATOM', '   CMOM(1)'
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
    do ihost = 1, 2
      nathost = nlbasis
      if (ihost==2) nathost = nrbasis
      write (1337, '(8X,A5,1X)', advance='no') chhost(ihost)
      ! ----------------------------------------------------------------------
      if (vacflag(ihost)) then
        do ih = 1, nlbasis
          do lm = 1, lmpotd
            cmomhost(lm, (ihost-1)*nlbasis+ih) = 0.e0_dp
            ! mapping the CMOMHOST array, ordering is important
          end do
        end do
        write (1337, '(A)') ' Vacuum setting    0.000'
        if (ihost==1) then
          write (1337, '(14X,24("-"))')
        else
          write (1337, '(8X,30("-"))')
        end if
        ! ----------------------------------------------------------------------
      else
        read (36+ihost, 110) naezl, lmpotl
        ! ......................................................................
        if (naezl/=nathost) then
          write (6, '(/,5X,2A)') 'ERROR: ', 'host not compatible with your input.'
          write (6, '(/,12X,A,I3,A,I3)') 'Charge moments tabulated for', naezl, ' host atoms, input NBASIS =', nathost
          stop '       < CMOMSREAD > '
        end if
        ! ......................................................................
        do ih = 1, naezl
          read (36+ihost, *) ihl
          if (ihl/=ih) then
            write (6, '(/,5X,2A,/)') 'ERROR reading host file', ' basis indexing wrong'
            stop '       < CMOMSREAD > '
          end if
          read (36+ihost, 100)(c00(lm), lm=1, lmpotl)
          ih1 = kaoez(1, naez+(ihost-1)*nlbasis+ih)

          do lm = 1, lmpotl
            cmomhost(lm, (ihost-1)*nlbasis+ih) = c00(lm)
            ! mapping the CMOMHOST array, ordering is important
          end do

          if (ih==1) then
            write (1337, '(1X,2I6,D12.4)') ih, ih1, c00(1)
          else
            write (1337, '(14X,2I6,D12.4)') ih, ih1, c00(1)
          end if
        end do
        ! ......................................................................
        if (ihost==1) then
          write (1337, '(14X,24("-"))')
        else
          write (1337, '(8X,30("-"))')
        end if
      end if
      ! ----------------------------------------------------------------------
    end do
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
    write (1337, *)

100 format (4d22.14)
110 format (5x, 2i6)
  end subroutine cmomsread

end module mod_cmomsread
