!------------------------------------------------------------------------------------
!> Summary: Reads in LDA+U arrays from formatted file `ldaupot`
!> Author: 
!> Reads in LDA+U arrays from formatted file `ldaupot`
!------------------------------------------------------------------------------------
module mod_startldau
  use :: mod_datatypes, only: dp
  private :: dp

contains

   !-------------------------------------------------------------------------------
   !> Summary: Reads in LDA+U arrays from formatted file `ldaupot`
   !> Author:
   !> Category: lda+u, input-output, initialization, KKRhost
   !> Deprecated: False 
   !> Reads in LDA+U arrays from formatted file `ldaupot`
   !-------------------------------------------------------------------------------
  subroutine startldau(itrunldau,idoldau,kreadldau,lopt,ueff,jeff,erefldau,natyp,   &
    nspin,wldau,uldau,phildau,irws,ntldau,itldau,irmd,natypd,nspind,mmaxd)

    use :: mod_readldaupot
    use :: mod_rinit
    use :: mod_cinit
    implicit none
    ! ..
    integer :: irmd, mmaxd, natypd, nspind, irws(natypd)
    ! ..
    ! .. Arguments ..
    integer :: itrunldau, idoldau, kreadldau, natyp, nspin, ntldau
    integer :: lopt(natypd), itldau(natypd)
    real (kind=dp) :: ueff(natypd), jeff(natypd), erefldau(natypd)
    real (kind=dp) :: wldau(mmaxd, mmaxd, nspind, natypd)
    real (kind=dp) :: uldau(mmaxd, mmaxd, mmaxd, mmaxd, natypd)
    complex (kind=dp) :: phildau(irmd, natypd)
    ! ..
    ! .. Locals ..
    integer :: i1, im1, im3, is, it, ll
    ! ..
    ! ----------------------------------------------------------------------
    itrunldau = 0
    idoldau = 1
    ntldau = 0
    do it = 1, natyp
      if (lopt(it)+1/=0) then
        ntldau = ntldau + 1
        itldau(ntldau) = it
      end if
    end do

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    write (1337, '(79("="),/,27X,A,/, 79("="),/)') 'LDA+U: starting parameters'
    write (1337, 100) natyp, ntldau
    write (1337, 110)
    write (1337, 120)
    do it = 1, ntldau
      i1 = itldau(it)
      write (1337, 130) i1, ueff(i1), jeff(i1), erefldau(i1)
    end do
    write (1337, 120)
    write (1337, *)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ! -> read in LDA+U from file if available (KREADLDAU=1)
    call rinit(mmaxd*mmaxd*nspind*natypd, wldau)
    call cinit(irmd*natypd, phildau)
    if (kreadldau==1) then
      write (1337, 140)
      call readldaupot(itrunldau,lopt,ueff,jeff,erefldau,natyp,wldau,uldau,phildau, &
        irws,ntldau,itldau,irmd,natypd,nspind,mmaxd)
    else
      write (1337, 150)
    end if

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    if (itrunldau/=0) then
      write (1337, 160) 'Coulomb matrix U(m1,m1,m3,m3)'
      do it = 1, ntldau
        i1 = itldau(it)
        ll = lopt(i1)
        ll = min(3, ll)
        write (1337, 170) i1
        do im1 = 1, 2*ll + 1
          write (1337, 180)(uldau(im1,im1,im3,im3,i1), im3=1, 2*ll+1)
        end do
        write (1337, *)
        if (it<ntldau) write (1337, 190)
      end do
      write (1337, 160) 'Interaction potential W(m1,m2)'
      do it = 1, ntldau
        i1 = itldau(it)
        ll = lopt(i1)
        ll = min(3, ll)
        do is = 1, nspin
          write (1337, 200) i1, is
          do im1 = 1, 2*ll + 1
            write (1337, 180)(wldau(im1,im3,is,i1), im3=1, 2*ll+1)
          end do
          write (1337, *)
        end do
        if (it<ntldau) write (1337, 190)
      end do
      write (1337, '(9X,60("-"))')
    else
      call rinit(mmaxd*mmaxd*mmaxd*mmaxd*natypd, uldau)
      call rinit(mmaxd*mmaxd*nspind*natypd, wldau)
      call cinit(irmd*natypd, phildau)
    end if
    write (1337, *)

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

100 format (6x, 'Number of atoms ', '  in the u.c. :', i4, /, 24x, 'using LDA+U :', i4, /)
110 format (9x, ' IT ', '   Ueff   ', '   Jeff   ', '   Eref   ', ' (Ry)')
120 format (9x, 40('-'))
130 format (9x, i3, 1x, 3f10.6)
140 format (9x, 'Reading in LDA+U potential information (file ldaupot)', /)
150 format (9x, 'LDA+U potential initialised (set to zero)')
160 format (9x, 60('-'), /, 9x, a, /, 9x, 60('-'), /)
170 format (9x, 'IT =', i3, /)
180 format (9x, 7f10.6)
190 format (11x, 58('~'))
200 format (9x, 'IT =', i3, ' ISPIN =', i2, /)
  end subroutine startldau

end module mod_startldau
