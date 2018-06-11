    Subroutine calc_orbitalmoment(lmax, lmsize, loperator)

      Use constants
      Use profiling
      Use mod_datatypes, Only: dp

      Implicit None
      Integer, Intent (In) :: lmax, lmsize
      Complex (Kind=dp), Dimension (lmsize, lmsize, 3), &
        Intent (Out) :: loperator
      Integer, Save :: first = 1
      Integer :: lval
      Complex (Kind=dp), Dimension (:, :, :), Allocatable :: lorbit_onel
      Integer :: lmmax, lstart, lstop, i_stat, i_all

      lmmax = (lmax+1)**2
      loperator = czero

      loperator(1, 1, 1) = czero
      loperator(1, 1, 2) = czero
      loperator(1, 1, 3) = czero

      Do lval = 1, lmax
        Allocate (lorbit_onel(2*lval+1,2*lval+1,3), Stat=i_stat)
        Call memocc(i_stat, product(shape(lorbit_onel))*kind(lorbit_onel), &
          'lorbit_onel', 'calc_orbitalmoment')
        lorbit_onel = czero
        Call calc_orbit_onel(lval, lorbit_onel)
        lstart = ((lval-1)+1)**2 + 1
        lstop = ((lval)+1)**2
        loperator(lstart:lstop, lstart:lstop, 1) = lorbit_onel(:, :, 1)
        loperator(lstart:lstop, lstart:lstop, 2) = lorbit_onel(:, :, 2)
        loperator(lstart:lstop, lstart:lstop, 3) = lorbit_onel(:, :, 3)
!
        i_all = -product(shape(lorbit_onel))*kind(lorbit_onel)
        Deallocate (lorbit_onel, Stat=i_stat)
        Call memocc(i_stat, i_all, 'lorbit_onel', 'calc_orbitalmoment')
      End Do
      If (lmsize/=lmmax) Then
        loperator(lmmax+1:lmsize, lmmax+1:lmsize, :) = loperator(:lmmax, &
          :lmmax, :)
      End If

! if (first==1) then
!  open(unit=423492157,file='out_Lx')
!  open(unit=423492158,file='out_Ly')
!  open(unit=423492159,file='out_Lz')
!  do ilm=1,lmsize
!    write(423492157,'(5000F)'),Loperator(ilm,:,1)
!    write(423492158,'(5000F)'),Loperator(ilm,:,2)
!    write(423492159,'(5000F)'),Loperator(ilm,:,3)
!  end do
!  close(423492157);close(423492158);close(423492159)
! end if

      first = 0
    End Subroutine


    Subroutine calc_orbit_onel(lval, lorbit_onel)

      Use constants
      Use mod_datatypes, Only: dp

      Implicit None
      Integer, Intent (In) :: lval
      Complex (Kind=dp), Dimension (2*lval+1, 2*lval+1, 3), &
        Intent (Out) :: lorbit_onel
      Complex (Kind=dp), Dimension (-lval:lval, -lval:lval) :: l_z
      Complex (Kind=dp), Dimension (-lval:lval, -lval:lval) :: l_x
      Complex (Kind=dp), Dimension (-lval:lval, -lval:lval) :: l_y
      Complex (Kind=dp), Dimension (-lval:lval, -lval:lval) :: l_dn
      Complex (Kind=dp), Dimension (-lval:lval, -lval:lval) :: l_up
      Integer :: i1
      Real (Kind=dp) :: lfac
!
      l_z = czero
      l_x = czero
      l_y = czero
      l_dn = czero
      l_up = czero
!
      lorbit_onel = czero
      Do i1 = -lval, lval
        l_z(-i1, i1) = -ci*i1
      End Do
!
      If (lval>0) Then
!
        lfac = sqrt(lval*(lval+1E0_dp))/sqrt(2E0_dp)
        l_dn(0, -1) = -ci*lfac
        l_dn(0, 1) = lfac
        l_dn(-1, 0) = ci*lfac
        l_dn(1, 0) = -lfac
!
        If (lval>1) Then
          Do i1 = 2, lval
            lfac = 0.5E0_dp*sqrt(lval*(lval+1E0_dp)-i1*(i1-1E0_dp))
            l_dn(-i1, -i1+1) = -lfac
            l_dn(-i1, i1-1) = ci*lfac
            l_dn(i1, -i1+1) = -ci*lfac
            l_dn(i1, i1-1) = -lfac
!
            lfac = 0.5E0_dp*sqrt(lval*(lval+1E0_dp)-(i1-1E0_dp)*i1)
            l_dn(-i1+1, -i1) = lfac
            l_dn(-i1+1, i1) = ci*lfac
            l_dn(i1-1, -i1) = -ci*lfac
            l_dn(i1-1, i1) = lfac
          End Do
        End If
      End If
!
      If (lval>0) Then
!
        lfac = sqrt(lval*(lval+1E0_dp))/sqrt(2E0_dp)
        l_up(0, -1) = -ci*lfac
        l_up(0, 1) = -lfac
        l_up(-1, 0) = ci*lfac
        l_up(1, 0) = lfac
!
        If (lval>1) Then
          Do i1 = 2, lval
            lfac = 0.5E0_dp*sqrt(lval*(lval+1E0_dp)-i1*(i1-1E0_dp))
            l_up(-i1, -i1+1) = lfac
            l_up(-i1, i1-1) = ci*lfac
            l_up(i1, -i1+1) = -ci*lfac
            l_up(i1, i1-1) = lfac
!
            lfac = 0.5E0_dp*sqrt(lval*(lval+1E0_dp)-(i1-1E0_dp)*i1)
            l_up(-i1+1, -i1) = -lfac
            l_up(-i1+1, i1) = ci*lfac
            l_up(i1-1, -i1) = -ci*lfac
            l_up(i1-1, i1) = -lfac
          End Do
        End If
      End If
!
      l_x = 0.5E0_dp*(l_up+l_dn)
      l_y = -0.5E0_dp*ci*(l_up-l_dn)
!
      lorbit_onel(:, :, 1) = l_x
      lorbit_onel(:, :, 2) = l_y
      lorbit_onel(:, :, 3) = l_z

    End Subroutine
