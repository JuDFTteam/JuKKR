  subroutine analytic_continuation(nw,nweb,nwef,onsite,struct)
! analytic continuation from line parallel to real axis
! structure of the data:
! nw + 1 energy panels
! first panel contains the original nescf points
! then next nw panels contain 2*(nescf+nweb+nwef) points
! the range is E-wmax to E+wmax, points in each panel are ordered like this
! for Eb: (Eb-w(i),Eb+w(i)), w(i) = wmax/nw
! for Ef: (Ef-w(i),Ef+w(i)), w(i) = wmax/nw
! the panels with nescf data don't enter the calculation here
  use global

  implicit none

! structure of the energy panels
  integer(kind=i4b), intent(in) :: nw, nweb, nwef
  logical,           intent(in) :: onsite, struct
! ----------------------------------------------------------------------
  complex(kind=c8b) :: gfold(nlmsb,nlmsb), gfnew(nlmsb,nlmsb), gf(nlmsb,nlmsb)
  complex(kind=c8b) :: de, eim, eold, enew
  integer(kind=i4b) :: ia, ia2, ja, ja2, iew, iew0, jw, ie


! ******************************
! Loop over blocks of the GF
  do ja2=1,nasusc2      ! atom j
    ja = iasusc2(ja2)
    do ia2=1,nasusc2    ! atom i
      ia = iasusc2(ia2)
! ******************************
!     Panels around bottom of the contour
      if (nweb > 0) then
      end if
!     Panels around Fermi energy
      if (nwef > 0) then
!       GF at Ef
        iew = nescf
        eold = esusc(iew)
        eim = cmplx(0.d0,-aimag(eold))
        write(*,'("eold,eim=",4es16.8)') eold, eim
!       Panels for E + w
        call projected_gf(iew,ia,ja,gfold,onsite,struct)
        if (lrot) call local_frame(ia,ja,magdir,gfold,gf)
!       Loop over subpanels
        do jw=1,nw
!         middle of the panel
          iew0 = nescf + 2*(jw-1)*(nescf+nweb+nwef) + 2*(nescf+nweb) + nwef
          do ie=1,nwef
            iew = iew0 + ie
            enew = esusc(iew)
            de = enew - eold
            write(*,'("enew,de=",i6,4es16.8)') iew, enew, de
!           GF at E + w
            call projected_gf(iew,ia,ja,gfnew,onsite,struct)
            if (lrot) call local_frame(ia,ja,magdir,gfnew,gf)
!           correction
            gf = gfnew
            gfnew = gfnew + (eim/de)*(gf-gfold)
            eold = enew
            gfold = gf
          end do
        end do
!       Panels for E - w
        call projected_gf(iew,ia,ja,gfold,onsite,struct)
        if (lrot) call local_frame(ia,ja,magdir,gfold,gf)
!       Loop over subpanels
        do jw=1,nw
!         middle of the panel
          iew0 = nescf + 2*(jw-1)*(nescf+nweb+nwef) + 2*(nescf+nweb) + nwef
          do ie=1,nwef
            iew = iew0 - ie
            enew = esusc(iew)
            de = enew - eold
!           GF at E + w
            call projected_gf(iew,ia,ja,gfnew,onsite,struct)
            if (lrot) call local_frame(ia,ja,magdir,gfnew,gf)
!           correction
            gf = gfnew
            gfnew = gfnew + (eim/de)*(gf-gfold)
            eold = enew
            gfold = gf
          end do
        end do
      end if
! ********
    end do
  end do
! ********
! All done!
  end subroutine analytic_continuation
