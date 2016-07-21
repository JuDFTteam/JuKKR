!*==wrmoms.f    processed by SPAG 6.05Rc at 15:39 on  1 Mar 2002
!     Write angular momentum resolved charges to stdout
  subroutine wrmoms(nspin, charge, i1, lmax, lmaxp1, first, last)
    implicit none
    integer, intent(in) :: nspin, i1, lmax, lmaxp1
    logical, intent(in) :: first, last
    double precision, intent(in) :: charge(0:lmaxp1,2)
    ! local variables
    character(len=5), parameter :: textns = ' ns =' ! non-spherical
    character(len=4), parameter :: textell(0:9) = [' s =',' p =',' d =',' f =',' g =',' h =',' i =',' j =',' k =',' l =']
    
    double precision :: muspin(0:lmaxp1+1), sumch(2)
    integer :: l

    muspin(:) = 0.d0
    sumch(:) = 0.d0
    
    if (nspin == 2) then
    
      if (first) then ! first atom
        write(6, '(/,78(1h#))')
        write(6, fmt='(15x,a)') "l-decomposed valence charges and magnetic moments"
        write(6, '(78(1h#))')
        write(6, fmt='(/,8x,a)') " ATOM        N_el spin dn  N_el spin up    m_spin"
        write(6, '(9x,49(1h=))')
      endif ! first atom

      do l = 0, lmaxp1
        sumch(1:2) = sumch(1:2) + charge(l,1:2)
        muspin(l) = charge(l,2) - charge(l,1)
        muspin(lmaxp1+1) = muspin(lmaxp1+1) + muspin(l)
      enddo ! l

      l = 0
        write(6, fmt='(i13,2x,a4,2(f14.8),2x,f8.4)') i1, textell(l), charge(l,1:2), muspin(l)
      do l = 1, lmax
        write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') textell(l), charge(l,1:2), muspin(l)
      enddo ! l
      l = lmaxp1
        write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') textns, charge(l,1:2), muspin(l)

      write(6, fmt='(16x,19(1h-),$)')
      write(6,'(23(1h-))')
      write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') ' TOT', sumch(1:2), muspin(lmaxp1+1)
      write(6,'(33x,f14.8)') sumch(1) + sumch(2)
      
      if (last) then ! last atom
        write(6,'(/,78(1h#))')
        write(6, *)
      else  ! last atom
        write(6, '(9x,49(1h=))') ! separation line
      endif ! last atom
    
    else

      if (first) then ! first atom
        write(6, '(/,44(1h#))')
        write(6, fmt='(8x,a)') "l-decomposed valence charges"
        write(6, '(44(1h#))')
        write(6, fmt='(/,8x," ATOM        N_el ")')
        write(6, '(9x,26(1h=))')
      endif ! first atom

      do l = 0, lmaxp1
        sumch(1) = sumch(1) + charge(l,1)
      enddo ! l

        write(6, fmt='(i13,2x,a4,f14.8)') i1, textell(0), charge(0,1)
      do l = 1, lmax
        write(6, fmt='(15x,a4,f14.8)') textell(l), charge(l,1)
      enddo ! l
        write(6, fmt='(15x,a4,f14.8)') textns, charge(lmaxp1,1)

      write(6, fmt='(16x,19(1h-))')
      write(6, fmt='(15x,a4,f14.8)') ' TOT', sumch(1)

      if (last) then ! last atom
        write(6,'(/,44(1h#))')
        write(6, *)
      else  ! last atom
        write(6, '(9x,26(1h=))') ! separation line
      endif ! last atom
    
    endif
    
  endsubroutine ! wrmoms

#if 0
  subroutine wrmoms(nspin, charge, i1, lmax, lmaxp1, first, last)
    implicit none
    integer, intent(in) :: nspin, i1, lmax, lmaxp1
    logical, intent(in) :: first, last
    double precision, intent(in) :: charge(0:lmaxp1,2)
    ! local variables
    character(len=5), parameter :: textns = ' ns =' ! non-spherical
    character(len=4), parameter :: textell(0:9) = [' s =',' p =',' d =',' f =',' g =',' h =',' i =',' j =',' k =',' l =']
    character(len=7), parameter :: textspin(0:2) = ['       ','spin dn','spin up']
    
    double precision :: muspin(0:lmaxp1+1), sumch(2)
    integer :: ispin, l

    if (first) then ! first atom

      write(6, *)
      if (nspin == 2) then
        write(6, '(78(1h#))')
        write(6, fmt='(15x,a)') "l-decomposed valence charges and magnetic moments"
        write(6, '(78(1h#))')
      else
        write(6, '(44(1h#))')
        write(6, fmt='(8x,a)') "l-decomposed valence charges"
        write(6, '(44(1h#))')
      endif
      
      write(6, *)

      write(6, fmt='(8x," ATOM      ",$)')
      do ispin = 1, nspin
        write(6, fmt='("  N_el ",a7,$)') textspin(ispin + nspin - 2)
      enddo ! ispin
      if (nspin == 2) write(6, fmt='("    m_spin",$)')
      write(6, *) ! finish the line

      if (nspin == 2) then 
        write(6, '(9x,49(1h=))')
      else
        write(6, '(9x,26(1h=))')
      endif

    endif ! first atom

    muspin(0:) = 0.d0
    sumch(1:2) = 0.d0
    do l = 0, lmaxp1
      do ispin = 1, nspin
        sumch(ispin) = sumch(ispin) + charge(l,ispin)
      enddo ! ispin
      muspin(l) = charge(l,2) - charge(l,1)
      muspin(lmaxp1+1) = muspin(lmaxp1+1) + muspin(l)
    enddo ! l

    do l = 0, 0
      if (nspin == 2) then
        write(6, fmt='(8x,i5,2x,a4,2(f14.8),2x,f8.4)') i1, textell(l), charge(l,1:nspin), muspin(l)
      else
        write(6, fmt='(10x,i3,2x,a4,f14.8)') i1, textell(l), charge(l,1)
      endif
    enddo ! l

    do l = 1, lmax
      if (nspin == 2) then
        write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') textell(l), charge(l,1:nspin), muspin(l)
      else
        write(6, fmt='(15x,a4,f14.8)') textell(l), charge(l,1)
      endif
    enddo ! l

    do l = lmaxp1, lmaxp1
      if (nspin == 2) then
        write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') textns, charge(l,1:nspin), muspin(l)
      else
        write(6, fmt='(15x,a4,f14.8)') textns, charge(l,1)
      endif
    enddo ! l

    write(6, fmt='(16x,19(1h-),$)')
    if (nspin == 2) then
      write(6,'(23(1h-))')
      write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') ' TOT', sumch(1:nspin), muspin(lmaxp1+1)
      write(6,'(33x,f14.8)') sumch(1) + sumch(2)
    else
      write(6, *)
      write(6, fmt='(15x,a4,f14.8)') ' TOT',sumch(1)
    endif

    if (.not. last) then ! not the last atom
      write(6, '(9x,26(1h=),$)')
      if (nspin == 2) write(6, '(23(1h=),$)')
    else
      ! last atom
      write(6, *)
      if (nspin == 2) then
        write(6,'(78(1h#))')
      else
        write(6,'(44(1h#))')
      endif
    endif ! last atom
    write(6, *)
    
  endsubroutine ! wrmoms

  subroutine wrmoms_nspin1(charge, i1, lmax, lmaxp1, first, last)
    implicit none
    integer, intent(in) :: i1, lmax, lmaxp1
    logical, intent(in) :: first, last
    double precision, intent(in) :: charge(0:)
    ! local variables
    character(len=5), parameter :: textns = ' ns =' ! non-spherical
    character(len=4), parameter :: textell(0:9) = [' s =',' p =',' d =',' f =',' g =',' h =',' i =',' j =',' k =',' l =']
    
    double precision :: sumch
    integer :: l

    if (first) then ! first atom
      write(6, *)
      write(6, '(44(1h#))')
      write(6, fmt='(8x,a)') "l-decomposed valence charges"
      write(6, '(44(1h#))')
      write(6, *)
      write(6, fmt='(8x," ATOM        N_el ")')
      write(6, '(9x,26(1h=))')
    endif ! first atom

    sumch = 0.d0
    do l = 0, lmaxp1
      sumch = sumch + charge(l)
    enddo ! l

      write(6, fmt='(i13,2x,a4,f14.8)') i1, textell(0), charge(0)
    do l = 1, lmax
      write(6, fmt='(15x,a4,f14.8)') textell(l), charge(l)
    enddo ! l
      write(6, fmt='(15x,a4,f14.8)') textns, charge(lmaxp1)

    write(6, fmt='(16x,19(1h-))')
    write(6, fmt='(15x,a4,f14.8)') ' TOT', sumch

    if (last) then ! last atom
      write(6, *)
      write(6,'(44(1h#))')
      write(6, *)
    else  ! last atom
      write(6, '(9x,26(1h=))') ! separation line
    endif ! last atom
    
  endsubroutine ! wrmoms
  
  subroutine wrmoms_nspin2(charge, i1, lmax, lmaxp1, first, last)
    implicit none
    integer, intent(in) :: i1, lmax, lmaxp1
    logical, intent(in) :: first, last
    double precision, intent(in) :: charge(0:lmaxp1,2)
    ! local variables
    character(len=5), parameter :: textns = ' ns =' ! non-spherical
    character(len=4), parameter :: textell(0:9) = [' s =',' p =',' d =',' f =',' g =',' h =',' i =',' j =',' k =',' l =']
    
    double precision :: muspin(0:lmaxp1+1), sumch(2)
    integer :: l

    if (first) then ! first atom
      write(6, *)
      write(6, '(78(1h#))')
      write(6, fmt='(15x,a)') "l-decomposed valence charges and magnetic moments"
      write(6, '(78(1h#))')
      write(6, *)
      write(6, fmt='(8x,a)') " ATOM        N_el spin dn  N_el spin up    m_spin"
      write(6, '(9x,49(1h=))')
    endif ! first atom

    muspin(0:) = 0.d0
    sumch(1:2) = 0.d0
    do l = 0, lmaxp1
      sumch(1:2) = sumch(1:2) + charge(l,1:2)
      muspin(l) = charge(l,2) - charge(l,1)
      muspin(lmaxp1+1) = muspin(lmaxp1+1) + muspin(l)
    enddo ! l

    l = 0
      write(6, fmt='(8x,i5,2x,a4,2(f14.8),2x,f8.4)') i1, textell(l), charge(l,1:2), muspin(l)
    do l = 1, lmax
      write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') textell(l), charge(l,1:2), muspin(l)
    enddo ! l
    l = lmaxp1
      write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') textns, charge(l,1:2), muspin(l)

    write(6, fmt='(16x,19(1h-),$)')
    write(6,'(23(1h-))')
    write(6, fmt='(15x,a4,2(f14.8),2x,f8.4)') ' TOT', sumch(1:2), muspin(lmaxp1+1)
    write(6,'(33x,f14.8)') sumch(1) + sumch(2)
    
    if (last) then ! last atom
      write(6, *)
      write(6,'(78(1h#))')
      write(6, *)
    else  ! last atom
      write(6, '(9x,49(1h=))') ! separation line
    endif ! last atom
    
  endsubroutine ! wrmoms
#endif
  