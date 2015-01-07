  subroutine order_lines(kpt2irr, nkpts_all, kpt2irr_ord, band_indices)
    integer,         intent(in)   :: kpt2irr(nkpts_all)
    integer,         intent(in)   :: nkpts_all
    integer,         allocatable  :: band_indices(:), kpt2irr_ord(:)
    logical                       :: beginning_of_band
    integer                       :: i1, i2, i3, i_ord, i_band

    i_ord =1
    i_band=1    
    !loop over the k-points to find the beginning of bands
    do i1=1,nkpts_all
      beginning_of_band=.true.

      ! check if i1 is the beginning of a band
      do i2=i1+1,nkpts_all
        if(kpt2irr(i2)==kpt2irr(i1))then
          beginning_of_band=.false.
          exit
        end if!kpt2irr(i2)==kpt2irr(i1)
      end do!i2=i1+1,nkpts_all

      ! check if i1 was not already treated (it could be the end of a band)
      if (i_ord>1)then
        do i3=1, i_ord-1
          if(kpt2irr_ord(i2)==kpt2irr(i1))then
            beginning_of_band=.false.
            exit
          end if!kpt2irr_ord(i2)==kpt2irr(i1)
        end do!i3=1, i_ord-1
      end if!i_ord>1

      !save all kpts from the band
      if(beginning_of_band==.true.)then
        call traceback_band(i_ord,i_band,i1,kpt2irr,nkpts_all,kpt2irr_ord,band_indices)
        i_band=i_band+1
      end if!end_of_band=.true.

    end do!i1=1,nkpts_all    
    
  end subroutine

  recursive subroutine traceback_band(i_ord, i_band, i, kpt2irr, nkpts_all, kpt2irr_ord, band_indices)
    integer,         intent(in)   :: kpt2irr(nkpts_all)
    integer,         intent(in)   :: i_ord, i_band, i, nkpts_all
    integer                       :: band_indices(nkpts_all), kpt2irr_ord(nkpts_all)
    integer                       :: i_sameline, j          
    logical                       :: end_of_band

    ! This routine calls itself recursively until the end of the band is reached 

    ! save the 2 kpts of the current line 
    kpt2irr_ord(i_ord)=kpt2irr(i)
    band_indices(i_ord)=i_band
    if(MOD(i, 2).ne.0)then ! i odd
      i_sameline = i+1
    else ! i even
      i_sameline = i-1
    end if
    band_indices(i_ord)=i_band
    kpt2irr_ord(i_ord)=kpt2irr(i)
    band_indices(i_ord+1)=i_band
    kpt2irr_ord(i_ord+1)=kpt2irr(i_sameline)

    end_of_band=.true.
    ! check if i_sameline is the end of a band
    do j=1,nkpts_all
      if(kpt2irr(j)==kpt2irr(i_sameline) .and. j.ne.i_sameline)then
        end_of_band=.false.
        exit
      end if!kpt2irr(j)==kpt2irr(i_sameline) .and. j.ne.i_sameline
    end do!j=1,nkpts_all

    if(end_of_band==.false.)then
      call traceback_band(i_ord+2, i_band, j, kpt2irr, nkpts_all, kpt2irr_ord, band_indices) 
    end if!end_of_band==.false. 
 
  end subroutine
