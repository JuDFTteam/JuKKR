  subroutine symmetrize_tgf()
! Enforces symmetry of collinear ASA SRA input
  use global

  implicit none

  complex(kind=c8b) :: tspdf(0:nlmax,nsmax)
  integer(kind=i4b) :: i2(2), ie, ialms, jalms, ia, ja, ilms, jlms, ilm, jlm, im, il, is, js


! **************
  do ie=1,nesusc
! **************

! ----------------------------------------------------------------------
! Symmetrize t-matrix
  do ia=1,nasusc
    tspdf = 0.d0
    do ilms=1,nlms
      i2 = i2lms(:,ilms)
      ilm = i2(1); is = i2(2)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      tspdf(il,is) = tspdf(il,is) + tmatrix(ilms,ilms,ia,ie)
    end do
    do is=1,nsmax
      do il=0,nlmax
        tspdf(il,is) = tspdf(il,is)/(2.d0*il+1.d0)
      end do
    end do
    tmatrix(:,:,ia,ie) = 0.d0
    do ilms=1,nlms
      i2 = i2lms(:,ilms)
      ilm = i2(1); is = i2(2)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      tmatrix(ilms,ilms,ia,ie) = tspdf(il,is)
    end do
  end do
! ----------------------------------------------------------------------
! Symmetrize structural GF
  do jalms=1,nalms
    i2 = i2alms(:,jalms)
    jlms = i2(1); ja = i2(2)
    i2 = i2lms(:,jlms)
    jlm = i2(1); js = i2(2)
    do ialms=1,nalms
      i2 = i2alms(:,ialms)
      ilms = i2(1); ia = i2(2)
      i2 = i2lms(:,jlms)
      ilm = i2(1); is = i2(2)
      if (is /= js) then
        gstruct(ialms,jalms,ie) = 0.d0
      else
        gstruct(ialms,jalms,ie) = 0.5d0*(gstruct(ialms,jalms,ie) + gstruct(jalms,ialms,ie))
        gstruct(jalms,ialms,ie) = gstruct(ialms,jalms,ie)
      end if
    end do
  end do
! ----------------------------------------------------------------------

! ******
  end do
! ******
  end subroutine symmetrize_tgf
