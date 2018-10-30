module mod_calctmat_bauernew_testtools

contains

subroutine switch_jlk(jlk)
implicit none
double complex :: jlk(:,:),temp
integer :: dim1,dim2,dim3,ir,lm

dim1=ubound(jlk,1)
dim2=ubound(jlk,2)
dim3=dim1/2
write(*,*) dim1,dim2

do ir=1,dim2
  write(345,'(5000E)') jlk(:,ir)
  do lm=1,dim3
!     print *, lm
!     print *, lm+dim3
    temp=jlk(lm+dim3,ir)
!     print *, lm+dim3,lm
    jlk(lm+dim3,ir) = jlk(lm,ir)
!     print *, lm
    jlk(lm,ir) = temp
  end do !lm
  write(345,'(5000E)') jlk(:,ir)
end do !dim2

end subroutine switch_jlk


subroutine switch_vll(vll)
implicit none
double complex :: vll(:,:,:),temp
integer :: dim1,dim2,dim3,dim1h,dim2h,ir,lm1,lm2

dim1=ubound(vll,1)
dim2=ubound(vll,2)
dim3=ubound(vll,3)
write(*,*) dim1,dim2,dim3
dim1h=dim1/2
dim2h=dim2/2

do ir=1,dim3
  do lm1=1,dim1h
    do lm2=1,dim2h

    temp=vll(lm2+dim2h,lm1+dim1h,ir)
    vll(lm2+dim2h,lm1+dim1h,ir) = vll(lm2,lm1,ir)
    vll(lm2,lm1,ir) = temp

    temp=vll(lm2+dim2h,lm1,ir)
    vll(lm2+dim2h,lm1,ir) = vll(lm2,lm1+dim1h,ir)
    vll(lm2,lm1+dim1h,ir) = temp

    end do !lm
  end do !lm
end do !dim2

end subroutine switch_vll







end module mod_calctmat_bauernew_testtools