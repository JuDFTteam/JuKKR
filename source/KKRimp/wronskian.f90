!------------------------------------------------------------------------------------
!> Summary: Calculation of the Wronskian of the wavefunctions for numerical tests
!> Author: David Bauer
!> Calculation of the Wronskian of the wavefunctions for numerical tests in the 
!> calculation of the t-matrix. The wronskian is calculated for both regular and
!> irregular solutions.
!------------------------------------------------------------------------------------
module mod_wronskian

contains
  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the Wronskian of the wavefunctions for numerical tests
  !> Author: David Bauer
  !> Category: numerical-tools, sanity-check, KKRimp
  !> Deprecated: False
  !> Calculation of the Wronskian of the wavefunctions for numerical tests in the 
  !> calculation of the t-matrix. The wronskian is calculated for both regular and
  !> irregular solutions.
  !-------------------------------------------------------------------------------
  subroutine calcwronskian(rll,sll,leftrll,leftsll,cellnew)
    use mod_mathtools, only: matmat1T, matmatT1, matvec_dzdz
    use type_cellnew, only:  cell_typenew
    use mod_cheb, only: getCLambdaCinv
    implicit none
  
    double complex ::  rll(:,:,:)
    double complex ::  sll(:,:,:)
    double complex ::  leftrll(:,:,:)
    double complex ::  leftsll(:,:,:)
    type(cell_typenew) :: cellnew
  
    double complex, allocatable ::  drlldr(:,:,:)
    double complex, allocatable ::  dslldr(:,:,:)
    double complex, allocatable ::  dleftrlldr(:,:,:)
    double complex, allocatable ::  dleftslldr(:,:,:)
    double complex, allocatable ::  fn(:),dfndr(:)
    double complex, allocatable ::  fn2(:),dfn2dr(:)
    double complex, allocatable ::  fn3(:),dfn3dr(:)
    double complex, allocatable ::  fn4(:),dfn4dr(:)
    double complex, allocatable ::  wronskian(:,:,:)
    double complex, allocatable ::  wronskian2(:,:,:)
  
    double complex             ::  CLambdaCinv(0:cellnew%ncheb,0:cellnew%ncheb)
    double precision           ::  CLambdaCinv_temp(0:cellnew%ncheb,0:cellnew%ncheb)
    double precision          :: widthfac
    integer :: ilm1,ilm2,ir,ipan,irstart,irstop
    integer :: lmsize,lmsize2,nrmax
  
    lmsize  = ubound(rll,2)
    lmsize2 = ubound(rll,1)
    nrmax   = ubound(rll,3)
  
    write(*,*) lmsize,lmsize2,nrmax
    allocate( drlldr(lmsize2,lmsize,nrmax),dslldr(lmsize2,lmsize,nrmax) )
  
    allocate( dleftrlldr(lmsize2,lmsize,nrmax),dleftslldr(lmsize2,lmsize,nrmax) )
  
    allocate( wronskian(lmsize,lmsize,nrmax) )
    allocate( wronskian2(lmsize2,lmsize2,nrmax) )
  
    allocate(fn(nrmax),fn2(nrmax),dfndr(nrmax),dfn2dr(nrmax))
    allocate(fn3(nrmax),fn4(nrmax),dfn3dr(nrmax),dfn4dr(nrmax))
  
    ! if (lmsize/=lmsize2) stop '[calcwronskian] error lmsize/=lmsize2'
  
  
    ! ############################################################3
    ! Differentiate the wave functions
    ! ############################################################3
  
    call getCLambdaCinv(cellnew%Ncheb,CLambdaCinv_temp)
    CLambdaCinv=CLambdaCinv_temp
    write(2213,'(50000E)') CLambdaCinv_temp
    write(2214,'(50000E)') CLambdaCinv
  
  
    do ilm1=1,lmsize2
      do ilm2=1,lmsize
        fn =rll(ilm1,ilm2,:)
        fn2=sll(ilm1,ilm2,:)
        fn3=leftrll(ilm1,ilm2,:)
        fn4=leftsll(ilm1,ilm2,:)
        do ipan=1,cellnew%npan_tot
          irstart=cellnew%ipan_intervall(ipan-1)+1
          irstop = cellnew%ipan_intervall(ipan)
          widthfac = 2.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
  
          dfndr(irstart:irstop) = matvec_dzdz(CLambdaCinv,fn(irstart:irstop))
          dfndr(irstart:irstop) = dfndr(irstart:irstop)*widthfac
          dfn2dr(irstart:irstop) = matvec_dzdz(CLambdaCinv,fn2(irstart:irstop))
          dfn2dr(irstart:irstop) = dfn2dr(irstart:irstop)*widthfac
          dfn3dr(irstart:irstop) = matvec_dzdz(CLambdaCinv,fn3(irstart:irstop))
          dfn3dr(irstart:irstop) = dfn3dr(irstart:irstop)*widthfac
          dfn4dr(irstart:irstop) = matvec_dzdz(CLambdaCinv,fn4(irstart:irstop))
          dfn4dr(irstart:irstop) = dfn4dr(irstart:irstop)*widthfac
        end do
        drlldr(ilm1,ilm2,:)=dfndr
        dslldr(ilm1,ilm2,:)=dfn2dr
        dleftrlldr(ilm1,ilm2,:)=dfn3dr
        dleftslldr(ilm1,ilm2,:)=dfn4dr
      end do 
    end do
  
    ! if (lmsize2==2*lmsize) then
    !   allocate(temp1(lmsize,lmsize))
    ! !   print *,lmsize
    ! !   stop 'sdf'
    ! 
    !   do ir=1,nrmax
    !     temp1=dleftslldr(:lmsize,:,ir)
    !     dleftslldr(:lmsize,:,ir)=dleftslldr(lmsize+1:,:,ir) !-leftsll(lmsize+1:,:,ir) !/cellnew%rmeshnew(ir)
    !     dleftslldr(lmsize+1:,:,ir)=temp1 ! +leftsll(:lmsize,:,ir)/cellnew%rmeshnew(ir))*1.0D0
    !     temp1=drlldr(:lmsize,:,ir)
    !     drlldr(:lmsize,:,ir)=drlldr(lmsize+1:,:,ir) !-rll(lmsize+1:,:,ir)!/cellnew%rmeshnew(ir)
    !     drlldr(lmsize+1:,:,ir)=temp1 ! +rll(:lmsize,:,ir)/cellnew%rmeshnew(ir))*1.0D0
    !   end do 
    ! end if
  
  
  
    do ir=1,nrmax
    !   write(2212,'(50000E)') rll(:,:,ir)
    !   write(2212,'(50000E)') dslldr(:,:,ir)
    !   write(2212,'(50000E)') sll(:,:,ir)
    !   write(2212,'(50000E)') drlldr(:,:,ir)
    !   write(*,*) 'q',ubound(dleftslldr),ubound(rll)
      wronskian(:,:,ir)= matmatT1( dleftslldr(:,:,ir),rll(:,:,ir)  )
      wronskian2(:,:,ir)= matmat1T( rll(:,:,ir),leftsll(:,:,ir)  )
    !   write(2212,'(50000E)') wronskian(:,:,ir)
    !   write(*,*) 'p',ubound(leftsll),ubound(drlldr)
      wronskian(:,:,ir)= wronskian(:,:,ir) - matmatT1(leftsll(:,:,ir),drlldr(:,:,ir))
      wronskian2(:,:,ir)= wronskian2(:,:,ir) - matmat1T(sll(:,:,ir),leftrll(:,:,ir))
  
    !   write(2212,'(50000E)') wronskian(:,:,ir)
    !   stop
    end do 
  
    ! write(2222,'(50000E)') rll(1,1,:)
    ! write(2223,'(50000E)') drlldr(1,1,:)
    ! write(2224,'(50000E)') sll(1,1,:)*drlldr(1,1,:)-dslldr(1,1,:)*rll(1,1,:)
    ! write(2225,'(50000E)') sll(1,1,:)*rll(lmsize+1,1,:)-sll(lmsize+1,1,:)*rll(1,1,:)
    ! write(2226,'(50000E)') rll(lmsize+1,1,:)
    ! write(2227,'(50000E)') drlldr(1,1,:)
  
  
    open(unit=3246762,file='test_wronskian')
    open(unit=3246763,file='test_wronskian2')
    ! open(unit=3246764,file='test_wronskian_sra')
  
    do ilm1=1,lmsize
      do ilm2=1, lmsize
        write(3246762,'(50000E)') wronskian(ilm2,ilm1,:)
      end do
    end do
  
    do ilm1=1,lmsize2
      do ilm2=1, lmsize2
        write(3246763,'(50000E)') wronskian2(ilm2,ilm1,:)
      end do
    end do
  
  
    ! if (lmsize2==2*lmsize) then
    !   do ilm1=1,lmsize
    !     do ilm2=1, lmsize
    !       write(3246764,'(50000E)') sll(ilm1,ilm2,:)*rll(lmsize+ilm1,ilm2,:)-sll(lmsize+ilm1,ilm2,:)*rll(ilm1,ilm2,:)
    !     end do
    !   end do
    ! end if
  
    close(3246762)
    close(3246763)
    close(3246764)
  
  end subroutine calcwronskian

end module mod_wronskian
