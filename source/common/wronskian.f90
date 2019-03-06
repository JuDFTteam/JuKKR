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
  !> irregular solutions and saved in the files `test_wronskian` and `test_wronskian2`
  !> The Wronskian relations are defined in 4.37 and 4.45 of the PhD thesis by 
  !> David Bauer (pp. 48).
  !-------------------------------------------------------------------------------
  subroutine calcwronskian(rll, sll, leftrll, leftsll, ncheb, npan_tot, ipan_intervall, rpan_intervall)
    
    use mod_datatypes, only: dp
    use mod_mathtools, only: matmat1T, matmatT1, matvec_dzdz
    use mod_cheb, only: getCLambdaCinv
    implicit none
    ! interface
    complex (kind=dp), intent(in) ::  rll(:,:,:) !! regular right wavefunction
    complex (kind=dp), intent(in) ::  sll(:,:,:) !! irregular right wavefunction
    complex (kind=dp), intent(in) ::  leftrll(:,:,:) !! regular left wavefunction
    complex (kind=dp), intent(in) ::  leftsll(:,:,:) !! irregular left wavefunction
    integer, intent(in) :: ncheb    !! number of Chebychev polynomials
    integer, intent(in) :: npan_tot !! number of panels
    integer, intent(in) :: ipan_intervall(:) !! index array for potential and shapefunction boundaries for each panel
    real (kind=dp), intent(in) :: rpan_intervall(:) !! radial values of panel boundaries 
    ! local arrays
    complex (kind=dp), allocatable ::  drlldr(:,:,:)     !! derivative of rll
    complex (kind=dp), allocatable ::  dslldr(:,:,:)     !! derivative of sll
    complex (kind=dp), allocatable ::  dleftrlldr(:,:,:) !! derivative of leftrll
    complex (kind=dp), allocatable ::  dleftslldr(:,:,:) !! derivative of leftsll
    complex (kind=dp), allocatable ::  fn(:),dfndr(:)    !! working arrays for radial derivatives
    complex (kind=dp), allocatable ::  fn2(:),dfn2dr(:)  !! working arrays for radial derivatives
    complex (kind=dp), allocatable ::  fn3(:),dfn3dr(:)  !! working arrays for radial derivatives
    complex (kind=dp), allocatable ::  fn4(:),dfn4dr(:)  !! working arrays for radial derivatives
    complex (kind=dp), allocatable ::  wronskian(:,:,:)  !! wronskian of the first kind (eq. 4.44, PhD Bauer)
    complex (kind=dp), allocatable ::  wronskian2(:,:,:) !! wronskian of the second kind (eq. 4.37, PhD Bauer)
    complex (kind=dp) :: CLambdaCinv(0:ncheb,0:ncheb) !! complex version of integration matrix in Chebychev mesh (imaginary part is zero)
    real (kind=dp) :: CLambdaCinv_temp(0:ncheb,0:ncheb) !! integration matrix in Chebychev mesh
    real (kind=dp)          :: widthfac !! width of panel = 1/2*(r(m)-r(m-1))
    ! working indices and matrix sizes
    integer :: ilm1,ilm2,ir,ipan,irstart,irstop
    integer :: lmsize,lmsize2,nrmax
  
    lmsize  = ubound(rll,2)
    lmsize2 = ubound(rll,1)
    nrmax   = ubound(rll,3)
  
    write(*,*) 'in calcwronskian:', lmsize,lmsize2,nrmax
    allocate( drlldr(lmsize2,lmsize,nrmax),dslldr(lmsize2,lmsize,nrmax) )
    allocate( dleftrlldr(lmsize2,lmsize,nrmax),dleftslldr(lmsize2,lmsize,nrmax) )
    allocate( wronskian(lmsize,lmsize,nrmax) )
    allocate( wronskian2(lmsize2,lmsize2,nrmax) )
    allocate( fn(nrmax),fn2(nrmax),dfndr(nrmax),dfn2dr(nrmax) )
    allocate( fn3(nrmax),fn4(nrmax),dfn3dr(nrmax),dfn4dr(nrmax) )
  
  
    ! ############################################################
    ! Differentiate the wave functions
    ! ############################################################

    ! construct Chebychev intetgration matrix
    call getCLambdaCinv(Ncheb, CLambdaCinv_temp(0:ncheb,0:ncheb))
    CLambdaCinv(0:ncheb,0:ncheb) = cmplx(CLambdaCinv_temp(0:ncheb,0:ncheb), 0.0_dp)
  
    ! perform radial derivative of rll, rllleft, sll and sllleft
    do ilm1=1,lmsize2
      do ilm2=1,lmsize
        fn =rll(ilm1,ilm2,:)
        fn2=sll(ilm1,ilm2,:)
        fn3=leftrll(ilm1,ilm2,:)
        fn4=leftsll(ilm1,ilm2,:)
        
        do ipan=1,npan_tot
          irstart=ipan_intervall(ipan-1)+1
          irstop = ipan_intervall(ipan)
          widthfac = 2.0D0/(rpan_intervall(ipan)-rpan_intervall(ipan-1))
  
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
  
    ! ############################################################
    ! Construct Wronskian
    ! ############################################################
    do ir=1,nrmax
      wronskian(:,:,ir)= matmatT1( dleftslldr(:,:,ir),rll(:,:,ir)  )
      wronskian2(:,:,ir)= matmat1T( rll(:,:,ir),leftsll(:,:,ir)  )
      wronskian(:,:,ir)= wronskian(:,:,ir) - matmatT1(leftsll(:,:,ir),drlldr(:,:,ir))
      wronskian2(:,:,ir)= wronskian2(:,:,ir) - matmat1T(sll(:,:,ir),leftrll(:,:,ir))
    end do 
  
    ! open output file for wronskian relations and write files out
    ! the files contail the left-hand site of eq. 4.44 and 4.37
    open(unit=3246762,file='test_wronskian')
    open(unit=3246763,file='test_wronskian2')
    write(3246762,'(a)') '# lm1, lm2, wronskian(lm1, lm2, 1:ir_max) (first kind)'
    write(3246763,'(a)') '# lm1, lm2, wronskian2(lm1, lm2, 1:ir_max) (second kind)'
  
    do ilm1=1, lmsize
      do ilm2=1, lmsize
        write(3246762,'(2i5,50000E)') ilm2, ilm1, wronskian(ilm2,ilm1,:)
      end do
    end do
  
    do ilm1=1, lmsize2
      do ilm2=1, lmsize2
        write(3246763,'(2i5,50000E)') ilm2, ilm1, wronskian2(ilm2,ilm1,:)
      end do
    end do
  
    close(3246762)
    close(3246763)
  
  end subroutine calcwronskian

end module mod_wronskian
