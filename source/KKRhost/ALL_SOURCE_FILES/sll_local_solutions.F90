module mod_sll_local_solutions

contains

subroutine sll_local_solutions(vll,tau,drpan2,csrc1,slc1sum, &
                         mihvy,mihvz,mijvy,mijvz, &
                         yif,zif, &
                         ncheb,ipan,lmsize,lmsize2,nrmax, &
                         nvec,jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor, &
                         cmodesll,LBESSEL,use_sratrick1)
      Use mod_datatypes, Only: dp
implicit none
      integer :: ncheb                               ! number of chebyshev nodes
      integer :: lmsize                              ! lm-components * nspin
      integer :: lmsize2                             ! lmsize * nvec
      integer :: nvec                                ! spinor integer
! nvec=1 non-rel, nvec=2 for sra and dirac
      integer :: nrmax                               ! total number of rad. mesh points

      integer :: LBESSEL, use_sratrick1      !  dimensions etc., needed only for host code interface


      complex (kind=dp),parameter:: cone=(1.0d0,0.0d0),czero=(0.0d0,0.0d0)

! running indices
      integer ivec, ivec2                            
      integer l1,l2,lm1,lm2,lm3
      integer info,icheb2,icheb,ipan,mn,nplm

! source terms
      complex (kind=dp) :: gmatprefactor               ! prefactor of green function
! non-rel: = kappa = sqrt e

      complex (kind=dp) :: HLK(LBESSEL,NRMAX), &
                        JLK(LBESSEL,NRMAX), &
                        HLK2(LBESSEL,NRMAX), &
                        JLK2(LBESSEL,NRMAX) 


      INTEGER JLK_INDEX(2*LMSIZE)


      character(len=1) :: cmodesll
! cmodesll ="1" : op( )=identity       for reg. solution
! cmodesll ="T" : op( )=transpose in L for reg. solution

      complex (kind=dp) :: vll(lmsize*nvec,lmsize*nvec,nrmax) ! potential term in 5.7

      complex (kind=dp) ::  &
                     mihvy(lmsize,lmsize),mihvz(lmsize,lmsize), &
                     mijvy(lmsize,lmsize),mijvz(lmsize,lmsize), &
                     yif(lmsize2,lmsize,0:ncheb), &       
                     zif(lmsize2,lmsize,0:ncheb)       
      complex (kind=dp) ::  &
                     srv(0:ncheb,lmsize2,0:ncheb,lmsize2), &
                     srv1(0:ncheb,lmsize,0:ncheb,lmsize), &
                     yill1(0:ncheb,lmsize,lmsize), zill1(0:ncheb,lmsize,lmsize), &
                     yill2(0:ncheb,lmsize,lmsize), zill2(0:ncheb,lmsize,lmsize), &
                     yill(0:ncheb,lmsize2,lmsize), zill(0:ncheb,lmsize2,lmsize), &
                     vjli(lmsize,lmsize2,0:ncheb), vhli(lmsize,lmsize2,0:ncheb), &
                     vjli_yill1(lmsize,lmsize), vhli_yill1(lmsize,lmsize), &
                     vjli_zill1(lmsize,lmsize), vhli_zill1(lmsize,lmsize), &
                     yill1temp(lmsize,lmsize), zill1temp(lmsize,lmsize)

      complex (kind=dp) ::  &
                     jlmkmn(0:ncheb,lmsize2,0:ncheb), &
                     hlmkmn(0:ncheb,lmsize2,0:ncheb)

! chebyshev arrays
      complex (kind=dp) zslc1sum(0:ncheb)
      real (kind=dp) drpan2
      real (kind=dp) &
                       csrc1(0:ncheb,0:ncheb), & ! Integration matrix from right ( C*S_R*C^-1 in eq. 5.54)
                       tau(0:ncheb), &    ! Radial mesh points
                       slc1sum(0:ncheb),taucsrcr,tau_icheb
      complex (kind=dp) :: gf_tau_icheb

      integer ipiv(0:ncheb,lmsize2)
      integer :: use_sratrick

if ( lmsize==1 ) then
  use_sratrick=0
else
  use_sratrick=use_sratrick1
end if

! initialization
  
  vhli=czero
  vjli=czero

  if (use_sratrick==0) then
    yill=czero
    zill=czero
  else
    yill1=czero
    zill1=czero
    yill2=czero
    zill2=czero
  end if

!---------------------------------------------------------------------
! 1. prepare VJLR, VNL, VHLR, which appear in the integrands
! TAU(K,IPAN) is used instead of TAU(K,IPAN)**2, which directly gives
! RLL(r) and SLL(r) multiplied with r. TAU is the radial mesh.
!
! 2. prepare the source terms YR, ZR, YI, ZI
! because of the conventions used by
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
! a factor sqrt(E) is included in the source terms
! this factor is removed by the definition of ZSLC1SUM given below
!
!vjlr = \kappa * J * V = \kappa * r * j *V
!vhlr = \kappa * H * V = \kappa * r * h *V
!
! i.e. prepare terms kappa*J*DV, kappa*H*DV appearing in 5.11, 5.12.

  do icheb = 0,ncheb
    mn = ipan*ncheb + ipan - icheb
    tau_icheb = tau(icheb)
    gf_tau_icheb = gmatprefactor*tau_icheb
    if     (cmodesll=='1') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjli(lm1,lm2+lmsize*(ivec2-1),icheb) = vjli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gf_tau_icheb*jlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
              vhli(lm1,lm2+lmsize*(ivec2-1),icheb) = vhli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gf_tau_icheb*hlk2(l1,mn)*vll(lm1+lmsize*(ivec-1),lm2+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmodesll=='T') then
      do ivec2=1,nvec
        do lm2 = 1,lmsize
          do ivec=1,nvec
            do lm1 = 1,lmsize
              l1 = jlk_index( lm1+lmsize*(ivec-1) )
              vjli(lm1,lm2+lmsize*(ivec2-1),icheb) = vjli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gf_tau_icheb*jlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
              vhli(lm1,lm2+lmsize*(ivec2-1),icheb) = vhli(lm1,lm2+lmsize*(ivec2-1),icheb) + &
                  gf_tau_icheb*hlk2(l1,mn)*vll(lm2+lmsize*(ivec-1),lm1+lmsize*(ivec2-1),mn)
            end do
          end do
        end do
      end do !nvec
    elseif (cmodesll=='0') then
              vjli(:,:,icheb) = czero
              vhli(:,:,icheb) = czero
    else
      stop '[rllsll] mode not known'
    end if

! calculation of the J (and H) matrix according to equation 5.69 (2nd eq.)
    if ( use_sratrick==0 ) then
      do ivec=1,nvec ! index for large/small component
        do lm1 = 1,lmsize
          l1 = jlk_index( lm1+lmsize*(ivec-1) )
          yill(icheb,lm1+lmsize*(ivec-1),lm1) =  tau_icheb*hlk(l1,mn)
          zill(icheb,lm1+lmsize*(ivec-1),lm1) =  tau_icheb*jlk(l1,mn)
        end do
      end do !ivec=1,nvec
    elseif ( use_sratrick==1 ) then
      do lm1 = 1,lmsize
        l1 = jlk_index( lm1 )
        l2 = jlk_index( lm1+lmsize )
        yill1(icheb,lm1,lm1) =  tau_icheb*hlk(l1,mn)
        zill1(icheb,lm1,lm1) =  tau_icheb*jlk(l1,mn)
        yill2(icheb,lm1,lm1) =  tau_icheb*hlk(l2,mn)
        zill2(icheb,lm1,lm1) =  tau_icheb*jlk(l2,mn)
      end do
    end if
  end do ! icheb

! calculation of A in 5.68
  if ( use_sratrick==0 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
         taucsrcr = tau(icheb)*csrc1(icheb,icheb2)*drpan2
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize2
          do ivec=1,nvec
            do lm3 = 1,lmsize
              lm1=lm3+(ivec-1)*lmsize
              l1 = jlk_index(lm1)
              srv(icheb,lm1,icheb2,lm2) = &
              taucsrcr*(-jlk(l1,mn)*vhli(lm3,lm2,icheb2) &
                        +hlk(l1,mn)*vjli(lm3,lm2,icheb2))
            end do
          end do
        end do
      end do
    end do
    do lm1 = 1,lmsize2
      do icheb = 0,ncheb
        srv(icheb,lm1,icheb,lm1) = srv(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do
  elseif  ( use_sratrick==1 ) then
    do icheb2 = 0,ncheb
      do icheb = 0,ncheb
         taucsrcr = tau(icheb)*csrc1(icheb,icheb2)*drpan2
         mn = ipan*ncheb + ipan - icheb
         do lm1 = 1,lmsize
            l1 = jlk_index(lm1)
            jlmkmn(icheb,lm1,icheb2) = taucsrcr*jlk(l1,mn)
            hlmkmn(icheb,lm1,icheb2) = taucsrcr*hlk(l1,mn)
         end do
      end do
    end do

        do lm2 = 1,lmsize
    do icheb2 = 0,ncheb
      call svpart(srv1(0,1,icheb2,lm2), &
                  jlmkmn(0,1,icheb2),hlmkmn(0,1,icheb2), &
                  vhli(1,lm2,icheb2),vjli(1,lm2,icheb2), &
                  ncheb,lmsize)
      end do
    end do
    do lm1 = 1,lmsize
      do icheb = 0,ncheb
        srv1(icheb,lm1,icheb,lm1) = srv1(icheb,lm1,icheb,lm1) + 1.d0
      end do
    end do

  else
    stop '[rllsll] error in inversion'
  end if

!-------------------------------------------------------
! determine the local solutions
! solve the equations SLV*YRLL=S and SLV*ZRLL=C
!                 and SRV*YILL=C and SRV*ZILL=S
! i.e., solve system A*U=J, see eq. 5.68.

  if ( use_sratrick==0 ) then
    nplm = (ncheb+1)*lmsize2

    if (cmodesll/='0') then
      call zgetrf(nplm,nplm,srv,nplm,ipiv,info)
      if (info/=0) stop 'rllsll: zgetrf'
      call zgetrs('n',nplm,lmsize,srv,nplm,ipiv,yill,nplm,info)
      call zgetrs('n',nplm,lmsize,srv,nplm,ipiv,zill,nplm,info)
    end if
  elseif ( use_sratrick==1 ) then
    nplm = (ncheb+1)*lmsize

    call zgetrf(nplm,nplm,srv1,nplm,ipiv,info)
    call zgetrs('n',nplm,lmsize,srv1,nplm,ipiv,yill1,nplm,info)
    call zgetrs('n',nplm,lmsize,srv1,nplm,ipiv,zill1,nplm,info)

    do icheb2 = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yill1temp(lm1,lm2) = yill1(icheb2,lm1,lm2)
          zill1temp(lm1,lm2) = zill1(icheb2,lm1,lm2)
        end do
      end do
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhli(1,1,icheb2), &
        lmsize,yill1temp,lmsize,czero,vhli_yill1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vhli(1,1,icheb2), &
        lmsize,zill1temp,lmsize,czero,vhli_zill1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjli(1,1,icheb2), &
        lmsize,yill1temp,lmsize,czero,vjli_yill1,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize,cone,vjli(1,1,icheb2), &
        lmsize,zill1temp,lmsize,czero,vjli_zill1,lmsize)

      do icheb = 0,ncheb
         taucsrcr =  tau(icheb)*csrc1(icheb,icheb2)*drpan2
        mn = ipan*ncheb + ipan - icheb
        do lm2 = 1,lmsize
            do lm3 = 1,lmsize
              lm1=lm3+lmsize
              l1 = jlk_index(lm1)

              yill2(icheb,lm3,lm2) = &
              yill2(icheb,lm3,lm2) + &
              taucsrcr*(jlk(l1,mn)*vhli_yill1(lm3,lm2) &
                       -hlk(l1,mn)*vjli_yill1(lm3,lm2))

              zill2(icheb,lm3,lm2) = &
              zill2(icheb,lm3,lm2) + &
              taucsrcr*(jlk(l1,mn)*vhli_zill1(lm3,lm2) &
                       -hlk(l1,mn)*vjli_zill1(lm3,lm2))

            end do
        end do
      end do
    end do

  else
    stop '[rllsll] error in inversion'
  end if

! Reorient indices for later use
  if ( use_sratrick==0 ) then
    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize2
          yif(lm1,lm2,icheb) = yill(icheb,lm1,lm2)
          zif(lm1,lm2,icheb) = zill(icheb,lm1,lm2)
        end do
      end do
    end do

  elseif ( use_sratrick==1 ) then

    do icheb = 0,ncheb
      do lm2 = 1,lmsize
        do lm1 = 1,lmsize
          yif(lm1,lm2,icheb) = yill1(icheb,lm1,lm2)
          zif(lm1,lm2,icheb) = zill1(icheb,lm1,lm2)
          yif(lm1+lmsize,lm2,icheb) = yill2(icheb,lm1,lm2)
          zif(lm1+lmsize,lm2,icheb) = zill2(icheb,lm1,lm2)
        end do
      end do
    end do

  end if

! Calculation of eq. 5.19-5.22

  do icheb = 0,ncheb
    zslc1sum(icheb) = slc1sum(icheb)*drpan2
  end do
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhli(1,1,0), &
        lmsize,yif(1,1,0),lmsize2,czero,mihvy,lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjli(1,1,0), &
        lmsize,yif(1,1,0),lmsize2,czero,mijvy,lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vhli(1,1,0), &
        lmsize,zif(1,1,0),lmsize2,czero,mihvz,lmsize)
  call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(0),vjli(1,1,0), &
        lmsize,zif(1,1,0),lmsize2,czero,mijvz,lmsize)
  do icheb = 1,ncheb
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhli(1,1,icheb), &
          lmsize,yif(1,1,icheb),lmsize2,cone,mihvy,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjli(1,1,icheb), &
          lmsize,yif(1,1,icheb),lmsize2,cone,mijvy,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vhli(1,1,icheb), &
          lmsize,zif(1,1,icheb),lmsize2,cone,mihvz,lmsize)
    call zgemm('n','n',lmsize,lmsize,lmsize2,zslc1sum(icheb),vjli(1,1,icheb), &
          lmsize,zif(1,1,icheb),lmsize2,cone,mijvz,lmsize)
  end do

end subroutine sll_local_solutions


subroutine svpart(srv1,jlmkmn,hlmkmn,vhli,vjli,ncheb,lmsize)
   ! this subroutine facilitates compile optimization by working with
   ! only two-dimensional arrays
      use mod_dataTypes, only: dp
implicit none
   integer :: ncheb,lmsize,icheb,lm1
   complex (kind=dp) :: srv1(0:ncheb,lmsize), &
                     jlmkmn(0:ncheb,lmsize),hlmkmn(0:ncheb,lmsize), &
                     vhli(lmsize),vjli(lmsize)
   
   do lm1 = 1,lmsize
      do icheb = 0,ncheb
         srv1(icheb,lm1) = &
          -jlmkmn(icheb,lm1)*vhli(lm1) &
          +hlmkmn(icheb,lm1)*vjli(lm1)
      end do
   end do
end subroutine svpart

end module mod_sll_local_solutions
