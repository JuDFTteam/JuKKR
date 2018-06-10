module mod_jijhelp

  implicit none
  private

  public :: set_jijcalc_flags, calc_dtmatjij

contains

  subroutine set_jijcalc_flags(t_dtmatjij, natypd, natomimpd, natomimp, &
    atomimp, iqat)

    use :: mod_types, only: t_inc, type_dtmatjijdij
    use :: mod_save_wavefun, only: t_wavefunctions

    implicit none
    type (type_dtmatjijdij), intent (inout) :: t_dtmatjij(t_inc%natyp)
    integer, intent (in) :: natypd, natomimpd, natomimp, atomimp(natomimpd), &
      iqat(natypd)

    integer :: i1, iq

    do iq = 1, t_inc%natyp
      do i1 = 1, natomimp
        if (iqat(iq)==atomimp(i1)) then
          t_dtmatjij(iq)%calculate = .true.
        end if
      end do !I1
    end do !IQ

! change saveing option for wavefunctions
    t_wavefunctions%save_rll = .true.
    t_wavefunctions%save_sll = .true.
    t_wavefunctions%save_rllleft = .true.
    t_wavefunctions%save_sllleft = .true.

!Test
!   write(*,*) '=========== TEST FOR TMAT-CALC============='
!   do IQ=1,t_inc%NATYP
!     write(*,'(A,I8,L3)') "atom",IQ, t_dtmatJij(IQ)%calculate
!   end do

  end subroutine



  subroutine calc_dtmatjij(lmmaxd, lmmaxso, lmpotd, ntotd, nrmaxd, nsra, &
    irmdnew, nspin, vins, rllleft, rll, rpan_intervall, ipan_intervall, &
    npan_tot, ncheb, cleb, icleb, iend, ncleb, rnew, dtmat)
! subroutine calc_dtmatJij(NTOTD,NRMAXD,NSRA,IRMDNEW,NSPIN,VINS,RLLLEFT,RLL,RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,NCLEB,RNEW,dtmat)

    implicit none
    double complex :: czero, cone
    parameter (czero=(0d0,0d0), cone=(1d0,0d0))

    integer, intent (in) :: lmmaxd, lmmaxso, lmpotd, nsra, irmdnew, nrmaxd, &
      nspin, iend, ncleb, ntotd !integer arguments that only define array sizes
!   integer, intent(in) :: NSRA,IRMDNEW,NRMAXD,NSPIN,IEND,NCLEB,NTOTD    !integer arguments that only define array sizes
    integer, intent (in) :: npan_tot, ncheb, ipan_intervall(0:ntotd), &
      icleb(ncleb, 4) !integer arguments
    double precision, intent (in) :: rpan_intervall(0:ntotd), &
      vins(irmdnew, lmpotd, nspin), cleb(*), rnew(nrmaxd)
    double complex, intent (in) :: rll(nsra*lmmaxso, lmmaxso, irmdnew), &
      rllleft(nsra*lmmaxso, lmmaxso, irmdnew)
    double complex, intent (out) :: dtmat(lmmaxso, lmmaxso, 3)

!..locals..
    integer :: ispin1, ispin2, jspin, ir, ishift1, ishift2, lm1, lm2
    double precision :: bins(irmdnew, lmpotd, 1)
    double complex :: bnspll0(lmmaxd, lmmaxd, irmdnew), sigma(2, 2, 3), &
      vnspll0(lmmaxso, lmmaxso), pnsil(lmmaxso, lmmaxso), &
      pnsir(lmmaxso, lmmaxso), cmattmp(lmmaxso, lmmaxso), &
      wr(lmmaxso, lmmaxso, irmdnew), ctmp(irmdnew)


    if (nspin/=2) stop 'calc_dtmatJij: NSPIN must be 2'

!   BINS(:,:,1) = (VINS(:,:,2)-VINS(:,:,1))/2 !Bxc-field  !CLEANUP: how to choose the sign here????
    bins(:, :, 1) = (vins(:,:,1)-vins(:,:,2))/2 !Bxc-field


!convert B_L into B_LL' by using the Gaunt coefficients
    bnspll0 = czero
    call vllmat(1, nrmaxd, irmdnew, lmmaxd, lmmaxd, bnspll0, bins, cleb, &
      icleb, iend, 1, 0d0, rnew, 0)

!get the pauli spin matrices
    call calc_sigma(sigma) !use this to perform automatically infinitesimal rotations


!loop over sigma_{x,y,z}
    do jspin = 1, 2 !x,y,z
      do ir = 1, irmdnew

!construct sigma*B_LL'(r) = VNSPLL0
        vnspll0 = czero
        do ispin1 = 1, nspin
          ishift1 = (ispin1-1)*lmmaxd
          do ispin2 = 1, nspin
            ishift2 = (ispin2-1)*lmmaxd
            vnspll0(ishift2+1:ishift2+lmmaxd, ishift1+1:ishift1+lmmaxd) &
              = sigma(ispin2, ispin1, jspin)*bnspll0(:, :, ir)
          end do !ispin2
        end do !ispin1

        pnsir(:, :) = rll(1:lmmaxso, :, ir)
        pnsil(:, :) = rllleft(1:lmmaxso, :, ir)

!calculate [Rleft * VNSPLL0 *Rright](r)
        call zgemm('N', 'N', lmmaxso, lmmaxso, lmmaxso, cone, vnspll0, &
          lmmaxso, pnsir, lmmaxso, czero, cmattmp, lmmaxso)

        call zgemm('T', 'N', lmmaxso, lmmaxso, lmmaxso, cone, pnsil, lmmaxso, &
          cmattmp, lmmaxso, czero, wr(:,:,ir), lmmaxso)

!add small component
        if (nsra==2) then
          pnsir(:, :) = rll(lmmaxso+1:2*lmmaxso, :, ir)
          pnsil(:, :) = -rllleft(lmmaxso+1:2*lmmaxso, :, ir)

!calculate [Rleft_small * VNSPLL0 *Rright_small](r) and add to large
!component
          call zgemm('N', 'N', lmmaxso, lmmaxso, lmmaxso, cone, vnspll0, &
            lmmaxso, pnsir, lmmaxso, czero, cmattmp, lmmaxso)

          call zgemm('T', 'N', lmmaxso, lmmaxso, lmmaxso, cone, pnsil, &
            lmmaxso, cmattmp, lmmaxso, cone, wr(:,:,ir), lmmaxso)
        end if !NSRA

      end do !ir

!perform radial integration for each matrix element dtmat(LM2,LM1) = \int dr {Rleft * VNSPLL0 *Rright}(r)
      do lm1 = 1, lmmaxso
        do lm2 = 1, lmmaxso
          ctmp = wr(lm2, lm1, :)
          call intcheb_cell(ctmp, dtmat(lm2,lm1,jspin), rpan_intervall, &
            ipan_intervall, npan_tot, ncheb, irmdnew)
        end do !LM2
      end do !LM1

    end do !jspin

  end subroutine


  subroutine calclambda(lambda, theta, phi)
    implicit none
    double complex :: lambda(2, 2, 3)
    double precision :: theta, phi
    double complex :: sigmatemp(2, 2, 3), sigma(2, 2, 3)
    integer :: ispin

    call calc_sigma(sigma)
    sigmatemp = sigma
    do ispin = 1, 3
      call rotatematrix(sigmatemp(:,:,ispin), theta, phi, 1, 1)
      lambda(:, :, ispin) = sigmatemp(:, :, ispin) !test
    end do !ispin
  end subroutine




  subroutine calc_sigma(sigma)
    implicit none
    double complex :: sigma(2, 2, 3)
    integer :: verbose
    character (len=*), parameter :: conventionmode = 'kkr'

    verbose = 0

    if (conventionmode=='normal') then
      sigma(1, 1, 1) = (0.0d0, 0.0d0)
      sigma(1, 2, 1) = (1.0d0, 0.0d0)
      sigma(2, 1, 1) = (1.0d0, 0.0d0)
      sigma(2, 2, 1) = (0.0d0, 0.0d0)

      sigma(1, 1, 2) = (0.0d0, 0.0d0)
      sigma(1, 2, 2) = (0.0d0, -1.0d0)
      sigma(2, 1, 2) = (0.0d0, 1.0d0)
      sigma(2, 2, 2) = (0.0d0, 0.0d0)

      sigma(1, 1, 3) = (1.0d0, 0.0d0)
      sigma(1, 2, 3) = (0.0d0, 0.0d0)
      sigma(2, 1, 3) = (0.0d0, 0.0d0)
      sigma(2, 2, 3) = (-1.0d0, 0.0d0)
    else if (conventionmode=='kkr') then
      sigma(1, 1, 1) = (0.0d0, 0.0d0)
      sigma(1, 2, 1) = (1.0d0, 0.0d0)
      sigma(2, 1, 1) = (1.0d0, 0.0d0)
      sigma(2, 2, 1) = (0.0d0, 0.0d0)

      sigma(1, 1, 2) = (0.0d0, 0.0d0)
      sigma(1, 2, 2) = (0.0d0, 1.0d0)
      sigma(2, 1, 2) = (0.0d0, -1.0d0)
      sigma(2, 2, 2) = (0.0d0, 0.0d0)

      sigma(1, 1, 3) = (-1.0d0, 0.0d0)
      sigma(1, 2, 3) = (0.0d0, 0.0d0)
      sigma(2, 1, 3) = (0.0d0, 0.0d0)
      sigma(2, 2, 3) = (1.0d0, 0.0d0)
    else
      stop '[calc_sigma] wrong mode'
    end if

    if (verbose==1) then
      write (*, *) '#################################'
      write (*, *) 'calculation OF Pauli matricies'
      write (*, *) '#################################'
      write (*, *) 'sigma_x'
      write (*, '(4f6.2)') sigma(1, 1, 1), sigma(1, 2, 1)
      write (*, '(4f6.2)') sigma(2, 1, 1), sigma(2, 2, 1)

      write (*, *) 'sigma_y'
      write (*, '(4f6.2)') sigma(1, 1, 2), sigma(1, 2, 2)
      write (*, '(4f6.2)') sigma(2, 1, 2), sigma(2, 2, 2)

      write (*, *) 'sigma_z'
      write (*, '(4f6.2)') sigma(1, 1, 3), sigma(1, 2, 3)
      write (*, '(4f6.2)') sigma(2, 1, 3), sigma(2, 2, 3)
      write (*, *) '#################################'
    end if

  end subroutine

end module
