module mod_jijhelp

  implicit none
  private

  public :: set_Jijcalc_flags, calc_dtmatJij
  
contains

  subroutine set_Jijcalc_flags(t_dtmatJij,NATYPD,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)

    use mod_types, only: t_inc, type_dtmatJijDij

    implicit none
    type(type_dtmatJijDij), intent(inout) :: t_dtmatJij(t_inc%NATYP)
    integer, intent(in) :: NATYPD,NATOMIMPD,NATOMIMP, ATOMIMP(NATOMIMPD),IQAT(NATYPD)

    integer :: I1, IQ

    DO IQ=1,t_inc%NATYP
     DO I1=1,NATOMIMP
       IF(IQAT(IQ)==ATOMIMP(I1))then
         t_dtmatJij(IQ)%calculate = .true.
       END IF
     END DO!I1
    END DO!IQ

    !Test
!   write(*,*) '=========== TEST FOR TMAT-CALC============='
!   do IQ=1,t_inc%NATYP
!     write(*,'(A,I8,L3)') "atom",IQ, t_dtmatJij(IQ)%calculate
!   end do

  end subroutine set_Jijcalc_flags



  subroutine calc_dtmatJij(LMAXD,LMMAXD,LMMAXSO,LMPOTD,NTOTD,NRMAXD,NSRA,IRMDNEW,NSPIN,VINS,RLLLEFT,RLL,RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,NCLEB,RNEW,dtmat)
! subroutine calc_dtmatJij(NTOTD,NRMAXD,NSRA,IRMDNEW,NSPIN,VINS,RLLLEFT,RLL,RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,NCLEB,RNEW,dtmat)

    implicit none
!   include 'inc.p'
!   INTEGER LMMAXD
!   PARAMETER (LMMAXD= (LMAXD+1)**2)
!   INTEGER LMMAXSO
!   PARAMETER (LMMAXSO=2*LMMAXD)
!   INTEGER LMPOTD
!   PARAMETER (LMPOTD= (LPOTD+1)**2)
    DOUBLE COMPLEX CZERO,CONE
    PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))

    integer, intent(in) :: LMAXD,LMMAXD,LMMAXSO,LMPOTD,NSRA,IRMDNEW,NRMAXD,NSPIN,IEND,NCLEB,NTOTD    !integer arguments that only define array sizes
!   integer, intent(in) :: NSRA,IRMDNEW,NRMAXD,NSPIN,IEND,NCLEB,NTOTD    !integer arguments that only define array sizes
    INTEGER, intent(in) :: NPAN_TOT,NCHEB,IPAN_INTERVALL(0:NTOTD), ICLEB(NCLEB,4) !integer arguments
    DOUBLE PRECISION, intent(in) :: RPAN_INTERVALL(0:NTOTD), VINS(IRMDNEW,LMPOTD,NSPIN), CLEB(*), RNEW(NRMAXD)
    DOUBLE COMPLEX,   intent(in) :: RLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW),&
                                  & RLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW)
    DOUBLE COMPLEX, intent(out) :: dtmat(LMMAXSO,LMMAXSO,3)

    !..locals..
    integer :: ispin1, ispin2, jspin, ir, ishift1, ishift2, lm1, lm2
    double precision :: BINS(IRMDNEW,LMPOTD,1)
    DOUBLE COMPLEX   :: BNSPLL0(LMMAXD,LMMAXD,IRMDNEW), sigma(2,2,3), VNSPLL0(LMMAXSO,LMMAXSO),&
                      & PNSIL(LMMAXSO,LMMAXSO), PNSIR(LMMAXSO,LMMAXSO), CMATTMP(LMMAXSO,LMMAXSO),&
                      & WR(LMMAXSO,LMMAXSO,IRMDNEW), CTMP(IRMDNEW)


    if(NSPIN/=2) stop 'calc_dtmatJij: NSPIN must be 2'
    
!   BINS(:,:,1) = (VINS(:,:,2)-VINS(:,:,1))/2 !Bxc-field  !CLEANUP: how to choose the sign here????
    BINS(:,:,1) = (VINS(:,:,1)-VINS(:,:,2))/2 !Bxc-field

!   write(65161,'(2ES25.16)') VINS

    !convert B_L into B_LL' by using the Gaunt coefficients
    BNSPLL0 = CZERO
    CALL VLLMAT( 1,NRMAXD,IRMDNEW,LMMAXD,LMMAXD,BNSPLL0,BINS, &
               & CLEB,ICLEB,IEND,1,0d0,RNEW,0   )

!   write(65162,'(2ES25.16)')  BNSPLL0
    !get the pauli spin matrices
    call calc_sigma(sigma)
!   write(65163,'(2ES25.16)') sigma


    !loop over sigma_{x,y,z}
    do jspin=1,3!x,y,z
     do ir=1,IRMDNEW

      !construct sigma*B_LL'(r) = VNSPLL0
      VNSPLL0 = CZERO
      do ispin1=1,NSPIN
        ishift1=(ispin1-1)*LMMAXD
        do ispin2=1,NSPIN
          ishift2=(ispin2-1)*LMMAXD
          VNSPLL0(ishift2+1:ishift2+LMMAXD,ishift1+1:ishift1+LMMAXD) = sigma(ispin2,ispin1,jspin)*BNSPLL0(:,:,ir)
        end do!ispin2
      end do!ispin1

     PNSIR(:,:)=RLL(1:LMMAXSO,:,ir)
!    if(jspin==1) write(65164,'(2ES25.16)') PNSIR
     PNSIL(:,:)=RLLLEFT(1:LMMAXSO,:,ir)
!    if(jspin==1) write(65165,'(2ES25.16)') PNSIL

     !calculate [Rleft * VNSPLL0 *Rright](r)
     CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VNSPLL0,&
     &             LMMAXSO,PNSIR,LMMAXSO,CZERO,CMATTMP,LMMAXSO)

     CALL ZGEMM('T','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSIL,&
     &             LMMAXSO,CMATTMP,LMMAXSO,CZERO,WR(:,:,IR),LMMAXSO)

     end do!ir

     !perform radial integration for each matrix element dtmat(LM2,LM1) = \int dr {Rleft * VNSPLL0 *Rright}(r)
     do LM1=1,LMMAXSO
      do LM2=1,LMMAXSO
        CTMP = WR(LM2,LM1,:)
        CALL INTCHEB_CELL(CTMP,dtmat(LM2,LM1,jspin),&
     &             RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
      end do!LM2
     end do!LM1

    end do!jspin

!   stop 'test-stop'

  end subroutine calc_dtmatJij



subroutine calc_sigma(sigma)
implicit none
double complex :: sigma(2,2,3)
integer        :: verbose
character(len=*),parameter :: conventionmode='kkr'
verbose=0

if (conventionmode=='normal') then
  sigma(1,1,1)=( 0.0D0, 0.0D0)
  sigma(1,2,1)=( 1.0D0, 0.0D0)
  sigma(2,1,1)=( 1.0D0, 0.0D0)
  sigma(2,2,1)=( 0.0D0, 0.0D0)
  
  sigma(1,1,2)=( 0.0D0, 0.0D0)
  sigma(1,2,2)=( 0.0D0,-1.0D0)
  sigma(2,1,2)=( 0.0D0, 1.0D0)
  sigma(2,2,2)=( 0.0D0, 0.0D0)
  
  sigma(1,1,3)=( 1.0D0, 0.0D0)
  sigma(1,2,3)=( 0.0D0, 0.0D0)
  sigma(2,1,3)=( 0.0D0, 0.0D0)
  sigma(2,2,3)=(-1.0D0, 0.0D0)
elseif (conventionmode=='kkr') then
  sigma(1,1,1)=( 0.0D0, 0.0D0)
  sigma(1,2,1)=( 1.0D0, 0.0D0)
  sigma(2,1,1)=( 1.0D0, 0.0D0)
  sigma(2,2,1)=( 0.0D0, 0.0D0)
  
  sigma(1,1,2)=( 0.0D0, 0.0D0)
  sigma(1,2,2)=( 0.0D0, 1.0D0)
  sigma(2,1,2)=( 0.0D0,-1.0D0)
  sigma(2,2,2)=( 0.0D0, 0.0D0)
  
  sigma(1,1,3)=(-1.0D0, 0.0D0)
  sigma(1,2,3)=( 0.0D0, 0.0D0)
  sigma(2,1,3)=( 0.0D0, 0.0D0)
  sigma(2,2,3)=( 1.0D0, 0.0D0)
else
  stop'[calc_sigma] wrong mode'
end if

if (verbose==1) then
  write(*,*) '#################################'
  write(*,*) 'calculation OF Pauli matricies'
  write(*,*) '#################################'
  write(*,*) 'sigma_x'
  write(*,'(4f6.2)') sigma(1,1,1), sigma(1,2,1)
  write(*,'(4f6.2)') sigma(2,1,1), sigma(2,2,1)
  
  write(*,*) 'sigma_y'
  write(*,'(4f6.2)') sigma(1,1,2), sigma(1,2,2)
  write(*,'(4f6.2)') sigma(2,1,2), sigma(2,2,2)
  
  write(*,*) 'sigma_z'
  write(*,'(4f6.2)') sigma(1,1,3), sigma(1,2,3)
  write(*,'(4f6.2)') sigma(2,1,3), sigma(2,2,3)
  write(*,*) '#################################'
end if

end subroutine calc_sigma

end module mod_jijhelp
