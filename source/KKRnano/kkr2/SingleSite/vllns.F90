  subroutine vllns(vnspll, vins, cleb, icleb, iend, lmax, irmd, irnsd, ncleb)
!-----------------------------------------------------------------------
!     calculates v_ll' from v_l.
!     to determine the non - spherical wavefunctions the potential
!         has to be lm1 and lm2 dependent . the potential is stored
!         only as lm dependent , therefore a transformation in the
!         following way has to be done :
!
!        vnsll(r,lm1,lm2)   =   {  c(lm1,lm2,lm3) *vins(r,lm3)  }
!                                  (summed over lm3 at the right site )
!        where c(lm1,lm2,lm3) are the gaunt coeffients .
!
!             (see notes by b.drittler)
!
!     attention : the gaunt coeffients are stored in an index array
!                  only for lm1.gt.lm2
!                 (see subroutine gaunt)
!
!                               b.drittler   july 1988
!-----------------------------------------------------------------------
!                          modified by r. zeller sep. 2000
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: lmax, irmd, irnsd, ncleb, iend
    double precision, intent(in) :: cleb(ncleb,2)
    double precision, intent(in) :: vins(irmd-irnsd:irmd,(2*lmax+1)**2)
    double precision, intent(inout) :: vnspll((lmax+1)**2,(lmax+1)**2,irmd-irnsd:irmd)
    integer, intent(in) :: icleb(ncleb,3)

    integer :: ir, j, lm1, lm2, lm3, lmmaxd, irmind
    
    irmind = irmd-irnsd
    lmmaxd = (lmax+1)**2

    do ir = irmind, irmd
      do lm1 = 1, lmmaxd
        do lm2 = 1, lm1
          vnspll(lm1,lm2,ir) = 0.d0 ! clear upper triangular matrix including diagonal
        enddo ! lm2
      enddo ! lm1
    enddo ! ir

    do j = 1, iend
      lm1 = icleb(j,1)
      lm2 = icleb(j,2)
      lm3 = icleb(j,3)
      do ir = irmind, irmd
        vnspll(lm1,lm2,ir) = vnspll(lm1,lm2,ir) + cleb(j,1)*vins(ir,lm3)
      enddo ! ir
    enddo ! j

!
!---> use symmetry of the gaunt coef.
!
    do ir = irmind, irmd
      do lm1 = 1, lmmaxd
        do lm2 = 1, lm1-1
          vnspll(lm2,lm1,ir) = vnspll(lm1,lm2,ir) ! copy lower trinagular matrix
        enddo ! lm2
        vnspll(lm1,lm1,ir) = vnspll(lm1,lm1,ir) + vins(ir,1) ! add diagonal terms
      enddo ! lm1
    enddo ! ir

! todo: check if an out loop over ir for all ops is better (will depend on ratio between iend and #ir)    
  endsubroutine
