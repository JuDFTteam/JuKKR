subroutine irwns(cr, dr, efac, qns, vnspll, icst, ipan, ircut, nsra,  &
                  pzlm, qzlm, pzekdr, qzekdr, cder, cmat, dder, dmat,  &
                  irmind, irmd, ipand, lmmaxd)
  implicit none
  !-----------------------------------------------------------------------
  !     determines the irregular non spherical wavefunctions in the n-th.
  !       born approximation ( n given by input parameter icst ) .
  !
  !
  !     using the wave functions pz and qz ( regular and irregular
  !       solution ) of the spherically averaged potential , the ir-
  !       regular wavefunction qns is determined by
  !
  !           qns(ir,lm1,lm2) = cr(ir,lm1,lm2)*pz(ir,l1)
  !
  !                                   + dr(ir,lm1,lm2)*qz(ir,l1)
  !
  !      the matrices cr and dr are determined by integral equations
  !        containing qns and only the non spherical contributions of
  !        the potential , stored in vinspll . these integral equations
  !        are solved iteratively with born approximation up to given n.
  !
  !     the original way of writing the cr and dr matrices in the equa-
  !        tion above caused numerical troubles . therefore here are used
  !        rescaled cr and dr matrices (compare subroutine wftsca):
  !
  !              ~
  !              cr(ir,lm1,lm2) = sqrt(e)**(l1+l2)
  !                             * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
  !
  !              ~
  !              dr(ir,lm1,lm2) = sqrt(e)**(l2-l1)
  !                             * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)
  !
  !     attention :  the sign of the dr matrix is changed to reduce the
  !     ===========  number of floating point operations
  !
  !     modified for the use of shape functions
  !
  !                              (see notes by b.drittler)
  !
  !                                b.drittler   mar.  1989
  !-----------------------------------------------------------------------
  !     modified by R. Zeller      Aug. 1994
  !-----------------------------------------------------------------------
  integer, intent(in) :: icst,ipan,ipand,irmd,irmind,lmmaxd,nsra
  double complex, intent(out)   :: cder(lmmaxd,lmmaxd,irmind:irmd), dder(lmmaxd,lmmaxd,irmind:irmd)
  double complex, intent(inout) :: cmat(lmmaxd,lmmaxd,irmind:irmd), dmat(lmmaxd,lmmaxd,irmind:irmd)
  double complex, intent(out)   :: cr(lmmaxd,lmmaxd), dr(lmmaxd,lmmaxd)
  double complex, intent(inout) :: qns(lmmaxd,lmmaxd,irmind:irmd,2)
  double complex, intent(in)    :: efac(lmmaxd)
  double complex, intent(in)    :: pzlm(lmmaxd,irmind:irmd,2)
  double complex, intent(in)    :: pzekdr(lmmaxd,irmind:irmd,2), qzekdr(lmmaxd,irmind:irmd,2)
  double complex, intent(in)    :: qzlm(lmmaxd,irmind:irmd,2)

  double precision, intent(in) :: vnspll(lmmaxd,lmmaxd,irmind:irmd)
  integer, intent(in) :: ircut(0:ipand)
  
  external :: csinwd, wfint, wfint0
  double complex, parameter :: CONE=(1.d0,0.d0)
  double complex :: efac2
  integer :: i, ir, irc1,j, lm
  
  irc1 = ircut(ipan)

  do i = 0, icst
    !---> set up integrands for i-th born approximation
    if (i == 0) then
      call wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
    else
      call wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
    endif ! first Born iteration

    !---> call integration subroutines
    call csinwd(cder,cmat,lmmaxd**2,irmind,irmd,ipan,ircut)
    call csinwd(dder,dmat,lmmaxd**2,irmind,irmd,ipan,ircut)
    do ir = irmind, irc1
      do lm = 1, lmmaxd
        dmat(lm,lm,ir) = dmat(lm,lm,ir) - CONE
      enddo ! lm
    enddo ! ir

    !---> calculate non sph. wft. in i-th born approximation
    do j = 1, nsra
      do ir = irmind, irc1
        do lm = 1, lmmaxd
          qns(:,lm,ir,j) = cmat(:,lm,ir)*pzlm(:,ir,j) - dmat(:,lm,ir)*qzlm(:,ir,j)
        enddo ! lm
      enddo ! ir
    enddo ! j
    
  enddo ! i

  !---> store c - and d - matrix
  cr(:,:) =  cmat(:,:,irmind)
  dr(:,:) = -dmat(:,:,irmind)

  !---> rescale with efac
  do j = 1, nsra
    do lm = 1, lmmaxd
      efac2 = 1.d0/efac(lm)
      do ir = irmind, irc1
        qns(:,lm,ir,j) = qns(:,lm,ir,j)*efac2
      enddo ! ir
    enddo ! lm
  enddo ! j

  endsubroutine
