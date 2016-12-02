  subroutine etotb1(ecou, epotin, espc, espv, exc, euldau, edcldau, ldau, kpre, lmax, lpot, lcoremax, nspin, i1, naez)
! ************************************************************************
!     calculate the total energy of the cluster .
!     gather all energy-parts which are calculated in different
!     subroutines .
!     since the program uses group theory only shell-indices
!     are used instead of atom-indices .
!
!                               b.drittler   may 1987
!
!     modified for supercells with nshell(i) atoms of type i in the
!     unit cell
!                               p.zahn       oct. 95
!
!
!-----------------------------------------------------------------------
!
!  main2
!    | 
!        +- results
!    |     |
!    |     +-etotb1
!
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: kpre, lmax, lpot, naez, nspin, lcoremax
    double precision, intent(in) :: ecou(0:lpot), epotin, espc(0:3,nspin), espv(0:lmax+1,nspin), exc(0:lpot), euldau, edcldau

    ! .. locals ..
    double precision :: bandet, ecous, edc, et, excs
    integer :: i1, ipot, is, ispin, l
    logical :: ldau
    character(len=*), parameter :: textl(0:6) = [' s =',' p =',' d =',' f =',' g =',' h =',' i =']
    character(len=*), parameter :: textns = ' ns ='
    character(len=*), parameter :: texts(3) = [' spin down   ',' spin  up    ',' paramagnetic']
!   logical, external :: test
    double precision, save :: etot, bandesum
    character(len=*), parameter :: FMTq = "(7X,a,2(a4,f15.8),/,(31X,2(a4,f15.8)))", &
            hLINE = "(5X,68('-'))", FMTe = "(5X,a,2(i3,1X,f15.8),/,(29X,2(i5,1X,f15.8)))"
  
    if (i1 == 1) then
      etot = 0.d0
      bandesum = 0.d0
      if (kpre == 1) write(6, fmt="(32('='),' TOTAL ENERGIES ',31('='),/)")
    endif

    if (kpre == 1) write(6, fmt="(3X,'Total energies atom ',i5,/,3x,23('-'))") i1

    edc = 0.d0
    et = 0.d0
    bandet = 0.d0
    ecous = 0.d0
    excs = 0.d0
    
    is = 0
    if (nspin == 1) is = is + 2
    do ispin = 1,nspin
      is = is + 1
      ipot = (i1-1)*nspin + ispin

      if (kpre == 1) then
          write(6, fmt="(5X,'single particle energies ',a13)") texts(is)
          write(6, fmt=FMTq) '  core   contribution : ', (textl(l), espc(l,ispin), l=0,lcoremax)
          write(6, fmt=FMTq) 'valence  contribution : ', (textl(l), espv(l,ispin), l=0,lmax)
          write(6, fmt="(7x,'                        ',a4,f15.8)") textns, espv(lmax+1,ispin)
      endif ! kpre

      do l = 0, lcoremax
        et = et + espc(l,ispin)
      enddo ! l

      do l = 0, lmax
          bandet = bandet + espv(l,ispin)
          et = et + espv(l,ispin)
      enddo ! l
      bandet = bandet + espv(lmax+1,ispin)
      et = et + espv(lmax+1,ispin)
    enddo ! ispin

    et = et + euldau
    bandet = bandet + euldau
    if ((ldau).and.(euldau /= 0.d0)) write(6, fmt="(5X,'LDA+U correction to the single particle energy   :',F16.8)") euldau
!
! --->  sum up coulomb and ex.-corel. contribution
!
    do l = 0, lpot
      ecous = ecous + ecou(l)
      excs = excs + exc(l)
    enddo ! l

    if (kpre == 1) then
      write(6, fmt=hLINE) ! horizontal line
      write(6, fmt="(5X,'total contribution of the single particle energies :',1X,f15.8)") et
      write(6, fmt="(5X,'                              band energy per atom :',1X,f15.10,/)") bandet
      write(6, fmt=FMTe) 'coulomb  contribution : ', (l, ecou(l), l=0,lpot)
      write(6, fmt=hLINE) ! horizontal line
      write(6, fmt="(5X,'tot. coulomb contribution : ',f15.8,/)") ecous
      write(6, fmt=FMTe) 'ex.-cor. contribution : ', (l,exc(l), l=0,lpot)
      write(6, fmt=hLINE) ! horizontal line
      write(6, fmt="(5X,'tot. ex.-cor. contribution : ',f15.8,/)") excs
      write(6, fmt="(5X,'eff. pot. contribution     : ',f15.8)") epotin
    endif ! kpre
    
    et = et + ecous + excs
    edc = edc + ecous + excs

    et = et + epotin - edcldau
    edc = edc + epotin - edcldau
    if ((ldau).and.(edcldau /= 0.d0)) write(6, fmt="(5X,'LDA+U double counting contribution               :',F16.8)") -edcldau

    if (kpre == 1) then
      write(6, fmt="(5X,'total double counting contribution                 :',1X,f15.8)") edc
      if ((ldau).and.(euldau - edcldau) /= 0.d0) write(6, fmt="(3X,'   including LDA+U correction :',F15.8)") euldau - edcldau
    endif ! kpre

    if (naez > 1) write(6, fmt="(3x,'Total contribution of atom',i5,' =',f15.8)") i1, et

    etot = etot + et
    bandesum = bandesum + bandet

    if (i1 == naez) write(6, fmt="(5X,'        sum of band energies :',1X,F20.10,/,3X,70('-'))") bandesum
!     if (i1 == naez) write(6, fmt="(/,3X,70('+'),/,15x,'TOTAL ENERGY in ryd. : ',f21.8,/15x,'                 eV  : ',F21.8,/,3x,70('+'))") etot, etot*13.605698066d0

  endsubroutine
