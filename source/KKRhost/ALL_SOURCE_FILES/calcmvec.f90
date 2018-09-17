module mod_calcmvec

contains

  subroutine calcmvec(nfilcbwf, splitss, iepath, nepath, irel, iprint, nt, nl, mezz, mezj, taut, tsst, iqat, nkmq, nkm, iecurr, netab, igrid, we, mvevdl0, mvevil, bmvevdl0, &
    bmvevil, r2drdi, jrws, imt, amemvec, ikmllim1, ikmllim2, imkmtab, ntmax, nlmax, nmuemax, nqmax, nkmmax, nmmax, nmvecmax, nrmax)
    ! ********************************************************************
    ! *                                                                  *
    ! *                                                                  *
    ! ********************************************************************
    use :: mod_types, only: t_inc
    use :: mod_datatypes, only: dp
    use :: mod_cintabr
    use :: mod_cinit
    implicit complex (kind=dp)(a-h, o-z)

    ! PARAMETER definitions

    complex (kind=dp) :: c0, c1
    parameter (c0=(0.0d0,0.0d0), c1=(1.0d0,0.0d0))
    real (kind=dp) :: pi
    parameter (pi=3.141592653589793238462643d0)
    complex (kind=dp) :: cpre
    parameter (cpre=-c1/pi)

    ! Dummy arguments

    complex (kind=dp) :: we
    integer :: iecurr, iepath, igrid, iprint, irel, nepath, netab, nfilcbwf, nkm, nkmmax, nl, nlmax, nmmax, nmuemax, nmvecmax, nqmax, nrmax, nt, ntmax
    logical :: splitss
    real (kind=dp) :: amemvec(nkmmax, nkmmax, 3, nmvecmax), r2drdi(nrmax, nmmax)
    complex (kind=dp) :: bmvevdl0(nlmax, ntmax, 3, nmvecmax), bmvevil(nlmax, ntmax, 3, nmvecmax), mezj(nkmmax, nkmmax, ntmax, nmvecmax), mezz(nkmmax, nkmmax, ntmax, nmvecmax), &
      mvevdl0(nlmax, ntmax, 3, nmvecmax), mvevil(nlmax, ntmax, 3, nmvecmax), taut(nkmmax, nkmmax, ntmax), tsst(nkmmax, nkmmax, ntmax)
    integer :: ikmllim1(nkmmax), ikmllim2(nkmmax), imkmtab(nkmmax), imt(ntmax), iqat(nqmax, ntmax), jrws(nmmax), nkmq(nqmax)

    ! Local variables

    ! F77--------------------------------------------------------------------
    ! ccc      COMPLEX*16 JF(NRMAX,2,NKMMAX),JG(NRMAX,2,NKMMAX),
    ! ccc     &           ZF(NRMAX,2,NKMMAX),ZG(NRMAX,2,NKMMAX),
    ! F77--------------------------------------------------------------------
    ! F90--------------------------------------------------------------------
    complex (kind=dp), allocatable :: jf(:, :, :), jg(:, :, :), zf(:, :, :), zg(:, :, :)
    ! F90--------------------------------------------------------------------
    complex (kind=dp) :: bmvevd(ntmax, 3, nmvecmax), bmvevdl(nlmax, 3, nmvecmax), bmvevdm(nlmax, nmuemax, 3, nmvecmax), bmvevi(ntmax, 3, nmvecmax), cwgt, &
      meirr(nkmmax, nkmmax, 3, nmvecmax), mereg(nkmmax, nkmmax, 3, nmvecmax), mvevd(ntmax, 3, nmvecmax), mvevdl(nlmax, 3, nmvecmax), mvevdm(nlmax, nmuemax, 3, nmvecmax), &
      mvevi(ntmax, 3, nmvecmax), w1(nkmmax, nkmmax), w2(nkmmax, nkmmax), w3(nkmmax, nkmmax), zfjf(2, 2), zfzf(2, 2), zgjg(2, 2), zgzg(2, 2)
    logical :: check
    integer :: i, i0, ikm, ikm1, ikm2, ikmcb(2, nkmmax), ikmt1, ikmt2, il, im, imkm1, imkm2, imv, imvec, ipol, it, iti, j, j1, j2, jtop, k, k1, k2, kapcb, l, li, lmax, m, mue, &
      muetop, n, nmvec, noswf, npol, nsol, nsolcb(nkmmax)
    real (kind=dp) :: mj, sum
    character (len=3) :: str3

    ! index 3:  ipol= 1,2,3  ==  (+),(-),(z)


    check = .true.
    check = .false.
    npol = 3
    nmvec = min(4, nmvecmax)

    if (iecurr==1 .and. iepath==1) then

      do imv = 1, nmvec
        do ipol = 1, npol
          do it = 1, nt
            do il = 1, nl
              mvevdl0(il, it, ipol, imv) = c0
              mvevil(il, it, ipol, imv) = c0
              bmvevdl0(il, it, ipol, imv) = c0
              bmvevil(il, it, ipol, imv) = c0
            end do
          end do
        end do
      end do

    end if

    noswf = nt*nkm
    nmvec = 3

    ! F90--------------------------------------------------------------------
    allocate (jf(nrmax,2,nkmmax), jg(nrmax,2,nkmmax), stat=it)
    if (it/=0) stop '      < CALCMVEC > : allocate JF/JG '
    allocate (zf(nrmax,2,nkmmax), zg(nrmax,2,nkmmax), stat=it)
    if (it/=0) stop '      < CALCMVEC > : allocate ZF/ZG '
    ! F90--------------------------------------------------------------------
    ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    do it = 1, nt

      m = nkmmax
      n = nkmq(iqat(1,it))
      im = imt(it)
      im = 1
      jtop = jrws(im)

      lmax = nl - 1

      call cinit(nkmmax*nkmmax*3*nmvecmax, mereg)
      call cinit(nkmmax*nkmmax*3*nmvecmax, meirr)

      ! =======================================================================
      ! calculate matrix elements
      ! =======================================================================

      do ikm = 1, n
        read (nfilcbwf, rec=ikm+(it-1)*nkm) iti, li, mj, nsol, str3, (kapcb, ikmcb(k,ikm), (zg(i,k,ikm),zf(i,k,ikm),i=1,jtop), k=1, nsol)
        if (it/=iti) stop ' IT(INI) <> IT  in <MENABIRR>'
        if (str3/='REG') stop 'WFT(INI) <> REG in <MENABIRR>'

        read (nfilcbwf, rec=ikm+(it-1)*nkm+noswf) iti, li, mj, nsol, str3, (kapcb, ikmcb(k,ikm), (jg(i,k,ikm),jf(i,k,ikm),i=1,jtop), k=1, nsol)
        nsolcb(ikm) = nsol
        if (it/=iti) stop ' IT(INI) <> IT  in <MENABIRR>'
        if (str3/='IRR') stop 'WFT(INI) <> IRR in <MENABIRR>'
        if (iprint>3 .and. (t_inc%i_write>0)) write (1337, *) iti, li, mj, nsol, str3, kapcb
      end do


      do ikmt2 = 1, n

        do ikmt1 = ikmllim1(ikmt2), ikmllim2(ikmt2)

          do imvec = 1, nmvec
            do k2 = 1, nsolcb(ikmt2)
              j2 = ikmcb(k2, ikmt2)
              do k1 = 1, nsolcb(ikmt1)
                j1 = ikmcb(k1, ikmt1)
                do ipol = 1, npol
                  if (abs(amemvec(j1,j2,ipol,imvec))>1e-8) go to 100
                end do
              end do
            end do
          end do
          ! -------------------------------------- all angular matrix elements =
          ! 0
          go to 110
          ! ---------------------------------- non-0 angular matrix elements
          ! found
          ! ------------------------------------- calculate radial matrix
          ! elements

100       continue
          call cintabr(zg(1,1,ikmt1), zg(1,1,ikmt2), zgzg, zf(1,1,ikmt1), zf(1,1,ikmt2), zfzf, r2drdi(1,im), nsolcb(ikmt1), nsolcb(ikmt2), jtop, nrmax)

          call cintabr(zg(1,1,ikmt1), jg(1,1,ikmt2), zgjg, zf(1,1,ikmt1), jf(1,1,ikmt2), zfjf, r2drdi(1,im), nsolcb(ikmt1), nsolcb(ikmt2), jtop, nrmax)

          ! -------------------------------------- calculate total matrix
          ! elements

          do k2 = 1, nsolcb(ikmt2)
            ikm2 = ikmcb(k2, ikmt2)
            imkm2 = imkmtab(ikm2)

            do k1 = 1, nsolcb(ikmt1)
              ikm1 = ikmcb(k1, ikmt1)
              imkm1 = imkmtab(ikm1)

              do imv = 1, nmvec
                do ipol = 1, npol
                  mereg(ikmt1, ikmt2, ipol, imv) = mereg(ikmt1, ikmt2, ipol, imv) + amemvec(ikm1, ikm2, ipol, imv)*zgzg(k1, k2)
                end do
              end do

            end do

            do imv = 1, nmvec
              do ipol = 1, npol
                mereg(ikmt2, ikmt2, ipol, imv) = mereg(ikmt2, ikmt2, ipol, imv) - amemvec(imkm2, imkm2, ipol, imv)*zfzf(k2, k2)
              end do
            end do

          end do

          if (ikmt1==ikmt2) then
            do k2 = 1, nsolcb(ikmt2)
              ikm2 = ikmcb(k2, ikmt2)
              imkm2 = imkmtab(ikm2)
              do k1 = 1, nsolcb(ikmt1)
                ikm1 = ikmcb(k1, ikmt1)
                imkm1 = imkmtab(ikm1)

                do imv = 1, nmvec
                  do ipol = 1, npol
                    meirr(ikmt1, ikmt2, ipol, imv) = meirr(ikmt1, ikmt2, ipol, imv) + amemvec(ikm1, ikm2, ipol, imv)*zgjg(k1, k2)
                  end do
                end do

              end do

              do imv = 1, nmvec
                do ipol = 1, npol
                  meirr(ikmt2, ikmt2, ipol, imv) = meirr(ikmt2, ikmt2, ipol, imv) - amemvec(imkm2, imkm2, ipol, imv)*zfjf(k2, k2)
                end do
              end do

            end do
          end if

110     end do

      end do

      if (check) then
        do i = 1, nkm
          do j = ikmllim1(i), ikmllim2(i)
            sum = abs(mezz(i,j,it,1)) + abs(mezj(i,j,it,1))
            sum = sum + abs(mereg(i,j,3,1)) + abs(meirr(i,j,3,1))
            if (sum>1d-8 .and. (t_inc%i_write>0)) then
              write (1337, *) ' spin '
              write (1337, '(2i3,2e17.8,2x,2e17.8)') i, j, mezz(i, j, it, 2), mezj(i, j, it, 2)
              write (1337, '(6x,2e17.8,2x,2e17.8)') mereg(i, j, 3, 1), meirr(i, j, 3, 1), (mezz(i,j,it,2)-mereg(i,j,3,1)), (mezj(i,j,it,2)-meirr(i,j,3,1))
              write (1337, *) ' orb '
              write (1337, '(2i3,2e17.8,2x,2e17.8)') i, j, mezz(i, j, it, 3), mezj(i, j, it, 3)
              write (1337, '(6x,2e17.8,2x,2e17.8)') mereg(i, j, 3, 2), meirr(i, j, 3, 2), (mezz(i,j,it,3)-mereg(i,j,3,2)), (mezj(i,j,it,3)-meirr(i,j,3,2))
            end if
          end do
        end do
      end if

      ! =======================================================================
      do imv = 1, nmvec
        do ipol = 1, npol

          mvevd(it, ipol, imv) = 0.0d0
          mvevi(it, ipol, imv) = 0.0d0

          if (.not. splitss) then
            bmvevd(it, ipol, imv) = 0.0d0
            bmvevi(it, ipol, imv) = 0.0d0
          end if

          call zgemm('N', 'N', n, n, n, cpre, mereg(1,1,ipol,imv), m, taut(1,1,it), m, c0, w1, m)
          cwgt = -1d0
          do j = 1, n
            call zaxpy(n, -cpre, meirr(1,j,ipol,imv), 1, w1(1,j), 1)
            call zcopy(n, taut(1,j,it), 1, w2(1,j), 1)
            call zaxpy(n, cwgt, tsst(1,j,it), 1, w2(1,j), 1)
          end do
          call zgemm('N', 'N', n, n, n, cpre, mereg(1,1,ipol,imv), m, w2, m, c0, w3, m)

          ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
          do l = 0, lmax
            il = l + 1

            if (irel<=1) then
              i0 = l*l
              do mue = 1, 2*l + 1
                mvevdm(il, mue, ipol, imv) = w1(i0+mue, i0+mue)
                bmvevdm(il, mue, ipol, imv) = w3(i0+mue, i0+mue)
              end do
            else
              i0 = 2*l*l + l + l
              muetop = 2*l + 2
              do mue = 1, muetop
                mvevdm(il, mue, ipol, imv) = w1(i0+mue, i0+mue)
                bmvevdm(il, mue, ipol, imv) = w3(i0+mue, i0+mue)
              end do
              i0 = 2*(l-1)*l + l + l
              do mue = 2, muetop - 1
                mvevdm(il, mue, ipol, imv) = mvevdm(il, mue, ipol, imv) + w1(i0+mue-1, i0+mue-1)
                bmvevdm(il, mue, ipol, imv) = bmvevdm(il, mue, ipol, imv) + w3(i0+mue-1, i0+mue-1)
              end do
            end if

            mvevdl(il, ipol, imv) = 0.0d0

            if (.not. splitss) bmvevdl(il, ipol, imv) = 0.0d0

            if (irel>1) then
              muetop = 2*l + 2
            else
              muetop = 2*l + 1
            end if

            ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            do mue = 1, muetop
              mvevdl(il, ipol, imv) = mvevdl(il, ipol, imv) + mvevdm(il, mue, ipol, imv)

              if (.not. splitss) bmvevdl(il, ipol, imv) = bmvevdl(il, ipol, imv) + bmvevdm(il, mue, ipol, imv)
            end do
            ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

            mvevil(il, it, ipol, imv) = mvevil(il, it, ipol, imv) + we*mvevdl(il, ipol, imv)
            mvevdl0(il, it, ipol, imv) = mvevdl(il, ipol, imv)

            if (.not. splitss) then
              bmvevil(il, it, ipol, imv) = bmvevil(il, it, ipol, imv) + we*bmvevdl(il, ipol, imv)
              bmvevdl0(il, it, ipol, imv) = bmvevdl(il, ipol, imv)
            end if

            if (igrid/=4 .or. iecurr<=netab) then

              mvevd(it, ipol, imv) = mvevd(it, ipol, imv) + mvevdl(il, ipol, imv)
              mvevi(it, ipol, imv) = mvevi(it, ipol, imv) + mvevil(il, it, ipol, imv)

              if (.not. splitss) then
                bmvevd(it, ipol, imv) = bmvevd(it, ipol, imv) + bmvevdl(il, ipol, imv)
                bmvevi(it, ipol, imv) = bmvevi(it, ipol, imv) + bmvevil(il, it, ipol, imv)
              end if
            end if

          end do
          ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
        end do
      end do

    end do
    ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    ! F90--------------------------------------------------------------------
    deallocate (jf, jg, zf, zg, stat=it)
    if (it/=0) stop '      < CALCMVEC > : deallocate JF/JG/ZF/ZG'
    ! F90--------------------------------------------------------------------
    if (splitss .and. ((iepath==1) .and. (iecurr==netab))) then
      do imv = 1, nmvec
        do ipol = 1, npol
          do it = 1, nt
            do il = 1, nl
              bmvevdl0(il, it, ipol, imv) = mvevdl0(il, it, ipol, imv)
            end do
          end do
        end do
      end do
    end if

    ! =======================================================================
    if ((igrid>=6) .or. (iecurr/=netab) .or. (iepath/=nepath)) return

    ! this part of the original Munich subroutine has been moved to
    !! mvecglobal >
    ! main2 --> tbkkr2 --> mvecglobal -- see makefile2

  end subroutine calcmvec

end module mod_calcmvec
