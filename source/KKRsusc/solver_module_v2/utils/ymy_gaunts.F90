  subroutine ymy_gaunts()
! Opens file with the LL mesh, generates spherical harmonics and Gaunts
! Makes use of global arrays ull, wll, cthetall, phill, ylm and gaunt
  use global

  implicit none

  real(kind=r8b)    :: fourpi = 1.d0!16.d0*atan(1.d0)
  integer(kind=i4b) :: ill, ngaunt(0:nlmax4)
  integer(kind=i4b) :: ngauntl(0:nlmax2,0:nlmax2,0:nlmax4)
  integer(kind=i4b) :: ilm, jlm, klm, im, il, jm, jl, km, kl, i2(2)
  real(kind=r8b)    :: norm
  real(kind=r8b), allocatable :: work(:)


  if (allocated(rgaunt)) return


! read Lebedev-Laikov mesh
!  open(file='lebedev_ascii.gga',unit=iofile)
!  read(iofile,*) nll, nlmaxll  ! I added the lmax in the file: nll=434 -> lmax=16
!  write(*,'(" Lebedev-Laikov mesh: nll, lmax=",2i4)') nll, nlmaxll
!  lmmaxll = (nlmaxll+1)**2
!  if (nlmaxll < 4*nlmax) stop 'LL mesh: nlmaxll < 4*nlmax!'
!  allocate(ull(3,nll),wll(nll),cthetall(nll),phill(nll))
!  allocate(ylm(nll,lmmax4),work(lmmax4))
!  allocate(rgaunt(lmmax2,lmmax2,lmmax4))
!  do ill=1,nll
!    read(iofile,*) ull(:,ill), wll(ill), cthetall(ill), phill(ill)
!    if (lhdio) write(iodb,'("ill=",i6,6es10.2)') ill, ull(:,ill), wll(ill), cthetall(ill), phill(ill)
!  end do
!  wll = fourpi*wll
!  write(*,'(" LL weights integrate to",f8.3)') sum(wll)
!  close(iofile)
! use built-in routines
  nlmaxll = 16
  lmmaxll = (nlmaxll+1)**2
  if (nlmaxll < 4*nlmax) stop 'LL mesh: nlmaxll < 4*nlmax!'
! this subroutine sets nll, allocates and fills in ull(3,nll) and wll(nll)
  call llsetup(nlmaxll)
  write(*,'(" Lebedev-Laikov mesh: nll, lmax=",2i4)') nll, nlmaxll
  allocate(ylm(nll,lmmax4),work(lmmax4))
  allocate(rgaunt(lmmax2,lmmax2,lmmax4))
  wll = fourpi*wll
  write(*,'(" LL weights integrate to",f8.3)') sum(wll)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! compute real spherical harmonics
!  open(file='ylm.dat',unit=iofile)
  do ill=1,nll
    call rymy(ull(:,ill),nlmax4,lmmax4,work)
!    call ymy2(cthetall(ill),phill(ill),nlmax2,lmmax2,work)
    ylm(ill,:) = work
!    write(iofile,'(1000es16.8)') ull(:,ill)*abs(work(16))
  end do
!  close(iofile)
! check various things
  do jlm=1,lmmax4
    norm = sum(wll*ylm(:,jlm))
    if (lhdio) write(iodb,'("<yilm>",2i4,es24.12)') i2lm(:,jlm), norm
    do ilm=1,lmmax4
      norm = sum(wll*ylm(:,jlm)*ylm(:,ilm))
      if (abs(norm) > ylmtol) then
        if (lhdio) write(iodb,'("<yilm*yjlm>",4i4,es24.12)') i2lm(:,ilm), i2lm(:,jlm), norm
      end if
    end do
  end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! now Gaunt numbers
!  open(file='gaunt_help.dat',unit=iofile,status='replace')
!  write(iofile,'("# klm jlm ilm  Gaunt number")')
  rgaunt = 0.d0; ngaunt = 0; ngauntl = 0
  do klm=1,lmmax4
    i2 = i2lm(:,klm)
    km = i2(1); kl = i2(2)
    do jlm=1,lmmax2 
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      do ilm=1,lmmax2 
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
!     -------------------------------------------------------------
!                          selection rules
!                 Homeier et al, Theochem 368, p31 (1996)
!     -------------------------------------------------------------
!     m selection rules
        if (abs(km) == abs(im+jm) .or. abs(km) == abs(im-jm)) then
!     parity selection rule
        if (mod(il+jl+kl,2) == 0) then
!     triangle inequality for l
        if (kl <= il + jl .and. kl >= max(abs(il-jl),min(abs(im+jm),abs(im-jm)))) then
          norm = sum(wll*ylm(:,ilm)*ylm(:,jlm)*ylm(:,klm))  ! this is where the work is
!         numerical selection rule
          if (abs(norm) > ylmtol) then
!            write(iodb,'("Gaunt: ",6i6,es8.1)') im, il, jm, jl, km, kl, norm
            rgaunt(ilm,jlm,klm) = norm
!            if (klm < 10 .and. ilm < 5 .and. jlm < 5) write(iofile,'(3i4,f16.8)') klm, jlm, ilm, real(norm)
            ngaunt(kl) = ngaunt(kl) + 1
            if (ngauntl(il,jl,kl) == 0) ngauntl(il,jl,kl) = 1
          end if
        end if
        end if
        end if
!     -------------------------------------------------------------
      end do
    end do
  end do
!  close(iofile)
  write(*,'(" Number of non-zero Gaunt numbers",i8,/)') sum(ngaunt)
  if (lhdio) write(iodb,'("Number of Gaunt number per l")')
  do kl=0,nlmax4
    if (lhdio) write(iodb,'(10i6)') kl, ngaunt(kl)
    do jl=0,nlmax2
      if (lhdio) write(iodb,'(10i2)') ngauntl(:,jl,kl)
    end do
  end do
! Test: m-sums
  do klm=1,lmmax4
    do il=0,nlmax2
      norm = 0.d0
      do im=-il,il
        ilm = lm2i(im,il)
        norm = norm + rgaunt(ilm,ilm,klm)
!        if (abs(rgaunt(ilm,ilm,klm)) > ylmtol) write(iodb,'("klm,ilm",2i4,f10.6)') klm, ilm, rgaunt(ilm,ilm,klm)
      end do
      if (abs(norm) > ylmtol .and. lhdio) write(iodb,'("gaunt m-sum: klm,il=",3i4,f10.6)') i2lm(:,klm), il, norm
    end do
  end do
! Test: Y_L1 Y_L1 = \sum_L3 C(L1,L2,L3) Y_L3
  do jlm=1,lmmax2
    do ilm=1,lmmax2
      do ill=1,nll
        norm = 0.d0
        do klm=1,lmmax4
          if (abs(rgaunt(ilm,jlm,klm)) > ylmtol) norm = norm + rgaunt(ilm,jlm,klm)*ylm(ill,klm)
        end do
        norm = ylm(ill,ilm)*ylm(ill,jlm) - norm
        if (abs(norm) > ylmtol .and. lhdio) write(iodb,'("ylm product to sum: ",4i4,es16.8)') ill, ilm, jlm, klm, norm
      end do
    end do
  end do
! All done!
  end subroutine ymy_gaunts
