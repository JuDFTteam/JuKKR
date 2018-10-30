  subroutine build_kxcalda(kxclm,kxcylm,kxcy00)
! assembles the xc ALDA kernel in the projection basis

  implicit none

  real(kind=r8b),    intent(in)  :: kxclm(nrmax,lmmax2,4,4,nasusc)
  complex(kind=c8b), intent(out) :: kxcylm(4,4,lmmax0,nasusc)
  complex(kind=c8b), intent(out) :: kxcy00(4,4,lmmax,lmmax,nasusc)
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: fourpi = 16.d0*atan(1.d0)
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: i2(2), i3(3), ie, ia, ja, jlm, ilm, j, i, nr, iq, jq
  real(kind=r8b)    :: dr(nrmax), maxnorm, dummy(5)
  real(kind=r8b)    :: kxc(nrmax), gaunti, gauntj, maxelem, doublegaunt
  complex(kind=c8b) :: work(nrmax), norm
  complex(kind=c8b) :: kerspin(nrmax,2,2,2,2)
  integer(kind=i4b) :: i1max, i2max, i3max, i4max

  kxcylm = 0.d0; kxcy00 = 0.d0; kernel = 0.d0
! use the ASA
! use the fact that the basis was constructed per l channel
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ia=1,nasusc    ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)
    maxnorm = 0.d0
!   xc kernel, from cartesian to spin labels
!   only spherical part
    do s4=1,nsmax
    do s3=1,nsmax
      do s2=1,nsmax
      do s1=1,nsmax
        kerspin(:,s1,s2,s3,s4) = 0.d0
        do j=1,4
        do i=1,4
          kerspin(1:nr,s1,s2,s3,s4) = kerspin(1:nr,s1,s2,s3,s4) + pc2s(s1,s2,i)*kxclm(1:nr,1,i,j,ia)*ds2c(j,s3,s4)
        end do
        end do
      end do
      end do
    end do
    end do
!    open(file="kerxc.scf",unit=iofile)
!    read(iofile,*)
!    do i=1,nr
!      read(iofile,*) dummy(1:5), kxc(i)
!    end do
!    kxc(1:nr) = kxc(1:nr)/sqrt(fourpi)
!    close(iofile)
!    where(rmesh(:,ia) < 0.1d0) kxc = 0.d0
!   kernel in cartesian labels
!   kernel from cartesian to spin labels
    do s4=1,nsmax
    do s3=1,nsmax
      do s2=1,nsmax
      do s1=1,nsmax
        if (sum(abs(kerspin(1:nr,s1,s2,s3,s4))) > atol) then
          do p4=1,iwsusc(2,s4,ia)
          do p3=1,iwsusc(2,s3,ia)
          do p2=1,iwsusc(2,s2,ia)
          do p1=1,iwsusc(2,s1,ia)
          work(1:nr) = phiref(1:nr,p1,2,s1,ia)*phiref(1:nr,p2,2,s2,ia)*phiref(1:nr,p3,2,s3,ia)*phiref(1:nr,p4,2,s4,ia)
          work(1:nr) = work(1:nr)*kerspin(1:nr,s1,s2,s3,s4)
          norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
          write(iodb,'(4es16.8)') br(nr,ia), nrv(nr,3,ia)
!          write(iodb,'(4es16.8)') phiref(nr,1,2,s1,ia), phiref(nr,1,2,s2,ia), phiref(nr,1,2,s3,ia), phiref(nr,1,2,s4,ia)
          write(iodb,'("U in eV:", 8i8,2f16.8)') s1, s2, s3, s4, p1, p2, p3, p4, norm*13.61d0
          end do
          end do
          end do
          end do
!          open(file="kerxc.dat",unit=400)
!          do i=1,nr
!            write(400,'(6es16.8)') rmesh(i,ia), phiref(i,1,2,1,ia), phiref(i,1,2,2,ia), br(i,ia), nrv(i,3,ia), kxc(i)
!          end do
!          close(400)
        end if
      end do
      end do
    end do
    end do
!   for my amusement: non-zero sums of gaunt products
    do lm4=5,9
    do lm3=5,9
      do lm2=5,9
      do lm1=5,9
        maxelem = 0.d0
        do ilm=1,lmmax2
          maxelem = maxelem + gaunt(lm1,lm2,ilm)*gaunt(lm3,lm4,ilm)
        end do
        if (abs(maxelem) > ylmtol) write(iodb,'(4i4,f10.3)') lm1, lm2, lm3, lm4, maxelem
      end do
      end do
    end do
    end do
    do jq=1+sum(nalmsbgf(1:ia-1)),sum(nalmsbgf(1:ia))
      i3 = i2almsbgf(:,jq)
      q3 = i3(1); q4 = i3(2)!; ia = i3(3)
      i3 = i2lmsb(:,q4,ia)
      p4 = i3(1); lm4 = i3(2); s4 = i3(3)
      i3 = i2lmsb(:,q3,ia)
      p3 = i3(1); lm3 = i3(2); s3 = i3(3)
      do iq=1+sum(nalmsbgf(1:ia-1)),sum(nalmsbgf(1:ia))
        i3 = i2almsbgf(:,iq)
        q1 = i3(1); q2 = i3(2)!; ia = i3(3)
        i3 = i2lmsb(:,q2,ia)
        p2 = i3(1); lm2 = i3(2); s2 = i3(3)
        i3 = i2lmsb(:,q1,ia)
        p1 = i3(1); lm1 = i3(2); s1 = i3(3)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       matrix elements of four basis functions
        work(1:nr) = phiref(1:nr,p1,i2lm(2,lm1),s1,ia)*phiref(1:nr,p2,i2lm(2,lm2),s2,ia)
        work(1:nr) = work(1:nr)*phiref(1:nr,p3,i2lm(2,lm3),s3,ia)*phiref(1:nr,p4,i2lm(2,lm4),s4,ia)
        work(1:nr) = work(1:nr)*kerspin(1:nr,s1,s2,s3,s4)
        norm = radint(nr,work,dr,,npanat(ia),ircutat(:,ia))
!        write(iodb,'(4i8,2es16.8)') q1, q2, q3, q4, norm
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          i1max = q1; i2max = q2; i3max = q3; i4max = q4
        end if
!        if (abs(norm) > susctol*maxelem) then
        doublegaunt = 0.d0
        do jlm=1,lmmax2
          gauntj = gaunt(lm3,lm4,jlm)
          if (abs(gauntj) > ylmtol) then
            gaunti = gaunt(lm1,lm2,jlm)
            if (abs(gaunti) > ylmtol) then
!              write(*,*) "jlm, ilm", jlm, ilm
              doublegaunt = doublegaunt + gaunti*gauntj
              if (jlm <= lmmax0) then
                if (lcartesian) then
                  do j=1,4
                    do i=1,4
                      kxcylm(i,j,jlm,ia) = kxcylm(i,j,jlm,ia) + norm*gaunti*gauntj*ps2c(i,s1,s2)*dc2s(s3,s4,j)
                    end do
                  end do
                else
                  j = is2i(s3,s4)
                    i = is2i(s1,s2)
                      kxcylm(i,j,jlm,ia) = kxcylm(i,j,jlm,ia) + norm*gaunti*gauntj
                end if
              end if
            end if
          end if
        end do
        kernel(iq,jq) = norm*doublegaunt
        if (abs(doublegaunt) < ylmtol) kernel(iq,jq) = 0.d0
!        end if
!       **** Spherical part ****
        if (lm1 == lm2 .and. lm3 == lm4) then
          if (lcartesian) then
            do j=1,4
              do i=1,4
                kxcy00(i,j,lm1,lm3,ia) = kxcy00(i,j,lm1,lm3,ia) + norm*ps2c(i,s1,s2)*dc2s(s3,s4,j)
              end do
            end do
          else
            j = is2i(s3,s4)
              i = is2i(s1,s2)
                kxcy00(i,j,lm1,lm3,ia) = kxcy00(i,j,lm1,lm3,ia) + norm
          end if
        else
!       TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
!          kxcsusc(iq,jq) = 0.d0
!       TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
        end if
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end do
    end do
    write(iodb,'(" xc kernel, ia=",i4)') ia
    write(iodb,'("kxc norm:",12i4,2es16.8)') i2lmsb(:,i1max,ia), i2lmsb(:,i2max,ia), i2lmsb(:,i3max,ia), i2lmsb(:,i4max,ia), maxnorm
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  kxcsusc = real(kxcsusc)
! All done
  end subroutine build_kxcalda

