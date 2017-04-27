  subroutine basis_change()
! Takes the KS susceptibility from the product to the density basis
  use global

  implicit none


  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0)
  integer(kind=i4b) :: nden0, nden1, ngf0, ngf1
  integer(kind=i4b) :: ia, ja, iq, iq1, jq, jq1, ib, jb, ilm, jlm, is, js, i4(4)
  real(kind=r8b)    :: start, finish
  complex(kind=c8b) :: density(ndensum), denlms(lmmax0,nsmax2,nasusc)

  call cpu_time(start)
! multiply the KS susc from the left by dengaunt
  denominator = 0.d0
  do jq1=1,ngfsum
    do ia=1,nasusc
      nden0 = sum(nalmsbden(1:ia-1))
      nden1 = sum(nalmsbden(1:ia))
      do iq=1+nden0,nden1
        ngf0 = sum(nalmsbgf(1:ia-1))
        ngf1 = sum(nalmsbgf(1:ia))
        do iq1=1+ngf0,ngf1
          if (abs(dengaunt(iq1-ngf0,iq-nden0,ia)) > ylmtol) then
            denominator(iq,jq1) = denominator(iq,jq1) + dengaunt(iq1-ngf0,iq-nden0,ia)*kssusc0(iq1,jq1)
          end if
        end do
      end do
    end do
  end do
! multiply the KS susc from the right by dengaunt
  kssusc0 = 0.d0
  do ia=1,nasusc
    nden0 = sum(nalmsbden(1:ia-1))
    nden1 = sum(nalmsbden(1:ia))
    do jq=1+nden0,nden1
      ngf0 = sum(nalmsbgf(1:ia-1))
      ngf1 = sum(nalmsbgf(1:ia))
      do jq1=1+ngf0,ngf1
        if (abs(dengaunt(jq1-ngf0,jq-nden0,ia)) > ylmtol) then
          do iq=1,ndensum
            kssusc0(iq,jq) = kssusc0(iq,jq) + denominator(iq,jq1)*dengaunt(jq1-ngf0,jq-nden0,ia)
          end do
        end if
      end do
    end do
  end do
! density = susceptibility * potential
  density = 0.d0
  call zgemv('N',ndensum,ndensum,cone,kssusc0,ngfsum,vlmsbden,1,czero,density,1)
  denlms = 0.d0
  do iq=1,ndensum
    i4 = i2almsbden(:,iq)
    ib = i4(1); ilm = i4(2); is = i4(3); ia = i4(4)
    denlms(ilm,is,ia) = denlms(ilm,is,ia) + density(iq)*suscnorm(iq)
  end do
  do ia=1,nasusc
    do ilm=1,lmmax0
      write(iodb,'("basis_change: denlms=",8f16.8)') denlms(ilm,:,ia)
    end do
  end do
  call cpu_time(finish)
  write(iodb,'("basis change: time=",f12.3," s")') finish - start
! All done!
  end subroutine basis_change

