! 13.10.95 ***************************************************************
subroutine ecoub(cmom, ecou, lmax, nspin, natyp, rho2ns, vm2z, z, r, drdi, &
  irws, kvmad, kshape, ircut, ipan, imaxsh, ifunm, ilm_map, ntcell, gsh, &
  thetas, lmsp, lpot)
  ! ************************************************************************

  ! attention : energy zero ---> electro static zero

  ! calculate the electrostatic potential-energies without the
  ! electron-nuclear interaction in the cell itself .
  ! the energy of the representive atom i is given by

  ! rc
  ! ecou(i) =  1/2 (  {  s dr' vm2z(r',lm,i)*rho2ns(r',lm,i,1) }
  ! 0

  ! -  z(i) * vmad ( ri )     )


  ! ( {..} = summed over lm )
  ! (see notes by b.drittler)
  ! vm2z is the coulomb potential of the atom without the nuclear
  ! potential of the atom
  ! rho2ns(...,1) is the real charge density times r**2

  ! both developed into spherical harmonics . (see deck rholm)

  ! z    is the nuclear charge of the atom

  ! vmad ( ri ) is a generalized madelung potential
  ! = 1/sqrt(4 pi) * vm2z(irws,1,is)
  ! - sqrt(4 pi) * 2 * cmom(1,ipot) / rws

  ! ( <..> = spherical averaged )

  ! attention : this subroutine has to be called before the
  ! exchange correlation potential is added to
  ! the potential vm2z .
  ! the energy calculated here is splitted into
  ! l-dependent parts to see the l -convergency .

  ! attention : in case of shape corrections the contribution of
  ! the coulomb potential the of the nucleus is
  ! analytically cancelled only in the muffin tin sphere
  ! in the interstial region it has to be taken into
  ! account ! see deck epotins

  ! modified for band structure code
  ! b.drittler   jan. 1990
  ! -----------------------------------------------------------------------
  use :: mod_types, only: t_inc
  use :: mod_datatypes, only: dp
  use global_variables
  implicit none
  ! ..
  ! .. Local Scalars ..
  integer :: lpot, kshape, kvmad, lmax, natyp, nspin
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: cmom(lmpotd, *), drdi(irmd, *), ecou(0:lpot, *), gsh(*), &
    r(irmd, *), rho2ns(irmd, lmpotd, natypd, *), thetas(irid, nfund, *), &
    vm2z(irmd, lmpotd, *), z(*)
  integer :: ifunm(natypd, *), ilm_map(ngshd, 3), imaxsh(0:lmpotd), ipan(*), &
    ircut(0:ipand, *), irws(*), ntcell(*), lmsp(natypd, *)
  ! ..
  ! .. External Subroutines ..
  real (kind=dp) :: rfpi, rhosp, sign, vm, vmad
  integer :: i, iatyp, icell, ifun, ipan1, ipot, ir, irc1, irh, irs1, ispin, &
    j, l, lm, lm2, m
  ! ..
  ! .. Intrinsic Functions ..
  real (kind=dp) :: er(irmd)
  ! ..

  external :: simp3, simpk


  intrinsic :: atan, sqrt
  ! --->   determine the right potential numbers - the coulomb potential
  rfpi = sqrt(16.d0*atan(1.0d0))
  ! is not spin dependend
  do iatyp = 1, natyp

    if (kshape/=0) then
      ipan1 = ipan(iatyp)
      icell = ntcell(iatyp)
      irs1 = ircut(1, iatyp)
      irc1 = ircut(ipan1, iatyp)
    else
      irs1 = irws(iatyp)
      irc1 = irs1
    end if




    ipot = iatyp*nspin

    do l = 0, lmax

      do i = 1, irc1
        er(i) = 0.0d0
      end do

      do ispin = 1, nspin

        if (ispin==nspin) then
          sign = 1.0d0
        else
          sign = -1.0d0
        end if
        ! --->           convolute with shape function
        do m = -l, l
          lm = l*l + l + m + 1

          do i = 1, irs1
            rhosp = (rho2ns(i,lm,iatyp,1)+sign*rho2ns(i,lm,iatyp,nspin))/4.0d0
            er(i) = er(i) + rhosp*vm2z(i, lm, ipot)
          end do

          if (kshape/=0) then

            ! --->                 remember that in the interstial -2z/r has
            ! to be taken into account
            do j = imaxsh(lm-1) + 1, imaxsh(lm)
              lm2 = ilm_map(j, 2)
              if (lmsp(icell,ilm_map(j,3))>0) then
                ifun = ifunm(icell, ilm_map(j,3))

                if (lm2==1) then
                  do ir = irs1 + 1, irc1
                    irh = ir - irs1
                    rhosp = (rho2ns(ir,lm,iatyp,1)+sign*rho2ns(ir,lm,iatyp, &
                      nspin))/2.0d0

                    ! (LM2.EQ.1)


                    er(ir) = er(ir) + rhosp*gsh(j)*thetas(irh, ifun, icell)*( &
                      vm2z(ir,1,ipot)/2.0d0-z(iatyp)/r(ir,iatyp)*rfpi)
                  end do
                  ! (LM2.EQ.1)
                else

                  do ir = irs1 + 1, irc1
                    irh = ir - irs1
                    rhosp = (rho2ns(ir,lm,iatyp,1)+sign*rho2ns(ir,lm,iatyp, &
                      nspin))/2.0d0
                    er(ir) = er(ir) + rhosp*gsh(j)*thetas(irh, ifun, icell)* &
                      vm2z(ir, lm2, ipot)/2.0d0
                  end do
                  ! (KSHAPE.NE.0)
                end if
              end if

            end do

          end if                   ! --->     now integrate

        end do

      end do


      ! --->   calculate the madelung potential
      if (kshape==0) then
        call simp3(er, ecou(l,iatyp), 1, irs1, drdi(1,iatyp))
      else
        call simpk(er, ecou(l,iatyp), ipan1, ircut(0,iatyp), drdi(1,iatyp))
      end if

    end do

    ! --->   add to ecou


    vmad = vm2z(irs1, 1, ipot)/rfpi - rfpi*2.0d0*cmom(1, iatyp)/r(irs1, iatyp)
    ! --->   option to calculate full generalized madelung potential
    ! rc
    ! vm(rn) = vmad +2*sqrt(4*pi)* s  dr*r*rho(lm=1,r)
    ecou(0, iatyp) = ecou(0, iatyp) - z(iatyp)*vmad/2.0d0
    ! 0



    ! atom nr. iatyp is the iatyp-th atom on the potential cards
    if (kvmad==1) then
      er(1) = 0.0d0
      ! e. g., in binary alloys iatyp=1 and iatyp=2 refer to host
      do i = 2, irs1
        er(i) = rho2ns(i, 1, iatyp, 1)/r(i, iatyp)
      end do

      call simp3(er, vm, 1, irs1, drdi(1,iatyp))
      vm = 2.0d0*rfpi*vm + vmad
      ! (KVMAD.EQ.1)



      if (t_inc%i_write>0) write (1337, fmt=110) iatyp, vmad
      if (t_inc%i_write>0) write (1337, fmt=100) iatyp, vm
    end if                         ! 13.10.95
                                   ! ***************************************************************
    ! ************************************************************************
  end do

  return
  ! attention : energy zero ---> electro static zero
100 format (10x, 'full generalized madelung pot. for atom', 1x, i3, 1x, ': ', &
    1p, d14.6)
110 format (10x, '     generalized madelung pot. for atom', 1x, i3, 1x, ': ', &
    1p, d14.6)
end subroutine ecoub
