!------------------------------------------------------------------------------------
!> Summary: Writes decimation potential-file `decimate.pot`
!> Author: V. Popescu
!> Writes decimation potential-file 'decimate.pot' to be later used 
!> for 2D systems with the `DECIMATE` option. Based on the host potentials, the 
!> single-site matrices of the host can be calculated directly on each particular 
!> energy-mesh
!------------------------------------------------------------------------------------
!> @note So far, only `SPHERICAL` case implemented
!> @endnote
!------------------------------------------------------------------------------------
module mod_outpothost
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Writes decimation potential-file `decimate.pot`
  !> Author: V. Popescu
  !> Category: potential, single-site, input-output, KKRhost
  !> Deprecated: False 
  !> Writes decimation potential-file 'decimate.pot' to be later used 
  !> for 2D systems with the `DECIMATE` option. Based on the host potentials, the 
  !> single-site matrices of the host can be calculated directly on each particular 
  !> energy-mesh
  !-------------------------------------------------------------------------------
  !> @note So far, only `SPHERICAL` case implemented
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine outpothost(alat,ins,krel,kmrot,nspin,naez,natyp,efermi,bravais,rbasis, &
    qmtet,qmphi,noq,kaoez,iqat,zat,conc,ipan,ircut,solver,soc,ctl,irws,rmt,rws,rr,  &
    drdi,visp,irshift,rmrel,drdirel,vtrel,btrel,lmaxd,natypd,naezd,ipand,irmd)

    use :: mod_version_info

    implicit none
    ! .. Input variables 
    integer, intent(in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent(in) :: irmd   !! Maximum number of radial points
    integer, intent(in) :: krel   ! Switch for non- (or scalar-) relativistic/relativistic (Dirac) program (0/1). Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent(in) :: naez   !! Number of atoms in unit cell
    integer, intent(in) :: lmaxd  !! Maximum l component in wave function expansion
    integer, intent(in) :: kmrot  !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
    integer, intent(in) :: naezd  !! Number of atoms in unit cell
    integer, intent(in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent(in) :: ipand
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent(in) :: natypd !! Number of kinds of atoms in unit cell
    real (kind=dp), intent(in) :: alat    !! Lattice constant in a.u.
    real (kind=dp), intent(in) :: efermi  !! Fermi energy
    character (len=10), intent(in) :: solver  !! Type of solver
    integer, dimension(naezd), intent(in) :: noq  !! Number of diff. atom types located
    integer, dimension(*), intent(in)     :: ipan !! Number of panels in non-MT-region
    integer, dimension(*), intent(in)     :: iqat !! The site on which an atom is located on a given site
    integer, dimension(*), intent(in)     :: irws !! R point at WS radius
    integer, dimension(*), intent(in)     :: irshift  !! shift of the REL radial mesh with respect no NREL
    integer, dimension(0:ipand, *), intent(in)  :: ircut  !! r points of panel borders
    integer, dimension(natypd, *), intent(in)   :: kaoez  !! Kind of atom at site in elem. cell
    real (kind=dp), dimension(*), intent(in) :: rmt !! Muffin-tin radius of true system
    real (kind=dp), dimension(*), intent(in) :: rws !! Wigner Seitz radius
    real (kind=dp), dimension(*), intent(in) :: zat !! Nuclear charge
    real (kind=dp), dimension(*), intent(in) :: conc  !! Concentration of a given atom
    real (kind=dp), dimension(*), intent(in) :: qmphi !! \( \phi \) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension(*), intent(in) :: qmtet !! \( \theta \) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension(irmd, *), intent(in)                :: rr !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension(krel*lmaxd+1, *), intent(in)        :: soc
    real (kind=dp), dimension(krel*lmaxd+1, *), intent(in)        :: ctl
    real (kind=dp), dimension(irmd, *), intent(in)                :: visp   !! Spherical part of the potential
    real (kind=dp), dimension(irmd, *), intent(in)                :: drdi   !! Derivative dr/di
    real (kind=dp), dimension(3, *), intent(in)                   :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
    real (kind=dp), dimension(3, 3), intent(in)                   :: bravais  !! Bravais lattice vectors
    real (kind=dp), dimension(irmd*krel+(1-krel), *), intent(in)  :: rmrel    !! radial mesh
    real (kind=dp), dimension(irmd*krel+(1-krel), *), intent(in)  :: vtrel    !! potential (spherical part)
    real (kind=dp), dimension(irmd*krel+(1-krel), *), intent(in)  :: btrel    !! magnetic field
    real (kind=dp), dimension(irmd*krel+(1-krel), *), intent(in)  :: drdirel  !! derivative of radial mesh
    ! .. Local variables
    integer :: i, iq, ir, is, int
    character (len=9), dimension(2) :: txtrel
    character (len=9), dimension(3) :: txtspin
    character (len=3), dimension (0:113) :: elemname
    ! ..
    data txtspin/'         ', 'spin UP  ', 'spin DOWN'/
    data txtrel/'(UP+DN)/2', '(UP-DN)/2'/
    ! .. 1     2     3     4    5     6     7     8     9     0
    data elemname/'Vac', 'H  ', 'He ', 'Li ', 'Be', 'B  ', 'C  ', 'N  ', 'O  ',     &
      'F  ', 'Ne ', 'Na ', 'Mg ', 'Al ', 'Si', 'P  ', 'S  ', 'Cl ', 'Ar ', 'K  ',   &
      'Ca ', 'Sc ', 'Ti ','V  ', 'Cr', 'Mn ', 'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ',    &
      'Ga ', 'Ge ', 'As ', 'Se', 'Br ', 'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ', 'Nb ',   &
      'Mo ', 'Tc ', 'Ru', 'Rh ', 'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ', 'Sb ', 'Te ',   &
      'I  ', 'Xe', 'Cs ', 'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ', 'Pm ', 'Sm ', 'Eu ',   &
      'Gd', 'Tb ', 'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ', 'Lu ', 'Hf ', 'Ta ', 'W ',    &
      'Re ', 'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ', 'Tl ', 'Pb ', 'Bi ', 'Po', 'At ',   &
      'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ', 'Pa ', 'U  ', 'Np ', 'Pu', 'Am ', 'Cm ',   &
      'Bk ', 'Cf ', 'Es ', 'Fm ', 'Md ', 'No ', 'Lr ', 'Rf', 'Db ', 'Sg ', 'Bh ',   &
      'Hs ', 'Mt ', 'Uun', 'Uuu', 'Uub', 'NoE'/
    ! ..
    write (1337, '(5X,A,A,/)') '< OUTPOTHOST > : ', 'creating decimate.pot file - host potential'

    open (37, file='decimate.pot', status='unknown')
    call version_print_header(37)
    write (37, fmt=*) 'Host structure and potential for decimation'
    write (37, fmt=100)
    write (37, fmt=140) krel, ins, nspin, kmrot
    write (37, fmt=150) naez, natyp, alat
    write (37, fmt=130) efermi
    write (37, fmt=110) bravais
    ! ----------------------------------------------------------------------
    ! here insert whatever is more needed for the structure (BZ,SYM etc),
    ! ref. system and so on for a full host calculation (1 iteration)
    ! ----------------------------------------------------------------------
    write (37, fmt=120)
    do iq = 1, naez
      write (37, fmt=160) iq, (rbasis(i,iq), i=1, 3)
    end do
    write (37, fmt=170)
    do iq = 1, naez
      write (37, fmt=180) iq, qmtet(iq), qmphi(iq), noq(iq), (kaoez(i,iq), i=1, noq(iq))
    end do
    if (krel==1) write (37, 210) solver
    write (37, 190)
    do i = 1, natyp
      write (37, fmt=200) i, zat(i), iqat(i), conc(i), irws(i), ipan(i), (ircut(iq,i), iq=0, ipan(i))
      if (krel==1) write (37, 220) soc(1, i), ctl(1, i)
    end do
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    do i = 1, natyp
      write (37, '(80("*"))')
      ir = int(zat(i))
      write (37, 230) i, elemname(ir), rmt(i), rws(i)
      if (krel==0) then
        write (37, '(4(A9,11X))') 'R MESH   ', 'DRDI     ', (txtspin(nspin+is-1), is=1, nspin)
        do ir = 1, irws(i)
          write (37, 240) rr(ir, i), drdi(ir, i), (visp(ir,(i-1)*natyp+is), is=1, nspin)
        end do
      else
        write (37, 250) irshift(i)
        write (37, '(4(A9,11X))') 'R MESH   ', 'DRDI     ', (txtrel(is), is=1, 2)
        do ir = 1, irws(i) - irshift(i)
          write (37, 240) rmrel(ir, i), drdirel(ir, i), vtrel(ir, i), btrel(ir, i)
        end do
      end if
    end do

    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    close (37)
    ! ..
100 format ('Vectors in lattice constant units')
110 format ('BRAVAIS ', /, 3f8.4, /, 3f8.4, /, 3f8.4)
120 format ('RBASIS')
130 format ('EFERMI=', f10.6)
140 format ('KREL  =', i3, ' INS  =', i3, ' NSPIN=', i3, ' KMROT=', i3)
150 format ('NAEZ  =', i3, ' NATYP=', i3, ' ALAT =', f12.8)
160 format ('SITE  :', i3, 3f12.8)
170 format ('Magnetisation angles, occupancies, types on each site')
180 format ('SITE  :', i3, ' THETA=', f9.4, ' PHI  =', f9.4, ' NOQ  =', i3, ' ITOQ :', 8i3)
190 format ('ATOMS')
200 format ('TYPE  :', i3, ' Z    =', f4.0, ' IQAT =', i3, ' CONC =', f7.4, /, 10x, ' IRWS =', i4, ' IPAN =', i3, ' IRCUT=', 6i4)
210 format ('SOLVER=', a10)
220 format (10x, ' SOC  =', f10.6, ' CTL  =', d13.6)
230 format ('ATOM  :', i3, 1x, a3, ': mesh and potential data', /, 'RMT   :', f12.8, /, 'RWS   :', f12.8)
240 format (1p, 4d20.12)
250 format ('ISHIFT:', i3)
  end subroutine outpothost

end module mod_outpothost
