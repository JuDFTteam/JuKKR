!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Runoptions to control code behaviour
!> Author: Bernd Zimmermann 
!> Runoptions to control code behaviour
!------------------------------------------------------------------------------------
module mod_runoptions

  implicit none

  public

  logical :: calc_DOS_Efermi = .false.                 !!calculate DOS at Fermi energy only (former: 'DOS-EF')
  logical :: calc_GF_Efermi = .false.                  !!calculation of cluster Green function at E Fermi (former: 'GF-EF')
  logical :: set_cheby_nospeedup = .false.             !!always calculate irregular solution in Chebychev solver (even if not needed) (former: 'norllsll')
  logical :: set_cheby_nosoc     = .false.             !!decouple matrices in Chebychev solver neglecting SOC (former: 'NOSOC')
  logical :: calc_complex_bandstructure = .false.      !!complex band structure (former: 'COMPLEX')
  logical :: calc_exchange_couplings = .false.         !!calculate magnetic exchange coupling parameters (former: 'XCPL')
  logical :: calc_exchange_couplings_energy = .false.  !!write energy-resolved Jij-files also if npol/=0 (former: 'Jijenerg')
  logical :: calc_gmat_lm_full = .false.               !!calculate all lm-lm components of systems greens function and store to file `gflle` (former: 'lmlm-dos')
  logical :: dirac_scale_SpeefOfLight = .false.        !!scale the speed of light for Dirac solver (former: 'CSCALE')
  logical :: disable_charge_neutrality = .false.       !!no charge neutrailty required: leaving Fermi level unaffected (former: 'no-neutr')
  logical :: disable_print_serialnumber = .false.      !!deactivate writing of serial number and version information to files (for backwards compatibility) (former: 'noserial')
  logical :: disable_reference_system = .false.        !!deactivate the tight-binding reference system (former: 'lrefsysf')
  logical :: disable_tmat_sratrick = .false.           !!deactivate SRATRICK in solver for t-matirx (former: 'nosph')
  logical :: fix_nonco_angles = .false.                !!fix direction of non-collinear magnetic moments (Chebychev solver) (former: 'FIXMOM')
  logical :: formatted_files = .false.                 !!write files ascii-format. only effective with some other write-options (former: 'fileverb')
  logical :: impurity_operator_only = .false.          !!only for `write_pkkr_operators`: disable costly recalculation of host operators (former: 'IMP_ONLY')
  logical :: modify_soc_Dirac = .false.                !!modify SOC for Dirac solver (former: 'SOC')
  logical :: no_madelung = .false.                     !!do not add some energy terms (coulomb, XC, eff. pot.) to total energy (former: 'NoMadel')
  logical :: print_Gij = .false.                       !!print cluster G_ij matrices to outfile (former: 'Gmatij')
  logical :: print_gmat = .false.                      !!print Gmat to outfile (former: 'Gmat')
  logical :: print_ickeck = .false.                    !!enable test-output of ICHECK matrix from gfmask (former: 'ICHECK')
  logical :: print_kmesh = .false.                     !!output of k-mesh (former: 'k-net')
  logical :: print_kpoints = .false.                   !!print k-points to outfile (former: 'BZKP')
  logical :: print_program_flow = .false.              !!monitor the program flow in some parts of the code (former: 'flow')
  logical :: print_radial_mesh = .false.               !!write mesh information to output (former: 'RMESH')
  logical :: print_refpot = .false.                    !!test output of refpot (former: 'REFPOT')
  logical :: print_tau_structure = .false.             !!write extensive information about k-mesh symmetrization and structure of site-diagonal tau matrices to output (former: 'TAUSTRUC')
  logical :: print_tmat = .false.                      !!print t-matrix to outfile (former: 'tmat')
  logical :: relax_SpinAngle_Dirac = .false.           !!relax the spin angle in a SCF calculation [only DIRAC mode] (former: 'ITERMDIR')
  logical :: search_Efermi = .false.                   !!modify convergence parameters to scan for fermi energy only (to reach charge neutrality). (former: 'SEARCHEF')
  logical :: set_gmat_to_zero = .false.                !!set GMAT=0 in evaluation of density (former: 'GMAT=0')
  logical :: set_empty_system = .false.                !!set potential and nuclear charge to zero (former: 'zeropot')
  logical :: set_kmesh_large = .false.                 !!set equal k-mesh (largest) for all energy points (former: 'fix mesh')
  logical :: set_kmesh_small = .false.                 !!set equal k-mesh (smallest) for all energy points (former: 'fix4mesh')
  logical :: set_tmat_noinversion = .false.            !!do not perform inversion to get msst = Delta t^-1, but msst = Delta t. (former: 'testgmat')
  logical :: simulate_asa = .false.                    !!set non-spherical potential to zero in full-potential calculation with Chebychev solver (former: 'simulasa')
  logical :: slow_mixing_Efermi = .false.              !!renormalize Fermi-energy shift by mixing factor during mixing (former: 'slow-neu')
  logical :: stop_1a = .false.                         !!stop after main1a (former: 'STOP1A')
  logical :: stop_1b = .false.                         !!stop after main1b (former: 'STOP1B')
  logical :: stop_1c = .false.                         !!stop after main1c (former: 'STOP1C')
  logical :: symmetrize_gmat = .false.                 !!use symmetrization [G(k) + G(-k)]/2 in k-point loop (former: 'symG(k)')
  logical :: symmetrize_potential_cubic= .false.       !!keep only symmetric part of potential (L=1,11,21,25,43,47). (former: 'potcubic')
  logical :: symmetrize_potential_madelung = .false.   !!symmetrize potential in consistency to madelung potential (former: 'potsymm')
  logical :: torque_operator_onlyMT = .false.          !!for torque operator: include only the part within the muffin tin (former: 'ONLYMT')
  logical :: torque_operator_onlySph = .false.         !!for torque operator: include only the spherically symmetric part (former: 'ONLYSPH')
  logical :: use_BdG = .false.                         !!use Bogoliubov-de-Gennes Formalism (former: 'useBdG')
  logical :: use_Chebychev_solver = .false.            !!use the Chebychev solver (former: 'NEWSOSOL')
  logical :: use_cond_LB = .false.                     !!perform calculation of conductance in Landauer-Büttiker formalism (former: 'CONDUCT')
  logical :: use_cont = .false.                        !!no usage of embedding points. NEMB is set to 0. (former: 'CONT')
  logical :: use_deci_onebulk= .false.                 !!in case of decimation: use same bulk on right and left. Speeds up calculations. (former: 'ONEBULK')
  logical :: use_decimation = .false.                  !!use Decimation technique for semi-infinite systems (former: 'DECIMATE')
  logical :: use_ewald_2d = .false.                    !!use 2D ewald sum instead of 3D sum (Attention: does not work always!) (former: 'ewald2d')
  logical :: use_full_BZ = .false.                     !!use full Brillouin zone, i.e. switch off symmetries for k-space integration (former: 'fullBZ')
  logical :: use_ldau = .false.                        !!use LDA+U as exchange-correlation potential (former: 'LDA+U')
  logical :: use_lloyd = .false.                       !!use Lloyds formula to correct finite angular momentum cutoff (former: 'LLOYD')
  logical :: use_qdos = .false.                        !!writes out qdos files for band structure calculations. (former: 'qdos')
  logical :: use_readcpa = .false.                     !!read cpa t-matrix from file (former: 'readcpa')
  logical :: use_rigid_Efermi = .false.                !!keep the Fermi energy fixed during self-consistency (former: 'rigid-ef')
  logical :: use_semicore = .false.                    !!use semicore contour (former: 'SEMICORE')
  logical :: use_spherical_potential_only = .false.    !!keeping only spherical component of potential (former: 'Vspher')
  logical :: use_virtual_atoms = .false.               !!add virtual atoms (former: 'VIRATOMS')
  logical :: write_BdG_tests = .false.                 !!test options for Bogouliubov-deGennes (former: 'BdG_dev')
  logical :: write_DOS = .false.                       !!write out DOS files in any case (also if npol!=0) (former: 'DOS')
  logical :: write_DOS_lm = .false.                    !!write out DOS files with decomposition into l and m components (former: 'lmdos')
  logical :: write_gmat_plain = .false.                !!write out Green function as plain text file (former: 'GPLAIN')
  logical :: write_green_host = .false.                !!write green function of the host to file `green_host` (former: 'WRTGREEN')
  logical :: write_green_imp = .false.                 !!write out impurity Green function to GMATLL_GES (former: 'GREENIMP')
  logical :: write_complex_qdos = .false.              !!write complex qdos to file (former: 'compqdos')
  logical :: write_cpa_projection_files = .false.      !!write CPA projectors to file (former: 'projfile')
  logical :: write_deci_pot = .false.                  !!write decimation-potential file (former: 'deci-pot')
  logical :: write_deci_tmat = .false.                 !!write t-matrix to file 'decifile' (former: 'deci-out')
  logical :: write_density_ascii = .false.             !!write density rho2ns to file densitydn.ascii (former: 'den-asci')
  logical :: write_energy_mesh = .false.               !!write out the energy mesh to file `emesh.scf` (former: 'EMESH')
  logical :: write_generalized_potential = .false.     !!write potential in general format. Usually prepares for running the VORONOI program. (former: 'GENPOT')
  logical :: write_gmat_file = .false.                 !!write GMAT to file (former: 'gmatfile')
  logical :: write_gref_file = .false.                 !!write GREF to file (former: 'greffile')
  logical :: write_gmat_ascii = .false.                !!write GMAT to formatted file `gmat.ascii` (former: 'gmatasci')
  logical :: write_kkrimp_input = .false.              !!write out files for KKRimp-code (former: 'KKRFLEX')
  logical :: write_kkrsusc_input = .false.             !!write out files for KKRsusc-code (former: 'KKRSUSC')
  logical :: write_kpts_file = .false.                 !!write and read k-mesh to/from file `kpoints` (former: 'kptsfile')
  logical :: write_lloyd_cdos_file = .false.           !!write Lloyd array to file  (former: 'wrtcdos')
  logical :: write_lloyd_dgref_file= .false.           !!write Lloyd array to file  (former: 'wrtdgref')
  logical :: write_lloyd_dtmat_file= .false.           !!write Lloyd array to file  (former: 'wrtdtmat')
  logical :: write_lloyd_files = .false.               !!write several Lloyd-arrays to files (former: 'llyfiles')
  logical :: write_lloyd_g0tr_file = .false.           !!write Lloyd array to file  (former: 'wrtgotr')
  logical :: write_lloyd_tralpha_file= .false.         !!write Lloyd array to file  (former: 'wrttral')
  logical :: write_madelung_file = .false.             !!write madelung summation to file 'abvmad.unformatted' instead of keeping it in memory (former: 'madelfil')
  logical :: write_pkkr_input = .false.                !!write out files for Pkkprime-code (former: 'FERMIOUT')
  logical :: write_pkkr_operators = .false.            !!for Fermi-surface output: calculate various operators in KKR basis. (former: 'OPERATOR')
  logical :: write_potential_tests = .false.           !!write potential at different steps in main2 to different files (former: 'vintrasp' and 'vpotout')
  logical :: write_rho2ns = .false.                    !!write array rho2ns into file out_rhoval (from main1c) and out_rhotot (from main2) (former: 'RHOVALTW' and 'RHOVALW')
  logical :: write_rhoq_input = .false.                !!write out files needed for rhoq module (Quasiparticle interference) (former: 'rhoqtest')
  logical :: write_tmat_file = .false.                 !!write t-matix to file (former: 'tmatfile')
  logical :: write_tb_coupling = .false.               !!write couplings in tight-binging reference system to file `couplings.dat` (former: 'godfrin')

  !some old run and test options have been removed:
  !  'atptshft': replaced by presence or absence of IVSHIFT in inputcard
  !  'RLL-SLL': now default. See inverse (new) option 'set_cheby_nospeedup'
  !  'full inv', 'SUPRCELL', 'godfrin': now controlled via keyword <INVMODE> = 0/2/3 (<INVMODE>=1 signals principal-layer technique; defaults are set so that <INVMODE>=0 for bulk and  <INVMODE>=1 for slab systems)
  !  'timings0', 'timings2', 'verbose1', 'verbose2': incorporated into <VERBOSITY> keyword: 0=old behaviour (default), 1=low, 2=medium, 3=high
  !  'MPIadapt' and 'MPIatom' and 'MPIenerg': combined into keyword <MPI_SCHEME>: 0 = best, 1 = atoms (default), 2 = energies 
  !  deleted because of no implementation (or little use): 'EigenV', 'SPARSE', 'WIRE', 'iso surf', 'wfct', 'EV', 'ND', 'WAIT'

  contains

  !-------------------------------------------------------------------------------
  !> Summary: Convenience function to set a runoption
  !> Author: Bernd Zimmermann
  !> Category: KKRhost
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Read and set runoptions (new style as <keyword>= value)
  !-------------------------------------------------------------------------------
  subroutine read_runoptions()

    implicit none

    call set_runoption(formatted_files               , '<formatted_files>'               , '<fileverb>')
    call set_runoption(use_ldau                      , '<use_ldau>'                      , '<LDA+U>'   )
    call set_runoption(set_kmesh_large               , '<set_kmesh_large>'               , '<fix mesh>')
    call set_runoption(write_madelung_file           , '<write_madelung_file>'           , '<madelfil>')
    call set_runoption(write_BdG_tests               , '<write_BdG_tests>'               , '<BdG_dev>' )
    call set_runoption(use_virtual_atoms             , '<use_virtual_atoms>'             , '<VIRATOMS>')
    call set_runoption(write_gmat_ascii              , '<write_gmat_ascii>'              , '<gmatasci>')
    call set_runoption(use_rigid_Efermi              , '<use_rigid_Efermi>'              , '<rigid-ef>')
    call set_runoption(use_Chebychev_solver          , '<use_Chebychev_solver>'          , '<NEWSOSOL>')
    call set_runoption(write_pkkr_input              , '<write_pkkr_input>'              , '<FERMIOUT>')
    call set_runoption(calc_complex_bandstructure    , '<calc_complex_bandstructure>'    , '<COMPLEX>' )
    call set_runoption(write_pkkr_operators          , '<write_pkkr_operators>'          , '<OPERATOR>')
    call set_runoption(print_ickeck                  , '<print_ickeck>'                  , '<ICHECK>'  )
    call set_runoption(print_Gij                     , '<print_Gij>'                     , '<Gmatij>'  )
    call set_runoption(modify_soc_Dirac              , '<modify_soc_Dirac>'              , '<SOC>'     )
    call set_runoption(write_lloyd_tralpha_file      , '<write_lloyd_tralpha_file>'      , '<wrttral>' )
    call set_runoption(write_lloyd_cdos_file         , '<write_lloyd_cdos_file>'         , '<wrtcdos>' )
    call set_runoption(calc_gmat_lm_full             , '<calc_gmat_lm_full>'             , '<lmlm-dos>')
    call set_runoption(simulate_asa                  , '<simulate_asa>'                  , '<simulasa>')
    call set_runoption(use_readcpa                   , '<use_readcpa>'                   , '<readcpa>' )
    call set_runoption(print_kpoints                 , '<print_kpoints>'                 , '<BZKP>'    )
    call set_runoption(use_cont                      , '<use_cont>'                      , '<CONT>'    )
    call set_runoption(print_tmat                    , '<print_tmat>'                    , '<tmat>'    )
    call set_runoption(use_BdG                       , '<use_BdG>'                       , '<useBdG>'  )
    call set_runoption(disable_reference_system      , '<disable_reference_system>'      , '<lrefsysf>')
    call set_runoption(relax_SpinAngle_Dirac         , '<relax_SpinAngle_Dirac>'         , '<ITERMDIR>')
    call set_runoption(write_rho2ns                  , '<write_rho2ns>'                  , '<RHOVALTW>', '<RHOVALW>' )
    call set_runoption(write_generalized_potential   , '<write_generalized_potential>'   , '<GENPOT>'  )
    call set_runoption(print_gmat                    , '<print_gmat>'                    , '<Gmat>'    )
    call set_runoption(write_lloyd_dtmat_file        , '<write_lloyd_dtmat_file>'        , '<wrtdtmat>')
    call set_runoption(write_deci_pot                , '<write_deci_pot>'                , '<deci-pot>')
    call set_runoption(torque_operator_onlySph       , '<torque_operator_onlySph>'       , '<ONLYSPH>' )
    call set_runoption(print_tau_structure           , '<print_tau_structure>'           , '<TAUSTRUC>')
    call set_runoption(write_complex_qdos            , '<write_complex_qdos>'            , '<compqdos>')
    call set_runoption(use_full_BZ                   , '<use_full_BZ>'                   , '<fullBZ>'  )
    call set_runoption(use_semicore                  , '<use_semicore>'                  , '<SEMICORE>')
    call set_runoption(set_gmat_to_zero              , '<set_gmat_to_zero>'              , '<GMAT=0>'  )
    call set_runoption(use_decimation                , '<use_decimation>'                , '<DECIMATE>')
    call set_runoption(write_kkrimp_input            , '<write_kkrimp_input>'            , '<KKRFLEX>' )
    call set_runoption(write_green_imp               , '<write_green_imp>'               , '<GREENIMP>')
    call set_runoption(use_cond_LB                   , '<use_cond_LB>'                   , '<CONDUCT>' )
    call set_runoption(write_potential_tests         , '<write_potential_tests>'         , '<vintrasp>', '<vpotout>' )
    call set_runoption(disable_print_serialnumber    , '<disable_print_serialnumber>'    , '<noserial>')
    call set_runoption(symmetrize_gmat               , '<symmetrize_gmat>'               , '<symG(k)>' )
    call set_runoption(write_density_ascii           , '<write_density_ascii>'           , '<den-asci>')
    call set_runoption(write_rhoq_input              , '<write_rhoq_input>'              , '<rhoqtest>')
    call set_runoption(use_ewald_2d                  , '<use_ewald_2d>'                  , '<ewald2d>' )
    call set_runoption(write_DOS                     , '<write_DOS>'                     , '<DOS>'     )
    call set_runoption(write_energy_mesh             , '<write_energy_mesh>'             , '<EMESH>'   )
    call set_runoption(dirac_scale_SpeefOfLight      , '<dirac_scale_SpeefOfLight>'      , '<CSCALE>'  )
    call set_runoption(write_kpts_file               , '<write_kpts_file>'               , '<kptsfile>')
    call set_runoption(slow_mixing_Efermi            , '<slow_mixing_Efermi>'            , '<slow-neu>')
    call set_runoption(use_deci_onebulk              , '<use_deci_onebulk>'              , '<ONEBULK>' )
    call set_runoption(write_gmat_plain              , '<write_gmat_plain>'              , '<GPLAIN>'  )
    call set_runoption(disable_charge_neutrality     , '<disable_charge_neutrality>'     , '<no-neutr>')
    call set_runoption(stop_1c                       , '<stop_1c>'                       , '<STOP1C>'  )
    call set_runoption(stop_1b                       , '<stop_1b>'                       , '<STOP1B>'  )
    call set_runoption(stop_1a                       , '<stop_1a>'                       , '<STOP1A>'  )
    call set_runoption(print_kmesh                   , '<print_kmesh>'                   , '<k-net>'   )
    call set_runoption(fix_nonco_angles              , '<fix_nonco_angles>'              , '<FIXMOM>'  )
    call set_runoption(set_kmesh_small               , '<set_kmesh_small>'               , '<fix4mesh>')
    call set_runoption(disable_tmat_sratrick         , '<disable_tmat_sratrick>'         , '<nosph>'   )
    call set_runoption(print_refpot                  , '<print_refpot>'                  , '<REFPOT>'  )
    call set_runoption(write_lloyd_files             , '<write_lloyd_files>'             , '<llyfiles>')
    call set_runoption(symmetrize_potential_madelung , '<symmetrize_potential_madelung>' , '<potsymm>' )
    call set_runoption(set_cheby_nosoc               , '<set_cheby_nosoc>'               , '<NOSOC>'   )
    call set_runoption(set_empty_system              , '<set_empty_system>'              , '<zeropot>' )
    call set_runoption(write_tmat_file               , '<write_tmat_file>'               , '<tmatfile>')
    call set_runoption(write_tb_coupling             , '<write_tb_coupling>'             , '<godfrin>' )
    call set_runoption(symmetrize_potential_cubic    , '<symmetrize_potential_cubic>'    , '<potcubic>')
    call set_runoption(print_radial_mesh             , '<print_radial_mesh>'             , '<RMESH>'   )
    call set_runoption(impurity_operator_only        , '<impurity_operator_only>'        , '<IMP_ONLY>')
    call set_runoption(write_gref_file               , '<write_gref_file>'               , '<greffile>')
    call set_runoption(write_green_host              , '<write_green_host>'              , '<WRTGREEN>')
    call set_runoption(no_madelung                   , '<no_madelung>'                   , '<NoMadel>' )
    call set_runoption(set_tmat_noinversion          , '<set_tmat_noinversion>'          , '<testgmat>')
    call set_runoption(calc_exchange_couplings_energy, '<calc_exchange_couplings_energy>', '<Jijenerg>')
    call set_runoption(calc_GF_Efermi                , '<calc_GF_Efermi>'                , '<GF-EF>'   )
    call set_runoption(write_cpa_projection_files    , '<write_cpa_projection_files>'    , '<projfile>')
    call set_runoption(write_kkrsusc_input           , '<write_kkrsusc_input>'           , '<KKRSUSC>' )
    call set_runoption(calc_exchange_couplings       , '<calc_exchange_couplings>'       , '<XCPL>'    )
    call set_runoption(print_program_flow            , '<print_program_flow>'            , '<flow>'    )
    call set_runoption(use_spherical_potential_only  , '<use_spherical_potential_only>'  , '<Vspher>'  )
    call set_runoption(search_Efermi                 , '<search_Efermi>'                 , '<SEARCHEF>')
    call set_runoption(set_cheby_nospeedup           , '<set_cheby_nospeedup>'           , '<norllsll>')
    call set_runoption(write_lloyd_g0tr_file         , '<write_lloyd_g0tr_file>'         , '<wrtgotr>' )
    call set_runoption(write_deci_tmat               , '<write_deci_tmat>'               , '<deci-out>')
    call set_runoption(calc_DOS_Efermi               , '<calc_DOS_Efermi>'               , '<DOS-EF>'  )
    call set_runoption(use_qdos                      , '<use_qdos>'                      , '<qdos>'    )
    call set_runoption(write_lloyd_dgref_file        , '<write_lloyd_dgref_file>'        , '<wrtdgref>')
    call set_runoption(torque_operator_onlyMT        , '<torque_operator_onlyMT>'        , '<ONLYMT>'  )
    call set_runoption(write_DOS_lm                  , '<write_DOS_lm>'                  , '<lmdos>'   )
    call set_runoption(use_lloyd                     , '<use_lloyd>'                     , '<LLOYD>'   )
    call set_runoption(write_gmat_file               , '<write_gmat_file>'               , '<gmatfile>')

    write (555,*) 'in read_runoptions:', calc_exchange_couplings

  end subroutine read_runoptions


  !-------------------------------------------------------------------------------
  !> Summary: Convenience function to set a runoption
  !> Author: Bernd Zimmermann
  !> Category: KKRhost
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Convenience function to set a runoptions. New style as <Keyword>= Value is used.
  !> In addition to a default keyword, up to two alternative versions can be given to set the runoption.
  !-------------------------------------------------------------------------------
  subroutine set_runoption(runop,key,key_alt1,key_alt2)
    use mod_ioinput, only: ioinput, convert_to_uppercase

    implicit none

    logical,          intent(inout)        :: runop    !!runoption to be set
    character(len=*), intent(in)           :: key      !!corresponding keyword in inputcard
    character(len=*), intent(in), optional :: key_alt1 !!alternative keyword in inputcard
    character(len=*), intent(in), optional :: key_alt2 !!alternative keyword in inputcard

    character(len=:), allocatable :: uio
    character(len=:), allocatable :: key_upper !! keyword 
    integer :: ier
    integer, parameter :: ifile = 7

    key_upper = convert_to_uppercase(trim(adjustl(key)))

    call ioinput(key, uio, 1, ifile, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) runop
      write (111, *) key_upper // '= ', runop

    else if(present(key_alt1)) then
      call ioinput(key_alt1, uio, 1, ifile, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) runop
        write (111, *) key_upper // '= ', runop

      else if(present(key_alt2)) then
        call ioinput(key_alt2, uio, 1, ifile, ier)
        if (ier==0) then
          read (unit=uio, fmt=*) runop
          write (111, *) key_upper // '= ', runop

        else
          write (111, *) 'Default ' // key_upper // '= ', runop
        end if
      else
        write (111, *) 'Default ' // key_upper // '= ', runop
      end if
    else
      write (111, *) 'Default ' // key_upper // '= ', runop
    end if

  end subroutine set_runoption


  !-------------------------------------------------------------------------------
  !> Summary: Sets the runoptions in case the old style (fixed format) is used 
  !> Author: Bernd Zimmermann
  !> Category: KKRhost
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Sets the runoptions in case the old style (fixed format) is used.
  !> New code behaviour! Case-inseneitivity is introduced.
  !-------------------------------------------------------------------------------
  subroutine set_old_runoption(keyword_in,invmod,verbosity,MPI_scheme)

    use :: mod_ioinput, only: convert_to_uppercase

    implicit none

    character (len=8), intent(in) :: keyword_in
    integer, intent(inout) :: invmod
    integer, intent(inout) :: verbosity
    integer, intent(inout) :: MPI_Scheme

    character (len=8) :: keyword

    !left-adjust and make uppercase for case-insensitive treatment
    keyword = convert_to_uppercase(ADJUSTL(keyword_in))

    if (keyword == '        ') then
      !do nothing

    !now follows a (long) list of straight forward replacements
    else if (keyword == 'FILEVERB') then
      formatted_files = .true.
      write (1337, *) "    Enable runoption 'formatted_files'"
    else if (keyword == 'LDA+U   ') then
      use_ldau = .true.
      write (1337, *) "    Enable runoption 'use_ldau'"
    else if (keyword == 'FIX MESH') then
      set_kmesh_large = .true.
      write (1337, *) "    Enable runoption 'set_kmesh_large'"
    else if (keyword == 'MADELFIL') then
      write_madelung_file = .true.
      write (1337, *) "    Enable runoption 'write_madelung_file'"
    else if (keyword == 'BDG_DEV ') then
      write_BdG_tests = .true.
      write (1337, *) "    Enable runoption 'write_BdG_tests'"
    else if (keyword == 'VIRATOMS') then
      use_virtual_atoms = .true.
      write (1337, *) "    Enable runoption 'use_virtual_atoms'"
    else if (keyword == 'GMATASCI') then
      write_gmat_ascii = .true.
      write (1337, *) "    Enable runoption 'write_gmat_ascii'"
    else if (keyword == 'RIGID-EF') then
      use_rigid_Efermi = .true.
      write (1337, *) "    Enable runoption 'use_rigid_Efermi'"
    else if (keyword == 'NEWSOSOL') then
      use_Chebychev_solver = .true.
      write (1337, *) "    Enable runoption 'use_Chebychev_solver'"
    else if (keyword == 'FERMIOUT') then
      write_pkkr_input = .true.
      write (1337, *) "    Enable runoption 'write_pkkr_input'"
    else if (keyword == 'COMPLEX ') then
      calc_complex_bandstructure = .true.
      write (1337, *) "    Enable runoption 'calc_complex_bandstructure'"
    else if (keyword == 'OPERATOR') then
      write_pkkr_operators = .true.
      write (1337, *) "    Enable runoption 'write_pkkr_operators'"
    else if (keyword == 'ICHECK  ') then
      print_ickeck = .true.
      write (1337, *) "    Enable runoption 'print_ickeck'"
    else if (keyword == 'GMATIJ  ') then
      print_Gij = .true.
      write (1337, *) "    Enable runoption 'print_Gij'"
    else if (keyword == 'SOC     ') then
      modify_soc_Dirac = .true.
      write (1337, *) "    Enable runoption 'modify_soc_Dirac'"
    else if (keyword == 'WRTTRAL ') then
      write_lloyd_tralpha_file = .true.
      write (1337, *) "    Enable runoption 'write_lloyd_tralpha_file'"
    else if (keyword == 'WRTCDOS ') then
      write_lloyd_cdos_file = .true.
      write (1337, *) "    Enable runoption 'write_lloyd_cdos_file'"
    else if (keyword == 'LMLM-DOS') then
      calc_gmat_lm_full = .true.
      write (1337, *) "    Enable runoption 'calc_gmat_lm_full'"
    else if (keyword == 'SIMULASA') then
      simulate_asa = .true.
      write (1337, *) "    Enable runoption 'simulate_asa'"
    else if (keyword == 'READCPA ') then
      use_readcpa = .true.
      write (1337, *) "    Enable runoption 'use_readcpa'"
    else if (keyword == 'BZKP    ') then
      print_kpoints = .true.
      write (1337, *) "    Enable runoption 'print_kpoints'"
    else if (keyword == 'CONT    ') then
      use_cont = .true.
      write (1337, *) "    Enable runoption 'use_cont'"
    else if (keyword == 'TMAT    ') then
      print_tmat = .true.
      write (1337, *) "    Enable runoption 'print_tmat'"
    else if (keyword == 'USEBDG  ') then
      use_BdG = .true.
      write (1337, *) "    Enable runoption 'use_BdG'"
    else if (keyword == 'LREFSYSF') then
      disable_reference_system = .true.
      write (1337, *) "    Enable runoption 'disable_reference_system'"
    else if (keyword == 'ITERMDIR') then
      relax_SpinAngle_Dirac = .true.
      write (1337, *) "    Enable runoption 'relax_SpinAngle_Dirac'"
    else if (keyword == 'RHOVALTW' .or. keyword == 'RHOVALW ') then
      write_rho2ns = .true.
      write (1337, *) "    Enable runoption 'write_rho2ns'"
    else if (keyword == 'GENPOT  ') then
      write_generalized_potential = .true.
      write (1337, *) "    Enable runoption 'write_generalized_potential'"
    else if (keyword == 'GMAT    ') then
      print_gmat = .true.
      write (1337, *) "    Enable runoption 'print_gmat'"
    else if (keyword == 'WRTDTMAT') then
      write_lloyd_dtmat_file = .true.
      write (1337, *) "    Enable runoption 'write_lloyd_dtmat_file'"
    else if (keyword == 'DECI-POT') then
      write_deci_pot = .true.
      write (1337, *) "    Enable runoption 'write_deci_pot'"
    else if (keyword == 'ONLYSPH ') then
      torque_operator_onlySph = .true.
      write (1337, *) "    Enable runoption 'torque_operator_onlySph'"
    else if (keyword == 'TAUSTRUC') then
      print_tau_structure = .true.
      write (1337, *) "    Enable runoption 'print_tau_structure'"
    else if (keyword == 'COMPQDOS') then
      write_complex_qdos = .true.
      write (1337, *) "    Enable runoption 'write_complex_qdos'"
    else if (keyword == 'FULLBZ  ') then
      use_full_BZ = .true.
      write (1337, *) "    Enable runoption 'use_full_BZ'"
    else if (keyword == 'SEMICORE') then
      use_semicore = .true.
      write (1337, *) "    Enable runoption 'use_semicore'"
    else if (keyword == 'GMAT=0  ') then
      set_gmat_to_zero = .true.
      write (1337, *) "    Enable runoption 'set_gmat_to_zero'"
    else if (keyword == 'DECIMATE') then
      use_decimation = .true.
      write (1337, *) "    Enable runoption 'use_decimation'"
    else if (keyword == 'KKRFLEX ') then
      write_kkrimp_input = .true.
      write (1337, *) "    Enable runoption 'write_kkrimp_input'"
    else if (keyword == 'GREENIMP') then
      write_green_imp = .true.
      write (1337, *) "    Enable runoption 'write_green_imp'"
    else if (keyword == 'CONDUCT ') then
      use_cond_LB = .true.
      write (1337, *) "    Enable runoption 'use_cond_LB'"
    else if (keyword == 'VINTRASP' .or. keyword == 'VPOTOUT ') then
      write_potential_tests = .true.
      write (1337, *) "    Enable runoption 'write_potential_tests'"
    else if (keyword == 'NOSERIAL') then
      disable_print_serialnumber = .true.
      write (1337, *) "    Enable runoption 'disable_print_serialnumber'"
    else if (keyword == 'SYMG(K) ') then
      symmetrize_gmat = .true.
      write (1337, *) "    Enable runoption 'symmetrize_gmat'"
    else if (keyword == 'DEN-ASCI') then
      write_density_ascii = .true.
      write (1337, *) "    Enable runoption 'write_density_ascii'"
    else if (keyword == 'RHOQTEST') then
      write_rhoq_input = .true.
      write (1337, *) "    Enable runoption 'write_rhoq_input'"
    else if (keyword == 'EWALD2D ') then
      use_ewald_2d = .true.
      write (1337, *) "    Enable runoption 'use_ewald_2d'"
    else if (keyword == 'DOS     ') then
      write_DOS = .true.
      write (1337, *) "    Enable runoption 'write_DOS'"
    else if (keyword == 'EMESH   ') then
      write_energy_mesh = .true.
      write (1337, *) "    Enable runoption 'write_energy_mesh'"
    else if (keyword == 'CSCALE  ') then
      dirac_scale_SpeefOfLight = .true.
      write (1337, *) "    Enable runoption 'dirac_scale_SpeefOfLight'"
    else if (keyword == 'KPTSFILE') then
      write_kpts_file = .true.
      write (1337, *) "    Enable runoption 'write_kpts_file'"
    else if (keyword == 'SLOW-NEU') then
      slow_mixing_Efermi = .true.
      write (1337, *) "    Enable runoption 'slow_mixing_Efermi'"
    else if (keyword == 'ONEBULK ') then
      use_deci_onebulk = .true.
      write (1337, *) "    Enable runoption 'use_deci_onebulk'"
    else if (keyword == 'GPLAIN  ') then
      write_gmat_plain = .true.
      write (1337, *) "    Enable runoption 'write_gmat_plain'"
    else if (keyword == 'NO-NEUTR') then
      disable_charge_neutrality = .true.
      write (1337, *) "    Enable runoption 'disable_charge_neutrality'"
    else if (keyword == 'STOP1C  ') then
      stop_1c = .true.
      write (1337, *) "    Enable runoption 'stop_1c'"
    else if (keyword == 'STOP1B  ') then
      stop_1b = .true.
      write (1337, *) "    Enable runoption 'stop_1b'"
    else if (keyword == 'STOP1A  ') then
      stop_1a = .true.
      write (1337, *) "    Enable runoption 'stop_1a'"
    else if (keyword == 'K-NET   ') then
      print_kmesh = .true.
      write (1337, *) "    Enable runoption 'print_kmesh'"
    else if (keyword == 'FIXMOM  ') then
      fix_nonco_angles = .true.
      write (1337, *) "    Enable runoption 'fix_nonco_angles'"
    else if (keyword == 'FIX4MESH') then
      set_kmesh_small = .true.
      write (1337, *) "    Enable runoption 'set_kmesh_small'"
    else if (keyword == 'NOSPH   ') then
      disable_tmat_sratrick = .true.
      write (1337, *) "    Enable runoption 'disable_tmat_sratrick'"
    else if (keyword == 'REFPOT  ') then
      print_refpot = .true.
      write (1337, *) "    Enable runoption 'print_refpot'"
    else if (keyword == 'LLYFILES') then
      write_lloyd_files = .true.
      write (1337, *) "    Enable runoption 'write_lloyd_files'"
    else if (keyword == 'POTSYMM ') then
      symmetrize_potential_madelung = .true.
      write (1337, *) "    Enable runoption 'symmetrize_potential_madelung'"
    else if (keyword == 'NOSOC   ') then
      set_cheby_nosoc = .true.
      write (1337, *) "    Enable runoption 'set_cheby_nosoc'"
    else if (keyword == 'ZEROPOT ') then
      set_empty_system = .true.
      write (1337, *) "    Enable runoption 'set_empty_system'"
    else if (keyword == 'TMATFILE') then
      write_tmat_file = .true.
      write (1337, *) "    Enable runoption 'write_tmat_file'"
    else if (keyword == 'GODFRIN ') then
      write_tb_coupling = .true.
      write (1337, *) "    Enable runoption 'write_tb_coupling'"
    else if (keyword == 'POTCUBIC') then
      symmetrize_potential_cubic = .true.
      write (1337, *) "    Enable runoption 'symmetrize_potential_cubic'"
    else if (keyword == 'RMESH   ') then
      print_radial_mesh = .true.
      write (1337, *) "    Enable runoption 'print_radial_mesh'"
    else if (keyword == 'IMP_ONLY') then
      impurity_operator_only = .true.
      write (1337, *) "    Enable runoption 'impurity_operator_only'"
    else if (keyword == 'GREFFILE') then
      write_gref_file = .true.
      write (1337, *) "    Enable runoption 'write_gref_file'"
    else if (keyword == 'WRTGREEN') then
      write_green_host = .true.
      write (1337, *) "    Enable runoption 'write_green_host'"
    else if (keyword == 'NOMADEL ') then
      no_madelung = .true.
      write (1337, *) "    Enable runoption 'no_madelung'"
    else if (keyword == 'TESTGMAT') then
      set_tmat_noinversion = .true.
      write (1337, *) "    Enable runoption 'set_tmat_noinversion'"
    else if (keyword == 'JIJENERG') then
      calc_exchange_couplings_energy = .true.
      write (1337, *) "    Enable runoption 'calc_exchange_couplings_energy'"
    else if (keyword == 'GF-EF   ') then
      calc_GF_Efermi = .true.
      write (1337, *) "    Enable runoption 'calc_GF_Efermi'"
    else if (keyword == 'PROJFILE') then
      write_cpa_projection_files = .true.
      write (1337, *) "    Enable runoption 'write_cpa_projection_files'"
    else if (keyword == 'KKRSUSC ') then
      write_kkrsusc_input = .true.
      write (1337, *) "    Enable runoption 'write_kkrsusc_input'"
    else if (keyword == 'XCPL    ') then
      calc_exchange_couplings = .true.
      write (1337, *) "    Enable runoption 'calc_exchange_couplings'"
    else if (keyword == 'FLOW    ') then
      print_program_flow = .true.
      write (1337, *) "    Enable runoption 'print_program_flow'"
    else if (keyword == 'VSPHER  ') then
      use_spherical_potential_only = .true.
      write (1337, *) "    Enable runoption 'use_spherical_potential_only'"
    else if (keyword == 'SEARCHEF') then
      search_Efermi = .true.
      write (1337, *) "    Enable runoption 'search_Efermi'"
    else if (keyword == 'NORLLSLL') then
      set_cheby_nospeedup = .true.
      write (1337, *) "    Enable runoption 'set_cheby_nospeedup'"
    else if (keyword == 'WRTGOTR ') then
      write_lloyd_g0tr_file = .true.
      write (1337, *) "    Enable runoption 'write_lloyd_g0tr_file'"
    else if (keyword == 'DECI-OUT') then
      write_deci_tmat = .true.
      write (1337, *) "    Enable runoption 'write_deci_tmat'"
    else if (keyword == 'DOS-EF  ') then
      calc_DOS_Efermi = .true.
      write (1337, *) "    Enable runoption 'calc_DOS_Efermi'"
    else if (keyword == 'QDOS    ') then
      use_qdos = .true.
      write (1337, *) "    Enable runoption 'use_qdos'"
    else if (keyword == 'WRTDGREF') then
      write_lloyd_dgref_file = .true.
      write (1337, *) "    Enable runoption 'write_lloyd_dgref_file'"
    else if (keyword == 'ONLYMT  ') then
      torque_operator_onlyMT = .true.
      write (1337, *) "    Enable runoption 'torque_operator_onlyMT'"
    else if (keyword == 'LMDOS   ') then
      write_DOS_lm = .true.
      write (1337, *) "    Enable runoption 'write_DOS_lm'"
    else if (keyword == 'LLOYD   ') then
      use_lloyd = .true.
      write (1337, *) "    Enable runoption 'use_lloyd'"
    else if (keyword == 'GMATFILE') then
      write_gmat_file = .true.
      write (1337, *) "    Enable runoption 'write_gmat_file'"

    !treat special replacements
    else if (keyword == 'ATPTSHFT') then
      write (1337, *) "    No need to call runoption atptshft any more, giving keyword IVSHIFT is sufficient"
    else if (keyword == 'RLL-SLL ') then
      write (1337, *) "    Runoption RLL-SLL is not default (has been removed)"
    else if (keyword == 'FULL INV') then
      invmod = 0
      write (1337, *) "    Setting inversion mode to Full Inversion"
    else if (keyword == 'SUPRCELL') then
      invmod = 2
      write (1337, *) "    Setting inversion mode to Supercell mode"
    else if (keyword == 'GODFRIN ') then
      invmod = 3
      write (1337, *) "    Setting inversion mode to Godfrin module"
    else if (keyword == 'VERBOSE1') then
      verbosity = 2
      write (1337, *) "    Setting verbosity level to medium because of verbose1"
    else if (keyword == 'VERBOSE2') then
      verbosity = 3
      write (1337, *) "    Setting verbosity level to high because of verbose2"
    else if (keyword == 'TIMINGS0') then
      verbosity = 1
      write (1337, *) "    Setting verbosity level to low because of timings0"
    else if (keyword == 'TIMINGS2') then
      verbosity = 3
      write (1337, *) "    Setting verbosity level to high because of timings2"
    else if (keyword == 'MPIADAPT') then
      MPI_scheme = 0
      write (1337, *) "    Setting MPI_Scheme=0 because of MPIadapt"
    else if (keyword == 'MPIATOM ') then
      MPI_scheme = 1
      write (1337, *) "    Setting MPI_Scheme=1 because of MPIatom"
    else if (keyword == 'MPIENERG') then
      MPI_scheme = 2
      write (1337, *) "    Setting MPI_Scheme=2 because of MPIenerg"
    else if (keyword == 'EIGENV  ' .or. keyword == 'SPARSE  ' .or. keyword == 'WIRE    ' .or. keyword == 'ISO SURF' .or. keyword == 'EIGENV  ' .or. keyword == 'WFCT    ' .or. keyword == 'EV      ' .or. keyword == 'ND      ' .or. keyword == 'WAIT    ') then
      write (1337, *) "    ### Ignoring option '" // keyword_in // "' because it is not implemented any more. ###"

    !default case: give a warning
    else
      write (1337, *) "    ### WARNING ###: run- or testoption '" // keyword_in // "' not recognized."
      write (1337, *) "                     Is there a typo?"
    end if


  end subroutine set_old_runoption


#ifdef CPP_MPI


  !-------------------------------------------------------------------------------
  !> Summary: Broadcasts the Runoptions after readin
  !> Author: Bernd Zimmermann
  !> Category: KKRhost
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Broadcasts the Runoptions after readin
  !-------------------------------------------------------------------------------
  subroutine bcast_runoptions()

    use :: mpi
    use :: mod_mympi, only: master

    implicit none
    integer :: ierr

    !call mpi_bcast(t_params%npan_log, 1, mpi_integer, master, mpi_comm_world, ierr)
    !call mpi_bcast(t_params%npan_eq, 1, mpi_integer, master, mpi_comm_world, ierr)

  end subroutine bcast_runoptions
#endif CPP_MPI

end module mod_runoptions