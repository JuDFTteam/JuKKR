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
  logical :: calc_cheby_sll = .false.                  !!always calculate irregular solution in Chebychev solver (even if not needed) (former: 'norllsll')
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
  logical :: symmetrize_Gmat = .false.                 !!use symmetrization [G(k) + G(-k)]/2 in k-point loop (former: 'symG(k)')
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
  logical :: write_pkkr_operators= .false.             !!for Fermi-surface output: calculate various operators in KKR basis. (former: 'OPERATOR')
  logical :: write_potential_tests = .false.           !!write potential at different steps in main2 to different files (former: 'vintrasp' and 'vpotout')
  logical :: write_rho2ns = .false.                    !!write array rho2ns into file out_rhoval (from main1c) and out_rhotot (from main2) (former: 'RHOVALTW' and 'RHOVALW')
  logical :: write_rhoq_input = .false.                !!write out files needed for rhoq module (Quasiparticle interference) (former: 'rhoqtest')
  logical :: write_tmat_file = .false.                 !!write t-matix to file (former: 'tmatfile')


  deleted = 'EigenV', 'SPARSE', 'WIRE', 'iso surf', 'wfct', 'EV', 'ND', 'WAIT', 
  changed  = 'NEWSOSOL', 'RLL-SLL', 'SUPRCELL', 'full inv', 'godfrin', 'MPIadapt', 'MPIatom', 'MPIenerg', 'NOSOC', 'alt mix', 'atptshft', 'spec mix', 'timings0', 'timings2', 

#ifdef CPP_MPI

contains

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

    call mpi_bcast(t_params%npan_log, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npan_eq, 1, mpi_integer, master, mpi_comm_world, ierr)

  end subroutine bcast_runoptions
#endif CPP_MPI
