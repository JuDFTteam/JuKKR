SHELL = /bin/sh

# Compiler and options
FC = ifort
# Debugger flags 
# FFLAGS= -O0 -g -debug all -CB -check -traceback -heap-arrays -warn all
# FFLAGS= -O0 -g -debug all -CB -check uninit -traceback -heap-arrays -warn all -ftrapuv
# Usual flags  
FFLAGS= -O3 -traceback -heap-arrays -g
FDBFLAGS= -C -traceback -heap-arrays

# Clear suffix list
.SUFFIXES:
# Suffixes in use
.SUFFIXES: .F90 .f90 .f .o .mod

########################################################################
#                       Folders configuration
########################################################################
OBJ   = obj

MEMDIR   = memory
MEM      = $(OBJ)/global.o $(OBJ)/init_param.o $(OBJ)/finit_param.o \
           $(OBJ)/init_arrays.o $(OBJ)/finit_arrays.o $(OBJ)/init_basis.o \
           $(OBJ)/finit_basis.o $(OBJ)/init_gfpartsfit.o $(OBJ)/finit_gfpartsfit.o\
           $(OBJ)/init_arrays_tgmat.o  $(OBJ)/init_gfpartsfit_coeff_gfpq.o

BASISDIR = basis
BASIS    = $(OBJ)/save_rmesh.o $(OBJ)/read_rmesh.o $(OBJ)/ref_wfn.o \
           $(OBJ)/new_basis.o $(OBJ)/new_basis1.o $(OBJ)/new_basis2.o \
           $(OBJ)/new_basis3.o $(OBJ)/find_basis1.o $(OBJ)/find_basis2.o \
           $(OBJ)/out_overlap.o $(OBJ)/save_wfns.o $(OBJ)/read_wfns.o \
           $(OBJ)/overlaps_gf.o $(OBJ)/overlaps_susc2.o $(OBJ)/basis_change.o \
           $(OBJ)/gradient_susc_basis.o

UTILSDIR = utils
UTILS    = $(OBJ)/radint.o $(OBJ)/ymy_gaunts.o $(OBJ)/rymy.o \
           $(OBJ)/orb_moment.o $(OBJ)/dlasrt2_saved.o $(OBJ)/zratint.o \
           $(OBJ)/chebyweights.o $(OBJ)/llsetup.o $(OBJ)/LebedevLaikov.o \
           $(OBJ)/bessel_new.o $(OBJ)/free_gf.o $(OBJ)/ratval.o \
           $(OBJ)/ratroots.o $(OBJ)/least_squares.o $(OBJ)/ratval2.o \
           $(OBJ)/rotvec.o $(OBJ)/mod_derivative_panels.o $(OBJ)/ymy.o \
           $(OBJ)/spline.o $(OBJ)/mod_spline_panels2.o $(OBJ)/mod_spline_panels.o \
           $(OBJ)/splint.o $(OBJ)/zchop.o

PROJDIR  = projection
PROJ     = $(OBJ)/reg_coeffs.o $(OBJ)/irr_coeffs.o $(OBJ)/out_coeffs.o \
           $(OBJ)/in_coeffs.o $(OBJ)/save_coeffs.o $(OBJ)/out_gmat.o \
           $(OBJ)/out_tgmat.o $(OBJ)/in_tmat.o $(OBJ)/in_gmat.o \
           $(OBJ)/save_tmati.o $(OBJ)/save_gsij.o $(OBJ)/save_tmcoll.o \
           $(OBJ)/save_gscoll.o $(OBJ)/projected_gf.o $(OBJ)/baryint_gf2.o \
           $(OBJ)/ratint_gf2.o $(OBJ)/biratint_gf.o $(OBJ)/ratval_gf.o \
           $(OBJ)/ratfit_gf.o $(OBJ)/ratfit_gf2.o $(OBJ)/out_tgmat_soc.o \
           $(OBJ)/reg_coeffs_soc.o $(OBJ)/irr_coeffs_soc.o $(OBJ)/out_coeffs_soc.o \
           $(OBJ)/tmat_proj.o $(OBJ)/spinrot_scatt_sol.o $(OBJ)/in_coeffs_soc.o $(OBJ)/in_gmat_soc.o \
	   $(OBJ)/in_tmat_soc.o

GSDIR    = groundstate
GS       = $(OBJ)/save_rho2ns.o $(OBJ)/read_rho2ns.o $(OBJ)/update_rho2ns.o \
           $(OBJ)/update_espv.o $(OBJ)/save_rhomat.o $(OBJ)/read_rhomat.o \
           $(OBJ)/mix_rhomat.o $(OBJ)/groundstate_new.o $(OBJ)/get_gfsum.o \
           $(OBJ)/integrated_charge.o $(OBJ)/charge_density.o \
           $(OBJ)/density_of_states.o $(OBJ)/get_xc2.o $(OBJ)/vxc_vwn.o \
           $(OBJ)/density_matrix.o $(OBJ)/current_density.o \
           $(OBJ)/current_soc.o $(OBJ)/current_zeeman.o \
           $(OBJ)/current_interpolation.o $(OBJ)/rho_interpolation.o

CORRDIR  = corrections
CORR     = $(OBJ)/build_selfe.o $(OBJ)/build_zeeman.o $(OBJ)/build_vscfb.o \
           $(OBJ)/ms_from_bxc3.o $(OBJ)/build_vxcdiff.o $(OBJ)/build_vsocz.o \
           $(OBJ)/build_vsocpm.o $(OBJ)/build_vldaub.o $(OBJ)/build_vsoc.o \
           $(OBJ)/new_directions.o $(OBJ)/spin_directions.o $(OBJ)/spin_rotation.o \
           $(OBJ)/local_frame.o $(OBJ)/orb_rotation.o $(OBJ)/build_vrotb.o \
           $(OBJ)/spinrot_tgf.o $(OBJ)/spinrot_gf.o $(OBJ)/spinrot_gf_sph.o \
           $(OBJ)/spinrot_rhomat.o $(OBJ)/lippmann_schwinger.o $(OBJ)/structural_gf.o \
           $(OBJ)/torque_soc.o $(OBJ)/torque_soc2.o $(OBJ)/lichtenstein_jij.o \
           $(OBJ)/build_dti.o $(OBJ)/build_dti_v2.o $(OBJ)/anisotropic_jij.o \
           $(OBJ)/symmetrize_tgf.o $(OBJ)/pot_correction.o $(OBJ)/build_vsocb.o \
           $(OBJ)/spinrot_tgf_sph.o $(OBJ)/ms_from_bxc4.o $(OBJ)/build_vscfb2.o \
           $(OBJ)/build_bsoctwist.o $(OBJ)/build_vsocz2.o $(OBJ)/build_vsocpm2.o \
           $(OBJ)/build_vxcdiff2.o $(OBJ)/build_grad_mass.o

SUSCDIR  = susceptibility
SUSC     = $(OBJ)/dyn_susc_expansion.o $(OBJ)/dynamic_susc.o $(OBJ)/dynamic_susc2.o \
           $(OBJ)/static_susc.o $(OBJ)/static_susc2.o $(OBJ)/static_susc_more.o \
           $(OBJ)/static_susc_enhanced.o $(OBJ)/add_susc_block.o $(OBJ)/susc_denominator.o \
           $(OBJ)/build_kxcalda2.o $(OBJ)/kxc_sumrule.o $(OBJ)/build_khartree.o \
           $(OBJ)/analytic_continuation.o $(OBJ)/add_susc_block_one_sided.o \
           $(OBJ)/current_correlation.o $(OBJ)/current_correlation_interpolation.o \
           $(OBJ)/current_correlation_divergence_interpolation.o $(OBJ)/susc_interpolation.o \
           $(OBJ)/current_correlation_divergence.o $(OBJ)/build_twist.o $(OBJ)/kxc_sumrule2.o



IODIR    = io
IO       = $(OBJ)/new_input.o $(OBJ)/load_input.o $(OBJ)/find_keyinline.o \
           $(OBJ)/find_keyinfile.o $(OBJ)/set_basis_options.o \
           $(OBJ)/set_spinrot_options.o $(OBJ)/set_bfield_options.o $(OBJ)/set_soc_options.o \
           $(OBJ)/set_selfe_options.o $(OBJ)/set_dos_options.o $(OBJ)/set_jij_options.o \
           $(OBJ)/set_susc_options.o $(OBJ)/set_fit_options.o $(OBJ)/set_ldau_options.o \
           $(OBJ)/set_groundstate_options.o $(OBJ)/set_groups.o $(OBJ)/set_global_options.o \
           $(OBJ)/close_outsusc.o $(OBJ)/outsusc_out.o $(OBJ)/outsusc_in.o \
           $(OBJ)/outsusc_in2.o $(OBJ)/outsusc_in_soc.o

########################################################################
#                          Solver module     
########################################################################
solver: $(MEM) $(UTILS) $(BASIS) $(PROJ) $(GS) $(CORR) $(SUSC) $(IO)

########################################################################
#                              Global
########################################################################
$(OBJ)/global.o: $(MEMDIR)/global.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@ ; mv global.mod $(OBJ)/

########################################################################
#                          Memory module
########################################################################
$(OBJ)/init_param.o: $(MEMDIR)/init_param.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/finit_param.o: $(MEMDIR)/finit_param.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/init_arrays.o: $(MEMDIR)/init_arrays.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/finit_arrays.o: $(MEMDIR)/finit_arrays.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/init_basis.o: $(MEMDIR)/init_basis.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/finit_basis.o: $(MEMDIR)/finit_basis.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/init_gfpartsfit.o: $(MEMDIR)/init_gfpartsfit.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/finit_gfpartsfit.o: $(MEMDIR)/finit_gfpartsfit.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/init_arrays_tgmat.o: $(MEMDIR)/init_arrays_tgmat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/init_gfpartsfit_coeff_gfpq.o: $(MEMDIR)/init_gfpartsfit_coeff_gfpq.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                          Utils module
########################################################################
$(OBJ)/radint.o: $(UTILSDIR)/radint.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ymy_gaunts.o: $(UTILSDIR)/ymy_gaunts.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/rymy.o: $(UTILSDIR)/rymy.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ymy.o: $(UTILSDIR)/ymy.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/orb_moment.o: $(UTILSDIR)/orb_moment.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/dlasrt2_saved.o: $(UTILSDIR)/dlasrt2_saved.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/zratint.o: $(UTILSDIR)/zratint.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/chebyweights.o: $(UTILSDIR)/chebyweights.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/llsetup.o: $(UTILSDIR)/llsetup.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/LebedevLaikov.o: $(UTILSDIR)/LebedevLaikov.f
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/bessel_new.o: $(UTILSDIR)/bessel_new.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@ ; mv bessel_new.mod $(OBJ)/

$(OBJ)/free_gf.o: $(UTILSDIR)/free_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ratval.o: $(UTILSDIR)/ratval.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ratval2.o: $(UTILSDIR)/ratval2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ratroots.o: $(UTILSDIR)/ratroots.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/least_squares.o: $(UTILSDIR)/least_squares.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/rotvec.o: $(UTILSDIR)/rotvec.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/mod_derivative_panels.o: $(UTILSDIR)/mod_derivative_panels.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@; mv mod_derivative_panels.mod $(OBJ)/

$(OBJ)/spline.o: $(UTILSDIR)/spline.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/mod_spline_panels2.o: $(UTILSDIR)/mod_spline_panels2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@; mv mod_spline_panels2.mod $(OBJ)/

$(OBJ)/mod_spline_panels.o: $(UTILSDIR)/mod_spline_panels.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@; mv mod_spline_panels.mod $(OBJ)/

$(OBJ)/splint.o: $(UTILSDIR)/splint.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/zchop.o: $(UTILSDIR)/zchop.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                           Basis module
########################################################################
$(OBJ)/save_rmesh.o: $(BASISDIR)/save_rmesh.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/read_rmesh.o: $(BASISDIR)/read_rmesh.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ref_wfn.o: $(BASISDIR)/ref_wfn.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/new_basis.o: $(BASISDIR)/new_basis.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/new_basis1.o: $(BASISDIR)/new_basis1.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/new_basis2.o: $(BASISDIR)/new_basis2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/new_basis3.o: $(BASISDIR)/new_basis3.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/find_basis1.o: $(BASISDIR)/find_basis1.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/find_basis2.o: $(BASISDIR)/find_basis2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/out_overlap.o: $(BASISDIR)/out_overlap.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_wfns.o: $(BASISDIR)/save_wfns.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/read_wfns.o: $(BASISDIR)/read_wfns.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/overlaps_gf.o: $(BASISDIR)/overlaps_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/overlaps_susc2.o: $(BASISDIR)/overlaps_susc2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/basis_change.o: $(BASISDIR)/basis_change.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/gradient_susc_basis.o: $(BASISDIR)/gradient_susc_basis.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                        Projection module
########################################################################
$(OBJ)/reg_coeffs.o: $(PROJDIR)/reg_coeffs.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/irr_coeffs.o: $(PROJDIR)/irr_coeffs.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/out_coeffs.o: $(PROJDIR)/out_coeffs.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/reg_coeffs_ia.o: $(PROJDIR)/reg_coeffs_ia.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/irr_coeffs_ia.o: $(PROJDIR)/irr_coeffs_ia.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/out_coeffs_ia.o: $(PROJDIR)/out_coeffs_ia.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/in_coeffs.o: $(PROJDIR)/in_coeffs.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_coeffs.o: $(PROJDIR)/save_coeffs.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_coeffs_ia.o: $(PROJDIR)/save_coeffs_ia.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/out_gmat.o: $(PROJDIR)/out_gmat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/out_tgmat.o: $(PROJDIR)/out_tgmat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/in_tmat.o: $(PROJDIR)/in_tmat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/in_gmat.o: $(PROJDIR)/in_gmat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_tmati.o: $(PROJDIR)/save_tmati.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_gsij.o: $(PROJDIR)/save_gsij.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_tmcoll.o: $(PROJDIR)/save_tmcoll.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_gscoll.o: $(PROJDIR)/save_gscoll.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/projected_gf.o: $(PROJDIR)/projected_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/baryint_gf2.o: $(PROJDIR)/baryint_gf2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ratint_gf2.o: $(PROJDIR)/ratint_gf2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/biratint_gf.o: $(PROJDIR)/biratint_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ratval_gf.o: $(PROJDIR)/ratval_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ratfit_gf.o: $(PROJDIR)/ratfit_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ratfit_gf2.o: $(PROJDIR)/ratfit_gf2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/out_tgmat_soc.o: $(PROJDIR)/out_tgmat_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/reg_coeffs_soc.o: $(PROJDIR)/reg_coeffs_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/irr_coeffs_soc.o: $(PROJDIR)/irr_coeffs_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/out_coeffs_soc.o: $(PROJDIR)/out_coeffs_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/tmat_proj.o: $(PROJDIR)/tmat_proj.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spinrot_scatt_sol.o: $(PROJDIR)/spinrot_scatt_sol.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/in_coeffs_soc.o: $(PROJDIR)/in_coeffs_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/in_gmat_soc.o: $(PROJDIR)/in_gmat_soc.F90    
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/in_tmat_soc.o: $(PROJDIR)/in_tmat_soc.F90 
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                       Groundstate module
########################################################################
$(OBJ)/save_rho2ns.o: $(GSDIR)/save_rho2ns.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/read_rho2ns.o: $(GSDIR)/read_rho2ns.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/update_rho2ns.o: $(GSDIR)/update_rho2ns.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/update_espv.o: $(GSDIR)/update_espv.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/save_rhomat.o: $(GSDIR)/save_rhomat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/read_rhomat.o: $(GSDIR)/read_rhomat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/mix_rhomat.o: $(GSDIR)/mix_rhomat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/groundstate_new.o: $(GSDIR)/groundstate_new.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/get_gfsum.o: $(GSDIR)/get_gfsum.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/integrated_charge.o: $(GSDIR)/integrated_charge.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/charge_density.o: $(GSDIR)/charge_density.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_density.o: $(GSDIR)/current_density.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_soc.o: $(GSDIR)/current_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_zeeman.o: $(GSDIR)/current_zeeman.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_interpolation.o: $(GSDIR)/current_interpolation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/density_of_states.o: $(GSDIR)/density_of_states.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/get_xc2.o: $(GSDIR)/get_xc2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/vxc_vwn.o: $(GSDIR)/vxc_vwn.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/density_matrix.o: $(GSDIR)/density_matrix.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/rho_interpolation.o: $(GSDIR)/rho_interpolation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                       Corrections module
########################################################################
$(OBJ)/build_selfe.o: $(CORRDIR)/build_selfe.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_zeeman.o: $(CORRDIR)/build_zeeman.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vscfb.o: $(CORRDIR)/build_vscfb.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ms_from_bxc3.o: $(CORRDIR)/ms_from_bxc3.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/ms_from_bxc4.o: $(CORRDIR)/ms_from_bxc4.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vxcdiff.o: $(CORRDIR)/build_vxcdiff.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vsocz.o: $(CORRDIR)/build_vsocz.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vsocpm.o: $(CORRDIR)/build_vsocpm.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vldaub.o: $(CORRDIR)/build_vldaub.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vsocb.o: $(CORRDIR)/build_vsocb.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vsoc.o: $(CORRDIR)/build_vsoc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vscfb2.o: $(CORRDIR)/build_vscfb2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/new_directions.o: $(CORRDIR)/new_directions.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spin_directions.o: $(CORRDIR)/spin_directions.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spin_rotation.o: $(CORRDIR)/spin_rotation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/local_frame.o: $(CORRDIR)/local_frame.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/orb_rotation.o: $(CORRDIR)/orb_rotation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vrotb.o: $(CORRDIR)/build_vrotb.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spinrot_tgf.o: $(CORRDIR)/spinrot_tgf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spinrot_gf.o: $(CORRDIR)/spinrot_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spinrot_gf_sph.o: $(CORRDIR)/spinrot_gf_sph.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spinrot_rhomat.o: $(CORRDIR)/spinrot_rhomat.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/lippmann_schwinger.o: $(CORRDIR)/lippmann_schwinger.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/structural_gf.o: $(CORRDIR)/structural_gf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/torque_soc.o: $(CORRDIR)/torque_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/torque_soc2.o: $(CORRDIR)/torque_soc2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/lichtenstein_jij.o: $(CORRDIR)/lichtenstein_jij.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_dti.o: $(CORRDIR)/build_dti.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_dti_v2.o: $(CORRDIR)/build_dti_v2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/anisotropic_jij.o: $(CORRDIR)/anisotropic_jij.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/symmetrize_tgf.o: $(CORRDIR)/symmetrize_tgf.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/pot_correction.o: $(CORRDIR)/pot_correction.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/spinrot_tgf_sph.o: $(CORRDIR)/spinrot_tgf_sph.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_bsoctwist.o: $(CORRDIR)/build_bsoctwist.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vsocz2.o: $(CORRDIR)/build_vsocz2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vsocpm2.o: $(CORRDIR)/build_vsocpm2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_vxcdiff2.o: $(CORRDIR)/build_vxcdiff2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_grad_mass.o: $(CORRDIR)/build_grad_mass.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                     Susceptibility module
########################################################################
$(OBJ)/dyn_susc_expansion.o: $(SUSCDIR)/dyn_susc_expansion.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/dynamic_susc.o: $(SUSCDIR)/dynamic_susc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/dynamic_susc2.o: $(SUSCDIR)/dynamic_susc2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/static_susc.o: $(SUSCDIR)/static_susc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/static_susc2.o: $(SUSCDIR)/static_susc2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/static_susc_more.o: $(SUSCDIR)/static_susc_more.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/static_susc_enhanced.o: $(SUSCDIR)/static_susc_enhanced.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/add_susc_block.o: $(SUSCDIR)/add_susc_block.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/susc_denominator.o: $(SUSCDIR)/susc_denominator.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_kxcalda2.o: $(SUSCDIR)/build_kxcalda2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/kxc_sumrule.o: $(SUSCDIR)/kxc_sumrule.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_khartree.o: $(SUSCDIR)/build_khartree.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/analytic_continuation.o: $(SUSCDIR)/analytic_continuation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/add_susc_block_one_sided.o: $(SUSCDIR)/add_susc_block_one_sided.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_correlation.o: $(SUSCDIR)/current_correlation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_correlation_interpolation.o: $(SUSCDIR)/current_correlation_interpolation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_correlation_divergence.o: $(SUSCDIR)/current_correlation_divergence.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/current_correlation_divergence_interpolation.o: $(SUSCDIR)/current_correlation_divergence_interpolation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/susc_interpolation.o: $(SUSCDIR)/susc_interpolation.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/build_twist.o: $(SUSCDIR)/build_twist.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/kxc_sumrule2.o: $(SUSCDIR)/kxc_sumrule2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                            I/O module
########################################################################
$(OBJ)/load_input.o: $(IODIR)/load_input.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/find_keyinline.o: $(IODIR)/find_keyinline.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/find_keyinfile.o: $(IODIR)/find_keyinfile.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/new_input.o: $(IODIR)/new_input.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_basis_options.o: $(IODIR)/set_basis_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_spinrot_options.o: $(IODIR)/set_spinrot_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_bfield_options.o: $(IODIR)/set_bfield_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_selfe_options.o: $(IODIR)/set_selfe_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_dos_options.o: $(IODIR)/set_dos_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_jij_options.o: $(IODIR)/set_jij_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_soc_options.o: $(IODIR)/set_soc_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_susc_options.o: $(IODIR)/set_susc_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_fit_options.o: $(IODIR)/set_fit_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_ldau_options.o: $(IODIR)/set_ldau_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_groundstate_options.o: $(IODIR)/set_groundstate_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_groups.o: $(IODIR)/set_groups.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/set_global_options.o: $(IODIR)/set_global_options.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/close_outsusc.o: $(IODIR)/close_outsusc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/outsusc_out.o: $(IODIR)/outsusc_out.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/outsusc_in.o: $(IODIR)/outsusc_in.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/outsusc_in2.o: $(IODIR)/outsusc_in2.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

$(OBJ)/outsusc_in_soc.o: $(IODIR)/outsusc_in_soc.F90
	$(FC) $(FFLAGS) -c -I$(OBJ) $< -o $@

########################################################################
#                          Clean up flags
########################################################################
.PHONY: clean

clean:
	@rm -f $(OBJ)/*.mod $(OBJ)/*.o *genmod*
