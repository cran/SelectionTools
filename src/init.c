/* Native routine registration for SelectionTools. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <stdlib.h>

/* .C entry points */
extern void allele_freq(void *, void *, void *, void *, void *, void *, void *, void *);
extern void append_population(void *, void *);
extern void concat_population(void *, void *);
extern void copy_population(void *, void *, void *, void *);
extern void correct_missing(void *, void *, void *, void *, void *, void *, void *, void *);
extern void cross(void *, void *, void *, void *, void *, void *, void *, void *);
extern void dh(void *, void *, void *, void *);
extern void divide_population(void *, void *, void *);
extern void effect_weight_set_all(void *);
extern void effmap_remove(void *);
extern void effmap_remove_all(void);
extern void evaluate_genome(void *, void *, void *, void *, void *, void *);
extern void evaluate_genome_region(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void evaluate_ld(void *, void *, void *, void *, void *);
extern void evaluate_locus(void *, void *, void *, void *, void *);
extern void evaluate_mdp(void *, void *, void *, void *);
extern void evaluate_population(void *, void *);
extern void st_r_info(void *, void *);
extern void gen_dist(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void generate_effect_file(void *, void *);
extern void generate_map_file(void *, void *, void *);
extern void generate_population(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void genotype_population(void *);
extern void get_allele_number(void *, void *);
extern void get_allow_zeros(void *);
extern void get_genome_par(void *, void *);
extern void get_genome_par_b(void *);
extern void get_info_level(void *);
extern void get_map_a(void *, void *, void *, void *);
extern void get_map_b(void *, void *, void *, void *);
extern void get_mdp(void *);
extern void get_population(void *, void *, void *, void *, void *, void *);
extern void get_population_info(void *, void *, void *);
extern void get_population_size(void *, void *, void *);
extern void get_population_value(void *, void *, void *, void *, void *);
extern void gs_build_V_01_GV(void *, void *, void *, void *);
extern void gs_build_Z_01_GV(void *, void *, void *, void *);
extern void gs_build_esvs_01_GV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_calc_ld_01_GV(void *, void *, void *, void *, void *);
extern void gs_calc_q_01_GV(void *, void *, void *, void *, void *, void *);
extern void gs_calc_rf_01_GV(void *, void *, void *, void *, void *);
extern void gs_check_pvals_01_GV(void *, void *);
extern void gs_chrom_stats_01_GV(void *, void *, void *, void *);
extern void gs_copy_marker_data_01_GV(void *, void *, void *);
extern void gs_cross_eval_es_01_GV(void *, void *, void *);
extern void gs_cross_eval_gd_01_GV(void *, void *, void *);
extern void gs_cross_eval_gd_fct_01_GV(void *, void *, void *, void *);
extern void gs_cross_eval_ma_01_GV(void *, void *);
extern void gs_cross_eval_mi_01_GV(void *, void *);
extern void gs_cross_eval_mu_01_GV(void *, void *);
extern void gs_cross_eval_va_01_GV(void *, void *, void *, void *, void *);
extern void gs_cross_info_01_GV(void *, void *, void *, void *, void *);
extern void gs_cross_validation_02_GV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_def_hblocks_01_GV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_esteff_lsq_01_GV(void *, void *, void *, void *);
extern void gs_esteff_rmla_01_GV(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_esteff_rmlc_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_esteff_rmlr_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_esteff_rmlv_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_esteff_rr_01_GV(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_estimate_bv_01_GV(void *, void *, void *, void *, void *);
extern void gs_get_info_level_GV(void *, void *);
extern void gs_get_map_01_GV(void *, void *, void *);
extern void gs_get_num_threads(void *);
extern void gs_lambda_aov_01_GV(void *, void *, void *, void *, void *, void *);
extern void gs_lambda_const_01_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_const_02_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_emstep_01_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_reg_01_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_rmla_01_GV(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_lambda_rmla_02_GV(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_lambda_rmlc_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_lambda_rmlr_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_lambda_rrblup_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_lambda_var_01_GV(void *, void *, void *, void *, void *, void *);
extern void gs_marker_stats_01_GV(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_mme_coeff_01_GV(void *, void *, void *, void *);
extern void gs_mme_restcoeff_01_GV(void *, void *, void *);
extern void gs_mme_restrhs_01_GV(void *, void *, void *);
extern void gs_mme_rhs_01_GV(void *, void *, void *, void *);
extern void gs_mme_solve_01_GV(void *, void *, void *, void *);
extern void gs_mmet_coeff_01_GV(void *, void *, void *, void *);
extern void gs_mmet_rhs_01_GV(void *, void *, void *, void *);
extern void gs_mmet_solve_01_GV(void *, void *, void *, void *);
extern void gs_neg_effall_01_GV(void *, void *);
extern void gs_pos_effall_01_GV(void *, void *);
extern void gs_r_info(void *, void *);
extern void gs_read_map_01_GV(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_read_marker_data_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_read_marker_data_02_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_read_marker_data_03_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_read_marker_data_04_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_read_performance_data_01_GV(void *, void *, void *, void *, void *);
extern void gs_recode_hbc_01_GV(void *, void *, void *);
extern void gs_recode_hil_01_GV(void *, void *, void *, void *);
extern void gs_recode_ref_01_GV(void *, void *, void *, void *);
extern void gs_recode_ref_02_GV(void *, void *, void *, void *, void *);
extern void gs_recode_ref_name_01_GV(void *, void *, void *, void *);
extern void gs_restrict_marker_data_03_GV(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_return_effects_01_GV(void *, void *, void *);
extern void gs_return_performance_data_01_GV(void *, void *, void *, void *);
extern void gs_return_pvals_01_GV(void *, void *, void *);
extern void gs_set_all_info_levels(void *);
extern void gs_set_allele_codes_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_set_effects_01_GV(void *, void *, void *);
extern void gs_set_hblocks_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_set_info_level_GV(void *, void *);
extern void gs_set_num_threads(void *);
extern void gs_set_performance_data_01_GV(void *, void *, void *, void *);
extern void gs_single_marker_aov_01_GV(void *, void *, void *, void *, void *);
extern void gs_single_marker_reg_01_GV(void *, void *, void *, void *);
extern void gs_testeff_lsq_01_GV(void *, void *, void *, void *, void *);
extern void gs_testeff_rmlc_01_GV(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_testeff_rmlv_01_GV(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gs_vc_rrblup_01_GV(void *, void *, void *, void *, void *, void *, void *);
extern void gs_w_aov_01_GV(void *, void *, void *, void *, void *);
extern void gs_write_map_01_GV(void *, void *, void *);
extern void gs_write_marker_data_02_GV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void haplotype_calc_d(void *, void *, void *);
extern void haplotype_free(void);
extern void haplotype_init(void *, void *, void *, void *);
extern void haplotype_next(void *, void *, void *, void *, void *, void *, void *);
extern void il_eval_library_GV(void *, void *, void *, void *, void *, void *);
extern void il_eval_lines_GV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void il_get_population_size_GV(void *, void *);
extern void init_population(void *, void *, void *, void *, void *, void *, void *);
extern void linkage_map_load(void *, void *, void *, void *, void *);
extern void list_eff(void *, void *);
extern void list_eff_parameter(void *, void *);
extern void list_populations(void *, void *);
extern void list_populations_parameter(void *, void *);
extern void load_effmap(void *, void *);
extern void load_effmap_mod(void *, void *, void *);
extern void optimize_population(void *);
extern void phenotype_population(void *, void *, void *, void *);
extern void print_version_2(void *);
extern void r_info(void *, void *);
extern void remove_evaluate_population(void *);
extern void remove_genotype_population(void *);
extern void remove_map(void);
extern void remove_population(void *);
extern void rename_population(void *, void *);
extern void reset_mdp(void);
extern void resize_population(void *, void *);
extern void rng_choose(void *);
extern void rng_info(void);
extern void rng_init(void);
extern void rng_list(void);
extern void sample_population(void *, void *, void *, void *);
extern void select_all_n_best(void *, void *, void *, void *, void *, void *, void *, void *);
extern void set_allow_zeros(void *);
extern void set_co_freq(void *);
extern void set_eff_weight(void *, void *);
extern void set_genome_par(void *, void *, void *);
extern void set_info_level(void *);
extern void set_mdp(void *);
extern void set_nolociinit(void *);
extern void set_population_gvalue(void *, void *, void *, void *);
extern void set_population_info(void *, void *, void *);
extern void sort_population(void *, void *, void *, void *, void *);
extern void ssd_mating(void *, void *, void *, void *, void *);
extern void st_get_info_level(void *);
extern void st_get_simpop2_GV(void *, void *, void *, void *);
extern void st_get_simpop_GV(void *, void *, void *, void *);
extern void st_set_info_level(void *);
extern void st_set_sim_ef(void *, void *, void *);
extern void st_set_sim_gp(void *, void *);
extern void st_set_sim_mp(void *, void *);
extern void st_set_sim_pp(void *, void *, void *);
extern void st_start_timer(void *, void *);
extern void st_stop_timer(void *, void *);
extern void swap_population_name(void *, void *);

/* .Call entry points */
extern SEXP allele_evaluate(SEXP, SEXP, SEXP, SEXP);
extern SEXP define_effects_df(SEXP, SEXP);
extern SEXP define_map_df(SEXP);
extern SEXP evaluate_haplotype(SEXP, SEXP, SEXP, SEXP);
extern SEXP genotype_evaluate(SEXP, SEXP, SEXP);
extern SEXP gs_read_marker_data_01_DF(SEXP, SEXP);
extern SEXP haplotype_return_genotype(void);
extern SEXP haplotype_return_table(void);
extern SEXP return_population(SEXP, SEXP);
extern SEXP splitdot(SEXP);
extern SEXP splitdt(SEXP);
extern SEXP st_destroy_phase(SEXP);
extern SEXP st_get_V(SEXP);
extern SEXP st_get_Z(SEXP);
extern SEXP st_get_im(SEXP, SEXP, SEXP);
extern SEXP st_get_y(SEXP);
extern SEXP st_has_openmp(void);
extern SEXP st_load_performance_data(SEXP, SEXP, SEXP);
extern SEXP st_load_vcf_data(SEXP, SEXP);
extern SEXP st_markerdata_to_STvcf(SEXP);
extern SEXP st_mixed(SEXP);
extern SEXP st_mxd_interface(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP st_phase_03_step1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP st_read_vcf(SEXP, SEXP);
extern SEXP st_switch_error(SEXP, SEXP);
extern SEXP st_write_vcf(SEXP, SEXP);

/*
 * These R-visible names are present in the R API but do not have a compatible
 * native implementation in the current sources.  Keeping explicit adapters
 * makes registration complete and turns the pre-existing invalid calls into
 * deterministic R errors instead of unsafe or unrestricted symbol lookup.
 */
static void st_unavailable_gs_cross_info_gd_01_GV(void *a1, void *a2, void *a3)
{
    (void) a1;
    (void) a2;
    (void) a3;
    error("Native routine 'gs_cross_info_gd_01_GV' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_gs_lambda_rmlv_01_GV(void *a1, void *a2, void *a3, void *a4, void *a5, void *a6)
{
    (void) a1;
    (void) a2;
    (void) a3;
    (void) a4;
    (void) a5;
    (void) a6;
    error("Native routine 'gs_lambda_rmlv_01_GV' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_gs_mme_invcoeff_01_GV(void *a1, void *a2, void *a3)
{
    (void) a1;
    (void) a2;
    (void) a3;
    error("Native routine 'gs_mme_invcoeff_01_GV' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_gs_mmet_coeff_03_GV(void *a1, void *a2, void *a3, void *a4)
{
    (void) a1;
    (void) a2;
    (void) a3;
    (void) a4;
    error("Native routine 'gs_mmet_coeff_03_GV' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_gs_reset(void)
{
    error("Native routine 'gs_reset' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_gs_restrict_marker_data_01_GV(void *a1, void *a2, void *a3, void *a4, void *a5, void *a6, void *a7, void *a8)
{
    (void) a1;
    (void) a2;
    (void) a3;
    (void) a4;
    (void) a5;
    (void) a6;
    (void) a7;
    (void) a8;
    error("Native routine 'gs_restrict_marker_data_01_GV' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_gs_set_openblas_threads(void *a1)
{
    (void) a1;
    error("Native routine 'gs_set_openblas_threads' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_gs_test_mp(void *a1)
{
    (void) a1;
    error("Native routine 'gs_test_mp' is not implemented by this SelectionTools source tree");
}

static void st_unavailable_st_set_matrix_ops(void *a1, void *a2)
{
    (void) a1;
    (void) a2;
    error("Native routine 'st_set_matrix_ops' is not implemented by this SelectionTools source tree");
}

static const R_CMethodDef CEntries[] = {
    {"allele_freq", (DL_FUNC) &allele_freq, 8},
    {"append_population", (DL_FUNC) &append_population, 2},
    {"concat_population", (DL_FUNC) &concat_population, 2},
    {"copy_population", (DL_FUNC) &copy_population, 4},
    {"correct_missing", (DL_FUNC) &correct_missing, 8},
    {"cross", (DL_FUNC) &cross, 8},
    {"dh", (DL_FUNC) &dh, 4},
    {"divide_population", (DL_FUNC) &divide_population, 3},
    {"effect_weight_set_all", (DL_FUNC) &effect_weight_set_all, 1},
    {"effmap_remove", (DL_FUNC) &effmap_remove, 1},
    {"effmap_remove_all", (DL_FUNC) &effmap_remove_all, 0},
    {"evaluate_genome", (DL_FUNC) &evaluate_genome, 6},
    {"evaluate_genome_region", (DL_FUNC) &evaluate_genome_region, 9},
    {"evaluate_ld", (DL_FUNC) &evaluate_ld, 5},
    {"evaluate_locus", (DL_FUNC) &evaluate_locus, 5},
    {"evaluate_mdp", (DL_FUNC) &evaluate_mdp, 4},
    {"evaluate_population", (DL_FUNC) &evaluate_population, 2},
    {"g_r_info", (DL_FUNC) &st_r_info, 2},
    {"gen_dist", (DL_FUNC) &gen_dist, 17},
    {"generate_effect_file", (DL_FUNC) &generate_effect_file, 2},
    {"generate_map_file", (DL_FUNC) &generate_map_file, 3},
    {"generate_population", (DL_FUNC) &generate_population, 17},
    {"genotype_population", (DL_FUNC) &genotype_population, 1},
    {"get_allele_number", (DL_FUNC) &get_allele_number, 2},
    {"get_allow_zeros", (DL_FUNC) &get_allow_zeros, 1},
    {"get_genome_par", (DL_FUNC) &get_genome_par, 2},
    {"get_genome_par_b", (DL_FUNC) &get_genome_par_b, 1},
    {"get_info_level", (DL_FUNC) &get_info_level, 1},
    {"get_map_a", (DL_FUNC) &get_map_a, 4},
    {"get_map_b", (DL_FUNC) &get_map_b, 4},
    {"get_mdp", (DL_FUNC) &get_mdp, 1},
    {"get_population", (DL_FUNC) &get_population, 6},
    {"get_population_info", (DL_FUNC) &get_population_info, 3},
    {"get_population_size", (DL_FUNC) &get_population_size, 3},
    {"get_population_value", (DL_FUNC) &get_population_value, 5},
    {"gs_build_V_01_GV", (DL_FUNC) &gs_build_V_01_GV, 4},
    {"gs_build_Z_01_GV", (DL_FUNC) &gs_build_Z_01_GV, 4},
    {"gs_build_esvs_01_GV", (DL_FUNC) &gs_build_esvs_01_GV, 10},
    {"gs_calc_ld_01_GV", (DL_FUNC) &gs_calc_ld_01_GV, 5},
    {"gs_calc_q_01_GV", (DL_FUNC) &gs_calc_q_01_GV, 6},
    {"gs_calc_rf_01_GV", (DL_FUNC) &gs_calc_rf_01_GV, 5},
    {"gs_check_pvals_01_GV", (DL_FUNC) &gs_check_pvals_01_GV, 2},
    {"gs_chrom_stats_01_GV", (DL_FUNC) &gs_chrom_stats_01_GV, 4},
    {"gs_copy_marker_data_01_GV", (DL_FUNC) &gs_copy_marker_data_01_GV, 3},
    {"gs_cross_eval_es_01_GV", (DL_FUNC) &gs_cross_eval_es_01_GV, 3},
    {"gs_cross_eval_gd_01_GV", (DL_FUNC) &gs_cross_eval_gd_01_GV, 3},
    {"gs_cross_eval_gd_fct_01_GV", (DL_FUNC) &gs_cross_eval_gd_fct_01_GV, 4},
    {"gs_cross_eval_ma_01_GV", (DL_FUNC) &gs_cross_eval_ma_01_GV, 2},
    {"gs_cross_eval_mi_01_GV", (DL_FUNC) &gs_cross_eval_mi_01_GV, 2},
    {"gs_cross_eval_mu_01_GV", (DL_FUNC) &gs_cross_eval_mu_01_GV, 2},
    {"gs_cross_eval_va_01_GV", (DL_FUNC) &gs_cross_eval_va_01_GV, 5},
    {"gs_cross_info_01_GV", (DL_FUNC) &gs_cross_info_01_GV, 5},
    {"gs_cross_info_gd_01_GV", (DL_FUNC) &st_unavailable_gs_cross_info_gd_01_GV, 3},
    {"gs_cross_validation_02_GV", (DL_FUNC) &gs_cross_validation_02_GV, 12},
    {"gs_def_hblocks_01_GV", (DL_FUNC) &gs_def_hblocks_01_GV, 10},
    {"gs_esteff_lsq_01_GV", (DL_FUNC) &gs_esteff_lsq_01_GV, 4},
    {"gs_esteff_rmla_01_GV", (DL_FUNC) &gs_esteff_rmla_01_GV, 8},
    {"gs_esteff_rmlc_01_GV", (DL_FUNC) &gs_esteff_rmlc_01_GV, 7},
    {"gs_esteff_rmlr_01_GV", (DL_FUNC) &gs_esteff_rmlr_01_GV, 7},
    {"gs_esteff_rmlv_01_GV", (DL_FUNC) &gs_esteff_rmlv_01_GV, 7},
    {"gs_esteff_rr_01_GV", (DL_FUNC) &gs_esteff_rr_01_GV, 9},
    {"gs_estimate_bv_01_GV", (DL_FUNC) &gs_estimate_bv_01_GV, 5},
    {"gs_get_info_level_GV", (DL_FUNC) &gs_get_info_level_GV, 2},
    {"gs_get_map_01_GV", (DL_FUNC) &gs_get_map_01_GV, 3},
    {"gs_get_num_threads", (DL_FUNC) &gs_get_num_threads, 1},
    {"gs_lambda_aov_01_GV", (DL_FUNC) &gs_lambda_aov_01_GV, 6},
    {"gs_lambda_const_01_GV", (DL_FUNC) &gs_lambda_const_01_GV, 5},
    {"gs_lambda_const_02_GV", (DL_FUNC) &gs_lambda_const_02_GV, 5},
    {"gs_lambda_emstep_01_GV", (DL_FUNC) &gs_lambda_emstep_01_GV, 5},
    {"gs_lambda_reg_01_GV", (DL_FUNC) &gs_lambda_reg_01_GV, 5},
    {"gs_lambda_rmla_01_GV", (DL_FUNC) &gs_lambda_rmla_01_GV, 8},
    {"gs_lambda_rmla_02_GV", (DL_FUNC) &gs_lambda_rmla_02_GV, 8},
    {"gs_lambda_rmlc_01_GV", (DL_FUNC) &gs_lambda_rmlc_01_GV, 7},
    {"gs_lambda_rmlr_01_GV", (DL_FUNC) &gs_lambda_rmlr_01_GV, 7},
    {"gs_lambda_rmlv_01_GV", (DL_FUNC) &st_unavailable_gs_lambda_rmlv_01_GV, 6},
    {"gs_lambda_rrblup_01_GV", (DL_FUNC) &gs_lambda_rrblup_01_GV, 7},
    {"gs_lambda_var_01_GV", (DL_FUNC) &gs_lambda_var_01_GV, 6},
    {"gs_marker_stats_01_GV", (DL_FUNC) &gs_marker_stats_01_GV, 9},
    {"gs_mme_coeff_01_GV", (DL_FUNC) &gs_mme_coeff_01_GV, 4},
    {"gs_mme_invcoeff_01_GV", (DL_FUNC) &st_unavailable_gs_mme_invcoeff_01_GV, 3},
    {"gs_mme_restcoeff_01_GV", (DL_FUNC) &gs_mme_restcoeff_01_GV, 3},
    {"gs_mme_restrhs_01_GV", (DL_FUNC) &gs_mme_restrhs_01_GV, 3},
    {"gs_mme_rhs_01_GV", (DL_FUNC) &gs_mme_rhs_01_GV, 4},
    {"gs_mme_solve_01_GV", (DL_FUNC) &gs_mme_solve_01_GV, 4},
    {"gs_mmet_coeff_01_GV", (DL_FUNC) &gs_mmet_coeff_01_GV, 4},
    {"gs_mmet_coeff_03_GV", (DL_FUNC) &st_unavailable_gs_mmet_coeff_03_GV, 4},
    {"gs_mmet_rhs_01_GV", (DL_FUNC) &gs_mmet_rhs_01_GV, 4},
    {"gs_mmet_solve_01_GV", (DL_FUNC) &gs_mmet_solve_01_GV, 4},
    {"gs_neg_effall_01_GV", (DL_FUNC) &gs_neg_effall_01_GV, 2},
    {"gs_pos_effall_01_GV", (DL_FUNC) &gs_pos_effall_01_GV, 2},
    {"gs_r_info", (DL_FUNC) &gs_r_info, 2},
    {"gs_read_map_01_GV", (DL_FUNC) &gs_read_map_01_GV, 8},
    {"gs_read_marker_data_01_GV", (DL_FUNC) &gs_read_marker_data_01_GV, 7},
    {"gs_read_marker_data_02_GV", (DL_FUNC) &gs_read_marker_data_02_GV, 7},
    {"gs_read_marker_data_03_GV", (DL_FUNC) &gs_read_marker_data_03_GV, 7},
    {"gs_read_marker_data_04_GV", (DL_FUNC) &gs_read_marker_data_04_GV, 7},
    {"gs_read_performance_data_01_GV", (DL_FUNC) &gs_read_performance_data_01_GV, 5},
    {"gs_recode_hbc_01_GV", (DL_FUNC) &gs_recode_hbc_01_GV, 3},
    {"gs_recode_hil_01_GV", (DL_FUNC) &gs_recode_hil_01_GV, 4},
    {"gs_recode_ref_01_GV", (DL_FUNC) &gs_recode_ref_01_GV, 4},
    {"gs_recode_ref_02_GV", (DL_FUNC) &gs_recode_ref_02_GV, 5},
    {"gs_recode_ref_name_01_GV", (DL_FUNC) &gs_recode_ref_name_01_GV, 4},
    {"gs_reset", (DL_FUNC) &st_unavailable_gs_reset, 0},
    {"gs_restrict_marker_data_01_GV", (DL_FUNC) &st_unavailable_gs_restrict_marker_data_01_GV, 8},
    {"gs_restrict_marker_data_03_GV", (DL_FUNC) &gs_restrict_marker_data_03_GV, 8},
    {"gs_return_effects_01_GV", (DL_FUNC) &gs_return_effects_01_GV, 3},
    {"gs_return_performance_data_01_GV", (DL_FUNC) &gs_return_performance_data_01_GV, 4},
    {"gs_return_pvals_01_GV", (DL_FUNC) &gs_return_pvals_01_GV, 3},
    {"gs_set_all_info_levels", (DL_FUNC) &gs_set_all_info_levels, 1},
    {"gs_set_allele_codes_01_GV", (DL_FUNC) &gs_set_allele_codes_01_GV, 7},
    {"gs_set_effects_01_GV", (DL_FUNC) &gs_set_effects_01_GV, 3},
    {"gs_set_hblocks_01_GV", (DL_FUNC) &gs_set_hblocks_01_GV, 7},
    {"gs_set_info_level_GV", (DL_FUNC) &gs_set_info_level_GV, 2},
    {"gs_set_num_threads", (DL_FUNC) &gs_set_num_threads, 1},
    {"gs_set_openblas_threads", (DL_FUNC) &st_unavailable_gs_set_openblas_threads, 1},
    {"gs_set_performance_data_01_GV", (DL_FUNC) &gs_set_performance_data_01_GV, 4},
    {"gs_single_marker_aov_01_GV", (DL_FUNC) &gs_single_marker_aov_01_GV, 5},
    {"gs_single_marker_reg_01_GV", (DL_FUNC) &gs_single_marker_reg_01_GV, 4},
    {"gs_test_mp", (DL_FUNC) &st_unavailable_gs_test_mp, 1},
    {"gs_testeff_lsq_01_GV", (DL_FUNC) &gs_testeff_lsq_01_GV, 5},
    {"gs_testeff_rmlc_01_GV", (DL_FUNC) &gs_testeff_rmlc_01_GV, 8},
    {"gs_testeff_rmlv_01_GV", (DL_FUNC) &gs_testeff_rmlv_01_GV, 8},
    {"gs_vc_rrblup_01_GV", (DL_FUNC) &gs_vc_rrblup_01_GV, 7},
    {"gs_w_aov_01_GV", (DL_FUNC) &gs_w_aov_01_GV, 5},
    {"gs_write_map_01_GV", (DL_FUNC) &gs_write_map_01_GV, 3},
    {"gs_write_marker_data_02_GV", (DL_FUNC) &gs_write_marker_data_02_GV, 14},
    {"haplotype_calc_d", (DL_FUNC) &haplotype_calc_d, 3},
    {"haplotype_free", (DL_FUNC) &haplotype_free, 0},
    {"haplotype_init", (DL_FUNC) &haplotype_init, 4},
    {"haplotype_next", (DL_FUNC) &haplotype_next, 7},
    {"il_eval_library_GV", (DL_FUNC) &il_eval_library_GV, 6},
    {"il_eval_lines_GV", (DL_FUNC) &il_eval_lines_GV, 11},
    {"il_get_population_size_GV", (DL_FUNC) &il_get_population_size_GV, 2},
    {"init_population", (DL_FUNC) &init_population, 7},
    {"linkage_map_load", (DL_FUNC) &linkage_map_load, 5},
    {"list_eff", (DL_FUNC) &list_eff, 2},
    {"list_eff_parameter", (DL_FUNC) &list_eff_parameter, 2},
    {"list_populations", (DL_FUNC) &list_populations, 2},
    {"list_populations_parameter", (DL_FUNC) &list_populations_parameter, 2},
    {"load_effmap", (DL_FUNC) &load_effmap, 2},
    {"load_effmap_mod", (DL_FUNC) &load_effmap_mod, 3},
    {"optimize_population", (DL_FUNC) &optimize_population, 1},
    {"phenotype_population", (DL_FUNC) &phenotype_population, 4},
    {"print_version_2", (DL_FUNC) &print_version_2, 1},
    {"r_info", (DL_FUNC) &r_info, 2},
    {"remove_evaluate_population", (DL_FUNC) &remove_evaluate_population, 1},
    {"remove_genotype_population", (DL_FUNC) &remove_genotype_population, 1},
    {"remove_map", (DL_FUNC) &remove_map, 0},
    {"remove_population", (DL_FUNC) &remove_population, 1},
    {"rename_population", (DL_FUNC) &rename_population, 2},
    {"reset_mdp", (DL_FUNC) &reset_mdp, 0},
    {"resize_population", (DL_FUNC) &resize_population, 2},
    {"rng_choose", (DL_FUNC) &rng_choose, 1},
    {"rng_info", (DL_FUNC) &rng_info, 0},
    {"rng_init", (DL_FUNC) &rng_init, 0},
    {"rng_list", (DL_FUNC) &rng_list, 0},
    {"sample_population", (DL_FUNC) &sample_population, 4},
    {"select_all_n_best", (DL_FUNC) &select_all_n_best, 8},
    {"set_allow_zeros", (DL_FUNC) &set_allow_zeros, 1},
    {"set_co_freq", (DL_FUNC) &set_co_freq, 1},
    {"set_eff_weight", (DL_FUNC) &set_eff_weight, 2},
    {"set_genome_par", (DL_FUNC) &set_genome_par, 3},
    {"set_info_level", (DL_FUNC) &set_info_level, 1},
    {"set_mdp", (DL_FUNC) &set_mdp, 1},
    {"set_nolociinit", (DL_FUNC) &set_nolociinit, 1},
    {"set_population_gvalue", (DL_FUNC) &set_population_gvalue, 4},
    {"set_population_info", (DL_FUNC) &set_population_info, 3},
    {"sort_population", (DL_FUNC) &sort_population, 5},
    {"ssd_mating", (DL_FUNC) &ssd_mating, 5},
    {"st_get_info_level", (DL_FUNC) &st_get_info_level, 1},
    {"st_get_simpop2_GV", (DL_FUNC) &st_get_simpop2_GV, 4},
    {"st_get_simpop_GV", (DL_FUNC) &st_get_simpop_GV, 4},
    {"st_r_info", (DL_FUNC) &st_r_info, 2},
    {"st_set_info_level", (DL_FUNC) &st_set_info_level, 1},
    {"st_set_matrix_ops", (DL_FUNC) &st_unavailable_st_set_matrix_ops, 2},
    {"st_set_sim_ef", (DL_FUNC) &st_set_sim_ef, 3},
    {"st_set_sim_gp", (DL_FUNC) &st_set_sim_gp, 2},
    {"st_set_sim_mp", (DL_FUNC) &st_set_sim_mp, 2},
    {"st_set_sim_pp", (DL_FUNC) &st_set_sim_pp, 3},
    {"st_start_timer", (DL_FUNC) &st_start_timer, 2},
    {"st_stop_timer", (DL_FUNC) &st_stop_timer, 2},
    {"swap_population_name", (DL_FUNC) &swap_population_name, 2},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"allele_evaluate", (DL_FUNC) &allele_evaluate, 4},
    {"define_effects_df", (DL_FUNC) &define_effects_df, 2},
    {"define_map_df", (DL_FUNC) &define_map_df, 1},
    {"evaluate_haplotype", (DL_FUNC) &evaluate_haplotype, 4},
    {"genotype_evaluate", (DL_FUNC) &genotype_evaluate, 3},
    {"gs_read_marker_data_01_DF", (DL_FUNC) &gs_read_marker_data_01_DF, 2},
    {"haplotype_return_genotype", (DL_FUNC) &haplotype_return_genotype, 0},
    {"haplotype_return_table", (DL_FUNC) &haplotype_return_table, 0},
    {"return_population", (DL_FUNC) &return_population, 2},
    {"splitdot", (DL_FUNC) &splitdot, 1},
    {"splitdt", (DL_FUNC) &splitdt, 1},
    {"st_destroy_phase", (DL_FUNC) &st_destroy_phase, 1},
    {"st_get_V", (DL_FUNC) &st_get_V, 1},
    {"st_get_Z", (DL_FUNC) &st_get_Z, 1},
    {"st_get_im", (DL_FUNC) &st_get_im, 3},
    {"st_get_y", (DL_FUNC) &st_get_y, 1},
    {"st_has_openmp", (DL_FUNC) &st_has_openmp, 0},
    {"st_load_performance_data", (DL_FUNC) &st_load_performance_data, 3},
    {"st_load_vcf_data", (DL_FUNC) &st_load_vcf_data, 2},
    {"st_markerdata_to_STvcf", (DL_FUNC) &st_markerdata_to_STvcf, 1},
    {"st_mixed", (DL_FUNC) &st_mixed, 1},
    {"st_mxd_interface", (DL_FUNC) &st_mxd_interface, 8},
    {"st_phase_03_step1", (DL_FUNC) &st_phase_03_step1, 26},
    {"st_read_vcf", (DL_FUNC) &st_read_vcf, 2},
    {"st_switch_error", (DL_FUNC) &st_switch_error, 2},
    {"st_write_vcf", (DL_FUNC) &st_write_vcf, 2},
    {NULL, NULL, 0}
};

void attribute_visible
R_init_SelectionTools(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
