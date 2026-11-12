
#include <R.h>
#include "SelectionToolsAllele.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <stdlib.h>

extern SEXP st_phase(SEXP, SEXP, SEXP, SEXP, SEXP);

extern void allele_freq(void *, void *, void *, void *, void *, void *, void *,
                        void *);
extern void append_population(void *, void *);
extern void concat_population(void *, void *);
extern void copy_population(void *, void *, void *, void *);
extern void correct_missing(void *, void *, void *, void *, void *, void *,
                            void *, void *);
extern void cross(void *, void *, void *, void *, void *, void *, void *,
                  void *);
extern void dh(void *, void *, void *, void *);
extern void divide_population(void *, void *, void *);
extern void effect_weight_set_all(void *);
extern void effmap_remove(void *);
extern void effmap_remove_all(void);
extern void evaluate_genome(void *, void *, void *, void *, void *, void *);
extern void evaluate_genome_region(void *, void *, void *, void *, void *,
                                   void *, void *, void *, void *);
extern void evaluate_ld(void *, void *, void *, void *, void *);
extern void evaluate_locus(void *, void *, void *, void *, void *);
extern void evaluate_mdp(void *, void *, void *, void *);
extern void evaluate_population(void *, void *);
extern void st_r_info(void *, void *);
extern void gen_dist(void *, void *, void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *, void *, void *,
                     void *, void *, void *);
extern void generate_effect_file(void *, void *);
extern void generate_map_file(void *, void *, void *);
extern void generate_population(void *, void *, void *, void *, void *, void *,
                                void *, void *, void *, void *, void *, void *,
                                void *, void *, void *, void *, void *);
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
extern void gs_build_esvs_01_GV(void *, void *, void *, void *, void *, void *,
                                void *, void *, void *, void *);
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
extern void gs_cross_validation_02_GV(void *, void *, void *, void *, void *,
                                      void *, void *, void *, void *, void *,
                                      void *, void *);
extern void gs_def_hblocks_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *);
extern void gs_esteff_lsq_01_GV(void *, void *, void *, void *);
extern void gs_esteff_rmla_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *, void *);
extern void gs_esteff_rmlc_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *);
extern void gs_esteff_rmlr_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *);
extern void gs_esteff_rmlv_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *);
extern void gs_esteff_rr_01_GV(void *, void *, void *, void *, void *, void *,
                               void *, void *, void *);
extern void gs_estimate_bv_01_GV(void *, void *, void *, void *, void *);
extern void gs_get_info_level_GV(void *, void *);
extern void gs_get_map_01_GV(void *, void *, void *);
extern void gs_get_num_threads(void *);
extern void gs_lambda_aov_01_GV(void *, void *, void *, void *, void *, void *);
extern void gs_lambda_const_01_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_const_02_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_emstep_01_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_reg_01_GV(void *, void *, void *, void *, void *);
extern void gs_lambda_rmla_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *, void *);
extern void gs_lambda_rmla_02_GV(void *, void *, void *, void *, void *, void *,
                                 void *, void *);
extern void gs_lambda_rmlc_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *);
extern void gs_lambda_rmlr_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *);
extern void gs_lambda_rrblup_01_GV(void *, void *, void *, void *, void *,
                                   void *, void *);
extern void gs_lambda_var_01_GV(void *, void *, void *, void *, void *, void *);
extern void gs_marker_stats_01_GV(void *, void *, void *, void *, void *,
                                  void *, void *, void *, void *);
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
extern void gs_read_map_01_GV(void *, void *, void *, void *, void *, void *,
                              void *, void *);
extern void gs_read_marker_data_01_GV(void *, void *, void *, void *, void *,
                                      void *, void *);
extern void gs_read_marker_data_02_GV(void *, void *, void *, void *, void *,
                                      void *, void *);
extern void gs_read_marker_data_03_GV(void *, void *, void *, void *, void *,
                                      void *, void *);
extern void gs_read_marker_data_04_GV(void *, void *, void *, void *, void *,
                                      void *, void *);
extern void gs_read_performance_data_01_GV(void *, void *, void *, void *,
                                           void *);
extern void gs_recode_hbc_01_GV(void *, void *, void *);
extern void gs_recode_hil_01_GV(void *, void *, void *, void *);
extern void gs_recode_ref_01_GV(void *, void *, void *, void *);
extern void gs_recode_ref_02_GV(void *, void *, void *, void *, void *);
extern void gs_recode_ref_name_01_GV(void *, void *, void *, void *);
extern void gs_restrict_marker_data_03_GV(void *, void *, void *, void *,
                                          void *, void *, void *, void *);
extern void gs_return_effects_01_GV(void *, void *, void *);
extern void gs_return_performance_data_01_GV(void *, void *, void *, void *);
extern void gs_return_pvals_01_GV(void *, void *, void *);
extern void gs_set_all_info_levels(void *);
extern void gs_set_allele_codes_01_GV(void *, void *, void *, void *, void *,
                                      void *, void *);
extern void gs_set_effects_01_GV(void *, void *, void *);
extern void gs_set_hblocks_01_GV(void *, void *, void *, void *, void *, void *,
                                 void *);
extern void gs_set_info_level_GV(void *, void *);
extern void gs_set_num_threads(void *);
extern void gs_set_performance_data_01_GV(void *, void *, void *, void *);
extern void gs_single_marker_aov_01_GV(void *, void *, void *, void *, void *);
extern void gs_single_marker_reg_01_GV(void *, void *, void *, void *);
extern void gs_testeff_lsq_01_GV(void *, void *, void *, void *, void *);
extern void gs_testeff_rmlc_01_GV(void *, void *, void *, void *, void *,
                                  void *, void *, void *);
extern void gs_testeff_rmlv_01_GV(void *, void *, void *, void *, void *,
                                  void *, void *, void *);
extern void gs_vc_rrblup_01_GV(void *, void *, void *, void *, void *, void *,
                               void *);
extern void gs_w_aov_01_GV(void *, void *, void *, void *, void *);
extern void gs_write_map_01_GV(void *, void *, void *);
extern void gs_write_marker_data_02_GV(void *, void *, void *, void *, void *,
                                       void *, void *, void *, void *, void *,
                                       void *, void *, void *, void *);
extern void haplotype_calc_d(void *, void *, void *);
extern void haplotype_free(void);
extern void haplotype_init(void *, void *, void *, void *);
extern void haplotype_next(void *, void *, void *, void *, void *, void *,
                           void *);
extern void il_eval_library_GV(void *, void *, void *, void *, void *, void *);
extern void il_eval_lines_GV(void *, void *, void *, void *, void *, void *,
                             void *, void *, void *, void *, void *);
extern void il_get_population_size_GV(void *, void *);
extern void init_population(void *, void *, void *, void *, void *, void *,
                            void *);
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
extern void rng_init_seed(int *, int *, int *);
extern void rng_list(void);
extern void sample_population(void *, void *, void *, void *);
extern void select_all_n_best(void *, void *, void *, void *, void *, void *,
                              void *, void *);
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

extern SEXP st_lmm_frontend(SEXP);
extern SEXP st_lmm_capture_body(SEXP);
extern SEXP st_lmm_capture_warning(SEXP, SEXP);
extern SEXP st_population_pedigree(SEXP, SEXP);
extern SEXP st_marker_results(SEXP);
extern SEXP st_marker_predict(SEXP);
extern SEXP st_marker_phenotypes(SEXP);
extern SEXP st_simulation_identity(SEXP, SEXP, SEXP);
extern SEXP st_population_marker_map(SEXP);
extern SEXP st_extend_marker_map(SEXP);
extern SEXP st_population_cross(SEXP, SEXP, SEXP, SEXP);
extern SEXP allele_evaluate(SEXP, SEXP, SEXP, SEXP);
extern SEXP define_effects_df(SEXP, SEXP);
extern SEXP define_map_df(SEXP);
extern SEXP evaluate_haplotype(SEXP, SEXP, SEXP, SEXP);
extern SEXP genotype_evaluate(SEXP, SEXP, SEXP);
extern SEXP gs_read_marker_data_DF_list(SEXP, SEXP);
extern SEXP gs_read_marker_data_DF_matrix(SEXP, SEXP);
extern SEXP haplotype_return_genotype(void);
extern SEXP haplotype_return_table(void);
extern SEXP return_population(SEXP, SEXP);
extern SEXP splitdot(SEXP);
extern SEXP splitdt(SEXP);
extern SEXP st_destroy_phase(SEXP);
extern SEXP st_read_map_df(SEXP, SEXP, SEXP);
extern SEXP st_get_V(SEXP);
extern SEXP gs_make_A(SEXP, SEXP);
extern SEXP st_get_Z(SEXP);
extern SEXP st_get_Z_named(SEXP, SEXP);
extern SEXP st_get_im(SEXP, SEXP, SEXP);
extern SEXP st_get_y(SEXP);
extern SEXP st_get_y_named(SEXP, SEXP);
extern SEXP st_has_openmp(void);
extern SEXP st_thread_control(SEXP);
extern SEXP st_load_performance_data(SEXP, SEXP, SEXP);
extern SEXP st_load_vcf_data(SEXP, SEXP, SEXP);
extern SEXP st_markerdata_to_STvcf(SEXP);
extern SEXP st_mxd_mv(SEXP);
extern SEXP st_mxd_mv_interface(SEXP);
extern SEXP st_mxd_mv_block(SEXP);
extern SEXP st_lmm_compile(SEXP);
extern SEXP st_kernel_structure(SEXP);
extern SEXP st_get_A(SEXP, SEXP);
extern SEXP st_lmm(SEXP);
extern SEXP st_lsm(SEXP);
extern SEXP st_predict(SEXP);
extern SEXP st_predict_method(SEXP);
extern SEXP st_lsm_method(SEXP);
extern SEXP st_lmm_method(SEXP);
extern SEXP st_mxd_mv_linear(SEXP);
extern SEXP st_read_vcf(SEXP, SEXP, SEXP);
extern SEXP st_switch_error(SEXP, SEXP);
extern SEXP st_write_vcf(SEXP, SEXP);

static void st_unavailable_gs_cross_info_gd_01_GV(void *a1, void *a2, void *a3)
{
    (void)a1;
    (void)a2;
    (void)a3;
    error("Native routine 'gs_cross_info_gd_01_GV' is not implemented by this "
          "SelectionTools source tree");
}

static void st_unavailable_gs_lambda_rmlv_01_GV(void *a1, void *a2, void *a3,
                                                void *a4, void *a5, void *a6)
{
    (void)a1;
    (void)a2;
    (void)a3;
    (void)a4;
    (void)a5;
    (void)a6;
    error("Native routine 'gs_lambda_rmlv_01_GV' is not implemented by this "
          "SelectionTools source tree");
}

static void st_unavailable_gs_mme_invcoeff_01_GV(void *a1, void *a2, void *a3)
{
    (void)a1;
    (void)a2;
    (void)a3;
    error("Native routine 'gs_mme_invcoeff_01_GV' is not implemented by this "
          "SelectionTools source tree");
}

static void st_unavailable_gs_mmet_coeff_03_GV(void *a1, void *a2, void *a3,
                                               void *a4)
{
    (void)a1;
    (void)a2;
    (void)a3;
    (void)a4;
    error("Native routine 'gs_mmet_coeff_03_GV' is not implemented by this "
          "SelectionTools source tree");
}

static void st_unavailable_gs_reset(void)
{
    error("Native routine 'gs_reset' is not implemented by this SelectionTools "
          "source tree");
}

static void st_unavailable_gs_restrict_marker_data_01_GV(void *a1, void *a2,
                                                         void *a3, void *a4,
                                                         void *a5, void *a6,
                                                         void *a7, void *a8)
{
    (void)a1;
    (void)a2;
    (void)a3;
    (void)a4;
    (void)a5;
    (void)a6;
    (void)a7;
    (void)a8;
    error("Native routine 'gs_restrict_marker_data_01_GV' is not implemented "
          "by this SelectionTools source tree");
}

static void st_unavailable_gs_set_openblas_threads(void *a1)
{
    (void)a1;
    error("Native routine 'gs_set_openblas_threads' is not implemented by this "
          "SelectionTools source tree");
}

static void st_unavailable_gs_test_mp(void *a1)
{
    (void)a1;
    error("Native routine 'gs_test_mp' is not implemented by this "
          "SelectionTools source tree");
}

static void st_unavailable_st_set_matrix_ops(void *a1, void *a2)
{
    (void)a1;
    (void)a2;
    error("Native routine 'st_set_matrix_ops' is not implemented by this "
          "SelectionTools source tree");
}

typedef struct
{
    void *value[17];
} st_C_arguments;
static void st_C_body_allele_freq(void *payload)
{
    st_C_arguments *args = payload;
    allele_freq(args->value[0], args->value[1], args->value[2], args->value[3],
                args->value[4], args->value[5], args->value[6], args->value[7]);
}
static void st_C_interface_allele_freq(void *a0, void *a1, void *a2, void *a3,
                                       void *a4, void *a5, void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_allele_freq, &args);
}
static void st_C_body_append_population(void *payload)
{
    st_C_arguments *args = payload;
    append_population(args->value[0], args->value[1]);
}
static void st_C_interface_append_population(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_append_population, &args);
}
static void st_C_body_concat_population(void *payload)
{
    st_C_arguments *args = payload;
    concat_population(args->value[0], args->value[1]);
}
static void st_C_interface_concat_population(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_concat_population, &args);
}
static void st_C_body_copy_population(void *payload)
{
    st_C_arguments *args = payload;
    copy_population(args->value[0], args->value[1], args->value[2],
                    args->value[3]);
}
static void st_C_interface_copy_population(void *a0, void *a1, void *a2,
                                           void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_copy_population, &args);
}
static void st_C_body_correct_missing(void *payload)
{
    st_C_arguments *args = payload;
    correct_missing(args->value[0], args->value[1], args->value[2],
                    args->value[3], args->value[4], args->value[5],
                    args->value[6], args->value[7]);
}
static void st_C_interface_correct_missing(void *a0, void *a1, void *a2,
                                           void *a3, void *a4, void *a5,
                                           void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_correct_missing, &args);
}
static void st_C_body_cross(void *payload)
{
    st_C_arguments *args = payload;
    cross(args->value[0], args->value[1], args->value[2], args->value[3],
          args->value[4], args->value[5], args->value[6], args->value[7]);
}
static void st_C_interface_cross(void *a0, void *a1, void *a2, void *a3,
                                 void *a4, void *a5, void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_cross, &args);
}
static void st_C_body_dh(void *payload)
{
    st_C_arguments *args = payload;
    dh(args->value[0], args->value[1], args->value[2], args->value[3]);
}
static void st_C_interface_dh(void *a0, void *a1, void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_dh, &args);
}
static void st_C_body_divide_population(void *payload)
{
    st_C_arguments *args = payload;
    divide_population(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_divide_population(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_divide_population, &args);
}
static void st_C_body_effect_weight_set_all(void *payload)
{
    st_C_arguments *args = payload;
    effect_weight_set_all(args->value[0]);
}
static void st_C_interface_effect_weight_set_all(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_effect_weight_set_all, &args);
}
static void st_C_body_effmap_remove(void *payload)
{
    st_C_arguments *args = payload;
    effmap_remove(args->value[0]);
}
static void st_C_interface_effmap_remove(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_effmap_remove, &args);
}
static void st_C_body_effmap_remove_all(void *payload)
{
    (void)payload;
    effmap_remove_all();
}
static void st_C_interface_effmap_remove_all(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_effmap_remove_all, &args);
}
static void st_C_body_evaluate_genome(void *payload)
{
    st_C_arguments *args = payload;
    evaluate_genome(args->value[0], args->value[1], args->value[2],
                    args->value[3], args->value[4], args->value[5]);
}
static void st_C_interface_evaluate_genome(void *a0, void *a1, void *a2,
                                           void *a3, void *a4, void *a5)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5}};
    st_work_native_interface(st_C_body_evaluate_genome, &args);
}
static void st_C_body_evaluate_genome_region(void *payload)
{
    st_C_arguments *args = payload;
    evaluate_genome_region(args->value[0], args->value[1], args->value[2],
                           args->value[3], args->value[4], args->value[5],
                           args->value[6], args->value[7], args->value[8]);
}
static void st_C_interface_evaluate_genome_region(void *a0, void *a1, void *a2,
                                                  void *a3, void *a4, void *a5,
                                                  void *a6, void *a7, void *a8)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8}};
    st_work_native_interface(st_C_body_evaluate_genome_region, &args);
}
static void st_C_body_evaluate_ld(void *payload)
{
    st_C_arguments *args = payload;
    evaluate_ld(args->value[0], args->value[1], args->value[2], args->value[3],
                args->value[4]);
}
static void st_C_interface_evaluate_ld(void *a0, void *a1, void *a2, void *a3,
                                       void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_evaluate_ld, &args);
}
static void st_C_body_evaluate_locus(void *payload)
{
    st_C_arguments *args = payload;
    evaluate_locus(args->value[0], args->value[1], args->value[2],
                   args->value[3], args->value[4]);
}
static void st_C_interface_evaluate_locus(void *a0, void *a1, void *a2,
                                          void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_evaluate_locus, &args);
}
static void st_C_body_evaluate_mdp(void *payload)
{
    st_C_arguments *args = payload;
    evaluate_mdp(args->value[0], args->value[1], args->value[2],
                 args->value[3]);
}
static void st_C_interface_evaluate_mdp(void *a0, void *a1, void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_evaluate_mdp, &args);
}
static void st_C_body_evaluate_population(void *payload)
{
    st_C_arguments *args = payload;
    evaluate_population(args->value[0], args->value[1]);
}
static void st_C_interface_evaluate_population(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_evaluate_population, &args);
}
static void st_C_body_st_r_info(void *payload)
{
    st_C_arguments *args = payload;
    st_r_info(args->value[0], args->value[1]);
}
static void st_C_interface_st_r_info(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_st_r_info, &args);
}
static void st_C_body_gen_dist(void *payload)
{
    st_C_arguments *args = payload;
    gen_dist(args->value[0], args->value[1], args->value[2], args->value[3],
             args->value[4], args->value[5], args->value[6], args->value[7],
             args->value[8], args->value[9], args->value[10], args->value[11],
             args->value[12], args->value[13], args->value[14], args->value[15],
             args->value[16]);
}
static void st_C_interface_gen_dist(void *a0, void *a1, void *a2, void *a3,
                                    void *a4, void *a5, void *a6, void *a7,
                                    void *a8, void *a9, void *a10, void *a11,
                                    void *a12, void *a13, void *a14, void *a15,
                                    void *a16)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11,
                            a12, a13, a14, a15, a16}};
    st_work_native_interface(st_C_body_gen_dist, &args);
}
static void st_C_body_generate_effect_file(void *payload)
{
    st_C_arguments *args = payload;
    generate_effect_file(args->value[0], args->value[1]);
}
static void st_C_interface_generate_effect_file(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_generate_effect_file, &args);
}
static void st_C_body_generate_map_file(void *payload)
{
    st_C_arguments *args = payload;
    generate_map_file(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_generate_map_file(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_generate_map_file, &args);
}
static void st_C_body_generate_population(void *payload)
{
    st_C_arguments *args = payload;
    generate_population(args->value[0], args->value[1], args->value[2],
                        args->value[3], args->value[4], args->value[5],
                        args->value[6], args->value[7], args->value[8],
                        args->value[9], args->value[10], args->value[11],
                        args->value[12], args->value[13], args->value[14],
                        args->value[15], args->value[16]);
}
static void st_C_interface_generate_population(void *a0, void *a1, void *a2,
                                               void *a3, void *a4, void *a5,
                                               void *a6, void *a7, void *a8,
                                               void *a9, void *a10, void *a11,
                                               void *a12, void *a13, void *a14,
                                               void *a15, void *a16)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11,
                            a12, a13, a14, a15, a16}};
    st_work_native_interface(st_C_body_generate_population, &args);
}
static void st_C_body_genotype_population(void *payload)
{
    st_C_arguments *args = payload;
    genotype_population(args->value[0]);
}
static void st_C_interface_genotype_population(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_genotype_population, &args);
}
static void st_C_body_get_allele_number(void *payload)
{
    st_C_arguments *args = payload;
    get_allele_number(args->value[0], args->value[1]);
}
static void st_C_interface_get_allele_number(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_get_allele_number, &args);
}
static void st_C_body_get_allow_zeros(void *payload)
{
    st_C_arguments *args = payload;
    get_allow_zeros(args->value[0]);
}
static void st_C_interface_get_allow_zeros(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_get_allow_zeros, &args);
}
static void st_C_body_get_genome_par(void *payload)
{
    st_C_arguments *args = payload;
    get_genome_par(args->value[0], args->value[1]);
}
static void st_C_interface_get_genome_par(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_get_genome_par, &args);
}
static void st_C_body_get_genome_par_b(void *payload)
{
    st_C_arguments *args = payload;
    get_genome_par_b(args->value[0]);
}
static void st_C_interface_get_genome_par_b(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_get_genome_par_b, &args);
}
static void st_C_body_get_info_level(void *payload)
{
    st_C_arguments *args = payload;
    get_info_level(args->value[0]);
}
static void st_C_interface_get_info_level(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_get_info_level, &args);
}
static void st_C_body_get_map_a(void *payload)
{
    st_C_arguments *args = payload;
    get_map_a(args->value[0], args->value[1], args->value[2], args->value[3]);
}
static void st_C_interface_get_map_a(void *a0, void *a1, void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_get_map_a, &args);
}
static void st_C_body_get_map_b(void *payload)
{
    st_C_arguments *args = payload;
    get_map_b(args->value[0], args->value[1], args->value[2], args->value[3]);
}
static void st_C_interface_get_map_b(void *a0, void *a1, void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_get_map_b, &args);
}
static void st_C_body_get_mdp(void *payload)
{
    st_C_arguments *args = payload;
    get_mdp(args->value[0]);
}
static void st_C_interface_get_mdp(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_get_mdp, &args);
}
static void st_C_body_get_population(void *payload)
{
    st_C_arguments *args = payload;
    get_population(args->value[0], args->value[1], args->value[2],
                   args->value[3], args->value[4], args->value[5]);
}
static void st_C_interface_get_population(void *a0, void *a1, void *a2,
                                          void *a3, void *a4, void *a5)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5}};
    st_work_native_interface(st_C_body_get_population, &args);
}
static void st_C_body_get_population_info(void *payload)
{
    st_C_arguments *args = payload;
    get_population_info(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_get_population_info(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_get_population_info, &args);
}
static void st_C_body_get_population_size(void *payload)
{
    st_C_arguments *args = payload;
    get_population_size(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_get_population_size(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_get_population_size, &args);
}
static void st_C_body_get_population_value(void *payload)
{
    st_C_arguments *args = payload;
    get_population_value(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4]);
}
static void st_C_interface_get_population_value(void *a0, void *a1, void *a2,
                                                void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_get_population_value, &args);
}
static void st_C_body_gs_build_V_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_build_V_01_GV(args->value[0], args->value[1], args->value[2],
                     args->value[3]);
}
static void st_C_interface_gs_build_V_01_GV(void *a0, void *a1, void *a2,
                                            void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_build_V_01_GV, &args);
}
static void st_C_body_gs_build_Z_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_build_Z_01_GV(args->value[0], args->value[1], args->value[2],
                     args->value[3]);
}
static void st_C_interface_gs_build_Z_01_GV(void *a0, void *a1, void *a2,
                                            void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_build_Z_01_GV, &args);
}
static void st_C_body_gs_build_esvs_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_build_esvs_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3], args->value[4], args->value[5],
                        args->value[6], args->value[7], args->value[8],
                        args->value[9]);
}
static void st_C_interface_gs_build_esvs_01_GV(void *a0, void *a1, void *a2,
                                               void *a3, void *a4, void *a5,
                                               void *a6, void *a7, void *a8,
                                               void *a9)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8, a9}};
    st_work_native_interface(st_C_body_gs_build_esvs_01_GV, &args);
}
static void st_C_body_gs_calc_ld_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_calc_ld_01_GV(args->value[0], args->value[1], args->value[2],
                     args->value[3], args->value[4]);
}
static void st_C_interface_gs_calc_ld_01_GV(void *a0, void *a1, void *a2,
                                            void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_calc_ld_01_GV, &args);
}
static void st_C_body_gs_calc_q_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_calc_q_01_GV(args->value[0], args->value[1], args->value[2],
                    args->value[3], args->value[4], args->value[5]);
}
static void st_C_interface_gs_calc_q_01_GV(void *a0, void *a1, void *a2,
                                           void *a3, void *a4, void *a5)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5}};
    st_work_native_interface(st_C_body_gs_calc_q_01_GV, &args);
}
static void st_C_body_gs_calc_rf_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_calc_rf_01_GV(args->value[0], args->value[1], args->value[2],
                     args->value[3], args->value[4]);
}
static void st_C_interface_gs_calc_rf_01_GV(void *a0, void *a1, void *a2,
                                            void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_calc_rf_01_GV, &args);
}
static void st_C_body_gs_check_pvals_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_check_pvals_01_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_check_pvals_01_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_check_pvals_01_GV, &args);
}
static void st_C_body_gs_chrom_stats_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_chrom_stats_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3]);
}
static void st_C_interface_gs_chrom_stats_01_GV(void *a0, void *a1, void *a2,
                                                void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_chrom_stats_01_GV, &args);
}
static void st_C_body_gs_copy_marker_data_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_copy_marker_data_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_copy_marker_data_01_GV(void *a0, void *a1,
                                                     void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_copy_marker_data_01_GV, &args);
}
static void st_C_body_gs_cross_eval_es_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_eval_es_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_cross_eval_es_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_cross_eval_es_01_GV, &args);
}
static void st_C_body_gs_cross_eval_gd_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_eval_gd_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_cross_eval_gd_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_cross_eval_gd_01_GV, &args);
}
static void st_C_body_gs_cross_eval_gd_fct_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_eval_gd_fct_01_GV(args->value[0], args->value[1], args->value[2],
                               args->value[3]);
}
static void st_C_interface_gs_cross_eval_gd_fct_01_GV(void *a0, void *a1,
                                                      void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_cross_eval_gd_fct_01_GV, &args);
}
static void st_C_body_gs_cross_eval_ma_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_eval_ma_01_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_cross_eval_ma_01_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_cross_eval_ma_01_GV, &args);
}
static void st_C_body_gs_cross_eval_mi_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_eval_mi_01_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_cross_eval_mi_01_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_cross_eval_mi_01_GV, &args);
}
static void st_C_body_gs_cross_eval_mu_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_eval_mu_01_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_cross_eval_mu_01_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_cross_eval_mu_01_GV, &args);
}
static void st_C_body_gs_cross_eval_va_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_eval_va_01_GV(args->value[0], args->value[1], args->value[2],
                           args->value[3], args->value[4]);
}
static void st_C_interface_gs_cross_eval_va_01_GV(void *a0, void *a1, void *a2,
                                                  void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_cross_eval_va_01_GV, &args);
}
static void st_C_body_gs_cross_info_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_info_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3], args->value[4]);
}
static void st_C_interface_gs_cross_info_01_GV(void *a0, void *a1, void *a2,
                                               void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_cross_info_01_GV, &args);
}
static void st_C_body_st_unavailable_gs_cross_info_gd_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_gs_cross_info_gd_01_GV(args->value[0], args->value[1],
                                          args->value[2]);
}
static void st_C_interface_st_unavailable_gs_cross_info_gd_01_GV(void *a0,
                                                                 void *a1,
                                                                 void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_st_unavailable_gs_cross_info_gd_01_GV,
                             &args);
}
static void st_C_body_gs_cross_validation_02_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_cross_validation_02_GV(args->value[0], args->value[1], args->value[2],
                              args->value[3], args->value[4], args->value[5],
                              args->value[6], args->value[7], args->value[8],
                              args->value[9], args->value[10], args->value[11]);
}
static void st_C_interface_gs_cross_validation_02_GV(
    void *a0, void *a1, void *a2, void *a3, void *a4, void *a5, void *a6,
    void *a7, void *a8, void *a9, void *a10, void *a11)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11}};
    st_work_native_interface(st_C_body_gs_cross_validation_02_GV, &args);
}
static void st_C_body_gs_def_hblocks_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_def_hblocks_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6], args->value[7], args->value[8],
                         args->value[9]);
}
static void st_C_interface_gs_def_hblocks_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6, void *a7, void *a8,
                                                void *a9)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8, a9}};
    st_work_native_interface(st_C_body_gs_def_hblocks_01_GV, &args);
}
static void st_C_body_gs_esteff_lsq_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_esteff_lsq_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3]);
}
static void st_C_interface_gs_esteff_lsq_01_GV(void *a0, void *a1, void *a2,
                                               void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_esteff_lsq_01_GV, &args);
}
static void st_C_body_gs_esteff_rmla_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_esteff_rmla_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6], args->value[7]);
}
static void st_C_interface_gs_esteff_rmla_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_gs_esteff_rmla_01_GV, &args);
}
static void st_C_body_gs_esteff_rmlc_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_esteff_rmlc_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6]);
}
static void st_C_interface_gs_esteff_rmlc_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_esteff_rmlc_01_GV, &args);
}
static void st_C_body_gs_esteff_rmlr_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_esteff_rmlr_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6]);
}
static void st_C_interface_gs_esteff_rmlr_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_esteff_rmlr_01_GV, &args);
}
static void st_C_body_gs_esteff_rmlv_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_esteff_rmlv_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6]);
}
static void st_C_interface_gs_esteff_rmlv_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_esteff_rmlv_01_GV, &args);
}
static void st_C_body_gs_esteff_rr_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_esteff_rr_01_GV(args->value[0], args->value[1], args->value[2],
                       args->value[3], args->value[4], args->value[5],
                       args->value[6], args->value[7], args->value[8]);
}
static void st_C_interface_gs_esteff_rr_01_GV(void *a0, void *a1, void *a2,
                                              void *a3, void *a4, void *a5,
                                              void *a6, void *a7, void *a8)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8}};
    st_work_native_interface(st_C_body_gs_esteff_rr_01_GV, &args);
}
static void st_C_body_gs_estimate_bv_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_estimate_bv_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4]);
}
static void st_C_interface_gs_estimate_bv_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_estimate_bv_01_GV, &args);
}
static void st_C_body_gs_get_info_level_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_get_info_level_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_get_info_level_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_get_info_level_GV, &args);
}
static void st_C_body_gs_get_map_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_get_map_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_get_map_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_get_map_01_GV, &args);
}
static void st_C_body_gs_get_num_threads(void *payload)
{
    st_C_arguments *args = payload;
    gs_get_num_threads(args->value[0]);
}
static void st_C_interface_gs_get_num_threads(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_gs_get_num_threads, &args);
}
static void st_C_body_gs_lambda_aov_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_aov_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3], args->value[4], args->value[5]);
}
static void st_C_interface_gs_lambda_aov_01_GV(void *a0, void *a1, void *a2,
                                               void *a3, void *a4, void *a5)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5}};
    st_work_native_interface(st_C_body_gs_lambda_aov_01_GV, &args);
}
static void st_C_body_gs_lambda_const_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_const_01_GV(args->value[0], args->value[1], args->value[2],
                          args->value[3], args->value[4]);
}
static void st_C_interface_gs_lambda_const_01_GV(void *a0, void *a1, void *a2,
                                                 void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_lambda_const_01_GV, &args);
}
static void st_C_body_gs_lambda_const_02_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_const_02_GV(args->value[0], args->value[1], args->value[2],
                          args->value[3], args->value[4]);
}
static void st_C_interface_gs_lambda_const_02_GV(void *a0, void *a1, void *a2,
                                                 void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_lambda_const_02_GV, &args);
}
static void st_C_body_gs_lambda_emstep_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_emstep_01_GV(args->value[0], args->value[1], args->value[2],
                           args->value[3], args->value[4]);
}
static void st_C_interface_gs_lambda_emstep_01_GV(void *a0, void *a1, void *a2,
                                                  void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_lambda_emstep_01_GV, &args);
}
static void st_C_body_gs_lambda_reg_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_reg_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3], args->value[4]);
}
static void st_C_interface_gs_lambda_reg_01_GV(void *a0, void *a1, void *a2,
                                               void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_lambda_reg_01_GV, &args);
}
static void st_C_body_gs_lambda_rmla_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_rmla_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6], args->value[7]);
}
static void st_C_interface_gs_lambda_rmla_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_gs_lambda_rmla_01_GV, &args);
}
static void st_C_body_gs_lambda_rmla_02_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_rmla_02_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6], args->value[7]);
}
static void st_C_interface_gs_lambda_rmla_02_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_gs_lambda_rmla_02_GV, &args);
}
static void st_C_body_gs_lambda_rmlc_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_rmlc_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6]);
}
static void st_C_interface_gs_lambda_rmlc_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_lambda_rmlc_01_GV, &args);
}
static void st_C_body_gs_lambda_rmlr_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_rmlr_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6]);
}
static void st_C_interface_gs_lambda_rmlr_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_lambda_rmlr_01_GV, &args);
}
static void st_C_body_st_unavailable_gs_lambda_rmlv_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_gs_lambda_rmlv_01_GV(args->value[0], args->value[1],
                                        args->value[2], args->value[3],
                                        args->value[4], args->value[5]);
}
static void
st_C_interface_st_unavailable_gs_lambda_rmlv_01_GV(void *a0, void *a1, void *a2,
                                                   void *a3, void *a4, void *a5)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5}};
    st_work_native_interface(st_C_body_st_unavailable_gs_lambda_rmlv_01_GV,
                             &args);
}
static void st_C_body_gs_lambda_rrblup_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_rrblup_01_GV(args->value[0], args->value[1], args->value[2],
                           args->value[3], args->value[4], args->value[5],
                           args->value[6]);
}
static void st_C_interface_gs_lambda_rrblup_01_GV(void *a0, void *a1, void *a2,
                                                  void *a3, void *a4, void *a5,
                                                  void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_lambda_rrblup_01_GV, &args);
}
static void st_C_body_gs_lambda_var_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_lambda_var_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3], args->value[4], args->value[5]);
}
static void st_C_interface_gs_lambda_var_01_GV(void *a0, void *a1, void *a2,
                                               void *a3, void *a4, void *a5)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5}};
    st_work_native_interface(st_C_body_gs_lambda_var_01_GV, &args);
}
static void st_C_body_gs_marker_stats_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_marker_stats_01_GV(args->value[0], args->value[1], args->value[2],
                          args->value[3], args->value[4], args->value[5],
                          args->value[6], args->value[7], args->value[8]);
}
static void st_C_interface_gs_marker_stats_01_GV(void *a0, void *a1, void *a2,
                                                 void *a3, void *a4, void *a5,
                                                 void *a6, void *a7, void *a8)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8}};
    st_work_native_interface(st_C_body_gs_marker_stats_01_GV, &args);
}
static void st_C_body_gs_mme_coeff_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mme_coeff_01_GV(args->value[0], args->value[1], args->value[2],
                       args->value[3]);
}
static void st_C_interface_gs_mme_coeff_01_GV(void *a0, void *a1, void *a2,
                                              void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_mme_coeff_01_GV, &args);
}
static void st_C_body_st_unavailable_gs_mme_invcoeff_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_gs_mme_invcoeff_01_GV(args->value[0], args->value[1],
                                         args->value[2]);
}
static void st_C_interface_st_unavailable_gs_mme_invcoeff_01_GV(void *a0,
                                                                void *a1,
                                                                void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_st_unavailable_gs_mme_invcoeff_01_GV,
                             &args);
}
static void st_C_body_gs_mme_restcoeff_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mme_restcoeff_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_mme_restcoeff_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_mme_restcoeff_01_GV, &args);
}
static void st_C_body_gs_mme_restrhs_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mme_restrhs_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_mme_restrhs_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_mme_restrhs_01_GV, &args);
}
static void st_C_body_gs_mme_rhs_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mme_rhs_01_GV(args->value[0], args->value[1], args->value[2],
                     args->value[3]);
}
static void st_C_interface_gs_mme_rhs_01_GV(void *a0, void *a1, void *a2,
                                            void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_mme_rhs_01_GV, &args);
}
static void st_C_body_gs_mme_solve_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mme_solve_01_GV(args->value[0], args->value[1], args->value[2],
                       args->value[3]);
}
static void st_C_interface_gs_mme_solve_01_GV(void *a0, void *a1, void *a2,
                                              void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_mme_solve_01_GV, &args);
}
static void st_C_body_gs_mmet_coeff_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mmet_coeff_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3]);
}
static void st_C_interface_gs_mmet_coeff_01_GV(void *a0, void *a1, void *a2,
                                               void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_mmet_coeff_01_GV, &args);
}
static void st_C_body_st_unavailable_gs_mmet_coeff_03_GV(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_gs_mmet_coeff_03_GV(args->value[0], args->value[1],
                                       args->value[2], args->value[3]);
}
static void st_C_interface_st_unavailable_gs_mmet_coeff_03_GV(void *a0,
                                                              void *a1,
                                                              void *a2,
                                                              void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_st_unavailable_gs_mmet_coeff_03_GV,
                             &args);
}
static void st_C_body_gs_mmet_rhs_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mmet_rhs_01_GV(args->value[0], args->value[1], args->value[2],
                      args->value[3]);
}
static void st_C_interface_gs_mmet_rhs_01_GV(void *a0, void *a1, void *a2,
                                             void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_mmet_rhs_01_GV, &args);
}
static void st_C_body_gs_mmet_solve_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_mmet_solve_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3]);
}
static void st_C_interface_gs_mmet_solve_01_GV(void *a0, void *a1, void *a2,
                                               void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_mmet_solve_01_GV, &args);
}
static void st_C_body_gs_neg_effall_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_neg_effall_01_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_neg_effall_01_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_neg_effall_01_GV, &args);
}
static void st_C_body_gs_pos_effall_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_pos_effall_01_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_pos_effall_01_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_pos_effall_01_GV, &args);
}
static void st_C_body_gs_r_info(void *payload)
{
    st_C_arguments *args = payload;
    gs_r_info(args->value[0], args->value[1]);
}
static void st_C_interface_gs_r_info(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_r_info, &args);
}
static void st_C_body_gs_read_map_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_read_map_01_GV(args->value[0], args->value[1], args->value[2],
                      args->value[3], args->value[4], args->value[5],
                      args->value[6], args->value[7]);
}
static void st_C_interface_gs_read_map_01_GV(void *a0, void *a1, void *a2,
                                             void *a3, void *a4, void *a5,
                                             void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_gs_read_map_01_GV, &args);
}
static void st_C_body_gs_read_marker_data_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_read_marker_data_01_GV(args->value[0], args->value[1], args->value[2],
                              args->value[3], args->value[4], args->value[5],
                              args->value[6]);
}
static void st_C_interface_gs_read_marker_data_01_GV(void *a0, void *a1,
                                                     void *a2, void *a3,
                                                     void *a4, void *a5,
                                                     void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_read_marker_data_01_GV, &args);
}
static void st_C_body_gs_read_marker_data_02_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_read_marker_data_02_GV(args->value[0], args->value[1], args->value[2],
                              args->value[3], args->value[4], args->value[5],
                              args->value[6]);
}
static void st_C_interface_gs_read_marker_data_02_GV(void *a0, void *a1,
                                                     void *a2, void *a3,
                                                     void *a4, void *a5,
                                                     void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_read_marker_data_02_GV, &args);
}
static void st_C_body_gs_read_marker_data_03_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_read_marker_data_03_GV(args->value[0], args->value[1], args->value[2],
                              args->value[3], args->value[4], args->value[5],
                              args->value[6]);
}
static void st_C_interface_gs_read_marker_data_03_GV(void *a0, void *a1,
                                                     void *a2, void *a3,
                                                     void *a4, void *a5,
                                                     void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_read_marker_data_03_GV, &args);
}
static void st_C_body_gs_read_marker_data_04_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_read_marker_data_04_GV(args->value[0], args->value[1], args->value[2],
                              args->value[3], args->value[4], args->value[5],
                              args->value[6]);
}
static void st_C_interface_gs_read_marker_data_04_GV(void *a0, void *a1,
                                                     void *a2, void *a3,
                                                     void *a4, void *a5,
                                                     void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_read_marker_data_04_GV, &args);
}
static void st_C_body_gs_read_performance_data_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_read_performance_data_01_GV(args->value[0], args->value[1],
                                   args->value[2], args->value[3],
                                   args->value[4]);
}
static void st_C_interface_gs_read_performance_data_01_GV(void *a0, void *a1,
                                                          void *a2, void *a3,
                                                          void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_read_performance_data_01_GV, &args);
}
static void st_C_body_gs_recode_hbc_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_recode_hbc_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_recode_hbc_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_recode_hbc_01_GV, &args);
}
static void st_C_body_gs_recode_hil_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_recode_hil_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3]);
}
static void st_C_interface_gs_recode_hil_01_GV(void *a0, void *a1, void *a2,
                                               void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_recode_hil_01_GV, &args);
}
static void st_C_body_gs_recode_ref_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_recode_ref_01_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3]);
}
static void st_C_interface_gs_recode_ref_01_GV(void *a0, void *a1, void *a2,
                                               void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_recode_ref_01_GV, &args);
}
static void st_C_body_gs_recode_ref_02_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_recode_ref_02_GV(args->value[0], args->value[1], args->value[2],
                        args->value[3], args->value[4]);
}
static void st_C_interface_gs_recode_ref_02_GV(void *a0, void *a1, void *a2,
                                               void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_recode_ref_02_GV, &args);
}
static void st_C_body_gs_recode_ref_name_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_recode_ref_name_01_GV(args->value[0], args->value[1], args->value[2],
                             args->value[3]);
}
static void st_C_interface_gs_recode_ref_name_01_GV(void *a0, void *a1,
                                                    void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_recode_ref_name_01_GV, &args);
}
static void st_C_body_st_unavailable_gs_reset(void *payload)
{
    (void)payload;
    st_unavailable_gs_reset();
}
static void st_C_interface_st_unavailable_gs_reset(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_st_unavailable_gs_reset, &args);
}
static void
st_C_body_st_unavailable_gs_restrict_marker_data_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_gs_restrict_marker_data_01_GV(
        args->value[0], args->value[1], args->value[2], args->value[3],
        args->value[4], args->value[5], args->value[6], args->value[7]);
}
static void st_C_interface_st_unavailable_gs_restrict_marker_data_01_GV(
    void *a0, void *a1, void *a2, void *a3, void *a4, void *a5, void *a6,
    void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(
        st_C_body_st_unavailable_gs_restrict_marker_data_01_GV, &args);
}
static void st_C_body_gs_restrict_marker_data_03_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_restrict_marker_data_03_GV(
        args->value[0], args->value[1], args->value[2], args->value[3],
        args->value[4], args->value[5], args->value[6], args->value[7]);
}
static void st_C_interface_gs_restrict_marker_data_03_GV(void *a0, void *a1,
                                                         void *a2, void *a3,
                                                         void *a4, void *a5,
                                                         void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_gs_restrict_marker_data_03_GV, &args);
}
static void st_C_body_gs_return_effects_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_return_effects_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_return_effects_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_return_effects_01_GV, &args);
}
static void st_C_body_gs_return_performance_data_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_return_performance_data_01_GV(args->value[0], args->value[1],
                                     args->value[2], args->value[3]);
}
static void st_C_interface_gs_return_performance_data_01_GV(void *a0, void *a1,
                                                            void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_return_performance_data_01_GV, &args);
}
static void st_C_body_gs_return_pvals_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_return_pvals_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_return_pvals_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_return_pvals_01_GV, &args);
}
static void st_C_body_gs_set_all_info_levels(void *payload)
{
    st_C_arguments *args = payload;
    gs_set_all_info_levels(args->value[0]);
}
static void st_C_interface_gs_set_all_info_levels(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_gs_set_all_info_levels, &args);
}
static void st_C_body_gs_set_allele_codes_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_set_allele_codes_01_GV(args->value[0], args->value[1], args->value[2],
                              args->value[3], args->value[4], args->value[5],
                              args->value[6]);
}
static void st_C_interface_gs_set_allele_codes_01_GV(void *a0, void *a1,
                                                     void *a2, void *a3,
                                                     void *a4, void *a5,
                                                     void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_set_allele_codes_01_GV, &args);
}
static void st_C_body_gs_set_effects_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_set_effects_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_set_effects_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_set_effects_01_GV, &args);
}
static void st_C_body_gs_set_hblocks_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_set_hblocks_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6]);
}
static void st_C_interface_gs_set_hblocks_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4, void *a5,
                                                void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_set_hblocks_01_GV, &args);
}
static void st_C_body_gs_set_info_level_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_set_info_level_GV(args->value[0], args->value[1]);
}
static void st_C_interface_gs_set_info_level_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_gs_set_info_level_GV, &args);
}
static void st_C_body_gs_set_num_threads(void *payload)
{
    st_C_arguments *args = payload;
    gs_set_num_threads(args->value[0]);
}
static void st_C_interface_gs_set_num_threads(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_gs_set_num_threads, &args);
}
static void st_C_body_st_unavailable_gs_set_openblas_threads(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_gs_set_openblas_threads(args->value[0]);
}
static void st_C_interface_st_unavailable_gs_set_openblas_threads(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_st_unavailable_gs_set_openblas_threads,
                             &args);
}
static void st_C_body_gs_set_performance_data_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_set_performance_data_01_GV(args->value[0], args->value[1],
                                  args->value[2], args->value[3]);
}
static void st_C_interface_gs_set_performance_data_01_GV(void *a0, void *a1,
                                                         void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_set_performance_data_01_GV, &args);
}
static void st_C_body_gs_single_marker_aov_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_single_marker_aov_01_GV(args->value[0], args->value[1], args->value[2],
                               args->value[3], args->value[4]);
}
static void st_C_interface_gs_single_marker_aov_01_GV(void *a0, void *a1,
                                                      void *a2, void *a3,
                                                      void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_single_marker_aov_01_GV, &args);
}
static void st_C_body_gs_single_marker_reg_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_single_marker_reg_01_GV(args->value[0], args->value[1], args->value[2],
                               args->value[3]);
}
static void st_C_interface_gs_single_marker_reg_01_GV(void *a0, void *a1,
                                                      void *a2, void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_gs_single_marker_reg_01_GV, &args);
}
static void st_C_body_st_unavailable_gs_test_mp(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_gs_test_mp(args->value[0]);
}
static void st_C_interface_st_unavailable_gs_test_mp(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_st_unavailable_gs_test_mp, &args);
}
static void st_C_body_gs_testeff_lsq_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_testeff_lsq_01_GV(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4]);
}
static void st_C_interface_gs_testeff_lsq_01_GV(void *a0, void *a1, void *a2,
                                                void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_testeff_lsq_01_GV, &args);
}
static void st_C_body_gs_testeff_rmlc_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_testeff_rmlc_01_GV(args->value[0], args->value[1], args->value[2],
                          args->value[3], args->value[4], args->value[5],
                          args->value[6], args->value[7]);
}
static void st_C_interface_gs_testeff_rmlc_01_GV(void *a0, void *a1, void *a2,
                                                 void *a3, void *a4, void *a5,
                                                 void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_gs_testeff_rmlc_01_GV, &args);
}
static void st_C_body_gs_testeff_rmlv_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_testeff_rmlv_01_GV(args->value[0], args->value[1], args->value[2],
                          args->value[3], args->value[4], args->value[5],
                          args->value[6], args->value[7]);
}
static void st_C_interface_gs_testeff_rmlv_01_GV(void *a0, void *a1, void *a2,
                                                 void *a3, void *a4, void *a5,
                                                 void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_gs_testeff_rmlv_01_GV, &args);
}
static void st_C_body_gs_vc_rrblup_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_vc_rrblup_01_GV(args->value[0], args->value[1], args->value[2],
                       args->value[3], args->value[4], args->value[5],
                       args->value[6]);
}
static void st_C_interface_gs_vc_rrblup_01_GV(void *a0, void *a1, void *a2,
                                              void *a3, void *a4, void *a5,
                                              void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_gs_vc_rrblup_01_GV, &args);
}
static void st_C_body_gs_w_aov_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_w_aov_01_GV(args->value[0], args->value[1], args->value[2],
                   args->value[3], args->value[4]);
}
static void st_C_interface_gs_w_aov_01_GV(void *a0, void *a1, void *a2,
                                          void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_gs_w_aov_01_GV, &args);
}
static void st_C_body_gs_write_map_01_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_write_map_01_GV(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_gs_write_map_01_GV(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_gs_write_map_01_GV, &args);
}
static void st_C_body_gs_write_marker_data_02_GV(void *payload)
{
    st_C_arguments *args = payload;
    gs_write_marker_data_02_GV(args->value[0], args->value[1], args->value[2],
                               args->value[3], args->value[4], args->value[5],
                               args->value[6], args->value[7], args->value[8],
                               args->value[9], args->value[10], args->value[11],
                               args->value[12], args->value[13]);
}
static void st_C_interface_gs_write_marker_data_02_GV(
    void *a0, void *a1, void *a2, void *a3, void *a4, void *a5, void *a6,
    void *a7, void *a8, void *a9, void *a10, void *a11, void *a12, void *a13)
{
    st_C_arguments args = {
        {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13}};
    st_work_native_interface(st_C_body_gs_write_marker_data_02_GV, &args);
}
static void st_C_body_haplotype_calc_d(void *payload)
{
    st_C_arguments *args = payload;
    haplotype_calc_d(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_haplotype_calc_d(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_haplotype_calc_d, &args);
}
static void st_C_body_haplotype_free(void *payload)
{
    (void)payload;
    haplotype_free();
}
static void st_C_interface_haplotype_free(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_haplotype_free, &args);
}
static void st_C_body_haplotype_init(void *payload)
{
    st_C_arguments *args = payload;
    haplotype_init(args->value[0], args->value[1], args->value[2],
                   args->value[3]);
}
static void st_C_interface_haplotype_init(void *a0, void *a1, void *a2,
                                          void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_haplotype_init, &args);
}
static void st_C_body_haplotype_next(void *payload)
{
    st_C_arguments *args = payload;
    haplotype_next(args->value[0], args->value[1], args->value[2],
                   args->value[3], args->value[4], args->value[5],
                   args->value[6]);
}
static void st_C_interface_haplotype_next(void *a0, void *a1, void *a2,
                                          void *a3, void *a4, void *a5,
                                          void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_haplotype_next, &args);
}
static void st_C_body_il_eval_library_GV(void *payload)
{
    st_C_arguments *args = payload;
    il_eval_library_GV(args->value[0], args->value[1], args->value[2],
                       args->value[3], args->value[4], args->value[5]);
}
static void st_C_interface_il_eval_library_GV(void *a0, void *a1, void *a2,
                                              void *a3, void *a4, void *a5)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5}};
    st_work_native_interface(st_C_body_il_eval_library_GV, &args);
}
static void st_C_body_il_eval_lines_GV(void *payload)
{
    st_C_arguments *args = payload;
    il_eval_lines_GV(args->value[0], args->value[1], args->value[2],
                     args->value[3], args->value[4], args->value[5],
                     args->value[6], args->value[7], args->value[8],
                     args->value[9], args->value[10]);
}
static void st_C_interface_il_eval_lines_GV(void *a0, void *a1, void *a2,
                                            void *a3, void *a4, void *a5,
                                            void *a6, void *a7, void *a8,
                                            void *a9, void *a10)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10}};
    st_work_native_interface(st_C_body_il_eval_lines_GV, &args);
}
static void st_C_body_il_get_population_size_GV(void *payload)
{
    st_C_arguments *args = payload;
    il_get_population_size_GV(args->value[0], args->value[1]);
}
static void st_C_interface_il_get_population_size_GV(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_il_get_population_size_GV, &args);
}
static void st_C_body_init_population(void *payload)
{
    st_C_arguments *args = payload;
    init_population(args->value[0], args->value[1], args->value[2],
                    args->value[3], args->value[4], args->value[5],
                    args->value[6]);
}
static void st_C_interface_init_population(void *a0, void *a1, void *a2,
                                           void *a3, void *a4, void *a5,
                                           void *a6)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    st_work_native_interface(st_C_body_init_population, &args);
}
static void st_C_body_linkage_map_load(void *payload)
{
    st_C_arguments *args = payload;
    linkage_map_load(args->value[0], args->value[1], args->value[2],
                     args->value[3], args->value[4]);
}
static void st_C_interface_linkage_map_load(void *a0, void *a1, void *a2,
                                            void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_linkage_map_load, &args);
}
static void st_C_body_list_eff(void *payload)
{
    st_C_arguments *args = payload;
    list_eff(args->value[0], args->value[1]);
}
static void st_C_interface_list_eff(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_list_eff, &args);
}
static void st_C_body_list_eff_parameter(void *payload)
{
    st_C_arguments *args = payload;
    list_eff_parameter(args->value[0], args->value[1]);
}
static void st_C_interface_list_eff_parameter(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_list_eff_parameter, &args);
}
static void st_C_body_list_populations(void *payload)
{
    st_C_arguments *args = payload;
    list_populations(args->value[0], args->value[1]);
}
static void st_C_interface_list_populations(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_list_populations, &args);
}
static void st_C_body_list_populations_parameter(void *payload)
{
    st_C_arguments *args = payload;
    list_populations_parameter(args->value[0], args->value[1]);
}
static void st_C_interface_list_populations_parameter(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_list_populations_parameter, &args);
}
static void st_C_body_load_effmap(void *payload)
{
    st_C_arguments *args = payload;
    load_effmap(args->value[0], args->value[1]);
}
static void st_C_interface_load_effmap(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_load_effmap, &args);
}
static void st_C_body_load_effmap_mod(void *payload)
{
    st_C_arguments *args = payload;
    load_effmap_mod(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_load_effmap_mod(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_load_effmap_mod, &args);
}
static void st_C_body_optimize_population(void *payload)
{
    st_C_arguments *args = payload;
    optimize_population(args->value[0]);
}
static void st_C_interface_optimize_population(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_optimize_population, &args);
}
static void st_C_body_phenotype_population(void *payload)
{
    st_C_arguments *args = payload;
    phenotype_population(args->value[0], args->value[1], args->value[2],
                         args->value[3]);
}
static void st_C_interface_phenotype_population(void *a0, void *a1, void *a2,
                                                void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_phenotype_population, &args);
}
static void st_C_body_print_version_2(void *payload)
{
    st_C_arguments *args = payload;
    print_version_2(args->value[0]);
}
static void st_C_interface_print_version_2(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_print_version_2, &args);
}
static void st_C_body_r_info(void *payload)
{
    st_C_arguments *args = payload;
    r_info(args->value[0], args->value[1]);
}
static void st_C_interface_r_info(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_r_info, &args);
}
static void st_C_body_remove_evaluate_population(void *payload)
{
    st_C_arguments *args = payload;
    remove_evaluate_population(args->value[0]);
}
static void st_C_interface_remove_evaluate_population(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_remove_evaluate_population, &args);
}
static void st_C_body_remove_genotype_population(void *payload)
{
    st_C_arguments *args = payload;
    remove_genotype_population(args->value[0]);
}
static void st_C_interface_remove_genotype_population(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_remove_genotype_population, &args);
}
static void st_C_body_remove_map(void *payload)
{
    (void)payload;
    remove_map();
}
static void st_C_interface_remove_map(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_remove_map, &args);
}
static void st_C_body_remove_population(void *payload)
{
    st_C_arguments *args = payload;
    remove_population(args->value[0]);
}
static void st_C_interface_remove_population(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_remove_population, &args);
}
static void st_C_body_rename_population(void *payload)
{
    st_C_arguments *args = payload;
    rename_population(args->value[0], args->value[1]);
}
static void st_C_interface_rename_population(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_rename_population, &args);
}
static void st_C_body_reset_mdp(void *payload)
{
    (void)payload;
    reset_mdp();
}
static void st_C_interface_reset_mdp(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_reset_mdp, &args);
}
static void st_C_body_resize_population(void *payload)
{
    st_C_arguments *args = payload;
    resize_population(args->value[0], args->value[1]);
}
static void st_C_interface_resize_population(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_resize_population, &args);
}
static void st_C_body_rng_choose(void *payload)
{
    st_C_arguments *args = payload;
    rng_choose(args->value[0]);
}
static void st_C_interface_rng_choose(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_rng_choose, &args);
}
static void st_C_body_rng_info(void *payload)
{
    (void)payload;
    rng_info();
}
static void st_C_interface_rng_info(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_rng_info, &args);
}
static void st_C_body_rng_init(void *payload)
{
    (void)payload;
    rng_init();
}
static void st_C_interface_rng_init(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_rng_init, &args);
}
static void st_C_body_rng_init_seed(void *payload)
{
    st_C_arguments *args = payload;
    rng_init_seed(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_rng_init_seed(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_rng_init_seed, &args);
}
static void st_C_body_rng_list(void *payload)
{
    (void)payload;
    rng_list();
}
static void st_C_interface_rng_list(void)
{
    st_C_arguments args = {{NULL}};
    st_work_native_interface(st_C_body_rng_list, &args);
}
static void st_C_body_sample_population(void *payload)
{
    st_C_arguments *args = payload;
    sample_population(args->value[0], args->value[1], args->value[2],
                      args->value[3]);
}
static void st_C_interface_sample_population(void *a0, void *a1, void *a2,
                                             void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_sample_population, &args);
}
static void st_C_body_select_all_n_best(void *payload)
{
    st_C_arguments *args = payload;
    select_all_n_best(args->value[0], args->value[1], args->value[2],
                      args->value[3], args->value[4], args->value[5],
                      args->value[6], args->value[7]);
}
static void st_C_interface_select_all_n_best(void *a0, void *a1, void *a2,
                                             void *a3, void *a4, void *a5,
                                             void *a6, void *a7)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4, a5, a6, a7}};
    st_work_native_interface(st_C_body_select_all_n_best, &args);
}
static void st_C_body_set_allow_zeros(void *payload)
{
    st_C_arguments *args = payload;
    set_allow_zeros(args->value[0]);
}
static void st_C_interface_set_allow_zeros(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_set_allow_zeros, &args);
}
static void st_C_body_set_co_freq(void *payload)
{
    st_C_arguments *args = payload;
    set_co_freq(args->value[0]);
}
static void st_C_interface_set_co_freq(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_set_co_freq, &args);
}
static void st_C_body_set_eff_weight(void *payload)
{
    st_C_arguments *args = payload;
    set_eff_weight(args->value[0], args->value[1]);
}
static void st_C_interface_set_eff_weight(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_set_eff_weight, &args);
}
static void st_C_body_set_genome_par(void *payload)
{
    st_C_arguments *args = payload;
    set_genome_par(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_set_genome_par(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_set_genome_par, &args);
}
static void st_C_body_set_info_level(void *payload)
{
    st_C_arguments *args = payload;
    set_info_level(args->value[0]);
}
static void st_C_interface_set_info_level(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_set_info_level, &args);
}
static void st_C_body_set_mdp(void *payload)
{
    st_C_arguments *args = payload;
    set_mdp(args->value[0]);
}
static void st_C_interface_set_mdp(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_set_mdp, &args);
}
static void st_C_body_set_nolociinit(void *payload)
{
    st_C_arguments *args = payload;
    set_nolociinit(args->value[0]);
}
static void st_C_interface_set_nolociinit(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_set_nolociinit, &args);
}
static void st_C_body_set_population_gvalue(void *payload)
{
    st_C_arguments *args = payload;
    set_population_gvalue(args->value[0], args->value[1], args->value[2],
                          args->value[3]);
}
static void st_C_interface_set_population_gvalue(void *a0, void *a1, void *a2,
                                                 void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_set_population_gvalue, &args);
}
static void st_C_body_set_population_info(void *payload)
{
    st_C_arguments *args = payload;
    set_population_info(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_set_population_info(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_set_population_info, &args);
}
static void st_C_body_sort_population(void *payload)
{
    st_C_arguments *args = payload;
    sort_population(args->value[0], args->value[1], args->value[2],
                    args->value[3], args->value[4]);
}
static void st_C_interface_sort_population(void *a0, void *a1, void *a2,
                                           void *a3, void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_sort_population, &args);
}
static void st_C_body_ssd_mating(void *payload)
{
    st_C_arguments *args = payload;
    ssd_mating(args->value[0], args->value[1], args->value[2], args->value[3],
               args->value[4]);
}
static void st_C_interface_ssd_mating(void *a0, void *a1, void *a2, void *a3,
                                      void *a4)
{
    st_C_arguments args = {{a0, a1, a2, a3, a4}};
    st_work_native_interface(st_C_body_ssd_mating, &args);
}
static void st_C_body_st_get_info_level(void *payload)
{
    st_C_arguments *args = payload;
    st_get_info_level(args->value[0]);
}
static void st_C_interface_st_get_info_level(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_st_get_info_level, &args);
}
static void st_C_body_st_get_simpop2_GV(void *payload)
{
    st_C_arguments *args = payload;
    st_get_simpop2_GV(args->value[0], args->value[1], args->value[2],
                      args->value[3]);
}
static void st_C_interface_st_get_simpop2_GV(void *a0, void *a1, void *a2,
                                             void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_st_get_simpop2_GV, &args);
}
static void st_C_body_st_get_simpop_GV(void *payload)
{
    st_C_arguments *args = payload;
    st_get_simpop_GV(args->value[0], args->value[1], args->value[2],
                     args->value[3]);
}
static void st_C_interface_st_get_simpop_GV(void *a0, void *a1, void *a2,
                                            void *a3)
{
    st_C_arguments args = {{a0, a1, a2, a3}};
    st_work_native_interface(st_C_body_st_get_simpop_GV, &args);
}
static void st_C_body_st_set_info_level(void *payload)
{
    st_C_arguments *args = payload;
    st_set_info_level(args->value[0]);
}
static void st_C_interface_st_set_info_level(void *a0)
{
    st_C_arguments args = {{a0}};
    st_work_native_interface(st_C_body_st_set_info_level, &args);
}
static void st_C_body_st_unavailable_st_set_matrix_ops(void *payload)
{
    st_C_arguments *args = payload;
    st_unavailable_st_set_matrix_ops(args->value[0], args->value[1]);
}
static void st_C_interface_st_unavailable_st_set_matrix_ops(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_st_unavailable_st_set_matrix_ops, &args);
}
static void st_C_body_st_set_sim_ef(void *payload)
{
    st_C_arguments *args = payload;
    st_set_sim_ef(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_st_set_sim_ef(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_st_set_sim_ef, &args);
}
static void st_C_body_st_set_sim_gp(void *payload)
{
    st_C_arguments *args = payload;
    st_set_sim_gp(args->value[0], args->value[1]);
}
static void st_C_interface_st_set_sim_gp(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_st_set_sim_gp, &args);
}
static void st_C_body_st_set_sim_mp(void *payload)
{
    st_C_arguments *args = payload;
    st_set_sim_mp(args->value[0], args->value[1]);
}
static void st_C_interface_st_set_sim_mp(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_st_set_sim_mp, &args);
}
static void st_C_body_st_set_sim_pp(void *payload)
{
    st_C_arguments *args = payload;
    st_set_sim_pp(args->value[0], args->value[1], args->value[2]);
}
static void st_C_interface_st_set_sim_pp(void *a0, void *a1, void *a2)
{
    st_C_arguments args = {{a0, a1, a2}};
    st_work_native_interface(st_C_body_st_set_sim_pp, &args);
}
static void st_C_body_st_start_timer(void *payload)
{
    st_C_arguments *args = payload;
    st_start_timer(args->value[0], args->value[1]);
}
static void st_C_interface_st_start_timer(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_st_start_timer, &args);
}
static void st_C_body_st_stop_timer(void *payload)
{
    st_C_arguments *args = payload;
    st_stop_timer(args->value[0], args->value[1]);
}
static void st_C_interface_st_stop_timer(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_st_stop_timer, &args);
}
static void st_C_body_swap_population_name(void *payload)
{
    st_C_arguments *args = payload;
    swap_population_name(args->value[0], args->value[1]);
}
static void st_C_interface_swap_population_name(void *a0, void *a1)
{
    st_C_arguments args = {{a0, a1}};
    st_work_native_interface(st_C_body_swap_population_name, &args);
}

void st_legacy_C_call(const char *name, void **value, int count)
{
    st_C_arguments args;
    if (count < 0 || count > 17)
        error("Invalid interface argument count");
    for (int i = 0; i < 17; i++)
        args.value[i] = i < count ? value[i] : NULL;
    if (!strcmp(name, "allele_freq"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_allele_freq(&args);
        return;
    }
    if (!strcmp(name, "append_population"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_append_population(&args);
        return;
    }
    if (!strcmp(name, "concat_population"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_concat_population(&args);
        return;
    }
    if (!strcmp(name, "copy_population"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_copy_population(&args);
        return;
    }
    if (!strcmp(name, "correct_missing"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_correct_missing(&args);
        return;
    }
    if (!strcmp(name, "cross"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_cross(&args);
        return;
    }
    if (!strcmp(name, "dh"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_dh(&args);
        return;
    }
    if (!strcmp(name, "divide_population"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_divide_population(&args);
        return;
    }
    if (!strcmp(name, "effect_weight_set_all"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_effect_weight_set_all(&args);
        return;
    }
    if (!strcmp(name, "effmap_remove"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_effmap_remove(&args);
        return;
    }
    if (!strcmp(name, "effmap_remove_all"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_effmap_remove_all(&args);
        return;
    }
    if (!strcmp(name, "evaluate_genome"))
    {
        if (count != 6)
            error("Invalid interface argument count");
        st_C_body_evaluate_genome(&args);
        return;
    }
    if (!strcmp(name, "evaluate_genome_region"))
    {
        if (count != 9)
            error("Invalid interface argument count");
        st_C_body_evaluate_genome_region(&args);
        return;
    }
    if (!strcmp(name, "evaluate_ld"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_evaluate_ld(&args);
        return;
    }
    if (!strcmp(name, "evaluate_locus"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_evaluate_locus(&args);
        return;
    }
    if (!strcmp(name, "evaluate_mdp"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_evaluate_mdp(&args);
        return;
    }
    if (!strcmp(name, "evaluate_population"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_evaluate_population(&args);
        return;
    }
    if (!strcmp(name, "g_r_info"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_st_r_info(&args);
        return;
    }
    if (!strcmp(name, "gen_dist"))
    {
        if (count != 17)
            error("Invalid interface argument count");
        st_C_body_gen_dist(&args);
        return;
    }
    if (!strcmp(name, "generate_effect_file"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_generate_effect_file(&args);
        return;
    }
    if (!strcmp(name, "generate_map_file"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_generate_map_file(&args);
        return;
    }
    if (!strcmp(name, "generate_population"))
    {
        if (count != 17)
            error("Invalid interface argument count");
        st_C_body_generate_population(&args);
        return;
    }
    if (!strcmp(name, "genotype_population"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_genotype_population(&args);
        return;
    }
    if (!strcmp(name, "get_allele_number"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_get_allele_number(&args);
        return;
    }
    if (!strcmp(name, "get_allow_zeros"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_get_allow_zeros(&args);
        return;
    }
    if (!strcmp(name, "get_genome_par"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_get_genome_par(&args);
        return;
    }
    if (!strcmp(name, "get_genome_par_b"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_get_genome_par_b(&args);
        return;
    }
    if (!strcmp(name, "get_info_level"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_get_info_level(&args);
        return;
    }
    if (!strcmp(name, "get_map_a"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_get_map_a(&args);
        return;
    }
    if (!strcmp(name, "get_map_b"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_get_map_b(&args);
        return;
    }
    if (!strcmp(name, "get_mdp"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_get_mdp(&args);
        return;
    }
    if (!strcmp(name, "get_population"))
    {
        if (count != 6)
            error("Invalid interface argument count");
        st_C_body_get_population(&args);
        return;
    }
    if (!strcmp(name, "get_population_info"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_get_population_info(&args);
        return;
    }
    if (!strcmp(name, "get_population_size"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_get_population_size(&args);
        return;
    }
    if (!strcmp(name, "get_population_value"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_get_population_value(&args);
        return;
    }
    if (!strcmp(name, "gs_build_V_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_build_V_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_build_Z_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_build_Z_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_build_esvs_01_GV"))
    {
        if (count != 10)
            error("Invalid interface argument count");
        st_C_body_gs_build_esvs_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_calc_ld_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_calc_ld_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_calc_q_01_GV"))
    {
        if (count != 6)
            error("Invalid interface argument count");
        st_C_body_gs_calc_q_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_calc_rf_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_calc_rf_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_check_pvals_01_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_check_pvals_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_chrom_stats_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_chrom_stats_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_copy_marker_data_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_copy_marker_data_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_eval_es_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_cross_eval_es_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_eval_gd_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_cross_eval_gd_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_eval_gd_fct_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_cross_eval_gd_fct_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_eval_ma_01_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_cross_eval_ma_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_eval_mi_01_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_cross_eval_mi_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_eval_mu_01_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_cross_eval_mu_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_eval_va_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_cross_eval_va_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_info_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_cross_info_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_info_gd_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_cross_info_gd_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_cross_validation_02_GV"))
    {
        if (count != 12)
            error("Invalid interface argument count");
        st_C_body_gs_cross_validation_02_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_def_hblocks_01_GV"))
    {
        if (count != 10)
            error("Invalid interface argument count");
        st_C_body_gs_def_hblocks_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_esteff_lsq_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_esteff_lsq_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_esteff_rmla_01_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_gs_esteff_rmla_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_esteff_rmlc_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_esteff_rmlc_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_esteff_rmlr_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_esteff_rmlr_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_esteff_rmlv_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_esteff_rmlv_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_esteff_rr_01_GV"))
    {
        if (count != 9)
            error("Invalid interface argument count");
        st_C_body_gs_esteff_rr_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_estimate_bv_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_estimate_bv_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_get_info_level_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_get_info_level_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_get_map_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_get_map_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_get_num_threads"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_gs_get_num_threads(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_aov_01_GV"))
    {
        if (count != 6)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_aov_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_const_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_const_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_const_02_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_const_02_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_emstep_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_emstep_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_reg_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_reg_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_rmla_01_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_rmla_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_rmla_02_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_rmla_02_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_rmlc_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_rmlc_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_rmlr_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_rmlr_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_rmlv_01_GV"))
    {
        if (count != 6)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_lambda_rmlv_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_rrblup_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_rrblup_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_lambda_var_01_GV"))
    {
        if (count != 6)
            error("Invalid interface argument count");
        st_C_body_gs_lambda_var_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_marker_stats_01_GV"))
    {
        if (count != 9)
            error("Invalid interface argument count");
        st_C_body_gs_marker_stats_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mme_coeff_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_mme_coeff_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mme_invcoeff_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_mme_invcoeff_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mme_restcoeff_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_mme_restcoeff_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mme_restrhs_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_mme_restrhs_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mme_rhs_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_mme_rhs_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mme_solve_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_mme_solve_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mmet_coeff_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_mmet_coeff_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mmet_coeff_03_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_mmet_coeff_03_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mmet_rhs_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_mmet_rhs_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_mmet_solve_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_mmet_solve_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_neg_effall_01_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_neg_effall_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_pos_effall_01_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_pos_effall_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_r_info"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_r_info(&args);
        return;
    }
    if (!strcmp(name, "gs_read_map_01_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_gs_read_map_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_read_marker_data_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_read_marker_data_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_read_marker_data_02_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_read_marker_data_02_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_read_marker_data_03_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_read_marker_data_03_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_read_marker_data_04_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_read_marker_data_04_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_read_performance_data_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_read_performance_data_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_recode_hbc_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_recode_hbc_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_recode_hil_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_recode_hil_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_recode_ref_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_recode_ref_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_recode_ref_02_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_recode_ref_02_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_recode_ref_name_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_recode_ref_name_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_reset"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_reset(&args);
        return;
    }
    if (!strcmp(name, "gs_restrict_marker_data_01_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_restrict_marker_data_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_restrict_marker_data_03_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_gs_restrict_marker_data_03_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_return_effects_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_return_effects_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_return_performance_data_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_return_performance_data_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_return_pvals_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_return_pvals_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_set_all_info_levels"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_gs_set_all_info_levels(&args);
        return;
    }
    if (!strcmp(name, "gs_set_allele_codes_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_set_allele_codes_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_set_effects_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_set_effects_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_set_hblocks_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_set_hblocks_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_set_info_level_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_gs_set_info_level_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_set_num_threads"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_gs_set_num_threads(&args);
        return;
    }
    if (!strcmp(name, "gs_set_openblas_threads"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_set_openblas_threads(&args);
        return;
    }
    if (!strcmp(name, "gs_set_performance_data_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_set_performance_data_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_single_marker_aov_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_single_marker_aov_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_single_marker_reg_01_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_gs_single_marker_reg_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_test_mp"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_gs_test_mp(&args);
        return;
    }
    if (!strcmp(name, "gs_testeff_lsq_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_testeff_lsq_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_testeff_rmlc_01_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_gs_testeff_rmlc_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_testeff_rmlv_01_GV"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_gs_testeff_rmlv_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_vc_rrblup_01_GV"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_gs_vc_rrblup_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_w_aov_01_GV"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_gs_w_aov_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_write_map_01_GV"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_gs_write_map_01_GV(&args);
        return;
    }
    if (!strcmp(name, "gs_write_marker_data_02_GV"))
    {
        if (count != 14)
            error("Invalid interface argument count");
        st_C_body_gs_write_marker_data_02_GV(&args);
        return;
    }
    if (!strcmp(name, "haplotype_calc_d"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_haplotype_calc_d(&args);
        return;
    }
    if (!strcmp(name, "haplotype_free"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_haplotype_free(&args);
        return;
    }
    if (!strcmp(name, "haplotype_init"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_haplotype_init(&args);
        return;
    }
    if (!strcmp(name, "haplotype_next"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_haplotype_next(&args);
        return;
    }
    if (!strcmp(name, "il_eval_library_GV"))
    {
        if (count != 6)
            error("Invalid interface argument count");
        st_C_body_il_eval_library_GV(&args);
        return;
    }
    if (!strcmp(name, "il_eval_lines_GV"))
    {
        if (count != 11)
            error("Invalid interface argument count");
        st_C_body_il_eval_lines_GV(&args);
        return;
    }
    if (!strcmp(name, "il_get_population_size_GV"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_il_get_population_size_GV(&args);
        return;
    }
    if (!strcmp(name, "init_population"))
    {
        if (count != 7)
            error("Invalid interface argument count");
        st_C_body_init_population(&args);
        return;
    }
    if (!strcmp(name, "linkage_map_load"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_linkage_map_load(&args);
        return;
    }
    if (!strcmp(name, "list_eff"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_list_eff(&args);
        return;
    }
    if (!strcmp(name, "list_eff_parameter"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_list_eff_parameter(&args);
        return;
    }
    if (!strcmp(name, "list_populations"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_list_populations(&args);
        return;
    }
    if (!strcmp(name, "list_populations_parameter"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_list_populations_parameter(&args);
        return;
    }
    if (!strcmp(name, "load_effmap"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_load_effmap(&args);
        return;
    }
    if (!strcmp(name, "load_effmap_mod"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_load_effmap_mod(&args);
        return;
    }
    if (!strcmp(name, "optimize_population"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_optimize_population(&args);
        return;
    }
    if (!strcmp(name, "phenotype_population"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_phenotype_population(&args);
        return;
    }
    if (!strcmp(name, "print_version_2"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_print_version_2(&args);
        return;
    }
    if (!strcmp(name, "r_info"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_r_info(&args);
        return;
    }
    if (!strcmp(name, "remove_evaluate_population"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_remove_evaluate_population(&args);
        return;
    }
    if (!strcmp(name, "remove_genotype_population"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_remove_genotype_population(&args);
        return;
    }
    if (!strcmp(name, "remove_map"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_remove_map(&args);
        return;
    }
    if (!strcmp(name, "remove_population"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_remove_population(&args);
        return;
    }
    if (!strcmp(name, "rename_population"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_rename_population(&args);
        return;
    }
    if (!strcmp(name, "reset_mdp"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_reset_mdp(&args);
        return;
    }
    if (!strcmp(name, "resize_population"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_resize_population(&args);
        return;
    }
    if (!strcmp(name, "rng_choose"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_rng_choose(&args);
        return;
    }
    if (!strcmp(name, "rng_info"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_rng_info(&args);
        return;
    }
    if (!strcmp(name, "rng_init"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_rng_init(&args);
        return;
    }
    if (!strcmp(name, "rng_init_seed"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_rng_init_seed(&args);
        return;
    }
    if (!strcmp(name, "rng_list"))
    {
        if (count != 0)
            error("Invalid interface argument count");
        st_C_body_rng_list(&args);
        return;
    }
    if (!strcmp(name, "sample_population"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_sample_population(&args);
        return;
    }
    if (!strcmp(name, "select_all_n_best"))
    {
        if (count != 8)
            error("Invalid interface argument count");
        st_C_body_select_all_n_best(&args);
        return;
    }
    if (!strcmp(name, "set_allow_zeros"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_set_allow_zeros(&args);
        return;
    }
    if (!strcmp(name, "set_co_freq"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_set_co_freq(&args);
        return;
    }
    if (!strcmp(name, "set_eff_weight"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_set_eff_weight(&args);
        return;
    }
    if (!strcmp(name, "set_genome_par"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_set_genome_par(&args);
        return;
    }
    if (!strcmp(name, "set_info_level"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_set_info_level(&args);
        return;
    }
    if (!strcmp(name, "set_mdp"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_set_mdp(&args);
        return;
    }
    if (!strcmp(name, "set_nolociinit"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_set_nolociinit(&args);
        return;
    }
    if (!strcmp(name, "set_population_gvalue"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_set_population_gvalue(&args);
        return;
    }
    if (!strcmp(name, "set_population_info"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_set_population_info(&args);
        return;
    }
    if (!strcmp(name, "sort_population"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_sort_population(&args);
        return;
    }
    if (!strcmp(name, "ssd_mating"))
    {
        if (count != 5)
            error("Invalid interface argument count");
        st_C_body_ssd_mating(&args);
        return;
    }
    if (!strcmp(name, "st_get_info_level"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_st_get_info_level(&args);
        return;
    }
    if (!strcmp(name, "st_get_simpop2_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_st_get_simpop2_GV(&args);
        return;
    }
    if (!strcmp(name, "st_get_simpop_GV"))
    {
        if (count != 4)
            error("Invalid interface argument count");
        st_C_body_st_get_simpop_GV(&args);
        return;
    }
    if (!strcmp(name, "st_r_info"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_st_r_info(&args);
        return;
    }
    if (!strcmp(name, "st_set_info_level"))
    {
        if (count != 1)
            error("Invalid interface argument count");
        st_C_body_st_set_info_level(&args);
        return;
    }
    if (!strcmp(name, "st_set_matrix_ops"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_st_unavailable_st_set_matrix_ops(&args);
        return;
    }
    if (!strcmp(name, "st_set_sim_ef"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_st_set_sim_ef(&args);
        return;
    }
    if (!strcmp(name, "st_set_sim_gp"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_st_set_sim_gp(&args);
        return;
    }
    if (!strcmp(name, "st_set_sim_mp"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_st_set_sim_mp(&args);
        return;
    }
    if (!strcmp(name, "st_set_sim_pp"))
    {
        if (count != 3)
            error("Invalid interface argument count");
        st_C_body_st_set_sim_pp(&args);
        return;
    }
    if (!strcmp(name, "st_start_timer"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_st_start_timer(&args);
        return;
    }
    if (!strcmp(name, "st_stop_timer"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_st_stop_timer(&args);
        return;
    }
    if (!strcmp(name, "swap_population_name"))
    {
        if (count != 2)
            error("Invalid interface argument count");
        st_C_body_swap_population_name(&args);
        return;
    }
    error("Unknown legacy interface");
}

static const R_CMethodDef CEntries[] = {
    {"allele_freq", (DL_FUNC)(void (*)(void))&st_C_interface_allele_freq, 8, NULL},
    {"append_population", (DL_FUNC)(void (*)(void))&st_C_interface_append_population, 2, NULL},
    {"concat_population", (DL_FUNC)(void (*)(void))&st_C_interface_concat_population, 2, NULL},
    {"copy_population", (DL_FUNC)(void (*)(void))&st_C_interface_copy_population, 4, NULL},
    {"correct_missing", (DL_FUNC)(void (*)(void))&st_C_interface_correct_missing, 8, NULL},
    {"cross", (DL_FUNC)(void (*)(void))&st_C_interface_cross, 8, NULL},
    {"dh", (DL_FUNC)(void (*)(void))&st_C_interface_dh, 4, NULL},
    {"divide_population", (DL_FUNC)(void (*)(void))&st_C_interface_divide_population, 3, NULL},
    {"effect_weight_set_all", (DL_FUNC)(void (*)(void))&st_C_interface_effect_weight_set_all,
     1, NULL},
    {"effmap_remove", (DL_FUNC)(void (*)(void))&st_C_interface_effmap_remove, 1, NULL},
    {"effmap_remove_all", (DL_FUNC)(void (*)(void))&st_C_interface_effmap_remove_all, 0, NULL},
    {"evaluate_genome", (DL_FUNC)(void (*)(void))&st_C_interface_evaluate_genome, 6, NULL},
    {"evaluate_genome_region", (DL_FUNC)(void (*)(void))&st_C_interface_evaluate_genome_region,
     9, NULL},
    {"evaluate_ld", (DL_FUNC)(void (*)(void))&st_C_interface_evaluate_ld, 5, NULL},
    {"evaluate_locus", (DL_FUNC)(void (*)(void))&st_C_interface_evaluate_locus, 5, NULL},
    {"evaluate_mdp", (DL_FUNC)(void (*)(void))&st_C_interface_evaluate_mdp, 4, NULL},
    {"evaluate_population", (DL_FUNC)(void (*)(void))&st_C_interface_evaluate_population, 2, NULL},
    {"g_r_info", (DL_FUNC)(void (*)(void))&st_C_interface_st_r_info, 2, NULL},
    {"gen_dist", (DL_FUNC)(void (*)(void))&st_C_interface_gen_dist, 17, NULL},
    {"generate_effect_file", (DL_FUNC)(void (*)(void))&st_C_interface_generate_effect_file, 2, NULL},
    {"generate_map_file", (DL_FUNC)(void (*)(void))&st_C_interface_generate_map_file, 3, NULL},
    {"generate_population", (DL_FUNC)(void (*)(void))&st_C_interface_generate_population, 17, NULL},
    {"genotype_population", (DL_FUNC)(void (*)(void))&st_C_interface_genotype_population, 1, NULL},
    {"get_allele_number", (DL_FUNC)(void (*)(void))&st_C_interface_get_allele_number, 2, NULL},
    {"get_allow_zeros", (DL_FUNC)(void (*)(void))&st_C_interface_get_allow_zeros, 1, NULL},
    {"get_genome_par", (DL_FUNC)(void (*)(void))&st_C_interface_get_genome_par, 2, NULL},
    {"get_genome_par_b", (DL_FUNC)(void (*)(void))&st_C_interface_get_genome_par_b, 1, NULL},
    {"get_info_level", (DL_FUNC)(void (*)(void))&st_C_interface_get_info_level, 1, NULL},
    {"get_map_a", (DL_FUNC)(void (*)(void))&st_C_interface_get_map_a, 4, NULL},
    {"get_map_b", (DL_FUNC)(void (*)(void))&st_C_interface_get_map_b, 4, NULL},
    {"get_mdp", (DL_FUNC)(void (*)(void))&st_C_interface_get_mdp, 1, NULL},
    {"get_population", (DL_FUNC)(void (*)(void))&st_C_interface_get_population, 6, NULL},
    {"get_population_info", (DL_FUNC)(void (*)(void))&st_C_interface_get_population_info, 3, NULL},
    {"get_population_size", (DL_FUNC)(void (*)(void))&st_C_interface_get_population_size, 3, NULL},
    {"get_population_value", (DL_FUNC)(void (*)(void))&st_C_interface_get_population_value, 5, NULL},
    {"gs_build_V_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_build_V_01_GV, 4, NULL},
    {"gs_build_Z_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_build_Z_01_GV, 4, NULL},
    {"gs_build_esvs_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_build_esvs_01_GV, 10, NULL},
    {"gs_calc_ld_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_calc_ld_01_GV, 5, NULL},
    {"gs_calc_q_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_calc_q_01_GV, 6, NULL},
    {"gs_calc_rf_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_calc_rf_01_GV, 5, NULL},
    {"gs_check_pvals_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_check_pvals_01_GV, 2, NULL},
    {"gs_chrom_stats_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_chrom_stats_01_GV, 4, NULL},
    {"gs_copy_marker_data_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_copy_marker_data_01_GV, 3, NULL},
    {"gs_cross_eval_es_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_eval_es_01_GV,
     3, NULL},
    {"gs_cross_eval_gd_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_eval_gd_01_GV,
     3, NULL},
    {"gs_cross_eval_gd_fct_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_eval_gd_fct_01_GV, 4, NULL},
    {"gs_cross_eval_ma_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_eval_ma_01_GV,
     2, NULL},
    {"gs_cross_eval_mi_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_eval_mi_01_GV,
     2, NULL},
    {"gs_cross_eval_mu_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_eval_mu_01_GV,
     2, NULL},
    {"gs_cross_eval_va_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_eval_va_01_GV,
     5, NULL},
    {"gs_cross_info_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_info_01_GV, 5, NULL},
    {"gs_cross_info_gd_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_cross_info_gd_01_GV, 3, NULL},
    {"gs_cross_validation_02_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_cross_validation_02_GV, 12, NULL},
    {"gs_def_hblocks_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_def_hblocks_01_GV, 10, NULL},
    {"gs_esteff_lsq_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_esteff_lsq_01_GV, 4, NULL},
    {"gs_esteff_rmla_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_esteff_rmla_01_GV, 8, NULL},
    {"gs_esteff_rmlc_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_esteff_rmlc_01_GV, 7, NULL},
    {"gs_esteff_rmlr_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_esteff_rmlr_01_GV, 7, NULL},
    {"gs_esteff_rmlv_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_esteff_rmlv_01_GV, 7, NULL},
    {"gs_esteff_rr_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_esteff_rr_01_GV, 9, NULL},
    {"gs_estimate_bv_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_estimate_bv_01_GV, 5, NULL},
    {"gs_get_info_level_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_get_info_level_GV, 2, NULL},
    {"gs_get_map_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_get_map_01_GV, 3, NULL},
    {"gs_get_num_threads", (DL_FUNC)(void (*)(void))&st_C_interface_gs_get_num_threads, 1, NULL},
    {"gs_lambda_aov_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_aov_01_GV, 6, NULL},
    {"gs_lambda_const_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_const_01_GV,
     5, NULL},
    {"gs_lambda_const_02_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_const_02_GV,
     5, NULL},
    {"gs_lambda_emstep_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_emstep_01_GV,
     5, NULL},
    {"gs_lambda_reg_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_reg_01_GV, 5, NULL},
    {"gs_lambda_rmla_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_rmla_01_GV, 8, NULL},
    {"gs_lambda_rmla_02_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_rmla_02_GV, 8, NULL},
    {"gs_lambda_rmlc_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_rmlc_01_GV, 7, NULL},
    {"gs_lambda_rmlr_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_rmlr_01_GV, 7, NULL},
    {"gs_lambda_rmlv_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_lambda_rmlv_01_GV, 6, NULL},
    {"gs_lambda_rrblup_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_rrblup_01_GV,
     7, NULL},
    {"gs_lambda_var_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_lambda_var_01_GV, 6, NULL},
    {"gs_marker_stats_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_marker_stats_01_GV,
     9, NULL},
    {"gs_mme_coeff_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mme_coeff_01_GV, 4, NULL},
    {"gs_mme_invcoeff_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_mme_invcoeff_01_GV, 3, NULL},
    {"gs_mme_restcoeff_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mme_restcoeff_01_GV,
     3, NULL},
    {"gs_mme_restrhs_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mme_restrhs_01_GV, 3, NULL},
    {"gs_mme_rhs_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mme_rhs_01_GV, 4, NULL},
    {"gs_mme_solve_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mme_solve_01_GV, 4, NULL},
    {"gs_mmet_coeff_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mmet_coeff_01_GV, 4, NULL},
    {"gs_mmet_coeff_03_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_mmet_coeff_03_GV, 4, NULL},
    {"gs_mmet_rhs_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mmet_rhs_01_GV, 4, NULL},
    {"gs_mmet_solve_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_mmet_solve_01_GV, 4, NULL},
    {"gs_neg_effall_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_neg_effall_01_GV, 2, NULL},
    {"gs_pos_effall_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_pos_effall_01_GV, 2, NULL},
    {"gs_r_info", (DL_FUNC)(void (*)(void))&st_C_interface_gs_r_info, 2, NULL},
    {"gs_read_map_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_read_map_01_GV, 8, NULL},
    {"gs_read_marker_data_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_read_marker_data_01_GV, 7, NULL},
    {"gs_read_marker_data_02_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_read_marker_data_02_GV, 7, NULL},
    {"gs_read_marker_data_03_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_read_marker_data_03_GV, 7, NULL},
    {"gs_read_marker_data_04_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_read_marker_data_04_GV, 7, NULL},
    {"gs_read_performance_data_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_read_performance_data_01_GV, 5, NULL},
    {"gs_recode_hbc_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_recode_hbc_01_GV, 3, NULL},
    {"gs_recode_hil_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_recode_hil_01_GV, 4, NULL},
    {"gs_recode_ref_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_recode_ref_01_GV, 4, NULL},
    {"gs_recode_ref_02_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_recode_ref_02_GV, 5, NULL},
    {"gs_recode_ref_name_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_recode_ref_name_01_GV, 4, NULL},
    {"gs_reset", (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_reset, 0, NULL},
    {"gs_restrict_marker_data_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_restrict_marker_data_01_GV, 8, NULL},
    {"gs_restrict_marker_data_03_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_restrict_marker_data_03_GV, 8, NULL},
    {"gs_return_effects_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_return_effects_01_GV, 3, NULL},
    {"gs_return_performance_data_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_return_performance_data_01_GV, 4, NULL},
    {"gs_return_pvals_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_return_pvals_01_GV,
     3, NULL},
    {"gs_set_all_info_levels", (DL_FUNC)(void (*)(void))&st_C_interface_gs_set_all_info_levels,
     1, NULL},
    {"gs_set_allele_codes_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_set_allele_codes_01_GV, 7, NULL},
    {"gs_set_effects_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_set_effects_01_GV, 3, NULL},
    {"gs_set_hblocks_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_set_hblocks_01_GV, 7, NULL},
    {"gs_set_info_level_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_set_info_level_GV, 2, NULL},
    {"gs_set_num_threads", (DL_FUNC)(void (*)(void))&st_C_interface_gs_set_num_threads, 1, NULL},
    {"gs_set_openblas_threads",
     (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_set_openblas_threads, 1, NULL},
    {"gs_set_performance_data_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_set_performance_data_01_GV, 4, NULL},
    {"gs_single_marker_aov_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_single_marker_aov_01_GV, 5, NULL},
    {"gs_single_marker_reg_01_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_single_marker_reg_01_GV, 4, NULL},
    {"gs_test_mp", (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_gs_test_mp, 1, NULL},
    {"gs_testeff_lsq_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_testeff_lsq_01_GV, 5, NULL},
    {"gs_testeff_rmlc_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_testeff_rmlc_01_GV,
     8, NULL},
    {"gs_testeff_rmlv_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_testeff_rmlv_01_GV,
     8, NULL},
    {"gs_vc_rrblup_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_vc_rrblup_01_GV, 7, NULL},
    {"gs_w_aov_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_w_aov_01_GV, 5, NULL},
    {"gs_write_map_01_GV", (DL_FUNC)(void (*)(void))&st_C_interface_gs_write_map_01_GV, 3, NULL},
    {"gs_write_marker_data_02_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_gs_write_marker_data_02_GV, 14, NULL},
    {"haplotype_calc_d", (DL_FUNC)(void (*)(void))&st_C_interface_haplotype_calc_d, 3, NULL},
    {"haplotype_free", (DL_FUNC)(void (*)(void))&st_C_interface_haplotype_free, 0, NULL},
    {"haplotype_init", (DL_FUNC)(void (*)(void))&st_C_interface_haplotype_init, 4, NULL},
    {"haplotype_next", (DL_FUNC)(void (*)(void))&st_C_interface_haplotype_next, 7, NULL},
    {"il_eval_library_GV", (DL_FUNC)(void (*)(void))&st_C_interface_il_eval_library_GV, 6, NULL},
    {"il_eval_lines_GV", (DL_FUNC)(void (*)(void))&st_C_interface_il_eval_lines_GV, 11, NULL},
    {"il_get_population_size_GV",
     (DL_FUNC)(void (*)(void))&st_C_interface_il_get_population_size_GV, 2, NULL},
    {"init_population", (DL_FUNC)(void (*)(void))&st_C_interface_init_population, 7, NULL},
    {"linkage_map_load", (DL_FUNC)(void (*)(void))&st_C_interface_linkage_map_load, 5, NULL},
    {"list_eff", (DL_FUNC)(void (*)(void))&st_C_interface_list_eff, 2, NULL},
    {"list_eff_parameter", (DL_FUNC)(void (*)(void))&st_C_interface_list_eff_parameter, 2, NULL},
    {"list_populations", (DL_FUNC)(void (*)(void))&st_C_interface_list_populations, 2, NULL},
    {"list_populations_parameter",
     (DL_FUNC)(void (*)(void))&st_C_interface_list_populations_parameter, 2, NULL},
    {"load_effmap", (DL_FUNC)(void (*)(void))&st_C_interface_load_effmap, 2, NULL},
    {"load_effmap_mod", (DL_FUNC)(void (*)(void))&st_C_interface_load_effmap_mod, 3, NULL},
    {"optimize_population", (DL_FUNC)(void (*)(void))&st_C_interface_optimize_population, 1, NULL},
    {"phenotype_population", (DL_FUNC)(void (*)(void))&st_C_interface_phenotype_population, 4, NULL},
    {"print_version_2", (DL_FUNC)(void (*)(void))&st_C_interface_print_version_2, 1, NULL},
    {"r_info", (DL_FUNC)(void (*)(void))&st_C_interface_r_info, 2, NULL},
    {"remove_evaluate_population",
     (DL_FUNC)(void (*)(void))&st_C_interface_remove_evaluate_population, 1, NULL},
    {"remove_genotype_population",
     (DL_FUNC)(void (*)(void))&st_C_interface_remove_genotype_population, 1, NULL},
    {"remove_map", (DL_FUNC)(void (*)(void))&st_C_interface_remove_map, 0, NULL},
    {"remove_population", (DL_FUNC)(void (*)(void))&st_C_interface_remove_population, 1, NULL},
    {"rename_population", (DL_FUNC)(void (*)(void))&st_C_interface_rename_population, 2, NULL},
    {"reset_mdp", (DL_FUNC)(void (*)(void))&st_C_interface_reset_mdp, 0, NULL},
    {"resize_population", (DL_FUNC)(void (*)(void))&st_C_interface_resize_population, 2, NULL},
    {"rng_choose", (DL_FUNC)(void (*)(void))&st_C_interface_rng_choose, 1, NULL},
    {"rng_info", (DL_FUNC)(void (*)(void))&st_C_interface_rng_info, 0, NULL},
    {"rng_init", (DL_FUNC)(void (*)(void))&st_C_interface_rng_init, 0, NULL},
    {"rng_init_seed", (DL_FUNC)(void (*)(void))&st_C_interface_rng_init_seed, 3, NULL},
    {"rng_list", (DL_FUNC)(void (*)(void))&st_C_interface_rng_list, 0, NULL},
    {"sample_population", (DL_FUNC)(void (*)(void))&st_C_interface_sample_population, 4, NULL},
    {"select_all_n_best", (DL_FUNC)(void (*)(void))&st_C_interface_select_all_n_best, 8, NULL},
    {"set_allow_zeros", (DL_FUNC)(void (*)(void))&st_C_interface_set_allow_zeros, 1, NULL},
    {"set_co_freq", (DL_FUNC)(void (*)(void))&st_C_interface_set_co_freq, 1, NULL},
    {"set_eff_weight", (DL_FUNC)(void (*)(void))&st_C_interface_set_eff_weight, 2, NULL},
    {"set_genome_par", (DL_FUNC)(void (*)(void))&st_C_interface_set_genome_par, 3, NULL},
    {"set_info_level", (DL_FUNC)(void (*)(void))&st_C_interface_set_info_level, 1, NULL},
    {"set_mdp", (DL_FUNC)(void (*)(void))&st_C_interface_set_mdp, 1, NULL},
    {"set_nolociinit", (DL_FUNC)(void (*)(void))&st_C_interface_set_nolociinit, 1, NULL},
    {"set_population_gvalue", (DL_FUNC)(void (*)(void))&st_C_interface_set_population_gvalue,
     4, NULL},
    {"set_population_info", (DL_FUNC)(void (*)(void))&st_C_interface_set_population_info, 3, NULL},
    {"sort_population", (DL_FUNC)(void (*)(void))&st_C_interface_sort_population, 5, NULL},
    {"ssd_mating", (DL_FUNC)(void (*)(void))&st_C_interface_ssd_mating, 5, NULL},
    {"st_get_info_level", (DL_FUNC)(void (*)(void))&st_C_interface_st_get_info_level, 1, NULL},
    {"st_get_simpop2_GV", (DL_FUNC)(void (*)(void))&st_C_interface_st_get_simpop2_GV, 4, NULL},
    {"st_get_simpop_GV", (DL_FUNC)(void (*)(void))&st_C_interface_st_get_simpop_GV, 4, NULL},
    {"st_r_info", (DL_FUNC)(void (*)(void))&st_C_interface_st_r_info, 2, NULL},
    {"st_set_info_level", (DL_FUNC)(void (*)(void))&st_C_interface_st_set_info_level, 1, NULL},
    {"st_set_matrix_ops",
     (DL_FUNC)(void (*)(void))&st_C_interface_st_unavailable_st_set_matrix_ops, 2, NULL},
    {"st_set_sim_ef", (DL_FUNC)(void (*)(void))&st_C_interface_st_set_sim_ef, 3, NULL},
    {"st_set_sim_gp", (DL_FUNC)(void (*)(void))&st_C_interface_st_set_sim_gp, 2, NULL},
    {"st_set_sim_mp", (DL_FUNC)(void (*)(void))&st_C_interface_st_set_sim_mp, 2, NULL},
    {"st_set_sim_pp", (DL_FUNC)(void (*)(void))&st_C_interface_st_set_sim_pp, 3, NULL},
    {"st_start_timer", (DL_FUNC)(void (*)(void))&st_C_interface_st_start_timer, 2, NULL},
    {"st_stop_timer", (DL_FUNC)(void (*)(void))&st_C_interface_st_stop_timer, 2, NULL},
    {"swap_population_name", (DL_FUNC)(void (*)(void))&st_C_interface_swap_population_name, 2, NULL},
    {NULL, NULL, 0, NULL}};

extern SEXP st_allele_repack(SEXP, SEXP);
#ifdef ST_ALLELE_TESTING
extern SEXP st_allele_test(SEXP, SEXP, SEXP);
#endif
extern SEXP st_data_set_list(void);
extern SEXP st_hybrid_data_set(SEXP);
extern SEXP st_compact_data_set(SEXP);
extern SEXP st_data_set_remove(SEXP);
extern SEXP st_data_set_remove_all(void);
extern SEXP st_data_combine(SEXP, SEXP);
extern SEXP st_keep_intersection_markers(SEXP, SEXP);
extern SEXP st_enrichment(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP st_detect_core(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP st_legacy_sim(SEXP);
extern SEXP st_legacy_marker(SEXP);

typedef struct
{
    SEXP value[26];
} st_interface_arguments;
static SEXP st_interface_body_st_lmm_frontend(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lmm_frontend(args->value[0]);
}
static SEXP st_interface_st_lmm_frontend(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_lmm_frontend, &args);
}
static SEXP st_interface_body_st_lmm_capture_body(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lmm_capture_body(args->value[0]);
}
static SEXP st_interface_st_lmm_capture_body(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_lmm_capture_body, &args);
}
static SEXP st_interface_body_st_lmm_capture_warning(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lmm_capture_warning(args->value[0], args->value[1]);
}
static SEXP st_interface_st_lmm_capture_warning(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_lmm_capture_warning, &args);
}
static SEXP st_interface_body_st_population_marker_map(void *payload)
{
    st_interface_arguments *args = payload;
    return st_population_marker_map(args->value[0]);
}
static SEXP st_interface_st_population_marker_map(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_population_marker_map, &args);
}
static SEXP st_interface_body_st_extend_marker_map(void *payload)
{
    st_interface_arguments *args = payload;
    return st_extend_marker_map(args->value[0]);
}
static SEXP st_interface_st_extend_marker_map(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_extend_marker_map, &args);
}
static SEXP st_interface_body_st_population_cross(void *payload)
{
    st_interface_arguments *args = payload;
    return st_population_cross(args->value[0], args->value[1], args->value[2],
                               args->value[3]);
}
static SEXP st_interface_st_population_cross(SEXP a0, SEXP a1, SEXP a2, SEXP a3)
{
    st_interface_arguments args = {{a0, a1, a2, a3}};
    return st_work_interface(st_interface_body_st_population_cross, &args);
}
static SEXP st_interface_body_st_marker_phenotypes(void *payload)
{
    st_interface_arguments *args = payload;
    return st_marker_phenotypes(args->value[0]);
}
static SEXP st_interface_st_marker_phenotypes(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_marker_phenotypes, &args);
}
static SEXP st_interface_body_st_marker_predict(void *payload)
{
    st_interface_arguments *args = payload;
    return st_marker_predict(args->value[0]);
}
static SEXP st_interface_st_marker_predict(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_marker_predict, &args);
}
static SEXP st_interface_body_st_marker_results(void *payload)
{
    st_interface_arguments *args = payload;
    return st_marker_results(args->value[0]);
}
static SEXP st_interface_st_marker_results(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_marker_results, &args);
}
static SEXP st_interface_body_st_simulation_identity(void *payload)
{
    st_interface_arguments *args = payload;
    return st_simulation_identity(args->value[0], args->value[1],
                                  args->value[2]);
}
static SEXP st_interface_st_simulation_identity(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_st_simulation_identity, &args);
}
static SEXP st_interface_body_st_population_pedigree(void *payload)
{
    st_interface_arguments *args = payload;
    return st_population_pedigree(args->value[0], args->value[1]);
}
static SEXP st_interface_st_population_pedigree(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_population_pedigree, &args);
}
static SEXP st_interface_body_st_allele_repack(void *payload)
{
    st_interface_arguments *args = payload;
    return st_allele_repack(args->value[0], args->value[1]);
}
static SEXP st_interface_st_allele_repack(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_allele_repack, &args);
}
#ifdef ST_ALLELE_TESTING
static SEXP st_interface_body_st_allele_test(void *payload)
{
    st_interface_arguments *args = payload;
    return st_allele_test(args->value[0], args->value[1], args->value[2]);
}
static SEXP st_interface_st_allele_test(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_st_allele_test, &args);
}
#endif
static SEXP st_interface_body_st_phase(void *payload)
{
    st_interface_arguments *args = payload;
    return st_phase(args->value[0], args->value[1], args->value[2],
                    args->value[3], args->value[4]);
}
static SEXP st_interface_st_phase(SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4)
{
    st_interface_arguments args = {{a0, a1, a2, a3, a4}};
    return st_work_interface(st_interface_body_st_phase, &args);
}
static SEXP st_interface_body_st_data_set_list(void *payload)
{
    (void)payload;
    return st_data_set_list();
}
static SEXP st_interface_st_data_set_list(void)
{
    st_interface_arguments args = {{R_NilValue}};
    return st_work_interface(st_interface_body_st_data_set_list, &args);
}
static SEXP st_interface_body_st_hybrid_data_set(void *payload)
{
    st_interface_arguments *args = payload;
    return st_hybrid_data_set(args->value[0]);
}
static SEXP st_interface_st_hybrid_data_set(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_hybrid_data_set, &args);
}
static SEXP st_interface_body_st_compact_data_set(void *payload)
{
    st_interface_arguments *args = payload;
    return st_compact_data_set(args->value[0]);
}
static SEXP st_interface_st_compact_data_set(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_compact_data_set, &args);
}
static SEXP st_interface_body_st_data_set_remove(void *payload)
{
    st_interface_arguments *args = payload;
    return st_data_set_remove(args->value[0]);
}
static SEXP st_interface_st_data_set_remove(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_data_set_remove, &args);
}
static SEXP st_interface_body_st_data_set_remove_all(void *payload)
{
    (void)payload;
    return st_data_set_remove_all();
}
static SEXP st_interface_st_data_set_remove_all(void)
{
    st_interface_arguments args = {{R_NilValue}};
    return st_work_interface(st_interface_body_st_data_set_remove_all, &args);
}
static SEXP st_interface_body_st_data_combine(void *payload)
{
    st_interface_arguments *args = payload;
    return st_data_combine(args->value[0], args->value[1]);
}
static SEXP st_interface_st_data_combine(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_data_combine, &args);
}
static SEXP st_interface_body_st_keep_intersection_markers(void *payload)
{
    st_interface_arguments *args = payload;
    return st_keep_intersection_markers(args->value[0], args->value[1]);
}
static SEXP st_interface_st_keep_intersection_markers(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_keep_intersection_markers,
                             &args);
}
static SEXP st_interface_body_st_enrichment(void *payload)
{
    st_interface_arguments *args = payload;
    return st_enrichment(args->value[0], args->value[1], args->value[2],
                         args->value[3], args->value[4], args->value[5],
                         args->value[6]);
}
static SEXP st_interface_st_enrichment(SEXP a0, SEXP a1, SEXP a2, SEXP a3,
                                       SEXP a4, SEXP a5, SEXP a6)
{
    st_interface_arguments args = {{a0, a1, a2, a3, a4, a5, a6}};
    return st_work_interface(st_interface_body_st_enrichment, &args);
}
static SEXP st_interface_body_st_detect_core(void *payload)
{
    st_interface_arguments *args = payload;
    return st_detect_core(args->value[0], args->value[1], args->value[2],
                          args->value[3], args->value[4], args->value[5]);
}
static SEXP st_interface_st_detect_core(SEXP a0, SEXP a1, SEXP a2, SEXP a3,
                                        SEXP a4, SEXP a5)
{
    st_interface_arguments args = {{a0, a1, a2, a3, a4, a5}};
    return st_work_interface(st_interface_body_st_detect_core, &args);
}
static SEXP st_interface_body_allele_evaluate(void *payload)
{
    st_interface_arguments *args = payload;
    return allele_evaluate(args->value[0], args->value[1], args->value[2],
                           args->value[3]);
}
static SEXP st_interface_allele_evaluate(SEXP a0, SEXP a1, SEXP a2, SEXP a3)
{
    st_interface_arguments args = {{a0, a1, a2, a3}};
    return st_work_interface(st_interface_body_allele_evaluate, &args);
}
static SEXP st_interface_body_define_effects_df(void *payload)
{
    st_interface_arguments *args = payload;
    return define_effects_df(args->value[0], args->value[1]);
}
static SEXP st_interface_define_effects_df(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_define_effects_df, &args);
}
static SEXP st_interface_body_define_map_df(void *payload)
{
    st_interface_arguments *args = payload;
    return define_map_df(args->value[0]);
}
static SEXP st_interface_define_map_df(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_define_map_df, &args);
}
static SEXP st_interface_body_evaluate_haplotype(void *payload)
{
    st_interface_arguments *args = payload;
    return evaluate_haplotype(args->value[0], args->value[1], args->value[2],
                              args->value[3]);
}
static SEXP st_interface_evaluate_haplotype(SEXP a0, SEXP a1, SEXP a2, SEXP a3)
{
    st_interface_arguments args = {{a0, a1, a2, a3}};
    return st_work_interface(st_interface_body_evaluate_haplotype, &args);
}
static SEXP st_interface_body_genotype_evaluate(void *payload)
{
    st_interface_arguments *args = payload;
    return genotype_evaluate(args->value[0], args->value[1], args->value[2]);
}
static SEXP st_interface_genotype_evaluate(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_genotype_evaluate, &args);
}
static SEXP st_interface_body_gs_read_marker_data_DF_list(void *payload)
{
    st_interface_arguments *args = payload;
    return gs_read_marker_data_DF_list(args->value[0], args->value[1]);
}
static SEXP st_interface_gs_read_marker_data_DF_list(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_gs_read_marker_data_DF_list,
                             &args);
}
static SEXP st_interface_body_gs_read_marker_data_DF_matrix(void *payload)
{
    st_interface_arguments *args = payload;
    return gs_read_marker_data_DF_matrix(args->value[0], args->value[1]);
}
static SEXP st_interface_gs_read_marker_data_DF_matrix(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_gs_read_marker_data_DF_matrix,
                             &args);
}
static SEXP st_interface_body_haplotype_return_genotype(void *payload)
{
    (void)payload;
    return haplotype_return_genotype();
}
static SEXP st_interface_haplotype_return_genotype(void)
{
    st_interface_arguments args = {{R_NilValue}};
    return st_work_interface(st_interface_body_haplotype_return_genotype,
                             &args);
}
static SEXP st_interface_body_haplotype_return_table(void *payload)
{
    (void)payload;
    return haplotype_return_table();
}
static SEXP st_interface_haplotype_return_table(void)
{
    st_interface_arguments args = {{R_NilValue}};
    return st_work_interface(st_interface_body_haplotype_return_table, &args);
}
static SEXP st_interface_body_return_population(void *payload)
{
    st_interface_arguments *args = payload;
    return return_population(args->value[0], args->value[1]);
}
static SEXP st_interface_return_population(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_return_population, &args);
}
static SEXP st_interface_body_splitdot(void *payload)
{
    st_interface_arguments *args = payload;
    return splitdot(args->value[0]);
}
static SEXP st_interface_splitdot(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_splitdot, &args);
}
static SEXP st_interface_body_splitdt(void *payload)
{
    st_interface_arguments *args = payload;
    return splitdt(args->value[0]);
}
static SEXP st_interface_splitdt(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_splitdt, &args);
}
static SEXP st_interface_body_st_destroy_phase(void *payload)
{
    st_interface_arguments *args = payload;
    return st_destroy_phase(args->value[0]);
}
static SEXP st_interface_st_destroy_phase(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_destroy_phase, &args);
}
static SEXP st_interface_body_st_read_map_df(void *payload)
{
    st_interface_arguments *args = payload;
    return st_read_map_df(args->value[0], args->value[1], args->value[2]);
}
static SEXP st_interface_st_read_map_df(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_st_read_map_df, &args);
}
static SEXP st_interface_body_st_get_V(void *payload)
{
    st_interface_arguments *args = payload;
    return st_get_V(args->value[0]);
}
static SEXP st_interface_st_get_V(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_get_V, &args);
}
static SEXP st_interface_body_gs_make_A(void *payload)
{
    st_interface_arguments *args = payload;
    return gs_make_A(args->value[0], args->value[1]);
}
static SEXP st_interface_gs_make_A(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_gs_make_A, &args);
}
static SEXP st_interface_body_st_get_Z(void *payload)
{
    st_interface_arguments *args = payload;
    return st_get_Z(args->value[0]);
}
static SEXP st_interface_st_get_Z(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_get_Z, &args);
}
static SEXP st_interface_body_st_get_Z_named(void *payload)
{
    st_interface_arguments *args = payload;
    return st_get_Z_named(args->value[0], args->value[1]);
}
static SEXP st_interface_st_get_Z_named(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_get_Z_named, &args);
}
static SEXP st_interface_body_st_get_im(void *payload)
{
    st_interface_arguments *args = payload;
    return st_get_im(args->value[0], args->value[1], args->value[2]);
}
static SEXP st_interface_st_get_im(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_st_get_im, &args);
}
static SEXP st_interface_body_st_get_y(void *payload)
{
    st_interface_arguments *args = payload;
    return st_get_y(args->value[0]);
}
static SEXP st_interface_st_get_y(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_get_y, &args);
}
static SEXP st_interface_body_st_get_y_named(void *payload)
{
    st_interface_arguments *args = payload;
    return st_get_y_named(args->value[0], args->value[1]);
}
static SEXP st_interface_st_get_y_named(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_get_y_named, &args);
}
static SEXP st_interface_body_st_has_openmp(void *payload)
{
    (void)payload;
    return st_has_openmp();
}
static SEXP st_interface_st_has_openmp(void)
{
    st_interface_arguments args = {{R_NilValue}};
    return st_work_interface(st_interface_body_st_has_openmp, &args);
}
static SEXP st_interface_body_st_thread_control(void *payload)
{
    st_interface_arguments *args = payload;
    return st_thread_control(args->value[0]);
}
static SEXP st_interface_st_thread_control(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_thread_control, &args);
}
static SEXP st_interface_body_st_load_performance_data(void *payload)
{
    st_interface_arguments *args = payload;
    return st_load_performance_data(args->value[0], args->value[1],
                                    args->value[2]);
}
static SEXP st_interface_st_load_performance_data(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_st_load_performance_data, &args);
}
static SEXP st_interface_body_st_load_vcf_data(void *payload)
{
    st_interface_arguments *args = payload;
    return st_load_vcf_data(args->value[0], args->value[1], args->value[2]);
}
static SEXP st_interface_st_load_vcf_data(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_st_load_vcf_data, &args);
}
static SEXP st_interface_body_st_markerdata_to_STvcf(void *payload)
{
    st_interface_arguments *args = payload;
    return st_markerdata_to_STvcf(args->value[0]);
}
static SEXP st_interface_st_markerdata_to_STvcf(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_markerdata_to_STvcf, &args);
}
static SEXP st_interface_body_st_mxd_mv(void *payload)
{
    st_interface_arguments *args = payload;
    return st_mxd_mv(args->value[0]);
}
static SEXP st_interface_st_mxd_mv(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_mxd_mv, &args);
}
static SEXP st_interface_body_st_lmm(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lmm(args->value[0]);
}
static SEXP st_interface_st_lmm(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_lmm, &args);
}
static SEXP st_interface_body_st_lsm(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lsm(args->value[0]);
}
static SEXP st_interface_st_lsm(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_lsm, &args);
}
static SEXP st_interface_body_st_predict(void *payload)
{
    st_interface_arguments *args = payload;
    return st_predict(args->value[0]);
}
static SEXP st_interface_st_predict(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_predict, &args);
}
static SEXP st_interface_body_st_predict_method(void *payload)
{
    st_interface_arguments *args = payload;
    return st_predict_method(args->value[0]);
}
static SEXP st_interface_st_predict_method(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_predict_method, &args);
}
static SEXP st_interface_body_st_lsm_method(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lsm_method(args->value[0]);
}
static SEXP st_interface_st_lsm_method(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_lsm_method, &args);
}
static SEXP st_interface_body_st_lmm_method(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lmm_method(args->value[0]);
}
static SEXP st_interface_st_lmm_method(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_lmm_method, &args);
}
static SEXP st_interface_body_st_mxd_mv_linear(void *payload)
{
    st_interface_arguments *args = payload;
    return st_mxd_mv_linear(args->value[0]);
}
static SEXP st_interface_st_mxd_mv_linear(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_mxd_mv_linear, &args);
}
static SEXP st_interface_body_st_mxd_mv_interface(void *payload)
{
    st_interface_arguments *args = payload;
    return st_mxd_mv_interface(args->value[0]);
}
static SEXP st_interface_st_mxd_mv_interface(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_mxd_mv_interface, &args);
}
static SEXP st_interface_body_st_mxd_mv_block(void *payload)
{
    st_interface_arguments *args = payload;
    return st_mxd_mv_block(args->value[0]);
}
static SEXP st_interface_st_mxd_mv_block(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_mxd_mv_block, &args);
}
static SEXP st_interface_body_st_lmm_compile(void *payload)
{
    st_interface_arguments *args = payload;
    return st_lmm_compile(args->value[0]);
}
static SEXP st_interface_st_lmm_compile(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_lmm_compile, &args);
}
static SEXP st_interface_body_st_kernel_structure(void *payload)
{
    st_interface_arguments *args = payload;
    return st_kernel_structure(args->value[0]);
}
static SEXP st_interface_st_kernel_structure(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_kernel_structure, &args);
}
static SEXP st_interface_body_st_get_A(void *payload)
{
    st_interface_arguments *args = payload;
    return st_get_A(args->value[0], args->value[1]);
}
static SEXP st_interface_st_get_A(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_get_A, &args);
}
static SEXP st_interface_body_st_read_vcf(void *payload)
{
    st_interface_arguments *args = payload;
    return st_read_vcf(args->value[0], args->value[1], args->value[2]);
}
static SEXP st_interface_st_read_vcf(SEXP a0, SEXP a1, SEXP a2)
{
    st_interface_arguments args = {{a0, a1, a2}};
    return st_work_interface(st_interface_body_st_read_vcf, &args);
}
static SEXP st_interface_body_st_switch_error(void *payload)
{
    st_interface_arguments *args = payload;
    return st_switch_error(args->value[0], args->value[1]);
}
static SEXP st_interface_st_switch_error(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_switch_error, &args);
}
static SEXP st_interface_body_st_write_vcf(void *payload)
{
    st_interface_arguments *args = payload;
    return st_write_vcf(args->value[0], args->value[1]);
}
static SEXP st_interface_st_write_vcf(SEXP a0, SEXP a1)
{
    st_interface_arguments args = {{a0, a1}};
    return st_work_interface(st_interface_body_st_write_vcf, &args);
}

static SEXP st_interface_body_st_legacy_sim(void *payload)
{
    st_interface_arguments *args = payload;
    return st_legacy_sim(args->value[0]);
}
static SEXP st_interface_st_legacy_sim(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_legacy_sim, &args);
}

static SEXP st_interface_body_st_legacy_marker(void *payload)
{
    st_interface_arguments *args = payload;
    return st_legacy_marker(args->value[0]);
}
static SEXP st_interface_st_legacy_marker(SEXP a0)
{
    st_interface_arguments args = {{a0}};
    return st_work_interface(st_interface_body_st_legacy_marker, &args);
}

static const R_CallMethodDef CallEntries[] = {
    {"st_legacy_marker", (DL_FUNC)(void (*)(void))&st_interface_st_legacy_marker, 1},
    {"st_legacy_sim", (DL_FUNC)(void (*)(void))&st_interface_st_legacy_sim, 1},
    {"st_lmm_frontend", (DL_FUNC)(void (*)(void))&st_interface_st_lmm_frontend, 1},
    {"st_lmm_capture_body", (DL_FUNC)(void (*)(void))&st_interface_st_lmm_capture_body, 1},
    {"st_lmm_capture_warning", (DL_FUNC)(void (*)(void))&st_interface_st_lmm_capture_warning,
     2},
    {"st_population_marker_map",
     (DL_FUNC)(void (*)(void))&st_interface_st_population_marker_map, 1},
    {"st_extend_marker_map", (DL_FUNC)(void (*)(void))&st_interface_st_extend_marker_map, 1},
    {"st_population_cross", (DL_FUNC)(void (*)(void))&st_interface_st_population_cross, 4},
    {"st_marker_phenotypes", (DL_FUNC)(void (*)(void))&st_interface_st_marker_phenotypes, 1},
    {"st_marker_predict", (DL_FUNC)(void (*)(void))&st_interface_st_marker_predict, 1},
    {"st_marker_results", (DL_FUNC)(void (*)(void))&st_interface_st_marker_results, 1},
    {"st_simulation_identity", (DL_FUNC)(void (*)(void))&st_interface_st_simulation_identity,
     3},
    {"st_population_pedigree", (DL_FUNC)(void (*)(void))&st_interface_st_population_pedigree,
     2},
    {"st_allele_repack", (DL_FUNC)(void (*)(void))&st_interface_st_allele_repack, 2},
#ifdef ST_ALLELE_TESTING
    {"st_allele_test", (DL_FUNC)(void (*)(void))&st_interface_st_allele_test, 3},
#endif
    {"st_phase", (DL_FUNC)(void (*)(void))&st_interface_st_phase, 5},
    {"st_data_set_list", (DL_FUNC)(void (*)(void))&st_interface_st_data_set_list, 0},
    {"st_hybrid_data_set", (DL_FUNC)(void (*)(void))&st_interface_st_hybrid_data_set, 1},
    {"st_compact_data_set", (DL_FUNC)(void (*)(void))&st_interface_st_compact_data_set, 1},
    {"st_data_set_remove", (DL_FUNC)(void (*)(void))&st_interface_st_data_set_remove, 1},
    {"st_data_set_remove_all", (DL_FUNC)(void (*)(void))&st_interface_st_data_set_remove_all,
     0},
    {"st_data_combine", (DL_FUNC)(void (*)(void))&st_interface_st_data_combine, 2},
    {"st_keep_intersection_markers",
     (DL_FUNC)(void (*)(void))&st_interface_st_keep_intersection_markers, 2},
    {"st_enrichment", (DL_FUNC)(void (*)(void))&st_interface_st_enrichment, 7},
    {"st_detect_core", (DL_FUNC)(void (*)(void))&st_interface_st_detect_core, 6},
    {"allele_evaluate", (DL_FUNC)(void (*)(void))&st_interface_allele_evaluate, 4},
    {"define_effects_df", (DL_FUNC)(void (*)(void))&st_interface_define_effects_df, 2},
    {"define_map_df", (DL_FUNC)(void (*)(void))&st_interface_define_map_df, 1},
    {"evaluate_haplotype", (DL_FUNC)(void (*)(void))&st_interface_evaluate_haplotype, 4},
    {"genotype_evaluate", (DL_FUNC)(void (*)(void))&st_interface_genotype_evaluate, 3},
    {"gs_read_marker_data_DF_list",
     (DL_FUNC)(void (*)(void))&st_interface_gs_read_marker_data_DF_list, 2},
    {"gs_read_marker_data_DF_matrix",
     (DL_FUNC)(void (*)(void))&st_interface_gs_read_marker_data_DF_matrix, 2},
    {"haplotype_return_genotype",
     (DL_FUNC)(void (*)(void))&st_interface_haplotype_return_genotype, 0},
    {"haplotype_return_table", (DL_FUNC)(void (*)(void))&st_interface_haplotype_return_table,
     0},
    {"return_population", (DL_FUNC)(void (*)(void))&st_interface_return_population, 2},
    {"splitdot", (DL_FUNC)(void (*)(void))&st_interface_splitdot, 1},
    {"splitdt", (DL_FUNC)(void (*)(void))&st_interface_splitdt, 1},
    {"st_destroy_phase", (DL_FUNC)(void (*)(void))&st_interface_st_destroy_phase, 1},
    {"st_read_map_df", (DL_FUNC)(void (*)(void))&st_interface_st_read_map_df, 3},
    {"st_get_V", (DL_FUNC)(void (*)(void))&st_interface_st_get_V, 1},
    {"gs_make_A", (DL_FUNC)(void (*)(void))&st_interface_gs_make_A, 2},
    {"st_get_Z", (DL_FUNC)(void (*)(void))&st_interface_st_get_Z, 1},
    {"st_get_Z_named", (DL_FUNC)(void (*)(void))&st_interface_st_get_Z_named, 2},
    {"st_get_im", (DL_FUNC)(void (*)(void))&st_interface_st_get_im, 3},
    {"st_get_y", (DL_FUNC)(void (*)(void))&st_interface_st_get_y, 1},
    {"st_get_y_named", (DL_FUNC)(void (*)(void))&st_interface_st_get_y_named, 2},
    {"st_has_openmp", (DL_FUNC)(void (*)(void))&st_interface_st_has_openmp, 0},
    {"st_thread_control", (DL_FUNC)(void (*)(void))&st_interface_st_thread_control, 1},
    {"st_load_performance_data",
     (DL_FUNC)(void (*)(void))&st_interface_st_load_performance_data, 3},
    {"st_load_vcf_data", (DL_FUNC)(void (*)(void))&st_interface_st_load_vcf_data, 3},
    {"st_markerdata_to_STvcf", (DL_FUNC)(void (*)(void))&st_interface_st_markerdata_to_STvcf,
     1},
    {"st_mxd_mv", (DL_FUNC)(void (*)(void))&st_interface_st_mxd_mv, 1},
    {"st_lmm", (DL_FUNC)(void (*)(void))&st_interface_st_lmm, 1},
    {"st_lsm", (DL_FUNC)(void (*)(void))&st_interface_st_lsm, 1},
    {"st_predict", (DL_FUNC)(void (*)(void))&st_interface_st_predict, 1},
    {"st_predict_method", (DL_FUNC)(void (*)(void))&st_interface_st_predict_method, 1},
    {"st_lsm_method", (DL_FUNC)(void (*)(void))&st_interface_st_lsm_method, 1},
    {"st_lmm_method", (DL_FUNC)(void (*)(void))&st_interface_st_lmm_method, 1},
    {"st_mxd_mv_linear", (DL_FUNC)(void (*)(void))&st_interface_st_mxd_mv_linear, 1},
    {"st_mxd_mv_interface", (DL_FUNC)(void (*)(void))&st_interface_st_mxd_mv_interface, 1},
    {"st_mxd_mv_block", (DL_FUNC)(void (*)(void))&st_interface_st_mxd_mv_block, 1},
    {"st_lmm_compile", (DL_FUNC)(void (*)(void))&st_interface_st_lmm_compile, 1},
    {"st_kernel_structure", (DL_FUNC)(void (*)(void))&st_interface_st_kernel_structure, 1},
    {"st_get_A", (DL_FUNC)(void (*)(void))&st_interface_st_get_A, 2},
    {"st_read_vcf", (DL_FUNC)(void (*)(void))&st_interface_st_read_vcf, 3},
    {"st_switch_error", (DL_FUNC)(void (*)(void))&st_interface_st_switch_error, 2},
    {"st_write_vcf", (DL_FUNC)(void (*)(void))&st_interface_st_write_vcf, 2},
    {NULL, NULL, 0}};

void attribute_visible R_init_SelectionTools(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
