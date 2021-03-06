from pydnameth import DataType, Task, Method
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np


def get_method_metrics_keys(config):
    metrics = []

    if config.experiment.data in [DataType.betas,
                                  DataType.betas_adj,
                                  DataType.epimutations,
                                  DataType.entropy,
                                  DataType.residuals,
                                  DataType.resid_old,
                                  DataType.cells,
                                  DataType.genes,
                                  DataType.bop]:

        if config.experiment.task == Task.table:

            if config.experiment.method == Method.heteroskedasticity:

                metrics = [
                    'item',
                    'aux',
                    'type',
                    'bp_lm',
                    'bp_lm_pvalue',
                    'bp_lm_pvalue_fdr_bh',
                    'bp_lm_pvalue_bonferroni',
                    'bp_fvalue',
                    'bp_f_pvalue',
                    'bp_f_pvalue_fdr_bh',
                    'bp_f_pvalue_bonferroni',
                    'w_lm',
                    'w_lm_pvalue',
                    'w_lm_pvalue_fdr_bh',
                    'w_lm_pvalue_bonferroni',
                    'w_fvalue',
                    'w_f_pvalue',
                    'w_f_pvalue_fdr_bh',
                    'w_f_pvalue_bonferroni',
                    'gq_fvalue',
                    'gq_f_pvalue',
                    'gq_f_pvalue_fdr_bh',
                    'gq_f_pvalue_bonferroni',
                    'gq_type'
                ]

            elif config.experiment.method == Method.manova:

                metrics = [
                    'item',
                    'aux',
                    'class',
                    'genes',
                ]

                method_params = config.experiment.method_params

                for key, values in method_params.items():
                    for val in values:
                        metrics.append(f'{val}_p_value_wilks')
                        metrics.append(f'{val}_p_value_pillai_bartlett')
                        metrics.append(f'{val}_p_value_lawley_hotelling')
                        metrics.append(f'{val}_p_value_roy')

            elif config.experiment.method == Method.formula:

                metrics = [
                    'item',
                    'aux',
                    'R2',
                    'R2_adj',
                    'mean'
                ]

                method_params = config.experiment.method_params

                exog_dict = {}
                for key, values in method_params.items():
                    if key == 'cells':
                        for val in values:
                            if val in config.cells_dict:
                                exog_dict[val] = config.cells_dict[val]
                            else:
                                raise ValueError(f'Wrong cell type in formula: {val}')
                    if key == 'observables':
                        for val in values:
                            if val in config.observables_categorical_dict:
                                exog_dict[val] = config.observables_categorical_dict[val]
                            else:
                                raise ValueError(f'Wrong observable in formula: {val}')

                exog_keys = []
                for exog_type, exog_data in exog_dict.items():
                    if config.is_observables_categorical.get(exog_type, False):
                        exog_keys.append('C(' + exog_type + ')')
                    else:
                        exog_keys.append(exog_type)
                formula = 'cpg ~ ' + ' + '.join(exog_keys)

                exog_dict['cpg'] = np.random.rand(len(config.attributes_indexes))
                data_df = pd.DataFrame(exog_dict)
                reg_res = smf.ols(formula=formula, data=data_df).fit()
                params = dict(reg_res.params)

                for key in params:
                    metrics.append(key)
                    metrics.append(key + '_std')
                    metrics.append(key + '_p_value')

            elif config.experiment.method == Method.formula_new:

                metrics = [
                    'item',
                    'aux',
                    'R2',
                    'R2_adj',
                    'mean'
                ]

                method_params = config.experiment.method_params
                formula = method_params['formula']

                dict_global = {}
                dict_global.update(config.observables_dict.items())
                if len(config.cells_dict) > 0:
                    dict_global.update(config.cells_dict.items())

                dict_global['cpg'] = np.random.rand(len(config.attributes_indexes))

                data_df = pd.DataFrame(dict_global)
                reg_res = smf.ols(formula=formula, data=data_df).fit()
                params = dict(reg_res.params)

                for key in params:
                    metrics.append(key)
                    metrics.append(key + '_std')
                    metrics.append(key + '_p_value')

            elif config.experiment.method == Method.linreg:

                metrics = [
                    'item',
                    'aux',
                    'mean',
                    'R2',
                    'R2_adj',
                    'f_stat',
                    'prob(f_stat)',
                    'log_likelihood',
                    'AIC',
                    'BIC',
                    'omnibus',
                    'prob(omnibus)',
                    'skew',
                    'kurtosis',
                    'durbin_watson',
                    'jarque_bera',
                    'prob(jarque_bera)',
                    'cond_no',
                    'intercept',
                    'slope',
                    'intercept_std',
                    'slope_std',
                    'intercept_p_value',
                    'slope_p_value',
                    'normality_p_value_shapiro',
                    'normality_p_value_ks_wo_params',
                    'normality_p_value_ks_with_params',
                    'normality_p_value_dagostino'
                ]

            elif config.experiment.method == Method.ancova:

                metrics = [
                    'item',
                    'aux',
                    'R2',
                    'R2_adj',
                    'f_stat',
                    'prob(f_stat)',
                    'intercept',
                    'category',
                    'x',
                    'x:category',
                    'intercept_std',
                    'category_std',
                    'x_std',
                    'x:category_std',
                    'intercept_pval',
                    'category_pval',
                    'x_pval',
                    'x:category_pval',
                    'intercept_pval_fdr_bh',
                    'category_pval_fdr_bh',
                    'x_pval_fdr_bh',
                    'x:category_pval_fdr_bh',
                    'intercept_pval_fdr_bon',
                    'category_pval_fdr_bon',
                    'x_pval_fdr_bon',
                    'x:category_pval_fdr_bon',
                ]

            elif config.experiment.method == Method.oma:

                metrics = [
                    'item',
                    'aux',
                    'lin_lin_corr_coeff',
                    'lin_lin_p_value',
                    'lin_lin_p_value_fdr_bh',
                    'lin_lin_p_value_bonferroni',
                    'lin_log_corr_coeff',
                    'lin_log_p_value',
                    'lin_log_p_value_fdr_bh',
                    'lin_log_p_value_bonferroni',
                    'log_lin_corr_coeff',
                    'log_lin_p_value',
                    'log_lin_p_value_fdr_bh',
                    'log_lin_p_value_bonferroni',
                    'log_log_corr_coeff',
                    'log_log_p_value',
                    'log_log_p_value_fdr_bh',
                    'log_log_p_value_bonferroni',
                ]

            elif config.experiment.method == Method.pbc:

                metrics = [
                    'item',
                    'aux',
                    'pbc_corr_coeff',
                    'pbc_p_value',
                    'pbc_p_value_fdr_bh',
                    'pbc_p_value_bonferroni',
                    'anova_p_value',
                    'anova_p_value_fdr_bh',
                    'anova_p_value_bonferroni',
                    'kw_p_value',
                    'kw_p_value_fdr_bh',
                    'kw_p_value_bonferroni'
                ]

            elif config.experiment.method == Method.variance:

                metrics = [
                    'item',
                    'aux',

                    'best_R2',

                    'increasing_div',
                    'increasing_sub',

                    'increasing_type',

                    'box_b_best_type',
                    'box_b_best_R2',
                    'box_b_lin_lin_R2',
                    'box_b_lin_lin_intercept',
                    'box_b_lin_lin_slope',
                    'box_b_lin_lin_intercept_std',
                    'box_b_lin_lin_slope_std',
                    'box_b_lin_lin_intercept_p_value',
                    'box_b_lin_lin_slope_p_value',
                    'box_b_lin_log_R2',
                    'box_b_lin_log_intercept',
                    'box_b_lin_log_slope',
                    'box_b_lin_log_intercept_std',
                    'box_b_lin_log_slope_std',
                    'box_b_lin_log_intercept_p_value',
                    'box_b_lin_log_slope_p_value',
                    'box_b_log_log_R2',
                    'box_b_log_log_intercept',
                    'box_b_log_log_slope',
                    'box_b_log_log_intercept_std',
                    'box_b_log_log_slope_std',
                    'box_b_log_log_intercept_p_value',
                    'box_b_log_log_slope_p_value',

                    'box_m_best_type',
                    'box_m_best_R2',
                    'box_m_lin_lin_R2',
                    'box_m_lin_lin_intercept',
                    'box_m_lin_lin_slope',
                    'box_m_lin_lin_intercept_std',
                    'box_m_lin_lin_slope_std',
                    'box_m_lin_lin_intercept_p_value',
                    'box_m_lin_lin_slope_p_value',
                    'box_m_lin_log_R2',
                    'box_m_lin_log_intercept',
                    'box_m_lin_log_slope',
                    'box_m_lin_log_intercept_std',
                    'box_m_lin_log_slope_std',
                    'box_m_lin_log_intercept_p_value',
                    'box_m_lin_log_slope_p_value',
                    'box_m_log_log_R2',
                    'box_m_log_log_intercept',
                    'box_m_log_log_slope',
                    'box_m_log_log_intercept_std',
                    'box_m_log_log_slope_std',
                    'box_m_log_log_intercept_p_value',
                    'box_m_log_log_slope_p_value',

                    'box_t_best_type',
                    'box_t_best_R2',
                    'box_t_lin_lin_R2',
                    'box_t_lin_lin_intercept',
                    'box_t_lin_lin_slope',
                    'box_t_lin_lin_intercept_std',
                    'box_t_lin_lin_slope_std',
                    'box_t_lin_lin_intercept_p_value',
                    'box_t_lin_lin_slope_p_value',
                    'box_t_lin_log_R2',
                    'box_t_lin_log_intercept',
                    'box_t_lin_log_slope',
                    'box_t_lin_log_intercept_std',
                    'box_t_lin_log_slope_std',
                    'box_t_lin_log_intercept_p_value',
                    'box_t_lin_log_slope_p_value',
                    'box_t_log_log_R2',
                    'box_t_log_log_intercept',
                    'box_t_log_log_slope',
                    'box_t_log_log_intercept_std',
                    'box_t_log_log_slope_std',
                    'box_t_log_log_intercept_p_value',
                    'box_t_log_log_slope_p_value',
                ]

            elif config.experiment.method == Method.cluster:

                metrics = [
                    'item',
                    'aux',
                    'number_of_clusters',
                    'number_of_noise_points',
                    'percent_of_noise_points',
                ]

            elif config.experiment.method == Method.polygon:

                if config.experiment.method_params['method'] == Method.linreg:

                    metrics = [
                        'item',
                        'aux',
                        'area_intersection',
                        'slope_intersection',
                        'max_abs_slope'
                    ]

                elif config.experiment.method_params['method'] == Method.variance:

                    metrics = [
                        'item',
                        'aux',
                        'area_intersection',
                        'increasing',
                        'increasing_id',
                        'begin_rel',
                        'end_rel'
                    ]

            elif config.experiment.method == Method.special:

                metrics = [
                    'item'
                ]

            elif config.experiment.method == Method.z_test_linreg:

                metrics = [
                    'item',
                    'aux',
                    'z_value',
                    'p_value',
                    'abs_z_value'
                ]

            elif config.experiment.method == Method.aggregator:

                metrics = [
                    'item',
                    'aux'
                ]

        elif config.experiment.task == Task.clock:

            if config.experiment.method == Method.linreg:
                metrics = [
                    'item',
                    'aux',
                    'R2',
                    'r',
                    'evs',
                    'mae',
                    'rmse',
                ]

    metrics = [f'{m}_{config.hash[0:8]}' if m not in ['item', 'aux'] else m for m in metrics]

    return metrics
