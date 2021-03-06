import statsmodels.api as sm
import statsmodels.stats.api as sms
import numpy as np


def process_heteroscedasticity(x, y, metrics_dict, suffix):
    x_with_const = sm.add_constant(x)

    results = sm.OLS(y, x_with_const).fit()

    bp_lm, bp_lm_pvalue, bp_fvalue, bp_f_pvalue = sms.het_breuschpagan(results.resid, results.model.exog)
    w_lm, w_lm_pvalue, w_fvalue, w_f_pvalue = sms.het_white(results.resid, results.model.exog)
    gq_fvalue, gq_f_pvalue, gq_type = sms.het_goldfeldquandt(results.resid, results.model.exog)

    beg_lim, end_lim = np.percentile(x, [33, 67])
    beg_ids = []
    end_ids = []
    for t_id, t in enumerate(x):
        if t < beg_lim:
            beg_ids.append(t_id)
        elif t > end_lim:
            end_ids.append(t_id)

    beg_std = np.std(np.array(y)[np.array(beg_ids)])
    end_std = np.std(np.array(y)[np.array(end_ids)])

    if end_std > beg_std:
        type = 'increasing'
    else:
        type = 'decreasing'

    metrics_dict['type' + suffix].append(type)

    metrics_dict['bp_lm' + suffix].append(bp_lm)
    metrics_dict['bp_lm_pvalue' + suffix].append(bp_lm_pvalue)
    metrics_dict['bp_fvalue' + suffix].append(bp_fvalue)
    metrics_dict['bp_f_pvalue' + suffix].append(bp_f_pvalue)

    metrics_dict['w_lm' + suffix].append(w_lm)
    metrics_dict['w_lm_pvalue' + suffix].append(w_lm_pvalue)
    metrics_dict['w_fvalue' + suffix].append(w_fvalue)
    metrics_dict['w_f_pvalue' + suffix].append(w_f_pvalue)

    metrics_dict['gq_fvalue' + suffix].append(gq_fvalue)
    metrics_dict['gq_f_pvalue' + suffix].append(gq_f_pvalue)
    metrics_dict['gq_type' + suffix].append(gq_type)
