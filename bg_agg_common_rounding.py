import matplotlib

matplotlib.use('Agg')

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import constants
import seaborn as sns
from scipy.stats import zscore, rankdata, wilcoxon, ks_2samp
import multiprocessing

INTRA = 0
INTER = 1


def calc_emp_pval(cur_rv, cur_dist):
    cur_dist = np.array(cur_dist, np.float32)
    emp_pvals = []
    if type(cur_rv) == str:
        hg_pvals = np.array(cur_rv[1:-1].split(", "), dtype=np.float32)
    else:
        hg_pvals = np.array([cur_rv], dtype=np.float32)

    for hg_pval in hg_pvals:
        pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), hg_pval, side='left')
        emp_pval = pos / float(np.size(cur_dist))
        emp_pvals.append(emp_pval)

    return emp_pvals


def get_sample_of_permuted_scores(args):
    mat, aneu_types, chr_interaction, pval_method, output = args
    s_sample = pd.Series()
    for a in np.arange(1, 23):
        for arm_a in ['p', 'q']:
            for b in np.arange(a, 23):
                for arm_b in ['p', 'q']:
                    for aneu_type_a, aneu_type_b in aneu_types:

                        if (arm_a <= arm_b and b == a) or (chr_interaction and b == a) or (
                                not chr_interaction and b != a): continue
                        if "{}{}".format(a, arm_a) not in mat.columns: continue
                        if "{}{}".format(b, arm_b) not in mat.columns: continue

                        s_sample.loc[
                            "{}{}{}-{}{}{}".format(a, arm_a, aneu_type_a, b, arm_b, aneu_type_b)] = pval_method(mat, a,
                                                                                                                arm_a,
                                                                                                                aneu_type_a,
                                                                                                                b,
                                                                                                                arm_b,
                                                                                                                aneu_type_b)

    output.append(s_sample)

    try:
        print("done {} permutations".format(len(output)))
    except:
        pass


def generate_bg_dists(aneu_types, chr_interaction, file_name, n_start, n_end, original, pval_method, p_factor,
                      file_suffix, use_cache):
    try:
        int(file_suffix.split('_')[-1])
        cache_file_format = "dist_{}_{}.tsv".format('_'.join(file_suffix.split('_')[:-1]), file_name)
    except:
        cache_file_format = "dist_{}_{}.tsv".format('_'.join(file_suffix.split('_')[:]), file_name)

    if use_cache and os.path.exists(os.path.join(constants.CACHE_FOLDER, cache_file_format.format(file_name))):
        # print "cache_test: ", os.path.join(constants.CACHE_FOLDER, cache_file_format.format(file_name))
        return pd.read_csv(os.path.join(constants.CACHE_FOLDER, cache_file_format.format(file_name)), sep='\t',
                           index_col=0)

    params = []
    p = multiprocessing.Pool(p_factor)
    output = multiprocessing.Manager().list([])
    perm_file_format = "{}_perm_{}.tsv"
    for f_index in np.arange(n_start, n_end):
        mat = np.load(os.path.join(constants.CACHE_FOLDER, file_name, perm_file_format.format(file_name, f_index)))
        mat = pd.DataFrame(data=mat, columns=original.columns, index=original.index)
        params.append([mat, aneu_types, chr_interaction, pval_method, output])

    p.map(get_sample_of_permuted_scores, params)
    p.close()

    dist = pd.concat([s_sample.to_frame() for s_sample in output], axis=1)
    dist.to_csv(os.path.join(constants.CACHE_FOLDER, cache_file_format.format(file_name)), sep='\t')

    return dist


def plot_dists(original, bg, interactions=None, file_name="", p_factor=10, file_suffix=None, use_cache=1):
    if interactions is None:
        interactions = sorted(list(original.index))

    if len(interactions) > 9:
        print("splitting interactions (total: {})".format(len(interactions)))
        interaction_sets = [interactions[a * 9:min((a + 1) * 9, len(interactions))] for a in
                            np.arange(len(interactions) / 9 + 1)]
    else:
        interaction_sets = [interactions]

    df_emp_pvals = pd.DataFrame()
    params = []
    p = multiprocessing.Pool(p_factor)
    output = multiprocessing.Manager().dict()
    for i_sets, interactions in enumerate(interaction_sets):
        # params.append([bg, file_name, i_sets, interactions, original, file_suffix, output])
        plot_interaction_set([bg, file_name, i_sets, interactions, original, file_suffix, output])

    # p.map(plot_interaction_set, params)
    # p.close()

    print("done!")
    for k, v in dict(output).iteritems():
        df_emp_pvals.loc[k, "pval"] = v

    df_emp_pvals['pval'][df_emp_pvals['pval'] == 0] = 1.0 / bg.shape[1]
    df_emp_pvals['qval'] = fdrcorrection0(df_emp_pvals.loc[:, 'pval'])[1]
    df_emp_pvals['enrichment_score'] = -np.log10(df_emp_pvals.loc[:, 'pval'])
    df_emp_pvals['zscore'] = zscore(df_emp_pvals.loc[:, 'enrichment_score'])
    df_emp_pvals['rank'] = rankdata(df_emp_pvals.loc[:, 'pval'])
    df_emp_pvals.to_csv(
        os.path.join(constants.OUTPUT_FOLDER, "emp_pvals_summary_{}_{}.tsv".format(file_suffix, file_name)), sep='\t')

    return df_emp_pvals


def plot_interaction_set(args):
    bg, file_name, i_sets, interactions, original, file_suffix, output = args

    fig_dim = max(int(np.ceil(np.sqrt(len(interactions)))), 2)
    fig, axs = plt.subplots(fig_dim, fig_dim, figsize=(fig_dim * 8, fig_dim * 8))
    #bg.to_csv("test_bg_values_brca_round.tsv", sep="\t")
    #bg = bg.apply(lambda a: np.round(a, -np.log(a)))#applying rounding
    for i, interaction in enumerate(interactions):
        bg.loc[interaction,:] = bg.loc[interaction,:].apply(lambda a: np.round(a, int(np.ceil(-np.log(a)))))  # applying rounding
        enrichment_score = -np.log10(np.round(original.loc[interaction, "pval"], int(np.ceil(-np.log(original.loc[interaction, "pval"])))))
        #enrichment_score = -np.log10(original.loc[interaction, "pval"])
        cur_bg = -np.log10(bg.loc[interaction, :])
        cur_bg[np.isinf(cur_bg)] = 100  # correct inf values
        p = calc_emp_pval(cur_rv=enrichment_score, cur_dist=cur_bg)[0]
        output[interaction] = p
        #plot_dist(enrichment_score, cur_bg, interaction, p, ax=axs[i / fig_dim][i % fig_dim])
    try:
        os.makedirs(os.path.join(constants.OUTPUT_FOLDER, "figures", file_name))
    except:
        pass

    #plt.savefig(os.path.join(constants.OUTPUT_FOLDER, "figures", file_name,
                             #("null_dist_{}_{}_{}.png").format(file_suffix, file_name, i_sets)))
    try:
        print("done {} interactions".format(len(output)))
    except:
        pass


def plot_dist(enrichment_score, bg, interaction, p, ax=None):
    a = bg
    b = [enrichment_score]
    bins = np.histogram(np.hstack((a, b)), bins=20)[1]
    sns.distplot(a, bins=bins, norm_hist=False, kde=False, ax=ax)
    sns.distplot(b, bins=bins, norm_hist=False, kde=False, ax=ax)
    ax.set_xlabel('enrichment scores')
    ax.set_ylabel('# of events')

    ax.annotate(str('real score'), xy=(enrichment_score, 1), xytext=(enrichment_score, 10), color='red',
                xycoords='data',
                textcoords='offset points',
                bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="->"))
    ax.set_yscale('log')
    ax.set_title("events of {}\nempirical p-value: {} (co-occurring score: {})".format(interaction, "%.2E" % p,
                                                                                       "%.2f" % enrichment_score))


