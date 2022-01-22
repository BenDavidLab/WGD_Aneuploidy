import matplotlib

#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
#matplotlib.use('pdf')

import numpy as np
import pandas as pd
import os

from scipy.stats import hypergeom, binom_test, ttest_rel
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import seaborn as sns
from scipy.stats import zscore, rankdata, wilcoxon, ks_2samp
import multiprocessing
import scipy.stats
import constants
import argparse
from statannot import add_stat_annotation
import statistics
import re


def data_exploration(ccle_mat):

    print("Entered")
    wgd_plus = {}
    wgd_minus = {}

    for cell_line in ccle_mat.CCLE_ID:
        #print(cell_line)
        tumor_type = cell_line[cell_line.index("_") + 1:]
        if not ccle_mat[(ccle_mat.CCLE_ID == cell_line) & (ccle_mat.ploidy < 2.5)].empty:
            if tumor_type in wgd_minus.keys():
                wgd_minus[tumor_type] += [cell_line]
            else:
                wgd_minus[tumor_type] = [cell_line]
        elif not ccle_mat[(ccle_mat.CCLE_ID == cell_line) & (ccle_mat.ploidy > 3)].empty:
            if tumor_type in wgd_plus.keys():
                wgd_plus[tumor_type] += [cell_line]
            else:
                wgd_plus[tumor_type] = [cell_line]
        else:
            pass

    all_types = []
    for cell_line in ccle_mat.CCLE_ID:
        tumor_type = cell_line[cell_line.index("_") + 1:]
        all_types.append(tumor_type)

    print(len(set(all_types)))

    count = 0
    types = []

    for k, v in wgd_plus.items():
        #all_types.append(k)
        if len(v) > 5 and len(wgd_minus[k]) > 5:
            count += 1
            print(k, len(v), len(wgd_minus[k]))
            types.append(k)

    print(all_types)
    return types


##===Whole-chromosome aneuploidy fraction analysis===
def plot_boxplot(graph_df, fname):
    """Plots grouped box-plot with t-test to indicate tumors with statistically significant
    difference in their fractions"""
    plt.figure(figsize=(30, 10))

    ax = sns.boxplot(x="Type", y="Fraction", data=graph_df, hue="WGD status", palette = "Set1", showfliers=False)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    ax.set_ylim((0,1.1))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    pairs = [((x,"WGD+"),(x,"WGD-")) for x in types]

    add_stat_annotation(ax, x="Type", y="Fraction", data=graph_df, hue="WGD status",box_pairs=pairs,test='t-test_ind', text_format='star', loc='inside', verbose=0)

    plt.show()
    #plt.savefig("{}.png".format(fname))
    #plt.clf()


def analyze_whole_chr_anp_frac_samplewise(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Fraction", "WGD status", "Type"])
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    for type in types:

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)
        #type_df = input_df[type in input_df.CCLE_ID]
        #wgd_plus_df = type_df[type_df.Genome_doublings > 0].iloc[:,12:51]
        #wgd_minus_df = type_df[type_df.Genome_doublings == 0].iloc[:,12:51]
        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:,4:]
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:,4:]


        wgd_plus_df_non_acro = wgd_plus_df.drop(["13q", "14q", "15q", "21q", "22q"], axis=1)
        wgd_minus_df_non_acro = wgd_minus_df.drop(["13q", "14q", "15q", "21q", "22q"], axis=1)
        #print(wgd_plus_df_non_acro)

        whole_chr_anp_frac_wgd_plus = []
        whole_chr_anp_frac_wgd_minus = []




        for i in range(0, len(wgd_plus_df_non_acro.index)):
            total_wgd_plus = 0
            arm_level = 0
            for j in range(0, 34, 2):

                total_wgd_plus += np.sum(np.logical_and(wgd_plus_df_non_acro.iloc[i,j] == 1, wgd_plus_df_non_acro.iloc[i,j+1] == 1)) \
                + np.sum(np.logical_and(wgd_plus_df_non_acro.iloc[i,j] == -1, wgd_plus_df_non_acro.iloc[i,j+1] == -1))

                if not pd.isnull(wgd_plus_df_non_acro.iloc[i,j]) and not pd.isnull(wgd_plus_df_non_acro.iloc[i,j+1]):

                    if wgd_plus_df_non_acro.iloc[i,j] != 0 and wgd_plus_df_non_acro.iloc[i,j+1] == 0:
                        arm_level += 1

                    elif wgd_plus_df_non_acro.iloc[i,j] == 0 and wgd_plus_df_non_acro.iloc[i,j+1] != 0:
                        arm_level += 1

                    elif wgd_plus_df_non_acro.iloc[i,j] != 0 and wgd_plus_df_non_acro.iloc[i,j+1] != 0 and \
                            wgd_plus_df_non_acro.iloc[i,j] != wgd_plus_df_non_acro.iloc[i,j+1]: #one gain, one loss
                        arm_level += 2


                else:

                    if not pd.isnull(wgd_plus_df_non_acro.iloc[i,j]) and wgd_plus_df_non_acro.iloc[i,j] != 0:
                        arm_level += 1

                    elif not pd.isnull(wgd_plus_df_non_acro.iloc[i,j+1]) and wgd_plus_df_non_acro.iloc[i,j+1] != 0:
                        arm_level += 1

                #print("j:", j, arm_level, total_wgd_plus, wgd_plus_df_non_acro.iloc[i,j], wgd_plus_df_non_acro.iloc[i,j+1])

            #print(type, arm_level, total_wgd_plus, total_wgd_plus / (total_wgd_plus + arm_level)) #, (wgd_plus_df.CCLE_ID)[i]
            value_series = [[type, arm_level, total_wgd_plus, total_wgd_plus / (total_wgd_plus + arm_level)]]

            #value_series.index = ["Type", "Arm_anp", "Whole_chr_anp", "wca_frac"]
            wgd_plus_fracs = wgd_plus_fracs.append(value_series, ignore_index=True)
            #print(wgd_plus_fracs)

            whole_chr_anp_frac_wgd_plus.append(total_wgd_plus / (total_wgd_plus + arm_level))

        for i in range(0, len(wgd_minus_df_non_acro.index)):
            total_wgd_minus = 0
            arm_level = 0
            for j in range(0, 34, 2):
                total_wgd_minus += np.sum(np.logical_and(wgd_minus_df_non_acro.iloc[i,j] == 1, wgd_minus_df_non_acro.iloc[i,j+1] == 1)) \
                + np.sum(np.logical_and(wgd_minus_df_non_acro.iloc[i,j] == -1, wgd_minus_df_non_acro.iloc[i,j+1] == -1))

                if not pd.isnull(wgd_minus_df_non_acro.iloc[i,j]) and not pd.isnull(wgd_minus_df_non_acro.iloc[i,j+1]):

                    if wgd_minus_df_non_acro.iloc[i,j] != 0 and wgd_minus_df_non_acro.iloc[i,j+1] == 0:
                        arm_level += 1

                    elif wgd_minus_df_non_acro.iloc[i,j] == 0 and wgd_minus_df_non_acro.iloc[i,j+1] != 0:
                        arm_level += 1

                    elif wgd_minus_df_non_acro.iloc[i,j] != 0 and wgd_minus_df_non_acro.iloc[i,j+1] != 0 and \
                            wgd_minus_df_non_acro.iloc[i,j] != wgd_minus_df_non_acro.iloc[i,j+1]: #one gain, one loss
                        arm_level += 2


                else:

                    if not pd.isnull(wgd_minus_df_non_acro.iloc[i,j]) and wgd_minus_df_non_acro.iloc[i,j] != 0:
                        arm_level += 1

                    elif not pd.isnull(wgd_minus_df_non_acro.iloc[i,j+1]) and wgd_minus_df_non_acro.iloc[i,j+1] != 0:
                        arm_level += 1

            value_series = [[type, arm_level, total_wgd_minus, total_wgd_minus / (total_wgd_minus + arm_level)]]

            wgd_minus_fracs = wgd_minus_fracs.append(value_series, ignore_index=True)
            #wgd_minus_fracs.columns = ["Type", "Arm_anp", "Whole_chr_anp", "wca_frac"]
            whole_chr_anp_frac_wgd_minus.append(total_wgd_minus/(total_wgd_minus+arm_level))

        whole_chr_anp_frac_sample = whole_chr_anp_frac_wgd_plus + whole_chr_anp_frac_wgd_minus

        whole_chr_anp_frac_sample = pd.Series(whole_chr_anp_frac_sample)
        wgd_status = ["WGD+"]*len(wgd_plus_df_non_acro.index) + ["WGD-"]*len(wgd_minus_df_non_acro.index)
        wgd_status = pd.Series(wgd_status)
        type_vec = pd.Series([type]*(len(wgd_plus_df_non_acro.index) + len(wgd_minus_df_non_acro.index)))

        whole_chr_anp_frac_sample = pd.concat([whole_chr_anp_frac_sample, wgd_status, type_vec], axis=1)

        whole_chr_anp_frac_sample.columns = ["Fraction", "WGD status", "Type"]
        pan_cancer_df = pd.concat([pan_cancer_df, whole_chr_anp_frac_sample], axis=0)


    #type_wise_df.to_csv(os.path.join(src_folder, "Type_wise_df.tsv"), sep = "\t")
    plot_boxplot(pan_cancer_df, "Pan-cancer whole-chr aneuploidy fraction (sample-wise)")
    pan_cancer_df.to_csv("C:\\Users\\kiwii\\Documents\\Kavya\\Uri\\Work\\Tasks\\CCLE analysis\\Prelim\\wca_samplewise_no-ke97.tsv", sep = "\t", index=0)
    wgd_plus_fracs.columns = ["Type", "Arm_anp", "Whole_chr_anp", "wca_frac"]
    wgd_minus_fracs.columns = ["Type", "Arm_anp", "Whole_chr_anp", "wca_frac"]

    analyze_whole_chr_anp_frac_samplewise_all(wgd_plus_fracs, wgd_minus_fracs)

def analyze_whole_chr_anp_frac_samplewise_all(wgd_plus_fracs, wgd_minus_fracs):

    wgd_status = ["WGD+"]*len(wgd_plus_fracs.index) + ["WGD-"]*len(wgd_minus_fracs.index)
    print(wgd_plus_fracs.head())
    print(wgd_minus_fracs.head())
    df_all = pd.concat([wgd_plus_fracs, wgd_minus_fracs], axis = 0, ignore_index=True)
    #df_all = pd.concat([df_all, pd.Series(wgd_status)], axis = 1)
    df_all["wgd_status"] = wgd_status
    print(df_all.head())
    df_all.to_csv(os.path.join(src_folder, "Type_wise_df_wca_frac.tsv"), sep = "\t", index = 0)

    plt.figure(figsize=(30, 10))
    ax = sns.boxplot(x="wgd_status", y="wca_frac", data=df_all, hue="wgd_status", palette="Set1", showfliers=False)

    ax.set_ylim((0, 1.1))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    pairs = [('WGD+', 'WGD-')]

    add_stat_annotation(ax, x="wgd_status", y="wca_frac", data=df_all,  box_pairs= pairs,test='t-test_ind',
                        text_format='star', loc='inside', verbose=0)

    plt.show()


def analyze_aneuploidy_score(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Fraction", "WGD status", "Type"])
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    score_wgd_plus = []
    score_wgd_minus = []
    types_vector_wgd_plus = []
    types_vector_wgd_minus = []

    for type in types:

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)

        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:, 4:]
        print(wgd_minus_df.columns)
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:, 4:]



        total_wgd_plus = []
        total_wgd_minus = []

        for i in range(0, len(wgd_plus_df.index)):
            total_wgd_plus.append(sum(wgd_plus_df.iloc[i,:] != 0))

        for i in range(0, len(wgd_minus_df.index)):
            total_wgd_minus.append(sum(wgd_minus_df.iloc[i,:] != 0))

        score_wgd_plus.extend(total_wgd_plus)
        score_wgd_minus.extend(total_wgd_minus)
        types_vector_wgd_plus.extend([type]*len(total_wgd_plus))
        types_vector_wgd_minus.extend([type] * len(total_wgd_minus))

    print(score_wgd_minus, score_wgd_plus)


    all_scores = score_wgd_plus + score_wgd_minus
    wgd_status = (["WGD+"] * len(score_wgd_plus)) + (["WGD-"] * len(score_wgd_minus))
    type_vector = types_vector_wgd_plus + types_vector_wgd_minus
    df_scores = pd.concat([pd.Series(all_scores), pd.Series(wgd_status), pd.Series(type_vector)], axis= 1, ignore_index= True)

    df_scores.columns = ["Score", "WGD status", "Type"]
    print(df_scores)

    df_scores.to_csv(
        "C:\\Users\\kiwii\\Documents\\Kavya\\Uri\\Work\\Tasks\\CCLE analysis\\Prelim\\anp_scores_samplewise_no-ke97.tsv", sep="\t",
        index=0)

    plt.figure(figsize=(30, 10))
    fig = sns.boxplot(x="Type", y="Score", data=df_scores, hue="WGD status", palette="Set1", showfliers=False)

    for tick in fig.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in fig.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    #ax.set(ylim = (0,1000))
    axes = fig.axes
    #axes.set_ylim(0, 1500)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    pairs = [((x,"WGD+"),(x,"WGD-")) for x in types]
    #pairs = [(("WGD+", x), ("WGD-", x)) for x in types]

    add_stat_annotation(fig, x="Type", y="Score", data=df_scores, hue="WGD status" ,box_pairs= pairs,test='t-test_ind',
                        text_format='star', loc='inside', verbose=0)

    plt.show()

def compute_wgd_plus_frac(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Type", "WGD+ count", "WGD- count", "WGD+ fraction"])
    pan_cancer_df = pd.DataFrame()
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    score_wgd_plus = []
    score_wgd_minus = []
    types_vector_wgd_plus = []
    types_vector_wgd_minus = []

    for type in types:

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)

        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:, 4:]
        print(wgd_minus_df.columns)
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:, 4:]

        pan_cancer_df = pan_cancer_df.append({"Type:":type, "WGD+ count": len(wgd_plus_df.index), "WGD- count":len(wgd_minus_df.index), "WGD+ fraction":len(wgd_plus_df.index)/(len(wgd_plus_df.index)+len(wgd_minus_df.index))}, ignore_index=True)

    #pan_cancer_df.columns = ["Type", "WGD+ count", "WGD- count", "WGD+ fraction"]

    pan_cancer_df.to_csv(
    "C:\\Users\\kiwii\\OneDrive\\Documents\\Kavya\\Uri\\Work\\Tasks\\CCLE analysis\\Prelim\\wgd_plus_frac_no-ke97.tsv", sep="\t",
    index=0)



def analyze_aneuploidy_score_experimental(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Fraction", "WGD status", "Type"])
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    score_wgd_plus = []
    score_wgd_minus = []
    types_vector_wgd_plus = []
    types_vector_wgd_minus = []

    for type in types:

        print(type)

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)

        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:, 4:]
        print(wgd_minus_df.columns)
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:, 4:]



        total_wgd_plus = []
        total_wgd_minus = []

        flag = 0

        for j in range(0, len(wgd_plus_df.index)):
        #for j in range(0, 1):
            #print("j: ", j)
            count_anp_sample = 0
            i = 0
            while i < 39:
            #total_wgd_plus.append(sum(wgd_plus_df.iloc[i,:] != 0))
                #print("i:", i)
                if i != 38:
                    if re.findall(r'\d+', wgd_plus_df.columns[i]) == re.findall(r'\d+', wgd_plus_df.columns[i+1]):
                        if (wgd_plus_df.iloc[j,i] == wgd_plus_df.iloc[j,i+1]) and (wgd_plus_df.iloc[j,i] != 0) and (wgd_plus_df.iloc[j,i] != None) and (wgd_plus_df.iloc[j,i+1] != None):
                            count_anp_sample += 1
                        else:
                            if wgd_plus_df.iloc[j,i] != None:
                                count_anp_sample += abs(wgd_plus_df.iloc[j,i])

                            if wgd_plus_df.iloc[j,i+1] != None:
                                count_anp_sample += abs(wgd_plus_df.iloc[j,i+1])
                        flag = 2
                    else: #acrocentric chromosomes
                        if i in [24, 25, 26, 37, 38]:
                            if wgd_plus_df.iloc[j,i] != None:
                                count_anp_sample += abs(wgd_plus_df.iloc[j,i])
                            flag = 1
                    #print(wgd_plus_df.columns[i], count_anp_sample, flag)
                    if flag == 2:
                        i += 2
                    else:
                        i += 1
                elif i == 38: #if i = 38 ie we are at the last column
                    if wgd_plus_df.iloc[j, i] != None:
                        count_anp_sample += abs(wgd_plus_df.iloc[j, i])
                    flag = 1
                    i += 1
                    #print(wgd_plus_df.columns[i], count_anp_sample, flag)


                #print(wgd_plus_df.columns[i], count_anp_sample)
            total_wgd_plus.append(count_anp_sample)


        flag = 0

        for j in range(0, len(wgd_minus_df.index)):
            #print("j: ", j)
            count_anp_sample = 0
            i = 0
            while i < 39:
            #total_wgd_plus.append(sum(wgd_plus_df.iloc[i,:] != 0))
                #print("i: ", i)
                if i != 38:
                    if re.findall(r'\d+', wgd_minus_df.columns[i]) == re.findall(r'\d+', wgd_minus_df.columns[i+1]):
                        if (wgd_minus_df.iloc[j,i] == wgd_minus_df.iloc[j,i+1]) and (wgd_minus_df.iloc[j,i] != 0) and (wgd_minus_df.iloc[j,i] != None) and (wgd_minus_df.iloc[j,i+1] != None):
                            count_anp_sample += 1
                        else:
                            if wgd_minus_df.iloc[j,i] != None:
                                count_anp_sample += abs(wgd_minus_df.iloc[j,i])

                            if wgd_minus_df.iloc[j,i+1] != None:
                                count_anp_sample += abs(wgd_minus_df.iloc[j,i+1])
                        flag = 2
                    else: #acrocentric chromosomes
                        if i in [24, 25, 26, 37, 38]:
                            if wgd_minus_df.iloc[j,i] != None:
                                count_anp_sample += abs(wgd_minus_df.iloc[j,i])
                            flag = 1
                    if flag == 2:
                        i += 2
                    else:
                        i += 1

                elif i == 38:  # if i = 38 ie we are at the last column
                    if wgd_minus_df.iloc[j, i] != None:
                        count_anp_sample += abs(wgd_minus_df.iloc[j, i])
                    flag = 1
                    i += 1
                    #print(wgd_plus_df.columns[i], count_anp_sample, flag)


            total_wgd_minus.append(count_anp_sample)


        # for i in range(0, len(wgd_minus_df.index)):
        #     total_wgd_minus.append(sum(wgd_minus_df.iloc[i,:] != 0))

        score_wgd_plus.extend(total_wgd_plus)
        score_wgd_minus.extend(total_wgd_minus)
        types_vector_wgd_plus.extend([type]*len(total_wgd_plus))
        types_vector_wgd_minus.extend([type] * len(total_wgd_minus))

    print(score_wgd_minus, score_wgd_plus)


    all_scores = score_wgd_plus + score_wgd_minus
    wgd_status = (["WGD+"] * len(score_wgd_plus)) + (["WGD-"] * len(score_wgd_minus))
    type_vector = types_vector_wgd_plus + types_vector_wgd_minus
    df_scores = pd.concat([pd.Series(all_scores), pd.Series(wgd_status), pd.Series(type_vector)], axis= 1, ignore_index= True)

    df_scores.columns = ["Score", "WGD status", "Type"]
    print(df_scores)

    df_scores.to_csv(
        "C:\\Users\\kiwii\\OneDrive\\Documents\\Kavya\\Uri\\Work\\Tasks\\CCLE analysis\\Prelim\\anp_scores_samplewise_no-ke97_wca.tsv", sep="\t",
        index=0)

    plt.figure(figsize=(30, 10))
    fig = sns.boxplot(x="Type", y="Score", data=df_scores, hue="WGD status", palette="Set1", showfliers=False)

    for tick in fig.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in fig.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    #ax.set(ylim = (0,1000))
    axes = fig.axes
    #axes.set_ylim(0, 1500)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    pairs = [((x,"WGD+"),(x,"WGD-")) for x in types]
    #pairs = [(("WGD+", x), ("WGD-", x)) for x in types]

    add_stat_annotation(fig, x="Type", y="Score", data=df_scores, hue="WGD status" ,box_pairs= pairs,test='t-test_ind',
                        text_format='star', loc='inside', verbose=0)

    plt.show()







def analyze_aneuploidy_score_all_new(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Fraction", "WGD status", "Type"])
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    score_wgd_plus = []
    score_wgd_minus = []
    types_vector_wgd_plus = []
    types_vector_wgd_minus = []

    for type in types:

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)

        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:, 3:]
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:, 3:]



        total_wgd_plus = []
        total_wgd_minus = []

        for i in range(0, len(wgd_plus_df.index)):
            total_wgd_plus.append(sum(wgd_plus_df.iloc[i,:] != 0))

        for i in range(0, len(wgd_minus_df.index)):
            total_wgd_minus.append(sum(wgd_minus_df.iloc[i,:] != 0))

        score_wgd_plus.extend(total_wgd_plus)
        score_wgd_minus.extend(total_wgd_minus)
        types_vector_wgd_plus.extend([type]*len(total_wgd_plus))
        types_vector_wgd_minus.extend([type] * len(total_wgd_minus))

    print(score_wgd_minus, score_wgd_plus)


    all_scores = score_wgd_plus + score_wgd_minus
    wgd_status = (["WGD+"] * len(score_wgd_plus)) + (["WGD-"] * len(score_wgd_minus))
    type_vector = types_vector_wgd_plus + types_vector_wgd_minus
    df_scores = pd.concat([pd.Series(all_scores), pd.Series(wgd_status), pd.Series(type_vector)], axis= 1, ignore_index= True)

    df_scores.columns = ["Score", "WGD status", "Type"]
    print(df_scores)

    plt.figure(figsize=(30, 10))
    fig = sns.boxplot(x="WGD status", y="Score", data=df_scores, hue="WGD status", palette="Set1", showfliers=False)

    for tick in fig.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in fig.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    #ax.set(ylim = (0,1000))
    axes = fig.axes
    #axes.set_ylim(0, 1500)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    pairs = [("WGD+", "WGD-")]
    #pairs = [(("WGD+", x), ("WGD-", x)) for x in types]

    add_stat_annotation(fig, x="WGD status", y="Score", data=df_scores,box_pairs= pairs,test='t-test_ind',
                        text_format='star', loc='inside', verbose=0)

    plt.show()

def analyze_aneuploidy_score_all(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Fraction", "WGD status", "Type"])
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    score_wgd_plus = []
    score_wgd_minus = []

    for type in types:

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)

        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:, 3:]
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:, 3:]



        total_wgd_plus = 0
        total_wgd_minus = 0

        for i in range(0, len(wgd_plus_df.index)):
            total_wgd_plus += sum(wgd_plus_df.iloc[i,:] != 0)

        for i in range(0, len(wgd_minus_df.index)):
            total_wgd_minus += sum(wgd_minus_df.iloc[i,:] != 0)

        score_wgd_plus.append(total_wgd_plus)
        score_wgd_minus.append(total_wgd_minus)

    print(score_wgd_minus, score_wgd_plus)


    all_scores = score_wgd_plus + score_wgd_minus
    wgd_status = (["WGD+"] * len(score_wgd_plus)) + (["WGD-"] * len(score_wgd_minus))
    df_scores = pd.concat([pd.Series(all_scores), pd.Series(wgd_status)], axis= 1, ignore_index= True)

    df_scores.columns = ["Score", "WGD status"]
    print(df_scores)

    plt.figure(figsize=(15, 5))
    fig = sns.boxplot(x="WGD status", y="Score", data=df_scores, hue="WGD status", palette="Set1", showfliers=False)

    for tick in fig.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in fig.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    #ax.set(ylim = (0,1000))
    axes = fig.axes
    axes.set_ylim(0, 1500)
    axes.set_xticklabels(fig.get_xticklabels(), rotation=90)
    pairs = [('WGD+', 'WGD-')]

    add_stat_annotation(fig, x="WGD status", y="Score", data=df_scores,  box_pairs= pairs,test='t-test_ind',
                        text_format='star', loc='inside', verbose=0)

    plt.show()

def analyze_loss_frac(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Fraction", "WGD status", "Type"])
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    fractions_wgd_plus = []
    fractions_wgd_minus = []
    types_vector_wgd_plus = []
    types_vector_wgd_minus = []

    for type in types:

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)

        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:, 3:]
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:, 3:]



        frac_wgd_plus = 0
        frac_wgd_minus = 0

        for i in range(0, len(wgd_plus_df.index)):
            loss = sum(wgd_plus_df.iloc[i,:] == -1)
            gain = sum(wgd_plus_df.iloc[i, :] == 1)
            if loss+gain != 0:
                frac_wgd_plus = loss/(loss+gain)
            else:
                frac_wgd_plus = 0
            fractions_wgd_plus.append(frac_wgd_plus)

        for i in range(0, len(wgd_minus_df.index)):
            loss = sum(wgd_minus_df.iloc[i,:] == -1)
            gain = sum(wgd_minus_df.iloc[i, :] == 1)
            if loss+gain != 0:
                frac_wgd_minus = loss/(loss+gain)
            else:
                frac_wgd_minus = 0
            fractions_wgd_minus.append(frac_wgd_minus)

        types_vector_wgd_plus.extend([type]*len(wgd_plus_df.index))
        types_vector_wgd_minus.extend([type] * len(wgd_minus_df.index))

    print(fractions_wgd_minus, fractions_wgd_plus)


    all_scores = fractions_wgd_plus + fractions_wgd_minus
    wgd_status = (["WGD+"] * len(fractions_wgd_plus)) + (["WGD-"] * len(fractions_wgd_minus))
    type_vector = types_vector_wgd_plus + types_vector_wgd_minus
    print(len(all_scores), len(wgd_status), len(type_vector))
    df_scores = pd.concat([pd.Series(all_scores), pd.Series(wgd_status), pd.Series(type_vector)], axis= 1, ignore_index= True)

    df_scores.columns = ["Score", "WGD status", "Type"]
    print(df_scores)

    df_scores.to_csv(
        "C:\\Users\\kiwii\\Documents\\Kavya\\Uri\\Work\\Tasks\\CCLE analysis\\Prelim\\loss_frac.tsv", sep="\t",
        index=0)

    plt.figure(figsize=(30, 10))
    fig = sns.boxplot(x="Type", y="Score", data=df_scores, hue="WGD status", palette="Set1", showfliers=False)

    for tick in fig.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in fig.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    #ax.set(ylim = (0,1000))
    axes = fig.axes
    axes.set_ylim(0, 1.1)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    pairs = [((x,"WGD+"),(x,"WGD-")) for x in types]
    #pairs = [(("WGD+", x), ("WGD-", x)) for x in types]

    add_stat_annotation(fig, x="Type", y="Score", data=df_scores, hue="WGD status" ,box_pairs= pairs,test='t-test_ind',
                        text_format='star', loc='inside', verbose=0)

    plt.show()

def analyze_loss_frac_all(types):
    """Generates analysis of fraction of whole chromosome aneuploidies, relative to fractions
    of whole-chromosome and arm-level aneuploidies
    Calculates a fraction corresponding to each sample"""
    pan_cancer_df = pd.DataFrame(columns=["Fraction", "WGD status", "Type"])
    type_wise_df = pd.DataFrame()

    wgd_plus_fracs = pd.DataFrame((pd.np.empty((0, 4))))
    wgd_minus_fracs = pd.DataFrame((pd.np.empty((0, 4))))

    fractions_wgd_plus = []
    fractions_wgd_minus = []
    types_vector_wgd_plus = []
    types_vector_wgd_minus = []

    for type in types:

        type_df = input_df[input_df.CCLE_ID.str.contains(type)]
        type_wise_df = pd.concat([type_wise_df, type_df], axis=0)

        wgd_minus_df = type_df[type_df.ploidy < 2.5].iloc[:, 3:]
        wgd_plus_df = type_df[type_df.ploidy > 3].iloc[:, 3:]



        frac_wgd_plus = 0
        frac_wgd_minus = 0

        for i in range(0, len(wgd_plus_df.index)):
            loss = sum(wgd_plus_df.iloc[i,:] == -1)
            gain = sum(wgd_plus_df.iloc[i, :] == 1)
            if loss+gain != 0:
                frac_wgd_plus = loss/(loss+gain)
            else:
                frac_wgd_plus = 0
            fractions_wgd_plus.append(frac_wgd_plus)

        for i in range(0, len(wgd_minus_df.index)):
            loss = sum(wgd_minus_df.iloc[i,:] == -1)
            gain = sum(wgd_minus_df.iloc[i, :] == 1)
            if loss+gain != 0:
                frac_wgd_minus = loss/(loss+gain)
            else:
                frac_wgd_minus = 0
            fractions_wgd_minus.append(frac_wgd_minus)

        types_vector_wgd_plus.extend([type]*len(wgd_plus_df.index))
        types_vector_wgd_minus.extend([type] * len(wgd_minus_df.index))

    print(fractions_wgd_minus, fractions_wgd_plus)
    print(statistics.mean(fractions_wgd_plus), statistics.mean(fractions_wgd_minus))


    all_scores = fractions_wgd_plus + fractions_wgd_minus
    wgd_status = (["WGD+"] * len(fractions_wgd_plus)) + (["WGD-"] * len(fractions_wgd_minus))
    type_vector = types_vector_wgd_plus + types_vector_wgd_minus
    print(len(all_scores), len(wgd_status), len(type_vector))
    df_scores = pd.concat([pd.Series(all_scores), pd.Series(wgd_status), pd.Series(type_vector)], axis= 1, ignore_index= True)

    df_scores.columns = ["Score", "WGD status", "Type"]
    print(df_scores)

    plt.figure(figsize=(30, 10))
    fig = sns.boxplot(x="WGD status", y="Score", data=df_scores, hue="WGD status", palette="Set1", showfliers=False)

    for tick in fig.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in fig.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    #ax.set(ylim = (0,1000))
    axes = fig.axes
    axes.set_ylim(0, 1.1)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    pairs = [("WGD+", "WGD-")]
    #pairs = [(("WGD+", x), ("WGD-", x)) for x in types]

    add_stat_annotation(fig, x="WGD status", y="Score", data=df_scores,box_pairs= pairs,test='t-test_ind',
                        text_format='star', loc='inside', verbose=0)

    plt.show()

if __name__ == "__main__":
    global src_folder
    src_folder = "C:\\Users\\kiwii\\OneDrive\\Documents\\Kavya\\Uri\\Work\\Tasks\\CCLE analysis\\Prelim"
    #input_df = pd.read_csv(os.path.join(src_folder, "Table_S1_arm_calls.txt"), sep = "\t")
    input_df = pd.read_csv(os.path.join(src_folder, "Type_wise_df_no-ke97.tsv"), sep="\t")
    types = data_exploration(input_df)
    print(types)
    #analyze_whole_chr_anp_frac_samplewise(types)
    #analyze_aneuploidy_score_all_new(types)
    #analyze_loss_frac_all(types)
    #analyze_loss_frac(types)
    #analyze_aneuploidy_score(types)
    analyze_aneuploidy_score_experimental(types)
    #compute_wgd_plus_frac(types)
