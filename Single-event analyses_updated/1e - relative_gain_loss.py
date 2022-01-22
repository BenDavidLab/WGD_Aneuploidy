import pandas as pd
import os
import constants
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import hypergeom

global input_df, input_df_pm
input_df = pd.read_csv(os.path.join("C:\\Users\\kiwii\\PycharmProjects\\pan_cancer_aneuploidy\\datasets\\Taylor_et_al._Arm-Level_WGD_TCGA_data.txt"), sep='\t', index_col=0, header=0)
input_df_pm = pd.read_csv(os.path.join("C:\\Users\\kiwii\\PycharmProjects\\pan_cancer_aneuploidy\\datasets\\Taylor_et_al._Arm-Level_WGD_TCGA_plus_minus_data.txt"), sep='\t', index_col=0, header=0)

#==============Intra-tumor trend analysis==================

##===Absolute and relative gain and loss fraction===

def plot_grouped_bars_gain(wgd_plus_frac, wgd_minus_frac, fname, xlabels, y_lim):
    """Plots grouped bar chart of gain or loss fractions
    of WGD+ and WGD- side-by-side for ease of comparison"""
    pos_wgd_minus = np.arange(len(wgd_minus_frac))
    pos_wgd_plus = [x + constants.b_width for x in pos_wgd_minus]

    plt.figure(figsize=(12, 4))
    plt.bar(pos_wgd_minus, wgd_minus_frac, color=  '#fa9c9c', width=constants.b_width)
    plt.bar(pos_wgd_plus, wgd_plus_frac, color=  '#950101', width=constants.b_width)
    plt.xlabel('Chromosome arm',fontsize=16.5)
    #plt.ylabel('Relative prevalence out of all chromosomal gains')
    plt.ylabel('Relative aneuploidy\nprevalence (gains)', fontsize=16.5)
    plt.xticks([x_pos + constants.b_width for x_pos in range(len(pos_wgd_plus))], xlabels,rotation = 90,fontsize=16.5)
    plt.yticks(fontsize=16.5)
    plt.tick_params(axis="x", which="both", bottom=False, top=False)
    plt.ylim((0, y_lim))
    #plt.title(fname.replace("_", " "))
    ##plt.title((fname.split("_"))[0])
    plt.legend(["WGD-", "WGD+"],fontsize=16.5)
    plt.savefig("{}_ax.pdf".format(fname), bbox_inches='tight', dpi=300)
    plt.clf()

def plot_grouped_bars_loss(wgd_plus_frac, wgd_minus_frac, fname, xlabels, y_lim):
    """Plots grouped bar chart of gain or loss fractions
    of WGD+ and WGD- side-by-side for ease of comparison"""
    pos_wgd_minus = np.arange(len(wgd_minus_frac))
    pos_wgd_plus = [x + constants.b_width for x in pos_wgd_minus]

    #fig = plt.figure(figsize=(12, 6))
    fig, ax = plt.subplots(1,1,figsize=(12, 4))
    ax.set_facecolor('white')

    plt.bar(pos_wgd_minus, wgd_minus_frac, color=  '#9cb3fa', width=constants.b_width)
    plt.bar(pos_wgd_plus, wgd_plus_frac, color=  '#00269c', width=constants.b_width)
    plt.xlabel('Chromosome arm',fontsize=16.5)
    #plt.ylabel('Relative prevalence out of all chromosomal losses')
    plt.ylabel('Relative aneuploidy\nprevalence (losses)', fontsize=16.5)
    plt.xticks([x_pos + constants.b_width for x_pos in range(len(pos_wgd_plus))], xlabels,rotation = 90,fontsize=16.5)
    plt.tick_params(axis="x", which="both", bottom=False, top=False)
    plt.yticks(fontsize=16.5)
    plt.ylim((0, y_lim))
    #plt.title(fname.replace("_", " "))
    ##plt.title((fname.split("_"))[0])
    plt.legend(["WGD-", "WGD+"],fontsize=16.5)
    #plt.savefig("{}_ax.svg".format(fname), bbox_inches='tight', format='svg', dpi=300)
    plt.savefig("{}_ax.pdf".format(fname), bbox_inches='tight', dpi=300)
    plt.clf()

def analyze_absolute_frac(type):
    """Analyzes absolute gain and loss fractions for all chromosome arms
    of WGD+ and WGD- for the given tumor type
    Absolute fraction refers to the fraction of aneuploidies across all samples,
    that each chromosome accounts for"""
    type_df = input_df[input_df.Type == type]
    wgd_plus_df = type_df[type_df.Genome_doublings > 0]
    wgd_minus_df = type_df[type_df.Genome_doublings == 0]

    wgd_plus_abs_gain_frac = []
    wgd_minus_abs_gain_frac = []

    wgd_plus_abs_loss_frac = []
    wgd_minus_abs_loss_frac = []

    for chr_arm in type_df.columns[12:51]:
        if constants.GAIN: #if analysis for absolute gain frac needs to be done (can be set in the 'constants.py' module)
            wgd_plus_arm_gain_frac = len(wgd_plus_df[wgd_plus_df[chr_arm] == 1].index)/len(wgd_plus_df.index)
            wgd_minus_arm_gain_frac = len(wgd_minus_df[wgd_minus_df[chr_arm] == 1].index)/len(wgd_minus_df.index)
            wgd_plus_abs_gain_frac.append(wgd_plus_arm_gain_frac)
            wgd_minus_abs_gain_frac.append(wgd_minus_arm_gain_frac)

            file_name_gain = "{}_abs_gain_frac".format(type)

        if constants.LOSS: #if analysis for absolute loss frac needs to be done (can be set in the 'constants.py' module)
            wgd_plus_arm_loss_frac = len(wgd_plus_df[wgd_plus_df[chr_arm] == -1].index)/len(wgd_plus_df.index)
            wgd_minus_arm_loss_frac = len(wgd_minus_df[wgd_minus_df[chr_arm] == -1].index)/len(wgd_minus_df.index)
            wgd_plus_abs_loss_frac.append(wgd_plus_arm_loss_frac)
            wgd_minus_abs_loss_frac.append(wgd_minus_arm_loss_frac)

            file_name_loss = "{}_abs_loss_frac".format(type)

    plot_grouped_bars_gain(wgd_plus_abs_gain_frac, wgd_minus_abs_gain_frac, file_name_gain, type_df.columns[12:51], 0.25)
    plot_grouped_bars_loss(wgd_plus_abs_loss_frac, wgd_minus_abs_loss_frac, file_name_loss, type_df.columns[12:51], 0.25)

def analyze_relative_frac(type):
    """Analyzes relative gain and loss fractions for all chromosome arms
    of WGD+ and WGD- for the given tumor type
    Relative fraction refers to the fraction of aneuploidies across all chromosomes,
    that each chromosome accounts for"""
    type_df = input_df[input_df.Type == type]
    wgd_plus_df = type_df[type_df.Genome_doublings > 0]
    wgd_minus_df = type_df[type_df.Genome_doublings == 0]

    wgd_plus_abs_gain_frac = []
    wgd_minus_abs_gain_frac = []

    wgd_plus_abs_loss_frac = []
    wgd_minus_abs_loss_frac = []

    total_chr_loss_wgd_minus = sum([len(wgd_minus_df[wgd_minus_df[chr_arm] == -1].index) for chr_arm in wgd_plus_df.columns[12:51]])
    total_chr_loss_wgd_plus = sum([len(wgd_plus_df[wgd_plus_df[chr_arm] == -1].index) for chr_arm in wgd_plus_df.columns[12:51]])

    total_chr_gain_wgd_minus = sum([len(wgd_minus_df[wgd_minus_df[chr_arm] == 1].index) for chr_arm in wgd_plus_df.columns[12:51]])
    total_chr_gain_wgd_plus = sum([len(wgd_plus_df[wgd_plus_df[chr_arm] == 1].index) for chr_arm in wgd_plus_df.columns[12:51]])

    for chr_arm in type_df.columns[12:51]:
        if constants.GAIN: #if analysis for relative gain frac needs to be done (can be set in the 'constants.py' module)
            wgd_plus_arm_gain_frac = len(wgd_plus_df[wgd_plus_df[chr_arm] == 1].index)/total_chr_gain_wgd_plus
            wgd_minus_arm_gain_frac = len(wgd_minus_df[wgd_minus_df[chr_arm] == 1].index)/total_chr_gain_wgd_minus
            wgd_plus_abs_gain_frac.append(wgd_plus_arm_gain_frac)
            wgd_minus_abs_gain_frac.append(wgd_minus_arm_gain_frac)

            file_name_gain = "{}_relative_gain_fraction".format(type)

        if constants.LOSS: #if analysis for relative loss frac needs to be done (can be set in the 'constants.py' module)
            wgd_plus_arm_loss_frac = len(wgd_plus_df[wgd_plus_df[chr_arm] == -1].index)/total_chr_loss_wgd_plus
            wgd_minus_arm_loss_frac = len(wgd_minus_df[wgd_minus_df[chr_arm] == -1].index)/total_chr_loss_wgd_minus
            wgd_plus_abs_loss_frac.append(wgd_plus_arm_loss_frac)
            wgd_minus_abs_loss_frac.append(wgd_minus_arm_loss_frac)

            file_name_loss = "{}_relative_loss_fraction".format(type)

    plot_grouped_bars_gain(wgd_plus_abs_gain_frac, wgd_minus_abs_gain_frac, file_name_gain, type_df.columns[12:51], 0.25)
    plot_grouped_bars_loss(wgd_plus_abs_loss_frac, wgd_minus_abs_loss_frac, file_name_loss, type_df.columns[12:51], 0.25)

analyze_relative_frac("UCEC")
