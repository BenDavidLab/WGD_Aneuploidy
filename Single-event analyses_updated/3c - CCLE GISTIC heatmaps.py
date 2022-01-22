import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import seaborn as sns

def read_data(DIR):
    "Read relevant data from gistic tables to make graph dataframes to send to heatmap function"

    gain_score_wgd_plus = pd.DataFrame()
    gain_score_wgd_minus = pd.DataFrame()

    loss_score_wgd_plus = pd.DataFrame()
    loss_score_wgd_minus = pd.DataFrame()

    gain_qval_wgd_plus = pd.DataFrame()
    gain_qval_wgd_minus = pd.DataFrame()

    loss_qval_wgd_plus = pd.DataFrame()
    loss_qval_wgd_minus = pd.DataFrame()

    tumor_types_order = ["LUNG", "OVARY", "SKIN", "BREAST", "URINARY_TRACT", "PANCREAS", "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "CENTRAL_NERVOUS_SYSTEM","BONE",
              "UPPER_AERODIGESTIVE_TRACT", "KIDNEY", "LARGE_INTESTINE", "LIVER", "STOMACH"]

    # tumor_types_order = ["COAD","ESCA",
    #                      "CESC", "HNSC",
    #                      "BLCA",
    #                      "GBM",
    #                      "ACC"]

    tumor_types_order.sort()  #we want alphabetical

    for tumor_type in tumor_types_order:
        wgd_minus = pd.read_csv(os.path.join(DIR, "{}_seg_file_hg19_2_wgd_minus.tsv//broad_significance_results.txt".format(tumor_type)), sep = "\t", header = 0) ###fill_in
        wgd_plus = pd.read_csv(os.path.join(DIR, "{}_seg_file_hg19_2_wgd_plus.tsv//broad_significance_results.txt".format(tumor_type)), sep = "\t", header = 0)  ###fill_in

        arm_names = wgd_plus['Arm']


        amp_score_wgd_plus = wgd_plus['Amp frequency']
        amp_score_wgd_minus = wgd_minus['Amp frequency']
        del_score_wgd_plus = wgd_plus['Del frequency']
        del_score_wgd_minus = wgd_minus['Del frequency']

        amp_qval_wgd_plus = wgd_plus['Amp q-value']
        amp_qval_wgd_minus = wgd_minus['Amp q-value']
        del_qval_wgd_plus = wgd_plus['Del q-value']
        del_qval_wgd_minus = wgd_minus['Del q-value']

        gain_score_wgd_minus = pd.concat([gain_score_wgd_minus, amp_score_wgd_minus], axis=1)
        gain_score_wgd_plus = pd.concat([gain_score_wgd_plus, amp_score_wgd_plus], axis=1)

        loss_score_wgd_minus = pd.concat([loss_score_wgd_minus, del_score_wgd_minus], axis=1)
        loss_score_wgd_plus = pd.concat([loss_score_wgd_plus, del_score_wgd_plus], axis=1)

        gain_qval_wgd_minus = pd.concat([gain_qval_wgd_minus, amp_qval_wgd_minus], axis=1)
        gain_qval_wgd_plus = pd.concat([gain_qval_wgd_plus, amp_qval_wgd_plus], axis=1)

        loss_qval_wgd_minus = pd.concat([loss_qval_wgd_minus, del_qval_wgd_minus], axis=1)
        loss_qval_wgd_plus = pd.concat([loss_qval_wgd_plus, del_qval_wgd_plus], axis=1)

    #plot_heatmap(amp_qval_wgd_minus, amp_qval_wgd_plus, amp_score_wgd_plus, "WGD+", arm_names)

    gain_qval_wgd_minus.columns = tumor_types_order
    gain_qval_wgd_plus.columns = tumor_types_order
    loss_qval_wgd_minus.columns = tumor_types_order
    loss_qval_wgd_plus.columns = tumor_types_order

    gain_score_wgd_minus.columns = tumor_types_order
    gain_score_wgd_plus.columns = tumor_types_order
    loss_score_wgd_minus.columns = tumor_types_order
    loss_score_wgd_plus.columns = tumor_types_order

    gain_qval_wgd_minus.index = arm_names
    gain_qval_wgd_plus.index = arm_names
    loss_qval_wgd_minus.index = arm_names
    loss_qval_wgd_plus.index = arm_names

    gain_score_wgd_minus.index = arm_names
    gain_score_wgd_plus.index = arm_names
    loss_score_wgd_minus.index = arm_names
    loss_score_wgd_plus.index = arm_names

    gain_qval_dfs = {"WGD+":gain_qval_wgd_plus, "WGD-": gain_qval_wgd_minus}
    loss_qval_dfs = {"WGD+": loss_qval_wgd_plus, "WGD-": loss_qval_wgd_minus}

    if HGLT_FLAG == 1:
        highlight_gain_wgd_minus, highlight_gain_wgd_plus = highlight_diff_cells_excl(gain_qval_dfs, tumor_types_order, THRESH)
        highlight_loss_wgd_minus, highlight_loss_wgd_plus = highlight_diff_cells_excl(loss_qval_dfs, tumor_types_order, THRESH)
    else:
        highlight_gain_wgd_minus, highlight_gain_wgd_plus = highlight_diff_cells_incl(gain_qval_dfs, tumor_types_order, THRESH)
        highlight_loss_wgd_minus, highlight_loss_wgd_plus = highlight_diff_cells_incl(loss_qval_dfs, tumor_types_order, THRESH)

    if UNIFIED == 1:
        #print('WGD+')
        unified_wgd_plus = make_unified_df(gain_qval_wgd_plus, loss_qval_wgd_plus, gain_score_wgd_plus, loss_score_wgd_plus, THRESH)
        #print("WGD-")
        unified_wgd_minus = make_unified_df(gain_qval_wgd_minus, loss_qval_wgd_minus, gain_score_wgd_minus, loss_score_wgd_minus,THRESH)
        #unified_wgd_plus.to_csv("unified_df_wgd_plus.tsv", sep = "\t")
        unified_dfs = {"WGD+": unified_wgd_plus, "WGD-": unified_wgd_minus}
        #highlight_wgd_minus, highlight_wgd_plus = find_unified_highlights_excl_freq(unified_dfs, gain_qval_wgd_minus, gain_qval_wgd_plus, loss_qval_wgd_minus, loss_qval_wgd_plus, tumor_types_order, THRESH)

        plot_heatmap_pan_cancer_anp_diff(unified_wgd_minus, "try_0.9_ccle_gistic_unified_wgd_minus_all_nohglt_freq_{}".format(RECUR), [], -1,
                                         "coolwarm")
        plot_heatmap_pan_cancer_anp_diff(unified_wgd_plus, "try_0.9_ccle_gistic_unified_wgd_plus_all_nohglt_freq_{}".format(RECUR), [], -1, "coolwarm")

    #gain_score_wgd_plus.to_csv("ccle_Gain_score_wgd_plus_all.tsv", sep = "\t")
    plot_heatmap_pan_cancer_anp_diff(gain_score_wgd_minus, "gistic_ccle_0.9_wgd_minus_gain_all", highlight_gain_wgd_minus, 0, "Reds")
    plot_heatmap_pan_cancer_anp_diff(gain_score_wgd_plus, "gistic_ccle_0.9_wgd_plus_gain_all", highlight_gain_wgd_plus, 0, "Reds")

    plot_heatmap_pan_cancer_anp_diff(loss_score_wgd_minus, "gistic_ccle_0.9_wgd_minus_loss_all", highlight_loss_wgd_minus, 0, "Blues")
    plot_heatmap_pan_cancer_anp_diff(loss_score_wgd_plus, "gistic_ccle_0.9_wgd_plus_loss_all", highlight_loss_wgd_plus, 0, "Blues")

def make_unified_df(gain_qval_df, loss_qval_df, gain_score_df, loss_score_df, THRESH):
    "Combine gain and loss dfs into a single df where we retain the more sig event in that arm"
    gain_qval_df.drop("21p", axis=0, inplace=True)
    gain_score_df.drop("21p", axis=0, inplace=True)
    loss_qval_df.drop("21p", axis=0, inplace=True)
    loss_score_df.drop("21p", axis=0, inplace=True)

    print(loss_qval_df.index, len(loss_qval_df.index))
    data = np.array([np.arange(39)] * 14).T
    unified_df = pd.DataFrame(data, index = loss_qval_df.index, columns = loss_qval_df.columns)
    #unified_df.drop("21p", axis=0, inplace=True)
    #print(len(unified_df.index))
    print(gain_qval_df.shape, loss_qval_df.shape, gain_score_df.shape, loss_score_df.shape)
    print(gain_qval_df.index)

    for i, col in enumerate(gain_qval_df.columns):
        #print("i: ", i, " ", col)
        for j, (x,y) in enumerate(zip(gain_qval_df.loc[:, col], loss_qval_df.loc[:, col])):
            #if ((x > 0 and y < 0) or  (x < 0 and y > 0)) and ((abs(x) < THRESH) or (abs(y) < THRESH)):
            #print("j: ", j, x, y, loss_score_df.iloc[j, i], gain_score_df.iloc[j, i])
            #print(unified_df.iloc[j, i])
            if x < y:
                #print(gain_score_df.iloc[j, i], loss_score_df.iloc[j, i], gain_score_df.iloc[j, i] - loss_score_df.iloc[j, i])
                unified_df.iloc[j, i] = (gain_score_df.iloc[j, i] - loss_score_df.iloc[j, i])
            else:
                if j == 7:
                    print(gain_score_df.iloc[j, i], loss_score_df.iloc[j, i],
                      -1*(loss_score_df.iloc[j, i]) - (gain_score_df.iloc[j, i]))
                #print(j, i)
                unified_df.iloc[j, i] = -1*(loss_score_df.iloc[j, i] - gain_score_df.iloc[j, i])
    return unified_df

def plot_heatmap_pan_cancer_anp_diff(graph_df, fname, highlight_cells, vmin, colour): #vmin says if unicolour or two-colour, colour = "vlag" for twocolour
    """Plots a heatmap showing pan-cancer trends in chromosome gains and losses
    and highlights events with different occurrence patterns ion the two groups"""
    #g = sns.clustermap(graph_df, cmap="vlag", vmin=-1, vmax=1, row_cluster=0, col_cluster=0)
    print(highlight_cells)
    sns.set(font_scale=0.7)
    #graph_df.values.astype(np.float64)))
    for y in graph_df.columns:
        if graph_df[y].dtype == object:
            graph_df[y] = pd.to_numeric(graph_df[y])

    graph_df.columns = ["BONE", "BREAST", "CNS", "HEMATO", "KIDNEY", "COLON", "LIVER", "LUNG","OVARY", "PANCREAS", "SKIN", "STOMACH", "UPPER\nAIRWAY", "URINARY\nTRACT"]
    #graph_df.set_index('Unnamed: 0', inplace=True)
    cmap = mcolors.LinearSegmentedColormap.from_list("n",
                                                     ["#003e92","#ffffff" ,"#ed0000"])
    sns.set(font_scale=0.6)
    g = sns.heatmap(graph_df, cmap=cmap, vmin=vmin, vmax=1, cbar_kws={'label': 'Gain frequency -\nloss frequency'})
    #ax = g.ax_heatmap --> use if you are using sns.clustermap
    ax = g
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    #ax = g.ax_heatmap


    cells = []

    tumors = ["LUNG", "OVARY", "SKIN", "BREAST", "URINARY_TRACT", "PANCREAS", "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "CENTRAL_NERVOUS_SYSTEM","BONE",
              "UPPER_AERODIGESTIVE_TRACT", "KIDNEY", "LARGE_INTESTINE", "LIVER", "STOMACH"]

    # tumors = ["COAD","ESCA",
    #                      "CESC", "HNSC",
    #                      "BLCA",
    #                      "GBM",
    #                      "ACC"]

    tumors.sort()

    # for i in range(0, len(highlight_cells.index)):
    #     #cells.append((highlight_cells.iloc[i, 1], highlight_cells.iloc[i, 0]))
    #     #cells.append((highlight_cells.iloc[i, 3], highlight_cells.iloc[i, 4]))
    #     patch = (tumors.index((highlight_cells.iloc[i, 5]).upper()), highlight_cells.iloc[i, 3])
    #     cells.append(patch)

    #print(cells[1:10])


    # for cell in highlight_cells:
    #
    #     ax.add_patch(Rectangle(cell, 1, 1, fill=False, edgecolor='black', lw=1.5))
    plt.savefig("{}.pdf".format(fname), bbox_inches='tight')
    plt.clf()

def highlight_diff_cells_excl(dfs, tumor_types_order, th):
    """Exclusive version: highlights only those which are sig in one but not other"""
    """Returns a list of tuples of coordinates pf events in the heatmap that show
    different tendencies of occurrence in the two groups"""
    for df in dfs.values():
        ##This ordering is based on alphabetical order
        #print(df.columns)
        df = df.loc[:, tumor_types_order]

    highlight_wgd_plus = []
    highlight_wgd_minus = []
    for i, col in enumerate(dfs["WGD+"].columns):
        for j, (x, y) in enumerate(zip(dfs["WGD+"].loc[:, col], dfs["WGD-"].loc[:, col])):
            ##if event is gained in WGD+ and neutral/lost in WGD- or vice versa,
            ##add its coordinates to the 'highlight' list
            if (x >= th and y < th): ##exclusivity condition
                highlight_wgd_minus.append((i, j))
            elif (x < th and y >= th):
                highlight_wgd_plus.append((i, j))
            else:
                pass
    return highlight_wgd_minus, highlight_wgd_plus


def highlight_diff_cells_incl(dfs, tumor_types_order, th):
    """Inclusive version: highlights all those which are sig in a group - may contain overlaps with other group"""
    """Returns a list of tuples of coordinates pf events in the heatmap that show
    different tendencies of occurrence in the two groups"""
    for df in dfs.values():
        ##This ordering is based on alphabetical order
        df = df.loc[:, tumor_types_order]

    highlight_wgd_plus = []
    highlight_wgd_minus = []
    for i, col in enumerate(dfs["WGD+"].columns):
        for j, (x, y) in enumerate(zip(dfs["WGD+"].loc[:, col], dfs["WGD-"].loc[:, col])):
            ##if event is gained in WGD+ and neutral/lost in WGD- or vice versa,
            ##add its coordinates to the 'highlight' list
            if (y < th): ##inclusivity condition
                highlight_wgd_minus.append((i, j))
            if (x < th):
                highlight_wgd_plus.append((i, j))
            else:
                pass
    return highlight_wgd_minus, highlight_wgd_plus

def find_unified_highlights_excl(dfs, gain_wgd_minus, gain_wgd_plus, loss_wgd_minus, loss_wgd_plus, tumor_types_order, th):
    """Exclusive version: highlights only those which are sig in one but not other"""
    """Returns a list of tuples of coordinates pf events in the heatmap that show
    different tendencies of occurrence in the two groups"""
    for df in dfs.values():
        ##This ordering is based on alphabetical order
        #print(df.columns)
        df = df.loc[:, tumor_types_order]

    highlight_wgd_plus = []
    highlight_wgd_minus = []
    for i, col in enumerate(dfs["WGD+"].columns):

        for j, (x, y) in enumerate(zip(dfs["WGD+"].loc[:, col], dfs["WGD-"].loc[:, col])):
            ##if event is gained in WGD+ and neutral/lost in WGD- or vice versa,
            ##add its coordinates to the 'highlight' list
            if (x > 0): ##exclusivity condition
                if (gain_wgd_plus.iloc[j, i] < th) and (gain_wgd_minus.iloc[j, i] >= th):
                    highlight_wgd_plus.append((i, j))
            if (x < 0):  ##exclusivity condition
                if (loss_wgd_plus.iloc[j, i] < th) and (loss_wgd_minus.iloc[j, i] >= th):
                    highlight_wgd_plus.append((i, j))
            if (y > 0): ##exclusivity condition
                if (gain_wgd_plus.iloc[j, i] >= th) and (gain_wgd_minus.iloc[j, i] < th):
                    highlight_wgd_minus.append((i, j))
            if (y < 0):  ##exclusivity condition
                if (loss_wgd_plus.iloc[j, i] >= th) and (loss_wgd_minus.iloc[j, i] < th):
                    highlight_wgd_minus.append((i, j))
            else:
                pass
    return highlight_wgd_minus, highlight_wgd_plus

# def find_highlights(grp, opp_grp, arm_names, th):
#     "Return the events to be highlighted"
#     grp = pd.concat([arm_names, grp], axis = 1)
#     opp_grp = pd.concat([arm_names, opp_grp], axis=1)
#
#     grp.columns = ["arms", "qvals"]
#     opp_grp.columns = ["arms", "qvals"]
#
#     sig_grp = grp.loc[grp["qvals"] < th, "arms"]
#     sig_opp_grp = opp_grp.loc[opp_grp["qvals"] < th, "arms"]
#
#     return set(sig_grp).difference(set(sig_opp_grp))
def find_unified_highlights_excl_freq(dfs, gain_wgd_minus, gain_wgd_plus, loss_wgd_minus, loss_wgd_plus, tumor_types_order, th):
    """Exclusive version: highlights only those which are sig in one but not other"""
    """Returns a list of tuples of coordinates pf events in the heatmap that show
    different tendencies of occurrence in the two groups"""
    for df in dfs.values():
        ##This ordering is based on alphabetical order
        #print(df.columns)
        df = df.loc[:, tumor_types_order]

    highlight_wgd_plus = []
    highlight_wgd_minus = []
    for i, col in enumerate(dfs["WGD+"].columns):
        print(col)
        for j, (x, y) in enumerate(zip(dfs["WGD+"].loc[:, col], dfs["WGD-"].loc[:, col])):
            ##if event is gained in WGD+ and neutral/lost in WGD- or vice versa,
            ##add its coordinates to the 'highlight' list
            if (x > 0): ##exclusivity condition
                if (gain_wgd_plus.iloc[j, i] < th) and (gain_wgd_minus.iloc[j, i] >= th):
                    if check_frequency((gain_wgd_plus.index)[j], "WGD+", 1, col) >= RECUR: #'1' indicates gain
                        highlight_wgd_plus.append((i, j))
            if (x < 0):  ##exclusivity condition
                if (loss_wgd_plus.iloc[j, i] < th) and (loss_wgd_minus.iloc[j, i] >= th):
                    if check_frequency((loss_wgd_plus.index)[j], "WGD+", -1, col) >= RECUR:  # '1' indicates gain
                        highlight_wgd_plus.append((i, j))
            if (y > 0): ##exclusivity condition
                if (gain_wgd_plus.iloc[j, i] >= th) and (gain_wgd_minus.iloc[j, i] < th):
                    if check_frequency((gain_wgd_minus.index)[j], "WGD-", 1, col) >= RECUR:  # '1' indicates gain
                        highlight_wgd_minus.append((i, j))
            if (y < 0):  ##exclusivity condition
                if (loss_wgd_plus.iloc[j, i] >= th) and (loss_wgd_minus.iloc[j, i] < th):
                    if check_frequency((loss_wgd_minus.index)[j], "WGD-", -1, col) >= RECUR:  # '1' indicates gain
                        highlight_wgd_minus.append((i, j))
            else:
                pass
    return highlight_wgd_minus, highlight_wgd_plus

def find_unified_highlights_incl_freq(dfs, gain_wgd_minus, gain_wgd_plus, loss_wgd_minus, loss_wgd_plus, tumor_types_order, th):
    """Exclusive version: highlights only those which are sig in one but not other"""
    """Returns a list of tuples of coordinates pf events in the heatmap that show
    different tendencies of occurrence in the two groups"""
    for df in dfs.values():
        ##This ordering is based on alphabetical order
        #print(df.columns)
        df = df.loc[:, tumor_types_order]

    highlight_wgd_plus = []
    highlight_wgd_minus = []
    for i, col in enumerate(dfs["WGD+"].columns):
        print(col)
        for j, (x, y) in enumerate(zip(dfs["WGD+"].loc[:, col], dfs["WGD-"].loc[:, col])):
            ##if event is gained in WGD+ and neutral/lost in WGD- or vice versa,
            ##add its coordinates to the 'highlight' list
            if (x > 0): ##exclusivity condition
                if (gain_wgd_plus.iloc[j, i] < th):
                    if check_frequency((gain_wgd_plus.index)[j], "WGD+", 1, col) >= RECUR: #'1' indicates gain
                        highlight_wgd_plus.append((i, j))
            if (x < 0):  ##exclusivity condition
                if (loss_wgd_plus.iloc[j, i] < th):
                    if check_frequency((loss_wgd_plus.index)[j], "WGD+", -1, col) >= RECUR:  # '1' indicates gain
                        highlight_wgd_plus.append((i, j))
            if (y > 0): ##exclusivity condition
                if (gain_wgd_minus.iloc[j, i] < th):
                    if check_frequency((gain_wgd_minus.index)[j], "WGD-", 1, col) >= RECUR:  # '1' indicates gain
                        highlight_wgd_minus.append((i, j))
            if (y < 0):  ##exclusivity condition
                if (loss_wgd_minus.iloc[j, i] < th):
                    if check_frequency((loss_wgd_minus.index)[j], "WGD-", -1, col) >= RECUR:  # '1' indicates gain
                        highlight_wgd_minus.append((i, j))
            else:
                pass
    return highlight_wgd_minus, highlight_wgd_plus

def check_frequency(arm, group, anp, tumor_type):
    "Checks if frequency of event in Taylor data is at least RECUR"
    taylor = pd.read_csv(os.path.join(TAYLOR_DIR, "Type_wise_df_no-ke97.tsv"), sep = "\t", header = 0)
    #print(tumor_type)
    #taylor_tumor_type = taylor[tumor_type in taylor["CCLE_ID"]]
    taylor_tumor_type = taylor[taylor['CCLE_ID'].str.contains(tumor_type)]
    if group == "WGD+":
        taylor_df = taylor_tumor_type[taylor_tumor_type["ploidy"] > 3]
    else:
        taylor_df = taylor_tumor_type[taylor_tumor_type["ploidy"] < 2.5]

    arm_col_num = pd.Index(taylor_df.columns).get_loc(arm)
    count_event = len(taylor_df[taylor_df.iloc[:,arm_col_num] == anp])
    n_sample = len(taylor_df.iloc[:,arm_col_num])

    if tumor_type == "OVARY":
        print(arm, count_event/float(n_sample))

    return count_event/float(n_sample)


if __name__ == "__main__":
    #DIR = "C:\\Users\\kiwii\\Documents\\Kavya\\Uri\\Work\\Tasks\\Gistic2\\kavya_cancer_1st.tar\\kavya_cancer_1st\\kavya\\"
    #DIR = "C:\\Users\\kiwii\\Documents\\Kavya\\Uri\\Work\\Tasks\\Gistic2\\CCLE\\broad_peaks_ccle\\"
    DIR = "C:\\Users\\kiwii\\OneDrive\\Documents\\Kavya\\Uri\\Work\\Cancer Research followup\\GISTIC analysis\\kavya_ccle_ci-90.tar\\kavya_ccle_ci-90\\kavya_ccle\\"
    global HGLT_FLAG, THRESH, UNIFIED, RECUR, TAYLOR_DIR
    HGLT_FLAG = 1 #1 = excl
    THRESH = 0.05
    UNIFIED = 1
    RECUR = 0.1
    TAYLOR_DIR = "C:\\Users\\kiwii\\OneDrive\\Documents\\Kavya\\Uri\\Work\\Tasks\\CCLE analysis\\Prelim\\"
    read_data(DIR)
