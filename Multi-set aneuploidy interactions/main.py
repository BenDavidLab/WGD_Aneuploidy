import math

import numpy
import pandas
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import patches
from scipy import stats

def combine_file(WGD_plus,WGD_minus,WGD_plus_count,WGD_minus_count,folder_path):
    interesting_list =[]
    qval = WGD_plus["Q.value"].apply(func=lambda x: -1*math.log10(x))
    WGD_plus["-log(q)"] = qval
    qval = WGD_minus["Q.value"].apply(func=lambda x: -1 * math.log10(x))
    WGD_minus["-log(q)"] = qval
    combined_files = WGD_plus.merge(WGD_minus,how= "inner",on="Intersections",suffixes=("_wgd(+)","_wgd(-)"))
    combined_files = combined_files[["Intersections","Observed.Overlap_wgd(+)","Observed.Overlap_wgd(-)","Degree_wgd(+)","FE_wgd(+)","Q.value_wgd(+)","-log(q)_wgd(+)","FE_wgd(-)","Q.value_wgd(-)","-log(q)_wgd(-)"]].copy()
    intersection_list = combined_files["Intersections"]
    for intersection in intersection_list:
        logy = float(combined_files[combined_files["Intersections"]==intersection]["-log(q)_wgd(+)"].values[0])
        logx = float(combined_files[combined_files["Intersections"]==intersection]["-log(q)_wgd(-)"].values[0])
        fex = float(combined_files[combined_files["Intersections"]==intersection]["FE_wgd(-)"].values[0])
        fey= float(combined_files[combined_files["Intersections"]==intersection]["FE_wgd(+)"].values[0])
        obsy = combined_files[combined_files["Intersections"]==intersection]["Observed.Overlap_wgd(+)"].values[0]
        obsx = combined_files[combined_files["Intersections"]==intersection]["Observed.Overlap_wgd(-)"].values[0]
        obsx_perc = (WGD_minus_count/100)*2
        obsy_prec = (WGD_plus_count/100)*2
        calc = (((logx>=1.3) and (logy <=0.0222763947111523)) or ((logx<=0.0222763947111523) and (logy>=1.3)))
        if ((fex>=1 and (obsx <=2 or obsx<= obsx_perc)) or (fey>=1 and (obsy<=2 or obsy<=obsy_prec))):
            interesting = "Not_Interesting"
        elif calc and(((fex>1.1) and (fey<0.9)) or ((fey>1.1) and (fex<0.9))):
            interesting = "Very_Interesting"
        elif calc:
            interesting = "Interesting"
        else:
            interesting = "Not_Interesting"
        interesting_list.append(interesting)
    combined_files["interesting_events"] = interesting_list
    combined_files.to_csv(folder_path, index=False)
    return combined_files



tumor_list = ["ACC","BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","UCEC"]
tumor_list1 = ["BRCA_DUMMY"]
for tumor in tumor_list:
    plus_path= r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\filtered_reccuring_arm_calls\super_exact_results_q_val_all_events_n_2/" +tumor+r"_WGD+ _super_results.csv"
    minus_path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\filtered_reccuring_arm_calls\super_exact_results_q_val_all_events_n_2/" +tumor+r"_WGD- _super_results.csv"
    samples_path = r"D:\user\Desktop\kavya_projct\samples_per_type.csv"
    samples_df = pd.read_csv(samples_path)
    WGD_plus= pd.read_csv(plus_path)
    WGD_plus_dummy = WGD_plus[WGD_plus["Degree"]==2]
    WGD_minus = pd.read_csv(minus_path)
    WGD_minus_dummy = WGD_minus[WGD_minus["Degree"]==2]
    WGD_plus_count = samples_df[samples_df["Tumor type"]==tumor]["WGD_plus"].values[0]
    WGD_minus_count =samples_df[samples_df["Tumor type"]==tumor]["WGD_minus"].values[0]
    folder_path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\for_graphs\all_events_n_2/"+tumor+"_all_events_combined.csv"
    combined_data = combine_file(WGD_plus,WGD_minus,WGD_plus_count,WGD_minus_count,folder_path)
    print(tumor)