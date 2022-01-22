import numpy as np
import pandas as pd
import re
from itertools import chain, combinations


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) 1,2,3)"
    s = list(iterable)
    iterable_list = chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    output = []
    for x in iterable_list:
        if len(x) >=2:
            output.append(x)
    return output



def is_relevant_wgd_plus(Intersection,df):
    cell_value = Intersection.split(" & ")
    power_set = powerset(cell_value)
    main_wgd_plus_qval = float(df[df["Intersections"] == Intersection]["Q.value_wgd(+)"].values[0])
    qval_list = []
    for set in power_set:
        right_format = " & ".join(set)
        qval = float(df[df["Intersections"]== right_format]["Q.value_wgd(+)"].values[0])
        inner_list = [right_format, qval]
        if main_wgd_plus_qval <= qval:
            relevance = True
            qval_list.append(inner_list)
        else:
            qval_list = ["Not Relevant"]
            relevance = False
            return qval_list,relevance
    return qval_list,relevance

def is_relevant_wgd_minus(Intersection,df):
    cell_value = Intersection.split(" & ")
    power_set = powerset(cell_value)
    main_wgd_minus_qval = float(df[df["Intersections"] == Intersection]["Q.value_wgd(-)"].values[0])
    qval_list= []
    for set in power_set:
        right_format = " & ".join(set)
        qval = float(df[df["Intersections"]== right_format]["Q.value_wgd(-)"].values[0])
        inner_list = [right_format, qval]
        qval_list.append(inner_list)
        if main_wgd_minus_qval <= qval:
            relevance = True
        else:
            qval_list = ["Not Relevant"]
            relevance = False
            return qval_list,relevance
    return qval_list,relevance






if __name__ == "__main__":
    TUMOR_LIST = ["ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LGG", "LIHC", "LUAD", "LUSC",
                  "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]
    TUMOR_LIST1 = ["BRCA"]
    FOLDER_PATH = r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\for_graphs\both_wgd"
    for tumor in TUMOR_LIST:
        raw_data  = pd.read_csv(FOLDER_PATH+"/"+tumor+"both_wgd_combined.csv")
        interesting_events = raw_data[raw_data["interesting_events" ]== "Very_Interesting"]["Intersections"].to_list()
        intersection_list = raw_data["Intersections"].to_list()
        powerset_relevant_list = []
        qval_list = []
        for event in intersection_list:
            if event in interesting_events:
                if raw_data[raw_data["Intersections"] == event]["Q.value_wgd(+)"].values[0] <= 0.05:
                    wgd_qval_list,wgd_powerset_relevant= is_relevant_wgd_plus(event,raw_data)
                else:
                    wgd_qval_list,wgd_powerset_relevant= is_relevant_wgd_minus(event,raw_data)
                if wgd_powerset_relevant==True:
                    powerset_relevant = "Yes"
                    qval_list.append(wgd_qval_list)


                else:
                    powerset_relevant = "No"
                    qval_list.append("Not Relevant")


            else:
                powerset_relevant = "No"
                qval_list.append("Not Relevant")

            powerset_relevant_list.append(powerset_relevant)
        print(tumor)
        raw_data["Powerset_relevant"] = powerset_relevant_list
        raw_data["wgd_Qval_list"] = qval_list
        final_folder = r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\finaldata\new_power_set"
        final_path = final_folder+"/"+tumor+"_final_data_fe0.9_qo.95_2%.csv"
        raw_data.to_csv(final_path,index=False)





