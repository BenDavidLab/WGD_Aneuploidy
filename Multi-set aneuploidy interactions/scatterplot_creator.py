import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches,markers,axes
from pylab import MaxNLocator
import colors_code as cl


def legend_creator(dict):
    tumor_list = list(dict.keys())
    color_list = list(dict.values())
    patches_list= []
    for color in color_list:
        patch = patches.Circle((0,0),0.1,color = color)
        patches_list.append(patch)
    return patches_list,tumor_list


def get_color(intersection,df,dict):
    type = df[df["Intersections"]==intersection]["tumor_type"].values[0]
    color = dict.get(type)
    return color



TUMOR_LIST = ["ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LGG", "LIHC", "LUAD", "LUSC",
              "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC"]
TUMOR_LIST1 = ["BRCA"]
color_db = pd.read_csv(r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\colors_rgb.csv")
chartreuse3  = '#%02x%02x%02x' %(102,205,0)
mediumpurple3 = '#%02x%02x%02x' %(137,104,205)
violetred1 = '#%02x%02x%02x' %(255,62,150)
turquoise3 = '#%02x%02x%02x' %(0,197,205)
blue4 = '#%02x%02x%02x' %(0,0,139)
slateblue2 = '#%02x%02x%02x' %(122,103,238)
gray75 = '#%02x%02x%02x' %(191,191,191)
firebrick4 = '#%02x%02x%02x' %(139,26,26)
chocolate1 = '#%02x%02x%02x' %(255,127,36)
gray49 = '#%02x%02x%02x' %(125,125,125)
springgreen4 = '#%02x%02x%02x' %(0,139,69)
bisque = '#%02x%02x%02x' %(255,228,196)
lavenderblush2 = '#%02x%02x%02x' %(238,224,229)
cyan4 = '#%02x%02x%02x' %(0,139,139)
purple3 = '#%02x%02x%02x' %(125,38,205)
gray32 = '#%02x%02x%02x' %(82,82,82)
deepskyblue1 = '#%02x%02x%02x' %(0,191,255)
blue = '#%02x%02x%02x' %(0,0,255)
gold = '#%02x%02x%02x' %(255,215,0)
gray14 = '#%02x%02x%02x' %(36,36,36)
lightslateblue = '#%02x%02x%02x' %(132,112,255)
salmon = '#%02x%02x%02x' %(250,128,114)


###scatter plot###
plot_data = pd.read_csv(r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\finaldata\all_tumor_types_interesting_events_1.1_10.csv")
for i in [2,3,4,5]:
    plot_data_filtered = plot_data[plot_data["Degree_wgd(+)"]== i]
    paper_df = plot_data_filtered[plot_data_filtered["paper"]=="Yes"]
    powerset_df = plot_data_filtered[(plot_data_filtered["Powerset_relevant"]=="Yes") & (plot_data_filtered["paper"]=="No") ]
    interesting_df = plot_data_filtered[(plot_data_filtered["Powerset_relevant"]=="No") & (plot_data_filtered["paper"]=="No")]
    x_paper = paper_df["-log(q)_wgd(-)"]
    y_paper=  paper_df["-log(q)_wgd(+)"]
    x_powerset =powerset_df["-log(q)_wgd(-)"]
    y_powerset = powerset_df["-log(q)_wgd(+)"]
    x_interesting = interesting_df["-log(q)_wgd(-)"]
    y_interesting = interesting_df["-log(q)_wgd(+)"]
    color_dict = {"ACC":chartreuse3,"BLCA":mediumpurple3,"BRCA":violetred1,"CESC":turquoise3,"COAD":blue4,"ESCA":slateblue2,"GBM":gray75,"HNSC":firebrick4,"KIRC":chocolate1,"LGG":gray49,"LIHC":springgreen4,"LUAD":bisque,"LUSC":lavenderblush2,"OV":cyan4,"PAAD":purple3,"PCPG":gray32,"PRAD":deepskyblue1,"READ":blue,"SARC":gold,"SKCM":gray14,"STAD":lightslateblue,"UCEC":salmon}
    color_list = [get_color(x,paper_df,color_dict) for x in paper_df["Intersections"]]
    scatterplot = plt.scatter (x_paper,y_paper,c =color_list,s=15,edgecolors="black",marker = "X",alpha=0.8)
    color_list = [get_color(x,powerset_df,color_dict) for x in powerset_df["Intersections"]]
    plt.scatter(x_powerset, y_powerset, c=color_list, s=15, edgecolors="black", marker="^",alpha=0.8)
    color_list = [get_color(x,interesting_df,color_dict) for x in interesting_df["Intersections"]]
    plt.scatter(x_interesting, y_interesting, c=color_list, s=15, edgecolors="black", marker=".",alpha=0.8)
    handles,lables = legend_creator(color_dict)
    plt.xlabel("-log(q)_wgd(-)")
    plt.ylabel("-log(q)_wgd(+)")
    plt.title("Interesting Events - all tumor types n="+str(i))
    plt.ylim(-0.05,2)
    color_legend = plt.legend(handles,lables,loc = "upper right" ,ncol=2,borderaxespad=0.,bbox_to_anchor=(0.99, 0.99), prop={'size': 8})
    plt.legend(["Showen in paper","Powerset significant","interesting events"],loc = "center right",borderaxespad=0,bbox_to_anchor=(0.99,0.4),prop={'size': 8})
    plt.gca().add_artist(color_legend)

    scatter_pdf_path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\finaldata\scatter_plots_1.1_10+2%/interesting_events_all_tumor_types_1.1_10+2%_n"+str(i)+".pdf"
    plt.savefig(scatter_pdf_path)
    plt.clf()


    events_list = plot_data["Intersections"].to_list()
    x_dict = {"2":0,"3":0,"4":0,"5":0}
    x_values = {"2":0,"3":0,"4":0,"5":0}
    for tumor in TUMOR_LIST:
        bar_data = plot_data[plot_data["tumor_type"]== tumor]
        #x_dict = {"2":len(bar_data[bar_data["Degree_wgd(+)"]==2]),"3":(len(bar_data[bar_data["Degree_wgd(+)"]==3])),"4":(len(bar_data[bar_data["Degree_wgd(+)"]==4])),"5":(len(bar_data[bar_data["Degree_wgd(+)"]==5]))}
        bar1 = {"2":len(bar_data[bar_data["paper"]=="Yes"]),"3":0,"4":0,"5":0}
        bar2 = {"2":len(bar_data[(bar_data["Powerset_relevant"]=="Yes") & (bar_data["Degree_wgd(+)"]==2)])-len(bar_data[bar_data["paper"]=="Yes"]),"3":len(bar_data[(bar_data["Degree_wgd(+)"]==3) &(bar_data["Powerset_relevant"]=="Yes")]),"4":len(bar_data[(bar_data["Degree_wgd(+)"]==4) &(bar_data["Powerset_relevant"]=="Yes")]),"5":len(bar_data[(bar_data["Degree_wgd(+)"]==5) &(bar_data["Powerset_relevant"]=="Yes")])}
        #bar3 = {"2":x_dict.get("2")-(bar1.get("2")+bar2.get("2")),"3":x_dict.get("3")-(bar1.get("3")+bar2.get("3")),"4":x_dict.get("4")-(bar1.get("4")+bar2.get("4")),"5":x_dict.get("5")-(bar1.get("5")+bar2.get("5"))}
        #barplot1 = plt.bar(list(bar1.keys()),list(bar1.values()),bottom=list(bar2.values()),width=0.5,color = "green")
        #bars= np.add(list(bar1.values()),list(bar2.values())).tolist()
        #barplot3  =plt.bar(list(bar3.keys()),list(bar3.values()),bottom= bars,width=0.5,color = "grey")

        barplot2= plt.bar(list(bar2.keys()),list(bar2.values()),bottom=list(x_values.values()),width=0.5,edgecolor = "black",color = color_dict.get(tumor),alpha = 0.75)
        for i in [2,3,4,5]:
            x_values[str(i)] += bar2[str(i)]
    for tumor in TUMOR_LIST:
        bar_data = plot_data[plot_data["tumor_type"]== tumor]
        bar1 = {"2": len(bar_data[bar_data["paper"] == "Yes"]), "3": 0, "4": 0, "5": 0}
        barplot2 = plt.bar(list(bar1.keys()), list(bar1.values()), bottom=list(x_values.values()), width=0.5,color=color_dict.get(tumor),edgecolor = "black",alpha = 0.75,hatch = "/")
        for i in [2,3,4,5]:
            x_values[str(i)] += bar1[str(i)]
    color_legend = plt.legend(handles,lables,loc = "upper right" ,ncol=2,borderaxespad=0,bbox_to_anchor=(0.99, 0.99), prop={'size': 8})
    plt.legend([patches.Rectangle((0,0),1,1,facecolor = "grey",hatch = "///",edgecolor="black")],["Showen in paper"],borderaxespad=0,bbox_to_anchor=(0.99,0.45),prop={'size': 9},loc = "center right")
    plt.gca().add_artist(color_legend)
    plt.xlabel("Number of arms events")
    plt.ylabel("Number of events ")
    plt.title("Events per arm")
    plt.ylim(0,max(list(x_values.values()))+1)
        #ticker_spacing = 1
        #plt.gca().yaxis.set_major_locator(plt.MultipleLocator(1))
    bar_pdf_path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\finaldata\bar_plots_1.1_10+2%/all_types_1.1_10+2%.pdf"
    plt.savefig(bar_pdf_path)
    plt.clf()

bar1 = {"2":len(plot_data[plot_data["paper"]=="Yes"]),"3":0,"4":0,"5":0}
bar2 = {"2": len(plot_data[(plot_data["Powerset_relevant"] == "Yes") & (plot_data["Degree_wgd(+)"] == 2)]) - len(plot_data[plot_data["paper"] == "Yes"]),"3": len(plot_data[(plot_data["Degree_wgd(+)"] == 3) & (plot_data["Powerset_relevant"] == "Yes")]),"4": len(plot_data[(plot_data["Degree_wgd(+)"] == 4) & (plot_data["Powerset_relevant"] == "Yes")]),"5": len(plot_data[(plot_data["Degree_wgd(+)"] == 5) & (plot_data["Powerset_relevant"] == "Yes")])}
barplot2 = plt.bar(list(bar2.keys()),list(bar2.values()),width=0.5,color = "grey")
barplot1 = plt.bar(list(bar1.keys()),list(bar1.values()),bottom=list(bar2.values()),width=0.5,color = "red")
plt.legend([(patches.Rectangle((0,0),1,1,facecolor="grey",edgecolor="grey")),(patches.Rectangle((0,0),1,1,facecolor="red",edgecolor="red"))],["Powerset significant events","Showen in paper"],loc = "upper right" ,ncol=2,borderaxespad=0,bbox_to_anchor=(0.99, 0.99), prop={'size': 8})
plt.xlabel("Number of arms involved in each event(N)")
plt.ylabel("Number of events ")
plt.title("Sum of events per N")
plt.ylim(0,max(list(x_values.values()))+1)
bar_pdf_path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Phd\kavya project\superexacttest\finaldata\bar_plots_1.1_10+2%/all_types_second_version_1.1_10+2%.pdf"
plt.savefig(bar_pdf_path)
plt.clf()