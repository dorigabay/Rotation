import matplotlib.pyplot as plt
import glob
import re
import numpy
import pandas as pd
import numpy as np
import os
import seaborn as sns
import shutil
from sklearn.neighbors import KernelDensity

# tips = sns.load_dataset("tips")
# print(tips)
# X= np.array(tips["total_bill"]).reshape(1, -1)
# kernel = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(X)
# hist = np.histogram(X,bins=100)
# print(hist)
# print(kernel
#       )
def flaten_data(data_frame):
    df = data_frame.drop([data_frame.columns[0], "time"], 1)
    RG_flatten = np.array([])
    resi_flatten = np.array([])
    FRET_flatten = np.array([])
    new_df = pd.DataFrame()
    for m in df:
        if "Rg_repeat" in m:
            RG_flatten = np.append(RG_flatten, df[m])
        elif "resi_distance" in m:
            resi_flatten = np.append(resi_flatten, df[m])
        elif "FRET" in m:
            FRET_flatten = np.append(FRET_flatten, df[m])
    new_df["Rg ($\AA$)"] = RG_flatten
    new_df["Residues Lables Distances ($\AA$)"] = resi_flatten
    new_df["smFRET"] = FRET_flatten
    return new_df

# def distribution_figure_2dimensions(dict_data,figure_type, config=None,title=None):
#     figure_types = {"Rg":"Rg ($\AA$)","residues":"Residues Lables Distances ($\AA$)","FRET":"smFRET"}
#     figure_titles = {"Rg":"Rg distribution","residues":"end to end residues distances","FRET":"smFRET"}
#     colors = {"segA":"red","segB":"green","segC":"blue"}
#     fig, ax = plt.subplots(1, 1, figsize=(16, 16))
#     print(dict_data)
#     for dim2,stylee in zip(dict_data,["-","-."]):
#         for name in dict_data[dim2].keys():
#             hist = np.histogram(dict_data[dim2][name][figure_types[figure_type]],bins=80,density=True,range=(0,80))
#             ax.plot(hist[1][:-1],hist[0],stylee,label=dim2+"_"+name,color = colors[name],lw=2)
#             ax.legend(prop={"size": 12})
#             ax.set_xlim((10,30))
#             if config != None:
#                 ax.set_xlim(config[dim2]["xlim"][figure_type])
#                 ax.set_ylim(config[dim2]["ylim"][figure_type])
#     fig.suptitle(title + " - "+figure_titles[figure_type], fontsize=16)
#     fig.tight_layout()
#     return fig, figure_type


def distribution_figure_2dimensions(dict_data,figure_type, config=None,title=None):
    figure_types = {"Rg":"Rg ($\AA$)","residues":"Residues Lables Distances ($\AA$)","FRET":"smFRET"}
    figure_titles = {"Rg":"Rg distribution","residues":"end to end residues distances","FRET":"smFRET"}
    dimension2 = [x for x in dict_data.keys()]
    fig, ax = plt.subplots(len(dimension2), 1, figsize=(16, 16))
    for col, dim2 in zip(ax, dict_data):
        col.set_title(dim2, rotation=0, size='large')
        for name in dict_data[dim2].keys():
            sns.kdeplot(ax=col, data=dict_data[dim2][name], label=name, linewidth=2.5, x=figure_types[figure_type])
            col.legend(prop={"size": 12})
            if config != None:
                col.set_xlim(config[dim2]["xlim"][figure_type])
                col.set_ylim(config[dim2]["ylim"][figure_type])
    fig.suptitle(title + " - "+figure_titles[figure_type], fontsize=16)
    fig.tight_layout()
    return fig, figure_type


# def residue_figure_2dimensions(dict_data, config=None,title=None):
#     dimension2 = [x for x in dict_data.keys()]
#     fig, ax = plt.subplots(len(dimension2), 1, figsize=(16, 16))
#     for col, dim2 in zip(ax, dict_data):
#         col.set_title(dim2, rotation=0, size='large')
#         for name in dict_data[dim2].keys():
#             sns.kdeplot(ax=col, data=dict_data[dim2][name], label=name, linewidth=2.5,x="Residues Lables Distances ($\AA$)")
#             col.legend(prop={"size": 12})
#             if config != None:
#                 col.set_xlim(config[dim2]["xlim"]["Rg"])
#                 col.set_ylim(config[dim2]["ylim"]["Rg"])
#     fig.suptitle(title + " - end to end residues distances", fontsize=16)
#     fig.tight_layout()
#     return fig, "Residue"


def Rg_figure_3dimensions(dict_data, config=None,title=None):
    dimension3 = [x for x in dict_data.keys()]
    dimension2 = dict_data[dimension3[0]].keys()
    fig, ax = plt.subplots(len(dimension3), len(dimension2), figsize=(16, 16))
    for row,dim3 in zip(ax,dict_data):
        row[0].set_ylabel(dim3, size='large')
        for col,dim2 in zip(row,dict_data[dim3]):
            col.set_title(dim2,rotation=0, size='large')
            for name in dict_data[dim3][dim2].keys():
                sns.kdeplot(ax=col, data=dict_data[dim3][dim2][name], label=name, linewidth=2.5, x="Rg ($\AA$)")
                col.legend(prop={"size": 12})
                if config != None:
                    col.set_xlim(config[dim2]["xlim"]["Rg"])
                    col.set_ylim(config[dim2]["ylim"]["Rg"])
    fig.suptitle(title+" - Rg distribution", fontsize=16)
    fig.tight_layout()
    return fig,"Rg"


def residue_figure_3dimensions(dict_data, config=None,title=None):
    dimension3 = [x for x in dict_data.keys()]
    dimension2 = dict_data[dimension3[0]].keys()
    fig, ax = plt.subplots(len(dimension3), len(dimension2), figsize=(16, 16))
    for row, dim3 in zip(ax, dict_data):
        row[0].set_ylabel(dim3, size='large')
        for col, dim2 in zip(row, dict_data[dim3]):
            col.set_title(dim2,rotation=0, size='large')
            for name in dict_data[dim3][dim2].keys():
                sns.kdeplot(ax=col, data=dict_data[dim3][dim2][name], label=name, linewidth=2.5,x="Residues Lables Distances ($\AA$)")
                col.set_title(dim2, size=12)
                col.legend(prop={"size": 12})
                if config != None:
                    col.set_xlim(config[dim2]["xlim"]["residue"])
                    col.set_ylim(config[dim2]["ylim"]["residue"])
    fig.suptitle(title+" - end to end residues distances", fontsize=16)
    fig.tight_layout()
    return fig,"Residue"


def reorganize_data(on_the_plot,directory):
    maindir = re.search(".*/(.*/)",directory).group(1)
    output_dir = "/home_e/dor/PycharmProjects/pythonProject/"+maindir
    directories_toiter = [directory]
    while directories_toiter:
        dir = directories_toiter.pop()
        dir_name = re.search(".*/(.*)/",dir).group(1)
        if on_the_plot in dir_name:
            directories_toiter1 = [dir]
            while directories_toiter1:
                dir1 = directories_toiter1.pop()
                if glob.glob(dir1+"*.csv"):
                    dir_out = dir1.replace(dir_name+"/","").replace(directory,output_dir)
                    os.makedirs(dir_out,exist_ok=True)
                    for file in glob.glob(dir1+"*.csv"):
                        shutil.copyfile(file,dir_out+dir_name+".csv")
                else:
                    for subdir1 in [x for x in os.walk(dir1)][0][1]:
                        directories_toiter1.append(dir1 + subdir1 + "/")
        else:
            for subdir in [x for x in os.walk(dir)][0][1]:
                directories_toiter.append(dir + subdir + "/")
    return "/home_e/dor/PycharmProjects/pythonProject/"+maindir


def putData_inDict_4dimensions(directory,on_the_plot):
    if on_the_plot != None:
        directory_organized = reorganize_data(on_the_plot, directory)
    else: directory_organized = directory
    directories_toiter = [directory_organized]
    complete_data = {}
    dim4,dim3,dim2,dim1 = "","","",""
    while directories_toiter:
        dir = directories_toiter.pop()
        if glob.glob(dir + "*.csv"):
            dim4,dim3,dim2 = re.search(directory_organized+"(.*)/",dir).group(1).split("/")
            dict_data = {}
            for file in glob.glob(dir + "*.csv"):
                dim1 = re.search(".*/(.*).csv", file).group(1)
                faltten_df = flaten_data(pd.read_csv(file))
                dict_data[dim1] = faltten_df
            if dim4 not in complete_data:
                complete_data[dim4] = {dim3:{dim2:dict_data}}
            elif dim3 not in complete_data[dim4]:
                complete_data[dim4][dim3] = {dim2:dict_data}
            else:
                complete_data[dim4][dim3][dim2] = dict_data
        else:
            for subdir in [x for x in os.walk(dir)][0][1]:
                directories_toiter.append(dir + subdir + "/")
    if on_the_plot != None:
        shutil.rmtree(directory_organized)
    return complete_data


def putData_inDict_2dimensions(directory,on_the_plot=None):
    if on_the_plot!=None:
        directory_organized = reorganize_data(on_the_plot, directory)
        print("done")
    else: directory_organized = directory
    directories_toiter = [directory_organized]
    complete_data = {}
    dim2, dim1 = "", ""
    while directories_toiter:
        dir = directories_toiter.pop()
        if glob.glob(dir + "*.csv"):
            dim2 = re.search(directory_organized+"(.*)/",dir).group(1).split("/")[0]
            dict_data = {}
            for file in glob.glob(dir + "*.csv"):
                dim1 = re.search(".*/(.*).csv", file).group(1)
                faltten_df = flaten_data(pd.read_csv(file))
                dict_data[dim1] = faltten_df
            if dim2 not in complete_data:
                complete_data[dim2] = dict_data
        else:
            for subdir in [x for x in os.walk(dir)][0][1]:
                directories_toiter.append(dir + subdir + "/")
    if on_the_plot != None:
        shutil.rmtree(directory_organized)
    return complete_data


def figures_4dimensions(directory,figure_type,on_the_plot=None,configuration=None,title=None):
    ouput_dir = directory+"graphs/"
    complete_data = putData_inDict_4dimensions(directory, on_the_plot)
    for dimesion4 in complete_data:
        if configuration != None:
            config = configuration[dimesion4]
        else:
            config = None
        for type in figure_type:
            figure,analysis_type = type(complete_data[dimesion4], config=config,title=title+" - "+dimesion4)
            os.makedirs(ouput_dir, exist_ok=True)
            figure.savefig(ouput_dir + dimesion4 +analysis_type+"_graph.png")
        print("saved")


def figures_2dimensions(directory,figure_type,on_the_plot=None,configuration=None,title=None):
    ouput_dir = directory+"graphs/"
    complete_data = putData_inDict_2dimensions(directory, on_the_plot)
    for type in figure_type:
        # if type == "smFRET":

        figure,analysis_type = distribution_figure_2dimensions(complete_data,type, config=configuration,title=title)
        os.makedirs(ouput_dir, exist_ok=True)
        figure.savefig(ouput_dir +analysis_type+"_graph.png")
        print("saved")


# def put_allData_together(directory):
#     directories_toiter = [directory]
#     while dir_to_iter:
#         dir = directories_toiter.pop()
#         if glob.glob(dir+"*.csv"):
#             if
#             titles = re.search(directory_organized + "(.*)/", dir).group(1).split("/")
#
#

configuration = {
    "temp03" : {"eps02":{"xlim":{"Rg":(0,60),"residue":(0,300)},
                       "ylim":{"Rg":(0,0.25),"residue":(0,0.03)}},
              "eps04":{"xlim":{"Rg":(8,30),"residue":(0,100)},
                       "ylim":{"Rg":(0,1.3),"residue":(0,0.05)}},
              "eps03":{"xlim":{"Rg":(8,40),"residue":(0,150)},
                       "ylim":{"Rg":(0,0.8),"residue":(0,0.05)}}},

    "temp04" : {"eps02":{"xlim":{"Rg":(0,60),"residue":(0,300)},
                       "ylim":{"Rg":(0,0.25),"residue":(0,0.03)}},
              "eps04":{"xlim":{"Rg":(8,30),"residue":(0,100)},
                       "ylim":{"Rg":(0,1.3),"residue":(0,0.05)}},
              "eps03":{"xlim":{"Rg":(8,40),"residue":(0,150)},
                       "ylim":{"Rg":(0,0.8),"residue":(0,0.05)}}},

    "temp05": {"eps02": {"xlim": {"Rg": (0, 60), "residue": (0, 300)},
                         "ylim": {"Rg": (0, 0.25), "residue": (0, 0.03)}},
               "eps04": {"xlim": {"Rg": (8, 30), "residue": (0, 100)},
                         "ylim": {"Rg": (0, 1.3), "residue": (0, 0.05)}},
               "eps03": {"xlim": {"Rg": (8, 40), "residue": (0, 150)},
                         "ylim": {"Rg": (0, 0.8), "residue": (0, 0.05)}}},
}

# figures_4dimensions("/home_e/dor/PycharmProjects/pythonProject/graphs/coarse/diff_temp_withAromatic_constant_bondAngles/","seg",[Rg_figure_3dimensions,residue_figure_3dimensions],configuration,title="With aromatic residues")
# figures_4dimensions("/home_e/dor/PycharmProjects/pythonProject/graphs/coarse/diff_temp_noAromatic_constantBondAngles/","seg",[Rg_figure_3dimensions,residue_figure_3dimensions],configuration,title="Without aromatic residues")
# figures_4dimensions("/home_e/dor/PycharmProjects/pythonProject/graphs/coarse/diff_temp_noAromatic_constantBondAngles/","eps",[Rg_figure_3dimensions,residue_figure_3dimensions],title="Without aromatic residues")
# figures_4dimensions("/home_e/dor/PycharmProjects/pythonProject/graphs/coarse/diff_temp_withAromatic_constant_bondAngles/","eps",[Rg_figure_3dimensions,residue_figure_3dimensions],title="Without aromatic residues")

# figures_2dimensions("/home_e/dor/PycharmProjects/pythonProject/graphs/coarse/withVSwithout_proline/","seg",[Rg_figure_2dimensions,residue_figure_2dimensions],title="The dihedral effect where there is only repulsion, bonds, and hookean")
# figures_2dimensions("/home_e/dor/PycharmProjects/pythonProject/graphs/atomistic/charmm_ff/histidine_nonprotonated/",
#                     ["Rg","residues"],on_the_plot="histidine_nonprotonated" ,title="Gromacs simulations under different salt conditions")
figures_2dimensions("/home_e/dor/PycharmProjects/pythonProject/graphs/new_experiments/nothing_butproline_strongAngles/",["Rg","residues"],on_the_plot="seg",title="babybe")