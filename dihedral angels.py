import MDAnalysis as MDA
import MDAnalysis.analysis.dihedrals as di
import MDAnalysis.lib.distances as dis
import numpy as np
import matplotlib.pyplot as plt
import math
import glob
import pandas as pd
import re
import seaborn as sns
import os


def calc_bonds(config_file,universe):
    config = open(config_file,"r").read()
    rows = config.split("\n")
    df = pd.DataFrame()
    names = {}
    for row in rows:
        all_angels = np.zeros(len(universe.trajectory))
        if "BondGroup" in row:
            location = re.search(r"BondGroup: ([A-Z]{3}), Location: (\d.*)", row).group(2)
            group_name = re.search(r"BondGroup: ([A-Z]{3}), Location: (\d.*)", row).group(1)
            names[location]=group_name
            atom_group = universe.select_atoms("name CA and resid "+str(location))
            idx=0
            for steps in universe.trajectory:
                pos = atom_group.positions
                all_angels[idx] = dis.calc_angles(pos[0,:],pos[1,:],pos[2,:])
                idx+=1
            df[location]=all_angels
        else:
            continue
    return df ,names

def calc_dihedrals(config_file,universe):
    config = open(config_file, "r").read()
    rows = config.split("\n")
    atom_group = []
    group_location = []
    names ={}
    for r in rows:
        if "DihedralGroup" in r:
            location = re.search(r"DihedralGroup: ([A-Z]{4}), Location: (\d.*)", r).group(2)
            group = re.search(r"DihedralGroup: ([A-Z]{4}), Location: (\d.*)", r).group(1)
            atom_group.append(universe.select_atoms("name CA and resid " + str(location)))
            group_location.append(location)
            names[location] = group
        else:
            continue
    for_run=di.Dihedral(atom_group).run().angles
    radians = (((for_run/180)+1)/2)
    radians = np.where(radians<.5,radians+.5,radians-.5)
    radians=radians*2*math.pi
    df = pd.DataFrame()
    df[group_location]=radians
    return df , names


def pdbAnalysis (directory):
    all_pdbs_files = glob.glob(directory+"/*.pdb")
    config = directory+"/config.txt"
    df_dihedral = pd.DataFrame()
    df_bonds = pd.DataFrame()
    dihedral_names = {}
    bonds_names = {}
    first = True
    for i in range(len(all_pdbs_files)):
        print(all_pdbs_files[i])
        universe = MDA.Universe(all_pdbs_files[i])
        if first:
            df_dihedral,dihedral_names = calc_dihedrals(config,universe)
            df_bonds,bonds_names = calc_bonds(config,universe)
            first =False
        else:
            df_dihedral = pd.concat([df_dihedral.reset_index(drop=True),calc_dihedrals(config,universe)[0].reset_index(drop=True)],axis=0)
            df_bonds = pd.concat([df_bonds.reset_index(drop=True),calc_bonds(config,universe)[0].reset_index(drop=True)],axis=0)
    pd.DataFrame.to_csv(df_dihedral,directory+"/all dihedral angels flatten.csv")
    pd.DataFrame.to_csv(df_bonds,directory+"/all bond angels flatten.csv")
    return df_dihedral,df_bonds,dihedral_names,bonds_names

def dist_repeats (directory,title,xlabel):
    '''
    get a dictionary of dataframes with their names, compute the distribution for each one and return
    :param data:
    :return: return overlayers of histogram plots
    '''
    df_dihedral,df_bonds,dihedral_names,bonds_names = pdbAnalysis(directory)
    plt.clf()
    plt.plot()
    sns.kdeplot(data=df_dihedral)
    plt.title(title)
    plt.ylabel("Distribution")
    plt.xlabel(xlabel)
    plt.legend(labels=dihedral_names.values())
    plt.ylim(0,0.4)
    plt.savefig("/home_e/dor/PycharmProjects/pythonProject/graphs/proline/" + str(title)
                + " " + str(xlabel) + " dihedral distribution.jpg")
    plt.clf()
    plt.plot()
    sns.kdeplot(data=df_bonds)
    plt.title(title)
    plt.ylabel("Distribution")
    plt.xlabel(xlabel)
    plt.legend(labels=bonds_names.values())
    plt.ylim(0,0.5)
    plt.savefig("/home_e/dor/PycharmProjects/pythonProject/graphs/" + str(title)
                + " " + str(xlabel) + " bonds distribution.jpg")


# def config_file(directory):
#     seqs = {"seqA": "GSGGGSGGGSG",
#             "seqB": "GSGGPSGGGSG",
#             "seqC": "GSGPPSGGGSG",
#             "seqD": "GSGGPSPGGSG",
#             "seqE": "GSGGPSGPGSG",
#             "seqF": "GSGGPPPGGSG",
#             "seqG": "GSGPSPGPGSG",
#             "seqH": "GSGGPPSPGSG"}
#     for seq in seqs:
#         config = open(directory + seq + "/output/Traj/config.txt", "w+")
#         matches = re.finditer(r'(?=(.{3}))', seqs[seq])
#         results = [(str(match.group(1)), int(match.start() + 1)) for match in matches]
#         for i in results:
#             config.write("BondGroup: " + i[0] + ", Location: " + str(i[1]) + "-" + str(i[1] + 2) + "\n")
#         if seq == "seqA":
#              text = """DihedralGroup: GSGG, Location: 1-4
# DihedralGroup: SGGG, Location: 2-5
# DihedralGroup: GGGS, Location: 3-6
# DihedralGroup: GGSG, Location: 4-7
# DihedralGroup: GSGG, Location: 5-8
# DihedralGroup: SGGG, Location: 6-9
# DihedralGroup: GGGS, Location: 7-10
# DihedralGroup: GGSG, Location: 8-11
#              """
#              config.write(text)
#         else:
#             matches = re.finditer(r'(?=(.{4}))', seqs[seq])
#             results = [(str(match.group(1)), int(match.start() + 1)) for match in matches]
#             for i in results:
#                 if "P" in i[0]:
#                     config.write("DihedralGroup: " + i[0] + ", Location: " + str(i[1]) + "-" + str(i[1] + 3) + "\n")
#                 else:
#                     continue

def config_file(directory,sequences_dict):
    seqs = sequences_dict
    for seq in seqs:
        config = open(directory +"/config.txt", "w+")
        matches = re.finditer(r'(?=(.{3}))', seqs[seq])
        results = [(str(match.group(1)), int(match.start() + 1)) for match in matches]
        for i in results:
            config.write("BondGroup: " + i[0] + ", Location: " + str(i[1]) + "-" + str(i[1] + 2) + "\n")
            matches = re.finditer(r'(?=(.{4}))', seqs[seq])
            results = [(str(match.group(1)), int(match.start() + 1)) for match in matches]
            for i in results:
                if "P" in i[0]:
                    config.write("DihedralGroup: " + i[0] + ", Location: " + str(i[1]) + "-" + str(i[1] + 3) + "\n")
                else:
                    continue


def raw_iter_dir(directory):
    config_file(directory)
    for subdirs in [x for x in os.walk(directory)][0][1]:
            dist_repeats(directory+subdirs+"/output/Traj","Distribution "+subdirs,"Angle")

segment_seqs = {"segA":"CNTGNKSAFPVRFHPHLQPPHHHQNATPSPAAFINNNTAANGSSAGSAWLFPAPATHNIQDEILGSEKAKSQQQEQQDPLEKQQLSPSPGQEAGILPETEKAKSEENQGDNSSENGNGKEKIRIESPVLTGFDYQEATGLGTSTQPLTSSC",
                "segB":"CSSLTGFSNWSAAIAPSSSTIINEDASFFHQGGVPAASANNGALLFQNFPHHVSPGFGGSFSPQIGPLSQHHPHHPHFQHHHSQHQQQRRSPASPHPPPFTHRNAAFNQLPHLANNLNKPPSPWSSYQSPSPTPSSSWSPGGGGYGGWGGSQGRDHRRC",
                "segC":"CLNGGITPLNSISPLKKNFASNHIQLQKYARPSSAFAPKSWMEDSLNRADNIFPFPDRPRTFDMHSLESSLIDIMRAENDTIKARTYGRRRGQSSLFPMEDGFLDDGRGDQPLHSGLGSPHSFSHQNC"}

config_file("/trajectories/dor/coarsed_grain/nothing_noproline/segA/output/Traj/",segment_seqs)
dist_repeats("/trajectories/dor/coarsed_grain/nothing_noproline/segA/output/Traj/","Distribution","Angle")


# def dist_repeats_prepared (dihedral_csv,bonds_csv,config_file,title+
