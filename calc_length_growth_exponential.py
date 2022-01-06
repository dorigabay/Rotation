import matplotlib.pyplot as plt
import glob
import re
import numpy
import pandas as pd
import numpy as np
import os
import seaborn as sns
import MDAnalysis as MDA
import math
import sys

def distance_2dots (point1,point2):
    if point1.ndim == 2:
        dis = math.sqrt((point1[0,0]-point2[0,0])**2+(point1[0,1]-point2[:,1])**2+(point1[0,2]-point2[0,2])**2)
        return dis
    if point1.ndim == 1:
        dis = math.sqrt((point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2)
        return dis


def Rg_per_frame (atom_group):
    center = atom_group.center_of_mass()
    all_atom = atom_group
    radi_sum = 0
    for i in all_atom.positions:
        radi_sum+=distance_2dots(center,i)
    return radi_sum/len(all_atom)


def calc_LSE (directory,universe,peptide_length,start_step,end_step):
    universe_CA = universe.select_atoms("name CA")
    N=peptide_length
    n_start = 3
    Rg_n_mean_STD = np.zeros((N-n_start+1,3))
    for n in range(n_start,N+1):
        Rg_collect_time = np.zeros((N-n+1, end_step-start_step))
        for step, t in zip(universe.trajectory[start_step:end_step], range(start_step, end_step+1)):
            Rg_collect = np.zeros(N-n+1)
            for i,j in zip(range(1,N-n+2),range(n,N+1)):
                Rg_collect[i-1] = Rg_per_frame(universe_CA.select_atoms("resid "+str(i)+"-"+str(j)))
            Rg_collect_time[:,(t-start_step)] =Rg_collect
        Rg_n_mean_STD[n-n_start,0]=n
        Rg_n_mean_STD[n-n_start,1]=np.average(Rg_collect_time.ravel()[np.flatnonzero(Rg_collect_time)])
        Rg_n_mean_STD[n-n_start,2]=np.std(Rg_collect_time.ravel()[np.flatnonzero(Rg_collect_time)])
        print("n: "+str(n)+ "is complete")
    pd.DataFrame.to_csv(pd.DataFrame(Rg_n_mean_STD, columns=["n", "Rg", "STD"]),directory+"_LSE.csv")
    return Rg_n_mean_STD


def iter_dir(directory,segment):
    segments = {"segA":151,"segB":159,"segC":128}
    simu_length = {"start":(1000,5001),"continue":(0,5001)}
    all_continue_xtc_files = glob.glob(directory + "/*continue.xtc")
    all_start_xtc_files = glob.glob(directory + "/*start.xtc")
    all_topo_files = glob.glob(directory + "/*.gro")
    if len(all_topo_files)==0:
        all_topo_files = glob.glob(directory + "/*.pdb")
    for i in range(len(all_start_xtc_files)):
        universe = MDA.Universe(all_topo_files[0], all_start_xtc_files[i])
        print("started:   "+all_start_xtc_files[i])
        calc_LSE(all_start_xtc_files[i],universe,segments[segment],simu_length["start"][0],simu_length["start"][1])
    for i in range(len(all_continue_xtc_files)):
        universe = MDA.Universe(all_topo_files[0], all_continue_xtc_files[i])
        calc_LSE(all_continue_xtc_files[i],universe,segments[segment],simu_length["continue"][0],simu_length["continue"][1])

if len(sys.argv)<3:
    exit("Not enough arguments")

iter_dir(sys.argv[1],sys.argv[2])
