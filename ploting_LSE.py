import glob
import re
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from lmfit import Model


def func(x,b):#,a):
    # return a*(x**b)
    return 2*(x**b)


def calc_Rsquared(x,y,parameters):
    y_mean = np.mean(y)
    total_sum_of_squares = np.sum((y-y_mean)**2)
    residuals_sum_of_squares = np.sum((y-func(x,*parameters))**2)
    Rsquared = 1-(residuals_sum_of_squares/total_sum_of_squares)
    return Rsquared


def fit_a_curve(train_data):
    polymer_length = train_data["n"]
    Rg = train_data["Rg"]
    STD = train_data["STD"]
    parameters, pcov = curve_fit(func, polymer_length, Rg, sigma=STD)
    perr = (np.sqrt(np.diag(pcov)))
    Rsquared_train = calc_Rsquared(polymer_length,Rg,parameters)
    return parameters,perr,Rsquared_train


def split_data(directory):
    all_LSE_files = glob.glob(directory + "*_LSE.csv")
    train_files = all_LSE_files[:-1]
    test_data = pd.read_csv(all_LSE_files[-1])
    train_data = pd.DataFrame()
    name = re.search(".*/(seg.).*", train_files[0]).group(1)
    for file in train_files:
        LSE = pd.read_csv(file)
        train_data = pd.concat([train_data, LSE], axis=0, ignore_index=True)
    return train_data,test_data,name

def plot_LSE(directory):
    train_data, test_data, name = split_data(directory)
    parameters, perr, Rsquared_train = fit_a_curve(train_data)
    polymer_length = test_data["n"]
    Rg = test_data["Rg"]
    STD = test_data["STD"]
    fig, ax1 = plt.subplots(1, 1)
    ax1.plot(polymer_length, Rg, color="blue", label="Rg for polymer length")
    ax1.fill_between(polymer_length, Rg - STD, Rg + STD, color="blue", alpha=0.2)
    y_expected = func(polymer_length, *parameters)
    Rsquared_test = calc_Rsquared(polymer_length, Rg, parameters)
    R_train ="R^2(train): "+str(round(Rsquared_train,2))
    R_test ="R^2(test): "+str(round(Rsquared_test,2))
    a="a: 2"
    # a = "a: "+str(round(parameters[1],2))+"  STD: "+str(round(perr[1],3))
    b = "b: "+str(round(parameters[0],2))+"  STD: "+str(round(perr[0],3))
    ax1.plot(polymer_length,y_expected,color="red",label= "fitted function\n"+a+"\n"+b+"\n"+R_train+"\n"+R_test)
    ax1.set_title("Rg as a function of polymer length, with the best fiiting curve")
    ax1.legend()
    # ax2.fill_between(polymer_length, Rg - STD, Rg + STD, color="green", alpha=0.2)
    # ax2.loglog(polymer_length,Rg,color = "green",label="log log")
    # ax2.plot(polymer_length, y_expected, color="red", label="fitted function")
    # plt.style.use("seaborn")
    ax1.set_ylabel("Rg ($\AA$")
    ax1.set_xlabel("polymer length")
    plt.legend()
    plt.savefig("/home_e/dor/PycharmProjects/pythonProject/graphs/atomistic/LSE_charmm2/LSE_a2_"+name+".png")


dir = {"segment A":"/trajectories/dor/atomistic/charm_ff/his_nonprotonated/salt125mM/segA/files_for_analysis/",
       "segment B":"/trajectories/dor/atomistic/charm_ff/his_nonprotonated/salt125mM/segB/files_for_analysis/LSE/",
       "segment C":"/trajectories/dor/atomistic/charm_ff/his_nonprotonated/salt125mM/segC/files_for_analysis/LSE/"}


for i in dir:
    plot_LSE(dir[i])