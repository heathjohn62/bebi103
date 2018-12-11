import numpy as np
import pandas as pd
import pystan
import sys
import bebi103
import bokeh

data = pd.read_csv("./9.1_data_wrangled.csv")
data = data.drop(columns = ["Unnamed: 0"])

mcmc_exp_2 = bebi103.stan.StanModel(file='./9.1_mcmc_exp_2.stan')

# I want to take random divsions. I chose 2 and 41, one from each cell
toy_data = pd.concat([data[data["div_num"]==2], data[data["div_num"]==41]])
# I need the division number to be [1, 2], and this can be changed with quick arithmetic. 
toy_data.loc[:,"div_num"] = (toy_data.copy()["div_num"].values % 2) + 1

df_toy_exp_2 = None
toy_samples_exp_2 = None
try:
    df_toy_exp_2 = pd.read_csv("./9.1_toy_mcmc_exp_2.csv")
    df_toy_exp_2 = df_toy_exp_2.drop(columns = ['Unnamed: 0'])
except FileNotFoundError:
    toy_exp_dict = {"N":len(toy_data.index),
                    "num_divisions": 2,
                    "num_cells": 2,
                    "div_num": toy_data["div_num"].values.astype(int),
                    "bac_id_per_div_num": [1,2],
                    "bac_id": toy_data["bac_id"].values.astype(int),
                    "time": toy_data["time (min)"].values,
                    "areas": toy_data['area (sq µm)'].values}
    toy_samples_exp_2 = mcmc_exp_2.sampling(data=toy_exp_dict, 
                                            control=dict(adapt_delta = .99,
                                                       max_treedepth = 13),
                                            warmup=2000, 
                                            iter=6000, 
                                            thin=2)
   
    df_toy_exp_2 = bebi103.stan.to_dataframe(toy_samples_exp_2, 
                                             inc_warmup = False, 
                                             diagnostics=False)
    df_toy_exp_2.to_csv("./9.1_toy_mcmc_exp_2.csv")
    
    # Check mcmc results
    orig_stdout = sys.stdout
    f = open('./9.1_mcmc_toy_exp_2_log.txt', 'a')
    sys.stdout = f
    bebi103.stan.check_all_diagnostics(toy_samples_exp_2)
    sys.stdout = orig_stdout
    f.close()