import numpy as np
import pandas as pd
import pystan
import sys
import bebi103
import bokeh

data = pd.read_csv("./9.1_data_wrangled.csv")
data = data.drop(columns = ["Unnamed: 0"])

mcmc_exp = bebi103.stan.StanModel(file='./9.1_mcmc_exp.stan')

# let's just see how many divisions we have in total. 
# We will need to feed this parameter into stan
div_num = data["div_num"].values
N = len(div_num)
num_divisions = len(np.unique(div_num))

# I need to make an array that specifies the bacterium for each division
unique_div_nums = np.unique(div_num)
bac_id_per_div_num = np.zeros(len(unique_div_nums))
for i in range(0, len(unique_div_nums)):
    div = unique_div_nums[i]
    sub_df = data[data["div_num"]==div]
    bac_id_per_div_num[i] = sub_df.iloc[0]["bac_id"]

df_exp = None
samples_exp = None
try:
    df_exp= pd.read_csv("./9.1_mcmc_exp.csv")
    df_exp = df_exp.drop(columns = ['Unnamed: 0'])
except FileNotFoundError:
    exp_dict = {"N":N,
                "num_divisions": num_divisions,
                "num_cells": 2,
                "div_num": data["div_num"].values.astype(int),
                "bac_id_per_div_num": bac_id_per_div_num.astype(int),
                "bac_id": data["bac_id"].values.astype(int),
                "time": data["time (min)"].values,
                "areas": data['area (sq µm)'].values}
    samples_exp = mcmc_exp.sampling(data=exp_dict, 
                                    control=dict(adapt_delta = .99,
                                                 max_treedepth = 13),
                                    warmup=2000, 
                                    iter=6000, 
                                    thin=1)
    
    df_exp = bebi103.stan.to_dataframe(samples_exp, 
                                       inc_warmup = False, 
                                       diagnostics=False)
    df_exp.to_csv("./9.1_mcmc_exp.csv")

    #  Check mcmc results
    orig_stdout = sys.stdout
    f = open('./9.1_mcmc_exp_log.txt', 'a')
    sys.stdout = f
    print("Log file for Exponential MCMC")
    print("John Heath")
    bebi103.stan.check_all_diagnostics(samples_exp)
    sys.stdout = orig_stdout
    f.close()
