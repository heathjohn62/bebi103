import numpy as np
import pandas as pd
import pystan
import sys
import bebi103
import bokeh

samples_exp, model = bebi103.stan.pickle_load_samples("./9.1_mcmc_exp_2.pkl")

print("Extracting to Dataframe!")
df_exp = bebi103.stan.to_dataframe(samples_exp, 
                                    inc_warmup = False, 
                                    diagnostics=False)

print("Saving dataframe!")
df_exp.to_csv("./9.1_mcmc_exp_2.csv")

print("Checking Diagnostics!")
#  Check mcmc results
orig_stdout = sys.stdout
f = open('./9.1_mcmc_exp_2_log.txt', 'a')
sys.stdout = f
print("Log file for Exponential MCMC")
bebi103.stan.check_all_diagnostics(samples_exp)

# Computing statistics for model comparison.
waic_results = bebi103.stan.waic(samples_exp, log_likelihood='log_lik')
loo_results = bebi103.stan.loo(samples_exp, log_likelihood='log_lik')

print('Exponential Model Effectiveness Statistics:')
print(waic_results, end='\n\n')
print(loo_results, end='\n\n\n\n')

sys.stdout = orig_stdout
f.close()
print("Done!")
