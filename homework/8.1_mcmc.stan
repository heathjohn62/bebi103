data 
{
    int N;
    int num_exps;
    int experiment[N];
    real fluorescence[N];
}

parameters 
{
    real<lower=0> theta_H;
    real<lower=0> sigma_H;
    real<lower=0> sigma_P;
    vector[num_exps] theta_tilde;
}

transformed parameters
{
    vector[num_exps] theta;
    
    theta = theta_H + theta_tilde * sigma_P;
}

model
{
    // hyperparameters
    theta_H ~ lognormal(0.25, 0.5); 
    sigma_H ~ normal(0, 0.5);
    sigma_P ~ normal(0, 0.5);
 

    // direct parameters
    theta_tilde ~ normal(0,1);
    
    // data!
    fluorescence ~ normal(theta[experiment], sigma_H);
}

// posterior predictive check
generated quantities 
{
    real fluorescence_ppc[N];
    for (i in 1:N) 
    {
        fluorescence_ppc[i] = normal_rng(theta[experiment[i]], sigma_H);
    }
}