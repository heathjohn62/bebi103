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
    real<lower=0> sigma_P_tilde;
    vector[num_exps] theta;
}

transformed parameters
{
    real<lower=0> sigma_P = sigma_P_tilde * .5;
}

model
{
    // hyperparameters
    theta_H ~ lognormal(0.25, 0.5); 
    sigma_H ~ normal(0, 0.5);
    sigma_P_tilde ~ normal(0, 1);

    // direct parameters
    for (i in 1:num_exps)
    { 
        theta[i] ~ normal(theta_H, sigma_P);
    }

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