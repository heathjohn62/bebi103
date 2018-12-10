data {
    int N;
    int num_exps;
    int experiment[N];
}

generated quantities {
    real<lower=0> theta_H;
    real<lower=0> sigma_H;
    real<lower=0> sigma_P;
    vector[num_exps] theta;
    vector[N] fluorescence;
    
    // Draw parameters
    theta_H = lognormal_rng(0.25, 0.5); 
    // Here, fabs is used to draw from the half-normal distribution.
    sigma_H = fabs(normal_rng(0, .5));
    sigma_P = fabs(normal_rng(0, .5));
    
    for (i in 1:num_exps)
    { 
        theta[i] = normal_rng(theta_H, sigma_P);
    }
    for (i in 1:N) 
    {
        fluorescence[i] = normal_rng(theta[experiment[i]], sigma_H);
    }
}