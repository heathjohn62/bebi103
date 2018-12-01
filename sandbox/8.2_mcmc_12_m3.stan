data{
    // Model 3
    int N;
    real uM_12[N];
}

parameters{
    real alpha;
    real sigma;
}

model{
    // Priors
    alpha ~ normal(1, 0.1);
    sigma ~ normal(400, 100);

    // Likelihood
    uM_12 ~ weibull(alpha, sigma);
}

generated quantities{
    // Parameters
    real log_like[N];
    real uM_12_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
        uM_12_ppc[i] = weibull_rng(alpha, sigma);
    }
    
    // Pointwise likelihood   
    for (i in 1:N) {
        log_like[i] = weibull_lpdf(uM_12[i] | alpha, sigma);
    }
}
