data{
    // Model 2
    int N;
    real uM_12[N];
}

parameters{
    real alpha_;
    real beta_;
}

model{
    // Priors
    alpha_ ~ normal(10, 3);
    beta_ ~ normal(10, 3);

    // Likelihood
    uM_12 ~ gamma(alpha_, beta_);
}

generated quantities{
    // Parameters
    real log_like[N];
    real uM_12_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
        uM_12_ppc[i] = gamma_rng(alpha_, beta_);
    }
    
    // Pointwise likelihood   
    for (i in 1:N) {
        log_like[i] = gamma_lpdf(uM_12[i] | alpha_, beta_);
    }
}
