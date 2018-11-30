data{
    // Model 1
    int N;
    real uM_12[N];
}

parameters{
    real tao;
}

transformed parameters{
    real beta_;
    beta_ = 1.0 / tao;
}

model{
    // Priors
    tao ~ normal(700, 100);

    // Likelihood
    uM_12 ~ exponential(beta_);
}

generated quantities{
    // Parameters
    real log_like[N];
    real uM_12_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
        uM_12_ppc[i] = exponential_rng(beta_);
    }
    
    // Pointwise likelihood   
    for (i in 1:N) {
        log_like[i] = exponential_lpdf(uM_12[i] | beta_);
    }
}