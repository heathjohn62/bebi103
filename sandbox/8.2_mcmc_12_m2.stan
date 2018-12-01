data{
    // Model 2
    int N;
    real uM_12[N];
}

parameters{
    real alpha_;
    real tao;
}

transformed parameters{
    real beta_;
    beta_ = 1.0 / tao;
}

model{
    // Priors
    alpha_ ~ normal(10, 3);
    tao ~ normal(150, 50);

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
