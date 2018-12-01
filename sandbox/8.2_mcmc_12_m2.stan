data{
    // Model 2
    int N;
    real time[N];
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
    time ~ gamma(alpha_, beta_);
}

generated quantities{
    // Parameters
    real log_like[N];
    real time_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
        time_ppc[i] = gamma_rng(alpha_, beta_);
    }
    
    // Pointwise likelihood   
    for (i in 1:N) {
        log_like[i] = gamma_lpdf(time[i] | alpha_, beta_);
    }
}
