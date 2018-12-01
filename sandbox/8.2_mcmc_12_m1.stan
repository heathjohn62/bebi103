data{
    // Model 1
    int N;
    real time[N];
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
    tao ~ normal(400, 100);

    // Likelihood
    time ~ exponential(beta_);
}

generated quantities{
    // Parameters
    real log_like[N];
    real time_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
        time_ppc[i] = exponential_rng(beta_);
    }
    
    // Pointwise likelihood   
    for (i in 1:N) {
        log_like[i] = exponential_lpdf(time[i] | beta_);
    }
}