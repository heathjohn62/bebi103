data{
    // Model 3
    int N;
    real time[N];
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
    time ~ weibull(alpha, sigma);
}

generated quantities{
    // Parameters
    real log_like[N];
    real time_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
        time_ppc[i] = weibull_rng(alpha, sigma);
    }
    
    // Pointwise likelihood   
    for (i in 1:N) {
        log_like[i] = weibull_lpdf(time[i] | alpha, sigma);
    }
}
