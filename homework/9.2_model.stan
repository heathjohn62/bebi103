data{
    int N;
    real cb[N];
}

parameters{
    real Kd;
    real f0;
    real fq;
    real<lower=0> sigma;
    real fl[N];
}

transformed parameters{
    real F[N];
    real ca;
    ca = 50;
    for (i in 1:N) {
        F[i] = f0 * ca - (2 * (f0 - fq) * ca * cb[i]) / (Kd + ca + cb[i] + sqrt(square(Kd + ca + cb[i])));
    }
}
model{
    // priors
    Kd ~ lognormal(1, 3);
    f0 ~ normal(500000, 50000);
    fq ~ normal(50000, 5000);
    sigma ~ normal(0, 5000);
    
    // likelihood - oh the problem is that this should be F[i]!
    fl ~ normal(F, sigma);
}

generated quantities{
    // Parameters
    real log_like[N];
    real fl_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
        fl_ppc[i] = normal_rng(F[i], sigma);
    }
    
    // Pointwise likelihood   
    for (i in 1:N) {
        log_like[i] = normal_lpdf(fl[i] | F, sigma);
    }
}
