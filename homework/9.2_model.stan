data{
    int N;
    real cb[N];
    real fl[N];
    
}

parameters{
    real Delta_G;
    real f0;
    real fq;

    real<lower=0> sigma;
}

transformed parameters{
    real Kd;
    real F[N];
    
    Kd = exp(Delta_G);
    for (i in 1:N) {
        F[i] = f0 * 50 - (2 * (f0 - fq) * 50 * cb[i]) / (Kd + 50 + cb[i] + sqrt(square(Kd + 50 + cb[i])));
    }
}
model{
    // priors
    Delta_G ~ normal(0, 1);
    f0 ~ normal(10000, 1000);
    fq ~ normal(5000, 500);
    sigma ~ normal(0, 5000);
    
    // likelihood 
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
