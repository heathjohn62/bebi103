data{
    int N;
    real cb[N];
    real fl[N];
    
}

parameters{
    real Delta_G;
    real f0;
    real fq;
    positive_ordered[2] sigma;
    real<lower=0, upper=1> w[N];
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
    
    for (i in 1:N) {
    target += log_mix(w[i],
                      normal_lpdf(fl[i] | F[i], sigma[1]),
                      normal_lpdf(fl[i] | F[i], sigma[2]));
  }
}

generated quantities{
    // Parameters
    real fl_ppc[N];
    
    // Draw from likelihood for post check
    for (i in 1:N) {
    if (uniform_rng(0.0, 1.0) < w[i]) {
      fl_ppc[i] = normal_rng(F[i], sigma[1]);
    }
    else {
      fl_ppc[i] = normal_rng(F[i], sigma[2]);
    }    
  }
    
}
