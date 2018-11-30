data{
    int N;
    real uM_12[N];
}

parameters{
    real alpha_;
    real beta_;
}

transformed parameters{
    real alpha_beta;
    alpha_beta = normal(700, 100)
    beta = alpha_ / beta_;
}

model{
    // Priors
    tao ~ normal(700, 100);
    
    alpha ~ geometric(0.1);

    // Likelihood
    uM_12 ~ gamma(beta_);
}