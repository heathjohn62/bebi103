data{
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