data{
    // model 3
    int N;
    real uM_12[N];
}

parameters{
    real alpha_;
    real_beta_;
}

model{
    // Priors
    alpha_ ~ normal(3, 0.05);
    sigma ~ normal(10, 3);

    // Likelihood
    uM_12 ~ weibull(alpha_, sigma_);
}
