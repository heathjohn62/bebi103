data{
    int N;
    int WT;
    int trials;
}

parameters{
    real theta;
}

model{
    // Priors
    theta ~ beta(1.0, 8.0);

    // Likelihood
    WT ~ binomial(trials, theta);
}