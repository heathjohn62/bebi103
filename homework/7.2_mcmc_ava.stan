data{
    int N;
    int AVA;
    int trials;
}

parameters{
    real theta;
}

model{
    // Priors
    theta ~ beta(5.0, 6.0);

    // Likelihood
    AVA ~ binomial(trials, theta);
}