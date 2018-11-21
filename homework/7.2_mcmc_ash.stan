data{
    int N;
    int ASH;
    int trials;
}

parameters{
    real theta;
}

model{
    // Priors
    theta ~ beta(2.0, 7.0);

    // Likelihood
    ASH ~ binomial(trials, theta);
}