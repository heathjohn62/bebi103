data{
    int N;
    int WT;
    int ASH;
    int AVA;
    int trials[N];
}

parameters{
    real theta_1;
    real theta_2;
    real theta_3;
}

model{
    // Priors
    theta_1 ~ beta(1, 8);
    theta_2 ~ beta(2, 7);
    theta_3 ~ beta(5, 6);

    // Likelihood
    WT ~ binomial(trials[1], theta_1);
    ASH ~ binomial(trials[2], theta_2);
    AVA ~ binomial(trials[3], theta_3);
}