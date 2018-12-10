data{
    int N;
    real cb[N];
}

generated quantities{
    // Parameters
    real Kd;
    real f0;
    real fq;
    real sigma;
    real ca;
    real F[N];
    real likelihood[N];
       
    Kd = lognormal_rng(1, 3);
    f0 = normal_rng(10000, 1000);
    fq = normal_rng(1000, 100);
    sigma = fabs(normal_rng(0, 5000));

    ca = 50;
    
    // Data
    for (i in 1:N) {
        F[i] = f0 * ca - (2 * (f0 - fq) * ca * cb[i]) / (Kd + ca + cb[i] + sqrt(square(Kd + ca + cb[i]) - 4 * ca * cb[i]));
        likelihood[i] = normal_rng(F[i], sigma);
    }
}