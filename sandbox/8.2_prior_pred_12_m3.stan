data{
    // Model 3
    int N;
}

generated quantities{
    // Parameters
    real uM_12[N];
    real alpha;
    real sigma;

    // Likelihood  
    alpha = normal_rng(1, 0.1);
    sigma = normal_rng(400, 100);
    
    // Data
    for (i in 1:N) {
        uM_12[i] = weibull_rng(alpha, sigma);
    }
}