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
    alpha = normal_rng(3, 0.05);
    sigma = normal_rng(10, 3);
    
    // Data
    for (i in 1:N) {
        uM_12[i] = weibull_rng(alpha, sigma);
    }
}