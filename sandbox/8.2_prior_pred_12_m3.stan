data{
    // Model 3
    int N;
}

generated quantities{
    // Parameters
    real time[N];
    real alpha;
    real sigma;

    // Likelihood  
    alpha = normal_rng(1, 0.1);
    sigma = normal_rng(400, 100);
    
    // Data
    for (i in 1:N) {
        time[i] = weibull_rng(alpha, sigma);
    }
}