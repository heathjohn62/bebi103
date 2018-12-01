data{
    // Model 2
    int N;
}

generated quantities{
    // Parameters
    real uM_12[N];
    real beta_;
    real alpha;
       
    alpha = normal_rng(10, 3);
    beta_ = normal_rng(10, 3);
    
    // Data
    for (i in 1:N) {
        uM_12[i] = gamma_rng(alpha, beta_);
    }
}