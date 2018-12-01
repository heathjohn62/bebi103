data{
    // Model 2
    int N;
}

generated quantities{
    // Parameters
    real time[N];
    real beta_;
    real alpha;
       
    alpha = normal_rng(10, 3);
    beta_ = 1 / normal_rng(150, 50);
    
    // Data
    for (i in 1:N) {
        time[i] = gamma_rng(alpha, beta_);
    }
}