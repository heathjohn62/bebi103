data{
    // Model 1
    int N;
}

generated quantities{
    // Parameters
    real uM_12[N];
    real beta_;
       
    beta_ = 1.0 / normal_rng(400, 100);
    
    // Data
    for (i in 1:N) {
        uM_12[i] = exponential_rng(beta_);
    }
}
