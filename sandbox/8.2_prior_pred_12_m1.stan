data{
    int N;
}

generated quantities{
    // Parameters
    real uM_12[N];
    real beta_;
       
    beta_ = 1.0/(20.0 * lognormal_rng(0.1, 0.7));
    
    // Data
    for (i in 1:N) {
        uM_12[i] = exponential_rng(beta_);
    }
}
