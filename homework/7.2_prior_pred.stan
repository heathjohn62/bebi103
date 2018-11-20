data{
    int N;
    real alpha[N];
    real beta_[N];
    int trials[N];
}

generated quantities{
    // Parameters
    real theta[N];
    real rev_prob[N];
    for (i in 1:N) {
        theta[i] = beta_rng(alpha[i], beta_[i]);
    }

    // Data
    for (i in 1:N) {
        rev_prob[i] = binomial_rng(trials[i], theta[i]) / trials[i];
    }
}