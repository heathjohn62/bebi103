data {
  int N;
}


generated quantities {
  int n[N];

  real theta = lognormal_rng(1.0, 8.0);
  
  for (i in 1:N) {
    n[i] = binomial_rng(126, theta);
  }
}
