data{
    int N;
    real cb[N];
}

generated quantities{
    // Parameters
    real Delta_G;
    real<lower=0> Kd;
    real<lower=0> f0;
    real<lower=0> fq;
    real<lower=0> sigma_g;

    real<lower=0> ca;
    vector[N] F;
    vector[N] fl;
       
    Delta_G = normal_rng(0,1); 
    Kd = exp(Delta_G);
    f0 = normal_rng(10000, 1000);
    fq = normal_rng(5000, 500);
    sigma_g = fabs(normal_rng(0, 5000));

    ca = 50;
    
    // Data
    for (i in 1:N) {
        F[i] = f0 * ca - (2 * (f0 - fq) * ca * cb[i]) / (Kd + ca + cb[i] + sqrt(square(Kd + ca + cb[i]) - 4 * ca * cb[i]));
        fl[i] = normal_rng(F[i], sigma_g);
    }
}