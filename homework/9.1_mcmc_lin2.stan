data
{
    int N; // number of total data points combined
    int num_divisions; // number of growth events
    int num_cells; // number of cells measured
    int div_num[N]; // Denotes growth period per data point
    int bac_id_per_div_num[num_divisions]; // bac_id per division number
    int bac_id[N]; // The bacterium being studied at each data point
    real time[N]; // Time elapsed since last division
    real areas[N]; // Cell area over time
}

parameters
{
    // top-layer hyperparameters
    real<lower=0> a0;
    real<lower=0> b;
    real<lower=0> sigma_a_0;
    real<lower=0> sigma_b_0_tilde;

    // bottom layer parameters
    vector<lower=0>[num_divisions] a0_2_tilde; // bottom layer growth coefficient 
    vector<lower=0>[num_divisions] b_2_tilde;  // bottom layer growth rate
    vector<lower=0>[num_cells] sigma_2; // variation in area measurement (unique per cell)
    vector<lower=0>[num_cells] nu_2_tilde;  // tail-weight in area measurement (unique per cell)

    // mid-layer parameters
    vector<lower=0>[num_cells] a0_1_tilde; // mid-layer growth coefficient
    vector<lower=0>[num_cells] sigma_a_1_tilde; 
    vector<lower=0>[num_cells] b_1_tilde; // mid-layer growth rate
    vector<lower=0>[num_cells] sigma_b_1_tilde;
}

transformed parameters
{
    real<lower=0> sigma_b_0;
    vector[N] model_pred; // expected value for cell area

    // mid-layer
    vector<lower=0>[num_cells] a0_1; // mid-layer growth coefficient
    vector<lower=0>[num_cells] sigma_a_1; 
    vector<lower=0>[num_cells] b_1; // mid-layer growth rate
    vector<lower=0>[num_cells] sigma_b_1;

    // bottom layer
    vector<lower=0>[num_divisions] a0_2; // bottom layer growth coefficient 
    vector<lower=0>[num_divisions] b_2;  // bottom layer growth rate
    vector<lower=1>[num_cells] nu_2;

    for (c in 1:num_cells)
    {
        nu_2[c] = 1 + 10 * nu_2_tilde[c];
        a0_1[c] = a0 + sigma_a_0 * a0_1_tilde[c];
        sigma_a_1[c] = 0.005 * sigma_a_1_tilde[c];
        b_1[c] = b + sigma_b_0 * b_1_tilde[c];
        sigma_b_1[c] = 0.005 * sigma_b_1_tilde[c];
    }

    for (i in 1:num_divisions)
    {
        a0_2[i] = a0_1[bac_id_per_div_num[i]] 
                + sigma_a_1[bac_id_per_div_num[i]] 
                * a0_2_tilde[i];
        b_2[i] = b_1[bac_id_per_div_num[i]] 
                    + sigma_b_1[bac_id_per_div_num[i]]
                    * b_2_tilde[i];
    }
    
    // Evaluate transformed params
    for (n in 1:N)
    {
        model_pred[n] = a0_2[div_num[n]] + b_2[div_num[n]] * time[n];
    }

}

model
{
    // priors on hyperparameters
    a0 ~ normal(600, 75);
    b ~ normal(7, 3); 
    
    sigma_a_0 ~ normal(0, 50);
    sigma_b_0_tilde ~ normal(0, 1.0);
    
    sigma_2 ~ normal(0, 20);
    nu_2_tilde ~ lognormal(0.1, 1.0);
    
    a0_1_tilde ~ normal(0, 1);
    sigma_a_1_tilde ~ normal(0, 1);
    b_1_tilde ~ normal(0, 1);
    sigma_b_1_tilde ~ normal(0, 1.0);
    a0_2_tilde ~ normal(0, 1);
    b_2_tilde ~ normal(0, 1);
    // likelihood
    for (n in 1:N)
    {
        areas[n] ~ student_t(nu_2[bac_id[n]], 
                             model_pred[n], 
                             sigma_2[bac_id[n]]);
    }
}


generated quantities
{
    real areas_ppc[N];
    // perform posterior predictive check
    for (n in 1:N)
    {
        areas_ppc[n] = student_t_rng(nu_2[bac_id[n]], 
                                  model_pred[n], 
                                  sigma_2[bac_id[n]]);
    }
}