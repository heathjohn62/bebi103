data
{
    int N; // number of total data points combined
    
    int num_divisions; // number of growth events
    int num_cells; // number of cells measured
    
    int div_num[N]; // Denotes growth period per data point
    int bac_id_per_div_num[num_divisions]; // bac_id per division number
    
    int bac_id[N]; // The bacterium being studied at each data point
    real time[N]; // Time elapsed since last division
}

generated quantities
{
    // PPC output
    vector[N] areas; // cell area over time
    vector[N] model_pred; // expected value for cell area

    // bottom layer parameters
    vector<lower=0>[num_divisions] a0_2; // bottom layer growth coefficient 
    vector[num_divisions] b_2;  // bottom layer b
    vector<lower=0>[num_cells] sigma_2; // variation in area measurement (unique per cell)
    vector<lower=1>[num_cells] nu_2;  // tail-weight in area measurement (unique per cell)

    // mid-layer parameters
    vector[num_cells] a0_1; // mid-layer growth coefficient
    vector<lower=0>[num_cells] sigma_a_1; // deviation in a from layer 1 to layer 2
    vector<lower=0>[num_cells] b_1; // mid-layer b
    vector<lower=0>[num_cells] sigma_b_1; // deviation in k from layer 1 to layer 2

    // top-layer hyperparameters
    real a0;
    real<lower=0> b;
    real<lower=0> sigma_a_0;
    real<lower=0> sigma_b_0;

    // We will generate this model from the top down
    // Level 0 hyperparameters
    a0 = normal_rng(600, 100);
    b = fabs(normal_rng(7, 3));
    sigma_a_0 = fabs(normal_rng(0, 50));
    sigma_b_0 = fabs(normal_rng(0, 0.25));

    // mid and low-level parameters
    for (c in 1:num_cells)
    {
        a0_1[c] = a0 + sigma_a_0 * normal_rng(0, 1);
        b_1[c] = b + sigma_b_0 * normal_rng(0, 1);
        
        sigma_2[c] = fabs(normal_rng(20, 1));
        nu_2[c] = 1 + 10 * lognormal_rng(0.1, 1.0);
        
        sigma_a_1[c] = 0.005 * 100 * fabs(normal_rng(0, 1));
        
        sigma_b_1[c] = 0.005 * fabs(normal_rng(0, .1));
    }

    for (i in 1:num_divisions)
    {
        a0_2[i] = a0_1[bac_id_per_div_num[i]] 
                  + sigma_a_1[bac_id_per_div_num[i]] 
                  * normal_rng(0, 1);
        b_2[i] = b_1[bac_id_per_div_num[i]] 
                 + sigma_b_1[bac_id_per_div_num[i]]
                 * normal_rng(0, 1);
    }

    for (n in 1:N)
    {
        model_pred[n] = a0_2[div_num[n]] + b_2[div_num[n]] * time[n];
        areas[n] = student_t_rng(nu_2[bac_id[n]], 
                                 model_pred[n], 
                                 sigma_2[bac_id[n]]);
    }
}