data{
    int N; // number of total data points combined
    
    int num_divisions; // number of growth events
    int num_cells; // number of cells measured
    
    int div_num[N]; // Denotes growth period per data point
    int bac_id_per_div_num[num_divisions]; // bac_id per division number
    
    int bac_id[N]; // The bacterium being studied at each data point
    
    // actual measurments 
    real areas[N];
    real time[N]; // Time elapsed since last division

}

parameters{

    // top-layer hyperparameters
    real a0;
    real<lower=0> b;
    real<lower=0> sigma_a_0;
    real<lower=0> sigma_b_0;
    
    // bottom layer parameters
    //vector<lower=0>[num_divisions] a0_2; // bottom layer growth coefficient 
    //vector[num_divisions] b_2;  // bottom layer b
    // vector<lower=0>[num_cells] sigma_2; // variation in area measurement (unique per cell)
    //vector<lower=1>[num_cells] nu_2;  // tail-weight in area measurement (unique per cell)
    
    vector[num_cells] a0_1_til;
    vector[num_cells] sigma_2_til;
    vector[num_cells] nu_2_til;
    vector[num_cells] b_1_til;
}

transformed parameters{

    // mid-layer parameters
    // vector[num_cells] a0_1; // mid-layer growth coefficient
    // vector<lower=0>[num_cells] sigma_a_1; // deviation in a from layer 1 to layer 2
    // vector<lower=0>[num_cells] b_1; // mid-layer b
    // vector<lower=0>[num_cells] sigma_b_1; // deviation in k from layer 1 to layer 2
    
    vector[num_cells] a0_1 = a0 + sigma_a_0 * a0_1_til;
    vector[num_cells] b_1 = b + sigma_b_0 * a0_1_til;
        
    vector[num_cells] sigma_2 = fabs(sigma_2_til);
    vector[num_cells] nu_2 = 1 + 10 * nu_2_til;
        
    vector[num_cells] sigma_a_1 = 0.005 * 100 * fabs(a0_1_til);
    vector[num_cells] sigma_b_1 = 0.005 * fabs(b_1_til);
    
    vector[N] model_pred; // expected value for cell area
    
    // ----
    
     vector[num_divisions] a0_2 = a0_1[bac_id_per_div_num] 
                  + sigma_a_1[bac_id_per_div_num] 
                  * a0_1_til;
     vector[num_divisions] b_2 = b_1[bac_id_per_div_num] 
                 + sigma_b_1[bac_id_per_div_num] * a0_1_til;
                 
    // area before we add the student t error to it
    for (n in 1:N)
    {
        model_pred[n] = a0_2[div_num[n]] + b_2[div_num[n]] * time[n];
    }
}

model{

    a0_1_til ~ normal(0, 1);
    sigma_2_til ~ normal(20, 1);
    nu_2_til ~ lognormal(0.1, 1.0);
    b_1_til ~ normal(0, .1);

    // Level 0 hyperparameters
    a0 ~ normal(600, 100);
    b ~ fabsnormal(7, 3); // should be fabs
    sigma_a_0 ~ normal(0, 50); // should be fabs
    sigma_b_0 ~ normal(0, 0.25); // should be fabs
    
    for (n in 1:N){
        areas[n] ~ student_t(nu_2[bac_id[n]], 
                                 model_pred[n], 
                                 sigma_2[bac_id[n]])
    }
}