
data {
  int n_time;
  vector[n_time] clim_means;
}

generated quantities {
  real y_sim[n_time]; 
  real alpha         = normal_rng(0, 1);
  real beta          = normal_rng(0, 1);
  real<lower=0> y_sd = student_t_rng(3, 0, 1);
  
  y_sim = normal_rng(alpha + beta * clim_means, y_sd);
  
}
