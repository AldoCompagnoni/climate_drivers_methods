
data {
  int<lower = 0> N;
  vector[N] x;
}
generated quantities {
  real alpha = normal_rng(0, 1);
  real beta = normal_rng(0, 1);
  int y_sim[N] = poisson_log_rng(alpha + beta * x);
}
