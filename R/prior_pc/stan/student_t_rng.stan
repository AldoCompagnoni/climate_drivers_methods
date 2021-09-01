
data {
  int N;
  real gamma_shape;
}

generated quantities {
  
  real y_sim[N] = student_t_rng(0,1,);
  
}
