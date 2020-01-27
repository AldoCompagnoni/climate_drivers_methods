
data {
  int<lower=0> n_time; // number of data points, length(y)
  int<lower=0> n_lag; // number of monthly lags
  vector[n_time] y;    // response
  matrix[n_time,n_lag] clim;  // matrix of climate covariates
}

parameters {
  
  real alpha;             //intercept
	real<lower = 0> sigma2; //error variance
	vector[n_lag] beta;     // regression parameters
	//hyperparameters prior
  real<lower = 0> lambda; //penalty parameter
  
}

transformed parameters {
  
  real<lower = 0> tau2; //prior variance
	real<lower = 0> y_sd; //error sd (for function "normal")
	vector[n_time] yhat; //mean normal model
	
	yhat    = alpha + clim * beta;
	tau2    = sigma2/lambda;
	y_sd    = sqrt(sigma2);
  
}

model {
  
  // prior regression coefficients: ridge
	beta   ~ normal(0, sqrt(tau2));
	lambda ~ cauchy(0, 1);
	
  //priors nuisance parameters: uniform on log(sigma^2) & mu
	target += -2 * log(y_sd); 
	
	// priors
  sigma2  ~ inv_gamma(0.01,0.01); // this is a variance
  
  y ~ normal(yhat, y_sd);
  
}

//   // hyper-parameters to weight climate effects
//   z ~ normal(0, 1);
//   sigma_beta ~ normal(0, 3);
//   
//   // parameters of data model
//   alpha ~ normal(0, 5);
//   y_sd ~ gamma(1, 1);
//   
//   y ~ normal(yhat, y_sd);
// }

// generated quantities {
//   vector[n_time] log_lik;
// 
//   for (n in 1:n_time) {
//     log_lik[n] = normal_lpdf(y[n] | yhat[n], y_sd);
//   }
// }
