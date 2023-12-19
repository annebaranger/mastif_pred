data {
  int<lower=1> N;
  int<lower=1> NX;
  int<lower=1> S;
  int<lower=1,upper=S> species[N];
  vector [N] ISP_meas;
  matrix [N,NX] X;
  vector [N] s_ISP;
  real <lower=0> prior_mean;
  real <lower=0> prior_sd;
  }

parameters {
  vector <lower=0> [N] ISP_real;
  vector[NX - 1] beta_raw;
  vector[S] alpha;
  real<lower=0> sigma;
  real<lower=0> sigma_s;
}

transformed parameters {
  vector[NX] beta;
  beta[1] = 0;  // Fixing the first value of beta to 0
  beta[2:NX] = beta_raw;  // The rest of the beta values
}

model {
  vector[N] mu;
  mu = X * beta + alpha[species];
  // vector [N] ISP_real;
  // Likelihood
  ISP_meas ~ normal(ISP_real, s_ISP);
  (log(ISP_real)-prior_mean)/prior_sd ~ normal(mu, sigma);
  // ISP_real ~ normal(mu, sigma);

  // Priors
  beta_raw ~ normal(0, 10);
  alpha ~ normal(0, sigma_s);
  sigma ~ gamma(2, 1);
  sigma_s ~ gamma(2, 1);
}
