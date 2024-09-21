data {
  int<lower=1> N;
  int<lower=1> NX;
  vector [N] ISP;
  matrix [N,NX] X;
  }

parameters {
  vector[NX] beta;
  // vector[NX - 1] beta_raw;
  real<lower=0> sigma;
}
// 
// transformed parameters {
//   vector[NX] beta;
//   beta[1] = 0;  // Fixing the first value of beta to 0
//   beta[2:NX] = beta_raw;  // The rest of the beta values
// }

model {
  vector[N] mu;
  mu = X * beta ;
  // Likelihood
  ISP ~ lognormal(mu, sigma);

  // Priors
  beta~ normal(0, 5);
  sigma ~ gamma(2, 1.5);
}
