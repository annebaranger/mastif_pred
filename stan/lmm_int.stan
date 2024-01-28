data {
  int<lower=1> N;
  int<lower=1> NX;
  int<lower=1> S;
  int<lower=1,upper=S> species[N];
  vector [N] ISP;
  matrix [N,NX] X;
  // vector [N] sISP;
  }

// parameters {
//   vector[NX] beta;
//   vector[S] alpha;
//   real<lower=0> sigma;
//   real<lower=0> sigma_s;
//   }
parameters {
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

// model {
//   vector[N] mu;
//   mu = X*beta+alpha[species];
//   //Likelihood
//   ISP~normal(mu,sigma);
//   
//   //Priors
//   beta ~ normal(0,1000);
//   alpha~normal(0,sigma_s);
//   sigma~cauchy(0,5);
//   sigma_s~cauchy(0,5);
//   }
model {
  vector[N] mu;
  mu = X * beta + alpha[species];
  // Likelihood
  ISP ~ normal(mu, sigma);
  
  // Priors
  beta_raw ~ normal(0, 1000);
  alpha ~ normal(0, sigma_s);
  sigma ~ cauchy(0, 5);
  sigma_s ~ cauchy(0, 5);
}
