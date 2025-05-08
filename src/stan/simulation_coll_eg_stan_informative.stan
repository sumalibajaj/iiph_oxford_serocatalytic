data{
  int<lower=0> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] y;
}
parameters{
  real alpha;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}
model{
  y ~ normal(alpha + beta1 * x1 + beta2 * x2, sigma);
  
  // add prior to make informative posteriors- change/remove & see what happens
  beta1 ~ normal(3, 0.01);
  // beta1 ~ normal(3, 100);
}

generated quantities{
  vector[N] y_rep;  // Posterior predictive samples
  for (n in 1:N){
    y_rep[n] = normal_rng(alpha + beta1 * x1[n] + beta2 * x2[n], sigma);
  }
}
