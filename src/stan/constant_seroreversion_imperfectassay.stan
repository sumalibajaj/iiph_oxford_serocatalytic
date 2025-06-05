functions {
  vector prob_infected_constant_model(
    array[] int ages,
    int n_ages,
    real foi,
    real mu
  ) {
    vector[n_ages] prob_infected;
    for (i in 1:n_ages) {
        prob_infected[i] = foi / (foi + mu) * (1 - exp(-(foi + mu) * ages[i]));  
    }
    return prob_infected;
  }
}

data {
  int<lower=0> n_observations; // each element is an ind or age group or subgroup
  array[n_observations] int n_seropositive;
  array[n_observations] int n_sample;
  array[n_observations] int age_group; // mid point (say) of age group
  real<lower=0> rate_seroreversion; // measured in units , per year (i.e. typical lifetime of detectability is 1 / this)
  real<lower=0, upper=1> sensitivity; // prob(test + | +), lies between 0 and 1
  real<lower=0, upper=1> specificity; // prob(test - | -), lies between 0 and 1
  
  // for simulating posterior predictive observations
  int<lower=1> age_max;
  array[age_max] int ages; // vector of integer ages from 1 to age_max
  int<lower=1> n_ages_fine;
  array[n_ages_fine] real ages_fine; // vector of fine ages to plot seroprev later
}

parameters {
   real<lower=0> foi;
}

transformed parameters {
  vector[n_observations] prob_infected;
  prob_infected = prob_infected_constant_model(age_group,n_observations,foi,rate_seroreversion);
}

model {
  for(i in 1:n_observations){
    real prob_positive = prob_infected[i] * sensitivity + (1 - prob_infected[i]) * (1-specificity);
    n_seropositive[i] ~ binomial(n_sample[i], prob_positive); // likelihood
  }
  
  foi ~ uniform(0, 10); // force of infection prior

}

generated quantities{
  vector[n_observations] log_likelihood;
  for(i in 1:n_observations){
    real prob_positive = prob_infected[i] * sensitivity + (1 - prob_infected[i]) * (1-specificity);
    log_likelihood[i] = binomial_lpmf(n_seropositive[i] | n_sample[i], prob_positive);
  }

  vector[age_max] prob_infected_expanded;
  vector[age_max] foi_expanded;
  for(i in 1:age_max) {
    foi_expanded[i] = foi;
  }
  // calculate posterior for ages 1 to age_max even if not in data-for plotting
	prob_infected_expanded = prob_infected_constant_model(ages, age_max, foi, rate_seroreversion);
	
	vector[n_ages_fine] seroprev_fine;
  for (i in 1:n_ages_fine) {
    real p_inf = foi / (foi + rate_seroreversion) * (1 - exp(-(foi + rate_seroreversion) * ages_fine[i]));
    seroprev_fine[i] = p_inf * sensitivity + (1 - p_inf) * (1 - specificity); // Adjusted for test performance
  }
}
