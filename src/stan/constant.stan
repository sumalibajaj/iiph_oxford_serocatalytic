functions {
  vector prob_infected_constant_model(
    array[] int ages,
    int n_ages,
    real foi
  ) {
    vector[n_ages] prob_infected;
    for (i in 1:n_ages) {
        prob_infected[i] = 1 - exp(-foi * ages[i]); 
    }
    return prob_infected;
  }
}

data {
  int<lower=0> n_observations; // each element is an ind or age group or subgroup
  int<lower=1> age_max;
  array[age_max] int ages; // vector of integer ages from 1 to age_max
  int<lower=1> n_ages_fine;
  array[n_ages_fine] real ages_fine; // vector of fine ages to plot seroprev later
  array[n_observations] int n_seropositive;
  array[n_observations] int n_sample;
  array[n_observations] int age_group; // mid point (say) of age group
}

parameters {
   real<lower=0> foi;
}

transformed parameters {
  vector[n_observations] prob_infected;
  prob_infected = prob_infected_constant_model(age_group,n_observations,foi);
}

model {
  for(i in 1:n_observations){
     n_seropositive[i] ~ binomial(n_sample[i], prob_infected[i]); // likelihood
  }
  
    foi ~ uniform(0, 10); // force of infection prior

}

generated quantities{
  vector[n_observations] log_likelihood;
  for(i in 1:n_observations){
    log_likelihood[i] = binomial_lpmf(n_seropositive[i] | n_sample[i], prob_infected[i]);
  }

  vector[age_max] prob_infected_expanded;
  vector[age_max] foi_expanded;
  for(i in 1:age_max) {
    foi_expanded[i] = foi;
  }
  // calculate posterior for ages 1 to age_max even if not in data-for plotting
	prob_infected_expanded = prob_infected_constant_model(ages,age_max,foi);
	
	vector[n_ages_fine] seroprev_fine;
  for (i in 1:n_ages_fine) {
    seroprev_fine[i] = 1 - exp(-foi * ages_fine[i]);
  }
}
