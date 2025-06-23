# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Define parameters
current_year <- 2025
max_age <- 80
# mu <- 0  # seroreversion rate
# sensitivity <- 1
# specificity <- 1
mu <- 1/3  # seroreversion rate
sensitivity <- 0.84
specificity <- 0.90
n_sample <- 100

# Define age-varying FOI: piecewise constant  
# Define FOI by age
ages <- 0:(max_age)

foi_values <- numeric(length(ages))
for (i in 1:length(ages)) {
  age <- ages[i]
  if (age <= 4) {
    foi_values[i] <- 0.01
  } else if (age <= 14) {
    foi_values[i] <- 0.06
  } else if (age <= 25) {
    foi_values[i] <- 0.04
  } else if (age <= 40) {
    foi_values[i] <- 0.02
  } else {
    foi_values[i] <- 0.01
  }
}

foi_by_age <- data.frame(age = ages, foi = foi_values)

write.csv(foi_by_age, "data/raw/foi_age_seroreversion_imperfectassay.csv", row.names = FALSE)

# Visualize FOI by age
ggplot(foi_by_age, aes(x = ages, y = foi)) +
  geom_step() +
  geom_point() +
  labs(y = "FOI",
       x = "Age") +
  theme_minimal()

# Create survey features
survey_features <- data.frame(
  age_min = seq(1, max_age),
  age_max = seq(1, max_age),
  n_sample = rep(n_sample, max_age)
)

# Function to compute true probability of being seropositive using FOI by age
ProbSeropos_serorev <- function(age, foi_by_age, mu) {
  prob <- 0
  for (j in 1:age) {
    foi <- foi_by_age$foi[foi_by_age$age == j]
    lambda_over_both <- foi / (foi + mu)
    prob <- lambda_over_both + exp(-(foi + mu)) * (prob - lambda_over_both)
  }
  return(prob)
}

# Simulate data

# Simulate seropositivity data
age <- (survey_features$age_min + survey_features$age_max) / 2
n_age_groups <- nrow(survey_features)
true_seropos <- numeric(n_age_groups) # we will fill this
obs_seropos <- numeric(n_age_groups) # we will fill this

for (i in 1:n_age_groups) {
  prob_true <- ProbSeropos_serorev(age[i], foi_by_age, mu)
  
  n_true_pos <- rbinom(1, n_sample, prob_true)
  n_true_neg <- n_sample - n_true_pos
  
  n_test_pos <- rbinom(1, n_true_pos, sensitivity) +
    rbinom(1, n_true_neg, 1 - specificity)
  
  true_seropos[i] <- n_true_pos
  obs_seropos[i] <- n_test_pos
}

# Assemble the final simulated dataset
serosurvey <- data.frame(
  survey_year = current_year,
  age_min = survey_features$age_min,
  age_max = survey_features$age_max,
  n_sample = survey_features$n_sample,
  n_seropositive = obs_seropos,
  model_type = "Age-varying FOI with seroreversion and imperfect assay"
)

# Visualise
ggplot(serosurvey, aes(x = age_min, y = n_seropositive/n_sample)) +
  geom_line() +
  geom_point() +
  labs(title = "Observed Seroprevalence by Age",
       y = "Observed Seroprevalence",
       x = "Age") +
  theme_minimal()

# Write to CSV
write.csv(serosurvey, "data/raw/serodata_age_seroreversion_imperfectassay.csv", row.names = FALSE)
