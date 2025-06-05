# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Define parameters
max_age <- 80
foi <- 0.05  # constant force of infection
mu <- 1/3   # seroreversion rate (e.g., 1/3 years)
sensitivity <- 0.84
specificity <- 0.90
n_sample <- 100

# Define FOI data frame
foi_constant <- data.frame(
  age = seq(1, max_age, 1),
  foi = rep(foi, max_age)
)

# Create survey features (1-year age bands)
survey_features <- data.frame(
  age_min = seq(1, max_age, 1),
  age_max = seq(1, max_age, 1),
  n_sample = rep(n_sample, max_age)
)

# Function to compute true probability of being seropositive (with seroreversion)
ProbSeropos_serorev <- function(age, foi, mu) {
  (foi / (foi + mu)) * (1 - exp(-(foi + mu) * age))
}

# Simulate data
age <- (survey_features$age_min + survey_features$age_max) / 2
n_age_groups <- nrow(survey_features)
true_seropos <- numeric(n_age_groups) # we will fill this
obs_seropos <- numeric(n_age_groups) # we will fill this

for (i in 1:n_age_groups) {
  # 1. True probability of being seropositive
  prob_true <- ProbSeropos_serorev(age[i], foi, mu)
  
  # 2. Simulate number truly seropositive
  n_true_pos <- rbinom(1, n_sample, prob_true)
  n_true_neg <- n_sample - n_true_pos
  
  # 3. Simulate test outcomes based on test characteristics
  n_test_pos <- rbinom(1, n_true_pos, sensitivity) + 
    rbinom(1, n_true_neg, 1 - specificity)
  
  true_seropos[i] <- n_true_pos
  obs_seropos[i] <- n_test_pos
}

# Assemble the final simulated dataset
serosurvey <- data.frame(
  age_min = survey_features$age_min,
  age_max = survey_features$age_max,
  n_sample = survey_features$n_sample,
  n_seropositive = obs_seropos,
  model_type = "Constant FOI with seroreversion and imperfect assay"
)

# Write to CSV (overwrite previous if needed)
write.csv(serosurvey, "data/raw/serodata_constant_seroreversion_imperfectassay.csv", row.names = FALSE)
