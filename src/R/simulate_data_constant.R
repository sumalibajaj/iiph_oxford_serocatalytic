# Run these two lines to install serofoi package once and then comment it out
# if(!require("remotes")) install.packages("remotes")
# remotes::install_github("epiverse-trace/serofoi")

library(serofoi)
library(ggplot2)
library(dplyr)
library(tidyr)

# Constant FoI
max_age <- 80
foi = 0.05
foi_constant <- data.frame(
  age = seq(1, max_age, 1),
  foi = rep(foi, max_age)
)

ggplot(foi_constant, aes(x = age, y = foi)) +
  geom_line() +
  theme_bw()

# Assume all individual one-year group ages are sampled
# and the sample size in each age is the same: n = 15
n_sample <- 15
survey_features <- data.frame(
  age_min = seq(1, max_age, 1),
  age_max = seq(1, max_age, 1),
  n_sample = rep(n_sample, max_age))

################################################################################
## From serofoi: using the existing function in the serofoi package
################################################################################
# Simulate serosurvey (without covariate)
serosurvey_constant <- simulate_serosurvey(
  model = "age",
  foi = foi_constant,
  survey_features = survey_features
) |>
  mutate(model_type = "constant FoI")


################################################################################
# Manual: Calculate the probability of being seropositive manually
################################################################################
ProbSeropos <- function(age, foi) {
  1 - exp(-foi * age)
}

age <- (survey_features$age_min + survey_features$age_max)/2
prob_seropos <- c()
Y_true_serostatus <- c()

for (i in 1:nrow(survey_features)) {
  prob_seropos[i] <- ProbSeropos(age[i], foi)
  # Y(a) ~ Bernoulli(X(a))
  Y_true_serostatus[i] <- rbinom(1, n_sample, prob_seropos[i])
} 

serosurvey_constant_manual <- data.frame(
  age_min = survey_features$age_min,
  age_max = survey_features$age_max,
  n_sample = survey_features$n_sample,
  n_seropositive = Y_true_serostatus,
  model_type = "Constant FOI"
)

write.csv(serosurvey_constant, "data/raw/serodata_no_sex.csv", row.names = FALSE)
