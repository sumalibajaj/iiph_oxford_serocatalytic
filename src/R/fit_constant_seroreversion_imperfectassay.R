library(dplyr)
library(ggplot2)
library(rstan)

# Read the simulated serosurvey
serosurvey <- read.csv("data/raw/serodata_constant_seroreversion_imperfectassay.csv")

# Add observed seroprevalence and uncertainty in the data (frequentist and Bayesian)
# Get data ready for stan model
serosurvey <- serosurvey %>%
  mutate(
    prev = n_seropositive/n_sample,
    sd = sqrt(prev*(1-prev)/n_sample)) %>% 
  mutate(
    a = 1 + n_seropositive,
    b = 1 + n_sample - n_seropositive) %>% 
  mutate(
    lower=qbeta(0.025, a, b),
    middle=qbeta(0.5, a, b),
    upper=qbeta(0.975, a, b)
  ) %>%
  mutate(age_group = floor(age_min + age_max)/2)

# Plot serosurvey
ggplot() +
  geom_pointrange(data = serosurvey, aes(x = (age_min + age_max)/2,
                                         ymin=prev-1.96*sd, y=prev, ymax=prev+1.96*sd),
                  colour="black") +
  labs(title = "Observed estimates of seroprevalence",
       x = "Age",
       y = "Seroprevalence") +
  theme_bw()

ggplot() +
  geom_pointrange(data = serosurvey, aes(x = (age_min + age_max)/2,
                                         ymin=lower, y=middle, ymax=upper),
                  colour="black") +
  labs(title = "Observed estimates of seroprevalence",
       x = "Age",
       y = "Seroprevalence") +
  theme_bw()


# to plot posterior sero prevalence as a curve and compare to real data
ages_fine <- seq(0, max(serosurvey$age_max), by = 0.1)  # or your desired range
n_ages_fine <- length(ages_fine)

# ------------------------------------------------------------------------------
# LETS TRY FITTING A CONSTANT FOI ONLY MODEL
# ------------------------------------------------------------------------------

stan_data <- list(n_observations = nrow(serosurvey),
                  n_seropositive = serosurvey$n_seropositive,
                  n_sample = serosurvey$n_sample,
                  age_group = serosurvey$age_group,
                  age_max = max(serosurvey$age_max),
                  ages = seq(1:max(serosurvey$age_max)),
                  ages_fine = ages_fine,
                  n_ages_fine = n_ages_fine)


# compile the stan model for constant foi
stanmodel <- stan_model("src/stan/constant.stan")

# Compile and fit the model
fit <- sampling(stanmodel, data = stan_data,
                iter = 2000, chains = 4)
print(fit, pars = "foi")

# View the posterior distribution of foi
hist(rstan::extract(fit, "foi")[[1]])

# Extract posterior predictive samples of seroprevalence for the fine ages
seroprev_fine <- rstan::extract(fit, pars = "seroprev_fine")$seroprev_fine  # matrix:iterations × N

# Plot observed vs posterior predictive
seroprev_fine_summary <- apply(seroprev_fine, 2, function(x){
  c(median = median(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
})

# Convert to a data frame for plotting
seroprev_fine_df <- as.data.frame(t(seroprev_fine_summary))
seroprev_fine_df$age <- ages_fine

# Plot and compare observed and posterior seroprevalence estimates
ggplot() +
  geom_pointrange(data = serosurvey, aes(x = age_group,
                                         ymin=lower, y=middle, ymax=upper),
                  colour="black") +
  geom_ribbon(data = seroprev_fine_df,
              aes(x = age, ymin = `lower.2.5%`, ymax = `upper.97.5%`),
              fill = "lightblue", alpha = 0.5) +
  geom_line(data = seroprev_fine_df,
            aes(x = age, y = median),
            color = "blue") +
  labs(title = "Observed and posterior estimates of seroprevalence",
       x = "Age",
       y = "Seroprevalence") +
  theme_minimal()

# ------------------------------------------------------------------------------
# With this model we estimated a much lower FOI, and the posterior fit is also poor
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# NOW LETS FIT A CONSTANT FOI WITH SEROREVERSION AND AN IMPERFECT ASSAY
# ------------------------------------------------------------------------------

rate_seroreversion <- 1/3
sensitivity <- 0.84
specificity <- 0.90

# Prepare data for Stan
stan_data <- list(
  n_observations = nrow(serosurvey),
  n_seropositive = serosurvey$n_seropositive,
  n_sample = serosurvey$n_sample,
  age_group = serosurvey$age_group,
  rate_seroreversion = rate_seroreversion,
  sensitivity = sensitivity,
  specificity = specificity,
  age_max = max(serosurvey$age_max),
  ages = seq(1, max(serosurvey$age_max)),
  ages_fine = ages_fine,
  n_ages_fine = n_ages_fine
)

# Compile and fit the updated Stan model
stanmodel <- stan_model("src/stan/constant_seroreversion_imperfectassay.stan")
fit <- sampling(stanmodel, data = stan_data, iter = 2000, chains = 4)

# Print posterior summaries
print(fit, pars = "foi")
hist(rstan::extract(fit, "foi")[[1]])

# Extract posterior predictive samples of seroprevalence for the fine ages
seroprev_fine <- rstan::extract(fit, pars = "seroprev_fine")$seroprev_fine

# Summarize posterior for plotting
seroprev_fine_summary <- apply(seroprev_fine, 2, function(x) {
  c(median = median(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
})
seroprev_fine_df <- as.data.frame(t(seroprev_fine_summary))
seroprev_fine_df$age <- ages_fine

# Plot observed vs posterior predictive seroprevalence
ggplot() +
  geom_pointrange(data = serosurvey, aes(x = age_group, ymin = lower, y = middle, ymax = upper),
                  colour = "black") +
  geom_ribbon(data = seroprev_fine_df,
              aes(x = age, ymin = `lower.2.5%`, ymax = `upper.97.5%`),
              fill = "lightblue", alpha = 0.5) +
  geom_line(data = seroprev_fine_df,
            aes(x = age, y = median),
            color = "blue") +
  labs(title = "Observed and posterior estimates of seroprevalence",
       x = "Age",
       y = "Seroprevalence") +
  theme_minimal()
