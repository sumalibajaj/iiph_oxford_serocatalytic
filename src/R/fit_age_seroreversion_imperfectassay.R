library(dplyr)
library(ggplot2)
library(rstan)

# Read the simulated serosurvey
serosurvey <- read.csv("data/raw/serodata_age_seroreversion_imperfectassay.csv")

# Add observed seroprevalence and uncertainty in the data (frequentist and Bayesian)
# Get data ready for stan model
serosurvey <- serosurvey %>%
  mutate(
    prev = n_seropositive / n_sample,
    sd = sqrt(prev * (1 - prev) / n_sample),
    a = 1 + n_seropositive,
    b = 1 + n_sample - n_seropositive,
    lower = qbeta(0.025, a, b),
    middle = qbeta(0.5, a, b),
    upper = qbeta(0.975, a, b),
    age_group = floor((age_min + age_max) / 2)
  )

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

# Create FOI indices - you cab modify this as you wish- kept it as equal chunks
max_age <- 80
window_size = 5
foi_index <- rep(1:(max_age/window_size), each = window_size)
n_foi <- length(unique(foi_index))

# Seroreversion and assay parameters
rate_seroreversion <- 1/3
sensitivity <- 0.84
specificity <- 0.90

# Define fine resolution ages for plotting
ages_fine <- seq(0, max(serosurvey$age_max), by = 0.1)
n_ages_fine <- length(ages_fine)

# Prepare data for Stan
stan_data <- list(
  n_observations = nrow(serosurvey),
  n_seropositive = serosurvey$n_seropositive,
  n_sample = serosurvey$n_sample,
  age_group = serosurvey$age_group,
  rate_seroreversion = rate_seroreversion,
  sensitivity = sensitivity,
  specificity = specificity,
  age_max = max(serosurvey$age_group),
  ages = seq(1, max(serosurvey$age_group)),
  ages_fine = ages_fine,
  n_ages_fine = n_ages_fine,
  foi_index = foi_index
)

# Compile and fit the updated Stan model
stanmodel <- stan_model("src/stan/age_seroreversion_imperfectassay.stan")
fit <- sampling(stanmodel, data = stan_data, iter = 2000, chains = 4)

print(fit, pars = "foi_vector")

# Extract posterior samples of foi_vector (dimension: iterations x n_foi)
foi_samples <- rstan::extract(fit, pars = "foi_vector")$foi_vector
# foi_samples is a matrix: rows = posterior draws, columns = time intervals

# foi_index: vector of length n_years, each value in 1:4
n_iter <- nrow(foi_samples)
n_years <- length(foi_index)

# Create a matrix to hold expanded FOI for all years & iterations
expanded_foi_samples <- matrix(NA, nrow = n_iter, ncol = n_years)
for (iter in 1:n_iter) {
  # For iteration i, get the 4 FOI values
  foi_iter <- foi_samples[iter, ]
  # Map each year to FOI chunk using foi_index
  expanded_foi_samples[iter, ] <- foi_iter[foi_index]
}

# Now compute summaries for each calendar year
foi_summary <- data.frame(
  age = seq(1, max_age),
  mean = apply(expanded_foi_samples, 2, mean),
  lower = apply(expanded_foi_samples, 2, quantile, probs = 0.025),
  upper = apply(expanded_foi_samples, 2, quantile, probs = 0.975)
)

# Read the time varying foi used to simulate the serosurvey
foi_obs <- read.csv("data/raw/foi_age_seroreversion_imperfectassay.csv")

ggplot(foi_summary, aes(x = age)) +
  geom_line(aes(y = mean), color = "#fca311") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#fca311") +
  geom_line(data = foi_obs, aes(x = age, y = foi)) +
  labs(title = "Observed (black) and posterior age-varying FOI",
       x = "Age",
       y = "FOI") +
  ylim(c(0,0.3)) +
  theme_minimal()

# Extract posterior predictive samples of seroprevalence
seroprev_fine <- rstan::extract(fit, pars = "seroprevalence")$seroprevalence  # matrix:iterations × N

# Plot observed vs posterior predictive
seroprev_fine_summary <- apply(seroprev_fine, 2, function(x){
  c(median = median(x), 
    lower = quantile(x, 0.025), 
    upper = quantile(x, 0.975))
})

# Convert to a data frame for plotting
seroprev_fine_df <- as.data.frame(t(seroprev_fine_summary))
seroprev_fine_df$age <- serosurvey$age_group

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


