library(rstan)
library(ggplot2)
# set.seed(123)

# Simulate data for simple linear regression
N <- 100  # Number of observations
x1 <- rnorm(N, mean = 5, sd = 2)
true_intercept <- 2
true_slope1 <- 3
sigma <- 1  # Noise standard deviation
y <- true_intercept + true_slope1 * x1 + rnorm(N, mean = 0, sd = sigma)
hist(x1)
plot(x1, y)

data <- list(N = N, x = x1, y = y)

# Classical (or frequentist) linear regression
m1 <- lm(y ~ x1)
summary(m1)

# Stan model for simple linear regression
stanmodel <- stan_model("src/stan/simulation_eg_stan.stan")

# Compile and fit the model
fit <- sampling(stanmodel, data = data, 
                iter = 2000, chains = 4)
print(fit, c("alpha", "beta", "sigma"))

# ------------------------------------------------------------------------------

# Case with collinearity
x2 <- x1 + rnorm(N, mean = 0, sd = 0.01)  # Highly correlated with x1
true_intercept <- 2 # same as before
true_slope1 <- 3 # same as before
true_slope2 <- 10 # new slope for new variable
sigma <- 1  # Noise standard deviation, same as before
y <- true_intercept + true_slope1 * x1 + true_slope2 * x2 + rnorm(N, mean = 0, sd = sigma)

data_collinear <- list(N = N, x1 = x1, x2 = x2, y = y)
hist(x2)
plot(x1, y)
points(x2, y, col= "blue")
plot(x1, x2)

# Classical (or frequentist) linear regression
m2 <- lm(y ~ x1 + x2)
summary(m2) # do you see anything surprising/unstable? hint- estimate, sd

stanmodel_collinear <- stan_model("src/stan/simulation_coll_eg_stan.stan")

# Fit collinear model
fit_collinear <- sampling(stanmodel_collinear, data = data_collinear, 
                          iter = 2000, chains = 4)
print(fit_collinear, pars = c("alpha", "beta1", "beta2", "sigma"))

# Extract posterior predictive samples
y_rep <- rstan::extract(fit_collinear, pars = "y_rep")$y_rep  # matrix:iterations × N

# Plot observed vs posterior predictive
y_pred_summary <- apply(y_rep, 2, function(x){
  c(median = median(x), 
    lower = quantile(x, 0.025), 
    upper = quantile(x, 0.975))
})

# Convert to a data frame for plotting
y_pred_df <- as.data.frame(t(y_pred_summary))
y_pred_df$observed <- y

# Why do the posterior predictive values look alright? hint- combination ok
ggplot(y_pred_df, aes(x = x1)) +
  geom_ribbon(aes(ymin = `lower.2.5%`, ymax = `upper.97.5%`), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = median), color = "blue") +
  geom_point(aes(y = observed), color = "black") +
  labs(title = "Posterior Predictive Check",
       x = "Observation",
       y = "Outcome") +
  theme_minimal()

ggplot(y_pred_df, aes(x = x2)) +
  geom_ribbon(aes(ymin = `lower.2.5%`, ymax = `upper.97.5%`), fill = "lightgreen", alpha = 0.5) +
  geom_line(aes(y = median), color = "darkgreen") +
  geom_point(aes(y = observed), color = "black") +
  labs(title = "Posterior Predictive Check",
       x = "Observation",
       y = "Outcome") +
  theme_minimal()

# Do the same but with the informative model
stanmodel_collinear_inf <- stan_model("src/stan/simulation_coll_eg_stan_informative.stan")

# Fit collinear model
fit_collinear_inf <- sampling(stanmodel_collinear_inf, data = data_collinear, 
                              iter = 2000, chains = 4)
print(fit_collinear_inf, pars = c("alpha", "beta1", "beta2", "sigma"))



