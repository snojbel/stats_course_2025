# Sampling the imaginary

library(cmdstanr)
library(rethinking)

# Grid approximation from globe
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(x = 6, size = 9, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)

plot(posterior)

# Sampling from a grid-approximate posterior. 

samples <- sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)

dens(samples)

# Sampling to summarize

# Intervals of defined boundries

sum(posterior[p_grid<0.5]) # sums the posterior for all values below 0.5. The probability of below 0.5
# This does not work well when many parameter values come into play -> solution is sampling

sum(samples < 0.5)/1e4   # Divide by number of samples
sum(samples > 0.5 & samples < 0.75)/1e4

# Intervals of defined probability mass

quantile(samples, 0.8) # lower 80 percent of data
quantile(samples, c(0.1, 0.9)) # Shows the middle 80 percent quantile. 

# Three waters three tosses

p_grid <- seq(0, 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(3, 3, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)

plot(p_grid, posterior)  # Highly skewed

samples <- sample(p_grid, 1e4, replace = T, prob = posterior)

quantile(samples, c(0.1, 0.9))

dens(samples)
abline(v = 0.5665666)
abline(v = 0.9739)

PI(samples, prob = 0.5)   #middle 50% of data
HPDI(samples, prob = 0.5)

# Point estimates

# Maximum a posteriori (MAP)

p_grid[which.max(posterior)]

chainmode(samples, adj = 0.01) # from samples

# Mean and median?
mean(samples)
median(samples)

# Loss function

sum(posterior*abs(0.5-p_grid))   # Loss at parameter value 0.5

loss <- sapply(p_grid, function(d) sum(posterior*abs(d-p_grid))) # Calculating over all parameter values

plot(p_grid, loss)

p_grid[which.min(loss)]  # Value that minimizes loss function.

# Different loss functions will lead to different point estimates. This was the absolute loss function
# and will result in the median. 

# Sampling to simulate prediction

# Dummy data for globe model

dbinom(0:2, size = 2, prob = 0.7) # asks how probable it is to see 0, 1, 2 tosses be water if 0.7 of globe
                                  # is water and 2 tosses are made.

# Data simulation via rbinom, random draw from binomial distribution
rbinom(1, size = 2, prob = 0.7)

dummy_data <- rbinom(1e5, size = 20, prob = 0.7)
table(dummy_data)/1e5
20*0.7
simplehist(dummy_data, xlab = "Dummy water count")

# We are uncertain about both our observations and the true value of the parameter 
# -> Posterior Predicitive Distribution

w <- rbinom(1e4, size = 9, prob = samples)

simplehist(w)

# Practice

P_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep (1, 1000)
likelihood <- dbinom(6, 9, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
set.seed(100)

samples <- sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)

# 3E1-7

# Intervals of defined boundries
sum(posterior[p_grid<0.2]) #or
sum(samples < 0.2) / 1e4

sum(posterior[p_grid>0.8])
sum(samples>0.8) / 1e4

sum(posterior[p_grid>0.2 & p_grid<0.8])
sum(samples>0.2 & samples <0.8) / 1e4

# Intervals of defined mass

quantile(samples, prob = 0.2)


quantile(samples, prob = 0.8)

interval1 <- HPDI(samples, prob = 0.66)

quantile(samples, c(0.17, 1-0.17))
interval2 <- PI(samples, prob = 0.66)

width1 <- interval1[2] - interval1[1]
width2 <- interval2[2] - interval2[1]

cbind(width1, width2)

# 3M 1-5

p_grid <- seq(0, 1, length.out(1000))
prior <- rep(1, 1000)
likelihood <- dbinom(8, 15, prob = p_grid)
posterior <- prior*likelihood
posterior <- posterior/sum(posterior)

plot(p_grid, posterior)

samples <- sample(p_grid, size = 1e4, prob = posterior, replace = TRUE)

interval <- HPDI(samples, prob = 0.9)
abline(v = interval)


Post_pred_dist <- rbinom(1e4, size = 15, prob = samples)
simplehist(Post_pred_dist)

sum(Post_pred_dist==8)/1e4

check <- rbinom(1e4, size = 9, prob = samples)
simplehist(check)
sum(check==6)/1e4

# New prior

p_grid <- seq(0, 1, length.out(1000))
prior <- c(rep(0, 1000/2), rep(1, 1000/2))
likelihood <- dbinom(8, 15, prob = p_grid)
posterior <- prior*likelihood
posterior <- posterior/sum(posterior)

plot(p_grid, posterior)

samples <- sample(p_grid, size = 1e4, prob = posterior, replace = TRUE)

interval <- HPDI(samples, prob = 0.9)
abline(v = interval)


Post_pred_dist <- rbinom(1e4, size = 15, prob = samples)
simplehist(Post_pred_dist)

sum(Post_pred_dist==8)/1e4

check <- rbinom(1e4, size = 9, prob = samples)
simplehist(check)
sum(check==6)/1e4


# 3H 1-5
data(homeworkch3)

trials <- 2 * length(birth1)
boys <- sum(birth1) + sum(birth2)

p_grid <- seq(0, 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(boys, trials, prob = p_grid)
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)

plot(p_grid, posterior, type = "l")
abline(v = 0.5, lty = 2)
p_grid[which.max(posterior)]

# Intervals

samples <- sample(p_grid, 1e4, replace = TRUE, prob = posterior)

HPDI(samples, prob = 0.5)
HPDI(samples, prob = 0.89)
HPDI(samples, prob = 0.9)


check <- rbinom(1e4, size = 200, prob = samples)
par(mfrow=c(1,1))
simplehist(check)
abline(v = 111, col = "hotpink")
dens(check)

par(mfrow=c(1,2))
dens(check/200)
plot(p_grid, posterior, type = "l")

check2 <- rbinom(1e4, size = 100, prob = samples)
dens(check2)
abline(v = sum(birth1), col = "hotpink")

# Checking independence. 

number <- sum(birth1==0)

check3 <- rbinom(1e4, size = 49, prob = samples)

simplehist(check3)
sum(birth2[birth1==0])
abline(v = sum(birth2[which(birth1==0)]), col = "hotpink")
