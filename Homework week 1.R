
# Homework week 1
# RAG experiment dice roll 123: loss, 456: win. 171 participants, 111 claim prizes. 
# Compute posterior distribution for the proportion of participants who were honest

library(rethinking)

# (1) Grid, (2) Prior, (3) Likelihood, (4) Posterior, (5) Stand. posterior
p_grid <- seq(0, 1, length.out = 1000)

prior <- rep(1, 1000)

likelihood <- dbinom(111, 171, prob = p_grid)

posterior <- likelihood*prior

posterior <- posterior/sum(posterior)

plot(p_grid, posterior, type = "l")  # Percent people who claimed a prize
abline(v = 0.5, col = "hotpink")     # Expected number of people to get a prize
# Clearly a few liars
# Since the expected number of prizes collected is 0.5, the proportion of
#liars can be observed by taking (p_grid - 0.5)/0.5, so prportion honest is simply
# 1 - liars
p_grid_honest <- (1-(p_grid-0.5)/0.5)
plot(p_grid_honest, posterior, type = "l", xlim = c(0,1))
# Does this account for observational uncertainty?

# Checking method 1
# Creating samples

samples <- sample(p_grid, 1e4, replace = T, prob =posterior)
samples <- samples*171
dens(samples)

# We know the probability of winning is 0.5 percent
expected <- rbinom(1e4, 171, prob = 0.5)
dens(expected)

honest <- expected/samples    # Proportion of honest to liars

par(mfrow = c(1, 2))
dens(honest, xlim = c(0, 1.1))
plot(p_grid_honest, posterior, type = "l", xlim = c(0, 1.1))
# seems legit

par(mfrow = c(1, 1))


# Compute posterior predictive distribution for the next 10 participants

samples <- sample(p_grid, 1e4, replace = T, prob = posterior)

posterior_pred_dist <- rbinom(1e4, 10, prob = samples)

simplehist(posterior_pred_dist, xlab = "Predicted number of prize claimers")


