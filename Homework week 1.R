
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


# Solutions from McElreath

compute_posterior <- function(Y, N, grid){
  ways <- sapply(grid, function(q) q^Y * (1 - q)^N)
  post <- ways/sum(ways)
  data.frame(grid, ways, post = round(post, 3))
}

# Posterior probability for not claiming prize
posterior <- compute_posterior(171-111, 111, seq(0, 1, length.out = 100))
plot(posterior$grid, posterior$post, type = "l")

# or

posterior <- dbeta(seq(0, 1, length.out = 100), 171-111+1, 111+1)   
plot(seq(1, 0, length.out = 100), posterior, type = "l")

# Find cheaters

compute_posterior <- function(Y, N, grid){
  ways <- sapply(grid, function(q) 
    {p <- q + (1 - q)*0.5; return(p^Y * (1 - p)^N)}) # p is proportion who claim prize
  post <- ways/sum(ways)
  data.frame(grid, ways, post = round(post, 3))
}


posterior2 <- compute_posterior(111, 171-111, seq(0, 1, len = 100))
plot(posterior2$grid, posterior2$post, type = "l")


p_samples <- rbeta(1e4, 111, 171-111)
# next ten rolls
Y_sim <- rbinom(1e4, size = 10, p = p_samples)

plot(table(Y_sim))
