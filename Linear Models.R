
# Linear Models

library(rethinking)

# Imagine 16 people on a football fields midline, they flip a coin and each time 
# they walk either a little backwards or a little forwards depending on the flip
# whats the distance traveled by each individual from the midline

pos <- replicate(1000, sum(runif(16, min = -1, max = 1)))

hist(pos)
plot(density(pos))

# Dozen loci with multiplicative but small and different change in growth

prod(1 + runif(12, 0, 0.1)) # takes the product of 12 numbers between 1 and 1.1

growth <- replicate(10000, prod(1 + runif(12, 0, 0.1)))

plot(density(growth))
hist(growth)
dens(growth, col = "cyan4", norm.comp = T)

# This becomes normal because the multiplicate effects are all small!

big <- replicate(10000, prod(1 + runif(12, 0, 0.5)))
small <- replicate(10000, prod(1 + runif(12, 0, 0.01)))
plot(density(big))
plot(density(small))

# So large deviates multiplied are not normal, but the logged values are

log.big <- replicate(1000, log(prod(1 + runif(12, 0, 0.5))))
dens(log.big, col ="mediumorchid", norm.comp = T)


# A Gaussian model for height

data(Howell1) # Real data from rethinking package
d <- Howell1
str(d) # display structure of r object

d2 <- d[d$age>=18, ]  # Remove any non adults 
str(d2)

dens(d2$height)

# We shall assume height is normally distributed. 

# Priors

# Mean prior: Normal(178, 20)
curve(dnorm(x, 178, 20), from = 100, to = 250)
  # Carries some information, but not alot. Assumes that most people are between 140 and 220. 

# Standard deviation prior
curve(dunif(x, min = 0, max = 50), from = 0, to = 60)

# Sampling from prior mu and sigma to check what they mean

sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)

dens(prior_h)

# Grid approximation
mu.list <- seq(from = 140, to = 160, length.out = 200)
sigma.list <- seq(from = 4, to = 9, length.out = 200)
post <- expand.grid(mu = mu.list, sigma = sigma.list) # Creates a dataframe with all possible 
                                                      # combinations of the two variables
post$LL <- sapply(1:nrow(post), function(i) sum(dnorm(
                  d2$height,
                  mean = post$mu[i],
                  sd = post$sigma[i],
                  log = TRUE
)))

post$prod <- post$LL + dnorm(post$mu, 178, 20, TRUE) + dunif(post$sigma, 0, 50, TRUE)
post$prob <- exp(post$prod - max(post$prod))

contour_xyz(post$mu, post$sigma, post$prob, xlim = c(150, 160), ylim = c(6, 9))
image_xyz(post$mu, post$sigma, post$prob)


# Sampling from distribution

sample.rows <- sample(1:nrow(post), 1e4, replace = TRUE, prob = post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]

plot(sample.mu, sample.sigma, cex = 0.5, pch = 16, col = col.alpha(rangi2,0.1))

dens(sample.mu)
dens(sample.sigma)

HPDI(sample.mu)
HPDI(sample.sigma)
