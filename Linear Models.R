
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


# Exploring non normality of sigma

d3 <- sample(d2$height, 20)

mu.list <- seq(from = 150, to = 170, length.out = 200)
sigma.list  <- seq(from = 4, to = 20, length.out = 200)
post2 <- expand.grid(mu = mu.list, sigma = sigma.list)
post2$LL <- sapply(1:nrow(post2), function(i)
                  sum(dnorm(d3, mean = post2$mu[i], sd = post2$sigma[i], log = TRUE)))
post2$prod <- post2$LL + dnorm(post2$mu, 178, 20, TRUE) + dunif(post2$sigma, 0, 50, TRUE)
post2$prob <- exp(post2$prod-max(post2$prod))

sample2.rows <- sample(1:nrow(post2), 1e4, replace = TRUE, prob = post2$prob)
sample2.mu <- post2$mu[sample2.rows]
sample2.sigma <- post2$sigma[sample2.rows]

plot(sample2.mu, sample2.sigma, cex = 0.5, col = rgb(0.5, 0.2, 0.5, alpha = 0.1), xlab = "mu",
     ylab = "sigma", pch = 16)

dens(sample2.sigma, norm.comp = T)


# Quadratic approximation 

# using MAP command from rethinking package. Needs the ingredients from model definition:

flist <- alist(
         height ~ dnorm(mu, sigma),
         mu ~ dnorm(178, 20),
         sigma ~ dunif(0, 50)
        )
m4.1 <- map(flist, data = d2)    # data fitting to model
precis(m4.1)                     # look at map


# More informative prior:

m4.2 <- map(
            alist(
                  height ~ dnorm(mu, sigma),
                  mu ~ dnorm(178, 0.1),
                  sigma ~ dunif(0, 50)
            ),
            data = d2
            )

precis(m4.2)
precis(m4.1)

# Sampling from MAP

vcov(m4.1) # Variance covariance matrix

diag(vcov(m4.1))
cov2cor(vcov(m4.1))

# Sample vectors from a multidimensional gaussian instead of just single values from gaussian. 

post <- extract.samples(m4.1, n = 1e4)  # Uses mvnorm under the hood
head(post)

precis(post)
plot(post, pch = 16, col = rgb(0.5, 0.2, 0.5, alpha = 0.1), cex = 0.5)

# Adding a predictor!

plot(d2$height, d2$weight)  # Seem to be related

m4.3 <- map(alist(
                  height ~ dnorm(mu, sigma),
                  mu <- a + b*weight,
                  a ~ dnorm(156, 100),
                  b ~ dnorm(0, 10),
                  sigma ~ dunif(0, 50)
                  ), data = d2 )

# Table of estimates

precis(m4.3)

cov2cor(vcov(m4.3))
precis(m4.3, corr = TRUE)    # to see how parameters are correlated

# Large correlations like this can be problematic in larger models. 
# so we can solve this by: x. Centering

d2$weight.c <- d2$weight - mean(d2$weight)   # centers data around 0
mean(d2$weight.c)   # should be close to zero

# refitting model

m4.4 <- map(alist( height ~ dnorm(mu, sigma),
                   mu <- a + b*weight.c,
                   a ~ dnorm(178, 100),
                   b ~ dnorm(0, 10),
                   sigma ~ dunif(0, 50)
                  ), data = d2)
precis(m4.4)
cov2cor(vcov(m4.4))
precis(m4.4, corr = TRUE)

# Plotting!

plot(height ~ weight, data = d2, cex = 0.5, pch = 16, col = "cyan4")
abline(a = coef(m4.3)["a"], b = coef(m4.3)["b"], lty = 2)

# above plot doesn't coomunicate uncertainty very well. 
# many different combinations of a and b could be likely, or only a few could
# the above graph does not communicate which.

post <- extract.samples(m4.3)

post[1:5,]

# Lets only use some of the data to ease learning

N <- 100    # change N to see how different amount of data results in different levels of spread in our regression lines
dN <- d2[1:N, ]
mN <- map(alist(
                height ~ dnorm(mu, sigma),
                mu <- a + b*weight,
                a ~ dnorm(178, 100),
                b ~ dnorm(0, 10),
                sigma ~ dunif(0, 50))
          , data = dN)

# extract 20 samples from posterior
post <- extract.samples(mN, n = 20)

# display raw data and sample size
plot(dN$weight, dN$height, 
     xlim = range(dN$weight), ylim = range(dN$height), 
     col = rangi2, xlab = "Weight", ylab = "Height",
     pch = 16)
mtext(concat("N = ", N))

# Plot in sample lines
for (i in 1:20){
  abline(a = post$a[i], b = post$b[i], col = col.alpha("black", 0.3))}

par(mfrow = c(3, 3))
par(mfrow = c(1, 1))

# Regression intervals and contours to show uncertainty. 

# Look at one single weight value

mu_at_50 <- post$a + post$b * 50
dens(mu_at_50, xlab = "mu ´weight = 50")    

HPDI(mu_at_50)

# This need to be computed for every weight value in the model.

mu <- link(m4.3)    # Provides 1000 samples of mean height for each case in the data, i.e. person
                    # each column is a different persons
str(mu)

# But we dont want it done for each person but rather each weight.
# so we define sequence of weights to compute predictions for
# these values will be each column¨

weight.seq <- seq(from = 25, to = 70, by = 1)

# use link to compute mu for each sample from posterior and for each weight in weight.seq

mu <- link(m4.3, data = data.frame(weight = weight.seq) )
str(mu)

colnames(mu) <- weight.seq
# use type = "n" to hide raw data
plot(height ~ weight, d2, type = "n")

for(i in 1:100){
  points(weight.seq, mu[i, ], pch = 16, col = col.alpha(rangi2, 0.1))
}

# summarize distribution of mu


mu.mean <- apply(mu, 2, mean)             # read as compute the mean of each column (dimensions 2) in matrix mu
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.89)

# plot raw data 
plot(height ~ weight, data = d2, col = col.alpha(rangi2, 0.4), pch = 16)

# plot MAP line, mean mu for each weight
lines(weight.seq, mu.mean)

# plot shaded region for 89% HPDI interval
shade(mu.HPDI, weight.seq)


# Prediction intervals
# So far we've only expressed our uncertainty in the mean. Not our uncertainty in height
# which also incorporates the standard deviation

# Lets simulate heights drawn from a Gaussian distribution with the correct mean height
# and standard deviation from a specific weight

sim.height <- sim(m4.3, data = list(weight = weight.seq))
str(sim.height)

# Summarize

height.PI <- apply(sim.height, 2, PI, prob = 0.89)

# all together now

plot(height ~ weight, d2, col = col.alpha(rangi2, 0.3), cex = 1, pch = 16)
# MAP line
lines(weight.seq, mu.mean)
# HPDI region of mu
shade(mu.HPDI, weight.seq) 
# PI region for simulated heights
shade(height.PI, weight.seq)

# Polynomial Regression

str(d)

plot(height ~ weight, data = d)  # Its curved when the yung uns aren't removed

# standardize predictor (weight)

d$weight.s <- (d$weight - mean(d$weight))/sd(d$weight) # Makes it easier to compare influence of predictors
                                                       # and also solves possible numerical glitches

                                                       # sets weight.s mean to zero and sd to 1, no info lost
par(mfrow = c(1, 2))
plot(height ~ weight.s, data = d)                      # same shape but different scale on x axis.


# Model

d$weight.s2 <- d$weight.s^2

m4.5 <- map(alist(
                  height ~ dnorm(mu, sigma),
                  mu <-  a + b1 * weight.s + b2 * weight.s2,
                  a ~ dnorm(178, 100),
                  b1 ~ dnorm(0, 10),
                  b2 ~ dnorm(0, 10), 
                  sigma ~ dunif(0, 50)),
            data = d)

precis(m4.5) # hard to interpret

# Plotting

weight.seq <- seq(from = -2.2, to = 2, length.out = 30)
pred_dat <- list(weight.s = weight.seq, weight.s2 = weight.seq^2)
mu <- link(m4.5, data = pred_dat)
colnames(mu) <- weight.seq
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = 0.89)
sim.height <- sim(m4.5, data = pred_dat)
height.PI <- apply(sim.height, 2, PI, prob = 0.89)

par(mfrow = c(1,1))

plot(height ~ weight.s, data = d, col = col.alpha(rangi2, 0.2), pch = 16)
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)



# Practice

# 4M 
# Simulate heights from prior

sim <- function(N){
                    sigma <- runif(N, 0, 10)
                    mu <- rnorm(N, 0, 10)
                    H <- rnorm(N, mu, sigma)
                    return(H)
}

sim(100)

# translate into map

M42 <- map(alist(
                  h ~ dnorm(mu, sigma)
                  mu ~ dnorm(0, 10)
                  sigma ~ dunif(0, 10)
                  ), data = X)


# 4H

# Simulate heights for specific weights

miss.data <- c(46.95, 43.72, 64.78, 32.59, 54.63)
miss.data <- (miss.data - mean(d$weight))/sd(d$weight)
miss.data.2 <- miss.data^2

# use link to compute mu for each sample from posterior and for each weight in weight.seq

mu <- link(m4.5, data = data.frame(weight.s = miss.data, weight.s2 = miss.data.2))
str(mu)

colnames(mu) <- miss.data

expect.2 <- apply(mu, 2, mean)
interval.2 <- apply(mu, 2, PI, prob = 0.89)

expect
expect.2
interval
interval.2

# Only yung uns

d4 <- d[d$age<18,]

M42 <- map(alist(
                 height ~ dnorm(mu, sigma),
                 mu <-  a + b * weight,
                 a ~ dnorm(170, 20),
                 b ~ dnorm(0, 10),
                 sigma ~ dunif(0, 50)
                 ), data = d4)

precis(M42)

plot(height ~ weight, data = d4, col = rgb(0.3, 0.1, 0.5, alpha = 0.2), pch = 16)
abline(a = 58, b = 2.7)

range(d4$weight)

weight.grid <- seq(from = 4, to = 45, by = 0.5)
sample.means.per.weight <- link(M42, data = data.frame(weight = weight.grid))
colnames(sample.means.per.weight) <- weight.grid

sample.mu <- apply(sample.means.per.weight, 2, mean)
sample.HPDI <- apply(sample.means.per.weight, 2, HPDI, prob = 0.89)

plot(height ~ weight, data = d4, col = rgb(0.3, 0.1, 0.5, alpha = 0.2), pch = 16, cex = 0.5)
lines(weight.grid, sample.mu)
shade(sample.HPDI, weight.grid)

samples <- extract.samples(M42)

height.interval <- sapply( weight.grid , function(z)
  HPDI( rnorm(10000, samples$a + samples$b*z, post$sigma) , 0.89) )
colnames(height.interval) <- weight.grid

shade(height.interval, weight.grid)


# Lets log it instead!

M.log <- map(alist(
                   height ~ dnorm(mu, sigma),
                   mu <- a + b * log(weight),
                   a ~ dnorm(178, 100),
                   b ~ dnorm(0, 100),
                   sigma ~ dunif(0, 50)
                   ), data = d)

plot(height ~ weight, data = d,
     col = col.alpha(rangi2, 0.4))
range(d$weight)
weight.grid <- seq(4, 63, by = 0.5)

samples <- extract.samples(M.log)
sample.means <- sapply(weight.grid, function(z)
  mean(samples$a + samples$b*log(z)))
sample.means.CI <- sapply(weight.grid, function(z)
  HPDI(samples$a + samples$b*log(z)))
sample.height.CI <- sapply(weight.grid, function(z)
  HPDI(rnorm(1e4, samples$a + samples$b*log(z), samples$sigma)))

lines(weight.grid, sample.means)
shade(sample.means.CI, weight.grid)
shade(sample.height.CI, weight.grid)
