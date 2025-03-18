# Overfittning, Regularization and Information Criteria
# Overfitting and underfitting are two monsters you must always navigate.
# Overfitting is often dealt with using regularizing prior (penalized likelihood)
# Underfitting often uses information criteria and cross validation.

# Models that are good at prediction are not always good causal models, often these
# features must be weiged against each other. 


# All data sets have regular and irregular features, our goal is to learn the regular features. 

library(rethinking)

# Hominin brain size -----------------------------------------------------------

sppnames <- c("afarensis", "africanus", "habilis", "boisei", "rudolfensis", "ergaster", "sapiens")
brainvolcc <- c(438, 452, 612, 521, 752, 871, 1350)
masskg <- c(37.0, 35.5, 34.5, 41.5, 55.5, 61.1, 53.5 )
d <- data.frame(species = sppnames, brain = brainvolcc, mass = masskg)

# Linear model

m6.1 <- lm(brain ~ mass, data = d)

1-var(resid(m6.1))/var(d$brain)    # R^2 score
summary(m6.1)

# Higher degree polynomials
m6.2 <- lm(brain ~ mass + I(mass^2), data = d)
m6.3 <- lm(brain ~ mass + I(mass^2) + I(mass^3), data = d)
m6.4 <- lm(brain ~ mass + I(mass^2) + I(mass^3) + I(mass^4), data = d)
m6.5 <- lm(brain ~ mass + I(mass^2) + I(mass^3) + I(mass^4) + I(mass^5), data = d)
m6.6 <- lm(brain ~ mass + I(mass^2) + I(mass^3) + I(mass^4) + I(mass^5) + I(mass^6), data = d)

summary(m6.6)

# R^2 increases with each added polynomial, until it reaches 1 at sixth degree polynomial
# because the number of parameters equals the number of data points. 
# The predictions get increasingly more ridiculus and would not be able to accomodate any
# new data point. 

# But one can also underfit models

m6.7 <- lm(brain ~ 1, data = d) # just intercept

# One can check to see how sensitive a model is to the sample

par(mfrow = c(1, 2))

plot(brain ~ mass, data = d, col = "slateblue")

# first order poly
for (i in 1:nrow(d)) {
  d.new <- d[-i, ]  # removes the ith row from data
  m0 <- lm(brain ~ mass, data = d.new)
  abline(m0, col = col.alpha("black", 0.5))
}

# fifth order

plot(brain ~ mass, data = d, col = "slateblue", ylim = c(-500, 2500))

for (i in 1:nrow(d)) {
  d.new <- d[-i, ]  # removes the ith row from data
  m0 <- lm(brain ~ mass + I(mass^2) + I(mass^3) + I(mass^4) + I(mass^5), data = d.new)
  curve(expr = cbind(1, poly(x, degree = 5, raw = TRUE)) %*% m0$coefficients,
        add = TRUE, col = col.alpha("black", 0.5))
}

# Under overfitting dichotomy is often described as "bias-variance trade-off"

# Log-likelihood ---------------------------------------------------------------
# Measures models deviance

m6.1 <- lm(brain ~ mass, data = d)
(-2) * logLik(m6.1)


# WAIC

data(cars)
range(cars$dist)
mean(cars$dist)


m <- quap(alist(dist ~ dnorm(mu, sigma),
                mu <- a + b*speed,
                a ~ dnorm(0, 100),
                b ~ dnorm(0, 10),
                sigma ~ dexp(1)),
          data = cars)

set.seed(94)
post <- extract.samples(m, n = 1e3)

n_samples <- 1000
logprob <- sapply(1:n_samples,
                  function(s){
                    mu <- post$a[s] + post$b[s] * cars$speed
                    dnorm(cars$dist, mu, post$sigma[s], log = TRUE)
                  })
# Obs is rows, samples in columns

# Now to compute lppd, the Bayesian deviance, we average the samples in each row, take
# the log, and add all of the logs together. However, to do this with precision, we need to do all of the
# averaging on the log scale. This is made easy with a function log_sum_exp, which computes the log
# of a sum of exponentiated terms. Then we can just subtract the log of the number of samples. This
# computes the log of the average.

n_cases <- nrow(cars)

lppd <- sapply(1:n_cases, function(i) log_sum_exp(logprob[i,] - log(n_samples)))  # because - log is same as dividing


pWAIC <- sapply(1:n_cases, function(i) var(logprob[i, ])) # Penatly term

-2*(sum(lppd)-sum(pWAIC))


# There will be simulation variance, because of how the samples are drawn from 
# the quap fit. But that variance remains much smaller than the standard error
# of WAIC itself. You can compute the standard error by computing the square root 
# of number of cases multiplied by the variance over the individual observation 
# terms in WAIC:

waic_vec <- -2*( lppd - pWAIC )
sqrt( n_cases*var(waic_vec) )


# Lets look at the marriage divorce example again


data("WaffleDivorce")

d <- WaffleDivorce
d$A <- standardize(d$MedianAgeMarriage)
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)

m5.1 <- quap(alist(D ~ dnorm(mu, sigma),
                   mu <- a + bA * A,
                   a ~ dnorm(0, 0.2),
                   bA ~ dnorm(0, 0.5),
                   sigma ~ dexp(1)),
             data = d)

m5.2 <- quap(alist(D ~ dnorm(mu, sigma),
                   mu <- a + bM * M,
                   a ~ dnorm(0, 0.2),
                   bM ~ dnorm(0, 0.5),
                   sigma ~ dexp(1)),
             data = d)

m5.3 <- quap(alist(D ~ dnorm(mu, sigma),
                   mu <- a + bA * A + bM * M,
                   a ~ dnorm(0, 0.2),
                   bA ~ dnorm(0, 0.5),
                   bM ~ dnorm(0, 0.5),
                   sigma ~ dexp(1)),
             data = d)

precis(m5.3)

set.seed(24071847)
compare( m5.1 , m5.2 , m5.3 , func=PSIS )
plot(compare( m5.1 , m5.2 , m5.3 , func=PSIS ))

# Gives warning that there are some points that have very high K-values They make
# it hard to estimate out-of-sample predictive accuracy. Why? 
# Because any new sample is unlikely to to contain these same outliers

# WAIC can also measure these, but will not warn us about them directly. 


set.seed(24071847)
PSIS_m5.3 <- PSIS(m5.3,pointwise=TRUE)
set.seed(24071847)
WAIC_m5.3 <- WAIC(m5.3,pointwise=TRUE)
plot( PSIS_m5.3$k , WAIC_m5.3$penalty , xlab="PSIS Pareto k" ,
      ylab="WAIC penalty" , col=rangi2 , lwd=2 )

# Robust Regression ----------------------------------------------------------

# Gaussian distributions are easily surprised. Student-t's less so (cause thicker tails)
# (The Student-t distribution arises from a mixture of Gaussian distributions with different variances.)

m5.3t <- quap(alist(D ~ dstudent(2, mu, sigma),   # Uses a t-dis with 2 df
                   mu <- a + bA * A + bM * M,
                   a ~ dnorm(0, 0.2),
                   bA ~ dnorm(0, 0.5),
                   bM ~ dnorm(0, 0.5),
                   sigma ~ dexp(1)),
             data = d)

PSIS(m5.3t)
# No warning

precis(m5.3t)

# Now shows stronger association, because Idaho was pulling in the other direction before

# Practice ---------------------------------------------------------------------

# 7E ------

# 1
# The three motivating criteria for information entropy: Information is the reduction in
# uncertainty when learning an outcome. (1)This uncertainty measure should be continuous.
# and (2) The uncertainty should be bigger with more outcomes. (3) and the uncertainty should
# be additive. 

# 2
# Entropy of a 0.7 weighted coin.:

p <- c( 0.7 , 0.3 )
-sum( p*log(p) )


# 3
# Entropy of: “1” 20%, “2” 25%, ”3” 25%, and ”4” 30% of the time.

p <- c(0.2, 0.25, 0.25, 0.3)
-sum(p*log(p))  

# 4 
# Entropy of dice with 1, 2, 3 equally often.

p <- c(1/3, 1/3, 1/3)
-sum(p*log(p))  

# 7M ----------

# 1
# AIC uses MAP estimates and flat priors whereas WAIC can use entire posterior distribution
# and will provide an estimate of standard error WAIC is more general

# 2 
# Difference between model selection and model comparison? Model selection refers to using 
# information critera to choose between different causal models and infer causal relationships
# this is inappropriate. Model comparison is when you look at different models and their fit
# to gain an understanding of how different parameters contribute and relate to your variables
# and possible choosing between different equivalent causal models. 

# 3
# Information critera are closely tied to the number of observeations and only their relative
# value compared to eachother holds any valuable information. IC values will always be higher
# for models with more observations. 

# 4 
# What happens to effective number of parameters as priors get tighter for PSIS and WAIC?

# Experiment: 

N <- 1000
cats <- rbern(N, prob = 0.5)   # red = 1
crazy <- rnorm(N, mean = cats, sd = 1)  # would imply that red cats are on average 1 crazier than non red cats

d <- data.frame(red = as.numeric(cats), crazy = crazy)
d$crazy <- standardize(d$crazy)


mcat <- quap(alist(crazy ~ dnorm(mu, sigma),
                   mu <- a + bR * red ,
                   a ~ dnorm(0, 0.05),
                   bR ~ dnorm(0.5, 0.01),
                   sigma ~ dexp(1)),
                   data = d)
precis(mcat)

PSIS(mcat)
WAIC(mcat)

# The effective number of parameters becomes smaller. 

# 5
# Informative priors reduce overfitting because it makes the model more sceptical of 
# extreme values which will otherwise have a lot of leverage. Makes it so the model 
# learns from regular features and doesn't get as influenced by irregular features. 

#6
# Overly informative priors causes underfitting because it doesn't allow the model
# to learn anything from the data. If enough data is provided the overly informative prior
# can be overcome. 

# 7H -----------------

# 1 
# Fit curve and compare models

data(Laffer)

d <- Laffer
plot(tax_revenue ~ tax_rate, data = d)
d$tax_rate.s <- standardize(d$tax_rate)
d$tax_revenue.s <- standardize(d$tax_revenue)

m7H.1 <- quap(alist(tax_revenue.s ~ dnorm(mu, sigma),
                    mu <- a + bR * tax_rate.s,
                    a ~ dnorm(0, 1),
                    bR ~ dnorm(0, 1),
                    sigma ~ dexp(1)),
              data = d)
precis(m7H.1)

m7H.1t <- quap(alist(tax_revenue.s ~ dstudent(2, mu, sigma),
                    mu <- a + bR * tax_rate.s,
                    a ~ dnorm(0, 1),
                    bR ~ dnorm(0, 1),
                    sigma ~ dexp(1)),
              data = d)
precis(m7H.1t)

d$tax_rate.2 <- (d$tax_rate.s)^2

m7H.2 <- quap(alist(tax_revenue.s ~ dnorm(mu, sigma),
                     mu <- a + bR * tax_rate.s + bR2 * tax_rate.2,
                     a ~ dnorm(0, 1),
                     bR ~ dnorm(0, 1),
                     bR2 ~ dnorm(0, 1),
                     sigma ~ dexp(1)),
               data = d)
precis(m7H.2)

m7H.2t <- quap(alist(tax_revenue.s ~ dstudent(2, mu, sigma),
                    mu <- a + bR * tax_rate.s + bR2 * tax_rate.2,
                    a ~ dnorm(0, 1),
                    bR ~ dnorm(0, 1),
                    bR2 ~ dnorm(0, 1),
                    sigma ~ dexp(1)),
              data = d)

precis(m7H.2t)

plot(compare(m7H.1, m7H.1t, m7H.2, m7H.2t, func = PSIS))

# 2
# Look at the effect of outliers. 

PSIS_72 <- PSIS(m7H.2, pointwise = T)
WAIC_72 <- WAIC(m7H.2,pointwise=TRUE)
plot( PSIS_72$k , WAIC_72$penalty , xlab="PSIS Pareto k" ,
      ylab="WAIC penalty" , col=col.alpha(rangi2, 0.3) , lwd=2, pch = 16 )

PSIS_71 <- PSIS(m7H.1, pointwise = T)
WAIC_71 <- WAIC(m7H.1,pointwise=TRUE)
points(PSIS_71$k , WAIC_71$penalty , xlab="PSIS Pareto k" ,
      ylab="WAIC penalty" , col=col.alpha("red", 0.3) , lwd=2, pch = 16)

PSIS_71t <- PSIS(m7H.1t, pointwise = T)
WAIC_71t <- WAIC(m7H.1t,pointwise=TRUE)
points( PSIS_71t$k , WAIC_71t$penalty , xlab="PSIS Pareto k" ,
      ylab="WAIC penalty" , col=col.alpha("green", 0.3) , lwd=2 , pch = 16)

PSIS_72t <- PSIS(m7H.2t, pointwise = T)
WAIC_72t <- WAIC(m7H.2t,pointwise=TRUE)
points( PSIS_72t$k , WAIC_72t$penalty , xlab="PSIS Pareto k" ,
        ylab="WAIC penalty" , col=col.alpha("pink", 0.3) , lwd=2 , pch = 16)

# Allowing the tax revenue to be student t distributed ment the model had a much better fit
# and reduced the problems arising from extreme values. 


# Counterfactual plots

tax.rate_seq <- seq(-3, 3, length.out = 30)
tax.rate_seq.2 <- tax.rate_seq^2

par(mfrow = c(1,2))

plot(tax_revenue.s ~ tax_rate.s, data = d, pch = 16, col = col.alpha("cyan4", 0.5),
     xlab = "Standardized tax rate", ylab = "Standardized tax revenue")

# 72
mu <- link(m7H.2, data = data.frame(tax_rate.s = tax.rate_seq, 
                               tax_rate.2 = tax.rate_seq.2), 
           n = 1e4)

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = 0.97)
lines(tax.rate_seq, mu.mean)
shade(mu.PI, tax.rate_seq)

Rev.sim <- sim(m7H.2, data = data.frame(tax_rate.s = tax.rate_seq, 
                                        tax_rate.2 = tax.rate_seq.2), n = 1e4)
Rev.PI <- apply(Rev.sim, 2, PI)

shade(Rev.PI, tax.rate_seq, col = col.alpha("grey", 0.3))

#72t
plot(tax_revenue.s ~ tax_rate.s, data = d, pch = 16, col = col.alpha("cyan4", 0.5),
     xlab = "Standardized tax rate", ylab = "Standardized tax revenue")
mu <- link(m7H.2t, data = data.frame(tax_rate.s = tax.rate_seq, 
                                    tax_rate.2 = tax.rate_seq.2), 
           n = 1e4)

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = 0.97)
lines(tax.rate_seq, mu.mean)
shade(mu.PI, tax.rate_seq)

Rev.sim <- sim(m7H.2t, data = data.frame(tax_rate.s = tax.rate_seq, 
                                        tax_rate.2 = tax.rate_seq.2), n = 1e4)
Rev.PI <- apply(Rev.sim, 2, PI)

shade(Rev.PI, tax.rate_seq, col = col.alpha("grey", 0.3))

# 3
# Islands with different proportion of birds. Lets caluclate entropy and K-L divergence.

p1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
p2 <- c(0.8, 0.1, 0.05, 0.025, 0.025)
p3 <- c(0.05, 0.15, 0.7, 0.05, 0.05)


H1 <- -sum(p1 * log(p1))
H2 <- -sum(p2 * log(p2))
H3 <- -sum(p3 * log(p3))
# H is the amount of uncertainty contained in a prob. distribution 
# high means a lot of uncertainty

KL12 <- sum(p1 * (log(p1) - log(p2)))  # the average difference in log probability between p1 and p2
KL13 <- sum(p1 * (log(p1) - log(p3)))
KL21 <- sum(p2 * (log(p2) - log(p1)))
KL23 <- sum(p2 * (log(p2) - log(p3)))
KL31 <- sum(p3 * (log(p3) - log(p1)))
KL32 <- sum(p3 * (log(p3) - log(p2)))
  # The additional uncertainty introduced by using p2 to predict p1 is KL1
  # Divergence is defined as the additional entropy induced by using q to predict p. 
  # So it’s just the difference between H(p), the actual entropy of events, and 
  # H(p, q): DKL(p, q) = H(p, q) − H(p)
  # So divergence really is measuring how far q is from the target p, in units of entropy

# Which island is best predictor? Check which has the smallest KL divergence

Total.KL <- c(KL12, KL13, KL21, KL23, KL31, KL32)
Total.KL.combined <- c(KL12 + KL13, KL21 + KL23, KL31 + KL32)

# Island number 1 is the best at prediciting the other islands. Because it actually would
# predict the same species fairly often, in comparison to the other two. As they would most
# often predict there very high proportion species. Which is very low in the others



# 5

data(foxes)
d <- foxes
str(d)
d$weight.s <- standardize(d$weight)
d$groupsize.s <- standardize(d$groupsize)
d$avgfood.s <- standardize(d$avgfood)
d$area.s <- standardize(d$area)

m5.1 <- quap(alist(weight ~ dnorm(mu, sigma),
                   mu <- a + bF*avgfood + bG*groupsize + bA*area,
                   a ~ dnorm(0, 1),
                   bF ~ dnorm(0, 1),
                   bG ~ dnorm(0, 1), 
                   bA ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data = d)

m5.2 <- quap(alist(weight ~ dnorm(mu, sigma),
                   mu <- a + bF*avgfood + bG*groupsize,
                   a ~ dnorm(0, 1),
                   bF ~ dnorm(0, 1),
                   bG ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data = d)
m5.3 <- quap(alist(weight ~ dnorm(mu, sigma),
                   mu <- a + bG*groupsize + bA*area,
                   a ~ dnorm(0, 1),
                   bG ~ dnorm(0, 1), 
                   bA ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data = d)

m5.4 <- quap(alist(weight ~ dnorm(mu, sigma),
                   mu <- a + bF*avgfood,
                   a ~ dnorm(0, 1),
                   bF ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data = d)

m5.5 <- quap(alist(weight ~ dnorm(mu, sigma),
                   mu <- a + bA*area,
                   a ~ dnorm(0, 1),
                   bA ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data = d)

compare(m5.1, m5.2, m5.3, m5.4, mx5.5, func = WAIC)
