# Adventure in covariance
library(rethinking)

# It is likely that more than one affect varies in a data set.
# For example different cafes will have different waiting times, but
# different times of day will also have different waiting times. 
# But different cafes will also differ in how much they differ.

# Any batch of parameters with index values should be pooled: Varying effects strategy.
# But cafes usually co-vary in their intercept and slopes.

# Varying slopes ---------------------------------------------------------------

# Simulate data:

a <- 3.5        # Average morning wait time
b <- (-1)       # Average difference afternoon wait time
sigma_a <- 1    # sd in average wait times
sigma_b <- 0.5  # sd in difference in wait times
rho <- -(0.7)   # correlation between wait time and difference in wait time (intercept and slopes)

Mu <- c(a, b)

cob_ab <- sigma_a * sigma_b * rho 
Sigma <- matrix(c(sigma_a^2, cob_ab, cob_ab, sigma_b^2), ncol = 2)

# Other way of defining it covariance matrix:

sigmas <- c(sigma_a, sigma_b)               # Standard deviation
Rho <- matrix(c(1, rho, rho, 1), ncol = 2)  # Correlation matrix

# now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 20

library(MASS)
set.seed(5)
vary_effects <- mvrnorm(N_cafes, Mu, Sigma)
colnames(vary_effects) <- c("Wait time", "Change in Afternoon wait time")
rownames(vary_effects) <- c(paste0("Cafe ", c(1:nrow(vary_effects))))

a_cafe <- vary_effects[, 1]
b_cafe <- vary_effects[, 2]

plot( a_cafe , b_cafe , col=rangi2 ,
      xlab="intercepts (a_cafe)" , ylab="slopes (b_cafe)",
      pch = 16 )
# overlay population distribution
library(ellipse)
for ( l in c(0.1,0.3,0.5,0.8,0.99) ){
  lines(ellipse(Sigma,centre=Mu,level=l),col=col.alpha("black",0.2))}
