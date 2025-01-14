
# First install Stan c++ toolchain: https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows

# Install Rstan interface: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# we recommend running this in a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

library(cmdstanr)
cmdstanr::install_cmdstan()


install.packages(c("coda","mvtnorm","devtools","dagitty"))
library(devtools)
devtools::install_github("rmcelreath/rethinking")

#### Globe tossing example: Trying to deduce the proportion of water on a globe 
#    by throwing it in the air and noting where your finger is when you catch it
#    on water = W, on land = L.

# Each toss is essentially a binomial trial: 
dbinom(6, 9, prob = 0.5)

# Grid approximation

# 1) Define grid

p_grid <- seq(from = 0, to = 1, length.out = 20)

# 2) Compute value of prior at each parameter value in grid

prior <- rep(1, times = 20)

  # Alt priors
  prior <- ifelse( p_grid < 0.5 , 0 , 1 )
  prior <- exp( -5*abs( p_grid - 0.5 ) )

# 3) Compute likelihood at each parameter value in grid of 6 water in 9 trials

likelihood <- dbinom(3, 4, prob = p_grid)

# 4) Compute unstandardized posterior

unstan_posterior <- likelihood*prior

# 5) Standerdize posterior

posterior <- unstan_posterior/(sum(unstan_posterior))

# Plot it
plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

# Quadratic estimation

library(rethinking)
globe.qa <- quap(
  alist(
    W ~ dbinom( W+L ,p) , # binomial likelihood
    p ~ dunif(0,1) # uniform prior
  ) ,
  data=list(W=6,L=3) )
# display summary of quadratic approximation
precis( globe.qa )

# Analytical calculation:

# analytical calculation
W <- 6
L <- 3
curve( dbeta( x , W+1 , L+1 ) , from=0 , to=1 )
# quadratic approximation
curve( dnorm( x , 0.67 , 0.16 ) , lty=2 , add=TRUE )


# Markov chain Monte Carlo

n_samples <- 1000
p <- rep( NA , n_samples )
p[1] <- 0.5
W <- 6
L <- 3
for ( i in 2:n_samples ) {
  p_new <- rnorm( 1 , p[i-1] , 0.1 )
  if ( p_new < 0 ) p_new <- abs( p_new )
  if ( p_new > 1 ) p_new <- 2 - p_new
  q0 <- dbinom( W , W+L , p[i-1] )
  q1 <- dbinom( W , W+L , p_new )
  p[i] <- ifelse( runif(1) < q1/q0 , p_new , p[i-1] )
}

dens( p , xlim=c(0,1) )
curve( dbeta( x , W+1 , L+1 ) , lty=2 , add=TRUE )
