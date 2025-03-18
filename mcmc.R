# Markov Chain Monte Carlo

library(rethinking)
# " A marvelous example of Fortuna Minerva's cooperation. "

# King Markov's island

num_weeks <- 1e5

positions <- rep(0, num_weeks)
current <- 10

for( i in 1:num_weeks){
  # current position
  positions[i] <- current
  # proposal island
  proposal <- current + sample(c(-1, 1), size = 1)
  # move around island 
  if(proposal > 10) proposal <- 1
  if(proposal < 1) proposal <- 10
  
  # move
  prob_move <- proposal/current
  current <- ifelse(runif(1) < prob_move, proposal, current)
}

plot( table( positions ), xlab = "Islands", ylab = "Weeks Spent" )

# All he needs to know is the size of the current island and size of the proposal island. 
# He can be ignorant of all else. 

# This is a version of a metropolis Algorithm, requires symmetrical proposals
# We want the ability of using asymmentrical proposals and more accuratly finding the
# posterior in fewer samples (efficiency) -> Gibbs sampling.
# Uses smart jumps to suggests proposals for the joint posterior distribution

# Problems: Need conjugate priors, might not want that. Inefficient with large number of
# parameters because they tend to get "stuck". " They dont know where they're going"

# Getting stuck with high parameter values: Concentration of Measure'
# Because: Most of the probability mass of a high parameter distribution is far away from the mode

D <- 10000
Tr <- 1e2
Y <- rmvnorm(Tr, rep(0, D), diag(D))
rad_dist <- function(Y) sqrt(Y^2)

Rd <- sapply( 1:T , function(i) rad_dist( Y[i,] ) )
dens( Rd )

# Need algorithms that focus on entire distribtution at once. -> Hamilitonian Monte Carlo


# Hamiltonian Monte Carlo ------------------------------
# Needs 2(4) things: (1) Negative-log probability of the data's given position. 
# (2) The gradiants around the position, given by derivatives. 
# (3) Step size, (4) Leap frog length

# Simulation two parameter example
# Suppose the data are 100 x values and 100 y values sampled from Normal(0, 1) 

# U is negative log probability
U <- function(q, a = 0, b = 1, k = 0, d = 1){
  muy <- q[1]
  mux <- q[2]
  U <- sum(dnorm(y, muy, 1, log = TRUE)) + sum(dnorm(x, mux, 1, log = TRUE)) +
            dnorm(muy, a, b, log = TRUE) + dnorm(mux, k, d, log = TRUE)
  return(-U)
}

# Gradient : Partial derivates for each parameter

U_gradient <- function(q, a = 0, b = 1, k = 0, d = 1){
  muy <- q[1]
  mux <- q[2]
  G1 <- sum(y - muy) + (a - muy/b^2)  # dU/dmuy
  G2 <- sum(x - mux) + (k - mux/b^2)  # dU/dmux
  return(c(-G1, -G2))
}

# test. data

set.seed(7)
y <- rnorm(50)
x <- rnorm(50)

x <- as.numeric(scale(x))
y <- as.numeric(scale(y))

# Fancy figure!

library(shape)
Q <- list()
Q$q <- c(-0.1, 0.2)
pr <- 0.3
plot(NULL, xlab = "mux", ylab = "muy", xlim = c(-pr, pr), ylim = c(-pr, pr))

step <-  0.03
L <- 11
n_samples <- 4

path_col <- col.alpha("black", 0.4)
points(Q$q[1], Q$q[2], pch = 4, col = "black")

for(i in 1:n_samples){
  Q <- HMC2(U, U_gradient, step, L, Q$q)
  if(n_samples < 10){
    for ( j in 1:L ) {
      K0 <- sum(Q$ptraj[j,]^2)/2 # kinetic energy
      lines( Q$traj[j:(j+1),1] , Q$traj[j:(j+1),2] , col=path_col , lwd=1+2*K0 )
    }
    points( Q$traj[1:L+1,] , pch=16 , col="white" , cex=0.35 )
    Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] ,
            arr.length=0.35 , arr.adj = 0.7 )
    text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=0.8 , pos=4 , offset=0.4 )
  }
  points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 16 , 1 ) ,
          col=ifelse( abs(Q$dH)>0.1 , "red" , "black" ) )
}


# Ulam -----------------------

data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )


m8.3 <- quap(
              alist(
                log_gdp_std ~ dnorm( mu , sigma ) ,
                mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
                a[cid] ~ dnorm( 1 , 0.1 ) ,
                b[cid] ~ dnorm( 0 , 0.3 ) ,
                sigma ~ dexp( 1 )
              ) , data=dd )
precis( m8.3 , depth=2 )


dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)
str(dat_slim)
# list allow elements to be of different lengths.


# Hamiltonian monte carlo posterior
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat_slim , chains=1 )

precis(m9.1, depth = 2)

m9.1 <- ulam(alist(
                log_gdp_std ~ dnorm( mu , sigma ) ,
                mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
                a[cid] ~ dnorm( 1 , 0.1 ) ,
                b[cid] ~ dnorm( 0 , 0.3 ) ,
                sigma ~ dexp( 1 )
              ), data=dat_slim , chains=4 , cores=4 )


precis(m9.1, depth = 2)
show(m9.1)
pairs(m9.1)
traceplot(m9.1)
trankplot(m9.1, n_cols = 2)


# Dangerous flat priors

y <- c(-1,1)
set.seed(11)
m9.2 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- alpha ,
    alpha ~ dnorm( 0 , 1000 ) ,
    sigma ~ dexp( 0.0001 )
  ) , data=list(y=y) , chains=3 )

precis(m9.2 )
traceplot(m9.2)
trankplot(m9.2)


# Lets make them sliiightly more informative

set.seed(11)
m9.3 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- alpha ,
    alpha ~ dnorm( 1 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=list(y=y) , chains=3 )

precis( m9.3 )
traceplot(m9.3)

# Non identifiable parameters

set.seed(41)
y <- rnorm( 100 , mean=0 , sd=1 )

set.seed(384)
m9.4 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- a1 + a2 ,
    a1 ~ dnorm( 0 , 1000 ),
    a2 ~ dnorm( 0 , 1000 ),
    sigma ~ dexp( 1 )
  ) , data=list(y=y) , chains=3 )

precis( m9.4 )
traceplot(m9.4)

# More informative priors

m9.5 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- a1 + a2 ,
    a1 ~ dnorm( 0 , 10 ),
    a2 ~ dnorm( 0 , 10 ),
    sigma ~ dexp( 1 )
  ) , data=list(y=y) , chains=3 )
precis( m9.5 )

show(m9.4)
show(m9.5)

# Much better estimate of their sum and also much faster
# Often, a model that is very slow to sample is under-identified.
# Just a little prior information telling the model 
# “none of these parameters can be 30 million” often helps. 

# Practice ---------------------------------------------------------------------


# 9M1

data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )

m9.1 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dexp( 1 )
), data=dat_slim , chains=4 , cores=4 )

precis(m9.1, depth = 2)
traceplot(m9.1)
post <- extract.samples(m9.1)
dens(post$a)
dens(post$b)

m9M.1 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 0 , 10 ) ,
  b[cid] ~ dnorm( 0 , 10 ) ,
  sigma ~ dexp( 1 )
), data=dat_slim , chains=4 , cores=4 )

precis(m9H.1, depth = 2)
traceplot(m9H.1)
post <- extract.samples(m9H.1)
dens(post$a)
dens(post$b)

y <- rcauchy(1e4, 0, 5)
mu <- sapply(1:length(y), function(i) sum(y[1:i]/i))
plot(mu, type = "l")


# 9M2

m9M.2 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dexp( 1 )
), data=dat_slim , chains=4 , cores=4 )

m9M.3 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dexp( 5 )
), data=dat_slim , chains=4 , cores=4 )

m9M.4 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dexp( 10 )
), data=dat_slim , chains=4 , cores=4 )

precis(m9M.2, depth = 2)
precis(m9M.3, depth = 2)
precis(m9M.4, depth = 2)

m9M.5 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dcauchy(0, 2)
), data=dat_slim , chains=4 , cores=4 )

m9M.6 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dcauchy(0, 0.5)
), data=dat_slim , chains=4 , cores=4 )

m9M.7 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dcauchy(0, 0.1 )
), data=dat_slim , chains=4 , cores=4 )

precis(m9M.5, depth = 2)
precis(m9M.6, depth = 2)
precis(m9M.7, depth = 2)

# 9M3

m9M.8 <- ulam(alist(
  log_gdp_std ~ dnorm( mu , sigma ) ,
  mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
  a[cid] ~ dnorm( 1 , 0.1 ) ,
  b[cid] ~ dnorm( 0 , 0.3 ) ,
  sigma ~ dcauchy(0, 2)
), data=dat_slim , chains=4 , cores=4, iter = 400, warmup = 50)

traceplot(m9M.8)
precis(m9M.8, depth = 2)
