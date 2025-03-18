# Models with memory

# "Cafes do differ, but they are also alike"
# "Anterograde amnesia is bad for learning about the world"
# Models that learn (at the same time) information about a cluster and information
# about the population of clusters are more effective. 
# So our models tracks the parameter for each cluster as well as the atleast two paramteres
# to describe the population: an average and a standard deviation. 
# These population parameters are used as the prior for all the cafes. 

library(rethinking)

# Tadpoles ---------------------------------------------------------------------

data("reedfrogs")
d <- reedfrogs

# Varying intercept model: Allowing each tank to have a separate intercept but
# also the variation of this intercept between tanks.

d$tank <- c(1:nrow(d))

dat <- list(S = d$surv,
            N = d$density,
            tank = d$tank)

# Unique intercept for each tank: 
m13.1 <- ulam(alist(
  S ~ dbinom(N, p),
  logit(p) <- a[tank],
  a[tank] ~ dnorm(0, 1.5)
), data = dat, chains = 4, cores = 4, log_lik = TRUE)

precis(m13.1, depth = 2)


# Not lets make it multilevel with hyperparameters


m13.2 <- ulam(alist(
  S ~ dbinom(N, p),
  logit(p) <- a[tank],
  a[tank] ~ dnorm(a_bar, sigma),
  a_bar ~ dnorm(0, 1.5),
  sigma ~ dexp(1)
), data = dat, chains = 4, cores = 4, log_lik = TRUE)

precis(m13.2, depth = 2)

compare(m13.1, m13.2)  # The model with more parameters actually have effectively less parameters :0
trankplot(m13.2)

# Plotting

post <- extract.samples(m13.2)

# Compute median intercept for each tank and transform to probability
d$propsurv.est <- logistic(apply(post$a, 2, mean))


plot( d$propsurv , ylim=c(0,1) , pch=16 , xaxt="n" ,
      xlab="tank" , ylab="proportion survival" , col=rangi2 )
axis( 1 , at=c(1,16,32,48) , labels=c(1,16,32,48) )

points(d$propsurv.est, col = rgb(0.1, 0.1, 0.1, 0.4), size = 3, pch = 16)

abline(h = mean(inv_logit(post$a_bar)), lty = 2)

abline( v=16.5 , lwd=0.5 )
abline( v=32.5 , lwd=0.5 )
text( 8 , 0 , "small tanks" )
text( 16+8 , 0 , "medium tanks" )
text( 32+8 , 0 , "large tanks" )

# See shrinkage (regression to mean) and a result of regularization (not bad)



# show first 100 populations in the posterior
plot( NULL , xlim=c(-3,4) , ylim=c(0,0.35) ,
      xlab="log-odds survive" , ylab="Density" )

for ( i in 1:100 ){
  curve( dnorm(x,post$a_bar[i],post$sigma[i]) , add=TRUE ,
         col=col.alpha("black",0.2) )}
# sample 8000 imaginary tanks from the posterior distribution
sim_tanks <- rnorm( 8000 , post$a_bar , post$sigma )
# transform to probability and visualize
dens( inv_logit(sim_tanks) , lwd=2 , adj=0.1 , xlab = "Probability survive")

# Over- and underfitting trade-off ---------------------------------------------

# Complete pooling vs No pooling vs Partial pooling= underfitting, overfitting, just right.
# Complete pooling assumes variation among ponds is zero.
# No pooling assumes the variation among ponds is infinite.
# Partial pooling assumes there is some variation, and estimates this. 

# Lets simulate some data to prove this to ourselves!

a_bar <- 1.5
sigma <- 1.5
nponds <- 60
Ni <- as.integer(rep(c(5, 10, 25, 35), each = 15)) # Density (initial number of individuals in pond)
a_pond <- rnorm(nponds, mean = a_bar, sd = sigma)  # Simulate log odds probability of survival

dsim <- data.frame(pond=1:nponds , Ni=Ni , true_a=a_pond )

# Simulate survival:

dsim$Si <- rbinom(n = nponds, size = Ni, prob = logistic(dsim$true_a))

# Comparison time:
# No pooling estimate:

dsim$p_nopool <- dsim$Si/dsim$Ni

# Partial pooling:

dat <- list(Si = dsim$Si, Ni = dsim$Ni, pond = dsim$pond)

m13.3 <- ulam(alist(
  Si ~ dbinom(Ni, p),
  logit(p) <-  a_pond[pond],
  a_pond[pond] ~ dnorm(a_bar, sigma),
  a_bar ~ dnorm(0, 1.5),
  sigma ~ dexp(1)
), data = dat, chains = 4, cores = 4)

traceplot(m13.3)
precis(m13.3, depth = 2)

post <- extract.samples(m13.3)
dsim$p_partpool <- apply(inv_logit(post$a_pond), 2, mean)



# Complete pooling:

m13.4 <- quap(alist(
  Si ~ dbinom(Ni, p),
  logit(p) <- a,
  a ~ dnorm(0, 1.5)
), data = dsim)

precis(m13.4)
post <- extract.samples(m13.4)

dsim$p_compool <- mean(inv_logit(post$a))

# True probabilities

dsim$p_true <- inv_logit( dsim$true_a )

# Differences:

nopool_error <- abs(dsim$p_nopool - dsim$p_true)
partpool_error <- abs(dsim$p_partpool - dsim$p_true)
compool_error <- abs(dsim$p_compool - dsim$p_true)

plot(1:60, nopool_error, xlab = "Pond", ylab = "Absolute Error", col = rangi2, pch = 16)
points(1:60, partpool_error, pch = 16, col = "pink")
points(1:60, compool_error, pch = 16, col = "palegreen")
abline(h = mean(nopool_error), col = rangi2, lty = 2)
abline(h = mean(partpool_error), col = "pink", lty = 2)
abline(h = mean(compool_error), col = "palegreen", lty = 2)

nopool_avg <- aggregate(nopool_error,list(dsim$Ni),mean)
partpool_avg <- aggregate(partpool_error,list(dsim$Ni),mean)

# at the higher density ponds (i.e more data) there is less difference between the methods
# both estimates are better in bigger ponds, and become approxiamtly even. 
# Ponds with extreme estimates will have less error in the nopool, because they are pulled
# towards the mean in the part pool. 



# Several clusters: Multilevel chimpanzee --------------------------------------
# Partial pooling over several clusters.

data(chimpanzees)
d <- chimpanzees

d$treatment <- 1 + d$prosoc_left + 2*d$condition

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  block_id = d$block,
  treatment = as.integer(d$treatment) )
set.seed(13)


m13.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + g[block_id] + b[treatment] ,
    b[treatment] ~ dnorm( 0 , 0.5 ),
    ## adaptive priors
    a[actor] ~ dnorm( a_bar , sigma_a ),
    g[block_id] ~ dnorm( 0 , sigma_g ),
    ## hyper-priors
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )



traceplot(m13.4)
trankplot(m13.4)

precis(m13.4, depth = 2)

# Low number of samples and high r-hat are proof of inefficient sampling.

# Only one multilevel: (ignores block)

set.seed(14)
m13.5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    b[treatment] ~ dnorm( 0 , 0.5 ),
    a[actor] ~ dnorm( a_bar , sigma_a ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )
precis(m13.5, depth = 2)
compare(m13.4, m13.5)

# More partial pooling: (on treatment effects)
set.seed(15)

m13.6 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + g[block_id] + b[treatment] ,
    b[treatment] ~ dnorm( 0 , sigma_b ),
    a[actor] ~ dnorm( a_bar , sigma_a ),
    g[block_id] ~ dnorm( 0 , sigma_g ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    sigma_b ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )

traceplot(m13.6)
trankplot(m13.6)
precis(m13.6, depth = 2)

coeftab( m13.4 , m13.6 )   # because there is a lot of data for each treatment this wont change the result much


# Divergent transitions --------------------------------------------------------
# When certain parts of the posterior are very steep the physics simulation behind
# MCMC has a hard time exploring it. Causing divergent transitions. Happens because
# the simulation is not actually continuos but discrete. 

# Example: The devils funnel

m13.7 <- ulam(alist(
  v ~ dnorm(0, 3),
  x ~ dnorm(0, exp(v))
), data = list(N = 1), chains = 4)
precis(m13.7)
traceplot(m13.7)

# Solution: Non centered parameterization:

m13.7nc <- ulam(
                 alist(
                   v ~ normal(0,3),
                   z ~ normal(0,1),
                   gq> real[1]:x <<- z*exp(v)
                 ), data=list(N=1) , chains=4 )
precis( m13.7nc )

# Also fix by changing stans acceptance rate, higher values usually leads to smaller
# step size, which can help solve divergent transition. (standard = 0.95)

m13.4b <- ulam(m13.4, chains = 4, cores = 4, control = list(adapt_delta = 0.99))
divergent(m13.4)
divergent(m13.4b)

# Helped slightly

precis(m13.4b)  # But still inefficient, number of sampels is way lower than 2000



set.seed(13)
m13.4nc <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a_bar + z[actor]*sigma_a + # actor intercepts
      x[block_id]*sigma_g + # block intercepts
      b[treatment] ,
    b[treatment] ~ dnorm( 0 , 0.5 ),
    z[actor] ~ dnorm( 0 , 1 ),
    x[block_id] ~ dnorm( 0 , 1 ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    gq> vector[actor]:a <<- a_bar + z*sigma_a,
    gq> vector[block_id]:g <<- x*sigma_g
  ) , data=dat_list , chains=4 , cores=4 )

divergent(m13.4nc)  # no divergent transitions

precis_c <- precis( m13.4 , depth=2 )
precis_nc <- precis( m13.4nc , depth=2 )
pars <- c( paste("a[",1:7,"]",sep="") , paste("g[",1:6,"]",sep="") ,
           paste("b[",1:4,"]",sep="") , "a_bar" , "sigma_a" , "sigma_g" )
neff_table <- cbind( precis_c[pars,"ess_bulk"] , precis_nc[pars,"ess_bulk"] )
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centered)" , ylab="n_eff (non-centered)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# Points above the diagonal are parameters that had better sampling for the 
# non centered model. 

# Non centered parameters will usually perform better for clusters
# if there is low variation. And also when if you have a large number of units 
# inside a cluster, but not much data for each unit.

# Its not always a good idea to have non centered priors, and sometimes
# only a good idea for a few parameters. 


# Multilevel posterior predictions ---------------------------------------------

# We no longer expect the posterior predictive distribution to match the sample,
# because the goal of partial pooling is to shrink estimates towards the grand mean.
# Giving up in sample prediction to increase out of sample prediction.


chimp <- 2   # Only looking at chimp 2
d_pred <- list(
  actor = rep(chimp,4),
  treatment = 1:4,
  block_id = rep(1,4)
)
p <- link( m13.4 , data=d_pred )
p_mu <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI )


# manually: 
post <- extract.samples(m13.4)

p_link <- function( treatment , actor=1 , block_id=1 ) {
  logodds <- with( post ,
                   a[,actor] + g[,block_id] + b[,treatment] )
  return( inv_logit(logodds) )
}


p_raw <- sapply( 1:4 , function(i) p_link( i , actor=2 , block_id=1 ) )
p_mu <- apply( p_raw , 2 , mean )
p_ci <- apply( p_raw , 2 , PI )


# Constructing new actors:

# Average pull of average monket for each treatment:
p_link_abar <- function( treatment ) {
  logodds <- with( post , a_bar + b[,treatment] )
  return( inv_logit(logodds) )
}


post <- extract.samples(m13.4)
p_raw <- sapply( 1:4 , function(i) p_link_abar( i ) )
p_mu <- apply( p_raw , 2 , mean )
p_ci <- apply( p_raw , 2 , PI )
plot( NULL , xlab="treatment" , ylab="proportion pulled left" ,
      ylim=c(0,1) , xaxt="n" , xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("R/N","L/N","R/P","L/P") )
lines( 1:4 , p_mu )
shade( p_ci , 1:4)


# Simulate new chimpanzees:

a_sim <- with( post , rnorm( length(post$a_bar) , a_bar , sigma_a ) )  # New chimp with average left pulls
p_link_asim <- function( treatment ) {
  logodds <- with( post , a_sim + b[,treatment] )
  return( inv_logit(logodds) )
}
p_raw_asim <- sapply( 1:4 , function(i) p_link_asim( i ) )


plot( NULL , xlab="treatment" , ylab="proportion pulled left" ,
      ylim=c(0,1) , xaxt="n" , xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("R/N","L/N","R/P","L/P") )
for ( i in 1:100 ) lines( 1:4 , p_raw_asim[i,] , col=grau(0.25) , lwd=2 )
