
# God spiked the integers 
# The moving pieces in generalized linear models are not as transperent as the ones
# in linear models. It's hard to interpret the meaning of a specific parameter by itself.
# Looking at just one gear of the machine is not enough. 
library(rethinking)

# Prosocial chimpanzee ---------------------------------------------------

data("chimpanzees")
d <- chimpanzees
?chimpanzees

d$treatment <- 1 + d$prosoc_left + 2*d$condition

# 1) One food, No partner, 2) two food, no partner, 3) One food, partner, 4) Two food, partner

xtabs( ~ treatment + prosoc_left + condition , d )

# Choosing priors

m11.1 <- quap(alist(pulled_left ~ dbinom(1, p),
                    logit(p) <- a,
                    a ~ dnorm(0, 10)),
              data = d)


set.seed(1999)

prior <- extract.prior(m11.1, n = 1e4)
p <- inv_logit(prior$a)
dens(p, adj = 0.1)   # Because most of the probability mass will be outside -4 and 4 

m11.1.1 <- quap(alist(pulled_left ~ dbinom(1, p),
                    logit(p) <- a,
                    a ~ dnorm(0, 1.5)),
              data = d)

prior <- extract.prior(m11.1.1, n = 1e4)
p <- inv_logit(prior$a)
dens(p, adj = 0.1)   # Now is prety flat, much better. 

# Next prior test

m11.2 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 10 )
  ) , data=d )

set.seed(1999)
prior <- extract.prior( m11.2 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
dens( abs( p[,1] - p[,2] ) , adj=0.1 )   # once again, kinda shit

# More reasonable priors

m11.3 <- quap( 
               alist(
                 pulled_left ~ dbinom( 1 , p ) ,
                 logit(p) <- a + b[treatment] ,
                 a ~ dnorm( 0 , 1.5 ),
                 b[treatment] ~ dnorm( 0 , 0.5 )
               ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.3 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
mean( abs( p[,1] - p[,2] ) )
dens( abs( p[,1] - p[,2] ), xlim = c(0, 1) , adj = 0.1)
dens(p[,1])


# With MCMC

# trim data

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = d$treatment
)

m11.4 <- ulam( 
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list, chains = 4, log_lik = TRUE)

precis( m11.4 , depth=2 )    # Hard to interpret 

post <- extract.samples(m11.4)
p_left <- inv_logit( post$a )
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )   # Preferance for pulling levers per monkey

labs <- c("R/N","L/N","R/P","L/P")
plot( precis( m11.4 , depth=2 , pars="b" ) , labels=labs )


# Differences between food on right partner/no partner and food on left partner no partner
diffs <- list(
               db13 = post$b[,1] - post$b[,3],    # Positive would mean they are less likely to pull left when there is a partner and food on right
               db24 = post$b[,2] - post$b[,4] )   # Here is should be negative if prosocial is real
plot( precis(diffs) )    # Not lots of evidence :() 


# Posterior Predictive check ------------------------


# Plotting observed

pl <- by( d$pulled_left , list( d$actor , d$treatment ) , mean )  # Each column is a 
#treatment. And the cells contain proportions of pulls that were of the left lever.
pl[4,]  # First actors proportion of left pulls per treatment

plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("actor ",j) , xpd=TRUE )

for ( j in (1:7)[-2] ) {
  lines( (j-1)*4+c(1,3) , pl[j,c(1,3)] , lwd=2 , col=rangi2 )
  lines( (j-1)*4+c(2,4) , pl[j,c(2,4)] , lwd=2 , col=rangi2 )
}

points( 1:28 , t(pl) , pch=16 , col="white" , cex=1.7 )
points( 1:28 , t(pl) , pch=c(1,1,16,16) , col=rangi2 , lwd=2 )

yoff <- 0.01
text( 1 , pl[1,1]-yoff , "R/N" , pos=1 , cex=0.8 )
text( 2 , pl[1,2]+yoff , "L/N" , pos=3 , cex=0.8 )
text( 3 , pl[1,3]-yoff , "R/P" , pos=1 , cex=0.8 )
text( 4 , pl[1,4]+yoff , "L/P" , pos=3 , cex=0.8 )
mtext( "observed proportions\n" )

# Drawing samples from posterior

dat <- list( actor=rep(1:7,each=4) , treatment=rep(1:4,times=7) )
p_post <- link( m11.4 , data=dat )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )


d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2


dat_list2 <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  side = d$side,
  cond = d$cond )
m11.5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + bs[side] + bc[cond] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    bs[side] ~ dnorm( 0 , 0.5 ),
    bc[cond] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list2 , chains=4 , log_lik=TRUE )

compare( m11.5 , m11.4 , func=PSIS )  # Not much difference

# Relative shark, absolute penguin

post <- extract.samples(m11.4)
mean( exp(post$b[,4]-post$b[,2]) )   # On average, the switch multiples the odds 
                                     # of pulling the left lever by 0.92, an 8% reduction in odds.

# Relative odds can be dangerous, because they don't inform us how big an effect actually is
# But conditionally very important. A penguin does not care about the absolute odds of dying of  a shark
# attack, but rather the relative odds of dying of a shark attack conditional on being in the ocean
# (and a delicious penguin).

# Aggregated binomial ----------------------------

data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2
d_aggregated <- aggregate(
  d$pulled_left ,
  list( treatment=d$treatment , actor=d$actor ,
        side=d$side , cond=d$cond ) ,
  sum )
colnames(d_aggregated)[5] <- "left_pulls"


dat <- with( d_aggregated , list(
  left_pulls = left_pulls,
  treatment = treatment,
  actor = actor,
  side = side,
  cond = cond ) )
m11.6 <- ulam(
  alist(
    left_pulls ~ dbinom( 18 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ) ,
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat , chains=4 , log_lik=TRUE )

precis(m11.6, depth = 2)
precis(m11.4, depth = 2)

# Very similar, makes sense.

compare(m11.4, m11.6, func = PSIS)

# Shows very different PSIS values
# Because more ways to see data in bernoulli from

# deviance of aggregated 6-in-9 
-2*dbinom(6,9,0.2,log=TRUE)
# deviance of dis-aggregated
-2*sum(dbern(c(1,1,1,1,1,1,0,0,0),0.2,log=TRUE))


# What’s the bottom line? If you want to calculate WAIC or PSIS, you should use a logistic
# regression data format, not an aggregated format. Otherwise you are implicitly assuming
# that only large chunks of the data are separable.  



# Aggregated binomial with unbalanced data, collage admission -------------

data("UCBadmit")
d <- UCBadmit
d

# Lets look for gender discrimination

dat_list <- list(
  admit = d$admit,
  applications = d$applications,
  gid = ifelse( d$applicant.gender=="male" , 1 , 2 )
)

m11.7 <- ulam(
  alist(admit ~ dbinom(applications, p),
        logit(p) <- a[gid],
        a[gid] ~ dnorm(0, 1.5)),
  data = dat_list, chains = 4) 
precis(m11.7, depth = 2)

# Lets compare differences

post <- extract.samples(m11.7)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis( list( diff_a=diff_a , diff_p=diff_p ) )

# Posterior prediction check

postcheck( m11.7 )
# draw lines connecting points from same dept
for ( i in 1:6 ) {
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
  text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}


dat_list$dept_id <- rep(1:6,each=2)
m11.8 <- ulam(
  alist(
    admit ~ dbinom( applications , p ) ,
    logit(p) <- a[gid] + delta[dept_id] ,
    a[gid] ~ dnorm( 0 , 1.5 ) ,
    delta[dept_id] ~ dnorm( 0 , 1.5 )
  ) , data=dat_list , chains=4 , iter=4000 )
precis( m11.8 , depth=2 )

# Comparing differences
post <- extract.samples(m11.8)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis( list( diff_a=diff_a , diff_p=diff_p ) )

# When sorted by department males seem to be (very slightly disadvantaged)

pg <- with( dat_list , sapply( 1:6 , function(k)
                               applications[dept_id==k]/sum(applications[dept_id==k]) ) )
rownames(pg) <- c("male","female")
colnames(pg) <- unique(d$dept)
round( pg , 2 )
# More males apply to A and B, which have high admittance rates

postcheck(m11.8)
pairs(m11.8)



# Poisson regression ---------------------------------------------------------
# No upper bound, very low probability many trials: Binomial -> Poisson

# Expected of binomial = N*p, variance = N*p(1-p), if p is very small this also becomes N*p

y <- rbinom(1e5, 1000, 1/1000)
c(mean(y), variance(y))

# Lets us model events where the the number of events is unknown or uncountably large. 


# The log link ensures that λi is always positive, which is required of the expected value of
# a count outcome. But as mentioned in the previous chapter, it also implies an exponential
# relationship between predictors and the expected value.

data("Kline")
d <- Kline
head(Kline)
d$P <- scale(log(d$population))
d$contact_id <- ifelse(d$contact=="high", 2, 1)


# Investigate priors on log scale

# A normal(0, 10) prior: (would become log normal)
curve( dlnorm( x , 0 , 10 ) , from=0 , to=100 , n=200 )
range(d$total_tools)   # Prior does not match data very well at all

a <- rnorm(1e4,0,10)
lambda <- exp(a)
mean( lambda )

# Which would mean a mean of 330468803110 tools per island 
# Playing with prior
curve( dlnorm( x , 3 , 0.5) , from=0 , to=100 , n=200 )

a <- rnorm(1e4,3,1)
lambda <- exp(a)
mean( lambda )

mean(d$total_tools)  # Dont actually look at this, we dont want our prior to fit out sample
# just live in the plausible outcome space

# For Poisson models, flat priors make no sense and can wreck Prague.

# Population prior

N <- 100
a <- rnorm( N , 3 , 0.5 )
b <- rnorm( N , 0 , 10 )
plot( NULL , xlim=c(-2,2) , ylim=c(0,100), xlab = "Log Population", ylab = "Total Tools" )
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=grau() ) # Horrible

N <- 100
a <- rnorm( N , 3 , 0.5 )
b <- rnorm( N , 0 , 0.2 )
plot( NULL , xlim=c(-2,2) , ylim=c(0,100), xlab = "Log Population", ylab = "Total Tools" )
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=grau() )

# Dampens prior to not believe explosive relationships 
# Strong relationships, possible but less likely


# Log population has a natural zero, lets not remove that. Its bad for thinking
range(d$population)

x_seq <- seq( from=log(100) , to=log(200000) , length.out=100 )
lambda <- sapply( x_seq , function(x) exp( a + b*x ) )
plot( NULL , xlim=range(x_seq) , ylim=c(0,500) , xlab="log population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( x_seq , lambda[i,] , col=grau() , lwd=1.5 )


# Natural scale

plot( NULL , xlim=range(exp(x_seq)) , ylim=c(0,500) , xlab="population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( exp(x_seq) , lambda[i,] , col=grau() , lwd=1.5 )


# Model

dat <- list( 
             T = d$total_tools ,
             P = d$P ,
             cid = d$contact_id )
# intercept only
m11.9 <- ulam(
  alist(
    T ~ dpois( lambda ),
    log(lambda) <- a,
    a ~ dnorm(3,0.5)
  ), data=dat , chains=4 , log_lik=TRUE )

precis(m11.9, depth = 2)


# interaction model
m11.10 <- ulam(
  alist(T ~ dpois( lambda ),
        log(lambda) <- a[cid] + b[cid]*P,
        a[cid] ~ dnorm( 3 , 0.5 ),
        b[cid] ~ dnorm( 0 , 0.2 )
  ), data=dat , chains=4 , log_lik=TRUE )


compare( m11.9 , m11.10 , func=PSIS )

k <- PSIS( m11.10 , pointwise=TRUE )$k
plot( dat$P , dat$T , xlab="log population (std)" , ylab="total tools" ,
      col=rangi2 , pch=ifelse( dat$cid==1 , 1 , 16 ) , lwd=2 ,
      ylim=c(0,75) , cex=1+normalize(k) )
# set up the horizontal axis values to compute predictions at
ns <- 100
P_seq <- seq( from=-1.4 , to=3 , length.out=ns )
# predictions for cid=1 (low contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )

lines( P_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE )
# predictions for cid=2 (high contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( P_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE )


plot( d$population , d$total_tools , xlab="population" , ylab="total tools" ,
      col=rangi2 , pch=ifelse( dat$cid==1 , 1 , 16 ) , lwd=2 ,
      ylim=c(0,75) , cex=1+normalize(k) )
ns <- 100
P_seq <- seq( from=-5 , to=3 , length.out=ns )
# 1.53 is sd of log(population)
# 9 is mean of log(population)
pop_seq <- exp( P_seq*1.53 + 9 )
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )

lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )


# A better model, scientific. Not as small world.

dat2 <- list( T=d$total_tools, P=d$population, cid=d$contact_id )
m11.11 <- ulam(
  alist(
    T ~ dpois( lambda ),
    lambda <- exp(a[cid])*P^b[cid]/g,
    a[cid] ~ dnorm(1,1),
    b[cid] ~ dexp(1),
    g ~ dexp(1)
  ), data=dat2 , chains=4 , log_lik=TRUE )


# Negative binomial/ Gamma poisson ---------------------------------------------
# Is poisson in disguise, or rather a mixture of poisson distributions. 

# We need to add a logarithm of exposure to the model!
# Lets make some data to analyze
num_days <- 30
y <- rpois( num_days , 1.5 )

num_weeks <- 4
y_new <- rpois(num_weeks, 0.5*7)

y_all <- c(y, y_new)
exposure <- c(rep(1, 30), rep(7, 4))
monastery <- c(rep(0, 30), rep(1, 4))

d <- data.frame(y = y_all, days = exposure, monastery = monastery)

# compute the offset
d$log_days <- log(d$days)

# fit the model
m11.12 <- quap(
  alist(
    y ~ dpois( lambda ),
    log(lambda) <- log_days + a + b*monastery,
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm( 0 , 1 )
  ), data=d )



post <- extract.samples( m11.12 )
lambda_old <- exp( post$a )
lambda_new <- exp( post$a + post$b )
precis( data.frame( lambda_old , lambda_new ) )


# Multinomial and categorical models -------------------------------------------

# Uses multinomial distribution. Need to build K-1 linear models (k being the number of
# possible outcomes/events we are trying to predict)

# Predictors matched to outcome (Careers choice example)

# Lets make data: Simulate career choices

N <- 500
income <- c(1, 2, 5)    # Expected income of each career
score <- 0.5*income     # Score for each career

# Convert scores to probabilties:
p <- softmax(score[1], score[2], score[3])    
# Formula for score 1: exp(score[1])/(sum(exp(score)))

# Simulate choice

career <- rep(NA, N)  # Vector to put choices in
set.seed(34302)
for(i in 1:N) career[i] <- sample(c(1:3), size = 1, prob = p)

plot(table(career))


# Stan model

# Same parameter linked to different outcomes 
code_m11.13 <- "
data{
int N; // number of individuals
int K; // number of possible careers
array[N] int career; // outcome
vector[K] career_income;
}
parameters{
vector[K-1] a; // intercepts
real<lower=0> b; // association of income with choice
}
model{
vector[K] p;
vector[K] s;
a ~ normal( 0 , 1 );
b ~ normal( 0 , 0.5 );
s[1] = a[1] + b*career_income[1];
s[2] = a[2] + b*career_income[2];
s[3] = 0; // pivot
p = softmax( s );
career ~ categorical( p );
}
"

dat_list <- list( N=N , K=3 , career=career , career_income=income )
m11.13 <- stan( model_code=code_m11.13 , data=dat_list , chains=4 )
precis( m11.13 , 2 )  # Cannot be interpreted by themselves easily, have to be converted


# Counterfactual 

post <- extract.samples(m11.13)

# logit scores

s1 <- with( post , a[,1] + b*income[1] )
s2_orig <- with( post , a[,2] + b*income[2] )
s2_new <- with( post , a[,2] + b*income[2]*2 )   # Counterfactual, double income


# compute probabilities for original and counterfactual
p_orig <- sapply( 1:length(post$b) , function(i)
  softmax( c(s1[i],s2_orig[i],0) ) )
p_new <- sapply( 1:length(post$b) , function(i)
  softmax( c(s1[i],s2_new[i],0) ) )
# summarize
p_diff <- p_new[2,] - p_orig[2,]
precis( p_diff )

apply(p_orig, 1, mean)
p

# Predictors matched to observations

N <- 500
family_income <- runif(N)
# Unique coefficient for each career (event)
b <- c(-2, 0, 2)     
career <- rep(NA, N)

for(i in 1:N){
  score <- 0.5*c(1:3) + b*family_income[i]
  p <- softmax(score[1], score[2], score[3])
  career[i] <- sample(1:3, size = 1, prob = p)
}


code_m11.14 <- "
    data{
      int N; // number of observations
      int K; // number of outcome values
      array[N] int career; // outcome
      array[N] real family_income;
    }
      parameters{
      vector[K-1] a; // intercepts
      vector[K-1] b; // coefficients on family income
    }
      model{
      vector[K] p;
      vector[K] s;
      a ~ normal(0,1.5);
      b ~ normal(0,1);
      for ( i in 1:N ) {
      for ( j in 1:(K-1) ) s[j] = a[j] + b[j]*family_income[i];
      s[K] = 0; // the pivot
      p = softmax( s );
      career[i] ~ categorical( p );
      }
    }
    "

dat_list <- list( N=N , K=3 , career=career , family_income=family_income )
m11.14 <- stan( model_code=code_m11.14 , data=dat_list , chains=4 )
precis( m11.14 , 2 )



# Multinomial in disguise as a serise of poissons ------------------------------


data(UCBadmit)
d <- UCBadmit

# binomial model of overall admission probability
m_binom <- quap(
  alist(
    admit ~ dbinom(applications,p),
    logit(p) <- a,
    a ~ dnorm( 0 , 1.5 )
  ), data=d )
precis(m_binom)

# Poisson model of overall admission rate and rejection rate
# 'reject' is a reserved word in Stan, cannot use as variable name
dat <- list( admit=d$admit , rej=d$reject )
m_pois <- ulam(
  alist(
    admit ~ dpois(lambda1),
    rej ~ dpois(lambda2),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    c(a1,a2) ~ dnorm(0,1.5)
  ), data=dat , chains=3 , cores=3 )
precis(m_pois)

inv_logit(coef(m_binom))  # Probability of being accepted

k <- coef(m_pois)
a1 <- k['a1']
a2 <- k['a2']
exp(a1)/(exp(a1)+exp(a2))    # Implied probability

# Almost exactly the same answer! :0


# Censoring and Survival -------------------------------------------------------

# Modeling rates in survival models. Rates -> displacement, real positive values
# Maximum entropy distributions are: Exponential (all we know is average displacement)
# or Gamma (fixed mean value and fixed mean magnitude (log))

# Censoring: When an event doesn't occur within our observational window
# You need to also account for the cats who haven’t yet been adopted. The cats who
# haven’t been adopted yet, but eventually will be adopted, clearly have longer waiting 
# times than the cats who have already been adopted.

# The chance of a cat being adopted after x amount of days is given by the cumulative probability
# distribution. So 1 - CDF is the probability of a cat not having been adopted by X amount of days
# = complementary cumulative probability distribution


data("AustinCats")
d <- AustinCats

d$adopt <- ifelse(d$out_event == "Adoption", 1L, 0L)
dat <- list(
  days_to_event = as.numeric(d$days_to_event),
  color_id = ifelse(d$color == "Black", 1L, 2L),
  adopted = d$adopt
)

m11.5 <- ulam(
  alist(
    days_to_event | adopted == 1 ~ exponential(lambda),
    days_to_event | adopted == 0 ~ custom(exponential_lccdf(!Y|lambda)),
    lambda <- 1.0/mu,
    log(mu) <- a[color_id],
    a[color_id] ~ dnorm(0, 1)
  ), data = dat, chains = 4, cores = 4)


traceplot(m11.5)
trankplot(m11.5)

precis(m11.5, 2)   # The black cats get adopted faster :(

post <- extract.samples(m11.5)
post$D <- exp(post$a)
precis(post, 2)


# Practice ---------------------------------------------------------------------


# 11H1

data("chimpanzees")
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = d$treatment
)

m11.4 <- ulam( 
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list, chains = 4, log_lik = TRUE)

precis( m11.4 , depth=2 )    # Hard to interpret 

post <- extract.samples(m11.4)
p_left <- inv_logit( post$a )
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )   # Preferance for pulling levers per monkey

labs <- c("R/N","L/N","R/P","L/P")
plot( precis( m11.4 , depth=2 , pars="b" ) , labels=labs )


# Differences between food on right partner/no partner and food on left partner no partner
diffs <- list(
  db13 = post$b[,1] - post$b[,3],    # Positive would mean they are less likely to pull left when there is a partner and food on right
  db24 = post$b[,2] - post$b[,4] )   # Here is should be negative if prosocial is real
plot( precis(diffs) )

mH1 <- quap( 
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list)

precis( mH1 , depth=2 )    # Hard to interpret 

post <- extract.samples(mH1)
p_left <- inv_logit( post$a )
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )   # Preferance for pulling levers per monkey

labs <- c("R/N","L/N","R/P","L/P")
plot( precis( mH1 , depth=2 , pars="b" ) , labels=labs )


# Differences between food on right partner/no partner and food on left partner no partner
diffs <- list(
  db13 = post$b[,1] - post$b[,3],    # Positive would mean they are less likely to pull left when there is a partner and food on right
  db24 = post$b[,2] - post$b[,4] )   # Here is should be negative if prosocial is real
plot( precis(diffs) )

# Pretty much the same. Because no interaction?

dat_list <- list(
  pulled_left = d$pulled_left,
  treatment = d$treatment
)

# H2
m11.3 <- ulam( 
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list, chains = 4, log_lik = TRUE )


compare(m11.3, m11.4, func = WAIC)
# Significantly better

# H3  Salmon pirates
library(MASS)
data("eagles")
?eagles
d <- eagles

dat_list <- list(
  success = d$y,
  attempts = d$n,
  P_size = ifelse( d$P == "L" , 2 , 1 ),   # Large/adult/ = 2
  P_age = ifelse(d$A == "A", 2, 1),
  V_size = ifelse(d$V == "L", 2, 1)
)

m11H.3 <- ulam(
  alist(success ~ dbinom(attempts, p),
        logit(p) <- a + bs[P_size] + ba[P_age] + bv[V_size] ,
        a ~ dnorm(0, 1.5),
        bs[P_size] ~ dnorm(0, 0.5), 
        ba[P_age] ~ dnorm(0, 0.5),
        bv[V_size] ~ dnorm(0, 0.5)), 
  data = dat_list, chains = 4, log_lik = TRUE) 

m11H.4 <- quap(
  alist(success ~ dbinom(attempts, p),
        logit(p) <- a + bs[P_size] + ba[P_age] + bv[V_size] ,
        a ~ dnorm(0, 1.5),
        bs[P_size] ~ dnorm(0, 0.5), 
        ba[P_age] ~ dnorm(0, 0.5),
        bv[V_size] ~ dnorm(0, 0.5)), 
  data = dat_list) 

precis(m11H.3, depth = 2)
plot(precis(m11H.4, depth = 2))



p_size <- inv_logit( post$bs[, 2] - post$bs[,1] )
p_age <- inv_logit( post$ba[, 2] - post$ba[,1] )
v_size <- inv_logit( post$bv[, 2] - post$bv[,1] )

dens(inv_logit(post$bs[,2]))   # Chance of success for large birds
dens(inv_logit(post$bs[,1]))
par(mfrow = c(1, 3))
dens(p_size, xlim = c(0, 1)) 
mtext("Effect of being large")
dens(p_age, xlim = c(0, 1)) 
mtext("Effect of Adult")
dens(v_size, xlim = c(0, 1))
mtext("Effect of attacking large victim")

# Attacking a large victim usually leads to failiure 


m11H.5 <- quap(
  alist(success ~ dbinom(attempts, p),
        logit(p) <- a + bs[P_size] * ba[P_age] + bv[V_size],
        a ~ dnorm(0, 1.5),
        bs[P_size] ~ dnorm(0, 0.5), 
        ba[P_age] ~ dnorm(0, 0.5),
        bv[V_size] ~ dnorm(0, 0.5)), 
  data = dat_list) 

precis(m11H.5, depth = 3)
plot(precis(m11H.4, depth = 2))



p_size <- inv_logit( post$bs[, 2] - post$bs[,1] )
p_age <- inv_logit( post$ba[, 2] - post$ba[,1] )
v_size <- inv_logit( post$bv[, 2] - post$bv[,1] )

dens(inv_logit(post$bs[,2]))   # Chance of success for large birds
dens(inv_logit(post$bs[,1]))
par(mfrow = c(1, 3))
dens(p_size, xlim = c(0, 1)) 
mtext("Effect of being large")
dens(p_age, xlim = c(0, 1)) 
mtext("Effect of Adult")
dens(v_size, xlim = c(0, 1))
mtext("Effect of attacking large victim")

# Attacking a large victim usually leads to failiure 




