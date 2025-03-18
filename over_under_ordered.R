
# Mosnters and mixtures

library(rethinking)


# Over dispersed counts ----------------------------------------------------------

# Processes are often variable mixtures which results in thicker tails -> Student -t 
# Over dispersion ->
# When counts arrise from a variety of processes 
# Dealing with over dispersion (Ideally you would find variables causing the dispersion and adding them to model)



# Beta binomial = assumes each binomial count observation has it's own probability of success

# The probabilities have a beta distribution: a probability distribution for probabilities. 

# Beta distribution has two parameters, average probability and a shape parameter, Theta.
# Theta = 2 means all proababilities between 0 and 1 are equally likely, above 2 makes them 
# more concentrated around the average and when its below 2 extreme probabilites closer to 
# one and zero are more likely than the mean probability. 

pbar <- 0.5
theta <- c(1:6)

par(mfrow = c(3, 2))
for(i in 1:length(theta)) {
  curve(dbeta2(x, pbar, theta[i]), from = 0, to = 1, 
        xlab = "Probability", ylab = "Density")
  abline(v = 0.5, lty = 2, col = "grey")
  mtext(paste0("theta = ", theta[i]))
}

par(mfrow = c(1, 1))
curve(dbeta2(x, 0.6, 4), from = 0, to = 1, 
      xlab = "Probability", ylab = "Density")
abline(v = 0.5, lty = 2, col = "grey")


# So we bind our predictor to phat, thereby changing the central tendancy of out distribution



data(UCBadmit)
d <- UCBadmit
d$gid <- ifelse( d$applicant.gender=="male" , 1L , 2L )

dat <- list( A=d$admit , N=d$applications , gid=d$gid )
m12.1 <- ulam(
  alist(
    A ~ dbetabinom( N , pbar , theta ),
    logit(pbar) <- a[gid],
    a[gid] ~ dnorm( 0 , 1.5 ),
    transpars> theta <<- phi + 2.0,   # Because we cant allow the theta to be less than 2
    phi ~ dexp(1)),                   # Transpars will make stan return it in the samples
  data=dat , chains=4 )


post <- extract.samples(m12.1)
post$da <- post$a[, 1] - post$a[, 2]
precis(post, depth = 3)

gid <- 2
curve(dbeta2(x, mean(logistic(post$a[,gid])), mean(post$theta)), from = 0, to = 1, 
      ylab = "Density", xlab = "Probability admit", lwd = 1.5, ylim = c(0, 3))

# Draw 50 posterior mean distributions
for(i in 1:50){
  p <- logistic(post$a[i, gid])
  theta <- post$theta[i]
  curve(dbeta2(x, p, theta), add = TRUE, col = rgb(0, 0, 0, 0.2))
}
mtext("Distribution of female admission rates")
# The model is no longer surprised when female admission is high in some departments

postcheck(m12.1)

# Negative binomial (Gamma Poisson) = Assumes each poisson count has its own rate.

# Estimates this rate from a gamma distribution, very useful, since it allows for changing the
# variance (esentially)
# Gamma distributions have two parameters, one for mean(lambda) and one for dispersion (phi).
# The variance of a gamma poisson is: lambda + lambda^2/phi
# Phi most be positive and with higher phi the distribution will be narrower 

data(Kline)
d <- Kline

d$P <- standardize(log(d$population))
d$contact_id <- ifelse(d$contact == "high", 2L, 1L)  # L indicates its an integer

dat2 <- list(
  T = d$total_tools,
  P = d$P,
  cid = d$contact_id
)

m12.2 <- ulam(
  alist(
    T ~ dgampois(lambda, phi),
    lambda <- exp(a[cid]) * P^b[cid] / g,
    a[cid] ~ dnorm(1, 1),
    b[cid] ~ dexp(1),
    g ~ dexp(1),
    phi ~ dexp(1)
  ), data = d, chains = 4, log_lik = TRUE)


# Zero inflated outcomes -------------------------------------------------------

# Mixtures of multiple processes, zeros can often arise is several different ways
# so not just one process. 

# Example: manuscript transcribing monks partaking in wine.
# Zeros can be generated from two processes (1) The monks spent the day drinking, (2)
# They worked but failed to finish any manuscript. 
# p = probability monks spend the day drinking, lambda = mean number of manuscripts completed
# when monks work.
# Drinking monks produce zero manuscrips, working monks produce a poisson number of manuscripts


# Generate data 

prob_drink <- 0.2  # 20 % of days
rate_work <- 1  # Average of 1 manuscript per day

# Sample one year of production
N <- 365

# simulate days monks drink
set.seed(365)
drink <- rbinom( N , 1 , prob_drink )
# simulate manuscripts completed
y <- (1-drink)*rpois( N , rate_work )

simplehist(y, xlab = "Manuscripts Completed", lwd = 4)

zeros_drink <- sum(drink)
zeros_work <- sum(y == 0 & drink == 0)
zeros_total <- sum(y == 0)
lines(c(0,0), c(zeros_work, zeros_total), lwd = 4, col = rangi2) # c(x, x), c(y, y)


m12.3 <- ulam( 
               alist(
                 y ~ dzipois( p , lambda ),
                 logit(p) <- ap,
                 log(lambda) <- al,
                 ap ~ dnorm( -1.5 , 1 ),    # Mostly between -2.5 and 0.5, meaning mostly below 0.5 percent in the logit scale
                 al ~ dnorm( 1 , 0.5 )
               ) , data=list(y=y) , chains=4 )
precis( m12.3 )

logistic(-1.29)
exp(0.01)

post <- extract.samples(m12.3)

mean(inv_logit(post$ap))  # Probability drink
mean(exp(post$al))        # Rate of transcription


# Ordered categorical outcomes -------------------------------------------------

# Treating categorical ordered outcomes as continuous measure is a bad idea. 
# a shift from 1 to 2 might be different from shifting between 5 and 6.

# Cumulative link function: The cumulative probability is the probability of a specific value
# and any value smaller than that. So the cumulative probability of 3 is the sum of the 
# probabilities of 3 and 2 and 1. 

# Linking a linear model to cumulative probability can guarantee a ordering of the
# outcomes 

# Moral reprehensibility

data("Trolley")
d <- Trolley

simplehist(d$response, xlim = c(1, 7), xlab = " Response")

# Our goal is to redescribe this on the log cumulative odds scale : Cumulative analog of
# logit link

# Proportion of each response variable:

pr_k <- table(d$response)/nrow(d)

cum_pr_k <- cumsum(pr_k)  # converts to cumulative
plot(1:7, cum_pr_k, type = "b", xlab = "Response", ylab = "Cumulative Proportion", ylim = c(0,1))

# Log cumulative odd:
logit <- function(x) log(x/(1-x)) # convenience function
round( lco <- logit( cum_pr_k ) , 2 )


# Incorportating no predictor variables:

m12.4 <- ulam(
               alist(
                 R ~ dordlogit( 0 , cutpoints ),
                 cutpoints ~ dnorm( 0 , 1.5 )
               ) , data=list( R=d$response ), chains=4 , cores=4 )

precis(m12.4, depth = 2)

round( inv_logit(coef(m12.4)) , 2 )

round( pk <- dordlogit( 1:7 , 0 , coef(m12.4) ) , 2 )
sum( pk*(1:7) )   # Average outcome variable


round( pk <- dordlogit( 1:7 , 0 , coef(m12.4)-0.5 ) , 2 )
sum(pk*(1:7))

# Adding predictors

dat <- list(
  R = d$response,
  A = d$action,
  I = d$intention,
  C = d$contact )
m12.5 <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints ),
    phi <- bA*A + bC*C + BI*I ,
    BI <- bI + bIA*A + bIC*C ,
    c(bA,bI,bC,bIA,bIC) ~ dnorm( 0 , 0.5 ),
    cutpoints ~ dnorm( 0 , 1.5 )
  ) , data=dat , chains=4 , cores=4 )
precis( m12.5 )

plot( precis(m12.5) , xlim=c(-1.4,0) )

plot( NULL , type="n" , xlab="intention" , ylab="probability" ,
      xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )


kA <- 0 # value for action 
kC <- 0 # value for contact
kI <- 0:1 # values of intention to calculate over
pdat <- data.frame(A=kA,C=kC,I=kI)
phi <- link( m12.5 , data=pdat )$phi

post <- extract.samples( m12.5 )
for ( s in 1:50 ) {
  pk <- pordlogit( 1:6 , phi[s,] , post$cutpoints[s,] )
  for ( i in 1:6 ) lines( kI , pk[,i] , col=grau(0.1) )
}


kA <- 0 # value for action
kC <- 1 # value for contact
kI <- 0:1 # values of intention to calculate over
pdat <- data.frame(A=kA,C=kC,I=kI)
s <- sim( m12.5 , data=pdat )
simplehist( s , xlab="response" )


# Ordered predictor variables --------------------------------------------------


edu_levels <- c( 6 , 1 , 8 , 4 , 7 , 2 , 5 , 3 )
d$edu_new <- edu_levels[ d$edu ]


library(gtools)
set.seed(1805)
delta <- rdirichlet( 10 , alpha=rep(2,7) )    # higher (2) means you expect them to be more equal
h <- 3
plot( NULL , xlim=c(1,7) , ylim=c(0,0.4) , xlab="index" , ylab="probability" )
for ( i in 1:nrow(delta) ) lines( 1:7 , delta[i,] , type="b" ,
                                  pch=ifelse(i==h,16,1) , lwd=ifelse(i==h,4,1.5) ,
                                  col=ifelse(i==h,"black",col.alpha("black",0.7)) )


dat <- list(
             R = d$response ,
             action = d$action,
             intention = d$intention,
             contact = d$contact,
             E = as.integer( d$edu_new ), # edu_new as an index
             alpha = rep( 2 , 7 ) ) # delta prior

m12.6 <- ulam(
  alist(
    R ~ ordered_logistic( phi , kappa ),
    phi <- bE*sum( delta_j[1:E] ) + bA*action + bI*intention + bC*contact,
    kappa ~ normal( 0 , 1.5 ),
    c(bA,bI,bC,bE) ~ normal( 0 , 1 ),
    vector[8]: delta_j <<- append_row( 0 , delta ),
    simplex[7]: delta ~ dirichlet( alpha )
  ), data=dat , chains=4 , cores=4 )


precis( m12.6 , depth=2 , omit="kappa" )


delta_labels <- c("Elem","MidSch","SHS","HSG","SCol","Bach","Mast","Grad") 
pairs( m12.6 , pars="delta" , labels=delta_labels )




dat$edu_norm <- normalize( d$edu_new ) 
m12.7 <- ulam(
  alist(R ~ ordered_logistic( mu , cutpoints ),
        mu <- bE*edu_norm + bA*action + bI*intention + bC*contact,
        c(bA,bI,bC,bE) ~ normal( 0 , 1 ),
        cutpoints ~ normal( 0 , 1.5 )
  ), data=dat , chains=4 , cores=4 )
precis( m12.7 )



# Practice ---------------------------------------------------------------------


# 12E:
#1 An ordered categorical variable is clearly ordered. Like Different levels of 
# temperature: Low, Medium, High. Where the effect of the variable would increase/decrease
# with the ordered. An unordered variable could be favorite ice cream flavors: Chocolate, Vanilla, Strawberry

#2 A log cumulative odds function. Uses the cumulative odds instead of just the odds.

#3 It will tend to underestimate the rate of events because a count distribution with 
# extra zeros added to it will have a lower mean. 

#4 Number of species of flowers in a 1x1 m square would probably be overdispersed?:
# arises do tue variation in underlying rates of different units. Number of ice creams
# sold at different shops aggregated (different shops will have different rates)
# A common cause of under disperion is autocorrelation.

# 12M:
#1 Compute log cumulative odds:

categories <- as.integer(c(1, 2, 3, 4))
ratings <- c(12, 36, 7, 41)
ratios <- ratings/sum(ratings)

cum_sums <- ratios[1]

for(i in 1:3){
  cum_sums[i+1] <- cum_sums[i] + ratios[i+1]
}

log_cum_odds <- log(cum_sums/(1-cum_sums))

#2

plot(1:4, cum_sums, type = "b", col = "black", lwd = 1.5, xlab = "Employee Rating",
     ylab = "Cumulative proportion", xaxt = "n", ylim = c(0, 1))
axis(1, c(1:4))
lines(1:4, cum_sums, type = "b", col = "black", lwd = 1.5)
for(i in 1:4) lines(c(i-0.01, i-0.01), c(0, cum_sums[i]), lwd = 3, col = "grey")
for(i in 1:4) lines(c(i+0.01, i+0.01), c(cum_sums[i], cum_sums[i] - ratios[i]), lwd = 3, col = "lightblue")



#3 Zero inlfated binomial 


# Generate data 
# On days you work there is a 0.8 percent chance you will succeed at finishing a task
# but on 20 percent of days you decide to just fuck it and drink instead. Then you will
# get nothing done. 

prob_drink <- 0.2  # 20 % of days
probability_success <- 0.8  # Average of 1 manuscript per day

# Sample one year of work
N <- 365

# simulate days I drink
drink <- rbinom( N , 1, prob_drink )
# simulate days with completed work
y <- (1-drink)*rbinom(N, 1, probability_success)

simplehist(y, xlab = "Work Completed", lwd = 4)

zeros_drink <- sum(drink)
zeros_work <- sum(y == 0 & drink == 0)
zeros_total <- sum(y == 0)
lines(c(0,0), c(zeros_work, zeros_total), lwd = 4, col = rangi2) # c(x, x), c(y, y)


m12M <- ulam( 
  alist(
    y ~ dzibinom(p_zero, N, p),
    logit(p_zero) <- ap,
    logit(p) <- al,
    ap ~ dnorm( -1.5 , 1 ),    # Mostly between -2.5 and 0.5, meaning mostly below 0.5 percent in the logit scale
    al ~ dnorm( 1.5 , 1 )
  ) , data=list(y=y, N = 365) , chains=4 )

traceplot(m12M)
trankplot(m12M)
precis( m12M )

logistic(-3.9)
logistic(-6.33)

post <- extract.samples(m12M)

mean(inv_logit(post$ap))        # Probability drink
mean(inv_logit(post$al))        # Probability of finishing work


#12H

#1
data(Hurricanes)
d <- Hurricanes
?Hurricanes

d$fem_std <- with(d, (femininity - mean(femininity))/sd(femininity))

data_list <- list(deaths = d$deaths, fem = d$fem_std)

m12H <- ulam( 
  alist(
    deaths ~ dpois(lambda),
    log(lambda) <- a + bf*fem,
    a ~ dnorm(0, 10),
    bf ~ dnorm( 0 , 1 )
  ), data=list(deaths = d$deaths, fem = d$fem_std), chains=4 )

traceplot(m12H)
trankplot(m12H)
precis(m12H)

m12H.2 <- ulam( 
  alist(
    deaths ~ dpois(lambda),
    log(lambda) <- a,
    a ~ dnorm(0, 10)
  ), data=list(deaths = d$deaths), chains=4)
traceplot(m12H.2)
precis(m12H.2)

plot(d$fem_std, d$deaths, pch = 16, col = rangi2, xlab = "Femininity", ylab = "Number of deaths")


fem_seq <- seq(-2, 2, length.out = 30)
lambda <- link(m12H, data = list(fem = fem_seq))

mu <- apply(lambda, 2, mean)
mu.PI <- apply(lambda, 2, PI)

lines(fem_seq, mu)
shade(mu.PI, fem_seq)

# compute sampling distribution
deaths_sim <- sim(m12H,data=list(fem = fem_seq))
deaths_sim.PI <- apply(deaths_sim,2,PI)
# superimpose sampling interval as dashed lines
lines( fem_seq , deaths_sim.PI[1,] , lty=2 )
lines( fem_seq , deaths_sim.PI[2,] , lty=2 )

# Data overdispersed, not a very good predictor. 


#2 Use gamma poisson (negative binomial)

m12H2 <- ulam(
  alist(
    deaths ~ dgampois(lambda, phi),
    log(lambda) <- a + bf*fem,
    a ~ dnorm(0, 10),
    bf ~ dnorm(0, 1),
    phi ~ dexp(1)
  ), data=list(deaths = d$deaths, fem = d$fem_std), chains = 4)

traceplot(m12H2)

precis(m12H2)  # relationship no longer as strong
plot(coeftab(m12H,m12H2))   # Much higher standard deviation, less sure about connection


# plot raw data
plot( d$fem_std , d$deaths , pch=16 ,
      col=rangi2 , xlab="femininity" , ylab="deaths" )
# compute model-based trend
pred_dat <- list( fem=seq(from=-2,to=1.5,length.out=30) )
lambda <- link(m12H2,data=pred_dat)
lambda.mu <- apply(lambda,2,mean)
lambda.PI <- apply(lambda,2,PI)
# superimpose trend
lines( pred_dat$fem , lambda.mu )
shade( lambda.PI , pred_dat$fem )
# compute sampling distribution
deaths_sim <- sim(m12H2,data=pred_dat)
deaths_sim.PI <- apply(deaths_sim,2,PI)
# superimpose sampling interval as dashed lines
lines( pred_dat$fem , deaths_sim.PI[1,] , lty=2 )
lines( pred_dat$fem , deaths_sim.PI[2,] , lty=2 )


#3 






