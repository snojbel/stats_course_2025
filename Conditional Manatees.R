# Coniditional Manatees

# "Every inference is conditional on something"

library(rethinking)
library(viridis)



# Africa is special ------------------------------------------------------------

data(rugged)
d <- rugged

# Wealth is often exponentially associated with parameters as having wealth generates
# more wealth. Meaning log(wealth) would be linearly associated.

d$log_gdp <- log(d$rgdppc_2000)

dd <- d[complete.cases(d$log_gdp), ]  # remove rows with NA in gdp



dd$log_gdp_stan <- dd$log_gdp/mean(dd$log_gdp) # rescaled as a proportion of the 
# international average. 1 means average, 0.8 means 80% of the average, 
# and 1.1 means 10% more than average.
dd$rugged_stan <- dd$rugged/max(dd$rugged)   # Goes from totally flat = 0 to maximum ruggedness = 1

m8.1 <- quap(alist(log_gdp_stan ~ dnorm(mu, sigma),
                   mu <- a + b*(rugged_stan - 0.215),  # 0.215 is mean of rugged_stan, makes it easier to choose priors
                   a ~ dnorm(1, 1),   # Alpha is log gdp when ruggedness is at its sample mean, meaning it should be close to 1 with the rescaling
                   b ~ dnorm(0, 1),
                   sigma ~ exp(1)), data = dd)  # 0 would indicate no positive or negative relationship

set.seed(7)
# first lets inspect our priors
prior <- extract.prior(m8.1)

plot(NULL, xlim = c(0, 1), ylim = c(0.5, 1.5), xlab = "Ruggedness", ylab = "Log GDP")
abline(h = min(dd$log_gdp_stan), lty = 2)
abline(h = max(dd$log_gdp_stan), lty = 2)

# draw lines from priod
rugged_seq <- seq(-0.1, 1.1, length.out = 30)
mu <- link(m8.1, post = prior, data = data.frame(rugged_stan = rugged_seq))
mu <- mu[[2]]    # extract mu
for(i in 1:50) lines(rugged_seq, mu[i,], col = col.alpha("black", 0.2))     
points(mean(dd$rugged_stan), 1, pch = 16, col = col.alpha("magenta", 0.4))

# Many of these lines are impossible considering the max and min gdp and, so we need tighter 
# prior on alpha
# First of all, lines should pass close to the point of mean ruggedness and mean log(GDP)
# So we need a tighter value of a
# An implausibly strong association between log gdp and ruggedness would
# be, for example, a line that goes from minimum ruggedness and extreme GDP on one end to
# maximum ruggedness and the opposite extreme of GDP on the other end. This would suggest
# all the variation in gdp is due to rugged terrain.
# The slope of such a line must be about 1.3 − 0.7 = 0.6, the difference
# between the maximum and minimum observed proportional log GDP. But very many lines
# in the prior have much more extreme slopes than this. Under the β ∼ Normal(0, 1) prior,
# more than half of all slopes will have absolute value greater than 0.6.

sum( abs(prior$b) > 0.6 ) / length(prior$b)  # proportion with ridiculous associations

m8.1 <- quap(alist(log_gdp_stan ~ dnorm(mu, sigma),
                   mu <- a + b*(rugged_stan - 0.215),  # 0.215 is mean of rugged_stan, makes it easier to choose priors
                   a ~ dnorm(1, 0.1),   # Alpha is log gdp when ruggedness is at its sample mean, meaning it should be close to 1 with the rescaling
                   b ~ dnorm(0, 0.3),
                   sigma ~ exp(1)), data = dd)  # 0 would indicate no positive or negative relationship

# first lets inspect our priors
prior <- extract.prior(m8.1)

plot(NULL, xlim = c(0, 1), ylim = c(0.5, 1.5), xlab = "Ruggedness", ylab = "Log GDP")
abline(h = min(dd$log_gdp_stan), lty = 2)
abline(h = max(dd$log_gdp_stan), lty = 2)

# draw lines from priod
rugged_seq <- seq(-0.1, 1.1, length.out = 30)
mu <- link(m8.1, post = prior, data = data.frame(rugged_stan = rugged_seq))
mu <- mu[[2]]    # extract mu
for(i in 1:50) lines(rugged_seq, mu[i,], col = col.alpha("black", 0.2))  
points(mean(dd$rugged_stan), 1, pch = 16, col = col.alpha("magenta", 0.4))

# Still has some impossible slopes, but much more plausible than before. Majority seem reasonable

precis(m8.1)
# No association , lets split it into africa and not africa

# First lets create index, 1 = africa, 2 = not africa
dd$cid <- ifelse(dd$cont_africa == 1, 1, 2)

m8.2 <- quap(alist(log_gdp_stan ~ dnorm(mu, sigma),
                   mu <- a[cid] + b*(rugged_stan - 0.215),
                   a[cid] ~ dnorm(1, 0.1),
                   b ~ dnorm(0, 0.3),
                   sigma ~ dexp(1))
             , data = dd)

precis(m8.2, depth = 2)

compare(m8.1, m8.2, func = WAIC)

# parameter 1, africa is reliable lower than parameter 2:

post <- extract.samples(m8.2)
post_mu <- post[[3]]
diff_al_a2 <- post_mu[, 1] - post_mu[, 2]
PI(diff_al_a2)


# Plotting posterior distribution

rugged_seq <- seq(-0.1, 1.1, length.out = 30)

mu_not_africa <- link(m8.2, 
                      data = data.frame(cid = 2, rugged_stan = rugged_seq))

mu_africa <- link(m8.2, 
                      data = data.frame(cid = 1, rugged_stan = rugged_seq))

mu.africa.mu <- apply(mu_africa, 2, mean)
mu.africa.ci <- apply(mu_africa, 2, PI, prob = 0.97)

mu.not.africa.mu <- apply(mu_not_africa, 2, mean)
mu.not.africa.ci <- apply(mu_not_africa, 2, PI, prob = 0.97)

plot(log_gdp_stan ~ rugged_stan, data = dd[dd$cid==1, ], 
     xlim = c(0, 1), ylim = c(0.5, 1.5), 
     col = col.alpha("magenta", 0.4),
     xlab = "Ruggedness", ylab = "Log GDP", pch = 16)
points(log_gdp_stan ~ rugged_stan, data = dd[dd$cid==2, ], col = col.alpha("cyan4", 0.4), pch = 16)
abline(h = min(dd$log_gdp_stan), lty = 2)
abline(h = max(dd$log_gdp_stan), lty = 2)
lines(rugged_seq, mu.not.africa.mu, col = "cyan4" )
lines(rugged_seq, mu.africa.mu, col = "magenta")
shade(mu.not.africa.ci, rugged_seq)
shade(mu.africa.ci, rugged_seq)

# This change in model will only display the overall lower wealth in african GDP but does
# not allow the slope of the association to change. So lets make slope conditional aswell

m8.3 <- quap(alist(log_gdp_stan ~ dnorm(mu, sigma),
                   mu <- a[cid] + b[cid]*(rugged_stan - 0.215),
                   a[cid] ~ dnorm(1, 0.1),
                   b[cid] ~ dnorm(0, 0.3),
                   sigma ~ dexp(1))
             , data = dd)

precis(m8.3, depth = 2)

compare( m8.1 , m8.2 , m8.3 , func=PSIS )
plot( PSIS( m8.3 , pointwise=TRUE )$k)  # shows some influential countries


# Plotting

par(mfrow = c(1, 2))

d.A1 <- dd[dd$cid == 1, ]

plot(d.A1$rugged_stan, d.A1$log_gdp_stan, pch = 16, col = rangi2,
     xlab = "Ruggedness ( Standardized)", ylab = "Log GDP (as proportion of mean)",
     xlim = c(0, 1))

mu <- link(m8.3, data = data.frame(cid = 1, rugged_stan = rugged_seq))
mu.mean <- apply(mu, 2, mean)
mu.CI.1 <- apply(mu, 2, PI, prob = 0.69)
mu.CI.2 <- apply(mu, 2, PI, prob = 0.89)
mu.CI.3 <- apply(mu, 2, PI, prob = 0.97)


shade(mu.CI.3, rugged_seq, col = col.alpha(mako(10)[9], 0.5))
shade(mu.CI.2, rugged_seq, col = col.alpha(mako(10)[8], 0.5))
shade(mu.CI.1, rugged_seq, col = col.alpha(mako(10)[7], 0.5))
lines(rugged_seq, mu.mean, col = col.alpha("black", 0.5), lwd = 2)

mtext("African Nations")



d.A2 <- dd[dd$cid == 2, ]

plot(d.A2$rugged_stan, d.A2$log_gdp_stan, pch = 16, col = rangi2,
     xlab = "Ruggedness ( Standardized)", ylab = "Log GDP (as proportion of mean)",
     xlim = c(0, 1))

mu <- link(m8.3, data = data.frame(cid = 2, rugged_stan = rugged_seq))
mu.mean <- apply(mu, 2, mean)
mu.CI.1 <- apply(mu, 2, PI, prob = 0.69)
mu.CI.2 <- apply(mu, 2, PI, prob = 0.89)
mu.CI.3 <- apply(mu, 2, PI, prob = 0.97)


shade(mu.CI.3, rugged_seq, col = col.alpha(mako(10)[9], 0.5))
shade(mu.CI.2, rugged_seq, col = col.alpha(mako(10)[8], 0.5))
shade(mu.CI.1, rugged_seq, col = col.alpha(mako(10)[7], 0.5))
lines(rugged_seq, mu.mean, col = col.alpha("black", 0.5), lwd = 2)

mtext("Non-African Nations")


compare(m8.1, m8.2, m8.3, func = PSIS) # Shows that model m8.3 does do better, but the error for the PSIS
                                 # values is as large as its standard deviation. Probably because of 
                                 # some very influential countries.
par(mfrow = c(1,1))
plot(PSIS(m8.3, pointwise = T)$k)   # A normally distributed variable might not be appropriate. 

# Our model does not distinguish between the statements: The association between log(gdp)
# and ruggedness depends on being in Afric and :
# The association of being in Africa with log GDP depends upon terrain ruggedness.

# Visualize this association:

rugged_seq <- seq(from = -0.2, to = 1.2, length.out = 30)

muA <- link(m8.3, data = data.frame(cid = 1, rugged_stan = rugged_seq))
muN <- link(m8.3, data = data.frame(cid = 2, rugged_stan = rugged_seq))

delta <- muA - muN

delta.mean <- apply(delta, 2, mean)

par(mfrow = c(1,1 ))
plot(rugged_seq, delta.mean, type = "l", xlab = "Rugged(standardized)", 
     ylab ="E(log(GDP)) of african country - E(log(GDP) non african country")
abline(h = 0, lty = 2)

# Any point above 1 is a nation where log(gdp) is predicted to be higher in africa than 
# in non african countries. When there is very high ruggedness african countries would
# be expected to have higher log(g)



# Continuous interactions ------------------------------------------------------
# When the effect of one parameters depends on another, continuously. 

data(tulips)
d <- tulips
str(d)

# blooms our outcome variable
# Water is soil moisture, (1) low, (3) high.
# shade is exposure, (1) High, (3) low
# Bed indicates clusters of plants 

# Sun and water will help plants grow. But one would be useless without the other, so we expect
# some interaction between the two.

# No interaction model

# preproccessing

d$blooms_std <- d$blooms/max(d$blooms)    # Makes it so 1 is maximum bloom, we don’t want to standardize blooms, because zero
                                          # is a meaningful boundary we want to preserve.
d$water_std <- d$water - mean(d$water)    # centering
d$shade_std <- d$shade - mean(d$shade)

a <- rnorm(1e4, 0.5, 1)
sum(a < 0 | a > 1) / length(a)  # A prior with these values would actually assign 62 % of its probability
                                # mass below zero and above 1.¨

a <- rnorm(1e4, 0.5, 0.25)
sum(a < 0 | a > 1) / length(a)  # Now its only 5 %, much better

a <- rnorm(1e4, 0, 0.25)
sum(a < -0.5 | a > 0.5) / length(a)  # So 95 % of slopes are within reasonable limits


m8.4 <- quap(alist(blooms_std ~ dnorm(mu, sigma),
                   mu <- a + bW * water_std + bS * shade_std,
                   a ~ dnorm(0.5, 0.25),
                   bW ~ dnorm(0, 0.25), 
                   bS ~ dnorm(0, 0.25),
                   sigma ~ dexp(1)),
             data = d)

# With interaction:

m8.5 <- quap(alist(blooms_std ~ dnorm(mu, sigma),
                   mu <- a + bW * water_std + bS * shade_std + bWS * water_std * shade_std,
                   a ~ dnorm(0.5, 0.25),
                   bW ~ dnorm(0, 0.25), 
                   bS ~ dnorm(0, 0.25),
                   bWS ~ dnorm(0, 0.25),
                   sigma ~ dexp(1)),
             data = d)

# Plotting interactions require us to plot the relationship at different levels
# : triptych plot!

par(mfrow = c(1,3))

for(s in -1:1){
  idx <- which(d$shade_std==s)
  plot(d$water_std[idx], d$blooms_std[idx], xlim = c(-1, 1), ylim = c(0, 1),
       xlab = "Water", ylab = "Blooms", pch = 16, col = col.alpha("cyan4", 0.5))
  mu <- link(m8.5, data = data.frame(shade_std = s, water_std = -1:1))
  for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",0.3) )
}

# Shows that with high shade the water is less effective or vice versa, relationship is both ways


# Analyzing our priors --------------------------------------------------------
set.seed(7)

prior <- extract.prior(m8.5)

for(s in -1:1){
  idx <- which(d$shade_std==s)
  plot(d$water_std[idx], d$blooms_std[idx], type = "n",  xlim = c(-1, 1), ylim = c(-1, 2),
       xlab = "Water", ylab = "Blooms", pch = 16, col = col.alpha("cyan4", 0.5))
  mu <- link(m8.5, post = prior, data = data.frame(shade_std = s, water_std = -1:1))
  for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",0.3) )
  lines(-1:1, mu[21, ], col = "black", lwd = 2)
  abline(h = 1, lty = 2)
  abline(h = 0, lty = 2)
}
# Not too bad, but also not great. Better than flat. 

#“These priors contain no bias towards positive or
# negative effects, and at the same time they very weakly bound the effects to 
# realistic ranges.”


# Practice ---------------------------------------------------------------------

# 8E1 Name possible variable with interaction
#1) Rise(bread) and yeast: Temperature
#2) Education and income: Country of origin
#3) Gasoline make car go: Model of car

# 8E2
# Which explanations invoke a interaction?
#  I think they all could, dont know that much about cars, but I would check all of them
# cause it doesnt sound unreasonable that they'd interact

# 8E3
# write linear models

#1) 
Onion cooking result ~ dnorm(mu, sigma)
mu <- a + bH * heat + bT * time + bHT * heat * time

#2)
car speed ~ dnorm(mu, sigma)
mu <- a + bC * cyl + bI * inj 

#3)
leftleaning ~ dnorm(mu, sigma)
mu <- a + bP * parents + bF * friends

# 8M3

year <- seq(1, 50, length.out = 50)
N_wolf <- 200
kW <- 400
rW <- 2.1
wolf_pop <- N_wolf
for(i in 1:49){
  wolf_pop[i+1] <- wolf_pop[i] + (rW * wolf_pop[i] * ((kW - wolf_pop[i])/kW))}

N_raven <- 50
kR <- 600
raven_pop <- N_raven

for(i in 1:49){
  r <- 2 + rnorm(1, mean = wolf_pop[i]/1000, sd = 0.25)
  
  raven_pop[i+1] <- raven_pop[i] + (r * raven_pop[i] * ((kR - raven_pop[i])/kR))
  }

plot(wolf_pop ~ year, type = "l")

# Is not very good to much cycles sorry future me


# 8H1-2

data(tulips)
d <- tulips
# Create ID

d$bid <- ifelse(d$bed == "a", 1,
                ifelse(d$bed == "b", 2, 3))
d$blooms_std <- d$blooms/max(d$blooms)
d$water_std <- d$water - mean(d$water)    # centering
d$shade_std <- d$shade - mean(d$shade)


m8H1.1 <- quap(alist(blooms_std ~ dnorm(mu, sigma),
                   mu <- a[bid] + bW * water_std + bS * shade_std + bWS * water_std * shade_std,
                   a[bid] ~ dnorm(0.5, 0.25),
                   bW ~ dnorm(0, 0.25), 
                   bS ~ dnorm(0, 0.25),
                   bWS ~ dnorm(0, 0.25),
                   sigma ~ dexp(1)),
             data = d)
precis(m8H1.1, depth = 2)
precis(m8.5, depth = 2)


compare(m8.5, m8H1.1, func = WAIC)
plot(compare(m8.5, m8H1.1, func = WAIC))
# Not much difference. as far as predictive ability goes. 

# 8H3

m8.3 <- quap(alist(log_gdp_stan ~ dnorm(mu, sigma),
                   mu <- a[cid] + b[cid]*(rugged_stan - 0.215),
                   a[cid] ~ dnorm(1, 0.1),
                   b[cid] ~ dnorm(0, 0.3),
                   sigma ~ dexp(1))
             , data = dd)

precis(m8.3, depth = 2)


psis <- PSIS(m8.3, pointwise = T)
waic <- WAIC(m8.3, pointwise = T)
plot(psis$k, waic$penalty, pch = 16, col = col.alpha("blue", 0.5))
influential <- which(psis$k > 0.25 & waic$penalty > 0.25)


dd[influential, ]$country
# all three of them are very rugged.

m8.3.2 <- quap(alist(log_gdp_stan ~ dstudent(2, mu, sigma),
                   mu <- a[cid] + b[cid]*(rugged_stan - 0.215),
                   a[cid] ~ dnorm(1, 0.1),
                   b[cid] ~ dnorm(0, 0.3),
                   sigma ~ dexp(1))
             , data = dd)
   

precis(m8.3.2, depth = 2)

psis <- PSIS(m8.3.2, pointwise = T)
waic <- WAIC(m8.3.2, pointwise = T)
points(psis$k, waic$penalty, pch = 16, col = col.alpha("magenta", 0.5))

# Makes the countries with very high influence have less soo, but still does not change the estimates to a large degree. 


# 8H4
data(nettle)
d <- nettle
str(d)
d$lang.per.cap <- with(d, num.lang / k.pop)
d$log.lang.per.cap <- log(d$lang.per.cap)
d$log.lang.per.cap <- (d$log.lang.per.cap - mean(d$log.lang.per.cap))/sd(d$log.lang.per.cap)
d$log.area <- log(d$area)
d$log.area_std <- with(d, (log.area-mean(log.area))/sd(log.area))
                       
d$mean.growing.season_std <- d$mean.growing.season/max(d$mean.growing.season)
range(d$mean.growing.season_std)
range(d$log.area_std)
range(d$log.lang.per.cap)

# a)

m8H4.1 <- quap(alist(log.lang.per.cap ~ dnorm(mu, sigma),
                     mu <- a + bG*(mean.growing.season_std - mean(mean.growing.season_std)) + bA*log.area_std,
                     a ~ dnorm(0, 0.5),
                     bG ~ dnorm(0, 0.5),
                     bA ~ dnorm(0, 0.5),
                     sigma ~ dexp(1)),
               data = d)
# inspect priors

prior <- extract.prior(m8H4.1)
growing_seq <- seq(-0.2, 1.2, length.out = 30)
area_seq <- seq(-4, 4, length.out = 30)

mu <- link(m8H4.1, post = prior, data = data.frame(mean.growing.season_std = growing_seq, log.area_std = area_seq))

plot(NULL, xlim = c(-0.3, 1.3), ylim = c(-4, 4), xlab = "Growing Season", ylab = "Log Language per Capita")
abline(h = min(d$log.lang.per.cap), lty = 2)
abline(h = max(d$log.lang.per.cap), lty = 2)
for(i in 1:50) lines(growing_seq, mu[i,], col = col.alpha("black", 0.2))     
points(mean(d$mean.growing.season_std), mean(d$log.lang.per.cap), pch = 16, col = col.alpha("magenta", 0.4))


# Not too bad, slopes probably many too steep.. 


precis(m8H4.1)
plot(PSIS(m8H4.1, pointwise = T)$k, WAIC(m8H4.1, pointwise = T)$penalty, pch = 16)

# some very influential points, try robust regression


m8H4.2 <- quap(alist(log.lang.per.cap ~ dstudent(2, mu, sigma),
                     mu <- a + bG*(mean.growing.season_std - mean(mean.growing.season_std)) + bA*log.area_std,
                     a ~ dnorm(0, 0.5),
                     bG ~ dnorm(0, 0.5),
                     bA ~ dnorm(0, 0.5),
                     sigma ~ dexp(1)),
               data = d)
points(PSIS(m8H4.2, pointwise = T)$k, WAIC(m8H4.2, pointwise = T)$penalty, col = col.alpha("magenta", 0.3), pch = 16)

precis(m8H4.2)

prior <- extract.prior(m8H4.2)
growing_seq <- seq(-0.2, 1.2, length.out = 30)
area_seq <- seq(-4, 4, length.out = 30)

mu <- link(m8H4.2, post = prior, data = data.frame(mean.growing.season_std = growing_seq, log.area_std = area_seq))

plot(NULL, xlim = c(-0.3, 1.3), ylim = c(-4, 4), xlab = "Growing Season", ylab = "Log Language per Capita")
abline(h = min(d$log.lang.per.cap), lty = 2)
abline(h = max(d$log.lang.per.cap), lty = 2)
for(i in 1:50) lines(growing_seq, mu[i,], col = col.alpha("black", 0.2))     
points(mean(d$mean.growing.season_std), mean(d$log.lang.per.cap), pch = 16, col = col.alpha("magenta", 0.4))

# much better

compare(m8H4.1, m8H4.2, func = PSIS)
plot(precis(m8H4.2))


growing_seq <- seq(-0.2, 1.2, length.out = 30)
area_seq <- seq(-4, 4, length.out = 30)

mu <- link(m8H4.2, data = data.frame(mean.growing.season_std = growing_seq, log.area_std = area_seq))

plot(NULL, xlim = c(-0.3, 1.3), ylim = c(-4, 4), xlab = "Growing Season", ylab = "Log Language per Capita")
abline(h = min(d$log.lang.per.cap), lty = 2)
abline(h = max(d$log.lang.per.cap), lty = 2)
for(i in 1:1000) lines(growing_seq, mu[i,], col = col.alpha("black", 0.2))     
points(mean(d$mean.growing.season_std), mean(d$log.lang.per.cap), pch = 16, col = col.alpha("magenta", 0.4))


# Its possible its positively associated, but not strong evidence for it

# b)

d$sd.growing.season_std <- d$sd.growing.season/max(d$sd.growing.season)


m8H4.3 <- quap(alist(log.lang.per.cap ~ dstudent(2, mu, sigma),
                     mu <- a + bS*(sd.growing.season_std - mean(sd.growing.season_std)) + bA*log.area_std,
                     a ~ dnorm(0, 0.5),
                     bS ~ dnorm(0, 0.5),
                     bA ~ dnorm(0, 0.5),
                     sigma ~ dexp(1)),
               data = d)


prior <- extract.prior(m8H4.3)
sd_seq <- seq(-0.2, 1.2, length.out = 30)
area_seq <- seq(-4, 4, length.out = 30)

mu <- link(m8H4.3, post = prior, data = data.frame(sd.growing.season_std = sd_seq, log.area_std = area_seq))

plot(NULL, xlim = c(-0.3, 1.3), ylim = c(-4, 4), xlab = "Growing Season", ylab = "Log Language per Capita")
abline(h = min(d$log.lang.per.cap), lty = 2)
abline(h = max(d$log.lang.per.cap), lty = 2)
for(i in 1:50) lines(growing_seq, mu[i,], col = col.alpha("black", 0.2))     
points(mean(d$mean.growing.season_std), mean(d$log.lang.per.cap), pch = 16, col = col.alpha("magenta", 0.4))


plot(precis(m8H4.3))
# c)


m8H4.5 <- quap(alist(log.lang.per.cap ~ dstudent(2, mu, sigma),
                     mu <- a + bS*(sd.growing.season_std - mean(sd.growing.season_std)) + 
                          bG*(mean.growing.season_std - mean(mean.growing.season_std)) + bA*log.area_std +
                           bGS*(mean.growing.season_std - mean(mean.growing.season_std))* (sd.growing.season_std - mean(sd.growing.season_std)),
                     a ~ dnorm(-5, 0.5),
                     bS ~ dnorm(0, 0.5),
                     bGS ~ dnorm(0, 0.5),
                     bG ~ dnorm(0, 0.5),
                     bA ~ dnorm(0, 0.5),
                     sigma ~ dexp(1)),
               data = d)
plot(precis(m8H4.5))


# : triptych plot!

par(mfrow = c(1,3))

seq(min(range(d$mean.growing.season_std)), max(range(d$mean.growing.season_std)), length.out = 4)

d$growing.id <- ifelse(d$mean.growing.season_std >= 0 & d$mean.growing.season_std < 0.333, 1, 
                       ifelse(d$mean.growing.season_std >= 0.333 & d$mean.growing.season_std < 0.666, 2, 3))

d$sd.id <- ifelse(d$sd.growing.season_std >= 0 & d$sd.growing.season_std < 0.333, 1, 
                       ifelse(d$sd.growing.season_std >= 0.333 & d$sd.growing.season_std < 0.666, 2, 3))

for(s in 1:3){
  idx <- which(d$growing.id==s)
  plot(d$sd.growing.season_std[idx], d$log.lang.per.cap[idx], xlim = c(0, 1), ylim = c(-4, 4),
       xlab = "SD growing season", ylab = "Language per capita", pch = 16, col = col.alpha("cyan4", 0.5))
  growing_seq <- seq(min(d$mean.growing.season_std[idx]), max(d$mean.growing.season_std[idx]), length.out = length(idx))
  sd_seq <- seq(0, 1, length.out = length(idx))
  area_seq <- seq(-4, 4, length.out = length(idx))
  mu <- link(m8H4.5, data = data.frame(mean.growing.season_std = growing_seq, sd.growing.season_std = sd_seq, log.area_std = area_seq))
  for ( i in 1:20 ) lines( sd_seq , mu[i,] , col=col.alpha("black",0.3) )
  mtext(s)
}



