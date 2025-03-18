# Homework week 3
library(rethinking)


data(foxes)
d <- foxes

#1 Causal effect of A

# Generative model
N <- 100
Area <- rnorm(N, mean = 0, sd = 1)
food <- rnorm(N, mean = Area*0.5, sd = 1)

test.d <- data.frame(Area = Area, Food = food)

# test of statistcal model on generative model. 

m1.test <- quap(alist(
  Food ~ dnorm(mu, sigma),
  mu <- Area * ba,
  ba ~ dnorm(0, 1),
  sigma ~ dexp(1)
), data = test.d)
precis(m1.test)

mu <- link(m1.test, data = data.frame(Area = seq(-2, 2, length.out = 30)))

mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

plot(Food ~ Area, data = test.d)
lines(seq(-2, 2, length.out = 30), mu_mean)
shade(mu_PI, seq(-2, 2, length.out = 30))

pred <- sim(m1.test, data = data.frame(Area = seq(-2, 2, length.out = 30)))
pred_PI <- apply(pred, 2, PI)
shade(pred_PI, seq(-2, 2, length.out = 30))

# "Real" data

d$avgfood_std <- with(d, (avgfood - mean(avgfood))/sd(avgfood))
d$area_std <- with(d, (area - mean(area))/sd(area))

m1.1 <- quap(alist(
  avgfood_std ~ dnorm(mu, sigma),
    mu <- area_std * ba,
    ba ~ dnorm(0, 1),
    sigma ~ dexp(1)
), data = d)

precis(m1.1)

area_seq <- seq(-2, 2, length.out = 30)

mu <- link(m1.1, data = data.frame(area_std = area_seq))

mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

plot(avgfood_std ~ area_std, data = d)
lines(area_seq, mu_mean)
shade(mu_PI, area_seq)

pred <- sim(m1.1, data = data.frame(area_std = area_seq))
pred_PI <- apply(pred, 2, PI)
shade(pred_PI, area_seq)

post <- extract.samples(m1.1)
dens(post)

# increasing the area by 1 std deviation will have an average increase of 0.88 sd food.

#2 

d$weight_std <- with(d, (weight-mean(weight)/sd(weight)))

m2.1 <- quap(alist(
  weight_std ~ dnorm(mu, sigma),
  mu <- avgfood_std*bf,
  bf ~ dnorm(0, 1),
  sigma ~ dexp(1)
), data = d)

precis(m2.1)  # The total effect of food is actually very small. 
post <- extract.samples(m2.1)
dens(post)

