# Homework week 4

library(rethinking)

#1 WAIC and PSIS on marriage. 
d <- sim_happiness( seed=1977 , N_years=1000 )
precis(d)

d2 <- d[ d$age>17 , ] # only adults
d2$A <- ( d2$age - 18 ) / ( 65 - 18 )  # 0 is now 18 and 1 is 65. 

d2$mid <- d2$married + 1

# With marriage (Collider)
m6.9 <- quap(
  alist(
    happiness ~ dnorm( mu , sigma ),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm( 0 , 1 ),
    bA ~ dnorm( 0 , 2 ),
    sigma ~ dexp(1)
  ) , data=d2 )
precis(m6.9,depth=2)

# Wihtout marriage
m6.10 <- quap( 
               alist(
                 happiness ~ dnorm( mu , sigma ),
                 mu <- a + bA*A,
                 a ~ dnorm( 0 , 1 ),
                 bA ~ dnorm( 0 , 2 ),
                 sigma ~ dexp(1)
               ) , data=d2 )
precis(m6.10)

compare(m6.9, m6.10, func = PSIS)
compare(m6.9, m6.10, func = WAIC)

# The information critera perfer the incorrect model which conditions on a collider. 

#2
data(foxes)
d <- foxes
d$group_std <- with(d, (groupsize - mean(groupsize))/sd(groupsize))
d$avgfood_std <- with(d, (avgfood - mean(avgfood))/sd(avgfood))
d$area_std <- with(d, (area - mean(area))/sd(area))
d$weight_std <- with(d, (weight - mean(weight))/sd(weight))


dat_list <- list(group = d$group_std,
                 food = d$avgfood_std,
                 area = d$area_std,
                 weight = d$weight_std)

m2.1 <- quap(alist(
  weight ~ dnorm(mu, sigma),
  mu <- a + bf * food,
  a ~ dnorm(0, 1),
  bf ~ dnorm(0, 0.5),
  sigma ~ dexp(1)
), data = dat_list)


m2.2 <- quap(alist(
  weight ~ dnorm(mu, sigma),
  mu <- a + bf * food + ba * area + bg * group,
  a ~ dnorm(0, 1)
  bf ~ dnorm(0, 0.5),
  ba ~ dnorm(0, 0.5),
  bg ~ dnorm(0, 0.5),
  sigma ~ dexp(1)
), data = dat_list)

precis(m2.1)
precis(m2.2)

compare(m2.1, m2.2, func = PSIS)

# The model that does the best according to infromation critera estimate direct effects, at best.


